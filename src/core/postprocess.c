/*
 * postprocess.c — per-frame density map post-processing.
 *
 * Extracted from render_image.c.  See postprocess.h for a description
 * of the two steps (noise floor + CIC hole inpainting).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "postprocess.h"

void postprocess_frame(float *data, int nx, int ny, float noise_floor_abs)
{
    int npix = nx * ny;

    /* -------------------------------------------------------------- */
    /* Step 1: Noise floor                                              */
    /*                                                                  */
    /* Zero pixels below NOISE_FLOOR_FRACTION of the peak value.       */
    /* Must be done after MPI_Reduce(SUM) — not per-task — because     */
    /* each task has a different local dmax, so per-task thresholds     */
    /* create visible stripe artefacts at slab boundaries.              */
    /* -------------------------------------------------------------- */
    float gdmax = 0.0f;
    for (int p = 0; p < npix; p++)
        if (data[p] > gdmax) gdmax = data[p];

    if (noise_floor_abs > 0.0f)
        for (int p = 0; p < npix; p++)
            if (data[p] < noise_floor_abs) data[p] = 0.0f;

    /* -------------------------------------------------------------- */
    /* Step 2: CIC zero-hole inpainting                                 */
    /*                                                                  */
    /* With ~0.03 particles/voxel on a 768³ grid the CIC 2×2×2 stencil */
    /* leaves contiguous empty regions several pixels wide after        */
    /* projection.  Iterative 3×3 propagation expands the filled        */
    /* frontier by one pixel per pass.  Only zero pixels are written;   */
    /* non-zero pixels are never modified, preserving dynamic range.    */
    /* -------------------------------------------------------------- */
    float *buf = (float *)malloc((size_t)npix * sizeof(float));
    if (!buf) {
        fprintf(stderr, "postprocess_frame: malloc failed for inpainting buffer\n");
        return;   /* skip inpainting; noise floor was still applied */
    }

    memcpy(buf, data, (size_t)npix * sizeof(float));

    int changed = 1;
    for (int pass = 0; pass < MAX_FILL_PASSES && changed; pass++) {
        changed = 0;
        for (int row = 1; row < ny - 1; row++) {
            for (int col = 1; col < nx - 1; col++) {
                int pidx = row * nx + col;
                if (buf[pidx] > 0.0f) continue;   /* already filled */

                float sum = 0.0f;
                int   nn  = 0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++) {
                        float v = buf[(row + dy) * nx + (col + dx)];
                        if (v > 0.0f) { sum += v; nn++; }
                    }

                if (nn > 0) {
                    data[pidx] = sum / (float)nn;
                    changed    = 1;
                }
            }
        }
        /* Sync buf so the next pass sees this pass's fills */
        memcpy(buf, data, (size_t)npix * sizeof(float));
    }

    free(buf);
}
