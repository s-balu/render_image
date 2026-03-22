/*
 * deposit_sph_2d.c — direct 2D SPH projection deposit
 *
 * Compiled as a separate translation unit so the compiler can inline
 * proj2_lookup() into the ix inner loop and auto-vectorise with SIMD.
 * When this code lives inside the large smooth_to_mesh.c, the compiler
 * makes conservative inlining decisions that prevent vectorisation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include "header.h"

#ifdef KERNEL_SMOOTHING

/* deposit_particle_t and deposit_sph_2d() are declared in header.h */

/* Local projected kernel table -- populated by deposit_sph_2d_init() */
#define D2D_NTAB 2048
static float d2d_proj[D2D_NTAB];
static int   d2d_ready = 0;

void deposit_sph_2d_init(void)
{
    if (d2d_ready) return;
    for (int i = 0; i < D2D_NTAB; i++) {
        float q2 = 4.0f * (float)i / (float)(D2D_NTAB - 1);
        d2d_proj[i] = get_projected_kernel_value(sqrtf(q2));
    }
    d2d_proj[D2D_NTAB - 1] = 0.0f;
    d2d_ready = 1;
}

/* Inline lookup -- compiler CAN inline this since it's in the same TU */
static inline float d2d_lookup(float q2)
{
    if (q2 >= 4.0f) return 0.0f;
    float fi = q2 * (float)(D2D_NTAB - 1) * 0.25f;
    int   i0 = (int)fi;
    return d2d_proj[i0] + (fi - (float)i0) * (d2d_proj[i0+1] - d2d_proj[i0]);
}

#ifndef MAX_H_DEPOSIT_PX
#define MAX_H_DEPOSIT_PX 8
#endif

long long deposit_sph_2d(long long n_packed,
                          const deposit_particle_t *particles,
                          int full_width, int height,
                          float *data2d)
{
    long long n_deposited = 0;
    memset(data2d, 0, (size_t)full_width * height * sizeof(float));

#ifdef ENABLE_OPENMP
    int NThreads = omp_get_max_threads();
    printf("deposit_sph_2d: %d threads, AoS packed\n", NThreads);
    fflush(stdout);
    #pragma omp parallel num_threads(NThreads) reduction(+:n_deposited)
    {
        int tid  = omp_get_thread_num();
        int y_lo = (long long)tid       * height / NThreads;
        int y_hi = (long long)(tid + 1) * height / NThreads;
        int n_rows = y_hi - y_lo;
        float *local = (float *)calloc((size_t)n_rows * full_width, sizeof(float));
        if (!local) goto d2d_done;
#else
    {
        int y_lo = 0, y_hi = height, n_rows = height;
        float *local = data2d;
#endif
        /* dx2 must hold up to 2*(2*MAX_H_DEPOSIT_PX)+1 entries.
         * At MAX_H_DEPOSIT_PX=41px: dnxy = ceil(2*41) = 82, nx <= 2*82+1 = 165.
         * Use 512 as a safe upper bound for any realistic cap value. */
        float dx2[512];

        for (long long ii = 0; ii < n_packed; ii++) {
            float px    = particles[ii].px;
            float py    = particles[ii].py;
            float dh_px = particles[ii].h;
            float sup   = 2.0f * dh_px;

            if (py + sup < (float)y_lo) continue;
            if (py - sup >= (float)y_hi) continue;
            if (px + sup < 0.0f || px - sup >= (float)full_width) continue;

            int dnxy = (int)ceilf(sup);
            float h_inv = 1.0f / dh_px;
            float norm  = h_inv * h_inv;

            int iy0 = (int)py - dnxy; if (iy0 < y_lo) iy0 = y_lo;
            int iy1 = (int)py + dnxy + 1; if (iy1 > y_hi) iy1 = y_hi;
            int ix0 = (int)px - dnxy; if (ix0 < 0) ix0 = 0;
            int ix1 = (int)px + dnxy + 1; if (ix1 > full_width) ix1 = full_width;
            int nx = ix1 - ix0; if (nx <= 0) continue;
            if (nx > 512) nx = 512;  /* guard: should never trigger with sane h cap */

            for (int k = 0; k < nx; k++) {
                float dx = (float)(ix0 + k) + 0.5f - px;
                dx2[k] = dx * dx * h_inv * h_inv;
            }

            for (int iy = iy0; iy < iy1; iy++) {
                float dy  = (float)iy + 0.5f - py;
                float dy2 = dy * dy * h_inv * h_inv;
                if (dy2 >= 4.0f) continue;
                float *row = local + (size_t)(iy - y_lo) * full_width + ix0;
                /* This ix loop is the hot path. With d2d_lookup inlined and
                 * dx2[] precomputed, the compiler should auto-vectorise with SIMD. */
                for (int k = 0; k < nx; k++) {
                    row[k] += d2d_lookup(dx2[k] + dy2) * norm;
                    n_deposited++;
                }
            }
        }

#ifdef ENABLE_OPENMP
        for (int r = y_lo; r < y_hi; r++) {
            float *src = local + (size_t)(r - y_lo) * full_width;
            float *dst = data2d + (size_t)r * full_width;
            for (int c = 0; c < full_width; c++) dst[c] += src[c];
        }
        free(local);
        d2d_done:;
    } /* end omp parallel */
#else
    }
#endif

    return n_deposited;
}

#endif /* KERNEL_SMOOTHING */
