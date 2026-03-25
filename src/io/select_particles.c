#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/*
 * select_particles: copy particles that are
 *   (a) inside the view volume, AND
 *   (b) of a type whose bit is set in ptype_mask.
 *
 * ptype_mask is a bitmask: bit i set => keep particle type i.
 * Pass -1 (all bits set) to keep all types regardless of species.
 *
 * The selection is done in-place: kept particles are compacted to the
 * front of the x/y/z/ptype arrays.  *NumPartKeep is set to the number
 * of particles that passed.
 */
void select_particles(float *x, float *y, float *z, int *ptype,
                      double BoxSize, long long NumPart,
                      float xc, float yc, float zc, float lbox,
                      int ptype_mask, long long *NumPartKeep)
{
    long long i;
    int j;
    int ix, iy, iz;
    float dx[3];

    fprintf(stdout, "BoxSize: %lf\n", BoxSize);
    fprintf(stdout, "Particle type mask: 0x%x\n", (unsigned)ptype_mask);
    fflush(stdout);

    *NumPartKeep = 0;
    
    fprintf(stdout, "NumPart: %llu\n", NumPart);
    fflush(stdout);

    for (i = 0; i < NumPart; i++) {
        /* Skip if this type is not requested */
        if (ptype_mask != -1 && !(ptype_mask & (1 << ptype[i])))
            continue;

        dx[0] = x[i] - xc;
        dx[1] = y[i] - yc;
        dx[2] = z[i] - zc;

        /* Apply periodic wrap */
        for (j = 0; j < 3; j++) {
            if (dx[j] >  0.5 * BoxSize) dx[j] -= BoxSize;
            if (dx[j] < -0.5 * BoxSize) dx[j] += BoxSize;
        }

        /* Keep only particles inside the view cube */
        ix = (int)(2.0f * dx[0] / lbox);
        iy = (int)(2.0f * dx[1] / lbox);
        iz = (int)(2.0f * dx[2] / lbox);

        if (ix == 0 && iy == 0 && iz == 0) {
            x    [*NumPartKeep] = x    [i];
            y    [*NumPartKeep] = y    [i];
            z    [*NumPartKeep] = z    [i];
            ptype[*NumPartKeep] = ptype[i];
            (*NumPartKeep)++;
        }
    }

    fprintf(stdout, "Selected %lld particles (mask 0x%x)\n",
            *NumPartKeep, (unsigned)ptype_mask);
    fflush(stdout);
}
