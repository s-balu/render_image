#include <stdlib.h>
#include <string.h>
#include "header.h"
#include <png.h>

/* Global variable definitions (declared extern in header.h) */
int ThisTask, NTask;
float *slab_x_lo = NULL;
float *slab_x_hi = NULL;
int SnapFormat;

/* ------------------------------------------------------------------ */
/* sim_info lifecycle helpers                                           */
/* ------------------------------------------------------------------ */

/*
 * sim_info_alloc: allocate the three dynamic arrays inside *h for
 * `num_types` particle types and zero-initialise them.
 * Call sim_info_free() when done.
 */
void sim_info_alloc(struct sim_info *h, int num_types)
{
    h->num_types = num_types;
    h->npart     = (int    *)calloc(num_types, sizeof(int));
    h->nall      = (int    *)calloc(num_types, sizeof(int));
    h->massarr   = (double *)calloc(num_types, sizeof(double));
}

void sim_info_free(struct sim_info *h)
{
    free(h->npart);
    free(h->nall);
    free(h->massarr);
    h->npart   = NULL;
    h->nall    = NULL;
    h->massarr = NULL;
    h->num_types = 0;
}
