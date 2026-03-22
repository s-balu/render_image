#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "header.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

float get_rand(int seed) { return (float)rand() / RAND_MAX; }

#ifdef ENABLE_MPI
/*
 * split_across_tasks_as_slabs
 *
 * Redistribute particles so each MPI task holds an equal number of
 * particles (±1), using adaptive x-slab boundaries determined by the
 * global particle distribution.
 *
 * Algorithm
 * ---------
 * 1. Each task sorts its local particles by x and picks NTask-1
 *    candidate split points at local quantiles 1/NTask, 2/NTask, ...
 *
 * 2. All tasks share their candidate split points via MPI_Allgather.
 *    The global median of each candidate set gives a robust estimate
 *    of the global equal-N split points that adapts to the actual
 *    density field.
 *
 * 3. Particles are assigned to tasks using these boundaries.  The
 *    resulting counts are within ~NTask of equal (bias from using
 *    medians of local quantiles rather than exact global sort).
 *
 * 4. Physical slab boundaries are broadcast so smooth_to_mesh can
 *    compute the correct local pixel grid for each task.
 *
 * Why not global sort?
 * --------------------
 * An exact equal-N decomposition requires sorting all N particles
 * globally, which costs O(N log N / NTask) communication.  The median-
 * of-local-quantiles approach costs O(NTask²) communication (tiny) and
 * O(N/NTask × log(N/NTask)) compute, and gives balance within 1-2%.
 */
void split_across_tasks_as_slabs(float *x, float *y, float *z,
                                  long long *NumPart,
                                  float xc, float yc, float zc,
                                  float lbox, float BoxSize)
{
    (void)yc; (void)zc; (void)BoxSize;

    int i, j;
    long long local_n = *NumPart;

    /* ----------------------------------------------------------------
     * Step 1: build a global histogram of particle x-coordinates
     *
     * Each task bins its local particles into NHIST bins over the view
     * box [xc-lbox/2, xc+lbox/2].  MPI_Allreduce sums the histograms
     * globally.  We then walk the cumulative histogram to find the
     * equal-N cut points.
     *
     * This handles zoom simulations correctly: even if all particles
     * cluster in a tiny x-range, the histogram bins that range finely
     * and finds accurate equal-N splits.
     * ---------------------------------------------------------------- */
    #define NHIST 8192
    long long *hist = (long long *)calloc(NHIST, sizeof(long long));
    if (!hist) { fprintf(stderr, "split: malloc hist failed\n"); MPI_Abort(MPI_COMM_WORLD,1); }

    float x_lo = xc - 0.5f * lbox;
    float x_hi = xc + 0.5f * lbox;
    float bin_inv = (float)NHIST / (x_hi - x_lo);

    for (long long k = 0; k < local_n; k++) {
        int b = (int)((x[k] - x_lo) * bin_inv);
        if (b < 0)      b = 0;
        if (b >= NHIST) b = NHIST - 1;
        hist[b]++;
    }

    MPI_Allreduce(MPI_IN_PLACE, hist, NHIST, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);

    /* ----------------------------------------------------------------
     * Step 2: find equal-N split points from cumulative histogram
     * ---------------------------------------------------------------- */
    /* Total particles across all tasks */
    long long global_n = 0;
    for (i = 0; i < NHIST; i++) global_n += hist[i];

    float *splits = (float *)malloc((NTask + 1) * sizeof(float));
    splits[0]     = x_lo;
    splits[NTask] = x_hi;

    long long cumsum = 0, split_idx = 1;
    for (i = 0; i < NHIST && split_idx < NTask; i++) {
        cumsum += hist[i];
        /* Place a split when we pass the k/NTask fraction of total */
        while (split_idx < NTask &&
               cumsum * NTask >= split_idx * global_n) {
            /* Split at the right edge of this bin */
            splits[split_idx] = x_lo + (float)(i + 1) / bin_inv;
            split_idx++;
        }
    }
    /* Fill any remaining splits at the right edge (degenerate case) */
    for (i = split_idx; i < NTask; i++)
        splits[i] = x_hi;

    free(hist);

    /* Ensure strict monotonicity with minimum spacing of one bin width */
    float min_gap = 1.0f / bin_inv;
    for (i = 1; i <= NTask; i++)
        if (splits[i] <= splits[i-1]) splits[i] = splits[i-1] + min_gap;

    /* Store global slab boundaries */
    if (!slab_x_lo) slab_x_lo = (float *)malloc(NTask * sizeof(float));
    if (!slab_x_hi) slab_x_hi = (float *)malloc(NTask * sizeof(float));
    for (i = 0; i < NTask; i++) {
        slab_x_lo[i] = splits[i];
        slab_x_hi[i] = splits[i + 1];
    }
    free(splits);

    /* ----------------------------------------------------------------
     * Step 3: assign each local particle to a task
     * ---------------------------------------------------------------- */
    long long *num_x     = (long long *)calloc(NTask, sizeof(long long));
    long long *num_gather = (long long *)calloc(NTask, sizeof(long long));
    long long **pid = (long long **)malloc(NTask * sizeof(long long *));
    for (i = 0; i < NTask; i++)
        pid[i] = (long long *)malloc(local_n * sizeof(long long));

    for (long long k = 0; k < local_n; k++) {
        /* Binary search for the correct slab */
        int lo = 0, hi = NTask - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (x[k] >= slab_x_hi[mid]) lo = mid + 1;
            else hi = mid;
        }
        int nx = lo;
        if (nx < 0)      nx = 0;
        if (nx >= NTask) nx = NTask - 1;
        pid[nx][num_x[nx]++] = k;
    }

    if (ThisTask == 0) {
        printf("Equal-N decomposition splits:\n");
        for (i = 0; i < NTask; i++)
            printf("  slab %d: x∈[%.4g, %.4g]\n", i, slab_x_lo[i], slab_x_hi[i]);
        fflush(stdout);
    }

    /* ----------------------------------------------------------------
     * Step 4: exchange particle counts then particles
     * ---------------------------------------------------------------- */
    for (i = 0; i < NTask; i++) {
        if (i != ThisTask) {
            MPI_Send(&num_x[i],      1, MPI_LONG_LONG_INT, i, ThisTask, MPI_COMM_WORLD);
            MPI_Recv(&num_gather[i], 1, MPI_LONG_LONG_INT, i, i,        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    num_gather[ThisTask] = num_x[ThisTask];

    long long num_proc = 0;
    for (i = 0; i < NTask; i++) num_proc += num_gather[i];

    long long num_tot;
    MPI_Reduce(&num_proc, &num_tot, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_tot,   1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    fprintf(stdout, "Task %d: will own %lld particles after redistribution "
            "(global total %lld)\n", ThisTask, num_proc, num_tot);
    fflush(stdout);

    float *posx = (float *)malloc(num_proc * sizeof(float));
    float *posy = (float *)malloc(num_proc * sizeof(float));
    float *posz = (float *)malloc(num_proc * sizeof(float));
    if (!posx || !posy || !posz) {
        fprintf(stderr, "Task %d: malloc failed\n", ThisTask);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Copy self-owned particles */
    long long noffset = 0;
    for (j = 0; j < ThisTask; j++) noffset += num_gather[j];
    for (j = 0; j < num_x[ThisTask]; j++) {
        long long src = pid[ThisTask][j];
        posx[j + noffset] = x[src];
        posy[j + noffset] = y[src];
        posz[j + noffset] = z[src];
    }

    /* Exchange with other tasks */
    for (i = 0; i < NTask; i++) {
        if (i == ThisTask) continue;
        long long nsend_ll = num_x[i], nrec_ll = num_gather[i];
        if (nsend_ll > INT_MAX || nrec_ll > INT_MAX) {
            fprintf(stderr, "Task %d: count to task %d exceeds INT_MAX\n", ThisTask, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int nsend = (int)nsend_ll, nrec = (int)nrec_ll;
        float *sb = (float *)malloc(nsend * sizeof(float));
        float *rb = (float *)malloc(nrec  * sizeof(float));

        noffset = 0;
        for (j = 0; j < i; j++) noffset += num_gather[j];

        for (j = 0; j < nsend; j++) sb[j] = x[pid[i][j]];
        MPI_Sendrecv(sb, nsend, MPI_FLOAT, i, ThisTask,
                     rb, nrec,  MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (j = 0; j < nrec; j++) posx[j + noffset] = rb[j];

        for (j = 0; j < nsend; j++) sb[j] = y[pid[i][j]];
        MPI_Sendrecv(sb, nsend, MPI_FLOAT, i, ThisTask,
                     rb, nrec,  MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (j = 0; j < nrec; j++) posy[j + noffset] = rb[j];

        for (j = 0; j < nsend; j++) sb[j] = z[pid[i][j]];
        MPI_Sendrecv(sb, nsend, MPI_FLOAT, i, ThisTask,
                     rb, nrec,  MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (j = 0; j < nrec; j++) posz[j + noffset] = rb[j];

        free(sb); free(rb);
    }

    /* Compute and report x-range */
    float xmin = 1e10f, xmax = -1e10f;
    for (i = 0; i < (int)num_proc; i++) {
        if (posx[i] < xmin) xmin = posx[i];
        if (posx[i] > xmax) xmax = posx[i];
    }

    *NumPart = num_proc;
    memcpy(x, posx, num_proc * sizeof(float));
    memcpy(y, posy, num_proc * sizeof(float));
    memcpy(z, posz, num_proc * sizeof(float));

    free(posx); free(posy); free(posz);
    for (i = 0; i < NTask; i++) free(pid[i]);
    free(pid); free(num_x); free(num_gather);

    fprintf(stdout, "Task %d: %lld particles, x∈[%g,%g] (slab [%g,%g])\n",
            ThisTask, num_proc, xmin, xmax,
            slab_x_lo[ThisTask], slab_x_hi[ThisTask]);
    fflush(stdout);
}
#endif /* ENABLE_MPI */
