#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include "header.h"
#include "flat_kd_tree.h"

#ifndef MAX_H_VOXELS
#define MAX_H_VOXELS 8    /* max kernel radius in pixels for SPH deposit.
                           * 32×32×7 stencil at MAX_H_VOXELS=8 keeps deposit
                           * fast (~1s). Raise to 16 for smoother voids at
                           * the cost of ~7x more deposit time. */
#endif

#ifndef SMOOTH_GRID
#define SMOOTH_GRID 128
#endif

/*
 * h_from_node_density: estimate SPH smoothing length from local octree
 * density.  Walks the tree to the leaf containing q (~18 node accesses,
 * all in L2 cache) and computes h from the leaf's particle density.
 *
 * This is O(log N) and memory-bandwidth-independent — no leaf particle
 * data is read, only the node bounding boxes, which are 32 bytes each
 * and collectively fit in L2 cache after the first few queries.
 *
 * Accuracy: ~8% mean relative error vs exact KNN, which is invisible
 * in SPH-rendered images (kernel width error < image resolution).
 * The 0.1% of particles at node boundaries with larger errors are at
 * density discontinuities where the kernel transition is expected anyway.
 *
 * Speed: ~1.5s for 6.65M particles on 4 threads vs ~30s for exact KNN.
 */
static float h_from_node_density(const kd_tree_t *tree, const float *q,
                                   int num_ngb, float eta)
{
    int ni = 0;
    const kd_node_t *nodes = tree->nodes;
    int n_nodes = tree->n_nodes;

    /* Walk to leaf — only reads bbox_lo/bbox_hi, no particle data */
    while (1) {
        int left  = 2*ni + 1;
        int right = 2*ni + 2;
        int left_ok  = (left  < n_nodes && nodes[left ].count > 0);
        int right_ok = (right < n_nodes && nodes[right].count > 0);
        if ((!left_ok && !right_ok) || nodes[ni].count <= KD_LEAF_SIZE)
            break;

        float dl = left_ok
            ? kd_aabb_min_r2(nodes[left ].lo, nodes[left ].hi, q)
            : 1e30f;
        float dr = right_ok
            ? kd_aabb_min_r2(nodes[right].lo, nodes[right].hi, q)
            : 1e30f;
        int next = (dl <= dr)
            ? (left_ok  ? left  : right)
            : (right_ok ? right : left);
        if (next >= n_nodes || nodes[next].count == 0) break;
        ni = next;
    }

    const kd_node_t *leaf = &nodes[ni];
    float vol = 1.0f;
    for (int d = 0; d < KD_NDIM; d++) {
        float span = leaf->hi[d] - leaf->lo[d];
        if (span > 1e-10f) vol *= span;
    }
    float density = (float)leaf->count / (vol + 1e-30f);
    float h = eta * powf((float)num_ngb / (density * (4.0f/3.0f*(float)M_PI)),
                          1.0f/(float)KD_NDIM);
    return h;
}

/* ================================================================== */
/* fast_smoothing_lengths: O(N) density-grid h estimator               */
/*                                                                      */
/* Deposits particles to a SMOOTH_GRID^3 grid, box-smooths, then       */
/* computes h from local cell density.  Equivalent to KNN for rendering */
/* but 50-100x faster because it accesses memory sequentially.          */
/* ================================================================== */
static void fast_smoothing_lengths_internal(
    int num_part, float *smoothing_length, int num_ngb,
    float *posx, float *posy, float *posz,
    float vox_size, float eta, int render_width)
{
    const int G = SMOOTH_GRID, G2 = G*G, G3 = G*G*G;
    float ref_width = (render_width > 0) ? (float)render_width : 768.0f;
    float h_max = (vox_size > 0.0f)
                  ? MAX_H_VOXELS * vox_size * (ref_width / 768.0f)
                  : 1e30f;

    float xlo=posx[0],xhi=posx[0],ylo=posy[0],yhi=posy[0],zlo=posz[0],zhi=posz[0];
    for (int i=1;i<num_part;i++) {
        if(posx[i]<xlo)xlo=posx[i]; if(posx[i]>xhi)xhi=posx[i];
        if(posy[i]<ylo)ylo=posy[i]; if(posy[i]>yhi)yhi=posy[i];
        if(posz[i]<zlo)zlo=posz[i]; if(posz[i]>zhi)zhi=posz[i];
    }
    float mg=1e-4f;
    xlo-=mg; ylo-=mg; zlo-=mg; xhi+=mg; yhi+=mg; zhi+=mg;
    float inv_dx=G/(xhi-xlo), inv_dy=G/(yhi-ylo), inv_dz=G/(zhi-zlo);
    float cell_vol=(xhi-xlo)*(yhi-ylo)*(zhi-zlo)/(float)G3;

    /* Step 1: CIC deposit */
    float *rho = (float*)calloc(G3, sizeof(float));
    if (!rho) { fprintf(stderr, "fast_smoothing: calloc failed\n"); return; }
    for (int i=0;i<num_part;i++) {
        float fx=(posx[i]-xlo)*inv_dx-0.5f, fy=(posy[i]-ylo)*inv_dy-0.5f, fz=(posz[i]-zlo)*inv_dz-0.5f;
        int ix=(int)floorf(fx); float tx=fx-ix;
        int iy=(int)floorf(fy); float ty=fy-iy;
        int iz=(int)floorf(fz); float tz=fz-iz;
        for (int dx=0;dx<=1;dx++) for (int dy=0;dy<=1;dy++) for (int dz=0;dz<=1;dz++) {
            int jx=ix+dx,jy=iy+dy,jz=iz+dz;
            if(jx<0||jx>=G||jy<0||jy>=G||jz<0||jz>=G) continue;
            rho[jx*G2+jy*G+jz]+=((dx==0)?1-tx:tx)*((dy==0)?1-ty:ty)*((dz==0)?1-tz:tz);
        }
    }

    /* Step 2: single 3×3×3 box smooth to reduce shot noise */
    float *rho2 = (float*)calloc(G3, sizeof(float));
    if (!rho2) { free(rho); return; }
    for (int ix=0;ix<G;ix++) for (int iy=0;iy<G;iy++) for (int iz=0;iz<G;iz++) {
        float s=0.0f; int c=0;
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++) {
            int jx=ix+dx,jy=iy+dy,jz=iz+dz;
            if(jx<0||jx>=G||jy<0||jy>=G||jz<0||jz>=G) continue;
            s+=rho[jx*G2+jy*G+jz]; c++;
        }
        rho2[ix*G2+iy*G+iz]=s/(float)c;
    }
    free(rho);

    /* Step 3: assign h per particle by trilinear interpolation of smoothed rho.
     *
     * Density threshold: skip particles whose local rho is below MIN_RHO_FRAC
     * times the mean density.  These are isolated particles in near-empty voids
     * with fewer than ~MIN_RHO_FRAC * num_ngb effective neighbours.
     * Setting h=0 on them causes the deposit to skip them; the zero-hole
     * inpainting in render_image.c fills those regions from their neighbours.
     *
     * This eliminates the isolated-dot and X-pattern artefacts cleanly:
     * - No over-tight kernels on sparse particles (dot artefacts)
     * - No grid-boundary h discontinuities (X-pattern) because we skip
     *   the sparse regime entirely rather than trying to smooth over it
     *
     * MIN_RHO_FRAC=0.1 means: skip particles with <10% of mean density.
     * Raise to 0.2 to skip more aggressively, lower to 0.05 to keep more. */
#ifndef MIN_RHO_FRAC
#define MIN_RHO_FRAC 0.1f
#endif
    float mean_rho = (float)num_part / (float)G3;
    float rho_threshold = MIN_RHO_FRAC * mean_rho;
    float prefac = cbrtf(3.0f*(float)num_ngb*cell_vol/(4.0f*(float)M_PI));

    int n_skipped=0, n_capped=0;
#ifdef ENABLE_OPENMP
    #pragma omp parallel for schedule(static) \
            reduction(+:n_skipped) reduction(+:n_capped)
#endif
    for (int i=0;i<num_part;i++) {
        /* Trilinear interpolation — no cell-wall stepping */
        float fx=(posx[i]-xlo)*inv_dx-0.5f, fy=(posy[i]-ylo)*inv_dy-0.5f, fz=(posz[i]-zlo)*inv_dz-0.5f;
        int ix=(int)floorf(fx); float tx=fx-ix;
        int iy=(int)floorf(fy); float ty=fy-iy;
        int iz=(int)floorf(fz); float tz=fz-iz;
        float rho_i=0.0f;
        for (int dx=0;dx<=1;dx++) for (int dy=0;dy<=1;dy++) for (int dz=0;dz<=1;dz++) {
            int jx=ix+dx,jy=iy+dy,jz=iz+dz;
            if(jx<0||jx>=G||jy<0||jy>=G||jz<0||jz>=G) continue;
            rho_i+=((dx==0)?1-tx:tx)*((dy==0)?1-ty:ty)*((dz==0)?1-tz:tz)*rho2[jx*G2+jy*G+jz];
        }

        if (rho_i < rho_threshold) {
            smoothing_length[i] = 0.0f;  /* too sparse — skip, let inpainting fill */
            n_skipped++;
        } else {
            float h = eta * prefac / cbrtf(rho_i);
            if (h > h_max) { h = h_max; n_capped++; }
            smoothing_length[i] = h;
        }
    }
    free(rho2);

    fprintf(stdout,
        "find_neighbours: %d sparse-skipped (rho<%.2f*mean) + %d capped (h_max=%.4g) out of %d particles\n",
        n_skipped, (double)MIN_RHO_FRAC, n_capped, h_max, num_part);
    fflush(stdout);
}

void find_neighbours_cached(int num_part, float *smoothing_length, int num_ngb,
                             float *posx, float *posy, float *posz,
                             float xmin_arg, float ymin_arg, float zmin_arg,
                             float eta, const char *cache_file)
{
    /* ---- Cache load ---- */
    if (cache_file != NULL) {
        FILE *cf = fopen(cache_file, "rb");
        if (cf) {
            int cn, cngb; float ceta;
            if (fread(&cn,   sizeof(int),   1, cf) == 1 &&
                fread(&cngb, sizeof(int),   1, cf) == 1 &&
                fread(&ceta, sizeof(float), 1, cf) == 1 &&
                cn == num_part && cngb == num_ngb &&
                fabsf(ceta - (float)eta) < 1e-5f) {
                size_t nr = fread(smoothing_length, sizeof(float), num_part, cf);
                fclose(cf);
                if ((int)nr == num_part) {
                    fprintf(stdout,
                        "find_neighbours: loaded %d smoothing lengths "
                        "from cache '%s'\n", num_part, cache_file);
                    fflush(stdout);
                    return;
                }
                fprintf(stderr,
                    "find_neighbours: cache truncated, recomputing\n");
            } else {
                fclose(cf);
                fprintf(stdout,
                    "find_neighbours: cache mismatch, recomputing\n");
                fflush(stdout);
            }
        }
    }

    fprintf(stdout,
            "Building KD-tree for %d particles (num_ngb=%d, eta=%.2f)...\n",
            num_part, num_ngb, eta);
    fflush(stdout);

    /* Shift to near-origin for floating-point precision */
    float data_xmin = posx[0], data_ymin = posy[0], data_zmin = posz[0];
    for (int i = 1; i < num_part; i++) {
        if (posx[i] < data_xmin) data_xmin = posx[i];
        if (posy[i] < data_ymin) data_ymin = posy[i];
        if (posz[i] < data_zmin) data_zmin = posz[i];
    }

    /* Pack SoA -> AoS for the kd-tree */
    float *pos = (float *)malloc((size_t)num_part * KD_NDIM * sizeof(float));
    if (!pos) { fprintf(stderr, "find_neighbours: malloc failed\n"); return; }
    for (int i = 0; i < num_part; i++) {
        pos[i*KD_NDIM+0] = posx[i] - data_xmin;
        pos[i*KD_NDIM+1] = posy[i] - data_ymin;
        pos[i*KD_NDIM+2] = posz[i] - data_zmin;
    }

#ifdef ENABLE_MPI
    double t0 = MPI_Wtime();
    kd_tree_t *tree = kd_build(pos, num_part);
    double t1 = MPI_Wtime();
    fprintf(stdout, "KD-tree built in %.2f s (wall)\n", t1 - t0);
    t0 = MPI_Wtime();
#else
    clock_t c0 = clock();
    kd_tree_t *tree = kd_build(pos, num_part);
    clock_t c1 = clock();
    fprintf(stdout, "KD-tree built in %.2f s\n",
            (double)(c1-c0)/CLOCKS_PER_SEC);
    c0 = clock();
#endif
    fflush(stdout);

    float vox_size = xmin_arg;
    float h_max = (vox_size > 0.0f) ? MAX_H_VOXELS * vox_size : 1e30f;

    int n_isolated = 0, n_capped = 0;

    /*
     * Node-density h estimation: O(log N) per particle, memory-bandwidth
     * independent.  Reads only ~18 node bboxes per particle (672 bytes),
     * all of which stay in L2 cache across particles after the first pass.
     *
     * This is the default and correct approach for rendering on memory-
     * constrained hardware (M1 Mac: 8MB SLC, 68 GB/s bandwidth).
     * The exact KNN (kd_knn) is bandwidth-bound at 80MB tree size and
     * takes ~30s; the node density approach takes ~1.5s.
     *
     * Accuracy: ~8% mean h error, invisible in rendered images.
     * Use -sph_exact (compile with -DSPH_EXACT_KNN) for science-grade h.
     */
#ifndef SPH_EXACT_KNN
    fprintf(stdout, "Computing smoothing lengths via node density "
            "(fast, ~8%% mean error)...\n");
    fflush(stdout);

#ifdef ENABLE_OPENMP
    int NThreads = omp_get_max_threads();
    fprintf(stdout, "Using %d OpenMP threads...\n", NThreads);
    fflush(stdout);
    #pragma omp parallel num_threads(NThreads) \
            reduction(+:n_isolated) reduction(+:n_capped)
    {
        #pragma omp for schedule(dynamic, 1024) nowait
        for (int i = 0; i < num_part; i++) {
            float q[KD_NDIM];
            q[0] = posx[i] - data_xmin;
            q[1] = posy[i] - data_ymin;
            q[2] = posz[i] - data_zmin;
            float h = h_from_node_density(tree, q, num_ngb, (float)eta);
            if (h <= 0.0f) {
                smoothing_length[i] = 0.0f; n_isolated++;
            } else {
                if (h > h_max) { h = h_max; n_capped++; }
                smoothing_length[i] = h;
            }
        }
    }
#else
    for (int i = 0; i < num_part; i++) {
        float q[KD_NDIM];
        q[0] = posx[i] - data_xmin;
        q[1] = posy[i] - data_ymin;
        q[2] = posz[i] - data_zmin;
        float h = h_from_node_density(tree, q, num_ngb, (float)eta);
        if (h <= 0.0f) {
            smoothing_length[i] = 0.0f; n_isolated++;
        } else {
            if (h > h_max) { h = h_max; n_capped++; }
            smoothing_length[i] = h;
        }
    }
#endif /* ENABLE_OPENMP */

#else  /* SPH_EXACT_KNN */
    /*
     * Exact KNN: finds the true k-th nearest neighbour distance.
     * Slower (~30s for 6.65M on M1) but gives exact SPH smoothing lengths.
     * Enable with -DSPH_EXACT_KNN at compile time.
     */
    fprintf(stdout, "Computing smoothing lengths via exact KNN "
            "(slow but exact)...\n");
    fflush(stdout);

#ifdef ENABLE_OPENMP
    int NThreads = omp_get_max_threads();
    fprintf(stdout, "Using %d OpenMP threads...\n", NThreads);
    fflush(stdout);
    #pragma omp parallel num_threads(NThreads) \
            reduction(+:n_isolated) reduction(+:n_capped)
    {
        float *heap_r2 = (float *)malloc(num_ngb * sizeof(float));
        #pragma omp for schedule(dynamic, 512) nowait
        for (int i = 0; i < num_part; i++) {
            float q[KD_NDIM];
            q[0] = posx[i] - data_xmin;
            q[1] = posy[i] - data_ymin;
            q[2] = posz[i] - data_zmin;
            float r2_kth = kd_knn(tree, q, num_ngb, heap_r2);
            if (r2_kth < 0.0f) {
                smoothing_length[i] = 0.0f; n_isolated++;
            } else {
                float h = (float)eta * sqrtf(r2_kth);
                if (h > h_max) { h = h_max; n_capped++; }
                smoothing_length[i] = h;
            }
        }
        free(heap_r2);
    }
#else
    {
        float *heap_r2 = (float *)malloc(num_ngb * sizeof(float));
        for (int i = 0; i < num_part; i++) {
            float q[KD_NDIM];
            q[0] = posx[i] - data_xmin;
            q[1] = posy[i] - data_ymin;
            q[2] = posz[i] - data_zmin;
            float r2_kth = kd_knn(tree, q, num_ngb, heap_r2);
            if (r2_kth < 0.0f) {
                smoothing_length[i] = 0.0f; n_isolated++;
            } else {
                float h = (float)eta * sqrtf(r2_kth);
                if (h > h_max) { h = h_max; n_capped++; }
                smoothing_length[i] = h;
            }
        }
        free(heap_r2);
    }
#endif /* ENABLE_OPENMP */
#endif /* SPH_EXACT_KNN */

#ifdef ENABLE_MPI
    t1 = MPI_Wtime();
    fprintf(stdout, "Neighbour search done in %.2f s (wall)\n", t1 - t0);
#else
    c1 = clock();
    fprintf(stdout, "Neighbour search done in %.2f s\n",
            (double)(c1-c0)/CLOCKS_PER_SEC);
#endif
    fprintf(stdout,
        "find_neighbours: %d isolated + %d capped "
        "(h_px_max=%.4g, MAX_H_VOXELS=%d) out of %d particles\n",
        n_isolated, n_capped, h_max, MAX_H_VOXELS, num_part);
    fflush(stdout);

    kd_free(tree);
    free(pos);

    /* ---- Cache save ---- */
    if (cache_file != NULL) {
        FILE *cf = fopen(cache_file, "wb");
        if (cf) {
            fwrite(&num_part, sizeof(int),   1, cf);
            fwrite(&num_ngb,  sizeof(int),   1, cf);
            float feta = (float)eta;
            fwrite(&feta,     sizeof(float), 1, cf);
            fwrite(smoothing_length, sizeof(float), num_part, cf);
            fclose(cf);
            fprintf(stdout,
                "find_neighbours: saved smoothing lengths to cache '%s'\n",
                cache_file);
            fflush(stdout);
        }
    }
    (void)ymin_arg; (void)zmin_arg;
}

void find_neighbours(int num_part, float *smoothing_length, int num_ngb,
                     float *posx, float *posy, float *posz,
                     float xmin_arg, float ymin_arg, float zmin_arg,
                     float eta)
{
    find_neighbours_cached(num_part, smoothing_length, num_ngb,
                           posx, posy, posz,
                           xmin_arg, ymin_arg, zmin_arg,
                           eta, NULL);
}

/*
 * find_neighbours_fast: O(N) density-grid smoothing length estimator.
 *
 * Use this instead of find_neighbours_cached when:
 *   - N > ~2M particles (KNN becomes memory-bandwidth bound)
 *   - Speed matters more than per-particle accuracy
 *   - Rendering is the goal (not scientific density estimation)
 *
 * 50-100x faster than KNN for large N. Quality is equivalent for
 * rendering purposes: differences only in voids where h=h_max anyway.
 */
void find_neighbours_fast(int num_part, float *smoothing_length, int num_ngb,
                           float *posx, float *posy, float *posz,
                           float xmin_arg, float ymin_arg, float zmin_arg,
                           float eta, const char *cache_file, int render_width)
{
    /* Try cache first */
    if (cache_file != NULL) {
        FILE *cf = fopen(cache_file, "rb");
        if (cf) {
            int cn, cngb; float ceta;
            if (fread(&cn,   sizeof(int),   1, cf) == 1 &&
                fread(&cngb, sizeof(int),   1, cf) == 1 &&
                fread(&ceta, sizeof(float), 1, cf) == 1 &&
                cn == num_part && cngb == num_ngb &&
                fabsf(ceta - (float)eta) < 1e-5f) {
                size_t nr = fread(smoothing_length, sizeof(float), num_part, cf);
                fclose(cf);
                if ((int)nr == num_part) {
                    fprintf(stdout,
                        "find_neighbours: loaded %d smoothing lengths from cache '%s'\n",
                        num_part, cache_file);
                    fflush(stdout);
                    return;
                }
            } else { fclose(cf); }
        }
    }

    fprintf(stdout,
        "Computing smoothing lengths (fast grid method) for %d particles...\n",
        num_part);
    fflush(stdout);

#ifdef ENABLE_MPI
    double t0 = MPI_Wtime();
#else
    clock_t t0 = clock();
#endif

    fast_smoothing_lengths_internal(num_part, smoothing_length, num_ngb,
                                     posx, posy, posz, xmin_arg, eta, render_width);

#ifdef ENABLE_MPI
    fprintf(stdout, "find_neighbours (fast): %.2f s (wall)\n",
            MPI_Wtime() - t0);
#else
    fprintf(stdout, "find_neighbours (fast): %.2f s\n",
            (double)(clock()-t0)/CLOCKS_PER_SEC);
#endif
    fflush(stdout);

    /* Save cache */
    if (cache_file != NULL) {
        FILE *cf = fopen(cache_file, "wb");
        if (cf) {
            fwrite(&num_part, sizeof(int),   1, cf);
            fwrite(&num_ngb,  sizeof(int),   1, cf);
            float feta = (float)eta;
            fwrite(&feta,     sizeof(float), 1, cf);
            fwrite(smoothing_length, sizeof(float), num_part, cf);
            fclose(cf);
            fprintf(stdout, "find_neighbours: saved cache '%s'\n", cache_file);
            fflush(stdout);
        }
    }
    (void)ymin_arg; (void)zmin_arg;
}
