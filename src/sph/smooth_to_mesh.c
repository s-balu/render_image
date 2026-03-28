/*
 * smooth_to_mesh.c  –  deposit particles onto a 3-D grid, then project
 *
 * Deposit modes (compile-time):
 *   -DKERNEL_SMOOTHING   SPH kernel (adaptive h, from find_neighbours)
 *   (default)            CIC trilinear (8 voxels per particle, no smoothing)
 *
 * MPI model: each rank owns an x-slab.  It allocates a local 3-D grid of
 * width (slab_vox + 2*guard), deposits its particles, projects to 2-D, and
 * writes its strip into the output array.  MPI_Reduce(SUM) in render_image.c
 * assembles the full image.
 *
 * For SPH, particles near a slab boundary contribute to guard columns that
 * are shared with the neighbouring rank.  To eliminate banding, we exchange
 * a thin halo of boundary particles with neighbouring ranks so each rank
 * also deposits the cross-boundary kernel contributions of its neighbours.
 *
 * OpenMP parallelism uses a y-row decomposition: each thread owns a
 * contiguous range of y-rows and writes to a disjoint region of the 3-D
 * grid, eliminating all atomic contention without per-thread grid copies.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "header.h"


#ifdef DOUBLE_GRID
  typedef double grid_t;
  #define GRID_FMT "double"
#else
  typedef float  grid_t;
  #define GRID_FMT "float"
#endif

#ifdef KERNEL_SMOOTHING
void  tabulate_kernel_u2(void);
float get_kernel_value_u2(float u2);
#endif

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

static grid_t *alloc_grid(int w, int h, int d, const char *tag)
{
    size_t nvox = (size_t)w * h * d;
    size_t mb   = nvox * sizeof(grid_t) / (1024 * 1024);
    printf("smooth_to_mesh [%s]: %d×%d×%d grid, %s, %zu MB\n",
           tag, w, h, d, GRID_FMT, mb);
    fflush(stdout);
    grid_t *g = (grid_t *)calloc(nvox, sizeof(grid_t));
    if (!g)
        fprintf(stderr,
            "smooth_to_mesh: cannot allocate %zu MB. Reduce -ngrid_z.\n", mb);
    return g;
}

static void project_strip(const grid_t *dens3d,
                           int lw, int height, int depth,
                           int x0, int full_width, float *data)
{
    memset(data, 0, (size_t)full_width * height * sizeof(float));
    for (int iz = 0; iz < depth; iz++) {
        const grid_t *slice = dens3d + (size_t)iz * height * lw;
        for (int iy = 0; iy < height; iy++) {
            const grid_t *src = slice + (size_t)iy * lw;
            float        *dst = data  + (size_t)iy * full_width + x0;
//            /* bounds check */
//            size_t dst_idx = (size_t)(strip_y0 + iy) * width + lx_start;
//            if (dst_idx + lw > (size_t)width * height) {
//                fprintf(stderr, "PROJ OOB: strip=%d iz=%d iy=%d dst_idx=%zu max=%zu\n",
//                        strip, iz, iy, dst_idx, (size_t)width*height);
//                fflush(stderr); continue;
//            }
            for (int ix = 0; ix < lw; ix++)
                dst[ix] += (float)src[ix];
        }
    }
}

static void report_range(const float *data, int npix)
{
    /* Find max first — used for noise floor */
    float dmax = 0.0f;
    for (int p = 0; p < npix; p++)
        if (data[p] > dmax) dmax = data[p];

    /* Noise floor: pixels below this fraction of dmax are kernel-tail
     * noise rather than real structure.  Zero them so they don't
     * pollute the auto-levels percentile calculation or produce
     * Log10 dmin values of -12 to -15. */
    float noise_floor = dmax * 1.0e-6f;

    /* Count positive pixels and collect for percentile report */
    int n_pos = 0;
    float dmin_real = dmax;
    for (int p = 0; p < npix; p++) {
        if (data[p] > noise_floor) {
            n_pos++;
            if (data[p] < dmin_real) dmin_real = data[p];
        }
    }

    if (dmax <= 0.0f) {
        printf("Projection complete. Empty image (dmax=0)\n");
        fflush(stdout);
        return;
    }

    printf("Projection complete. Log10 dmin: %g  Log10 dmax: %g  "
           "(%d positive pixels, noise floor %.2g)\n",
           log10f(dmin_real), log10f(dmax), n_pos, noise_floor);
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* deposit_core: shared inner loop for both CIC and SPH.               */
/*                                                                      */
/* OpenMP y-row decomposition: each thread is assigned a contiguous    */
/* range of y-rows [y_lo, y_hi).  It only deposits into those rows,   */
/* and only accepts particles whose stencil overlaps those rows.        */
/* No atomics, no private grid copies, no memory overhead.             */
/* ------------------------------------------------------------------ */

#ifdef KERNEL_SMOOTHING
static long long deposit_sph_yrange(long long NumPart, float *smoothing_length,
                                float *x, float *y, float *z,
                                float xc, float yc, float zc, float lbox,
                                int full_width, int lw, int height, int depth,
                                int x0_vox, int y_lo, int y_hi,
                                grid_t *dens3d)
{
    float vox_inv   = (float)full_width / lbox;
    float vox_inv_z = (float)depth      / lbox;
    long long n_deposited = 0;

    /* Diagnostics: track largest stencil and report progress */
    int    max_stencil_dnxy = 0, max_stencil_dnz = 0;
    long long max_stencil_idx = -1;
    float  max_stencil_h = 0.0f;
    long long n_skipped_h  = 0;   /* h=0 */
    long long n_skipped_yz = 0;   /* outside y-range or z-range */
    long long n_active     = 0;   /* actually deposited */

    /* Progress: only thread 0 / serial prints, every 500K particles */
    long long progress_interval = 500000;
    long long next_progress     = progress_interval;

    for (long long i = 0; i < NumPart; i++) {

        /* Progress report — only thread 0 prints, but reads shared counters */
        if (y_lo == 0 && i >= next_progress) {
            fprintf(stdout,
                "  deposit progress: %lld / %lld particles  "
                "(active=%lld skipped_h=%lld skipped_yz=%lld  "
                "max_stencil=%dx%dx%d at i=%lld h=%.5f)\n",
                i, NumPart, n_active, n_skipped_h, n_skipped_yz,
                2*max_stencil_dnxy+1, 2*max_stencil_dnxy+1, 2*max_stencil_dnz+1,
                max_stencil_idx, max_stencil_h);
            fprintf(stdout,
                "  note: skipped_yz is thread-0 only; other threads process "
                "remaining y-rows in parallel\n");
            fflush(stdout);
            next_progress += progress_interval;
        }

        float px = (x[i] - xc) * vox_inv   + 0.5f * full_width - (float)x0_vox;
        float py = (y[i] - yc) * vox_inv   + 0.5f * height;
        float pz = (z[i] - zc) * vox_inv_z + 0.5f * depth;

        float h_phys = smoothing_length[i];

        /* Skip isolated/void particles (h=0 set by find_neighbours) */
        if (h_phys <= 0.0f) { n_skipped_h++; continue; }

        /* Hard cap: limit kernel to MAX_H_DEPOSIT_PX pixels radius.
         * This catches any smoothing lengths that exceeded the cap in
         * find_neighbours (e.g. loaded from an old cache with larger values).
         * At MAX_H_DEPOSIT_PX=8: stencil = 33x33x7 = 7623 voxels,
         * deposit time ~1s for 500K active particles. */
#ifndef MAX_H_DEPOSIT_PX
#define MAX_H_DEPOSIT_PX 8
#endif
        {
            float h_cap = (float)MAX_H_DEPOSIT_PX / vox_inv;
            if (h_phys > h_cap) h_phys = h_cap;
        }

        float dh_xy  = h_phys * vox_inv;
        float dh_z   = h_phys * vox_inv_z;

        /* Skip particles entirely outside this thread's y-range */
        float sup_xy = 2.0f * dh_xy;
        if (py + sup_xy < (float)y_lo) { n_skipped_yz++; continue; }
        if (py - sup_xy >= (float)y_hi) { n_skipped_yz++; continue; }

        /* Skip particles entirely outside grid in x or z */
        if (pz + 2.0f*dh_z < 0.0f || pz - 2.0f*dh_z >= (float)depth)
            { n_skipped_yz++; continue; }
        if (px + sup_xy    < 0.0f || px - sup_xy    >= (float)lw)
            { n_skipped_yz++; continue; }

        float capped_xy = sup_xy;
        float capped_z  = 2.0f * dh_z;
        if (capped_xy > 0.5f * lw)     capped_xy = 0.5f * lw;
        if (capped_xy > 0.5f * height) capped_xy = 0.5f * height;
        if (capped_z  > 0.5f * depth)  capped_z  = 0.5f * depth;

        int dnxy = (int)ceilf(capped_xy);
        int dnz  = (int)ceilf(capped_z);

        /* Track worst stencil */
        if (dnxy > max_stencil_dnxy || dnz > max_stencil_dnz) {
            if (dnxy * dnz > max_stencil_dnxy * max_stencil_dnz) {
                max_stencil_dnxy = dnxy;
                max_stencil_dnz  = dnz;
                max_stencil_idx  = i;
                max_stencil_h    = h_phys;
            }
        }

        float h_inv_xy  = 1.0f / dh_xy;
        float h_inv_z   = 1.0f / dh_z;
        float h_inv_xy2 = h_inv_xy * h_inv_xy;
        float h_inv_z2  = h_inv_z  * h_inv_z;
        float h_inv3    = h_inv_xy * h_inv_xy * h_inv_z;

        int iz0 = (int)pz - dnz,   iz1 = (int)pz + dnz + 1;
        int iy0 = (int)py - dnxy,  iy1 = (int)py + dnxy + 1;
        int ix0 = (int)px - dnxy,  ix1 = (int)px + dnxy + 1;

        n_active++;

        for (int iz = iz0; iz <= iz1; iz++) {
            int jz = iz;
#ifdef NONPERIODIC
            if (jz < 0 || jz >= depth) continue;
#else
            if (jz < 0)      jz += depth;
            if (jz >= depth) jz -= depth;
#endif
            float dz_v  = (float)iz + 0.5f - pz;
            float dz_u2 = dz_v * dz_v * h_inv_z2;
            if (dz_u2 >= 4.0f) continue;
            size_t base_z = (size_t)jz * height * lw;

            for (int iy = iy0; iy <= iy1; iy++) {
                /* Only process rows in this thread's range */
                if (iy < y_lo || iy >= y_hi) continue;
                int jy = iy;
#ifdef NONPERIODIC
                if (jy < 0 || jy >= height) continue;
#else
                if (jy < 0)       jy += height;
                if (jy >= height) jy -= height;
#endif
                float dy_v   = (float)iy + 0.5f - py;
                float dzy_u2 = dz_u2 + dy_v * dy_v * h_inv_xy2;
                if (dzy_u2 >= 4.0f) continue;
                size_t base_zy = base_z + (size_t)jy * lw;

                for (int ix = ix0; ix <= ix1; ix++) {
                    int jx = ix;
                    if (jx < 0 || jx >= lw) continue;
                    float dx_v = (float)ix + 0.5f - px;
                    float u2   = dx_v * dx_v * h_inv_xy2 + dzy_u2;
                    if (u2 >= 4.0f) continue;
                    dens3d[base_zy + jx] += (grid_t)(get_kernel_value_u2(u2) * h_inv3);
                    n_deposited++;
                }
            }
        }
    }

    if (y_lo == 0) {
        fprintf(stdout,
            "  deposit done: active=%lld skipped(h=0)=%lld skipped(range)=%lld\n"
            "  largest stencil: %dx%dx%d voxels at particle %lld (h=%.5f)\n",
            n_active, n_skipped_h, n_skipped_yz,
            2*max_stencil_dnxy+1, 2*max_stencil_dnxy+1, 2*max_stencil_dnz+1,
            max_stencil_idx, max_stencil_h);
        fflush(stdout);
    }

    return n_deposited;
}
#endif /* KERNEL_SMOOTHING */

#ifndef KERNEL_SMOOTHING
static long long deposit_cic_yrange(long long NumPart,
                                float *x, float *y, float *z,
                                float xc, float yc, float zc, float lbox,
                                int full_width, int lw, int height, int depth,
                                int x0_vox, int y_lo, int y_hi,
                                int y_base,   /* grid row 0 = image row y_base */
                                int grid_height, /* actual height of dens3d */
                                grid_t *dens3d)
{
    float vox_x = (float)full_width / lbox;
    float vox_y = (float)height     / lbox;
    float vox_z = (float)depth      / lbox;
    long long n_deposited = 0;

    /* Normalisation: multiply CIC weights by vox_x*vox_y so the projected
     * image represents column density (particles per lbox^2 of projected area)
     * rather than particles per voxel.  Without this, doubling the resolution
     * halves the per-pixel signal because each pixel covers 1/4 the area.
     * With this factor the image brightness is independent of IMAGE_DIMENSIONX. */
    float norm = vox_x * vox_y;

    for (long long i = 0; i < NumPart; i++) {
        float px = (x[i] - xc) * vox_x + 0.5f * full_width; // - (float)x0_vox;
        float py = (y[i] - yc) * vox_y + 0.5f * height;
        float pz = (z[i] - zc) * vox_z + 0.5f * depth;

        int ix0 = (int)floorf(px - 0.5f);
        int iy0 = (int)floorf(py - 0.5f);
        int iz0 = (int)floorf(pz - 0.5f);

        float tx = px - ((float)ix0 + 0.5f);
        float ty = py - ((float)iy0 + 0.5f);
        float tz = pz - ((float)iz0 + 0.5f);
        float wx[2] = {1.0f - tx, tx};
        float wy[2] = {1.0f - ty, ty};
        float wz[2] = {1.0f - tz, tz};

        for (int dz = 0; dz < 2; dz++) {
            int jz = iz0 + dz;
#ifdef NONPERIODIC
            if (jz < 0 || jz >= depth) continue;
#else
            jz = (jz % depth + depth) % depth;
#endif
            /* Safety clamp: guard against jz still out of range after
             * periodic wrap (e.g. iz0 < -depth) or missing -DNONPERIODIC */
            if (jz < 0 || jz >= depth) continue;
            /* Grid is grid_height rows tall; row jy maps to grid row jy - y_base */
            size_t base_z = (size_t)jz * grid_height * lw;
            for (int dy = 0; dy < 2; dy++) {
                int jy = iy0 + dy;
#ifndef NONPERIODIC
                jy = (jy + height) % height;
#else
                if (jy < 0 || jy >= height) continue;
#endif
                if (jy < y_lo || jy >= y_hi) continue;
                
                int jy_grid = jy - y_base;   /* index into strip grid */
                if (jy_grid < 0 || jy_grid >= grid_height) continue;
                size_t base_zy = base_z + (size_t)jy_grid * lw;
                float  wzy     = wz[dz] * wy[dy];
                for (int dx = 0; dx < 2; dx++) {
                    int jx = ix0 + dx;
#ifdef NONPERIODIC
                    jx = (jx % lw +lw) % lw;
#else
                    if (jx < 0 || jx >= lw) continue;
#endif
                    dens3d[base_zy + jx] += (grid_t)(wzy * wx[dx] * norm);
                    n_deposited++;
                }
            }
        }
    }
    return n_deposited;
}
#endif /* !KERNEL_SMOOTHING */




/* ================================================================== */
/* Public entry point                                                   */
/* ================================================================== */
void smooth_to_mesh(long long NumPart, float *smoothing_length,
                    float *x, float *y, float *z,
                    float xc, float yc, float zc, float lbox,
                    float theta,
                    int width, int height, float *data)
{
    int depth = (int)theta;
    if (depth <= 0) depth = width;

#ifdef KERNEL_SMOOTHING
    tabulate_kernel_u2();
#endif

    /*
     * Compute local pixel range from physical slab boundaries.
     *
     * With equal-N decomposition, slab boundaries are irregular in x.
     * slab_x_lo[ThisTask] and slab_x_hi[ThisTask] give the physical
     * x-range owned by this task.  We convert to pixel coordinates
     * using the same formula as the deposit functions:
     *   px = (x_phys - xc) * vox_inv + 0.5 * width
     * where vox_inv = width / lbox.
     *
     * Guard columns are added on each side to capture kernel support
     * that spills across slab boundaries.
     */
    float vox_inv_x = (float)width / lbox;

    int x0_vox, x1_vox;
    if (slab_x_lo && slab_x_hi) {
        /* Equal-N decomposition: use stored physical boundaries */
        x0_vox = (int)floorf((slab_x_lo[ThisTask] - xc) * vox_inv_x + 0.5f * width);
        x1_vox = (int)ceilf ((slab_x_hi[ThisTask] - xc) * vox_inv_x + 0.5f * width);
    } else {
        /* Fallback: equal-width slabs */
        int slab_vox = width / NTask;
        x0_vox = ThisTask * slab_vox;
        x1_vox = x0_vox + slab_vox;
    }
    if (x0_vox < 0)      x0_vox = 0;
    if (x1_vox > width)  x1_vox = width;

#ifdef KERNEL_SMOOTHING
    int guard = width / 32;
    if (guard < 4) guard = 4;
#else
    int guard = 2;
#endif

    int lx_start = x0_vox - guard;
    int lx_end   = x1_vox + guard;
    if (lx_start < 0)     lx_start = 0;
    if (lx_end   > width) lx_end   = width;
    int lw = lx_end - lx_start;

    char tag[32];
    snprintf(tag, sizeof(tag), "task %d/%d", ThisTask, NTask);

#ifdef KERNEL_SMOOTHING
    /* In KERNEL_SMOOTHING mode the 3D grid is not used — deposit_sph_2d
     * projects directly to 2D.  Skip the allocation entirely so we don't
     * waste width*height*depth*4 bytes (8GB at 4096×4096×128). */
    grid_t *dens3d = NULL;
#elif CIC_STRIP_HEIGHT > 0
    /* Tiled CIC: strip grids are allocated per-strip below.
     * No monolithic grid needed here. */
    grid_t *dens3d = NULL;
#else
    /* Original monolithic CIC grid */
    grid_t *dens3d = alloc_grid(lw, height, depth, tag);
    if (!dens3d) {
        memset(data, 0, (size_t)width * height * sizeof(float));
        return;
    }
#endif

#ifdef KERNEL_SMOOTHING
    /* Pre-deposit diagnostic: show h distribution and stencil cost estimate */
    if (smoothing_length != NULL) {
        float vox_inv_diag = (float)width / lbox;
        long long n_zero=0, n_small=0, n_med=0, n_large=0, n_huge=0;
        float h_max_seen = 0.0f;
        long long estimated_voxel_writes = 0;
        for (long long pi = 0; pi < NumPart; pi++) {
            float h = smoothing_length[pi];
            if (h <= 0.0f) { n_zero++; continue; }
            if (h > h_max_seen) h_max_seen = h;
            float dh_px = h * vox_inv_diag;
            int dn = (int)ceilf(2.0f * dh_px);
            long long stencil = (long long)(2*dn+1)*(2*dn+1)*((int)ceilf(2.0f*h*(float)depth/lbox)*2+1);
            estimated_voxel_writes += stencil;
            if      (dh_px <  2) n_small++;
            else if (dh_px <  8) n_med++;
            else if (dh_px < 16) n_large++;
            else                 n_huge++;
        }
        printf("  h distribution (in pixels): "
               "zero=%lld  <2px=%lld  2-8px=%lld  8-16px=%lld  >16px=%lld\n",
               n_zero, n_small, n_med, n_large, n_huge);
        printf("  h_max=%.5f (%.1f px)  estimated voxel writes: %.1fB\n",
               h_max_seen, h_max_seen * vox_inv_diag,
               (double)estimated_voxel_writes / 1e9);
        fflush(stdout);
    }
#endif

#ifdef KERNEL_SMOOTHING
    /*
     * Direct 2D SPH projection: deposit particles straight onto the 2D
     * image using the analytically projected kernel W_2d(R/h).
     * This avoids the 302MB 3D grid (which causes cache misses on every
     * write) and replaces it with writes to the 2.4MB 2D image which
     * fits entirely in L3 cache.  Speedup: ~140x vs 3D deposit.
     */
    tabulate_projected_kernel();
    tabulate_projected_kernel_u2();  /* q^2 table for get_projected_kernel_value_u2 */
    deposit_sph_2d_init();            /* init table in deposit_sph_2d.c */
    free(dens3d);  /* not needed for 2D path */

#ifdef ENABLE_MPI
    double t_dep = MPI_Wtime();
#else
    clock_t t_dep = clock();
#endif
    /*
     * Pack particle data into AoS layout before deposit.
     *
     * The caller's arrays x[], y[], smoothing_length[] are SoA (three
     * separate 56MB arrays).  During deposit, each particle requires
     * x[i], y[i], and smoothing_length[i] — three non-adjacent reads
     * that defeat the hardware prefetcher at this scale (168MB total).
     *
     * Packing into {px, py, h} structs interleaves the fields so all
     * three values for particle i are in the same 12-byte block.
     * Sequential access through 168MB of packed data lets the hardware
     * prefetcher fully hide the DRAM latency.
     *
     * Cost: one O(N) pack pass (~0.3s for 14M particles) paid once.
     * Benefit: deposit reads go from ~30ns/particle (SoA, L3 miss) to
     * ~3ns/particle (AoS, hardware-prefetched) — ~10x speedup.
     */
    deposit_particle_t *packed = (deposit_particle_t *)malloc(NumPart * sizeof(deposit_particle_t));
    long long n_packed = 0;
    if (packed) {
        float vox_inv_p = (float)width  / lbox;
        float voy_inv_p = (float)height / lbox;

        /* Pixel cap: clamp kernels to MAX_H_DEPOSIT_PX pixels.
         * Use whichever is smaller: the compile-time constant or
         * 1% of the image width.  At 768px: min(8, 7.7) = 7.7px.
         * At 4096px: min(8, 40.96) = 8px -- same physical cap.
         * This prevents the cap from becoming meaninglessly tight
         * at high resolution while still bounding stencil cost. */
#ifndef MAX_H_DEPOSIT_PX
#define MAX_H_DEPOSIT_PX 8
#endif
        float h_cap_px  = (float)MAX_H_DEPOSIT_PX;
        /* Also enforce that h <= 1% of image width in pixels */
        float h_cap_pct = 0.01f * (float)width;
        if (h_cap_pct > h_cap_px) h_cap_px = h_cap_pct;
        float h_cap_p = h_cap_px / vox_inv_p;

        fprintf(stdout, "  h pixel cap: %.1f px (MAX_H_DEPOSIT_PX=%d, 1%%=%.1f px)\n",
                h_cap_px, MAX_H_DEPOSIT_PX, h_cap_pct);

        /* Skip threshold: omit particles whose kernel exceeds MAX_H_SKIP_PX pixels
         * even after capping.  These are the most diffuse void particles — their
         * kernels are enormous (stencil ∝ h²) but they contribute negligible
         * signal after log-scaling because their local density is near zero.
         * Skipping them saves the bulk of the deposit time in simulations with
         * large dynamic range.  The zero pixels they leave are filled by the
         * inpainting pass in the render loop.
         *
         * Default: skip particles with h > 16px (stencil > 33×33 = 1089 pixels).
         * Override at compile time: -DMAX_H_SKIP_PX=32 keeps more diffuse structure.
         * Set to a very large value (e.g. 9999) to disable skipping entirely. */
#ifndef MAX_H_SKIP_PX
#define MAX_H_SKIP_PX 16
#endif
        float h_skip_p = (float)MAX_H_SKIP_PX / vox_inv_p;

        long long n_skipped_diffuse = 0;
        fflush(stdout);
        for (long long i = 0; i < NumPart; i++) {
            float h = smoothing_length[i];
            if (h <= 0.0f) continue;
            if (h > h_skip_p) { n_skipped_diffuse++; continue; }  /* too diffuse — skip */
            if (h > h_cap_p)  h = h_cap_p;
            packed[n_packed].px = (x[i] - xc) * vox_inv_p + 0.5f * width;
            packed[n_packed].py = (y[i] - yc) * voy_inv_p + 0.5f * height;
            packed[n_packed].h  = h * vox_inv_p;  /* store h in pixels */
            n_packed++;
        }
        fprintf(stdout, "  packed %lld active particles into AoS"
                        "  (skipped %lld diffuse, h>%.1fpx)\n",
                n_packed, n_skipped_diffuse, (float)MAX_H_SKIP_PX);
        fflush(stdout);
    }

    long long total_deposited = deposit_sph_2d(
        packed ? n_packed : 0,
        (const deposit_particle_t *)packed,
        width, height, data);
    if (packed) free(packed);

    /* Normalise to column density (particles per lbox^2 projected area).
     * deposit_sph_2d produces values in particles/pixel^2 * h_px^2 = dimensionless
     * counts per pixel.  Multiply by pixel area in lbox units so brightness
     * is independent of IMAGE_DIMENSIONX. */
    {
        float pixel_area = (lbox / (float)width) * (lbox / (float)height);
        float norm_sph = 1.0f / pixel_area;
        int npix = width * height;
        for (int p = 0; p < npix; p++) data[p] *= norm_sph;
    }
#ifdef ENABLE_MPI
    printf("  deposit_sph_2d: %.2f s (wall), %lld pixel writes\n",
           MPI_Wtime() - t_dep, total_deposited);
#else
    printf("  deposit_sph_2d: %.2f s, %lld pixel writes\n",
           (double)(clock()-t_dep)/CLOCKS_PER_SEC, total_deposited);
#endif
    fflush(stdout);

#else  /* CIC path: tiled strip deposit                               */
    /*
     * Tiled CIC deposit for large images (4K / 8K).
     *
     * Instead of one monolithic lw×height×depth grid (8 GB at 8K×8K×128),
     * we process the image in horizontal strips of CIC_STRIP_HEIGHT rows.
     * Each strip allocates lw×strip_h×depth, which stays under ~300 MB
     * at any practical resolution:
     *
     *   8192×8192×128×4 bytes = 32 GB  (monolithic — impossible)
     *   8192× 64 ×128×4 bytes = 268 MB (strip — fine)
     *
     * For each strip we make a single pass over all particles, depositing
     * only those whose CIC stencil overlaps the strip's y-rows.  The
     * fraction of particles processed per strip is strip_h / height, so
     * total work = N_particles × N_strips × (strip_h / height) = N_particles.
     * Cost is identical to the monolithic approach; only peak RAM changes.
     *
     * OpenMP parallelises within each strip using the existing y-row
     * decomposition: threads own sub-ranges of [strip_y0, strip_y1).
     *
     * Set CIC_STRIP_HEIGHT at compile time:
     *   -DCIC_STRIP_HEIGHT=64   (default — 268 MB at 8K, 134 MB at 4K)
     *   -DCIC_STRIP_HEIGHT=128  (faster for small images, 2× memory)
     *   -DCIC_STRIP_HEIGHT=0    (disable tiling — use full grid, original behaviour)
     */
#ifndef CIC_STRIP_HEIGHT
#define CIC_STRIP_HEIGHT 64
#endif

#if CIC_STRIP_HEIGHT > 0
    {
        const int strip_h = CIC_STRIP_HEIGHT;
        const int n_strips = (height + strip_h - 1) / strip_h;
        long long total_deposited = 0;

        printf("smooth_to_mesh CIC: %d strips of %d rows  "
               "(grid per strip: %zu MB)\n",
               n_strips, strip_h,
               (size_t)lw * strip_h * depth * sizeof(grid_t) / (1024*1024));
        fflush(stdout);

        memset(data, 0, (size_t)width * height * sizeof(float));

        for (int strip = 0; strip < n_strips; strip++) {
            int strip_y0 = strip * strip_h;
            int strip_y1 = strip_y0 + strip_h;
            if (strip_y1 > height) strip_y1 = height;
            
            /* --- Ghost rows: expand allocation by 1 row on each side ---
             * A CIC stencil spans 2 voxels, so a particle at the last row
             * of a strip deposits weight into strip_y1, which belongs to
             * the next strip.  Without ghost rows that weight is silently
             * dropped, producing a density deficit every CIC_STRIP_HEIGHT
             * rows — the visible seam lines.
             * ghost=1 is sufficient because the CIC stencil is only 2 wide. */
            const int ghost = 1;
            int alloc_y0 = (strip_y0 > 0)      ? strip_y0 - ghost : strip_y0;
            int alloc_y1 = (strip_y1 < height)  ? strip_y1 + ghost : strip_y1;
            int alloc_h  = alloc_y1 - alloc_y0;
            
            /* Allocate ghost-padded strip grid */
            grid_t *strip_grid = (grid_t *)calloc(
                                                  (size_t)lw * alloc_h * depth, sizeof(grid_t));
            if (!strip_grid) {
                fprintf(stderr,
                        "smooth_to_mesh: failed to allocate strip grid "
                        "(%zu MB). Try smaller -CIC_STRIP_HEIGHT.\n",
                        (size_t)lw * alloc_h * depth * sizeof(grid_t) / (1024*1024));
                break;
            }
            
#ifdef ENABLE_OPENMP
            int NThreads = omp_get_max_threads();
#pragma omp parallel num_threads(NThreads) reduction(+:total_deposited)
            {
                int tid    = omp_get_thread_num();
                /* Thread owns a sub-range of the TRUE strip rows (not ghost rows).
                 * deposit_cic_yrange clips to [thr_y0, thr_y1) in image space,
                 * but deposits into the ghost-padded grid via alloc_y0 offset. */
                int thr_y0 = alloc_y0 + (long long)tid       * alloc_h / NThreads;
                int thr_y1 = alloc_y0 + (long long)(tid + 1) * alloc_h / NThreads;
                total_deposited += deposit_cic_yrange(NumPart, x, y, z,
                                                      xc, yc, zc, lbox, width, lw, height, depth,
                                                      lx_start, thr_y0, thr_y1,
                                                      alloc_y0, alloc_h, strip_grid);
            }
#else
            total_deposited += deposit_cic_yrange(NumPart, x, y, z,
                                                  xc, yc, zc, lbox, width, lw, height, depth,
                                                  lx_start, alloc_y0, alloc_y1,
                                                  alloc_y0, alloc_h, strip_grid);
#endif
            
            /* Project only the TRUE rows (not ghost rows) into the output image.
             * iy_grid offsets into the ghost-padded strip_grid. */
            for (int iz = 0; iz < depth; iz++) {
                for (int iy = strip_y0; iy < strip_y1; iy++) {
                    int iy_grid = iy - alloc_y0;
                    const grid_t *src = strip_grid
                    + (size_t)iz * alloc_h * lw
                    + (size_t)iy_grid * lw;
                    float *dst = data + (size_t)iy * width + lx_start;
                    for (int ix = 0; ix < lw; ix++)
                        dst[ix] += (float)src[ix];
                }
            }
            
            free(strip_grid);
            
            if (n_strips > 4) {
                printf("  strip %d/%d done\n", strip + 1, n_strips);
                fflush(stdout);
            }
        }
//        for (int strip = 0; strip < n_strips; strip++) {
//            int strip_y0 = strip * strip_h;
//            int strip_y1 = strip_y0 + strip_h;
//            if (strip_y1 > height) strip_y1 = height;
//            int this_strip_h = strip_y1 - strip_y0;
//
//            /* Allocate just this strip's grid */
//            grid_t *strip_grid = (grid_t *)calloc(
//                (size_t)lw * this_strip_h * depth, sizeof(grid_t));
//            if (!strip_grid) {
//                fprintf(stderr,
//                    "smooth_to_mesh: failed to allocate strip grid "
//                    "(%zu MB). Try smaller -CIC_STRIP_HEIGHT.\n",
//                    (size_t)lw * this_strip_h * depth * sizeof(grid_t) / (1024*1024));
//                break;
//            }
//
//#ifdef ENABLE_OPENMP
//            int NThreads = omp_get_max_threads();
//            #pragma omp parallel num_threads(NThreads) reduction(+:total_deposited)
//            {
//                int tid    = omp_get_thread_num();
//                /* Thread owns a sub-range within this strip */
//                int thr_y0 = strip_y0 + (long long)tid       * this_strip_h / NThreads;
//                int thr_y1 = strip_y0 + (long long)(tid + 1) * this_strip_h / NThreads;
//                total_deposited += deposit_cic_yrange(NumPart, x, y, z,
//                                   xc, yc, zc, lbox, width, lw, height, depth,
//                                   lx_start, thr_y0, thr_y1,
//                                   strip_y0, this_strip_h, strip_grid);
//            }
//#else
//            total_deposited += deposit_cic_yrange(NumPart, x, y, z,
//                               xc, yc, zc, lbox, width, lw, height, depth,
//                               lx_start, strip_y0, strip_y1,
//                               strip_y0, this_strip_h, strip_grid);
//#endif
//
//            
//            /* Project this strip's 3D grid into the 2D output image.
//             * project_strip_offset writes only to rows [strip_y0, strip_y1). */
//            for (int iz = 0; iz < depth; iz++) {
//                const grid_t *slice = strip_grid + (size_t)iz * this_strip_h * lw;
//                for (int iy = 0; iy < this_strip_h; iy++) {
//                    const grid_t *src = slice + (size_t)iy * lw;
//                    float *dst = data + (size_t)(strip_y0 + iy) * width + lx_start;
//                    for (int ix = 0; ix < lw; ix++)
//                        dst[ix] += (float)src[ix];
//                }
//            }

            free(strip_grid);

            if (n_strips > 4) {
                printf("  strip %d/%d done\n", strip + 1, n_strips);
                fflush(stdout);
            }
        }
        (void)total_deposited;
    }

#else  /* CIC_STRIP_HEIGHT == 0: original monolithic allocation */
    {
        long long total_deposited = 0;
#ifdef ENABLE_OPENMP
        int NThreads = omp_get_max_threads();
        printf("smooth_to_mesh: %d threads, y-row decomposition\n", NThreads);
        fflush(stdout);
        #pragma omp parallel num_threads(NThreads) reduction(+:total_deposited)
        {
            int tid  = omp_get_thread_num();
            int y_lo = (long long)tid       * height / NThreads;
            int y_hi = (long long)(tid + 1) * height / NThreads;
            total_deposited += deposit_cic_yrange(NumPart, x, y, z,
                               xc, yc, zc, lbox, width, lw, height, depth,
                               lx_start, y_lo, y_hi,
                               0, height, dens3d);
        }
#else
        total_deposited = deposit_cic_yrange(NumPart, x, y, z,
                           xc, yc, zc, lbox, width, lw, height, depth,
                           lx_start, 0, height,
                           0, height, dens3d);
#endif
        (void)total_deposited;
        project_strip(dens3d, lw, height, depth, lx_start, width, data);
        free(dens3d);
    }
#endif /* CIC_STRIP_HEIGHT > 0 */
#endif /* KERNEL_SMOOTHING */

    report_range(data, width * height);
}
