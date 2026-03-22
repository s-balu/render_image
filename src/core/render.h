#ifndef RENDER_IMAGE_RENDER_H
#define RENDER_IMAGE_RENDER_H

/*
 * render.h — Rendering pipeline function prototypes.
 *
 * Covers particle selection, density deposit, SPH smoothing lengths,
 * image output, and the direct 2D SPH deposit.
 *
 * Include this in files that participate in the render pipeline:
 *   render_image.c, smooth_to_mesh.c, write_to_ppm.c, deposit_sph_2d.c
 */

#include "types.h"      /* deposit_particle_t, bitmap_t, pixel_t */
#include "colormap.h"   /* render_config_t */

/* ------------------------------------------------------------------ */
/* Particle selection (select_particles.c)                              */
/* ------------------------------------------------------------------ */

/*
 * Keep particles inside the view volume whose type bit is set in ptype_mask.
 * Compacts survivors to the front of x/y/z/ptype in-place.
 * Sets *NumPartKeep to the number of particles that passed.
 * Pass ptype_mask = -1 to keep all types.
 */
void select_particles(float *x, float *y, float *z, int *ptype,
                      double BoxSize, long long NumPart,
                      float xc, float yc, float zc, float lbox,
                      int ptype_mask, long long *NumPartKeep);

/* ------------------------------------------------------------------ */
/* MPI slab decomposition (split_across_tasks.c)                        */
/* ------------------------------------------------------------------ */

void split_across_tasks_as_slabs(float *x, float *y, float *z,
                                  long long *NumPart,
                                  float xc, float yc, float zc,
                                  float lbox, float BoxSize);

/* ------------------------------------------------------------------ */
/* Smoothing lengths (find_neighbours.c)                                */
/* ------------------------------------------------------------------ */

/*
 * Exact KNN smoothing lengths using a flat kd-tree.
 * Results are optionally cached to disk; reloaded if num_part / num_ngb /
 * eta match the stored values.
 * xmin_arg: voxel size = lbox / IMAGE_DIMENSIONX (used for h pixel cap).
 */
void find_neighbours_cached(int num_part, float *smoothing_length,
                             int num_ngb,
                             float *posx, float *posy, float *posz,
                             float xmin_arg, float ymin_arg, float zmin_arg,
                             float sph_eta, const char *cache_file);

/*
 * O(N) density-grid smoothing length estimator.
 * ~0.9s for 14M particles vs ~160s for exact KNN.
 * render_width: IMAGE_DIMENSIONX — used to keep the physical h cap
 * constant regardless of output resolution.
 */
void find_neighbours_fast(int num_part, float *smoothing_length,
                           int num_ngb,
                           float *posx, float *posy, float *posz,
                           float xmin_arg, float ymin_arg, float zmin_arg,
                           float sph_eta, const char *cache_file,
                           int render_width);

/* ------------------------------------------------------------------ */
/* Density deposit and projection (smooth_to_mesh.c)                    */
/* ------------------------------------------------------------------ */

/*
 * Build the 2D density image from particle positions.
 *
 * With -DKERNEL_SMOOTHING: direct 2D SPH projection using the analytically
 * integrated cubic spline kernel.  smoothing_length must be non-NULL.
 *
 * Without -DKERNEL_SMOOTHING (CIC mode): cloud-in-cell deposit to a 3D
 * grid of depth (int)theta, then summed along z.  smoothing_length ignored.
 *
 * Output is written to data[width * height] in row-major order.
 */
void smooth_to_mesh(long long NumPart, float *smoothing_length,
                    float *x, float *y, float *z,
                    float xc, float yc, float zc, float lbox,
                    float theta,
                    int width, int height, float *data);

/* ------------------------------------------------------------------ */
/* Direct 2D SPH deposit (deposit_sph_2d.c)                             */
/* ------------------------------------------------------------------ */

/* Initialise the local projected kernel lookup table.
 * Must be called once before deposit_sph_2d(). */
void deposit_sph_2d_init(void);

/*
 * Deposit n_packed pre-packed particles directly onto the 2D image.
 * particles[i].h is the kernel radius in pixels.
 * Writes to data2d[height * width], zeroing it first.
 * Returns total number of pixel writes (diagnostic).
 */
long long deposit_sph_2d(long long n_packed,
                          const deposit_particle_t *particles,
                          int full_width, int height,
                          float *data2d);

/* ------------------------------------------------------------------ */
/* Image output (write_to_ppm.c)                                        */
/* ------------------------------------------------------------------ */

/* Simple greyscale PNG (no colormap). */
void write_to_png(char *image_file, int width, int height, float *data);

/* Full pipeline: colormap + opacity + background compositing. */
void write_to_png_ex(const char *image_file, int width, int height,
                     const float *data, const render_config_t *cfg);

/* PPM output (debug / fallback). */
void write_to_ppm(char *image_file, int width, int height,
                  int MaxColorComponentValue, float *data);

/* Low-level PNG helpers (used internally by write_to_ppm.c). */
pixel_t *pixel_at(bitmap_t *bitmap, int x, int y);
int      save_png_to_file(bitmap_t *bitmap, const char *path);

#endif /* RENDER_IMAGE_RENDER_H */
