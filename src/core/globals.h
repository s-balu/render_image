#ifndef RENDER_IMAGE_GLOBALS_H
#define RENDER_IMAGE_GLOBALS_H

/*
 * globals.h — extern declarations for process-global mutable state.
 *
 * These are defined once in header.c and shared across translation units.
 * Include this header only in files that actually read or write these
 * variables.  Files that don't need them (kernels.c, colormap.c,
 * make_tree.c, walk_tree.c, select_particles.c, write_to_ppm.c) should
 * not include it.
 *
 * Future direction: thread these through a render_context_t struct rather
 * than keeping them as globals, which would make the code testable and
 * eventually thread-safe.
 */

/* MPI rank and world size.  Both are 0 and 1 respectively in serial builds. */
extern int ThisTask;
extern int NTask;

/* Per-task x-slab boundaries in physical coordinates.
 * slab_x_lo[t] / slab_x_hi[t] are the physical x-range owned by task t.
 * Set by split_across_tasks_as_slabs(), read by smooth_to_mesh().
 * Allocated as NTask floats by split_across_tasks_as_slabs(). */
extern float *slab_x_lo;
extern float *slab_x_hi;

/* Snapshot format tag set by the I/O layer. */
extern int SnapFormat;

/* Number of particles on this MPI task (set during I/O). */
extern int NThisTask;

#endif /* RENDER_IMAGE_GLOBALS_H */
