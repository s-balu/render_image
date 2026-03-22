#ifndef RENDER_IMAGE_HEADER_H
#define RENDER_IMAGE_HEADER_H

/*
 * header.h — compatibility umbrella header.
 *
 * All declarations have been split into focused headers:
 *
 *   types.h    — POD structs and constants (sim_info, tree nodes, etc.)
 *   globals.h  — extern globals (ThisTask, NTask, slab_x_lo, etc.)
 *   tree.h     — octree build and walk prototypes
 *   kernels.h  — SPH kernel function prototypes
 *   io.h       — snapshot reader prototypes (HDF5 + Gadget binary)
 *   render.h   — rendering pipeline prototypes (deposit, output, etc.)
 *   colormap.h — render_config_t, palette and opacity enums
 *
 * This file includes all of them so existing source files that do
 * #include "header.h" continue to compile without modification.
 *
 * New code should include only the headers it actually needs:
 *
 *   kernels.c            -> types.h, kernels.h
 *   make_tree.c          -> types.h, tree.h
 *   walk_tree.c          -> types.h, tree.h, kernels.h
 *   find_neighbours.c    -> types.h, globals.h, tree.h, kernels.h
 *   io.c                 -> types.h, globals.h, io.h
 *   select_particles.c   -> types.h, render.h
 *   split_across_tasks.c -> types.h, globals.h, render.h
 *   smooth_to_mesh.c     -> types.h, globals.h, kernels.h, render.h
 *   deposit_sph_2d.c     -> types.h, kernels.h, render.h
 *   write_to_ppm.c       -> types.h, render.h, colormap.h
 *   render_image.c       -> types.h, globals.h, io.h, render.h, colormap.h
 *   header.c             -> types.h, globals.h, io.h
 */

#include "types.h"
#include "globals.h"
#include "tree.h"
#include "kernels.h"
#include "io.h"
#include "render.h"
#include "colormap.h"

#endif /* RENDER_IMAGE_HEADER_H */
