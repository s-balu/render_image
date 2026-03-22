#ifndef RENDER_IMAGE_TREE_H
#define RENDER_IMAGE_TREE_H

/*
 * tree.h — octree build and walk function prototypes.
 *
 * Struct definitions (tree_node, tree_node_nd, etc.) live in types.h.
 * Include this header in files that build or search the octree:
 *   make_tree.c, walk_tree.c, find_neighbours.c, kernels.c
 *
 * render_image.c and smooth_to_mesh.c do NOT need this header — they
 * call find_neighbours() which owns the tree internally.
 */

#include "types.h"

/* ------------------------------------------------------------------ */
/* Tree building (make_tree.c)                                          */
/* ------------------------------------------------------------------ */

struct tree_node    *add_to_node(struct tree_node *, float, float, float);
struct tree_node_nd *add_to_node_nd(struct tree_node_nd *, float *,
                                     float *, float *);
void make_tree(struct point *, int, struct tree_node_nd **);

/* ------------------------------------------------------------------ */
/* Tree walking (walk_tree.c)                                           */
/* ------------------------------------------------------------------ */

/* Flat-buffer neighbour search — no malloc per neighbour, SIMD-friendly */
void get_multiple_nodes_nd_flat(struct tree_node_nd *, const float *,
                                 float, ngb_buf_t *);

/* Potential/density estimation */
void get_distance_to_nth_nearest_neighbour(struct point *, int, int,
                                            struct tree_node_nd *);
void get_kernel_density_estimate(struct point *, int, struct tree_node_nd *);
void get_potential_estimate(struct point *, int, struct tree_node_nd *);

/* Internal walk functions (used within walk_tree.c) */
void sweep_over_nodes_nd(struct tree_node_nd *, float *, float,
                          int *, int *);
void add_node_index_nd(struct tree_node_nd *, int *);
void build_interaction_list(struct tree_node_nd *,
                             struct interaction_list **, int *);

/* Legacy 1-D tree (kept for compatibility) */
void get_node(struct tree_node *, float, int *, struct link_list **);
void get_multiple_nodes(struct tree_node *, float, float,
                         int *, struct link_list **);
void scan_nodes(struct tree_node *, int *);

/* Legacy linked-list interface */
void get_node_nd(struct tree_node_nd *, float *, int *, struct link_list **);
void get_multiple_nodes_nd(struct tree_node_nd *, float *, float,
                            int *, int *, struct link_list **);
void free_link_list(struct link_list *);

/* Utility */
void get_points(struct point *, int *, int);
int  cmpfunc(const void *, const void *);
float get_rand(int);

#endif /* RENDER_IMAGE_TREE_H */
