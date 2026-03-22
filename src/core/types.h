#ifndef RENDER_IMAGE_TYPES_H
#define RENDER_IMAGE_TYPES_H

/*
 * types.h — Plain-old-data types shared across the codebase.
 *
 * Rules for this file:
 *   - structs and typedefs only
 *   - #define constants only
 *   - NO function declarations
 *   - NO extern globals
 *   - NO includes of project headers (only standard library types if needed)
 *
 * Every other header may safely include this one.  Changing a struct here
 * will trigger recompilation of everything that includes it — keep structs
 * stable and add new fields at the end.
 */

#include <stdint.h>
#include <stddef.h>   /* size_t */

/* ------------------------------------------------------------------ */
/* Spatial constants                                                    */
/* ------------------------------------------------------------------ */

#define MAXNODE    8
#define NDIM       3
#define NCHILDREN  (1 << NDIM)   /* children per octree node: 2^NDIM = 8 */

/* Upper bound on particle types in the Gadget binary header block.
 * The actual number of types present is stored in sim_info.num_types. */
#define GADGET_MAX_TYPES  6

/* Maximum neighbours in a flat buffer.  Increase for very dense fields. */
#define MAX_NGP  200000

/* ------------------------------------------------------------------ */
/* Tree structures                                                      */
/* ------------------------------------------------------------------ */

struct link_list {
    float x[NDIM];
    struct link_list *ptr;
};

struct interaction_list {
    int index;
    struct tree_node_nd *node;
    struct interaction_list *left, *right;
};

struct tree_node {
    float x[MAXNODE];
    float xmin, xmax;
    int   num_members;
    int   i;
    int   split;
    struct tree_node *left, *right;
};

/* N-dimensional octree node.
 * children[k] uses a bitmask: bit i set => coord i is in the upper half-cell. */
struct tree_node_nd {
    float x[NDIM][MAXNODE];
    float xmin[NDIM], xmax[NDIM];
    float xmean[NDIM];
    int   num_members;
    int   index;
    int   split;
    struct tree_node_nd *children[NCHILDREN];
};

struct btree_node {
    float x;
    struct btree_node *left, *right;
};

/* ------------------------------------------------------------------ */
/* Particle / simulation types                                          */
/* ------------------------------------------------------------------ */

struct point {
    float pos[NDIM];
    float mass;
    float density;
    float dist_ngb;
};

/*
 * sim_info: everything read from a snapshot header.
 *
 * npart / nall / massarr are heap-allocated arrays of length num_types.
 * Call sim_info_alloc() after reading num_types, sim_info_free() when done.
 */
struct sim_info {
    int    num_types;      /* actual number of particle types in this file  */
    int   *npart;          /* [num_types] particle counts in this file      */
    int   *nall;           /* [num_types] total particle counts (all files) */
    double *massarr;       /* [num_types] table masses (0 = per-particle)   */
    double time;
    double redshift;
    int    NumFiles;
    double BoxSize;
    double Omega0, OmegaLambda, HubbleParam;
    int    SnapFormat;
};

/* ------------------------------------------------------------------ */
/* Neighbour buffer                                                     */
/* ------------------------------------------------------------------ */

/*
 * ngb_buf_t: flat heap-allocated buffer replacing malloc'd linked lists.
 *
 * SoA layout (separate x0/x1/x2) lets the distance loop auto-vectorise.
 * r2[] is filled during the tree walk so no second distance pass is needed.
 * Allocate once per thread with ngb_buf_alloc(num_ngb).
 */
typedef struct {
    float *x0;       /* x-coordinates of candidate neighbours */
    float *x1;       /* y-coordinates                         */
#if NDIM == 3
    float *x2;       /* z-coordinates                         */
#endif
    float *r2;       /* squared distances (filled in-place)   */
    int    count;
    int    capacity;
} ngb_buf_t;

/* ------------------------------------------------------------------ */
/* Image types                                                          */
/* ------------------------------------------------------------------ */

typedef struct { uint8_t red; uint8_t green; uint8_t blue; } pixel_t;
typedef struct { pixel_t *pixels; size_t width; size_t height; }  bitmap_t;

/* ------------------------------------------------------------------ */
/* SPH deposit types                                                    */
/* ------------------------------------------------------------------ */

/* Packed AoS particle record for deposit_sph_2d.
 * Stores pre-projected pixel coordinates and h in pixels so the
 * deposit loop reads all three values from the same cache line. */
typedef struct { float px; float py; float h; } deposit_particle_t;

#endif /* RENDER_IMAGE_TYPES_H */
