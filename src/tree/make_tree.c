#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

/* ------------------------------------------------------------------
 * octant_index: compute which child octant x belongs to, and narrow
 * the bounding box to that octant.  Returns a bitmask k where bit i
 * is set when coordinate i is in the upper half of the cell.
 * ------------------------------------------------------------------ */
static inline int octant_index(const float *x, const float *xmid,
                                float *node_xmin, float *node_xmax,
                                const float *xmin, const float *xmax)
{
    int k = 0;
    for (int i = 0; i < NDIM; i++) {
        node_xmin[i] = xmin[i];
        node_xmax[i] = xmax[i];
        if (x[i] > xmid[i]) {
            k += (1 << i);
            node_xmin[i] = xmid[i];
        } else {
            node_xmax[i] = xmid[i];
        }
    }
    return k;
}

/* ------------------------------------------------------------------
 * make_tree: entry point.  Computes bounding box then inserts every
 * point into the octree rooted at *root.
 * ------------------------------------------------------------------ */
void make_tree(struct point *P, int num_pts, struct tree_node_nd **root)
{
    int i, j;
    float x[NDIM];
    float xmin[NDIM], xmax[NDIM];
    struct tree_node_nd *root_nd;

    for (j = 0; j < NDIM; j++) {
        xmin[j] =  1.e10f;
        xmax[j] = -1.e10f;
    }

    for (i = 0; i < num_pts; i++)
        for (j = 0; j < NDIM; j++) {
            if (xmin[j] > P[i].pos[j]) xmin[j] = P[i].pos[j];
            if (xmax[j] < P[i].pos[j]) xmax[j] = P[i].pos[j];
        }

    for (j = 0; j < NDIM; j++)
        printf("Extent: xmin[%d]: %g, xmax[%d]: %g\n", j, xmin[j], j, xmax[j]);

    root_nd = NULL;
    for (i = 0; i < num_pts; i++) {
        memcpy(x, &(P[i].pos[0]), NDIM * sizeof(float));
        root_nd = add_to_node_nd(root_nd, x, xmin, xmax);
    }

    *root = root_nd;
}

/* ------------------------------------------------------------------
 * add_to_node: 1-D binary tree (kept for legacy callers).
 * ------------------------------------------------------------------ */
struct tree_node *add_to_node(struct tree_node *node, float x,
                               float xmin, float xmax)
{
    int i;
    float xmid = 0.5f * (xmin + xmax);

    if (node == NULL) {
        node = (struct tree_node *)malloc(sizeof(struct tree_node));
        node->num_members = 1;
        node->x[0] = x;
        node->xmin  = xmin;
        node->xmax  = xmax;
        node->split = 0;
        node->left  = node->right = NULL;
    } else {
        if (node->split == 0) {
            node->num_members++;
            if (node->num_members > MAXNODE) {
                node->split = 1;
                for (i = 0; i < MAXNODE; i++) {
                    if (node->x[i] <= xmid)
                        node->left  = add_to_node(node->left,  node->x[i], xmin, xmid);
                    else
                        node->right = add_to_node(node->right, node->x[i], xmid, xmax);
                }
                if (x <= xmid)
                    node->left  = add_to_node(node->left,  x, xmin, xmid);
                else
                    node->right = add_to_node(node->right, x, xmid, xmax);
            } else {
                node->x[node->num_members - 1] = x;
            }
        } else {
            node->num_members++;
            if (x <= xmid)
                node->left  = add_to_node(node->left,  x, xmin, xmid);
            else
                node->right = add_to_node(node->right, x, xmid, xmax);
        }
    }
    return node;
}

/* ------------------------------------------------------------------
 * add_to_node_nd: insert point x into the N-D octree.
 *
 * Changes vs. original:
 *  - children[k] replaces p0..p7 (single array, indexed by bitmask).
 *  - octant_index() eliminates the four copies of the bitmask logic.
 *  - xmean is updated with a correct incremental-mean formula.
 * ------------------------------------------------------------------ */
struct tree_node_nd *add_to_node_nd(struct tree_node_nd *node, float *x,
                                     float *xmin, float *xmax)
{
    int i, k;
    float xmid[NDIM];
    float y[NDIM];
    float node_xmin[NDIM], node_xmax[NDIM];

    for (i = 0; i < NDIM; i++)
        xmid[i] = 0.5f * (xmin[i] + xmax[i]);

#ifdef VERBOSE
    for (i = 0; i < NDIM; i++)
        printf("(xmin, xmax)[%d]: (%g,%g) x[%d]: %g\n",
               i, xmin[i], xmax[i], i, x[i]);
    printf("\n");
#endif

    if (node == NULL) {
        /* --- create new leaf node --- */
        node = (struct tree_node_nd *)malloc(sizeof(struct tree_node_nd));
        node->num_members = 1;
        for (i = 0; i < NDIM; i++) {
            node->x[i][0] = x[i];
            node->xmin[i] = xmin[i];
            node->xmax[i] = xmax[i];
            node->xmean[i] = x[i];
        }
        node->split = 0;
        for (i = 0; i < NCHILDREN; i++)
            node->children[i] = NULL;

    } else if (node->split == 0) {
        /* --- leaf node: try to add to bucket --- */
        node->num_members++;

        if (node->num_members > MAXNODE) {
            /* Bucket full: split and redistribute existing points */
            node->split = 1;

            for (i = 0; i < MAXNODE; i++) {
                /* extract the i-th stored point */
                for (int j = 0; j < NDIM; j++)
                    y[j] = node->x[j][i];

                k = octant_index(y, xmid, node_xmin, node_xmax, xmin, xmax);
                node->children[k] = add_to_node_nd(node->children[k],
                                                    y, node_xmin, node_xmax);
            }

            /* Now insert the new point */
            k = octant_index(x, xmid, node_xmin, node_xmax, xmin, xmax);
            node->children[k] = add_to_node_nd(node->children[k],
                                                x, node_xmin, node_xmax);
        } else {
            /* Still room in the bucket */
            int m = node->num_members - 1;
            for (i = 0; i < NDIM; i++) {
                node->x[i][m] = x[i];
                /* Correct incremental mean: mean_n = mean_{n-1} + (x - mean_{n-1})/n */
                node->xmean[i] += (x[i] - node->xmean[i]) / (float)node->num_members;
            }
        }

    } else {
        /* --- internal node: route to the correct child --- */
        node->num_members++;
        k = octant_index(x, xmid, node_xmin, node_xmax, xmin, xmax);
        node->children[k] = add_to_node_nd(node->children[k],
                                            x, node_xmin, node_xmax);
    }

    return node;
}

