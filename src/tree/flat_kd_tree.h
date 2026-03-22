#ifndef FLAT_KD_TREE_H
#define FLAT_KD_TREE_H

/*
 * flat_kd_tree.h  —  cache-friendly flat-array kd-tree, v2
 *
 * Design goals vs v1:
 *   - Zero per-node malloc during build (in-place partial sort)
 *   - Iterative KNN search (explicit stack, no recursion overhead)
 *   - Nodes packed tightly: 32 bytes each, fits ~500K nodes in 16MB L3
 *   - Children at 2i+1, 2i+2 — no pointer chasing
 *
 * Build: O(N log N), ~1-2s for 6.65M particles
 * Search: O(log N) average, ~30 cycles/particle after L3 warm-up
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifndef KD_LEAF_SIZE
#define KD_LEAF_SIZE 32
#endif
#ifndef KD_NDIM
#define KD_NDIM 3
#endif

/* 32-byte node — two nodes per cache line */
typedef struct {
    float lo[KD_NDIM];   /* AABB min  12 bytes */
    int   start;         /*            4 bytes */
    float hi[KD_NDIM];   /* AABB max  12 bytes */
    int   count;         /*            4 bytes */
} kd_node_t;             /*           32 bytes total */

typedef struct {
    kd_node_t *nodes;
    float     *pos;      /* AoS: [x0,y0,z0, x1,y1,z1, ...], sorted in-place */
    int        n_nodes;
    int        n_pts;
} kd_tree_t;

/* ------------------------------------------------------------------ */
/* In-place partial sort along one dimension using introselect          */
/* Rearranges rows of pos[start..start+count) so that rows [start..    */
/* start+count/2) have dim-coord <= median and rest have >= median.    */
/* ------------------------------------------------------------------ */
static void kd_nth_row(float *pos, int start, int count, int dim)
{
    int lo = 0, hi = count - 1, k = count / 2;
    while (lo < hi) {
        /* Median-of-3 pivot */
        int mid = lo + (hi - lo) / 2;
        float plo = pos[(start+lo) *KD_NDIM+dim];
        float pmid= pos[(start+mid)*KD_NDIM+dim];
        float phi = pos[(start+hi) *KD_NDIM+dim];
        float pivot;
        if (plo <= pmid && pmid <= phi) pivot = pmid;
        else if (phi <= pmid && pmid <= plo) pivot = pmid;
        else if (pmid <= plo && plo <= phi) pivot = plo;
        else if (phi <= plo && plo <= pmid) pivot = plo;
        else pivot = phi;

        int i = lo, j = hi;
        while (i <= j) {
            while (pos[(start+i)*KD_NDIM+dim] < pivot) i++;
            while (pos[(start+j)*KD_NDIM+dim] > pivot) j--;
            if (i <= j) {
                /* Swap full rows */
                float *ra = pos + (start+i)*KD_NDIM;
                float *rb = pos + (start+j)*KD_NDIM;
                for (int d = 0; d < KD_NDIM; d++) {
                    float t = ra[d]; ra[d] = rb[d]; rb[d] = t;
                }
                i++; j--;
            }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else break;
    }
}

static void kd_compute_bbox(const float *pos, int start, int count,
                              float *lo, float *hi)
{
    for (int d = 0; d < KD_NDIM; d++) lo[d] = hi[d] = pos[start*KD_NDIM+d];
    for (int i = start+1; i < start+count; i++)
        for (int d = 0; d < KD_NDIM; d++) {
            float v = pos[i*KD_NDIM+d];
            if (v < lo[d]) lo[d] = v;
            if (v > hi[d]) hi[d] = v;
        }
}

/* Iterative build using an explicit stack of (node_idx, start, count) */
typedef struct { int node; int start; int count; } build_task_t;

static kd_tree_t *kd_build(float *pos, int n_pts)
{
    kd_tree_t *tree = (kd_tree_t *)malloc(sizeof(kd_tree_t));
    tree->n_pts = n_pts;
    tree->pos   = pos;

    /* Compute max depth and allocate node array */
    int max_depth = 1;
    { int n = n_pts; while (n > KD_LEAF_SIZE) { n = (n+1)/2; max_depth++; } }
    tree->n_nodes = 1 << max_depth;  /* complete binary tree upper bound */
    tree->nodes = (kd_node_t *)calloc(tree->n_nodes, sizeof(kd_node_t));

    /* Build stack */
    int stack_cap = max_depth * 4 + 16;
    build_task_t *stack = (build_task_t *)malloc(stack_cap * sizeof(build_task_t));
    int sp = 0;
    stack[sp++] = (build_task_t){0, 0, n_pts};

    while (sp > 0) {
        build_task_t t = stack[--sp];
        int ni = t.node, s = t.start, cnt = t.count;
        kd_node_t *node = &tree->nodes[ni];
        node->start = s;
        node->count = cnt;
        kd_compute_bbox(pos, s, cnt, node->lo, node->hi);

        if (cnt <= KD_LEAF_SIZE) continue;  /* leaf */

        /* Choose longest axis */
        int best_dim = 0;
        float best_span = node->hi[0] - node->lo[0];
        for (int d = 1; d < KD_NDIM; d++) {
            float span = node->hi[d] - node->lo[d];
            if (span > best_span) { best_span = span; best_dim = d; }
        }

        int left_cnt  = cnt / 2;
        int right_cnt = cnt - left_cnt;
        kd_nth_row(pos, s, cnt, best_dim);

        int left_ni  = 2*ni + 1;
        int right_ni = 2*ni + 2;

        /* Grow stack if needed */
        if (sp + 2 >= stack_cap) {
            stack_cap *= 2;
            stack = (build_task_t *)realloc(stack, stack_cap * sizeof(build_task_t));
        }
        stack[sp++] = (build_task_t){right_ni, s + left_cnt, right_cnt};
        stack[sp++] = (build_task_t){left_ni,  s,             left_cnt};
    }

    free(stack);
    return tree;
}

static void kd_free(kd_tree_t *tree)
{
    if (!tree) return;
    free(tree->nodes);
    free(tree);
}

/* ------------------------------------------------------------------ */
/* Max-heap for k nearest                                               */
/* ------------------------------------------------------------------ */
static inline void kd_heap_push(float *r2, int *count, int k, float val)
{
    int n = *count;
    if (n < k) {
        r2[n] = val;
        (*count)++;
        /* Sift up */
        int i = n;
        while (i > 0) {
            int p = (i-1)>>1;
            if (r2[p] >= r2[i]) break;
            float t = r2[p]; r2[p] = r2[i]; r2[i] = t;
            i = p;
        }
    } else if (val < r2[0]) {
        r2[0] = val;
        /* Sift down */
        int i = 0;
        while (1) {
            int l = 2*i+1, r = 2*i+2, lg = i;
            if (l < k && r2[l] > r2[lg]) lg = l;
            if (r < k && r2[r] > r2[lg]) lg = r;
            if (lg == i) break;
            float t = r2[i]; r2[i] = r2[lg]; r2[lg] = t;
            i = lg;
        }
    }
}

static inline float kd_aabb_min_r2(const float *lo, const float *hi,
                                     const float *q)
{
    float r2 = 0.0f;
    for (int d = 0; d < KD_NDIM; d++) {
        float v = q[d];
        float diff = v < lo[d] ? lo[d]-v : v > hi[d] ? v-hi[d] : 0.0f;
        r2 += diff*diff;
    }
    return r2;
}

/* ------------------------------------------------------------------ */
/* Iterative KNN search — no recursion, explicit stack                  */
/* ------------------------------------------------------------------ */
#define KD_STACK_SIZE 128   /* depth 18 tree needs at most 36 entries */

static float kd_knn(const kd_tree_t *tree, const float *q,
                     int k, float *heap_r2)
{
    int heap_count = 0;
    float worst = 1e30f;

    /* Iterative traversal stack: stores node indices */
    int stack[KD_STACK_SIZE];
    int sp = 0;
    stack[sp++] = 0;

    const kd_node_t *nodes = tree->nodes;
    const float     *pos   = tree->pos;
    int              n_nodes = tree->n_nodes;

    while (sp > 0) {
        int ni = stack[--sp];
        if (ni >= n_nodes) continue;
        const kd_node_t *node = &nodes[ni];
        if (node->count == 0) continue;

        /* Prune */
        float min_r2 = kd_aabb_min_r2(node->lo, node->hi, q);
        if (heap_count == k && min_r2 >= worst) continue;

        int left  = 2*ni + 1;
        int right = 2*ni + 2;
        int is_leaf = (left >= n_nodes ||
                       nodes[left].count == 0) &&
                      (right >= n_nodes ||
                       nodes[right].count == 0);

        if (is_leaf || node->count <= KD_LEAF_SIZE) {
            /* Leaf: check all particles */
            const float *p = pos + (size_t)node->start * KD_NDIM;
            int cnt = node->count;
            for (int i = 0; i < cnt; i++, p += KD_NDIM) {
                float r2 = 0.0f;
                for (int d = 0; d < KD_NDIM; d++) {
                    float diff = q[d] - p[d]; r2 += diff*diff;
                }
                kd_heap_push(heap_r2, &heap_count, k, r2);
            }
            if (heap_count == k) worst = heap_r2[0];
        } else {
            /* Push farther child first (so nearer child is processed first) */
            float dl = (left  < n_nodes && nodes[left ].count > 0)
                        ? kd_aabb_min_r2(nodes[left ].lo, nodes[left ].hi, q)
                        : 1e30f;
            float dr = (right < n_nodes && nodes[right].count > 0)
                        ? kd_aabb_min_r2(nodes[right].lo, nodes[right].hi, q)
                        : 1e30f;
            if (sp + 2 >= KD_STACK_SIZE) {
                /* Stack overflow safety: fall back to nearer child only */
                int near = (dl <= dr) ? left : right;
                if (near < n_nodes && nodes[near].count > 0)
                    stack[sp++] = near;
            } else if (dl <= dr) {
                if (right < n_nodes && nodes[right].count > 0) stack[sp++] = right;
                if (left  < n_nodes && nodes[left ].count > 0) stack[sp++] = left;
            } else {
                if (left  < n_nodes && nodes[left ].count > 0) stack[sp++] = left;
                if (right < n_nodes && nodes[right].count > 0) stack[sp++] = right;
            }
        }
    }

    return (heap_count >= k) ? heap_r2[0] : -1.0f;
}

#endif /* FLAT_KD_TREE_H */
