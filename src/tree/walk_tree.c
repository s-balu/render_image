#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef ENABLE_OPENMP
#include "omp.h"
#endif

#include "header.h"

/* Forward declaration: sweep_over_nodes_nd is defined later in this file
 * but called by get_potential_estimate which appears first. */
void sweep_over_nodes_nd(struct tree_node_nd *node, float *x, float theta,
                          int *count, int *ncount);

#include "header.h"

/* ------------------------------------------------------------------ */
/* ngb_buf_t lifecycle                                                  */
/* ------------------------------------------------------------------ */

ngb_buf_t *ngb_buf_alloc(int num_ngb)
{
    /*
     * Size the buffer so the sphere filter (which limits actual fill to
     * ~2*num_ngb) has plenty of headroom without wasting cache.
     * 8*num_ngb bytes per float array × 4 arrays = 32*num_ngb bytes total.
     * For num_ngb=32: 32*32*4 = 4 KB — fits easily in L1.
     * For num_ngb=64: 8 KB.
     */
    int cap = 8 * num_ngb;
    if (cap < 64) cap = 64;   /* minimum sanity */

    ngb_buf_t *buf = (ngb_buf_t *)malloc(sizeof(ngb_buf_t));
    if (!buf) return NULL;

    buf->x0 = (float *)malloc(cap * sizeof(float));
    buf->x1 = (float *)malloc(cap * sizeof(float));
#if NDIM == 3
    buf->x2 = (float *)malloc(cap * sizeof(float));
#endif
    buf->r2 = (float *)malloc(cap * sizeof(float));

    if (!buf->x0 || !buf->x1 || !buf->r2
#if NDIM == 3
        || !buf->x2
#endif
    ) {
        ngb_buf_free(buf);
        return NULL;
    }

    buf->count    = 0;
    buf->capacity = cap;
    return buf;
}

void ngb_buf_free(ngb_buf_t *buf)
{
    if (!buf) return;
    free(buf->x0);
    free(buf->x1);
#if NDIM == 3
    free(buf->x2);
#endif
    free(buf->r2);
    free(buf);
}

/* ------------------------------------------------------------------ */
/* Comparator                                                           */
/* ------------------------------------------------------------------ */
int cmpfunc(const void *elem1, const void *elem2)
{
    float a = *(const float *)elem1;
    float b = *(const float *)elem2;
    return (a > b) - (a < b);
}

/* ------------------------------------------------------------------ */
/* free_link_list                                                       */
/* ------------------------------------------------------------------ */
void free_link_list(struct link_list *ngb)
{
    struct link_list *tmp;
    while (ngb != NULL) {
        tmp = ngb->ptr;
        free(ngb);
        ngb = tmp;
    }
}

/* ------------------------------------------------------------------ */
/* octant_index_walk                                                    */
/* ------------------------------------------------------------------ */
static inline int octant_index_walk(const float *x, const float *xmid,
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

/* ================================================================== */
/* get_multiple_nodes_nd_flat                                           */
/*                                                                      */
/* Sphere-filtered range search into a heap-allocated buffer.          */
/* The buffer is grown if it's full rather than silently clamping.     */
/* ================================================================== */
void get_multiple_nodes_nd_flat(struct tree_node_nd *node, const float *x,
                                  float search, ngb_buf_t *buf)
{
    if (node == NULL) return;

    float r2_limit = search * search;

    /* AABB rejection */
#if NDIM == 3
    if ((x[0] - search) >= node->xmax[0] || (x[0] + search) <= node->xmin[0] ||
        (x[1] - search) >= node->xmax[1] || (x[1] + search) <= node->xmin[1] ||
        (x[2] - search) >= node->xmax[2] || (x[2] + search) <= node->xmin[2])
        return;
#elif NDIM == 2
    if ((x[0] - search) >= node->xmax[0] || (x[0] + search) <= node->xmin[0] ||
        (x[1] - search) >= node->xmax[1] || (x[1] + search) <= node->xmin[1])
        return;
#else
    if ((x[0] - search) >= node->xmax[0] || (x[0] + search) <= node->xmin[0])
        return;
#endif

    if (node->split == 0) {
        int nm = node->num_members;
        int c  = buf->count;

        for (int i = 0; i < nm; i++) {
#if NDIM == 3
            float d0 = x[0] - node->x[0][i];
            float d1 = x[1] - node->x[1][i];
            float d2 = x[2] - node->x[2][i];
            float r2 = d0*d0 + d1*d1 + d2*d2;
#elif NDIM == 2
            float d0 = x[0] - node->x[0][i];
            float d1 = x[1] - node->x[1][i];
            float r2 = d0*d0 + d1*d1;
#else
            float d0 = x[0] - node->x[0][i];
            float r2 = d0*d0;
#endif
            if (r2 > r2_limit) continue;

            /* Grow buffer if needed — shouldn't happen often with correct base_f */
            if (c >= buf->capacity) {
                int newcap = buf->capacity * 2;
                float *nx0 = (float *)realloc(buf->x0, newcap * sizeof(float));
                float *nx1 = (float *)realloc(buf->x1, newcap * sizeof(float));
                float *nr2 = (float *)realloc(buf->r2, newcap * sizeof(float));
#if NDIM == 3
                float *nx2 = (float *)realloc(buf->x2, newcap * sizeof(float));
                if (!nx0 || !nx1 || !nx2 || !nr2) {
                    fprintf(stderr, "ngb_buf realloc failed\n"); break;
                }
                buf->x2 = nx2;
#else
                if (!nx0 || !nx1 || !nr2) {
                    fprintf(stderr, "ngb_buf realloc failed\n"); break;
                }
#endif
                buf->x0 = nx0; buf->x1 = nx1; buf->r2 = nr2;
                buf->capacity = newcap;
            }

            buf->x0[c] = node->x[0][i];
            buf->x1[c] = node->x[1][i];
#if NDIM == 3
            buf->x2[c] = node->x[2][i];
#endif
            buf->r2[c] = r2;
            c++;
        }
        buf->count = c;
    } else {
        for (int i = 0; i < NCHILDREN; i++)
            get_multiple_nodes_nd_flat(node->children[i], x, search, buf);
    }
}

/* ================================================================== */
/* nth_element_float — introselect                                      */
/* ================================================================== */
static void nth_element_float(float *arr, int n, int k)
{
    int left = 0, right = n - 1;
    int limit = n * 2;

    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid]   < arr[left])  { float t = arr[mid];   arr[mid]   = arr[left];  arr[left]  = t; }
        if (arr[right] < arr[left])  { float t = arr[right]; arr[right] = arr[left];  arr[left]  = t; }
        if (arr[right] < arr[mid])   { float t = arr[right]; arr[right] = arr[mid];   arr[mid]   = t; }

        float pivot = arr[mid];
        { float t = arr[mid]; arr[mid] = arr[right-1]; arr[right-1] = t; }

        int i = left, j = right - 1;
        while (1) {
            while (arr[++i] < pivot) {}
            while (arr[--j] > pivot) {}
            if (i >= j) break;
            float t = arr[i]; arr[i] = arr[j]; arr[j] = t;
        }
        { float t = arr[i]; arr[i] = arr[right]; arr[right] = t; }

        if (i <= k) left  = i + 1;
        if (i >= k) right = i - 1;

        if (--limit == 0) {
            for (int a = left + 1; a <= right; a++) {
                float v = arr[a]; int b = a;
                while (b > left && arr[b-1] > v) { arr[b] = arr[b-1]; b--; }
                arr[b] = v;
                if (b <= k && k <= a) return;
            }
            return;
        }
    }
}

/* ================================================================== */
/* get_distance_to_nth_nearest_neighbour                                */
/* ================================================================== */
void get_distance_to_nth_nearest_neighbour(struct point *P, int num_pts,
                                            int num_ngb,
                                            struct tree_node_nd *root)
{
    int i;

    /*
     * Estimate mean inter-particle separation.
     *
     * The root cell volume can be highly non-cubic (e.g. x-span 6, z-span 100)
     * when the view box doesn't span the full simulation depth.  Using the
     * full root volume would give a mean_sep 8-16x too large, causing the
     * initial search sphere to collect thousands of extra candidates.
     *
     * Instead we use the MINIMUM axis span to estimate the local density.
     * This gives a search radius appropriate for the densest dimension
     * (the x-slab), which is what determines the 32nd-nearest-neighbour
     * distance in practice.
     */
    float min_span = root->xmax[0] - root->xmin[0];
    for (i = 1; i < NDIM; i++) {
        float span = root->xmax[i] - root->xmin[i];
        if (span < min_span) min_span = span;
    }
    /* Local density estimate: num_pts particles in a cube of side min_span */
    float mean_sep = min_span * powf((float)num_pts, -1.0f / (float)NDIM);
#if NDIM == 2
    mean_sep *= (float)(M_PI / 4.0);
#elif NDIM == 3
    mean_sep *= (float)(M_PI / 6.0);
#endif

#if NDIM == 3
    float base_f = powf(3.0f * 2.0f * (float)num_ngb / (4.0f * (float)M_PI), 1.0f/3.0f);
#elif NDIM == 2
    float base_f = sqrtf(2.0f * (float)num_ngb / (float)M_PI);
#else
    float base_f = 2.0f * (float)num_ngb;
#endif

    const float WARMSTART_CAP = 4.0f;

    printf("Using a base search radius of %g (mean separation = %g)...\n",
           base_f * mean_sep, mean_sep);

#ifdef ENABLE_OPENMP
    int NThreads = omp_get_max_threads();
    printf("Using %d OpenMP threads for neighbour search...\n", NThreads);
    fflush(stdout);

#pragma omp parallel num_threads(NThreads) private(i)
    {
        int tid     = omp_get_thread_num();
        int i_start = (long long)tid       * num_pts / NThreads;
        int i_end   = (long long)(tid + 1) * num_pts / NThreads;
#else
    {
        int i_start = 0, i_end = num_pts;
#endif

        ngb_buf_t *buf = ngb_buf_alloc(num_ngb);
        if (!buf) {
            fprintf(stderr, "Failed to allocate neighbour buffer\n");
#ifdef ENABLE_OPENMP
            goto thread_done;
#else
            return;
#endif
        }

        float fsearch = base_f;

        for (i = i_start; i < i_end; i++) {
            float x0 = P[i].pos[0];
#if NDIM >= 2
            float x1 = P[i].pos[1];
#endif
#if NDIM == 3
            float x2 = P[i].pos[2];
            float xarr[3] = {x0, x1, x2};
#elif NDIM == 2
            float xarr[2] = {x0, x1};
#else
            float xarr[1] = {x0};
#endif

            /* MAX_EXPAND: particles needing more expansions than this
             * have no real neighbours in the view volume (e.g. z-boundary
             * particles from a deep periodic box projected onto a thin slice).
             * Mark them with dist_ngb = -1 so find_neighbours can exclude them. */
            #define MAX_EXPAND 4
            int attempts = 0;
            int isolated = 0;
            while (1) {
                buf->count = 0;
                get_multiple_nodes_nd_flat(root, xarr, fsearch * mean_sep, buf);
                if (buf->count > num_ngb) break;
                fsearch *= 1.5f;
                if (++attempts > MAX_EXPAND) {
                    isolated = 1;
                    break;
                }
            }

            if (isolated) {
                P[i].dist_ngb = -1.0f;
                /* Reset fsearch so the next particle starts fresh */
                fsearch = base_f;
                continue;
            }

            int nc = buf->count;
            int k  = (num_ngb < nc) ? num_ngb : nc - 1;

            if (nc > num_ngb)
                nth_element_float(buf->r2, nc, k);

            P[i].dist_ngb = sqrtf(buf->r2[k]);

            float next_f = 1.3f * P[i].dist_ngb / mean_sep;
            if (next_f < base_f * 0.5f)         next_f = base_f * 0.5f;
            if (next_f > base_f * WARMSTART_CAP) next_f = base_f * WARMSTART_CAP;
            fsearch = next_f;
        }

        ngb_buf_free(buf);
#ifdef ENABLE_OPENMP
        thread_done:;
    }
#else
    }
#endif
}

#ifdef KERNEL_SMOOTHING
/* ================================================================== */
/* get_kernel_density_estimate                                          */
/* ================================================================== */
void get_kernel_density_estimate(struct point *P, int num_pts,
                                  struct tree_node_nd *root)
{
    ngb_buf_t *buf = ngb_buf_alloc(64);

    for (int i = 0; i < num_pts; i++) {
        float x0 = P[i].pos[0];
#if NDIM >= 2
        float x1 = P[i].pos[1];
#endif
#if NDIM == 3
        float x2 = P[i].pos[2];
        float xarr[3] = {x0, x1, x2};
#elif NDIM == 2
        float xarr[2] = {x0, x1};
#else
        float xarr[1] = {x0};
#endif

        float h    = P[i].dist_ngb;
        float hfac = 1.0f / h;
#if NDIM == 2
        hfac /= h;
#elif NDIM == 3
        hfac /= h * h;
#endif

        buf->count = 0;
        get_multiple_nodes_nd_flat(root, xarr, h, buf);

        double density = 0.0;
        for (int j = 0; j < buf->count; j++) {
            float xvar = sqrtf(buf->r2[j]) / h;
            density += P[i].mass * get_kernel_value(xvar) * hfac;
        }
        P[i].density = (float)density;
    }
    ngb_buf_free(buf);
}
#endif

/* ================================================================== */
/* get_potential_estimate                                               */
/* ================================================================== */
void get_potential_estimate(struct point *P, int num_pts,
                             struct tree_node_nd *root)
{
    int count = 0, ncount = 0;
    float x[NDIM];
    memcpy(x, &(P[0].pos[0]), NDIM * sizeof(float));
    sweep_over_nodes_nd(root, x, 0.3f, &count, &ncount);
    printf("%d %d\n", count, ncount);
    (void)num_pts;
}

/* ================================================================== */
/* 1-D walkers                                                          */
/* ================================================================== */
void scan_nodes(struct tree_node *node, int *count)
{
    if (node != NULL) {
        scan_nodes(node->left, count); (*count)++; scan_nodes(node->right, count);
    }
}

void get_node(struct tree_node *node, float x, int *count, struct link_list **ngb)
{
    struct link_list *tmp;
    if (node == NULL) return;
    float xmid = 0.5f * (node->xmin + node->xmax);
    struct tree_node *child = (x <= xmid) ? node->left : node->right;
    if (child) {
        get_node(child, x, count, ngb);
    } else {
        *count += node->num_members;
        for (int i = 0; i < node->num_members; i++) {
            tmp = (struct link_list *)malloc(sizeof(struct link_list));
            tmp->x[0] = node->x[i]; tmp->ptr = *ngb; *ngb = tmp;
        }
    }
}

void get_multiple_nodes(struct tree_node *node, float x, float search,
                         int *count, struct link_list **ngb)
{
    struct link_list *tmp;
    if (node == NULL) return;
    if ((x - search) < node->xmax && (x + search) > node->xmin) {
        get_multiple_nodes(node->left, x, search, count, ngb);
        if (node->split == 0) {
            *count += node->num_members;
            for (int i = 0; i < node->num_members; i++) {
                tmp = (struct link_list *)malloc(sizeof(struct link_list));
                tmp->x[0] = node->x[i]; tmp->ptr = *ngb; *ngb = tmp;
            }
        }
        get_multiple_nodes(node->right, x, search, count, ngb);
    }
}

/* ================================================================== */
/* N-D utilities                                                        */
/* ================================================================== */
void add_node_index_nd(struct tree_node_nd *node, int *count)
{
    if (node == NULL) return;
    node->index = ++(*count);
    for (int i = 0; i < NCHILDREN; i++)
        add_node_index_nd(node->children[i], count);
}

void get_node_nd(struct tree_node_nd *node, float *x, int *count,
                  struct link_list **ngb)
{
    float xmid[NDIM], node_xmin[NDIM], node_xmax[NDIM];
    struct link_list *tmp;
    if (node == NULL) return;
    for (int i = 0; i < NDIM; i++) xmid[i] = 0.5f*(node->xmin[i]+node->xmax[i]);
    int k = octant_index_walk(x, xmid, node_xmin, node_xmax, node->xmin, node->xmax);
    if (node->children[k]) {
        get_node_nd(node->children[k], x, count, ngb);
    } else {
        *count += node->num_members;
        for (int i = 0; i < node->num_members; i++) {
            tmp = (struct link_list *)malloc(sizeof(struct link_list));
            for (int j = 0; j < NDIM; j++) tmp->x[j] = node->x[j][i];
            tmp->ptr = *ngb; *ngb = tmp;
        }
    }
}

void get_multiple_nodes_nd(struct tree_node_nd *node, float *x, float search,
                            int *count, int *ncount, struct link_list **ngb)
{
    struct link_list *tmp;
    if (node == NULL) return;
    for (int i = 0; i < NDIM; i++)
        if ((x[i]-search) >= node->xmax[i] || (x[i]+search) <= node->xmin[i]) return;
    if (node->split == 0) {
        (*ncount)++;
        *count += node->num_members;
        for (int i = 0; i < node->num_members; i++) {
            tmp = (struct link_list *)malloc(sizeof(struct link_list));
            for (int j = 0; j < NDIM; j++) tmp->x[j] = node->x[j][i];
            tmp->ptr = *ngb; *ngb = tmp;
        }
    } else {
        for (int i = 0; i < NCHILDREN; i++)
            get_multiple_nodes_nd(node->children[i], x, search, count, ncount, ngb);
    }
}

void build_interaction_list(struct tree_node_nd *node,
                             struct interaction_list **list, int *index)
{
    struct interaction_list *tmp;
    if (node == NULL) return;
    build_interaction_list(node->children[0], list, index);
    tmp = (struct interaction_list *)malloc(sizeof(struct interaction_list));
    tmp->index = (*index)++; tmp->node = node; tmp->left = tmp->right = NULL;
    if (!*list) { *list = tmp; } else { (*list)->left = tmp; tmp->right = *list; *list = tmp; }
    for (int i = 1; i < NCHILDREN; i++)
        build_interaction_list(node->children[i], list, index);
}

void sweep_over_nodes_nd(struct tree_node_nd *node, float *x, float theta,
                          int *count, int *ncount)
{
    float xmid[NDIM], xmean[NDIM], dist_to_cell, cell_size, cell_offset;
    if (node == NULL) return;
    cell_size = 1.0f; cell_offset = 0.0f;
    for (int i = 0; i < NDIM; i++) {
        xmid[i]  = 0.5f*(node->xmin[i]+node->xmax[i]);
        xmean[i] = node->xmean[i];
        cell_size   *= (node->xmax[i]-node->xmin[i]);
        cell_offset += (xmean[i]-xmid[i])*(xmean[i]-xmid[i]);
    }
    cell_size   = powf(cell_size, 1.0f/(float)NDIM);
    cell_offset = sqrtf(cell_offset);
    dist_to_cell = 0.0f;
    for (int i = 0; i < NDIM; i++)
        dist_to_cell += (x[i]-xmean[i])*(x[i]-xmean[i]);
    if (dist_to_cell < (cell_size/theta + cell_offset)) {
        for (int i = 0; i < NCHILDREN; i++)
            sweep_over_nodes_nd(node->children[i], x, theta, count, ncount);
        printf("%g %g %g %g %d\n", x[0],x[1],xmid[0],xmid[1],node->num_members);
        *count += node->num_members;
    }
}
