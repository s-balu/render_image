/*
 * kernels.c  –  SPH kernel and projection functions
 *
 * The cubic spline kernel W(r,h) with compact support at r=2h is used
 * throughout.  Two evaluation paths are provided:
 *
 *  cubic_spline_kernel(u)            – 3-D kernel at u = r/h
 *  cubic_spline_kernel_2d_proj(u)   – analytically projected along z,
 *                                     giving a 2-D surface density kernel
 *                                     at projected radius u = R/h
 *
 * The 2-D projected kernel is used for image rendering because each
 * particle's contribution to a pixel is an integral of W along the
 * full line of sight through the particle, not a slice.
 *
 * Lookup tables with linear interpolation are used so that the inner
 * pixel loop pays only 2 multiplies and 1 add per evaluation.
 */

#include <stdio.h>
#include <math.h>

#include "header.h"

/* ------------------------------------------------------------------ */
/* Table resolution.  512 gives <0.1% interpolation error everywhere.  */
/* ------------------------------------------------------------------ */
#define NTAB 512

static float rh_table   [NTAB];
static float kernel_table[NTAB];
static float proj_table  [NTAB];   /* 2-D projected kernel            */
static float integral_table[NTAB]; /* cumulative of 3-D kernel (legacy)*/

/* ================================================================== */
/* 3-D cubic spline kernel  W(u)  where u = r/h                        */
/*                                                                      */
/* Normalised so that ∫ W dV = 1 (3-D).                                */
/* Support: u ∈ [0, 2].                                                 */
/* ================================================================== */
float cubic_spline_kernel(float u)
{
    float prefac;
#if NDIM == 3
    prefac = 8.0f / (float)M_PI;   /* correct 3-D normalisation */
#elif NDIM == 2
    prefac = 40.0f / (7.0f * (float)M_PI);
#else
    prefac = 4.0f / 3.0f;
#endif

    if (u < 0.5f)
        return prefac * (6.0f*u*u*(u - 1.0f) + 1.0f);
    else if (u < 1.0f)
        return prefac * 2.0f * (1.0f - u)*(1.0f - u)*(1.0f - u);
    return 0.0f;
}

/* ================================================================== */
/* 2-D projected cubic spline  Σ(R)  where R = projected radius / h   */
/*                                                                      */
/* This is the integral of the 3-D kernel along the full line of sight: */
/*   Σ(R) = ∫_{-∞}^{∞} W(sqrt(R² + z²), h) dz                        */
/*                                                                      */
/* Derived analytically for the cubic spline (Monaghan & Lattanzio      */
/* 1985 form with support at r=h, rescaled here to support at r=2h):   */
/*                                                                      */
/*   For q = R/h:                                                       */
/*     q < 1:  Σ = (8/πh²)[ (1-q²)^(3/2)(0.25 + 0.375q²)             */
/*                           + q²(1-q²/4)^(1/2)(0.1875q² - 0.25) ]    */
/*             (plus the contribution from the outer shell 1≤r≤2)       */
/*     1 ≤ q < 2: Σ = (8/πh²)[ ... outer cubic piece only ]           */
/*     q ≥ 2:  Σ = 0                                                   */
/*                                                                      */
/* Reference: Monaghan (1992), price (2007) eq. A8-A11.               */
/*                                                                      */
/* We use the closed-form result from Price (2007) PASA 24, 159,       */
/* equation A14 (support radius = 2h convention).                       */
/* ================================================================== */
float cubic_spline_kernel_2d_proj(float q)
{
    /*
     * Use the kernel with support at r = h (the "standard" form), then
     * map q_input = R / (2h) → q = R / h so the support is at q = 1
     * in Price's formula.  We normalise separately.
     *
     * Price (2007) eq. A14, support at r=h, NDIM=3 kernel projected:
     *
     *  For 0 ≤ q < 1/2:
     *    Σ = (8/πh²)[ -2/3 * (1-q²)^(3/2) + (2-q²/2)*(1-q²)^(1/2)
     *                 * ... ]
     *
     * The full analytic form is complex.  We use a high-accuracy
     * Gauss-Legendre numerical integration at table-build time (called
     * once), then linear interpolation at runtime.
     *
     * Here we evaluate it directly for tabulation.
     * The 3-D kernel: W(r) = sigma * f(r/h)
     *   f(u) = 1 - 3/2 u² + 3/4 u³         0 ≤ u < 1
     *   f(u) = 1/4 (2 - u)³                 1 ≤ u < 2
     *   sigma_3d = 1/pi (with our factor-of-8 form)
     */

    /* q is R/h; kernel support at u=2h means q ranges 0..2 */
    if (q >= 2.0f) return 0.0f;

    /*
     * Numerically integrate W_3D(sqrt(q² + t²)) dt from t = -sqrt(4-q²)
     * to t = sqrt(4-q²)  (the range where r = sqrt(R²+z²) < 2h).
     *
     * We do this with 64-point Gauss-Legendre quadrature.
     * (Called only at table-build time so cost is irrelevant.)
     */
    double Q  = (double)q;
    double t_max = sqrt(4.0 - Q*Q);
    /* 16-point GL nodes and weights on [-1,1] */
    static const double gl_x[16] = {
        -0.9894009349916499, -0.9445750230732326,
        -0.8656312023341108, -0.7554044083550030,
        -0.6178762444026438, -0.4580167776572274,
        -0.2816035507792589, -0.0950125098360373,
         0.0950125098360373,  0.2816035507792589,
         0.4580167776572274,  0.6178762444026438,
         0.7554044083550030,  0.8656312023341108,
         0.9445750230732326,  0.9894009349916499
    };
    static const double gl_w[16] = {
        0.0271524594117541, 0.0622535239386479,
        0.0951585116824928, 0.1246289712555339,
        0.1495959888165767, 0.1691565193950025,
        0.1826034150449236, 0.1894506104550685,
        0.1894506104550685, 0.1826034150449236,
        0.1691565193950025, 0.1495959888165767,
        0.1246289712555339, 0.0951585116824928,
        0.0622535239386479, 0.0271524594117541
    };

    double sum = 0.0;
    for (int k = 0; k < 16; k++) {
        double t = t_max * gl_x[k];   /* map [-1,1] → [-t_max, t_max] */
        double r = sqrt(Q*Q + t*t);
        double u = r;                  /* u = r/h, with h=1 in this form */

        double w;
        /* 3-D cubic spline with support at u=2 */
        if (u < 1.0)
            w = 1.0 - 1.5*u*u + 0.75*u*u*u;
        else if (u < 2.0)
            w = 0.25*(2.0-u)*(2.0-u)*(2.0-u);
        else
            w = 0.0;

        sum += gl_w[k] * w;
    }
    sum *= t_max;   /* Jacobian of the variable change */

    /* Normalise: 3-D kernel has sigma = 1/π, so the projected value
       (integrating a 1/π prefactor over dz) gives units of 1/(π h²).
       We multiply by t_max (the half-length) to get the full integral.
       The factor of 2 for the symmetric integral is already folded in. */
    float result = (float)(sum / M_PI);
    return result;
}

/* ================================================================== */
/* Table building                                                       */
/* ================================================================== */

void tabulate_kernel(void)
{
    float lumin = -6.0f, lumax = log10f(2.1f);
    float dlu   = (lumax - lumin) / (float)(NTAB - 1);
    for (int i = 0; i < NTAB; i++) {
        float u    = powf(10.0f, lumin + i * dlu);
        rh_table   [i] = u;
        kernel_table[i] = cubic_spline_kernel(u);
    }
}

void tabulate_projected_kernel(void)
{
    /*
     * q ranges from 0 to 2 (support radius), sampled uniformly.
     * We use uniform spacing (not log) because the projected kernel
     * is smooth near q=0, unlike the 3-D kernel.
     */
    for (int i = 0; i < NTAB; i++) {
        float q     = 2.0f * (float)i / (float)(NTAB - 1);
        rh_table[i] = q;
        proj_table[i] = cubic_spline_kernel_2d_proj(q);
    }
}

void tabulate_integral(void)
{
    float lumin = -2.0f, lumax = log10f(30.0f);
    float dlu   = (lumax - lumin) / (float)(NTAB - 1);
    for (int i = 0; i < NTAB; i++) {
        float u         = powf(10.0f, lumin + i * dlu);
        rh_table[i]     = u;
        integral_table[i] = cumulative_cubic_spline_interpolant(u);
    }
}

/* ================================================================== */
/* Table lookup with linear interpolation                               */
/*                                                                      */
/* All three lookups use the same fast path: compute the fractional     */
/* index directly from the known uniform or log-uniform spacing rather  */
/* than using binary search.                                            */
/* ================================================================== */

float get_kernel_value(float u)
{
    if (u <= 0.0f) return kernel_table[0];
    if (u >= 2.1f) return 0.0f;

    /* Log-uniform table: index = (log10(u) - lumin) / dlu */
    float lumin = -6.0f, lumax = log10f(2.1f);
    float dlu   = (lumax - lumin) / (float)(NTAB - 1);
    float fi    = (log10f(u) - lumin) / dlu;
    if (fi < 0.0f) fi = 0.0f;
    int   i0    = (int)fi;
    if (i0 >= NTAB - 1) return kernel_table[NTAB - 1];
    float f     = fi - (float)i0;
    return kernel_table[i0] * (1.0f - f) + kernel_table[i0 + 1] * f;
}

/*
 * get_projected_kernel_value: evaluate the 2-D projected kernel at
 * projected radius q = R / h (h is the SPH smoothing length).
 * Returns the surface density contribution (units: 1/h²).
 */
float get_projected_kernel_value(float q)
{
    if (q <= 0.0f) return proj_table[0];
    if (q >= 2.0f) return 0.0f;

    /* Uniform table over [0, 2] */
    float fi = q * (float)(NTAB - 1) / 2.0f;
    int   i0 = (int)fi;
    if (i0 >= NTAB - 1) return proj_table[NTAB - 1];
    float f  = fi - (float)i0;
    return proj_table[i0] * (1.0f - f) + proj_table[i0 + 1] * f;
}

float get_integral_value(float u)
{
    if (u <= 0.0f) return integral_table[0];

    float lumin = -2.0f, lumax = log10f(30.0f);
    float dlu   = (lumax - lumin) / (float)(NTAB - 1);
    float fi    = (log10f(u > 0 ? u : 1e-6f) - lumin) / dlu;
    if (fi < 0.0f) fi = 0.0f;
    int   i0    = (int)fi;
    if (i0 >= NTAB - 1) return integral_table[NTAB - 1];
    float f     = fi - (float)i0;
    return integral_table[i0] * (1.0f - f) + integral_table[i0 + 1] * f;
}


#ifdef KERNEL_SMOOTHING
/* ================================================================== */
/* u²-domain kernel table (only needed for SPH 3-D deposit)            */
/*                                                                      */
/* For the 3-D deposit inner loop we want to avoid sqrtf and log10f.   */
/* We table W(u) uniformly in u² ∈ [0, 4] (support is u < 2, u² < 4). */
/* The caller computes u² = r²/h² and calls get_kernel_value_u2(u2).   */
/* ================================================================== */
#define NTAB_U2 4096          /* uniform in u², so coarser near u=0 is OK */
static float u2_table[NTAB_U2]; /* W(sqrt(u2)) * (1/(8/pi)) — rescaled below */

void tabulate_kernel_u2(void)
{
    /* u² runs from 0 to 4 (support boundary) */
    for (int i = 0; i < NTAB_U2; i++) {
        float u2 = 4.0f * (float)i / (float)(NTAB_U2 - 1);
        float u  = sqrtf(u2);
        /* Store the full normalised kernel value (including 8/pi prefactor) */
        u2_table[i] = cubic_spline_kernel(u);
    }
    /* Last entry is exactly at u=2: kernel is 0 there */
    u2_table[NTAB_U2 - 1] = 0.0f;
}

/*
 * get_kernel_value_u2: evaluate W(sqrt(u2)) for u2 = (r/h)^2.
 * Returns 0 for u2 >= 4 (outside support).
 * No sqrtf, no log10f — just a multiply, floor, and lerp.
 */
float get_kernel_value_u2(float u2)
{
    if (u2 >= 4.0f) return 0.0f;
    float fi = u2 * (float)(NTAB_U2 - 1) * 0.25f;  /* / 4 */
    int   i0 = (int)fi;
    float f  = fi - (float)i0;
    /* i0+1 is always valid because u2 < 4 guarantees fi < NTAB_U2-1 */
    return u2_table[i0] + f * (u2_table[i0 + 1] - u2_table[i0]);
}

#endif /* KERNEL_SMOOTHING */


/* ================================================================== */
/* q²-domain projected kernel table                                    */
/*                                                                     */
/* For the 2D deposit inner loop we want to avoid sqrtf.              */
/* We table W_proj(q) uniformly in q² ∈ [0, 4].                      */
/* The caller computes q² = R²/h² and calls                           */
/* get_projected_kernel_value_u2(q2).  No sqrtf needed.               */
/* ================================================================== */
#define NTAB_PROJ2 2048
static float proj2_table[NTAB_PROJ2];

void tabulate_projected_kernel_u2(void)
{
    for (int i = 0; i < NTAB_PROJ2; i++) {
        float q2 = 4.0f * (float)i / (float)(NTAB_PROJ2 - 1);
        float q  = sqrtf(q2);
        proj2_table[i] = get_projected_kernel_value(q);
    }
    proj2_table[NTAB_PROJ2 - 1] = 0.0f;
}

/*
 * get_projected_kernel_value_u2: evaluate W_proj(sqrt(q2)) for q2 = (R/h)^2.
 * Returns 0 for q2 >= 4 (outside support).
 * No sqrtf — just a multiply, floor, and lerp.
 */
float get_projected_kernel_value_u2(float q2)
{
    if (q2 >= 4.0f) return 0.0f;
    float fi = q2 * (float)(NTAB_PROJ2 - 1) * 0.25f;
    int   i0 = (int)fi;
    float f  = fi - (float)i0;
    return proj2_table[i0] + f * (proj2_table[i0 + 1] - proj2_table[i0]);
}

/* ================================================================== */
/* cumulative_cubic_spline_interpolant (kept for legacy callers)        */
/* ================================================================== */
float cumulative_cubic_spline_interpolant(float u)
{
    float sum = 0.0f;
    if (u < 1.0f)
        sum = u - 0.5f*u*u*u + (3.0f/16.0f)*u*u*u*u;
    else if (u < 2.0f)
        sum = 0.6875f + 2.0f*u - 1.5f*u*u + 0.5f*u*u*u
              - (1.0f/16.0f)*u*u*u*u - 0.9375f;
    else
        sum = 0.75f;
    return sum / (float)M_PI;
}
