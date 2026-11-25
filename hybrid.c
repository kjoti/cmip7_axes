/*
 * hybrid.c
 *
 * standard_hybrid_sigma: `p = a * p0 + b * ps`
 * reference: miroc/util/hybrid/mkzlhyb.F
 *
 * Six GTAXLOC.*:
 *   HETA[0-9]*:      eta (full level)
 *   HETA[0-9]*_A     FA  (full level)
 *   HETA[0-9]*_B     FB  (full level)
 *   HETA[0-9]*.M     eta (half level)
 *   HETA[0-9]*_A.M   FA  (half level)
 *   HETA[0-9]*_B.M   FB  (half level)
 */
#include <assert.h>
#include <math.h>
#include "cmip7_axes.h"

#define MAX_PLEVS 1000

const double CP = 1004.60;   /* specific heat of air */
const double RAIR = 287.04;  /* gas constant of air */
const double PREF0 = 1e5;    /* reference pressure [Pa] */
const double PINT = 35000.0; /* interface pressure for pure-P coord. [Pa] */


/*
 * `full`: size len
 * `half`: size len + 1
 */
static void
calc_full_from_half(double *full, const double *half, size_t len)
{
    double ak1 = 1.0 + RAIR / CP;
    double akiv = CP / RAIR;
    double z[MAX_PLEVS + 1];
    double dz, ds;
    size_t k;

    assert(len <= MAX_PLEVS);

    for (k = 0; k < len + 1; ++k) {
        z[k] = pow(half[k], ak1);
    }
    for (k = 0; k < len; ++k) {
        dz = z[k] - z[k + 1];
        ds = half[k] - half[k + 1];
        full[k] = pow(dz / (ak1 * ds), akiv);
    }

    assert(full[0] <= 1.0);
    for (k = 1; k < len; ++k) {
        /* Check descending order. */
        assert(full[k] < full[k - 1]);

        /* half[k] >= full[k] >= half[k+1] */
        assert(full[k] <= half[k]);
        assert(full[k] >= half[k + 1]);
    }
}


/*
 * plevs(k) = fa(k) + fb(k) * PREF0
 *
 * fa: [Pa]
 * fb: [1]
 */
static void
hybrid_sigma_coef(double *fa, double *fb,
                  const double *plevs, size_t len, double p_intr)
{
    double sigma;
    size_t k;

    assert(p_intr >= 0.);
    assert(p_intr < PREF0);
    for (k = 0; k < len; ++k) {
        sigma = (plevs[k] - p_intr) / (PREF0 - p_intr);
        if (sigma < 0.)
            sigma = 0.;

        fa[k] = plevs[k] - sigma * PREF0;
        fb[k] = sigma;
    }
}


static double
pressure_inter(const double *plevs, size_t len)
{
#if 0
    return PINT;
#else
    size_t k;

    for (k = 0; k < len; ++k) {
        if (plevs[k] <= PINT)
            return plevs[k];
    }
    return 0.;
#endif
}


/*
 * `HETA90.M` (half levels
 * (cf. include/zll90h.F)
 */
GT3_Dim *
alevs_heta90_half()
{
    double sigma[] = {
        1.0000     , 0.99500    , 0.98800    , 0.97800    , 0.96600    ,
        0.95100    , 0.93300    , 0.91200    , 0.88900    , 0.86300    ,
        0.83500    , 0.80500    , 0.76900    , 0.72500    , 0.67300    ,
        0.61300    , 0.55100    , 0.48800    , 0.43234    , 0.38328    ,
        0.34010    , 0.30217    , 0.26888    , 0.23969    , 0.21412    ,
        0.19171    , 0.17208    , 0.15488    , 0.13980    , 0.12657    ,
        0.11495    , 0.10474    , 0.95750E-01, 0.87826E-01, 0.80833E-01,
        0.74649E-01, 0.69122E-01, 0.64004E-01, 0.59266E-01, 0.54878E-01,
        0.50815E-01, 0.47052E-01, 0.43569E-01, 0.40339E-01, 0.37342E-01,
        0.34559E-01, 0.31971E-01, 0.29563E-01, 0.27321E-01, 0.25233E-01,
        0.23286E-01, 0.21471E-01, 0.19778E-01, 0.18198E-01, 0.16724E-01,
        0.15349E-01, 0.14066E-01, 0.12869E-01, 0.11753E-01, 0.10712E-01,
        0.97436E-02, 0.88418E-02, 0.80032E-02, 0.72244E-02, 0.65018E-02,
        0.58326E-02, 0.52137E-02, 0.46427E-02, 0.41168E-02, 0.36338E-02,
        0.31915E-02, 0.27878E-02, 0.24205E-02, 0.20879E-02, 0.17880E-02,
        0.15190E-02, 0.12791E-02, 0.10666E-02, 0.87986E-03, 0.71704E-03,
        0.57647E-03, 0.45645E-03, 0.35525E-03, 0.26415E-03, 0.18611E-03,
        0.12303E-03, 0.75427E-04, 0.42304E-04, 0.21356E-04, 0.95203E-05,
        0.0000,
    };
    size_t k;

    assert(sizeof sigma / sizeof(double) == 91);

    /* Check descending order. */
    assert(sigma[0] == 1.0);
    for (k = 1; k < 91; ++k) {
        assert(sigma[k] < sigma[k - 1]);
    }

    return new_GT3_Dim(91, "HETA90.M", "eta", "1", sigma);
}


/*
 * `HETA90` (full levels)
 */
GT3_Dim *
alevs_heta90()
{
    GT3_Dim *bnds = alevs_heta90_half();
    GT3_Dim *lev;
    double full[90];

    calc_full_from_half(full, bnds->values, 90);

    lev = new_GT3_Dim(90, "HETA90", "eta", "1", full);
    GT3_freeDim(bnds);
    return lev;
}


/*
 * `HETA40` in half levels (cf. include/zll40h.F)
 */
GT3_Dim *
alevs_heta40_half()
{
    double sigma[] = {
        1.0000     , 0.99500    , 0.98800    , 0.97800    , 0.96600    ,
        0.95100    , 0.93300    , 0.91200    , 0.88900    , 0.86300    ,
        0.83500    , 0.80500    , 0.76900    , 0.72500    , 0.67300    ,
        0.61300    , 0.53655    , 0.46698    , 0.40642    , 0.35372    ,
        0.30785    , 0.26793    , 0.23318    , 0.20295    , 0.17663    ,
        0.15372    , 0.13379    , 0.11644    , 0.10134    , 0.88200E-01,
        0.76763E-01, 0.66809E-01, 0.58145E-01, 0.50475E-01, 0.43462E-01,
        0.36917E-01, 0.30384E-01, 0.22927E-01, 0.14526E-01, 0.70000E-02,
        0.0000
    };
    int k;

    assert(sizeof sigma / sizeof(double) == 41);

    /* Check descending order. */
    assert(sigma[0] == 1.0);
    for (k = 1; k < 41; ++k) {
        assert(sigma[k] < sigma[k - 1]);
    }

    return new_GT3_Dim(41, "HETA40.M", "eta", "1", sigma);
}


/*
 * `HETA40` in full levels.
 */
GT3_Dim *
alevs_heta40()
{
    GT3_Dim *bnds = alevs_heta40_half();
    GT3_Dim *lev;
    double full[40];

    calc_full_from_half(full, bnds->values, 40);

    lev = new_GT3_Dim(40, "HETA40", "eta", "1", full);
    GT3_freeDim(bnds);
    return lev;
}


/*
 * Generate six GTAXLOC files (HETA90*).
 */
void
make_heta90()
{
    const size_t NLEVS = 90;
    double plev_half[91], plev[91];
    double a_half[91], b_half[91];
    double a[90], b[90];
    double p_intr;
    GT3_Dim *axes[6];
    size_t i;

    axes[0] = alevs_heta90_half(); /* HEATA90.m */
    axes[1] = alevs_heta90();      /* HEATA90 */

    for (i = 0; i < NLEVS + 1; ++i)
        plev_half[i] = axes[0]->values[i] * PREF0;
    for (i = 0; i < NLEVS; ++i)
        plev[i] = axes[1]->values[i] * PREF0;

    /* interface pressure [Pa] */
    p_intr = pressure_inter(plev_half, NLEVS + 1);
    printf("HETA90: PINT = %10.2f Pa\n", p_intr);

    /* half levels */
    hybrid_sigma_coef(a_half, b_half, plev_half, NLEVS + 1, p_intr);

    /* full levels */
    hybrid_sigma_coef(a, b, plev, NLEVS, p_intr);

    axes[2] = new_GT3_Dim(NLEVS,     "HETA90_A",   "FA", "Pa", a);
    axes[3] = new_GT3_Dim(NLEVS,     "HETA90_B",   "FB", "1",  b);
    axes[4] = new_GT3_Dim(NLEVS + 1, "HETA90_A.M", "FA", "Pa", a_half);
    axes[5] = new_GT3_Dim(NLEVS + 1, "HETA90_B.M", "FB", "1",  b_half);

    write_axes(axes, 6);
    for (i = 0; i < 6; ++i) {
        GT3_freeDim(axes[i]);
    }
}


/*
 * Generate six GTAXLOC files (HETA40*).
 */
void
make_heta40()
{
    const size_t NLEVS = 40;
    double plev_half[41], plev[41];
    double a_half[41], b_half[41];
    double a[40], b[40];
    double p_intr;
    GT3_Dim *axes[6];
    size_t i;

    axes[0] = alevs_heta40_half(); /* HEATA40.m */
    axes[1] = alevs_heta40();      /* HEATA40 */

    for (i = 0; i < NLEVS + 1; ++i)
        plev_half[i] = axes[0]->values[i] * PREF0;
    for (i = 0; i < NLEVS; ++i)
        plev[i] = axes[1]->values[i] * PREF0;

    /* interface pressure [Pa] */
    p_intr = pressure_inter(plev_half, NLEVS + 1);
    printf("HETA40: PINT = %10.2f Pa\n", p_intr);

    /* half levels */
    hybrid_sigma_coef(a_half, b_half, plev_half, NLEVS + 1, p_intr);

    /* full levels */
    hybrid_sigma_coef(a, b, plev, NLEVS, p_intr);

    axes[2] = new_GT3_Dim(NLEVS,     "HETA40_A",   "FA", "Pa", a);
    axes[3] = new_GT3_Dim(NLEVS,     "HETA40_B",   "FB", "1",  b);
    axes[4] = new_GT3_Dim(NLEVS + 1, "HETA40_A.M", "FA", "Pa", a_half);
    axes[5] = new_GT3_Dim(NLEVS + 1, "HETA40_B.M", "FB", "1",  b_half);

    write_axes(axes, 6);
    for (i = 0; i < 6; ++i) {
        GT3_freeDim(axes[i]);
    }
}
