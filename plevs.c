/*
 * plevs.c
 */
#include <assert.h>
#include "cmip7_axes.h"


static GT3_Dim *
plevs_cmip6p38(void)
{
    double plevs[] = {
        100000.0, 97500.0,  95000.0,  92500.0,  90000.0,
        87500.0,  85000.0,  82500.0,  80000.0,  77500.0,
        75000.0,  70000.0,  65000.0,  60000.0,  55000.0,
        50000.0,  45000.0,  40000.0,  35000.0,  30000.0,
        25000.0,  22500.0,  20000.0,  17500.0,  15000.0,
        12500.0,  10000.0,   7000.0,   5000.0,   3000.0,
        2000.0,   1000.0,     700.0,    500.0,    300.0,
        200.0,    100.0,
        40.0,
    };
    assert(sizeof plevs / sizeof(double) == 38);

    return new_GT3_Dim(38, "CMIP6P38", "Pressure", "Pa", plevs);
}


static GT3_Dim *
plevs_ecmanlp37(void)
{
    double plevs[] = {
        1000.0, 975.0,  950.0,  925.0,  900.0,
        875.0,  850.0,  825.0,  800.0,  775.0,
        750.0,  700.0,  650.0,  600.0,  550.0,
        500.0,  450.0,  400.0,  350.0,  300.0,
        250.0,  225.0,  200.0,  175.0,  150.0,
        125.0,  100.0,   70.0,   50.0,   30.0,
        20.0,   10.0,     7.0,    5.0,    3.0,
        2.0,    1.0,
    };
    assert(sizeof plevs / sizeof(double) == 37);

    return new_GT3_Dim(38, "ECMANLP37", "pressure", "hPa", plevs);
}



static GT3_Dim *
plevs_ecmanlp23(void)
{
    double plevs[] = {
        1000.0,  925.0,  850.0,  775.0,  700.0,
        600.0,   500.0,  400.0,  300.0,  250.0,
        200.0,   150.0,  100.0,   70.0,   50.0,
        30.0,     20.0,   10.0,    7.0,    5.0,
        3.0,       2.0,    1.0,
    };
    assert(sizeof plevs / sizeof(double) == 23);

    return new_GT3_Dim(23, "ECMANLP23", "Pressure", "hPa", plevs);
}


static GT3_Dim *
plevs_cmip6p7h(void)
{
    double plevs[] = {
        925.0,  850.0,  700.0,  600.0,  500.0,
        250.0,   50.0,
    };
    assert(sizeof plevs / sizeof(double) == 7);

    return new_GT3_Dim(7, "CMIP6P7h", "Pressure", "hPa", plevs);
}


static GT3_Dim *
plevs_cmip6p39(void)
{
    double plevs[] = {
        100000.0, 92500.0,  85000.0,  70000.0,  60000.0,
        50000.0,  40000.0,  30000.0,  25000.0,  20000.0,
        17000.0,  15000.0,  13000.0,  11500.0,  10000.0,
        9000.0,    8000.0,   7000.0,   5000.0,   3000.0,
        2000.0,    1500.0,   1000.0,    700.0,    500.0,
        300.0,      200.0,    150.0,    100.0,     70.0,
        50.0,        40.0,     30.0,     20.0,     15.0,
        10.0,         7.0,      5.0,      3.0,
    };
    assert(sizeof plevs / sizeof(double) == 39);

    return new_GT3_Dim(39, "CMIP6P39", "Pressure", "Pa", plevs);
}


static GT3_Dim *
plevs_cmip6p7c(void)
{
    double plevs[] = {
        1000.0,  850.0,  700.0,  500.0,  300.0,
        200.0,   100.0,
    };
    assert(sizeof plevs / sizeof(double) == 7);

    return new_GT3_Dim(7, "CMIP6P7c", "Pressure", "hPa", plevs);
}


void
make_plevs()
{
    axis_gen_func_t gen_axes[] = {
        plevs_cmip6p38,
        plevs_ecmanlp37,
        plevs_ecmanlp23,
        plevs_cmip6p7h,
        plevs_cmip6p39,
        plevs_cmip6p7c,
    };
    size_t nfuncs = sizeof gen_axes / sizeof(axis_gen_func_t);
    size_t i;

    for (i = 0; i < nfuncs; ++i) {
        GT3_Dim *dim = gen_axes[i]();

        if (dim) {
            printf("Generating %s ... ", dim->name);
            if (write_gtaxloc(dim) < 0)
                printf("failed\n");
            else
                printf("done\n");

            GT3_freeDim(dim);
        }
    }
}
