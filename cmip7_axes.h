/*
 * cmip7_axes.h
 */
#ifndef CMIP7_AXIS_H
#define CMIP7_AXIS_H

#include "gtool3.h"

typedef GT3_Dim* (*axis_gen_func_t)(void);

/* hybrid.c (HETA90) */
GT3_Dim *alevs_heta90_half(void);
GT3_Dim *alevs_heta90(void);
void make_heta90(void);

/* hybrid.c (HETA40) */
GT3_Dim *alevs_heta40_half(void);
GT3_Dim *alevs_heta40(void);
void make_heta40();

/* ocean.c */
void make_med52(void);

/* plevs.c */
void make_plevs(void);

/* utils.c */
GT3_Dim *new_GT3_Dim(size_t size,
                     const char *name,
                     const char *title,
                     const char *units,
                     const double *values);
int write_gtaxloc(const GT3_Dim *dim);
int write_axes(GT3_Dim **axes, size_t naxes);

#endif /* !CMIP7_AXIS_H */
