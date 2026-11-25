/*
 * utils.c
 */
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmip7_axes.h"


GT3_Dim *
new_GT3_Dim(size_t size,
            const char *name,
            const char *title,
            const char *units,
            const double *values)
{
    GT3_Dim *dim = NULL;
    size_t i;

    if ((dim = malloc(sizeof(GT3_Dim))) == NULL) {
        perror(NULL);
        return NULL;
    }
    dim->name = NULL;
    dim->values = NULL;
    dim->title = NULL;
    dim->unit = NULL;

    if ((dim->name = strdup(name)) == NULL
        || (dim->title = strdup(title)) == NULL
        || (dim->unit = strdup(units)) == NULL
        || (dim->values = malloc(size * sizeof(double))) == NULL) {
        perror(NULL);
        GT3_freeDim(dim);
        return NULL;
    }

    dim->len = size;
    dim->range[0] = -999.0;
    dim->range[1] = -999.0;
    dim->cyclic = 0;
    for (i = 0; i < size; ++i)
        dim->values[i] = values[i];

    return dim;
}


int
write_gtaxloc(const GT3_Dim *dim)
{
    char path[PATH_MAX + 1];
    FILE *fp;
    int rc = 0;

    snprintf(path, sizeof path, "GTAXLOC.%s", dim->name);
    if ((fp = fopen(path, "wb")) == NULL) {
        perror(path);
        return -1;
    }
    if (GT3_writeDimFile(fp, dim, "UR8") < 0) {
        GT3_printErrorMessages(stderr);
        rc = -1;
    }

    fclose(fp);
    return rc;
}


int
write_axes(GT3_Dim **axes, size_t naxes)
{
    int count = 0;
    size_t i;

    for (i = 0; i < naxes; ++i) {
        if (axes[i]) {
            printf("Generating %s ... ", axes[i]->name);
            if (write_gtaxloc(axes[i]) < 0)
                printf("failed\n");
            else {
                count++;
                printf("done\n");
            }
        }
    }
    return count;
}
