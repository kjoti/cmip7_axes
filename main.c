/*
 * main.c
 */
#include <stdio.h>
#include <stdlib.h>
#include "cmip7_axes.h"


int
main(int argc, char **argv)
{
    make_plevs();
    make_heta90();
    make_heta40();
    make_med52();
    return 0;
}
