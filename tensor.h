#ifndef TENSOR_H
#define TENSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "binary_cube_sequential.h"

typedef struct
{
    int numtensor;
    double* coord;
    double* inputtensor;
    double* fibers;

} tensorfield;

void sparse_readtensorfiles(char* tensorfile,tensorfield* T,int skip);
void fiberstotensors(tensorfield* T);
void simple_averagetensors(cube *c,tensorfield* T);

#endif
