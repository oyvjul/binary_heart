#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "binary_cube_sequential.h"

const double SIGMA_l =  1.21321;
const double SIGMA_t =  0.2121;

typedef struct {
    int numtensor;
    double* coord;
    double* inputtensor;
    double* fibers;

} tensorfield;

void sparse_readtensorfiles(char* tensorfile,tensorfield* T,int skip);
void fiberstotensors(tensorfield* T);
void simple_averagetensors(cube c,tensorfield* T);
