#ifndef READ_MESH_H
#define READ_MESH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "binary_cube_sequential.h"

typedef struct
{
    int* elements;
    int* neighbours;
    double* Gx;
    double* Gy;
    double* Gz;

    double* nodes;
    double* centroid;
    double* area;
    double* Nx;
    double* Ny;
    double* Nz;
    double* Sx;
    double* Sy;
    double* Sz;

    double* tensor;
    double* volume;
    double totalVolume;
    int numtet;
    int numnodes;
} meshdata;

/*double ***dallocate_3d(int x, int y, int z);
void dinit_3d(double*** matrix, int x, int y, int z);*/
void readmesh(char* infile, meshdata *M);
void compute_minmax(double *x_max, double *x_min, double *y_max, double *y_min, double *z_max, double *z_min, meshdata *m);
void calculate_centroid(meshdata *m);
double determinant(double *a, double *b, double *c, double *d);
int inside(int numtet, int *elements, double *nodes, double point_x, double point_y, double point_z);
void init_cube_grid(cube *c, meshdata *m, int x, int y, int z);

#endif
