#ifndef BINARY_CUBE_SEQUENTIAL_H
#define BINARY_CUBE_SEQUENTIAL_H

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
  double ***u_old;
  double ***u_new;
  double ***tensor_x0;
  double ***tensor_x1;
  double ***tensor_y0;
  double ***tensor_y1;
  double ***tensor_z0;
  double ***tensor_z1;
  double *grid_x;
  double *grid_y;
  double *grid_z;
}cube;

double ***dallocate_3d(int x, int y, int z);
void dinit_3d(double*** matrix, int x, int y, int z);
void init_meshdata(cube *c, int x, int y, int z);
void read_meshdata(double ***u_old, int x, int y, int z);
void init_tensor(cube c, int x, int y, int z);

#endif
