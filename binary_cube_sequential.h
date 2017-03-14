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

  double x_step;
  double y_step;
  double z_step;

  int x;
  int y;
  int z;
}cube;

double ***dallocate_3d(int x, int y, int z);
void dinit_3d(double*** matrix, int x, int y, int z);
void init_cubedata(cube *c, int x, int y, int z);
void read_cubemesh(cube *c, char *cubefile);
void init_tensor(cube c, int x, int y, int z);

#endif
