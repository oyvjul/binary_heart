#ifndef BINARY_CUBE_SEQUENTIAL_H
#define BINARY_CUBE_SEQUENTIAL_H

#include <stdio.h>
#include <stdlib.h>

double ***dallocate_3d(int x, int y, int z);
void dinit_3d(double*** matrix, int x, int y, int z);
void init_cubedata(cube *c, int x, int y, int z);
void read_cubemesh(cube *c, char *cubefile);
void init_tensor(cube c, int x, int y, int z);

#endif
