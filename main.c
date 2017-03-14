#include <stdio.h>
#include "tensor.h"
#include "binary_cube_sequential.h"
#include "read_mesh.h"

int main(int argc, char *argv[])
{
  int x = 128;
  int y = 128;
  int z = 128;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  int i, j, k;
  cube *c = (cube*)malloc(sizeof(cube));
  meshdata *m = (meshdata*)malloc(sizeof(meshdata));
  tensorfield *t = (tensorfield *)malloc(sizeof(tensorfield));

  /*init_meshdata(&c, x, y, z);
  read_cubemesh(c.u_old, x, y, z);
  sparse_readtensorfiles("mesh_new/3Dheart.1", T, 1000);
  fiberstotensors(T);
  simple_averagetensors(c, T, x, y, z);*/
  return 0;
}
