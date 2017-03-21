#include <stdio.h>
#include "tensor.h"
#include "binary_cube_sequential.h"
#include "read_mesh.h"
#include "io.h"
#include "diffusion.h"

int main(int argc, char *argv[])
{
  int x = 8;
  int y = 8;
  int z = 8;
  double count = 0.0;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  int i, j, k;
  cube *c = (cube*)malloc(sizeof(cube));
  meshdata *m = (meshdata*)malloc(sizeof(meshdata));
  tensorfield *t = (tensorfield *)malloc(sizeof(tensorfield));

  init_cubedata(c, x, y, z);
  init_cube_grid(c, m);
  read_cubemesh(c, "mesh_new/128x128x128.txt");
  //sparse_readtensorfiles("mesh_new/3Dheart.1", t, 1000);
  //fiberstotensors(t);
  //simple_averagetensors(c, t);

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        c->u_old[i][j][k] = count;
        count++;
      }
    }
  }

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        divergence_direction_x(c->u_old, 2.1, 2.1, 2.1, c->x_step, c->y_step, c->y_step, i, j, k);
      }
    }
  }

  /*write_binaryformat("mesh_new/test.tensor", c->u_old, x+3, y+3, z+3);
  read_binaryformat("mesh_new/test.tensor", &c->tensor_x0, &x+3, &y+3, &z+3);

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        printf("%f ", c->tensor_x0[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }*/

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        if(c->u_old[i][j][k] == 1)
        {
          count_inside++;
        }
        else
        {
          count_outside++;
        }
      }
    }
  }

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        if(c->tensor_x0[i][j][k] == 0)
        {
          count_outside_tensor++;
        }
        else
        {
          count_inside_tensor++;
        }
      }
    }
  }

  double l2_uold = 0;

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        l2_uold += c->u_old[i][j][k]*c->u_old[i][j][k];
      }
    }
  }

  double l2_tensor = 0;

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        l2_tensor += c->tensor_x0[i][j][k]*c->tensor_x0[i][j][k];
      }
    }
  }

  printf("DEBUG: \n");
  printf("U_OLD: inside points: %d \t outside points: %d sum check: %0.12f \n", count_inside, count_outside, l2_uold);
  printf("TENSOR: inside points: %d \t outside points: %d sum check: %0.12f \n", count_inside_tensor, count_outside_tensor, l2_tensor);*/

  return 0;
}
