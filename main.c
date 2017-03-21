#include <stdio.h>
#include "tensor.h"
#include "binary_cube_sequential.h"
#include "read_mesh.h"
#include "io.h"
#include "diffusion.h"
#include "omp.h"

int main(int argc, char *argv[])
{
  int x = 6;
  int y = 6;
  int z = 6;
  double count = 1.0;
  int i, j, k;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  cube *c = (cube*)malloc(sizeof(cube));
  meshdata *m = (meshdata*)malloc(sizeof(meshdata));
  tensorfield *t = (tensorfield *)malloc(sizeof(tensorfield));

  init_cubedata(c, x, y, z);
  init_cube_grid(c, m);



  //read_cubemesh(c, "mesh_new/128x128x128.txt");
  //sparse_readtensorfiles("mesh_new/3Dheart.1", t, 1000);
  //fiberstotensors(t);
  //simple_averagetensors(c, t);
  #ifdef GENERATE_AND_DUMP_TENSORS
    count_inside_tensor = 0;
    count_outside_tensor = 0;
    count_inside = 0;
    for(i = 1; i <= x+1; i++)
    {
      printf("%f ", c->grid_x[i]);
    }
    double start = omp_get_wtime();
    sparse_readtensorfiles("mesh_new/3Dheart.1", t, 1000);
    fiberstotensors(t);
    generate_tensor(c, t, m);
    printf("GENERATED GENERATE_TENSORS_ONLY\n");
    double end = omp_get_wtime();
    printf("it took : %0.12f \n", end-start);

    for(i = 1; i <= c->z+1; i++)
    {
      for(j = 1; j <= c->y+1; j++)
      {
        for(k = 1; k <= c->x+1; k++)
        {
          if(c->tensor_x0[i][j][k] == 0)
          {
            count_outside_tensor++;
          }
          else
          {
            count_inside_tensor++;
          }

          count_inside++;
        }
      }
    }

    double l2_tensor = 0;

    for(i = 1; i <= c->z+1; i++)
    {
      for(j = 1; j <= c->y+1; j++)
      {
        for(k = 1; k <= c->x+1; k++)
        {
          l2_tensor += c->tensor_x0[i][j][k]*c->tensor_x0[i][j][k];
        }
      }
    }

    printf("TENSOR: inside points: %d \t outside points: %d sum check: %d  total: %d\n", count_inside_tensor, count_outside_tensor, c->x*c->y*c->z, count_inside);

    /*for(i = 1; i <= z+1; i++)
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

    /*write_binaryformat("mesh_new/tensor_x0.tensor", c->tensor_x0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_x1.tensor", c->tensor_x1, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_y0.tensor", c->tensor_y0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_y1.tensor", c->tensor_y1, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_z0.tensor", c->tensor_z0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_z1.tensor", c->tensor_z1, x+3, y+3, z+3);*/
    //read_binaryformat("mesh_new/test.tensor", &c->tensor_x0, &x+3, &y+3, &z+3);
  #else

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

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        divergence_cell_direction_z(c->u_old, count, count, count, c->x_step, c->y_step, c->y_step, i, j, k);
      }
    }
  }*/
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
  #endif

  return 0;
}
