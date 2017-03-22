#include <stdio.h>
#include <stdlib.h>
#include "tensor.h"
#include "binary_cube_sequential.h"
#include "read_mesh.h"
#include "io.h"
#include "diffusion.h"
//#include "omp.h"

void init_sequential_data(cube *c, int x, int y, int z)
{
  c->tensor_x0 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_x0, x+3, y+3, z+3);

  c->tensor_x1 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_x1, x+3, y+3, z+3);

  c->tensor_y0 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_y0, x+3, y+3, z+3);

  c->tensor_y1 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_y1, x+3, y+3, z+3);

  c->tensor_z0 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_z0, x+3, y+3, z+3);

  c->tensor_z1 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_z1, x+3, y+3, z+3);

  c->x = x;
  c->y = y;
  c->z = z;
}

void free_sequential_data(double ***tensor)
{
  free(tensor[0][0]);
  free(tensor[0]);
  free(tensor);
}

void init_data(cube *c, int x, int y, int z)
{
  c->u_old = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->u_old, x+3, y+3, z+3);

  c->u_new = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->u_new, x+3, y+3, z+3);

  c->tensor_val_11 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_11, x+3, y+3, z+3);

  c->tensor_val_12 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_12, x+3, y+3, z+3);

  c->tensor_val_13 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_13, x+3, y+3, z+3);

  c->tensor_val_22 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_22, x+3, y+3, z+3);

  c->tensor_val_23 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_23, x+3, y+3, z+3);

  c->tensor_val_33 = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_val_33, x+3, y+3, z+3);

  /*c->grid_x = (double*)calloc(x+2, sizeof(double));
  c->grid_y = (double*)calloc(y+2, sizeof(double));
  c->grid_z = (double*)calloc(z+2, sizeof(double));*/
}

void decompose(int n, int dim, int coord, int* start, int* end)
{
  int length, rest;

  length = n/dim;
  rest = n%dim;
  *start = coord * length + (coord < rest ? coord : rest);
  *end = *start + length - (coord < rest ? 0 : 1);
}

int main(int argc, char *argv[])
{
  int x = 10;
  int y = 10;
  int z = 10;
  double count = 1.0;
  int i, j, k;
  int ii, jj, kk;
  int start_x, start_y, start_z;
  int size, rank;
  int coords[3];
  int periods[3];
  int dims[3];
  int nd;
  int procs_x, procs_y, procs_z;
  int left, right, up, down, z_up, z_down;
  int x0, x1, y0, y1, z0, z1;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  cube *c = (cube*)malloc(sizeof(cube));
  meshdata *m = (meshdata*)malloc(sizeof(meshdata));
  tensorfield *t = (tensorfield *)malloc(sizeof(tensorfield));

  init_sequential_data(c, x, y, z);
  init_cube_grid_mpi(c, m);

  /*read_binaryformat("mesh_new/tensor_x0.tensor", &c->tensor_x0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_x1.tensor", &c->tensor_x1, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_y0.tensor", &c->tensor_y0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_y1.tensor", &c->tensor_y1, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_z0.tensor", &c->tensor_z0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_z1.tensor", &c->tensor_z1, x+3, y+3, z+3);*/


  for(i = 1; i <= c->z+1; i++)
  {
    for(j = 1; j <= c->y+1; j++)
    {
      for(k = 1; k <= c->x+1; k++)
      {
        c->tensor_x0[i][j][k] = count;
        count++;
      }
    }
  }

  //init_data(c, x, y, z);

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
  }


  //read_cubemesh(c, "mesh_new/128x128x128.txt");
  //sparse_readtensorfiles("mesh_new/3Dheart.1", t, 1000);
  //fiberstotensors(t);
  //simple_averagetensors(c, t);
  #ifdef GENERATE_AND_DUMP_TENSORS
    count_inside_tensor = 0;
    count_outside_tensor = 0;
    count_inside = 0;
    double start = omp_get_wtime();
    sparse_readtensorfiles("mesh_new/3Dheart.1", t, 1000);
    fiberstotensors(t);
    generate_tensor(c, t, m);
    printf("GENERATE_TENSORS_ONLY\n");
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
  #endif

  return 0;
}
