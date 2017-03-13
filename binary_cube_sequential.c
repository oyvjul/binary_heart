#include <stdio.h>
#include <stdlib.h>
#include "binary_cube_sequential.h"
#include "tensor.h"

double ***dallocate_3d(int x, int y, int z)
{
  int i, j;
  double *storage = (double*)malloc(x * y * z * sizeof(*storage));
  double *alloc = storage;
  double ***matrix;
  matrix = (double***)malloc(z * sizeof(double**));

  for (i = 0; i < z; i++)
  {
    matrix[i] = (double**)malloc(y * sizeof(**matrix));

    for (j = 0; j < y; j++)
    {
      matrix[i][j] = alloc;
      alloc += x;
    }
  }

  return matrix;
}

void dinit_3d(double*** matrix, int x, int y, int z)
{
  int i, j, k;

  for(i = 0; i < z; i++)
  {
    for(j = 0; j < y; j++)
    {
      for(k = 0; k < x; k++)
      {
        matrix[i][j][k] = 0.0;
      }
    }
  }
}

void init_meshdata(cube *c, int x, int y, int z)
{
  c->u_old = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->u_old, x+3, y+3, z+3);

  c->u_new = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->u_new, x+3, y+3, z+3);

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
}

void read_meshdata(double ***u_old, int x, int y, int z)
{
  int i, j, k;

  FILE *fp = fopen ("mesh_new/128x128x128.txt","r");

  if(fp==NULL)
  {
      printf("\nFailure to open input file.\n");
      exit(0);
  }

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        fscanf(fp,"%lf ", &u_old[i][j][k]);
        //printf("%f ", tensor_x[i][j][k]);
      }
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }
}

void init_tensor(cube c, int x, int y, int z)
{
  int i, j, k;
  double delta_x, delta_y, delta_z;

  double x_step = (3 - 1)/(double)x;
  double y_step = (2 - 1)/(double)y;
  double z_step = (5 - 4)/(double)z;

  double *grid_x = (double*)calloc(x+2, sizeof(double));
  double *grid_y = (double*)calloc(y+2, sizeof(double));
  double *grid_z = (double*)calloc(z+2, sizeof(double));

  for(i = 1; i <= x+1; i++)
  {
    grid_x[i] = 1 + x_step*(i-1);
  }

  for(i = 1; i <= y+1; i++)
  {
    grid_y[i] = 1 + y_step*(i-1);
  }

  for(i = 1; i <= z+1; i++)
  {
    grid_z[i] = 4 + z_step*(i-1);
  }

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {

          if((c.u_old[i-1][j][k] != 0 && c.u_old[i+1][j][k] != 0) &&
          (c.u_old[i][j-1][k] != 0 && c.u_old[i][j+1][k] != 0) &&
          (c.u_old[i][j][k-1] != 0 && c.u_old[i][j][k+1] != 0) &&
          c.u_old[i][j][k] != 0)
          {
            c.tensor_x0[i][j][k] = 1.0;
            c.tensor_y0[i][j][k] = 1.0;
            c.tensor_z0[i][j][k] = 1.0;
          }
          else
          {
            c.tensor_x0[i][j][k] = 0.0;
            c.tensor_y0[i][j][k] = 0.0;
            c.tensor_z0[i][j][k] = 0.0;
          }
      }
    }
  }
}

/*void calculate_center(double ***u_old, double *center_x, double *center_y, double *center_z)
{
  int i, j, k;
  int a, b, c, d, e, f, g, h;
  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {

      }
    }
  }
}*/

int main(int argc, char *argv[])
{
  int x = 128;
  int y = 128;
  int z = 128;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  int i, j, k;
  cube c;
  tensorfield *T = (tensorfield *)malloc(sizeof(tensorfield));
  //char tensorfile[100] = "mesh_new/3Dheart.1";
  //char *tensorfiledsadas;



  init_meshdata(&c, x, y, z);
  read_meshdata(c.u_old, x, y, z);
  sparse_readtensorfiles("mesh_new/3Dheart.1", T, 1000);
  fiberstotensors(T);
  simple_averagetensors(c, T, x, y, z);
  //init_tensor(c.u_old, c.tensor_x1, c.tensor_y1, c.tensor_z1, x, y, z);
  //init_tensor(c, x, y, z);

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        if(c.u_old[i][j][k] == 1)
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

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        if(c.tensor_x1[i][j][k] == 1)
        {
          count_inside_tensor++;
        }
        else
        {
          count_outside_tensor++;
        }
      }
    }
  }*/

  double l2_uold = 0;

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        l2_uold += c.u_old[i][j][k]*c.u_old[i][j][k];
      }
    }
  }

  double l2_tensor = 0;

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        l2_tensor += c.tensor_x1[i][j][k]*c.tensor_x1[i][j][k];
      }
    }
  }*/

  printf("DEBUG: \n");
  printf("U_OLD: inside points: %d \t outside points: %d sum check: %0.12f \n", count_inside, count_outside, l2_uold);
  //printf("TENSOR: inside points: %d \t outside points: %d sum check: %0.12f \n", count_inside_tensor, count_outside_tensor, l2_tensor);
  return 0;
}
