#include <stdio.h>
#include "tensor.h"
//#include "binary_cube_sequential.h"
#include "read_mesh.h"
#include "io.h"
#include "diffusion.h"
#include "omp.h"
#include "math.h"

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

void init_sequential_data(cube *c, int x, int y, int z)
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

  c->grid_x = (double*)calloc(x+2, sizeof(double));
  c->grid_y = (double*)calloc(y+2, sizeof(double));
  c->grid_z = (double*)calloc(z+2, sizeof(double));

  c->center_x = (double*)calloc(x+2, sizeof(double));
  c->center_y = (double*)calloc(y+2, sizeof(double));
  c->center_z = (double*)calloc(z+2, sizeof(double));

  c->x = x;
  c->y = y;
  c->z = z;
}

int main(int argc, char *argv[])
{
  int x = 10;
  int y = 10;
  int z = 10;
  double count = 1.0;
  int i, j, k;
  int count_inside, count_outside;
  int count_inside_tensor, count_outside_tensor;
  cube *c = (cube*)malloc(sizeof(cube));
  meshdata *m = (meshdata*)malloc(sizeof(meshdata));
  tensorfield *t = (tensorfield *)malloc(sizeof(tensorfield));

  init_sequential_data(c, x, y, z);
  init_cube_grid(c, m);



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

    printf("TENSOR: inside points: %d \t outside points: %d sum check: %d  total: %d\n", count_inside_tensor, count_outside_tensor, c->x*c->y*c->z, count_inside);

    /*for(i = 1; i <= z; i++)
    {
      for(j = 1; j <= y; j++)
      {
        for(k = 1; k <= x; k++)
        {
          printf("%f ", c->tensor_x0[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }*/

    write_binaryformat("mesh_new/tensor_x0.tensor", c->tensor_x0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_x1.tensor", c->tensor_x1, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_y0.tensor", c->tensor_y0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_y1.tensor", c->tensor_y1, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_z0.tensor", c->tensor_z0, x+3, y+3, z+3);
    write_binaryformat("mesh_new/tensor_z1.tensor", c->tensor_z1, x+3, y+3, z+3);
    //read_binaryformat("mesh_new/test.tensor", &c->tensor_x0, &x+3, &y+3, &z+3);
  #else

  /*double W = 2.0;
  double pi =  3.14159265358979323846;
  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        c->u_old[i][j][k] = W;
        c->tensor_x0[i][j][k] = 0.1;
        count++;
      }
    }
  }*/

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
  //write_binaryformat("mesh_new/test.tensor", c->u_old, x+3, y+3, z+3);
  //read_binaryformat("mesh_new/tensor_x0.tensor", &c->tensor_x0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_x1.tensor", &c->tensor_x1, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_y0.tensor", &c->tensor_y0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_y1.tensor", &c->tensor_y1, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_z0.tensor", &c->tensor_z0, x+3, y+3, z+3);
  read_binaryformat("mesh_new/tensor_z1.tensor", &c->tensor_z1, x+3, y+3, z+3);

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
  int kx = 0.5;

  //exp(-t*3*W*W*pi*pi*k)*cos(pi*x*W)*cos(pi*y*W)*cos(pi*z*W);
  double start1 = omp_get_wtime();
  for(int l = 1; l <= 2; l++)
  {
  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        /*c->u_new[i][j][k] = divergence_cell_direction_x(c->u_old, c->tensor_x0[i][j][k], c->tensor_x1[i][j][k], c->tensor_y0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k)
                          + divergence_cell_direction_y(c->u_old, c->tensor_x1[i][j][k], c->tensor_y1[i][j][k], c->tensor_z0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k)
                          + divergence_cell_direction_x(c->u_old, c->tensor_y0[i][j][k], c->tensor_z0[i][j][k], c->tensor_z1[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k);*/

        /*c->u_new[i][j][k] = c->u_old[i][j][k] + (1*(divergence_cell_direction_x(c->u_old, c->tensor_x0[i][j][k], 0, 0, 1, 1, 1, i, j, k)
                          + divergence_cell_direction_y(c->u_old, 0, c->tensor_x0[i][j][k], 0, 1, 1, 1, i, j, k)
                          + divergence_cell_direction_x(c->u_old, 0, 0, c->tensor_x0[i][j][k], 1, 1, 1, i, j, k)));*/

        c->u_new[i][j][k] =  ((divergence_cell_direction_x(c->u_old, c->tensor_x0[i][j][k], 0, 0, c->x_step, c->y_step, c->y_step, i, j, k)
                            + divergence_cell_direction_y(c->u_old, 0, c->tensor_x0[i][j][k], 0, c->x_step, c->y_step, c->y_step, i, j, k)
                            + divergence_cell_direction_x(c->u_old, 0, 0, c->tensor_x0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k)));
      }
    }
  }

  double ***temp;
  temp = c->u_old;
  c->u_old = c->u_new;
  c->u_new = temp;
}


/*for(int l = 0; l < 40; l++)
{
for(i = 1; i <= z+1; i++)
{
  for(j = 1; j <= y+1; j++)
  {
    for(k = 1; k <= x+1; k++)
    {
      c->u_new[i][j][k] = divergence_cell_direction_x(c->u_old, c->tensor_x0[i][j][k], c->tensor_x1[i][j][k], c->tensor_y0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k);
    }

    for(k = 1; k <= x+1; k++)
    {
      c->u_new[i][j][k] += divergence_cell_direction_y(c->u_old, c->tensor_x0[i][j][k], c->tensor_x1[i][j][k], c->tensor_y0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k);
    }

    for(k = 1; k <= x+1; k++)
    {
      c->u_new[i][j][k] += divergence_cell_direction_z(c->u_old, c->tensor_x0[i][j][k], c->tensor_x1[i][j][k], c->tensor_y0[i][j][k], c->x_step, c->y_step, c->y_step, i, j, k);
    }
  }
}

double ***temp;
temp = c->u_old;
c->u_old = c->u_new;
c->u_new = temp;
}*/
  double ***debug;
  debug = dallocate_3d(c->x+3, c->y+3, c->z+3);
  dinit_3d(debug, c->x+3, c->y+3, c->z+3);
  //double pi =  3.14159265358979323846;

  double end1 = omp_get_wtime();
  printf("it took: %0.12f \n", end1-start1);
  double test;

  for(int l = 1; l <= 1; l++)
  {
    for(i = 1; i <= z+1; i++)
    {
      for(j = 1; j <= y+1; j++)
      {
        for(k = 1; k <= x+1; k++)
        {
          double x_coord = c->grid_x[k]+(c->x_step/2);
          double y_coord = c->grid_y[j]+(c->y_step/2);
          double z_coord = c->grid_z[i]+(c->z_step/2);
          debug[i][j][k] = exp(-0.1*3*W*W*pi*pi*0.2)*cos(pi*x_coord*W)*cos(pi*y_coord*W)*cos(pi*z_coord*W);
          //printf("%f ", test);
          //printf("%f ", c->u_old[i][j][k]);
        }
        //printf("\n");
      }
      //printf("\n");
    }
  }

  for(i = 0; i <= 3; i++)
  {
    for(j = 0; j <= y+2; j++)
    {
      for(k = 0; k <= x+2; k++)
      {
        //debug[i][j][k] = exp(-l*3*W*W*pi*pi*c->tensor_x0[i][j][k])*cos(pi*x*W)*cos(pi*y*W)*cos(pi*z*W);
        //printf("%f ", debug[i][j][k]);
        printf("%f ", c->u_old[i][j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }
    //p
  /*for(i = 1; i <= z+1; i++)
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
