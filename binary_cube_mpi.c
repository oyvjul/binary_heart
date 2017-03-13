#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct
{
  double ***u_old;
  double ***u_new;
  double ***tensor_x;
  double ***tensor_y;
  double ***tensor_z;
}cube;

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

  c->tensor_x = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_x, x+3, y+3, z+3);

  c->tensor_y = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_y, x+3, y+3, z+3);

  c->tensor_z = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(c->tensor_z, x+3, y+3, z+3);
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
  int x = 128;
  int y = 128;
  int z = 128;
  int nx, ny, nz;
  double ***heart_cube, ***u_old;
  double count = 1.0;

  MPI_Comm comm3d;
  MPI_Request req[12];
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  procs_x = 2;
  procs_y = 2;
  procs_z = 2;

  dims[0] = procs_x;
  dims[1] = procs_y;
  dims[2] = procs_z;

  periods[0] = 0;
  periods[1] = 0;
  periods[2] = 0;

  if(procs_x*procs_y*procs_z != size)
  {
    if(rank == 0)
      fprintf(stderr, "Product of grid block dimensions must match the number of processes\n");
  }

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm3d);
  MPI_Cart_get(comm3d, 3, dims, periods, coords);
  MPI_Cart_shift(comm3d, 0, 1, &left, &right);
  MPI_Cart_shift(comm3d, 1, 1, &up, &down);
  MPI_Cart_shift(comm3d, 2, 1, &z_up, &z_down);

  decompose(x+1, dims[0], coords[0], &x0, &x1);
  decompose(y+1, dims[1], coords[1], &y0, &y1);
  decompose(z+1, dims[2], coords[2], &z0, &z1);

  start_x = 1;
  start_y = 1;
  start_z = 1;

  if(x0 == 0)
    x0 = 1;

  if(y0 == 0)
    y0 = 1;

  if(z0 == 0)
    z0 = 1;

  if(z1 == x)
    z1++;

  if(y1 == y)
    y1++;

  if(x1 == z)
    x1++;

  nx = x1 - x0 +1;
  ny = y1 - y0 +1;
  nz = z1 - z0 +1;

  if(left >= 0)
  {
    x0--;
    start_x = 0;
  }

  if(right >= 0)
    x1++;

  if(up >= 0)
  {
    y0--;
    start_y = 0;
  }

  if(down >= 0)
    y1++;

  if(z_up >= 0)
  {
    z0--;
    start_z = 0;
  }

  if(z_down >= 0)
    z1++;

  heart_cube = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(heart_cube, x+3, y+3, z+3);
  u_old = dallocate_3d(x+3, y+3, z+3);
  dinit_3d(u_old, x+3, y+3, z+3);

  /*for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        u_old[i][j][k] = count;
        count++;
      }
    }
  }*/

  //FILE *fp = fopen ("mesh_new/10x10x10_test.txt","r");
  FILE *fp = fopen ("mesh_new/128x128x128.txt","r");

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        fscanf(fp,"%lf ", &heart_cube[i][j][k]);
      }
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }

  ii = start_z;
  for(i = z0; i <= z1; i++)
  {
    jj = start_y;
    for(j = y0; j <= y1; j++)
    {
      kk = start_x;
      for(k = x0; k <= x1; k++)
      {
        u_old[ii][jj][kk] = heart_cube[i][j][k];
        kk++;
      }
      jj++;
    }
    ii++;
  }

  double l2_uold = 0;
  double all_sum = 0;

  for(i = 1; i <= nz; i++)
  {
    for(j = 1; j <= ny; j++)
    {
      for(k = 1; k <= nx; k++)
      {
        l2_uold += u_old[i][j][k]*u_old[i][j][k];
      }
    }
  }

  MPI_Allreduce(&l2_uold, &all_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(rank == 0)
  {
    printf("Result: %0.15f \n", all_sum);
    //printf("Iter: %d, r: %d", max_time, r);
  }

  /*double l2_tensor = 0;

  for(i = 1; i <= z+1; i++)
  {
    for(j = 1; j <= y+1; j++)
    {
      for(k = 1; k <= x+1; k++)
      {
        l2_tensor += tensor_x[i][j][k]*tensor_x[i][j][k];
      }
    }
  }*/

  if(rank == 0)
  {
    printf("rank %d, x0,x1,y0,y1,z0,z1:(%d, %d, %d, %d, %d, %d) nx: %d ny: %d, nz: %d \n", rank, x0, x1, y0, y1, z0, z1, nx, ny, nz);
    printf("rank: %d, l,r,u,d,zup,zdown(%d, %d, %d, %d, %d, %d) \n", rank, left, right, up, down, z_up, z_down);

    /*for(i = z0; i <= z1; i++)
    {
      for(j = y0; j <= y1; j++)
      {
        for(k = x0; k <= x1; k++)
        {
          printf("%f ", heart_cube[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }*/

    /*for(i = start_z; i <= nz+1; i++)
    {
      for(j = start_y; j <= ny+1; j++)
      {
        for(k = start_x; k <= nx+1; k++)
        {
          printf("%f ", u_old[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }*/
    printf("\n");
  }


  MPI_Finalize();
  return 0;
}
