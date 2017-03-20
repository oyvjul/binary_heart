#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include "binary_cube_sequential.h"

int direction_x, direction_y, direction_z;

int main(int agrc, char *argv[])
{

  return 0;
}

//(i+1/2, j+1/2, k+1/2)
double cell_upper_right_one(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = (u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i+1][j+1][k] + u[i+1][j][k] -
                u[i][j+1][k+1] - u[i][j][k+1] - u[i][j+1][k] - u[i][j][k])/(2*delta_x);

  direction_y = (u[i+1][j+1][k+1] + u[i][j+1][k+1] + u[i+1][j+1][k] + u[i][j+1][k] -
                u[i+1][j][k+1] - u[i][j][k+1] - u[i+1][j][k] - u[i][j][k])/(2*delta_y);

  direction_z = (u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i][j+1][k+1] + u[i][j][k+1] -
                u[i+1][j+1][k] - u[i+1][j][k] - u[i][j+1][k] - u[i][j][k])/(2*delta_z);

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i+1/2, j+1/2, k-1/2)
double cell_upper_right(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = (u[i+1][j+1][k] + u[i+1][j][k] + u[i+1][j+1][k-1] + u[i+1][j][k-1] -
                u[i][j+1][k] - u[i][j][k] - u[i][j+1][k-1] - u[i][j][k-1])/(2*delta_x);

  direction_y = (u[i+1][j+1][k] + u[i][j+1][k] + u[i+1][j+1][k-1] + u[i][j+1][k-1] -
                u[i+1][j][k] - u[i][j][k] - u[i+1][j][k-1] - u[i][j][k-1])/(2*delta_y);

  direction_z = (u[i+1][j+1][k] + u[i+1][j][k] + u[i][j+1][k] + u[i][j][k] -
                u[i+1][j+1][k-1] - u[i+1][j][k-1] - u[i][j+1][k-1] - u[i][j][k-1])/(2*delta_z);

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i+1/2, j-1/2, k-1/2)
double cell_lower_right(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = (u[i+1][j][k] + u[i+1][j-1][k] + u[i+1][j][k-1] + u[i+1][j-1][k-1] -
                u[i][j][k] - u[i][j-1][k] - u[i][j][k-1] - u[i][j-1][k-1])/(2*delta_x);

  direction_y = u[i+1][j][k] + u[i][j][k] + u[i+1][j][k-1] + u[i][j][k-1] -
                u[i+1][j-1][k] - u[i][j-1][k] - u[i+1][j-1][k-1] - u[i][j-1][k-1];

  direction_z = u[i+1][j][k] + u[i+1][j-1][k] + u[i][j][k] + u[i][j-1][k] -
                u[i+1][j][k-1] - u[i+1][j-1][k-1] - u[i][j][k-1] - u[i][j-1][k-1];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i+1/2, j-1/2, k+1/2)
double cell_lower_right_one(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i+1][j][k] + u[i+1][j-1][k] -
                u[i][j][k+1] - u[i][j-1][k+1] - u[i][j][k] - u[i][j-1][k];

  direction_y = u[i+1][j][k+1] + u[i][j][k+1] + u[i+1][j][k] + u[i][j][k] -
                u[i+1][j-1][k+1] - u[i][j-1][k+1] - u[i+1][j-1][k] - u[i][j-1][k];

  direction_z = u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i][j][k+1] + u[i][j-1][k+1] -
                u[i+1][j][k] - u[i+1][j-1][k] - u[i][j][k] - u[i][j-1][k];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i-1/2, j+1/2, k+1/2)
double cell_upper_left_one(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = u[i][j+1][k+1] + u[i][j][k+1] + u[i][j+1][k] + u[i][j][k] -
                u[i-1][j+1][k+1] - u[i-1][j][k+1] - u[i-1][j+1][k] - u[i-1][j][k];

  direction_y = u[i][j+1][k+1] + u[i-1][j+1][k+1] + u[i][j+1][k] + u[i-1][j+1][k] -
                u[i][j][k+1] - u[i-1][j][k+1] - u[i][j][k] - u[i-1][j][k];

  direction_z = u[i][j+1][k+1] + u[i][j][k+1] + u[i-1][j+1][k+1] + u[i-1][j][k+1] -
                u[i][j+1][k] - u[i][j][k] - u[i-1][j+1][k] - u[i-1][j][k];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i-1/2, j+1/2, k-1/2)
double cell_upper_left(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = u[i][j+1][k] + u[i][j][k] + u[i][j+1][k-1] + u[i][j][k-1] -
                u[i-1][j+1][k] - u[i-1][j][k] - u[i-1][j+1][k-1] - u[i-1][j][k-1];

  direction_y = u[i][j+1][k] + u[i-1][j+1][k] + u[i][j+1][k-1] + u[i-1][j+1][k-1] -
                u[i][j][k] - u[i-1][j][k] - u[i][j][k-1] - u[i-1][j][k-1];

  direction_z = u[i][j+1][k] + u[i][j][k] + u[i-1][j+1][k] + u[i-1][j][k] -
                u[i][j+1][k-1] - u[i][j][k-1] - u[i-1][j+1][k-1] - u[i-1][j][k-1];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i-1/2, j-1/2, k+1/2)
double cell_lower_left_one(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = u[i][j][k+1] + u[i][j-1][k+1] + u[i][j][k] + u[i][j-1][k] -
                u[i-1][j][k+1] - u[i-1][j-1][k+1] - u[i-1][j-1][k] - u[i-1][j-1][k];

  direction_y = u[i][j][k+1] + u[i-1][j][k+1] + u[i][j][k] + u[i-1][j][k] -
                u[i][j-1][k+1] - u[i-1][j-1][k+1] - u[i][j-1][k] - u[i-1][j-1][k];

  direction_z = u[i][j][k+1] + u[i][j-1][k+1] + u[i-1][j][k+1] + u[i-1][j-1][k+1] -
                u[i][j][k] - u[i][j-1][k] - u[i-1][j][k] - u[i-1][j-1][k];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}

//(i-1/2, j-1/2, k-1/2)
double cell_lower_left(double ***u, double *** tensor_x, double ***tensor_y, double ***tensor_z, int i, int j, int k)
{
  direction_x = u[i][j][k] + u[i][j-1][k] + u[i][j][k-1] + u[i][j-1][k-1] -
                u[i-1][j][k] - u[i-1][j-1][k] - u[i-1][j][k-1] - u[i-1][j-1][k-1];

  direction_y = u[i][j][k] + u[i-1][j][k] + u[i][j][k-1] + u[i-1][j][k-1] -
                u[i][j-1][k] - u[i-1][j-1][k] - u[i][j-1][k-1] - u[i-1][j-1][k-1];

  direction_z = u[i][j][k] + u[i][j-1][k] + u[i-1][j][k] + u[i-1][j-1][k] -
                u[i][j][k-1] - u[i][j-1][k-1] - u[i-1][j][k-1] - u[i-1][j-1][k-1];

  return (tensor_x*direction_x) + (tensor_y*direction_y) + (tensor_z*direction_z);
}
