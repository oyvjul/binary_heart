#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include "binary_cube_sequential.h"

int main(int agrc, char *argv[])
{

  return 0;
}

double cell_upper_right_one(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i+1][j+1][k] + u[i+1][j][k] -
                u[i][j+1][k+1] - u[i][j][k+1] - u[i][j+1][k] - u[i][j][k];

  direction_y = u[i+1][j+1][k+1] + u[i][j+1][k+1] + u[i+1][j+1][k] + u[i][j+1][k] -
                u[i+1][j][k+1] - u[i][j][k+1] - u[i+1][j][k] - u[i][j][k];

  direction_z = u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i][j+1][k+1] + u[i][j][k+1] -
                u[i+1][j+1][k] - u[i+1][j][k] - u[i][j+1][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_upper_right(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i+1][j+1][k] + u[i+1][j][k] + u[i+1][j+1][k-1] + u[i+1][j][k-1] -
                u[i][j+1][k] - u[i][j][k] - u[i][j+1][k-1] - u[i][j][k-1];

  direction_y = u[i+1][j+1][k] + u[i][j+1][k] + u[i+1][j+1][k-1] + u[i][j+1][k-1] -
                u[i+1][j][k] - u[i][j][k] - u[i+1][j][k-1] - u[i][j][k-1];

  direction_z = u[i+1][j+1][k] + u[i+1][j][k] + u[i][j+1][k] + u[i][j][k] -
                u[i+1][j+1][k-1] - u[i+1][j][k-1] - u[i][j+1][k-1] - u[i][j][k-1];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_lower_right(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_y = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_z = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_lower_right_one(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i+1][j][k] + u[i+1][j-1][k] -
                u[i][j][k+1] - u[i][j-1][k+1] - u[i][j][k] - u[i][j-1][k];

  direction_y = u[i+1][j][k+1] + u[i][j][k+1] + u[i+1][j][k] + u[i][j][k] -
                u[i+1][j-1][k+1] - u[i][j-1][k+1] - u[i+1][j-1][k] - u[i][j-1][k];

  direction_z = u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i][j][k+1] + u[i][j-1][k+1] -
                u[i+1][j][k] - u[i+1][j-1][k] - u[i][j][k] - u[i][j-1][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_upper_left(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_y = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_z = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_upper_left_one(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_y = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_z = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_lower_left(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_y = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_z = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}

double cell_lower_left_one(double ***u, int i, int j, int k)
{
  int direction_x, direction_y, direction_z;

  direction_x = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_y = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  direction_z = u[i][j][k] + u[i][j][k] + u[i][j][k] + u[i][j][k] -
                u[i][j][k] - u[i][j][k] - u[i][j][k] - u[i][j][k];

  return (tensor_x0*direction_x) + (tensor_x1*direction_y) + (tensor_x1*direction_z);
}
