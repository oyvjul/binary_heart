#include "diffusion.h"
//#include <math.h>
//#include "binary_cube_sequential.h"

double cell_direction_x, cell_direction_y, cell_direction_z;
double upper_right_one, lower_right_one, upper_right, lower_right, upper_left_one, lower_left_one, upper_left, lower_left;

//(i+1/2, j+1/2, k+1/2)
double flux_upper_right_one(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i+1][j+1][k] + u[i+1][j][k] -
                u[i][j+1][k+1] - u[i][j][k+1] - u[i][j+1][k] - u[i][j][k])/(2*delta_x);

  cell_direction_y = (u[i+1][j+1][k+1] + u[i][j+1][k+1] + u[i+1][j+1][k] + u[i][j+1][k] -
                u[i+1][j][k+1] - u[i][j][k+1] - u[i+1][j][k] - u[i][j][k])/(2*delta_y);

  cell_direction_z = (u[i+1][j+1][k+1] + u[i+1][j][k+1] + u[i][j+1][k+1] + u[i][j][k+1] -
                u[i+1][j+1][k] - u[i+1][j][k] - u[i][j+1][k] - u[i][j][k])/(2*delta_z);

  return ((tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z));
}

//(i+1/2, j+1/2, k-1/2)//DONE
double flux_upper_right(double ***u, double tensor_x, double tensor_y, double tensor_z,
                        double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i+1][j+1][k] + u[i+1][j][k] + u[i+1][j+1][k-1] + u[i+1][j][k-1] -
                u[i][j+1][k] - u[i][j][k] - u[i][j+1][k-1] - u[i][j][k-1])/(2*delta_x);

  cell_direction_y = (u[i+1][j+1][k] + u[i][j+1][k] + u[i+1][j+1][k-1] + u[i][j+1][k-1] -
                u[i+1][j][k] - u[i][j][k] - u[i+1][j][k-1] - u[i][j][k-1])/(2*delta_y);

  cell_direction_z = (u[i+1][j+1][k] + u[i+1][j][k] + u[i][j+1][k] + u[i][j][k] -
                u[i+1][j+1][k-1] - u[i+1][j][k-1] - u[i][j+1][k-1] - u[i][j][k-1])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i+1/2, j-1/2, k-1/2)//DONE
double flux_lower_right(double ***u, double tensor_x, double tensor_y, double tensor_z,
                        double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i+1][j][k] + u[i+1][j-1][k] + u[i+1][j][k-1] + u[i+1][j-1][k-1] -
                u[i][j][k] - u[i][j-1][k] - u[i][j][k-1] - u[i][j-1][k-1])/(2*delta_x);

  cell_direction_y = (u[i+1][j][k] + u[i][j][k] + u[i+1][j][k-1] + u[i][j][k-1] -
                u[i+1][j-1][k] - u[i][j-1][k] - u[i+1][j-1][k-1] - u[i][j-1][k-1])/(2*delta_y);

  cell_direction_z = (u[i+1][j][k] + u[i+1][j-1][k] + u[i][j][k] + u[i][j-1][k] -
                u[i+1][j][k-1] - u[i+1][j-1][k-1] - u[i][j][k-1] - u[i][j-1][k-1])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i+1/2, j-1/2, k+1/2)//DONE
double flux_lower_right_one(double ***u, double tensor_x, double tensor_y, double tensor_z,
                            double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i+1][j][k] + u[i+1][j-1][k] -
                u[i][j][k+1] - u[i][j-1][k+1] - u[i][j][k] - u[i][j-1][k])/(2*delta_x);

  cell_direction_y = (u[i+1][j][k+1] + u[i][j][k+1] + u[i+1][j][k] + u[i][j][k] -
                u[i+1][j-1][k+1] - u[i][j-1][k+1] - u[i+1][j-1][k] - u[i][j-1][k])/(2*delta_y);

  cell_direction_z = (u[i+1][j][k+1] + u[i+1][j-1][k+1] + u[i][j][k+1] + u[i][j-1][k+1] -
                u[i+1][j][k] - u[i+1][j-1][k] - u[i][j][k] - u[i][j-1][k])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i-1/2, j+1/2, k+1/2)//DONE
double flux_upper_left_one(double ***u, double tensor_x, double tensor_y, double tensor_z,
                           double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i][j+1][k+1] + u[i][j][k+1] + u[i][j+1][k] + u[i][j][k] -
                u[i-1][j+1][k+1] - u[i-1][j][k+1] - u[i-1][j+1][k] - u[i-1][j][k])/(2*delta_x);

  cell_direction_y = (u[i][j+1][k+1] + u[i-1][j+1][k+1] + u[i][j+1][k] + u[i-1][j+1][k] -
                u[i][j][k+1] - u[i-1][j][k+1] - u[i][j][k] - u[i-1][j][k])/(2*delta_y);

  cell_direction_z = (u[i][j+1][k+1] + u[i][j][k+1] + u[i-1][j+1][k+1] + u[i-1][j][k+1] -
                u[i][j+1][k] - u[i][j][k] - u[i-1][j+1][k] - u[i-1][j][k])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i-1/2, j+1/2, k-1/2)//DONE
double flux_upper_left(double ***u, double tensor_x, double tensor_y, double tensor_z,
                       double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i][j+1][k] + u[i][j][k] + u[i][j+1][k-1] + u[i][j][k-1] -
                u[i-1][j+1][k] - u[i-1][j][k] - u[i-1][j+1][k-1] - u[i-1][j][k-1])/(2*delta_x);

  cell_direction_y = (u[i][j+1][k] + u[i-1][j+1][k] + u[i][j+1][k-1] + u[i-1][j+1][k-1] -
                u[i][j][k] - u[i-1][j][k] - u[i][j][k-1] - u[i-1][j][k-1])/(2*delta_y);

  cell_direction_z = (u[i][j+1][k] + u[i][j][k] + u[i-1][j+1][k] + u[i-1][j][k] -
                u[i][j+1][k-1] - u[i][j][k-1] - u[i-1][j+1][k-1] - u[i-1][j][k-1])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i-1/2, j-1/2, k+1/2)//DONE
double flux_lower_left_one(double ***u, double tensor_x, double tensor_y, double tensor_z,
                           double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i][j][k+1] + u[i][j-1][k+1] + u[i][j][k] + u[i][j-1][k] -
                u[i-1][j][k+1] - u[i-1][j-1][k+1] - u[i-1][j-1][k] - u[i-1][j-1][k])/(2*delta_x);

  cell_direction_y = (u[i][j][k+1] + u[i-1][j][k+1] + u[i][j][k] + u[i-1][j][k] -
                u[i][j-1][k+1] - u[i-1][j-1][k+1] - u[i][j-1][k] - u[i-1][j-1][k])/(2*delta_y);

  cell_direction_z = (u[i][j][k+1] + u[i][j-1][k+1] + u[i-1][j][k+1] + u[i-1][j-1][k+1] -
                u[i][j][k] - u[i][j-1][k] - u[i-1][j][k] - u[i-1][j-1][k])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

//(i-1/2, j-1/2, k-1/2)
double flux_lower_left(double ***u, double tensor_x, double tensor_y, double tensor_z,
                       double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  cell_direction_x = (u[i][j][k] + u[i][j-1][k] + u[i][j][k-1] + u[i][j-1][k-1] -
                u[i-1][j][k] - u[i-1][j-1][k] - u[i-1][j][k-1] - u[i-1][j-1][k-1])/(2*delta_x);

  cell_direction_y = (u[i][j][k] + u[i-1][j][k] + u[i][j][k-1] + u[i-1][j][k-1] -
                u[i][j-1][k] - u[i-1][j-1][k] - u[i][j-1][k-1] - u[i-1][j-1][k-1])/(2*delta_y);

  cell_direction_z = (u[i][j][k] + u[i][j-1][k] + u[i-1][j][k] + u[i-1][j-1][k] -
                u[i][j][k-1] - u[i][j-1][k-1] - u[i-1][j][k-1] - u[i-1][j-1][k-1])/(2*delta_z);

  return (tensor_x*cell_direction_x) + (tensor_y*cell_direction_y) + (tensor_z*cell_direction_z);
}

double divergence_cell_direction_x(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  upper_right_one = flux_upper_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right_one = flux_lower_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_right = flux_upper_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right = flux_lower_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left_one = flux_upper_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left_one = flux_lower_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left = flux_upper_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left = flux_lower_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);

  return (upper_right_one + lower_right_one + upper_right + lower_right - upper_left_one - lower_left_one - upper_left - lower_left)/(2*delta_x);
}

double divergence_cell_direction_y(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  upper_right_one = flux_upper_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left_one = flux_upper_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_right = flux_upper_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left = flux_upper_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right_one = flux_lower_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left_one = flux_lower_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right =  flux_lower_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left = flux_lower_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);

  return (upper_right_one + upper_left_one + upper_right + upper_left - lower_right_one - lower_left_one - lower_right - lower_left)/(2*delta_y);
}

double divergence_cell_direction_z(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k)
{
  upper_right_one = flux_upper_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right_one = flux_lower_right_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left_one = flux_upper_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left_one = flux_lower_left_one(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_right = flux_upper_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_right = flux_lower_right(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  upper_left = flux_upper_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);
  lower_left = flux_lower_left(u, tensor_x, tensor_y, tensor_z, delta_x, delta_y, delta_z, i, j, k);

  return (upper_right_one + lower_right_one + upper_left_one + lower_left_one - upper_right - lower_right - upper_left - lower_left)/(2*delta_z);
}
