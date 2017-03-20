#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <stdio.h>
#include <stdlib.h>

double divergence_direction_x(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k);
double divergence_direction_y(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k);
double divergence_direction_z(double ***u, double tensor_x, double tensor_y, double tensor_z, double delta_x, double delta_y, double delta_z, int i, int j, int k);


#endif
