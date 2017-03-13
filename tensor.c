#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tensor.h"
#include "binary_cube_sequential.h"

const double SIGMA_l =  1.21321;
const double SIGMA_t =  0.2121;

void sparse_readtensorfiles(char* tensorfile,tensorfield* T,int skip)
{
    char fname[200];
    char suffix[20];
    FILE *fpdata;
    double dump;
    int err;
    strcpy(fname,tensorfile);
    strcpy(suffix,".fibers");
    strcat(fname,suffix);
    fpdata = fopen(fname, "r");
    if(fpdata==NULL)
    {
        printf("\nFailure to open input file %s .\n",fname);
        exit(0);
    }
    err = fscanf(fpdata, "%d", &T->numtensor);
    int total=T->numtensor;
    T->numtensor/=skip;
    T->coord = (double*)malloc(3*T->numtensor*sizeof(double));
    T->fibers = (double*)malloc(3*T->numtensor*sizeof(double));
    int pos=0;
    for (int i = 0; i < T->numtensor; i++)
    {
        if(pos<total)
        for(int j=0;j<3;j++)
            err = fscanf(fpdata, "%lf", &T->fibers[i*3+j]);
        else
            break;
        for(int j=0;j<3;j++)
            err = fscanf(fpdata, "%lf", &dump);

    }
    fclose(fpdata);

    strcpy(fname,tensorfile);
    strcpy(suffix,".tensorcoord");
    strcat(fname,suffix);
    fpdata = fopen(fname, "r");
    if(fpdata==NULL)
    {
        printf("\nFailure to open input file %s .\n",fname);
        exit(0);
    }
    pos=0;
    for (int i = 0; i < T->numtensor; i++)
    {
        if(pos<total)
        for(int j=0;j<3;j++)
            err = fscanf(fpdata, "%lf", &T->coord[i*3+j]);
        else
            break;
        for(int j=0;j<3;j++)
            err = fscanf(fpdata, "%lf", &dump);
    }
    fclose(fpdata);

}



void fiberstotensors(tensorfield* T)
{
    //nonzeroes on main diagonal
    // +3: (1,2) , (2,1)
    // +4: (1,3) , (3,1)
    // +5: (2,3) , (3,2)

    //double start=omp_get_wtime();
    T->inputtensor = (double*)malloc(6*T->numtensor*sizeof(double));
    for (int i = 0; i < T->numtensor; i++)
    {
        for(int j=0;j<3;j++)
          T->inputtensor[i*6+j]=((SIGMA_l-SIGMA_t)*T->fibers[i*3+j]*T->fibers[i*3+j])+SIGMA_t;

      T->inputtensor[i*6+3]=((SIGMA_l-SIGMA_t)*T->fibers[i*3+0]*T->fibers[i*3+1]);
      T->inputtensor[i*6+4]=((SIGMA_l-SIGMA_t)*T->fibers[i*3+0]*T->fibers[i*3+2]);
      T->inputtensor[i*6+5]=((SIGMA_l-SIGMA_t)*T->fibers[i*3+1]*T->fibers[i*3+2]);
    }
    free(T->fibers);
    //double runtime=omp_get_wtime()-start;
    //printf("Fiberstotensors. Time taken: %lf \n",runtime);
}


/*void simple_averagetensors(Nodemeshdata* M,tensorfield* T)
{
    double lambda=-1.0;  //store negative lambda
    M->tensor = (double*)calloc(6*M->numtet,sizeof(double));
    //nonzeroes on main diagonal
    // +3: (1,2) , (2,1)
    // +4: (1,3) , (3,1)
    // +5: (2,3) , (3,2)
    double start=omp_get_wtime();

    #pragma omp parallel for
    for(int i=0; i<M->numtet; i++)
    {
        double norm=0.0;
        for(int j=0; j<T->numtensor; j++)
        {
            double d=
(M->centroid[i*3+0]-T->coord[j*3+0])*(M->centroid[i*3+0]-T->coord[j*3+0])+
(M->centroid[i*3+1]-T->coord[j*3+1])*(M->centroid[i*3+1]-T->coord[j*3+1])+
(M->centroid[i*3+2]-T->coord[j*3+2])*(M->centroid[i*3+2]-T->coord[j*3+2]);
            double e=exp(lambda*d);
            norm+=e;
            M->tensor[6*i+0]+=T->inputtensor[6*j+0]*e;
            M->tensor[6*i+1]+=T->inputtensor[6*j+1]*e;
            M->tensor[6*i+2]+=T->inputtensor[6*j+2]*e;
            M->tensor[6*i+3]+=T->inputtensor[6*j+3]*e;
            M->tensor[6*i+4]+=T->inputtensor[6*j+4]*e;
            M->tensor[6*i+5]+=T->inputtensor[6*j+5]*e;
        }

        M->tensor[6*i+0]/=norm;
        M->tensor[6*i+1]/=norm;
        M->tensor[6*i+2]/=norm;
        M->tensor[6*i+3]/=norm;
        M->tensor[6*i+4]/=norm;
        M->tensor[6*i+5]/=norm;

        //printf("%d : %e %e %e  %e %e %e \n",i,M->tensor[6*i+0],M->tensor[6*i+1],M->tensor[6*i+2],M->tensor[6*i+3],M->tensor[6*i+4],M->tensor[6*i+5]);
    }
    double runtime=omp_get_wtime()-start;
    printf("Averagetensors. Time taken: %lf \n",runtime);
}*/

void simple_averagetensors(cube c,tensorfield* T, int x, int y, int z)
{
    double lambda=-1.0;  //store negative lambda
    int i, j, k, l;
    double d, e, norm;
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
              norm = 0.0;
              //for(l = 0; l < T->numtensor; l++)
              for(l = 0; l < 10; l++)
              {

                d = ((grid_x[i]+(x_step/2))-T->coord[l*3+0])*((grid_x[i]+(x_step/2))-T->coord[l*3+0])+
                    ((grid_y[j+1])-T->coord[l*3+1])*((grid_y[j+1])-T->coord[l*3+1])+
                    ((grid_z[k]+(z_step/2))-T->coord[l*3+2])*((grid_z[k]+(z_step/2))-T->coord[l*3+2]);

                e = exp(lambda*d);
                norm+=e;
                c.tensor_x0[i][j][k]+=T->inputtensor[6*l+0]*e;
                c.tensor_x1[i][j][k]+=T->inputtensor[6*l+1]*e;
                c.tensor_y0[i][j][k]+=T->inputtensor[6*l+2]*e;
                c.tensor_y1[i][j][k]+=T->inputtensor[6*l+3]*e;
                c.tensor_z0[i][j][k]+=T->inputtensor[6*l+4]*e;
                c.tensor_z1[i][j][k]+=T->inputtensor[6*l+5]*e;
              }

              c.tensor_x0[i][j][k]/=norm;
              c.tensor_x1[i][j][k]/=norm;
              c.tensor_y0[i][j][k]/=norm;
              c.tensor_y1[i][j][k]/=norm;
              c.tensor_z0[i][j][k]/=norm;
              c.tensor_z1[i][j][k]/=norm;
            }
            else
            {
              c.tensor_x0[i][j][k] = 0.0;
              c.tensor_x1[i][j][k] = 0.0;
              c.tensor_y0[i][j][k] = 0.0;
              c.tensor_y1[i][j][k] = 0.0;
              c.tensor_z0[i][j][k] = 0.0;
              c.tensor_z1[i][j][k] = 0.0;
            }
        }
      }
    }
}

/*tensorfield *T = (tensorfield *)malloc(sizeof(tensorfield));

    sparse_readtensorfiles(tensorfile,T,1000);
    fiberstotensors(T);
    simple_averagetensors(M,T);*/
