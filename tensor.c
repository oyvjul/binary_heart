#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tensor.h"

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

void simple_averagetensors(cube c,tensorfield* T)
{
    double lambda=-1.0;  //store negative lambda
    //M->tensor = (double*)calloc(6*M->numtet,sizeof(double));
    //nonzeroes on main diagonal
    // +3: (1,2) , (2,1)
    // +4: (1,3) , (3,1)
    // +5: (2,3) , (3,2)
    //double start=omp_get_wtime();

    //#pragma omp parallel for
    /*for(int i=0; i<M->numtet; i++)
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
    }*/
    //double runtime=omp_get_wtime()-start;
    //printf("Averagetensors. Time taken: %lf \n",runtime);
}

/*tensorfield *T = (tensorfield *)malloc(sizeof(tensorfield));

    sparse_readtensorfiles(tensorfile,T,1000);
    fiberstotensors(T);
    simple_averagetensors(M,T);*/

int main(int argc, char *argv[])
{

  return 0;
}
