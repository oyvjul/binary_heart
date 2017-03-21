#include <stdio.h>
#include "io.h"

void read_binaryformat(char* filename, double ****matrix, int x, int y, int z)
{
    int i;
    FILE* fp = fopen (filename,"rb");
    /*fread (x, sizeof(int), 1, fp);
    fread (y, sizeof(int), 1, fp);
    fread (z, sizeof(int), 1, fp);*/
    fread ((*matrix)[0][0], sizeof(double), x*y*z, fp);
    fclose (fp);
}

void write_binaryformat(char* filename, double ***matrix, int x, int y, int z)
{
    FILE *fp = fopen (filename,"wb");
    /*fwrite (&x, sizeof(int), 1, fp);
    fwrite (&y, sizeof(int), 1, fp);
    fwrite (&z, sizeof(int), 1, fp);*/
    fwrite (matrix[0][0], sizeof(double), x*y*z, fp);
    fclose (fp);
}
