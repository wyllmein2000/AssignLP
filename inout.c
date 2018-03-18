//
//
#include <stdio.h>
#include <math.h>

void readdata (char *fn, float *p, int n) {
   FILE *fp = fopen(fn,"rb");
   for (int i=0; i<n; i++)
       fread (&p[i],sizeof(float),1,fp);
   fclose(fp);
}

void writedata (char *fn, float *p, int n) {
     FILE *fp = fopen(fn,"wb");
     for (int i=0; i<n; i++)
         fwrite(&p[i],sizeof(float),1,fp);
     fclose(fp);
}
