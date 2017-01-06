#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <gsl/gsl_linalg.h>

void findFactor(double *result,double *a_data,double *P,int row,int col,int indexI,int indexJ);

void findAveragePx(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s)
{
   int i,j,k,l,m,istart,iend,jstart,jend,kstart,kend,indexI,indexJ,indexK,cnt;
   double minXDomain,minYDomain,tmp;
   double weight,x,y,z,wxl,wxr,wyl,wyr,aveP;
   double *w,*h,*gridP,*lambda;

   Particle ***particle;
   ptclList *p;
   particle=D->particle;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   weight=1.0/((double)D->numberInCell)/D->resolX/D->resolY;
   minXDomain=(double)D->minXDomain;
   minYDomain=(double)D->minYDomain;

   switch (D->dimension) {
   //2D
   case 2 :
     lambda=(double *)malloc(16*sizeof(double ));
     w=(double *)malloc(4*sizeof(double ));
     gridP=(double *)malloc(4*sizeof(double ));
     h=(double *)malloc(4*sizeof(double ));

     k=0;
//i=82;
//j=79;
     for(i=istart; i<iend; i++)  
       for(j=jstart; j<jend; j++)  
       {
         for(l=0; l<16; l++)  lambda[l]=0.0;
         for(l=0; l<4; l++)  h[l]=0.0;

         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           x=p->x-minXDomain;         y=p->y-minYDomain;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<16; l++)  {
             indexI=l%4;
             indexJ=l/4;
             lambda[l]+=w[indexI]*w[indexJ];
           }
           cnt++;
           p=p->next;
         }
         gridP[0]=den1[i][j][k];
         gridP[1]=den2[i][j][k];
         gridP[2]=den3[i][j][k];
         gridP[3]=den4[i][j][k];

//         if(cnt>0)  findFactor(h,lambda,gridP,4,4,i,j); else;
         findFactor(h,lambda,gridP,4,4,i,j);

         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           aveP=0.0;
           x=p->x-minXDomain;         y=p->y-minYDomain;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<4; l++)  aveP+=h[l]*w[l];
           p->p1=aveP/weight;
           p=p->next;
           cnt++;
         }

       }	//End of i,j
     free(lambda);
     free(w);
     free(gridP);
     free(h);
     break;
   }
}


/*
void findAveragePy(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s)
{
   int i,j,k,l,istart,iend,jstart,jend,kstart,kend,indexI,indexJ,indexK,cnt;
   double weight,x,y,z,wxl,wxr,wyl,wyr,aveP;
   double *lambda,*w,*h,*gridP;
   Particle ***particle;
   ptclList *p;
   particle=D->particle;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->jstart;
   kend=D->jend;
   weight=1.0/((double)D->numberInCell)/D->resolX/D->resolY;

   switch (D->dimension) {
   //2D
   case 2 :
     lambda=(double *)malloc(16*sizeof(double ));
     w=(double *)malloc(4*sizeof(double ));
     gridP=(double *)malloc(4*sizeof(double ));
     h=(double *)malloc(4*sizeof(double ));

     k=0;
     for(i=istart; i<iend; i++)  
       for(j=jstart; j<jend; j++)  
       {
         for(l=0; l<16; l++)  lambda[l]=0.0;

         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           x=p->x;         y=p->y;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<16; l++)  {
             indexI=l%4;
             indexJ=l/4;
             lambda[l]+=w[indexI]*w[indexJ];
           }
           cnt++;
           p=p->next;
         }
         gridP[0]=den1[i][j][k];
         gridP[1]=den2[i][j][k];
         gridP[2]=den3[i][j][k];
         gridP[3]=den4[i][j][k];
         findFactor(h,lambda,gridP,4,4);

         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           aveP=0.0;
           x=p->x;         y=p->y;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<4; l++)  aveP+=h[0]*w[0];
           p->p2=aveP/weight;
           p=p->next;
         }


       }	//End of i,j

     free(lambda);
     free(w);
     free(gridP);
     free(h);
     break;
   }
}

void findAveragePz(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s)
{
   int i,j,k,l,istart,iend,jstart,jend,kstart,kend,indexI,indexJ,indexK,cnt;
   double weight,x,y,z,wxl,wxr,wyl,wyr,aveP;
   double *lambda,*w,*h,*gridP;
   Particle ***particle;
   ptclList *p;
   particle=D->particle;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->jstart;
   kend=D->jend;
   weight=1.0/((double)D->numberInCell)/D->resolX/D->resolY;

   switch (D->dimension) {
   //2D
   case 2 :
     lambda=(double *)malloc(16*sizeof(double ));
     w=(double *)malloc(4*sizeof(double ));
     gridP=(double *)malloc(4*sizeof(double ));
     h=(double *)malloc(4*sizeof(double ));

     k=0;
     for(i=istart; i<iend; i++)  
       for(j=jstart; j<jend; j++)  
       {
         for(l=0; l<16; l++)  lambda[l]=0.0;

         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           x=p->x;         y=p->y;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<16; l++)  {
             indexI=l%4;
             indexJ=l/4;
             lambda[l]+=w[indexI]*w[indexJ];
           }
           cnt++;
           p=p->next;
         }
         gridP[0]=den1[i][j][k];
         gridP[1]=den2[i][j][k];
         gridP[2]=den3[i][j][k];
         gridP[3]=den4[i][j][k];
         findFactor(h,lambda,gridP,4,4);

         p=particle[i][j][k].head[s]->pt;
         while(p)   {
           aveP=0.0;
           x=p->x;         y=p->y;
           x=x-(int)x;     y=y-(int)y;
           wxl=1.0-x;   wxr=1.0-wxl;
           wyl=1.0-y;   wyr=1.0-wyl;
           w[0]=wxl*wyl;
           w[1]=wxr*wyl;
           w[2]=wxl*wyr;
           w[3]=wxr*wyr;
           for(l=0; l<4; l++)  aveP+=h[0]*w[0];
           p->p3=aveP/weight;
           p=p->next;
         }


       }	//End of i,j

     free(lambda);
     free(w);
     free(gridP);
     free(h);
     break;
   }
}
*/
