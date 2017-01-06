#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <math.h>
#include <mpi.h>

void cleanSubBoundary(Domain *D)
{
   int i,j,k,s;
   Particle ***particle;
   particle=D->particle;

   if(D->dimension==2)  {
     k=0;
     for(i=0; i<D->reIend+1; i++)
       for(j=1; j<D->reJend+1; j++)  {
         for(s=0; s<D->nSpecies; s++)
           free(particle[i][j][k].head[s]);
         free(particle[i][j][k].head);
       }
     for(i=0; i<D->reNxSub+3; i++)  {
       for(j=1; j<D->reNySub+3; j++)  
         free(D->particle[i][j]);
       free(D->particle[i]);
     }
     free(particle);
   }
}

void deleteParticleMemory(Domain *D,int s)
{
   int i,j,k;
   Particle ***particle;
   particle=D->particle;
   ptclList *p,*tmp;

   //remove particles
   if(D->dimension==2)
   {
     k=0;
     for(i=0; i<D->reIend+1; i++)
       for(j=1; j<D->reJend+1; j++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           tmp=p->next;
           particle[i][j][k].head[s]->pt=tmp;
           p->next=NULL;
           free(p);
           p=particle[i][j][k].head[s]->pt;
         }
       }
   }	//End of dimension==2


}

void cleanBoundary(Domain *D)
{
   free(D->recv);
   free(D->minXSubList);
   free(D->maxXSubList);
   free(D->minYSubList);
   free(D->maxYSubList);
   free(D->minZSubList);
   free(D->maxZSubList);
   free(D->recvData);
   free(D->sendData);
   free(D->recvDataCnt);
   free(D->coreCnt);
   free(D->reMinXSubList);
   free(D->reMaxXSubList);
   free(D->reMinYSubList);
   free(D->reMaxYSubList);
   free(D->reMinZSubList);
   free(D->reMaxZSubList);
}

void deleteField(double ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)  
   {
     for(j=0; j<ny; j++)  
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}

void deleteFieldInt(int ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)
   {
     for(j=0; j<ny; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


