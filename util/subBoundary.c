#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void calParameter(int nx,int *istart,int *iend,int *saveNxSub,int rankX,int *biasX,int reIend,int reNxSub,int L);


void boundary(Domain *D)
{
   int subX,subY,subZ,remainX,remainY,remainZ;
   int minX,minY,minZ,maxX,maxY,maxZ,rank,rankX,rankY,rankZ,tmpX,tmpY,tmpZ;
   char fileName[100];
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //set memory
   D->recv=(int *)malloc(nTasks*sizeof(int ));
   D->minXSubList=(int *)malloc(nTasks*sizeof(int ));
   D->maxXSubList=(int *)malloc(nTasks*sizeof(int ));
   D->minYSubList=(int *)malloc(nTasks*sizeof(int ));
   D->maxYSubList=(int *)malloc(nTasks*sizeof(int ));
   D->minZSubList=(int *)malloc(nTasks*sizeof(int ));
   D->maxZSubList=(int *)malloc(nTasks*sizeof(int ));
   D->recvData=(double **)malloc(nTasks*sizeof(double *));
   D->sendData=(double **)malloc(nTasks*sizeof(double *));
   D->recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
   D->sharePNum=(int *)malloc(nTasks*sizeof(int ));
   D->coreCnt=(int *)malloc(nTasks*sizeof(int ));
   D->reMinXSubList=(int *)malloc(nTasks*sizeof(int ));
   D->reMaxXSubList=(int *)malloc(nTasks*sizeof(int ));
   D->reMinYSubList=(int *)malloc(nTasks*sizeof(int ));
   D->reMaxYSubList=(int *)malloc(nTasks*sizeof(int ));
   D->reMinZSubList=(int *)malloc(nTasks*sizeof(int ));
   D->reMaxZSubList=(int *)malloc(nTasks*sizeof(int ));

   //restore meta data
   sprintf(fileName,"dumpParticle%d.h5",D->step);
   if(myrank==0)  {
     restoreIntMeta(fileName,"/nSpecies",&D->nSpecies);
     restoreIntMeta(fileName,"/nx",&D->nx);
     restoreIntMeta(fileName,"/ny",&D->ny);
     restoreIntMeta(fileName,"/nz",&D->nz);
     restoreIntMeta(fileName,"/minXDomain",&D->minXDomain);
     restoreIntMeta(fileName,"/minYDomain",&D->minYDomain);
     restoreIntMeta(fileName,"/minZDomain",&D->minZDomain);
   }    else   ;
   MPI_Bcast(&D->nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minYDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minZDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   if(D->minXDomain>0)  D->minXDomain-=1;  else ;

   D->nxSub=D->nx/D->L;
   subX=D->nxSub;
   remainX=D->nx%D->L;
   minX=maxX=0;
   D->nySub=D->ny/D->M;
   subY=D->nySub;
   remainY=D->ny%D->M;
   minY=maxY=0;
   D->nzSub=D->nz/D->N;
   subZ=D->nzSub;
   remainZ=D->nz%D->N;
   minZ=maxZ=0;

   minX=maxX=D->minXDomain;
   for(rankX=0; rankX<D->L; rankX++)
   {
     if(rankX<remainX)   tmpX=subX+1;
     else                tmpX=subX;
     minX=maxX;
     maxX=minX+tmpX;

     minZ=maxZ=D->minZDomain;
     for(rankZ=0; rankZ<D->N; rankZ++)
     {
       if(rankZ<remainZ)   tmpZ=subZ+1;
       else                tmpZ=subZ;
       minZ=maxZ;
       maxZ=minZ+tmpZ;

       minY=maxY=D->minYDomain;
       for(rankY=0; rankY<D->M; rankY++)
       {
         if(rankY<remainY)   tmpY=subY+1;
         else                tmpY=subY;
         minY=maxY;
         maxY=minY+tmpY;

         rank=rankY+rankZ*D->M+rankX*(D->M*D->N);
         if(myrank==rank)
         {
           D->nxSub=tmpX;
           D->nySub=tmpY;
           D->nzSub=tmpZ;
           D->minXSub=minX;
           D->maxXSub=maxX;
           D->minYSub=minY;
           D->maxYSub=maxY;
           D->minZSub=minZ;
           D->maxZSub=maxZ;
         }
       }
     }
   }

   D->minX=D->minXDomain;
   D->maxX=(D->minXDomain+D->nx);
   D->minY=D->minYDomain;
   D->maxY=(D->minYDomain+D->ny);
   D->minZ=D->minZDomain;
   D->maxZ=(D->minZDomain+D->nz);

   //restore minSub, maxSub
   MPI_Gather(&D->minXSub,1,MPI_INT,D->minXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->minXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->maxXSub,1,MPI_INT,D->maxXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->maxXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->minYSub,1,MPI_INT,D->minYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->minYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->maxYSub,1,MPI_INT,D->maxYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->maxYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->minZSub,1,MPI_INT,D->minZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->minZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->maxZSub,1,MPI_INT,D->maxZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->maxZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);

   D->istart=2;
   D->iend=D->nxSub+2;
   D->jstart=0;
   D->jend=1;
   D->kstart=0;
   D->kend=1;
   if(D->dimension>1)  {
     D->jstart=2;
     D->jend=D->nySub+2;
   }
   if(D->dimension>2)  {
     D->kstart=2;
     D->kend=D->nzSub+2;
   }

   D->rankX=myrank/(D->M*D->N);
   D->rankZ=(myrank%(D->M*D->N))/D->M;
   D->rankY=(myrank%(D->M*D->N))%D->M;

   D->saveJstart=0;
   D->saveKstart=0;
   D->saveJend=1;
   D->saveKend=1;
   D->saveNySub=1;
   D->saveNzSub=1;
   calParameter(D->nx+5,&D->saveIstart,&D->saveIend,&D->saveNxSub,D->rankX,&D->biasX,D->iend,D->nxSub,D->L);
   if(D->dimension>1)
     calParameter(D->ny+5,&D->saveJstart,&D->saveJend,&D->saveNySub,D->rankY,&D->biasY,D->jend,D->nySub,D->M);
   if(D->dimension>2)
     calParameter(D->nz+5,&D->saveKstart,&D->saveKend,&D->saveNzSub,D->rankZ,&D->biasZ,D->kend,D->nzSub,D->N);
   MPI_Barrier(MPI_COMM_WORLD);

}

void resolHighBoundary(Domain *D)
{
   int i,j,k,s,nxSub1D,nySub2D,nzSub3D;
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   D->reNx=D->nx*D->resolX;
   D->reNy=D->ny*D->resolY;
   D->reNz=D->nz*D->resolZ;

   D->reNxSub=D->nxSub*D->resolX;
   D->reNySub=D->nySub*D->resolY;
   D->reNzSub=D->nzSub*D->resolZ;
   D->reIend=D->reNxSub+2;
   D->reJend=D->reKend=1;
   if(D->dimension>1)
     D->reJend=D->reNySub+2;
   if(D->dimension>2)
     D->reKend=D->reNzSub+2;

   //set domain range
   D->reMinX=D->minXDomain*D->resolX;
   D->reMaxX=(D->minXDomain+D->nx)*D->resolX;
   D->reMinY=D->minYDomain*D->resolY;
   D->reMaxY=(D->minYDomain+D->ny)*D->resolY;
   D->reMinZ=D->minZDomain*D->resolZ;
   D->reMaxZ=(D->minZDomain+D->nz)*D->resolZ;
   D->reMinXSub=D->minXSub*D->resolX;
   D->reMinYSub=D->minYSub*D->resolY;
   D->reMinZSub=D->minZSub*D->resolZ;
   D->reMaxXSub=D->maxXSub*D->resolX;
   D->reMaxYSub=D->maxYSub*D->resolY;
   D->reMaxZSub=D->maxZSub*D->resolZ;

   D->reSaveJstart=0;
   D->reSaveKstart=0;
   D->reSaveJend=1;
   D->reSaveKend=1;
   D->reSaveNySub=1;
   D->reSaveNzSub=1;
   calParameter(D->reNx+5,&D->reSaveIstart,&D->reSaveIend,&D->reSaveNxSub,D->rankX,&D->biasX,D->reIend,D->reNxSub,D->L);
   if(D->dimension>1)
     calParameter(D->reNy+5,&D->reSaveJstart,&D->reSaveJend,&D->reSaveNySub,D->rankY,&D->biasY,D->reJend,D->reNySub,D->M);
   if(D->dimension>2)
     calParameter(D->reNz+5,&D->reSaveKstart,&D->reSaveKend,&D->reSaveNzSub,D->rankZ,&D->biasZ,D->reKend,D->reNzSub,D->N);
   MPI_Barrier(MPI_COMM_WORLD);

   //restore reminSub, remaxSub
   MPI_Gather(&D->reMinXSub,1,MPI_INT,D->reMinXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMinXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->reMaxXSub,1,MPI_INT,D->reMaxXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMaxXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->reMinYSub,1,MPI_INT,D->reMinYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMinYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->reMaxYSub,1,MPI_INT,D->reMaxYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMaxYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->reMinZSub,1,MPI_INT,D->reMinZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMinZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Gather(&D->reMaxZSub,1,MPI_INT,D->reMaxZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->reMaxZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);

   //particle memory
   nxSub1D=D->reNxSub+3;
   nySub2D=1;
   nzSub3D=1;
   if(D->dimension>1) {
     nySub2D=D->reNySub+3;
   }
   if(D->dimension>2) {
     nzSub3D=D->reNzSub+3;
   }

   D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
   for(i=0; i<nxSub1D; i++) {
     D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
     for(j=0; j<nySub2D; j++)
       D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
   }
   if(D->dimension==2)
   {
     k=0;
     for(i=0; i<D->reIend+1; i++)  
       for(j=D->jstart-1; j<D->reJend+1; j++)
       {
         D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)
         {
           D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j][k].head[s]->pt = NULL;
         }
       }
   }


}

void calParameter(int nx,int *istart,int *iend,int *saveNxSub,int rankX,int *biasX,int reIend,int reNxSub,int L)
{
  if(L==1) {
    *istart=0;
    *iend=reIend+3;
    *biasX=0;
    *saveNxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0;
      *iend=reIend;
      *saveNxSub=reNxSub+2;
      *biasX=0;
    }  else if(rankX==L-1)  {
      *istart=2;
      *iend=reIend+3;
      *saveNxSub=reNxSub+3;
      *biasX=2;
    } else  {
      *istart=2;
      *iend=reIend;
      *saveNxSub=reNxSub;
      *biasX=2;
    }
  }
}



double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)
   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0.0;
       }
   return field;
}
  
int ***memoryAsignInt(int nx, int ny, int nz)
{
   int i,j,k;
   int ***field;

   field = (int ***)malloc((nx)*sizeof(int **));
   for(i=0; i<nx; i++)
   {
     field[i] = (int **)malloc((ny)*sizeof(int *));
     for(j=0; j<ny; j++)
       field[i][j] = (int *)malloc((nz)*sizeof(int ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0;
       }
   return field;
}

