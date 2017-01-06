#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void saveParticleHDF(Domain *D,int iteration,int s,double minPx)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub,cnt,totalCnt,start,index;
    int minXSub,minYSub,minZSub;
    double dx,dy,dz,lambda,tmpDouble;
    char name[100];
    double *saveDouble;
    int *saveInt,offset[2];
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;
    int myrank, nTasks;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    int recv[nTasks];

    hid_t file_id;
    herr_t status;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    lambda=D->lambda;
    minXSub=D->minXSub;
    minYSub=D->minYSub;
    minZSub=D->minZSub;


    sprintf(name,"%dParticle%d.h5",s,iteration);
    if(myrank==0)    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    switch(D->dimension) {
    //2D
    case 2:
      k=0;
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)    {
            if(p->p1>=minPx)  cnt++;
            else	;
            p=p->next;
          }
        }
      saveDouble = (double *)malloc(cnt*sizeof(double ));      
      saveInt = (int *)malloc(cnt*sizeof(int ));      
      MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

      start=0;
      for(i=0; i<myrank; i++)
        start+=recv[i];
      totalCnt=0;
      for(i=0; i<nTasks; i++)
        totalCnt+=recv[i];

      if(myrank==0)
        saveIntMeta(name,"totalCnt",&totalCnt,1);
      else 	;
      MPI_Barrier(MPI_COMM_WORLD);

      if(totalCnt>0)
      {
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              if(p->p1>=minPx)   {
                tmpDouble=((i-istart+minXSub)+p->x)*dx*lambda;
                saveDouble[index]=tmpDouble;
                index++;
              }      else 	;
              p=p->next;
            }
          }
        saveParticleComp_Double(saveDouble,name,"x",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)            {
              if(p->p1>=minPx)              {
                tmpDouble=((j-jstart+minYSub)+p->y)*dy*lambda;
                saveDouble[index]=tmpDouble;
                index++;
              }     else	;
              p=p->next;
            }
          } 
        saveParticleComp_Double(saveDouble,name,"y",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)         {
              if(p->p1>=minPx)              {
                saveDouble[index]=p->p1;
                index++;
              }      else	;
              p=p->next;
            }
          } 
        saveParticleComp_Double(saveDouble,name,"px",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)          {
              if(p->p1>=minPx)            {
                saveDouble[index]=p->p2;
                index++;
              }             else	;
              p=p->next;
            }
          } 
        saveParticleComp_Double(saveDouble,name,"py",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)          {
              if(p->p1>=minPx)              {
                saveDouble[index]=p->p3;
                index++;
              }            else	;
              p=p->next;
            }
          } 
        saveParticleComp_Double(saveDouble,name,"pz",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)          {
              if(p->p1>=minPx)            {
                saveInt[index]=p->index;
                index++;
              }            else	;
              p=p->next;
            }
          } 
        saveParticleComp_Int(saveInt,name,"index",totalCnt,cnt,start);

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)          {
              if(p->p1>=minPx)             {
                saveInt[index]=p->core;
                index++;
              }              else	;
              p=p->next;
            }
          } 
        saveParticleComp_Int(saveInt,name,"core",totalCnt,cnt,start);
      }	//End of totalCnt>0
      else	;

      free(saveDouble);
      free(saveInt);
      break;
/*
    //3D
    case 3:
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
                cnt++;
              else	;
              p=p->next;
            }
          }
      saveDouble = (double *)malloc(cnt*sizeof(double ));      
      saveInt = (int *)malloc(cnt*sizeof(int ));      
      MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      start=0;
      for(i=0; i<myrank; i++)
        start+=recv[i];
      totalCnt=0;
      for(i=0; i<nTasks; i++)
        totalCnt+=recv[i];

      if(myrank==0)
        saveIntMeta(name,"totalCnt",&totalCnt);
      else 	;
      MPI_Barrier(MPI_COMM_WORLD);

      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                tmpDouble=((i-istart+minXSub)+p->x)*dx*lambda;
                saveDouble[index]=tmpDouble;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"x",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                tmpDouble=((j-jstart+minYSub)+p->y)*dy*lambda;
                saveDouble[index]=tmpDouble;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"y",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                tmpDouble=((k-kstart+minZSub)+p->z)*dz*lambda;
                saveDouble[index]=tmpDouble;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"z",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                saveDouble[index]=p->p1;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"px",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                saveDouble[index]=p->p2;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"py",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                saveDouble[index]=p->p3;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Double(saveDouble,name,"pz",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                saveInt[index]=p->index;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Int(saveInt,name,"index",totalCnt,cnt,start);
      index=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(p->p1>=minPx)
              {
                saveInt[index]=p->core;
                index++;
              }
              else	;
              p=p->next;
            }
          } 
      MPI_Barrier(MPI_COMM_WORLD);
      saveParticleComp_Int(saveInt,name,"core",totalCnt,cnt,start);

      free(saveDouble);
      free(saveInt);
      break;
*/
    }		//End of switch(dimension....)
}
