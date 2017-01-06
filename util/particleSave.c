#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreDoubleData(char *fileName,char *dataName,double *data,int totalCnt,int dataNum,int offSet);
void restoreIntData(char *fileName,char *dataName,int *data,int totalCnt,int dataNum,int offSet);

void highParticleSave(Domain *D,char *fileName,int s)
{
   int i,j,k,n,sub,remain,tmpInt,rank,start,flag,unitX,unitY,unitZ;
   int startIndexX,startIndexY,startIndexZ,indexI,indexJ,indexK;
   int core,dataCnt;
   double x,y,z,xx,yy,zz,dx,dy;
   int reStart,*testList,index,outCnt;
   char dataName[100],name[100];
   FILE *out;
   Particle ***particle;
   ptclList *p;
   particle=D->particle;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   for(i=0; i<nTasks; i++)  {
     D->sharePNum[i]=0;
     D->recvData[i]=0;
     D->coreCnt[i]=0;
     D->recvDataCnt[i]=0;
   }
   
   sprintf(dataName,"%dParticle/totalCnt",s);
   if(myrank==0) 
     restoreIntMeta(fileName,dataName,&D->totalCnt);
   else    ;
   MPI_Bcast(&D->totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

   sub=D->totalCnt/nTasks;
   remain=D->totalCnt%nTasks;
   for(rank=0; rank<nTasks; rank++) {
     if(rank<remain)  tmpInt=sub+1;
     else             tmpInt=sub;
     if(myrank==rank)
       D->cntSub=tmpInt;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Gather(&D->cntSub,1,MPI_INT,D->recv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(D->recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

   D->dataX = (double *)malloc(D->cntSub*sizeof(double ));
   D->dataY = (double *)malloc(D->cntSub*sizeof(double ));
   D->dataIndex = (int *)malloc(D->cntSub*sizeof(int ));
   D->dataCores = (int *)malloc(D->cntSub*sizeof(int ));

   start=0;
   for(i=0; i<myrank; i++)
     start+=D->recv[i];

   if(D->totalCnt>0)
   {
     sprintf(dataName,"%dParticle/x",s);
     restoreDoubleData(fileName,dataName,D->dataX,D->totalCnt,D->cntSub,start);
     sprintf(dataName,"%dParticle/y",s);
     restoreDoubleData(fileName,dataName,D->dataY,D->totalCnt,D->cntSub,start);
//     sprintf(dataName,"%dParticle/px",s);
//     restoreDoubleData(fileName,dataName,D->dataPx,D->totalCnt,D->cntSub,start);
//     sprintf(dataName,"%dParticle/py",s);
//     restoreDoubleData(fileName,dataName,D->dataPy,D->totalCnt,D->cntSub,start);
//     sprintf(dataName,"%dParticle/pz",s);
//     restoreDoubleData(fileName,dataName,D->dataPz,D->totalCnt,D->cntSub,start);
     sprintf(dataName,"%dParticle/index",s);
     restoreIntData(fileName,dataName,D->dataIndex,D->totalCnt,D->cntSub,start);
     sprintf(dataName,"%dParticle/core",s);
     restoreIntData(fileName,dataName,D->dataCores,D->totalCnt,D->cntSub,start);

     dx=1.0/(double)(D->resolX);
     dy=1.0/(double)(D->resolY);
     n=D->resolX*D->resolY;
     D->reCntSub=n*D->cntSub;
     D->reTotalCnt=n*D->totalCnt;
     reStart=n*start;
     testList = (int *)malloc(D->reCntSub*sizeof(int ));
     outCnt=0;
     for(i=0; i<D->cntSub; i++)   {
       x=D->dataX[i];
       y=D->dataY[i];
       for(j=0; j<D->resolY; j++)
         for(k=0; k<D->resolX; k++) {
           xx=(x-0.5+dx*(0.5+k));
           yy=(y-0.5+dy*(0.5+j));
           if(D->minX<D->minX || xx>=D->maxX || yy<D->minY || yy>=D->maxY)  {
             testList[i*n+j*D->resolX+k]=1;
             outCnt++;
           }
           else   testList[i*n+j*D->resolX+k]=0;
         }
     }

     //recalculation of saving memory
     D->reCntSub=D->cntSub*n-outCnt;
     MPI_Gather(&D->reCntSub,1,MPI_INT,D->recv,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
     D->reTotalCnt=0;
     for(i=0; i<nTasks; i++)
       D->reTotalCnt+=D->recv[i];
     reStart=0;
     for(i=0; i<myrank; i++)
       reStart+=D->recv[i];
     MPI_Barrier(MPI_COMM_WORLD);

     D->reDataX = (double *)malloc(D->reCntSub*sizeof(double ));
     D->reDataY = (double *)malloc(D->reCntSub*sizeof(double ));
     D->reDataIndex = (int *)malloc(D->reCntSub*sizeof(int ));
     D->reDataCores = (int *)malloc(D->reCntSub*sizeof(int ));
     D->reDataCore = (int *)malloc(D->reCntSub*sizeof(int ));

     index=0;
     for(i=0; i<D->cntSub; i++)   {
       x=D->dataX[i];
       y=D->dataY[i];
       for(j=0; j<D->resolY; j++)
         for(k=0; k<D->resolX; k++)
           if(testList[i*n+j*D->resolX+k]==0)  {
             D->reDataX[index]=(x-0.5+dx*(0.5+k));
             D->reDataY[index]=(y-0.5+dy*(0.5+j));
             D->reDataIndex[index]=D->dataIndex[i];
             D->reDataCores[index]=D->dataCores[i];
             index++;
           }  else   ;
     }
     
     free(D->dataX);
     free(D->dataY);
     free(D->dataIndex);
     free(D->dataCores);
     free(testList);

     for(i=0; i<D->reCntSub; i++)
       D->reDataCore[i]=-1;

     unitX=D->nx/D->L+1;
     unitY=D->ny/D->M+1;
     for(i=0; i<D->reCntSub; i++)  {
       x=D->reDataX[i];         
       y=D->reDataY[i];          
       startIndexX=((int)(x-D->minXDomain))/unitX;
       startIndexY=((int)(y-D->minYDomain))/unitY;
       startIndexZ=0;
       rank=startIndexY+startIndexZ*D->M+startIndexX*D->M*D->N;
       flag=0;
       while(flag==0)  {
         if(D->minXSubList[rank]<=x && x<D->maxXSubList[rank] &&
            D->minYSubList[rank]<=y && y<D->maxYSubList[rank]) { 
           flag=1;           
           if(rank==myrank)  {
             indexI=(int)(x-D->minXSub)+2;
             indexJ=(int)(y-D->minYSub)+2;
             if(indexI>=D->istart && indexI<D->iend  &&
                indexJ>=D->jstart && indexJ<D->jend)  {
               p = (ptclList *)malloc(sizeof(ptclList));
               p->next = particle[indexI][indexJ][0].head[s]->pt;
               particle[indexI][indexJ][0].head[s]->pt=p;
            
               p->x=D->reDataX[i];
               p->y=D->reDataY[i];
               p->z=0.0;
               p->index=D->reDataIndex[i];
               p->core=D->reDataCores[i];
             } else        ;

           }  else ;
           D->reDataCore[i]=rank;
           D->sharePNum[rank]+=1;
         } 
         else      rank++  ;
       }
     }
     MPI_Barrier(MPI_COMM_WORLD);

     //set count '0' at myrank due to no reason for sharing
     for(i=0; i<nTasks; i++)
       if(myrank==i)  D->sharePNum[i]=0;
  
     for(i=0; i<nTasks; i++)   {          
       if(myrank!=i)    
         MPI_Send(&D->sharePNum[i],1,MPI_INT,i,myrank,MPI_COMM_WORLD);            }
     for(i=0; i<nTasks; i++)   {          
       if(myrank!=i)    {
         MPI_Recv(&D->recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
       }  else     ;
     }
     MPI_Barrier(MPI_COMM_WORLD);

     //memory for send and recving data
     dataCnt=4;
     for(i=0; i<nTasks; i++)   {          
       D->sendData[i]=(double *)malloc(D->sharePNum[i]*dataCnt*sizeof(double ));
       D->recvData[i]=(double *)malloc(D->recvDataCnt[i]*dataCnt*sizeof(double ));
     }        
  
     for(i=0; i<D->reCntSub; i++)  
     {          
       core=D->reDataCore[i];
       if(myrank!=core)  {
         n=D->coreCnt[core];
         D->sendData[core][n*dataCnt+0]=D->reDataX[i];
         D->sendData[core][n*dataCnt+1]=D->reDataY[i];
         D->sendData[core][n*dataCnt+2]=D->reDataIndex[i];
         D->sendData[core][n*dataCnt+3]=D->reDataCores[i];
         D->coreCnt[core]+=1;
       }   else    ;
     }

     for(i=0; i<nTasks; i++)
     {
       if(myrank==i)  {
         for(j=0; j<nTasks; j++)
           if(i!=j)
             MPI_Send(D->sendData[j],D->sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
       } 
       else  {
         MPI_Recv(D->recvData[i],D->recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
         for(j=0; j<D->recvDataCnt[i]; j++)  {
           x=D->recvData[i][j*dataCnt+0];
           y=D->recvData[i][j*dataCnt+1];
           indexI=(int)(x-D->minXSub)+2;
           indexJ=(int)(y-D->minYSub)+2;
           p = (ptclList *)malloc(sizeof(ptclList));
           p->next = particle[indexI][indexJ][0].head[s]->pt;
           particle[indexI][indexJ][0].head[s]->pt=p;
           p->x=x;
           p->y=y;
           p->z=0.0;
           p->index=D->recvData[i][j*dataCnt+2];
           p->core=D->recvData[i][j*dataCnt+3];
         }
       } 
       MPI_Barrier(MPI_COMM_WORLD);
     }

     free(D->reDataX);
     free(D->reDataY);
     free(D->reDataIndex);
     free(D->reDataCore);
     free(D->reDataCores);
     for(i=0; i<nTasks; i++)   {          
       free(D->sendData[i]);
       free(D->recvData[i]);
     }        

/*
sprintf(name,"particle%d",myrank);
out=fopen(name,"w");
k=0;
for(i=2; i<D->iend; i++)
  for(j=2; j<D->jend; j++)
  {
    for(s=0; s<D->nSpecies; s++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)  {
        fprintf(out,"%g %g\n",p->x,p->y);
        p=p->next;
      }
    }
  }
fclose(out);
*/



   }	//End of totalCnt>0

}

