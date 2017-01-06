#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreDoubleData(char *fileName,char *dataName,double *data,int totalCnt,int dataNum,int offSet);
void restoreIntData(char *fileName,char *dataName,int *data,int totalCnt,int dataNum,int offSet);

void restoreDump(Domain D,int iteration)
{
  char name[100];
  void restoreDumpHDF(Domain *D,int iteration);

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((Pukhov-1)*3+2) :
    if(D.saveDumpMode==HDF)  {
      restoreDumpHDF(&D,iteration);
    }
    break;
  }

}


void restoreDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,n,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int indexI,indexJ,indexK;
    int nxSub,nySub,nzSub,nSpecies,totalCnt,start,core,dataCnt;
    int minXDomain,minYDomain,minZDomain,L,M,N,flag,shareCnt;
    int startIndexX,startIndexY,startIndexZ,unitX,unitY,unitZ;
    int rankX,rankY,rankZ,tmp,shift,rank,cntSub,remain,sub;
    int *recv,*minXSub,*maxXSub,*minYSub,*maxYSub,*minZSub,*maxZSub;
    int *sharePNum,*recvDataCnt,*dataCore,*coreCnt,*dataIndex,*dataCores;
    double *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz;
    double **sendData,**recvData,*shareData;
    double x,y,z,px,py,pz;
    int offset[3];
    char name[100],name2[100];
    void restoreIntMeta();
    void restoreFieldComp();
    ptclList *p;
    Particle ***particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,filespace;
    hsize_t metaDim[1];
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
    L=D->L;  M=D->M;  N=D->N;

    recv=(int *)malloc(nTasks*sizeof(int ));
    minXSub=(int *)malloc(nTasks*sizeof(int ));
    maxXSub=(int *)malloc(nTasks*sizeof(int ));
    minYSub=(int *)malloc(nTasks*sizeof(int ));
    maxYSub=(int *)malloc(nTasks*sizeof(int ));
    minZSub=(int *)malloc(nTasks*sizeof(int ));
    maxZSub=(int *)malloc(nTasks*sizeof(int ));
    sharePNum=(int *)malloc(nTasks*sizeof(int ));
    recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    coreCnt=(int *)malloc(nTasks*sizeof(int ));

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(name,"/minXDomain",&D->minXDomain,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      ny=D->ny+5;      
      nx=D->nx+5;     
      istart=0;
      iend+=3; 
      jstart=0;
      jend+=3;
      nxSub+=5;
      nySub+=5;

      offset[0]=D->minXSub;
      offset[1]=D->minYSub-D->minYDomain;
      offset[2]=0;
      
      restoreFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      restoreFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

      //restore minSub, maxSub
      unitX=D->nx/D->L+1;
      unitY=D->ny/D->M+1;
      MPI_Gather(&D->minXSub,1,MPI_INT,minXSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&D->maxXSub,1,MPI_INT,maxXSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&D->minYSub,1,MPI_INT,minYSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&D->maxYSub,1,MPI_INT,maxYSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&D->minZSub,1,MPI_INT,minZSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minZSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&D->maxZSub,1,MPI_INT,maxZSub,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxZSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);

      //restore particle
      sprintf(name,"dumpParticle%d.h5",iteration);
      if(myrank==0)  {
        restoreIntMeta(name,"/nSpecies",&nSpecies,1);
      }    else	;
      MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

      for(s=0; s<nSpecies; s++)
      {
        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvData[i]=0;
          coreCnt[i]=0;
        }

        sprintf(name2,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(name,name2,&totalCnt,1);
        else	;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmp=sub+1;
          else		   tmp=sub;
          if(myrank==rank)
            cntSub=tmp;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        dataX = (double *)malloc(cntSub*sizeof(double ));
        dataY = (double *)malloc(cntSub*sizeof(double ));
        dataPx = (double *)malloc(cntSub*sizeof(double ));
        dataPy = (double *)malloc(cntSub*sizeof(double ));
        dataPz = (double *)malloc(cntSub*sizeof(double ));
        dataIndex = (int *)malloc(cntSub*sizeof(int ));
        dataCore = (int *)malloc(cntSub*sizeof(int ));
        dataCores = (int *)malloc(cntSub*sizeof(int ));

        sprintf(name2,"%dParticle/x",s);
        restoreDoubleData(name,name2,dataX,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/y",s);
        restoreDoubleData(name,name2,dataY,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/px",s);
        restoreDoubleData(name,name2,dataPx,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/py",s);
        restoreDoubleData(name,name2,dataPy,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/pz",s);
        restoreDoubleData(name,name2,dataPz,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/index",s);
        restoreIntData(name,name2,dataIndex,totalCnt,cntSub,start);
        sprintf(name2,"%dParticle/core",s);
        restoreIntData(name,name2,dataCores,totalCnt,cntSub,start);

        for(i=0; i<cntSub; i++)
          dataCore[i]=-1;

        for(i=0; i<cntSub; i++)  {
          x=dataX[i]-D->minXDomain;         
          y=dataY[i];          
          startIndexX=((int)x)/unitX;
          startIndexY=((int)(y-D->minYDomain))/unitY;
          startIndexZ=0;
          rank=startIndexY+startIndexZ*M+startIndexX*M*N;
          flag=0;
          while(flag==0)  {
            if(minXSub[rank]<=x && x<maxXSub[rank] &&
               minYSub[rank]<=y && y<maxYSub[rank]) { 
              flag=1;           
              if(rank==myrank)  {
                indexI=(int)(x-D->minXSub)+2;
                indexJ=(int)(y-D->minYSub)+2;
                p = (ptclList *)malloc(sizeof(ptclList));
                p->next = particle[indexI][indexJ][0].head[s]->pt;
                particle[indexI][indexJ][0].head[s]->pt=p;
            
                p->x=(x-D->minXSub)-(int)(x-D->minXSub);
                p->y=(y-D->minYSub)-(int)(y-D->minYSub);
                p->z=0.0;
                p->p1=dataPx[i];
                p->p2=dataPy[i];
                p->p3=dataPz[i];
                p->index=dataIndex[i];
                p->core=dataCores[i];
              }
              dataCore[i]=rank;
              sharePNum[rank]+=1;
            } 
            else 	rank++	;
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //set count '0' at myrank due to no reason for sharing
        for(i=0; i<nTasks; i++)
          if(myrank==i)  sharePNum[i]=0;

        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    
            MPI_Send(&sharePNum[i],1,MPI_INT,i,myrank,MPI_COMM_WORLD);         
        }
        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    {
            MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&mpi_status);
          }  else	;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //memory for send and recving data
        dataCnt=7;
        for(i=0; i<nTasks; i++)   {          
          sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
          recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
        }        

        for(i=0; i<cntSub; i++)  
        {          
          core=dataCore[i];
          x=dataX[i]-D->minXDomain;          
          y=dataY[i];
          if(myrank!=core)  {
            n=coreCnt[core];
            sendData[core][n*dataCnt+0]=dataX[i]-D->minXDomain;
            sendData[core][n*dataCnt+1]=dataY[i];
            sendData[core][n*dataCnt+2]=dataPx[i];
            sendData[core][n*dataCnt+3]=dataPy[i];
            sendData[core][n*dataCnt+4]=dataPz[i];
            sendData[core][n*dataCnt+5]=dataIndex[i];
            sendData[core][n*dataCnt+6]=dataCores[i];
            coreCnt[core]+=1;
          }	else	;
        }
//lala
        for(i=0; i<nTasks; i++)
        {
          if(myrank==i)  {
            for(j=0; j<nTasks; j++)
              if(i!=j)
                MPI_Send(sendData[j],sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
          } 
          else  {
            MPI_Recv(recvData[i],recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&mpi_status);
            for(j=0; j<recvDataCnt[i]; j++)  {
              x=recvData[i][j*dataCnt+0];
              y=recvData[i][j*dataCnt+1];
              indexI=(int)(x-D->minXSub)+2;
              indexJ=(int)(y-D->minYSub)+2;
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ][0].head[s]->pt;
              particle[indexI][indexJ][0].head[s]->pt=p;
            
              p->x=(x-D->minXSub)-(int)(x-D->minXSub);
              p->y=(y-D->minYSub)-(int)(y-D->minYSub);
              p->z=0.0;
              p->p1=recvData[i][j*dataCnt+2];
              p->p2=recvData[i][j*dataCnt+3];
              p->p3=recvData[i][j*dataCnt+4];
              p->index=recvData[i][j*dataCnt+5];
              p->core=recvData[i][j*dataCnt+6];
            }
          } 
          MPI_Barrier(MPI_COMM_WORLD);
        }

        free(dataX);
        free(dataY);
        free(dataPx);
        free(dataPy);
        free(dataPz);
        free(dataIndex);
        free(dataCore);
        free(dataCores);
        for(i=0; i<nTasks; i++)   {          
          free(recvData[i]);
          free(sendData[i]);
        }        

      }			//End of for(nSpecies)

    break;
/*
    //3D
    case (Split-1)*3+3:
      ny=D->ny+5;      
      nx=D->nx+5;     
      nz=D->nz+5;      
      istart=0;
      iend+=3; 
      jstart=0;
      jend+=3;
      kstart=0;
      kend+=3;
      nxSub+=5;
      nySub+=5;
      nzSub+=5;
      offset[1]=0;
      offset[0]=D->minYSub-D->minYDomain;
      offset[2]=D->minZSub-D->minZDomain;
      
      //save metaData for field
      sprintf(name,"dump%d.h5",iteration);
      if(myrank==0)
        restoreIntMeta(name,"/minXSub",&(D->minXSub));
      else	;
      MPI_Bcast(&(D->minXSub),1,MPI_INT,0,MPI_COMM_WORLD);
//here
      restoreField3D(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->ExC,name,"/ExC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->PrC,name,"/PrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->PlC,name,"/PlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->BxC,name,"/BxC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->SrC,name,"/SrC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->SlC,name,"/SlC",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JxOld,name,"/JxOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JyOld,name,"/JyOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreField3D(D->JzOld,name,"/JzOld",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      //restore metaData for particle
      if(myrank==0)
        restoreIntMeta(name,"/nSpecies",&nSpecies);
      else	;
      MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

      //restore particle
      for(s=0; s<nSpecies; s++)
      {
        sprintf(name2,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(name,name2,&totalCnt);
        else	;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/offSet",s);
        if(myrank==0)
          restoreIntArray(name,name2,offSet,nTasks);
        else	;
        MPI_Bcast(offSet,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        sprintf(name2,"%dParticle/counts",s);
        if(myrank==0)
          restoreIntArray(name,name2,counts,nTasks);
        else	;
        MPI_Bcast(counts,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        if(counts[myrank]>0)
        {
          dataI = (int *)malloc(counts[myrank]*sizeof(int ));
          dataJ = (int *)malloc(counts[myrank]*sizeof(int ));
          dataK = (int *)malloc(counts[myrank]*sizeof(int ));
          dataCore = (int *)malloc(counts[myrank]*sizeof(int ));
          dataIndex = (int *)malloc(counts[myrank]*sizeof(int ));
          dataX = (double *)malloc(counts[myrank]*sizeof(double ));
          dataY = (double *)malloc(counts[myrank]*sizeof(double ));
          dataZ = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPx = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPy = (double *)malloc(counts[myrank]*sizeof(double ));
          dataPz = (double *)malloc(counts[myrank]*sizeof(double ));
          sprintf(name2,"%dParticle/i",s);
          restoreIntData(name,name2,dataI,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/j",s);
          restoreIntData(name,name2,dataJ,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/k",s);
          restoreIntData(name,name2,dataK,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/core",s);
          restoreIntData(name,name2,dataCore,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/index",s);
          restoreIntData(name,name2,dataIndex,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/x",s);
          restoreDoubleData(name,name2,dataX,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/y",s);
          restoreDoubleData(name,name2,dataY,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/z",s);
          restoreDoubleData(name,name2,dataZ,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/px",s);
          restoreDoubleData(name,name2,dataPx,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/py",s);
          restoreDoubleData(name,name2,dataPy,totalCnt,counts[myrank],offSet[myrank]);
          sprintf(name2,"%dParticle/pz",s);
          restoreDoubleData(name,name2,dataPz,totalCnt,counts[myrank],offSet[myrank]);

          for(n=0; n<counts[myrank]; n++)
          {
            i=dataI[n];
            j=dataJ[n];
            k=dataK[n];
            p = (ptclList *)malloc(sizeof(ptclList));
            p->next = particle[i][j][k].head[s]->pt;
            particle[i][j][k].head[s]->pt=p;
            
            p->x=dataX[n];
            p->y=dataY[n];
            p->z=dataZ[n];
            p->p1=dataPx[n];
            p->p2=dataPy[n];
            p->p3=dataPz[n];
            p->index=dataIndex[n];
            p->core=dataCore[n];
          }

          free(dataI);
          free(dataJ);
          free(dataK);
          free(dataCore);
          free(dataIndex);
          free(dataX);
          free(dataY);
          free(dataZ);
          free(dataPx);
          free(dataPy);
          free(dataPz);
        }
      }			//End of for(nSpecies)
    break;
*/
    free(recv);
    free(minXSub);
    free(maxXSub);
    free(minYSub);
    free(maxYSub);
    free(minZSub);
    free(maxZSub);
    free(sharePNum);
    free(recvDataCnt);
    free(coreCnt);
    free(recvData);
    free(sendData);
  }		//End of switch(dimension....)
}

void restoreIntData(char *fileName,char *dataName,int *data,int totalCnt,int dataNum,int offSet)
{
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=totalCnt;     
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=dataNum;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dread(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreDoubleData(char *fileName,char *dataName,double *data,int totalCnt,int dataNum,int offSet)
{
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,group_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=totalCnt;     
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=dataNum;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreField2D(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int *offSet)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=nx;      
  dimsf[1]=ny;     
  filespace=H5Screate_simple(2,dimsf,NULL);

  count[1]=nySub;
  count[0]=nxSub;
  offset[1]=offSet[1];
  offset[0]=offSet[0];
  memspace=H5Screate_simple(2,count,NULL);

  field = (double *)malloc(nxSub*nySub*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(i=istart; i<iend; i++)
  {
    for(j=jstart; j<jend; j++)
      data[i][j][0]=field[start+j-jstart];
    start+=nySub;
  }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny;      
  dimsf[1]=nx;     
  dimsf[2]=nz;     
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=nySub;
  count[1]=nxSub;
  count[2]=nzSub;
  offset[0]=offSet[1];
  offset[1]=offSet[0];
  offset[2]=offSet[2];
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(j=jstart; j<jend; j++)
    for(i=istart; i<iend; i++)
    {
      for(k=kstart; k<kend; k++)
        data[i][j][k]=field[start+k-kstart];
      start+=nzSub;
    }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

