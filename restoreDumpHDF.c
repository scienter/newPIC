#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

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
    double *data;
    double **sendData,**recvData,*shareData;
    double x,y,z,px,py,pz;
    int offset[3];
    char name[100],dataName[100];
    void restoreIntMeta();
    void restoreFieldComp();
    ptclList *p;
    Particle ***particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

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
    //1D
    case (Split-1)*3+1:
      nx=D->nx+5; nxSub+=5; istart=0; iend+=3; 
      ny=1; nySub=1; jstart=0; jend=1;
      nz=1; nzSub=1; kstart=0; kend=1;
      offset[0]=D->minXSub;
      offset[1]=0;
      offset[2]=0;
      
      restoreFieldComp(D->Ex,name,"/Ex",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Pr,name,"/Pr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Pl,name,"/Pl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Bx,name,"/Bx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Sr,name,"/Sr",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Sl,name,"/Sl",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Jx,name,"/Jx",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Jy,name,"/Jy",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      restoreFieldComp(D->Jz,name,"/Jz",nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      break;

    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5; nxSub+=5; istart=0; iend+=3; 
      ny=D->ny+5; nySub+=5; jstart=0; jend+=3;
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
      break;
    }


    //particle resotre
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

    //restore minSub, maxSub
    unitX=D->nx/D->L+1;
    unitY=D->ny/D->M+1;
    unitZ=D->nz/D->N+1;
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

    int cntList[nSpecies];

    for(s=0; s<nSpecies; s++)  {
      sprintf(dataName,"%dtotalCnt",s);
      if(myrank==0)  
        restoreIntMeta(name,dataName,&cntList[s],1);
      else	;
      MPI_Bcast(&cntList[s],1,MPI_INT,0,MPI_COMM_WORLD);
    }

    //set HDF5 parameter
    hid_t file_id,dset_id,plist_id,attr_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
    herr_t ierr;

    for(s=0; s<nSpecies; s++)
    {
      totalCnt=cntList[s];
      if(totalCnt>0)
      {
        //open file
        plist_id=H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
        file_id=H5Fopen(name,H5F_ACC_RDWR,plist_id);
        H5Pclose(plist_id);

        //set dataset
        sprintf(dataName,"%d",s);
        dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
/*
      //attribute
      sprintf(dataName,"%dtotalCnt",s);
      attr_id=H5Aopen(dset_id,dataName,H5P_DEFAULT);
      H5Aread(attr_id,H5T_NATIVE_INT,&totalCnt);
      H5Aclose(attr_id);
      MPI_Barrier(MPI_COMM_WORLD);
*/
        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmp=sub+1;
          else	         tmp=sub;
          if(myrank==rank)
            cntSub=tmp;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvData[i]=0;
          coreCnt[i]=0;
        }

        data = (double *)malloc(cntSub*9*sizeof(double ));

        //file space
        dimsf[0]=totalCnt;
        dimsf[1]=9;
        filespace=H5Screate_simple(2,dimsf,NULL);

        //memory space
        dimsf[0]=cntSub;
        dimsf[1]=9;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cntSub;
        block[1]=9;
        offSet[0]=start;
        offSet[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);

        //hyperslab in memory space
        offSet[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);

        //read data
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Fclose(file_id);

        //for testing the core of each particle
        dataCore = (int *)malloc(cntSub*sizeof(int ));       
        for(i=0; i<cntSub; i++)
          dataCore[i]=-1;

        for(i=0; i<cntSub; i++)  {          
          x=data[i*9+0]-D->minXDomain;         
          y=data[i*9+1];          
          z=data[i*9+2];          
          startIndexX=((int)x)/unitX;
          startIndexY=((int)(y-D->minYDomain))/unitY;
          startIndexZ=((int)(z-D->minZDomain))/unitZ;
          rank=startIndexY+startIndexZ*M+startIndexX*M*N;
          flag=0;

          while(flag==0)  {
            if(minXSub[rank]<=x && x<maxXSub[rank] &&
               minYSub[rank]<=y && y<maxYSub[rank] && 
               minZSub[rank]<=z && z<maxZSub[rank]) { 
              flag=1;           
              if(rank==myrank)  {
                indexI=(int)(x-D->minXSub)+D->istart;
                indexJ=(int)(y-D->minYSub)+D->jstart;
                indexK=(int)(z-D->minZSub)+D->kstart;

                p = (ptclList *)malloc(sizeof(ptclList));
                p->next = particle[indexI][indexJ][indexK].head[s]->pt;
                particle[indexI][indexJ][indexK].head[s]->pt=p;
            
                p->x=(x-D->minXSub)-(int)(x-D->minXSub);
                p->y=(y-D->minYSub)-(int)(y-D->minYSub);
                p->z=(z-D->minZSub)-(int)(z-D->minZSub);
                p->p1=data[i*9+3];
                p->p2=data[i*9+4];
                p->p3=data[i*9+5];
                p->index=data[i*9+6];
                p->core=data[i*9+7];
                p->weight=data[i*9+8];
              }
              dataCore[i]=rank;
              sharePNum[rank]+=1;
            } 
            else 	rank++	;
          }	//End of while(flag=0)
        }

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

        //set memory for send and recving data
        dataCnt=9;
        for(i=0; i<nTasks; i++)   {          
          sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
          recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
        }        

        for(i=0; i<cntSub; i++)  
        {          
          core=dataCore[i];
          if(myrank!=core)  {
            n=coreCnt[core];
            sendData[core][n*dataCnt+0]=data[i*9+0]-D->minXDomain;
            sendData[core][n*dataCnt+1]=data[i*9+1];
            sendData[core][n*dataCnt+2]=data[i*9+2];
            sendData[core][n*dataCnt+3]=data[i*9+3];
            sendData[core][n*dataCnt+4]=data[i*9+4];
            sendData[core][n*dataCnt+5]=data[i*9+5];
            sendData[core][n*dataCnt+6]=data[i*9+6];
            sendData[core][n*dataCnt+7]=data[i*9+7];
            sendData[core][n*dataCnt+8]=data[i*9+8];
            coreCnt[core]+=1;
          }	else	;
        }

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
              z=recvData[i][j*dataCnt+2];
              indexI=(int)(x-D->minXSub)+D->istart;
              indexJ=(int)(y-D->minYSub)+D->jstart;
              indexK=(int)(z-D->minZSub)+D->kstart;
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ][indexK].head[s]->pt;
              particle[indexI][indexJ][indexK].head[s]->pt=p;
            
              p->x=(x-D->minXSub)-(int)(x-D->minXSub);
              p->y=(y-D->minYSub)-(int)(y-D->minYSub);
              p->z=(z-D->minZSub)-(int)(z-D->minZSub);
              p->p1=recvData[i][j*dataCnt+3];
              p->p2=recvData[i][j*dataCnt+4];
              p->p3=recvData[i][j*dataCnt+5];
              p->index=recvData[i][j*dataCnt+6];
              p->core=recvData[i][j*dataCnt+7];
              p->weight=recvData[i][j*dataCnt+8];
            }
          } 
          MPI_Barrier(MPI_COMM_WORLD);
        }


        free(data);
        free(dataCore);
        for(i=0; i<nTasks; i++)   {          
          free(recvData[i]);
          free(sendData[i]);
        }        

      }	else ;	//End of totalCnt>0
    } 		//End of nSpecies 

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
//lala

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

