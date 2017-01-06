#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreIntMeta(char *fileName,char *dataName,int *data);
void saveIntMeta(char *fileName,char *dataName,int *data);
void restoreDoubleData(char *fileName,char *dataName,double *data,int totalCnt,int dataNum,int offSet);
void restoreIntData(char *fileName,char *dataName,int *data,int totalCnt,int dataNum,int offSet);
void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void resolCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int rankX,int edge);
void resolCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int rankY,int edge);
void calBField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX);
void calParameter(int nx,int *istart,int *iend,int *saveNxSub,int rankX,int *biasX,int reIend,int reNxSub,int L);
double ***memoryAsign(int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);
int ***memoryAsignInt(int nx, int ny, int nz);
void deleteFieldInt(int ***field,int nx,int ny,int nz);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
double randomValue(double beta);
double calInter2D(double ***field,int i,int j,double x,double y);



int main(int argc, char *argv[])
{
    int dimension,L,M,N,resolX,resolY,resolZ,step;
    int i,j,k,s,l,m,n,nx,ny,nz,nSpecies,reNx,reNy,reNz,index,number;
    int totalCnt,remain,sub,cntSub,reCntSub,start,rank,cnt,dataCnt,core;
    int startIndexX,startIndexY,startIndexZ,flag;
    int reTotalCnt,reStart,tmpInt,unitX,unitY,unitZ;
    double x,y,z,dx,dy,dz,randTest,tmpDouble,sum;
    double minX,maxX,minY,maxY,minZ,maxZ,**recvData,**sendData;
    int *recv,*dataIndex,*dataCores,*dataCore,*reDataIndex,*reDataCores,*intData,*coreCnt,*sharePNum,*recvDataCnt;
    int ***reCnt;
    double *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz;
    double *reDataX,*reDataY,*reDataZ,*reDataPx,*reDataPy,*reDataPz;
    int minXDomain,minYDomain,minZDomain;
    int minXSub,maxXSub,minYSub,maxYSub,minZSub,maxZSub;
    int *minXSubList,*maxXSubList,*minYSubList,*maxYSubList,*minZSubList,*maxZSubList;
    int istart,iend,jstart,jend,kstart,kend,offset[3],offSet[3];
    int reNxSub,reNySub,reNzSub,reIend,reJend,reKend;
    int saveNxSub,saveNySub,saveNzSub,biasX,biasY,biasZ;
    int saveIstart,saveJstart,saveKstart,saveIend,saveJend,saveKend;
    int remainX,remainY,remainZ,subX,subY,subZ,tmpX,tmpY,tmpZ;
    int rankX,rankY,rankZ,nxSub,nySub,nzSub;
    double ***reField,***reField1,***fieldOld,***fieldNow,***fieldNext;
    double ***density,***reDensity;

    FILE *out;
    char fileName[100],dataName[100],outFile[100],fileName1[100];
    hid_t file_id,group_id;
    herr_t hdfstatus;

    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(argc < 9) 
    {  
      printf("mpirun -np N ./resolParticle dimension step L M N resolX resolY resolZ\n"); 
      exit(0); 
    }

    dimension=atoi(argv[1]);
    step=atoi(argv[2]);
    L=atoi(argv[3]);
    M=atoi(argv[4]);
    N=atoi(argv[5]);
    resolX=atoi(argv[6]);
    resolY=atoi(argv[7]);
    resolZ=atoi(argv[8]);
    if(L*M*N!=nTasks)  {
      printf("check nTasks!. Now L=%d,M=%d,N=%d\n",L,M,N);
      exit(0);
    }


    //save outfile
    sprintf(outFile,"redumpParticle%d.h5",step);
    if(myrank==0)
    {
      file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    //set memory
    recv=(int *)malloc(nTasks*sizeof(int ));
    minXSubList=(int *)malloc(nTasks*sizeof(int ));
    maxXSubList=(int *)malloc(nTasks*sizeof(int ));
    minYSubList=(int *)malloc(nTasks*sizeof(int ));
    maxYSubList=(int *)malloc(nTasks*sizeof(int ));
    minZSubList=(int *)malloc(nTasks*sizeof(int ));
    maxZSubList=(int *)malloc(nTasks*sizeof(int ));
    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
    sharePNum=(int *)malloc(nTasks*sizeof(int ));
    coreCnt=(int *)malloc(nTasks*sizeof(int ));

    //restore meta data
    sprintf(fileName,"dumpParticle%d.h5",step);
    if(myrank==0)  {
      restoreIntMeta(fileName,"/nSpecies",&nSpecies);
      restoreIntMeta(fileName,"/nx",&nx);
      restoreIntMeta(fileName,"/ny",&ny);
      restoreIntMeta(fileName,"/nz",&nz);
      restoreIntMeta(fileName,"/minXDomain",&minXDomain);
      restoreIntMeta(fileName,"/minYDomain",&minYDomain);
      restoreIntMeta(fileName,"/minZDomain",&minZDomain);
    }    else   ;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&minYDomain,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&minZDomain,1,MPI_INT,0,MPI_COMM_WORLD);
    if(minXDomain>0)  minXDomain-=1;  else	;

    nxSub=nx/L;
    subX=nxSub;
    remainX=nx%L;
    minX=maxX=0;
    nySub=ny/M;
    subY=nySub;
    remainY=ny%M;
    minY=maxY=0;
    nzSub=nz/N;
    subZ=nzSub;
    remainZ=nz%N;
    minZ=maxZ=0;

    minX=maxX=0;
    for(rankX=0; rankX<L; rankX++)
    {
      if(rankX<remainX)   tmpX=subX+1;
      else                tmpX=subX;
      minX=maxX;
      maxX=minX+tmpX;

      minZ=maxZ=minZDomain;
      for(rankZ=0; rankZ<N; rankZ++)
      {
        if(rankZ<remainZ)   tmpZ=subZ+1;
        else                tmpZ=subZ;
        minZ=maxZ;
        maxZ=minZ+tmpZ;

        minY=maxY=minYDomain;
        for(rankY=0; rankY<M; rankY++)
        {
          if(rankY<remainY)   tmpY=subY+1;
          else                tmpY=subY;
          minY=maxY;
          maxY=minY+tmpY;

          rank=rankY+rankZ*M+rankX*(M*N);
          if(myrank==rank)
          {
             nxSub=tmpX;
             nySub=tmpY;
             nzSub=tmpZ;
             minXSub=minX;
             maxXSub=maxX;
             minYSub=minY;
             maxYSub=maxY;
             minZSub=minZ;
             maxZSub=maxZ;
          }
        }
      }
    }
    istart=2;
    iend=nxSub+2;
    jstart=0;
    jend=1;
    kstart=0;
    kend=1;
    if(dimension>1)  {
      jstart=2;
      jend=nySub+2;
    }
    if(dimension>2)  {
      kstart=2;
      kend=nzSub+2;
    }

    rankX=myrank/(M*N);
    rankZ=(myrank%(M*N))/M;
    rankY=(myrank%(M*N))%M;

    //setting resolution change
    reNx=nx*resolX;
    reNy=ny*resolY;
    reNz=nz*resolZ;

    reNxSub=nxSub*resolX;
    reNySub=nySub*resolY;
    reNzSub=nzSub*resolZ;
    reIend=reNxSub+2;
    reJend=reKend=1;
    if(dimension>1)
      reJend=reNySub+2;
    if(dimension>2)
      reKend=reNzSub+2;

    //set domain range
    minX=minXDomain*resolX;
    maxX=(minXDomain+nx)*resolX;
    minY=minYDomain*resolY;
    maxY=(minYDomain+ny)*resolY;
    minZ=minZDomain*resolZ;
    maxZ=(minZDomain+ny)*resolZ;

    calParameter(reNx+5,&saveIstart,&saveIend,&saveNxSub,rankX,&biasX,reIend,reNxSub,L);
    saveJstart=saveKstart=0;
    saveKend=saveKend=saveNySub=saveNzSub=1;
    if(dimension>1)
      calParameter(reNy+5,&saveJstart,&saveJend,&saveNySub,rankY,&biasY,reJend,reNySub,M);
    if(dimension>2)
      calParameter(reNz+5,&saveKstart,&saveKend,&saveNzSub,rankZ,&biasZ,reKend,reNzSub,N);

    switch (dimension)  {
    case 2 :
      fieldOld=memoryAsign(nxSub+5,nySub+5,1);
      fieldNow=memoryAsign(nxSub+5,nySub+5,1);
      fieldNext=memoryAsign(nxSub+5,nySub+5,1);

      offset[0]=minXSub;
      offset[1]=minYSub-minYDomain;
      offset[2]=0;

      offSet[0]=(minXSub)*resolX+biasX;
      offSet[1]=(minYSub-minYDomain)*resolY+biasY;
      offSet[2]=0;

      //restore minSub, maxSub
      MPI_Gather(&minXSub,1,MPI_INT,minXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&maxXSub,1,MPI_INT,maxXSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxXSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&minYSub,1,MPI_INT,minYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&maxYSub,1,MPI_INT,maxYSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxYSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&minZSub,1,MPI_INT,minZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(minZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Gather(&maxZSub,1,MPI_INT,maxZSubList,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(maxZSubList,nTasks,MPI_INT,0,MPI_COMM_WORLD);

      for(s=0; s<nSpecies; s++)
      {

        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvData[i]=0;
          coreCnt[i]=0;
          recvDataCnt[i]=0;
        }

        sprintf(dataName,"%dParticle/totalCnt",s);
        if(myrank==0) 
          restoreIntMeta(fileName,dataName,&totalCnt);
        else    ;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmpInt=sub+1;
          else             tmpInt=sub;
          if(myrank==rank)
            cntSub=tmpInt;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        dataX = (double *)malloc(cntSub*sizeof(double ));
        dataY = (double *)malloc(cntSub*sizeof(double ));
        dataPx = (double *)malloc(cntSub*sizeof(double ));
        dataPy = (double *)malloc(cntSub*sizeof(double ));
        dataPz = (double *)malloc(cntSub*sizeof(double ));
        dataIndex = (int *)malloc(cntSub*sizeof(int ));
        dataCore = (int *)malloc(cntSub*sizeof(int ));
        dataCores = (int *)malloc(cntSub*sizeof(int ));

        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        if(totalCnt>0)
        {
          sprintf(dataName,"%dParticle/x",s);
          restoreDoubleData(fileName,dataName,dataX,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/y",s);
          restoreDoubleData(fileName,dataName,dataY,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/px",s);
          restoreDoubleData(fileName,dataName,dataPx,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/py",s);
          restoreDoubleData(fileName,dataName,dataPy,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/pz",s);
          restoreDoubleData(fileName,dataName,dataPz,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/index",s);
          restoreIntData(fileName,dataName,dataIndex,totalCnt,cntSub,start);
          sprintf(dataName,"%dParticle/core",s);
          restoreIntData(fileName,dataName,dataCores,totalCnt,cntSub,start);

          for(i=0; i<cntSub; i++) {
            dataCore[i]=-1;
          }

          cnt=0;
          unitX=nx/L+1;
          unitY=ny/M+1;
          for(i=0; i<cntSub; i++)  {
            x=dataX[i]-minXDomain;
            y=dataY[i]-minYDomain;
            startIndexX=((int)x)/unitX;
            startIndexY=((int)y)/unitY;
            startIndexZ=0;
            rank=startIndexY+startIndexZ*M+startIndexX*M*N;
            flag=0;
            x=dataX[i];
            y=dataY[i];
            while(flag==0 && rank<nTasks)  {
              if(minXSubList[rank]<=x && x<maxXSubList[rank] &&
                 minYSubList[rank]<=y && y<maxYSubList[rank]) {
                flag=1;
                if(rank==myrank)  cnt++;     else ;
                dataCore[i]=rank;
                sharePNum[rank]+=1;
              }
              else   {
                rank++  ;
              }
            }
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
              MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
            }  else     ;
          }
          MPI_Barrier(MPI_COMM_WORLD);

          reCntSub=cnt;
          for(i=0; i<nTasks; i++)
            reCntSub+=recvDataCnt[i];

          reDataX = (double *)malloc(reCntSub*sizeof(double ));
          reDataY = (double *)malloc(reCntSub*sizeof(double ));
          reDataPx = (double *)malloc(reCntSub*sizeof(double ));
          reDataPy = (double *)malloc(reCntSub*sizeof(double ));
          reDataPz = (double *)malloc(reCntSub*sizeof(double ));
          reDataIndex = (int *)malloc(reCntSub*sizeof(int ));
          reDataCores = (int *)malloc(reCntSub*sizeof(int ));

          cnt=0;
          for(i=0; i<cntSub; i++)  {
            x=dataX[i]-minXDomain;
            y=dataY[i]-minYDomain;
            startIndexX=((int)x)/unitX;
            startIndexY=((int)y)/unitY;
            startIndexZ=0;
            rank=startIndexY+startIndexZ*M+startIndexX*M*N;
            flag=0;
            x=dataX[i];
            y=dataY[i];
            while(flag==0 && rank<nTasks)  {
              if(minXSubList[rank]<=x && x<maxXSubList[rank] &&
                 minYSubList[rank]<=y && y<maxYSubList[rank]) {
                flag=1;
                if(rank==myrank)  {
                  reDataX[cnt]=dataX[i];
                  reDataY[cnt]=dataY[i];
                  reDataPx[cnt]=dataPx[i];
                  reDataPy[cnt]=dataPy[i];
                  reDataPz[cnt]=dataPz[i];
                  reDataIndex[cnt]=dataIndex[i];
                  reDataCores[cnt]=dataCores[i];
                  cnt++;                    
                }    else ;
              }
              else   {
                rank++  ;
              }
            }
          }	//End of for(i ~ cntSub)

          //memory for send and recving data
          dataCnt=7;
          for(i=0; i<nTasks; i++)   {          
            sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
            recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
          }        
          for(i=0; i<cntSub; i++)  
          {          
            core=dataCore[i];
            x=dataX[i];          
            y=dataY[i];
            if(myrank!=core)  {
              n=coreCnt[core];
              sendData[core][n*dataCnt+0]=dataX[i];
              sendData[core][n*dataCnt+1]=dataY[i];
              sendData[core][n*dataCnt+2]=dataPx[i];
              sendData[core][n*dataCnt+3]=dataPy[i];
              sendData[core][n*dataCnt+4]=dataPz[i];
              sendData[core][n*dataCnt+5]=dataIndex[i];
              sendData[core][n*dataCnt+6]=dataCores[i];
              coreCnt[core]+=1;
            }   else    ;
          }

          for(i=0; i<nTasks; i++)
          {
            if(myrank==i)  {
              for(j=0; j<nTasks; j++)
                if(i!=j)
                  MPI_Send(sendData[j],sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
            } 
            else  {
              MPI_Recv(recvData[i],recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
              for(j=0; j<recvDataCnt[i]; j++)  {
                reDataX[cnt]=recvData[i][j*dataCnt+0];
                reDataY[cnt]=recvData[i][j*dataCnt+1];
                reDataPx[cnt]=recvData[i][j*dataCnt+1];
                reDataPy[cnt]=recvData[i][j*dataCnt+1];
                reDataPz[cnt]=recvData[i][j*dataCnt+1];
                reDataIndex[cnt]=recvData[i][j*dataCnt+1];
                reDataCores[cnt]=recvData[i][j*dataCnt+1];
                cnt++;                    
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);
          }

sprintf(fileName,"particle%d",myrank);
out=fopen(fileName,"w");
for(i=0; i<cnt; i++)
   fprintf(out,"%g %g\n",reDataX[i],reDataY[i]);
fclose(out);

          for(i=0; i<nTasks; i++)   {          
            free(recvData[i]);
            free(sendData[i]);
          } 
          free(reDataX);
          free(reDataY);
          free(reDataPx);
          free(reDataPy);
          free(reDataPz);
          free(reDataIndex);
          free(reDataCores);
        }	//End of totalCnt>0

        free(dataX);
        free(dataY);
        free(dataPx);
        free(dataPy);
        free(dataPz);
        free(dataIndex);
        free(dataCore);
        free(dataCores);

      }	//End of nSpecies
//        deleteField(reField,reNxSub+5,nySub+5,1);
//        deleteField(reField1,reNxSub+5,reNySub+5,1);

   
      deleteField(fieldOld,nxSub+5,nySub+5,1);
      deleteField(fieldNow,nxSub+5,nySub+5,1);
      deleteField(fieldNext,nxSub+5,nySub+5,1);

      break;
    }

/*
    minXDomain*=resolX;
    minYDomain*=resolY;
    minZDomain*=resolZ;
    if(myrank==0)  {
      saveIntMeta(outFile,"/nSpecies",&nSpecies);
      saveIntMeta(outFile,"/nx",&reNx);
      saveIntMeta(outFile,"/ny",&reNy);
      saveIntMeta(outFile,"/nz",&reNz);
      saveIntMeta(outFile,"/minXDomain",&minXDomain);
      saveIntMeta(outFile,"/minYDomain",&minYDomain);
      saveIntMeta(outFile,"/minZDomain",&minZDomain);
    }    else   ;
*/ 
    free(recv);
    free(recvData);
    free(sendData);
    free(recvDataCnt);
    free(coreCnt);
    free(sharePNum);
    free(minXSubList);
    free(maxXSubList);
    free(minYSubList);
    free(maxYSubList);
    free(minZSubList);
    free(maxZSubList);
  
    MPI_Finalize();

    return 0;
}

void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
  memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
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

void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;
  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
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

void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
  int i,j,k,start;
  double *field;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
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

void calBField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX)
{
  int i,j,k;
  double y1,y2,y3,a,b,c,dx,x;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        y1=fieldOld[i][j][k];
        y2=fieldNow[i][j][k];
        y3=fieldNext[i][j][k];
        a=0.5*(y1+y3)-y2;
        b=2.0*y2-1.5*y1-0.5*y3;
        c=y1;
        dx=1.0/((double)resolX);
        x=0.5+dx*0.5;
        fieldNow[i][j][k]=(a*x*x+b*x+c);
      }
}

void resolCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int rankY,int edge)
{
  int jj,j1,i,j,k,n;
  double y1,y2,y3,a,b,c,dy,yy,y;

  dy=1.0/((double)resolY);
  if(edge==1)
  {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(n=0; n<resolY; n++)  {
            yy=dy*n;
            jj=((int)(0.5+yy))+j;
            y1=field[i][jj-1][k];
            y2=field[i][jj][k];
            y3=field[i][jj+1][k];
            a=0.5*(y1+y3)-y2;
            b=2.0*y2-1.5*y1-0.5*y3;
            c=y1;
            y=1.0+yy-(int)(0.5+yy);
            j1=(j-jstart)*resolY+jstart;
            reField[i][j1+n][k]=a*y*y+b*y+c;
          }
        }
  }
  else          //edge==0
  {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(n=0; n<resolY; n++)     {
            y=0.5+0.5*dy+dy*n;
            y1=field[i][j-1][k];
            y2=field[i][j][k];
            y3=field[i][j+1][k];
            a=0.5*(y1+y3)-y2;
            b=2.0*y2-1.5*y1-0.5*y3;
            c=y1;
            j1=(j-jstart)*resolY+jstart;
            reField[i][j1+n][k]=a*y*y+b*y+c;
          }
        }
  }             //End of edge==0

}

void resolCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int rankX,int edge)
{
                             int ii,i1,i,j,k,n;
  double y1,y2,y3,a,b,c,dx,xx,x;

  dx=1.0/((double)resolX);
  if(edge==1)
  {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(n=0; n<resolX; n++)  {
            xx=dx*n;
            ii=((int)(0.5+xx-dx*dx))+i;
            y1=field[ii-1][j][k];
            y2=field[ii][j][k];
            y3=field[ii+1][j][k];
            a=0.5*(y1+y3)-y2;
            b=2.0*y2-1.5*y1-0.5*y3;
            c=y1;
            x=1.0+xx-(int)(0.5+xx);
            i1=(i-istart)*resolX+istart;
            reField[i1+n][j][k]=a*x*x+b*x+c;
          }
        }
  }
  else          //edge==0
  {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(n=0; n<resolX; n++)     {
            x=0.5+0.5*dx+dx*n;
            y1=field[i-1][j][k];
            y2=field[i][j][k];
            y3=field[i+1][j][k];
            a=0.5*(y1+y3)-y2;
            b=2.0*y2-1.5*y1-0.5*y3;
            c=y1;
            i1=(i-istart)*resolX+istart;
            reField[i1+n][j][k]=a*x*x+b*x+c;
          }
        }
  }             //End of edge==0

}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=100, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}

double calInter2D(double ***field,int i,int j,double x,double y)
{
   double result;

   result=(1.0-x)*(1.0-y)*field[i][j][0]
         +(1.0-x)*     y *field[i][j+1][0]
         +     x *(1.0-y)*field[i+1][j][0]
         +     x *     y *field[i+1][j+1][0];
   return result;
}
