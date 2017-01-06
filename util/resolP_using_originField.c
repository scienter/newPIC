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
double calMomentum(double old2,double old1,double p,int resolX);

int main(int argc, char *argv[])
{
    int dimension,L,M,N,resolX,resolY,resolZ,step;
    int i,j,k,s,n,nx,ny,nz,nSpecies,reNx,reNy,reNz,index;
    int totalCnt,remain,sub,cntSub,start,rank,outCnt;
    int reTotalCnt,reCntSub,reStart,tmpInt,minXDomain,minYDomain,minZDomain;
    double x,y,z,dx,dy,dz,xx,yy,zz,px,py,pz,old1,old2,tmpDouble;
    double minX,maxX,minY,maxY,minZ,maxZ;
    int *recv,*dataIndex,*dataCores,*reDataIndex,*reDataCores,*intData;
    int *testList;
    double *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz,*doubleData;
    double *dataOld1Px,*dataOld1Py,*dataOld1Pz;
    double *dataOld2Px,*dataOld2Py,*dataOld2Pz;

    FILE *out;
    char fileName[100],dataName[100],outFile[100],fileName1[100];
    hid_t file_id,group_id;
    herr_t hdfstatus;

    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(argc < 8) 
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

    //set memory
    recv=(int *)malloc(nTasks*sizeof(int ));


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

    //save outfile
    sprintf(outFile,"redumpParticle%d.h5",step);
    if(myrank==0)
    {
      file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;
    reNx=nx;
    reNy=ny;
    reNz=nz;
    minX=minXDomain;
    maxX=(minXDomain+nx);
    minY=minYDomain;
    maxY=(minYDomain+ny);
    minZ=minZDomain;
    maxZ=(minZDomain+ny);

    switch (dimension)  {
    case 2 :

      if(myrank==0)
        saveIntMeta(outFile,"/nSpecies",&nSpecies);
      else	;

      for(s=0; s<nSpecies; s++)
      {
        if(myrank==0)    {
          file_id=H5Fopen(outFile,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(dataName,"%dParticle",s);
          group_id=H5Gcreate2(file_id,dataName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }  else         ;
        MPI_Barrier(MPI_COMM_WORLD);

        sprintf(dataName,"%dParticle/totalCnt",s);
        if(myrank==0)
          restoreIntMeta(fileName,dataName,&totalCnt);
        else    ;
        MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

        //calculating cntSub
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
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        dataX = (double *)malloc(cntSub*sizeof(double ));
        dataY = (double *)malloc(cntSub*sizeof(double ));
        dataPx = (double *)malloc(cntSub*sizeof(double ));
        dataPy = (double *)malloc(cntSub*sizeof(double ));
        dataPz = (double *)malloc(cntSub*sizeof(double ));
        dataIndex = (int *)malloc(cntSub*sizeof(int ));
        dataCores = (int *)malloc(cntSub*sizeof(int ));
        dataOld1Px = (double *)malloc(cntSub*sizeof(double ));
        dataOld1Py = (double *)malloc(cntSub*sizeof(double ));
        dataOld1Pz = (double *)malloc(cntSub*sizeof(double ));
        dataOld2Px = (double *)malloc(cntSub*sizeof(double ));
        dataOld2Py = (double *)malloc(cntSub*sizeof(double ));
        dataOld2Pz = (double *)malloc(cntSub*sizeof(double ));

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
        sprintf(dataName,"%dParticle/old1Px",s);
        restoreDoubleData(fileName,dataName,dataOld1Px,totalCnt,cntSub,start);
        sprintf(dataName,"%dParticle/old1Py",s);
        restoreDoubleData(fileName,dataName,dataOld1Py,totalCnt,cntSub,start);
        sprintf(dataName,"%dParticle/old1Pz",s);
        restoreDoubleData(fileName,dataName,dataOld1Pz,totalCnt,cntSub,start);
        sprintf(dataName,"%dParticle/old2Px",s);
        restoreDoubleData(fileName,dataName,dataOld2Px,totalCnt,cntSub,start);
        sprintf(dataName,"%dParticle/old2Py",s);
        restoreDoubleData(fileName,dataName,dataOld2Py,totalCnt,cntSub,start);
        sprintf(dataName,"%dParticle/old2Pz",s);
        restoreDoubleData(fileName,dataName,dataOld2Pz,totalCnt,cntSub,start);

        dx=1.0/(double)(resolX);
        dy=1.0/(double)(resolY);
        n=resolX*resolY;
        reCntSub=n*cntSub;       
        reTotalCnt=n*totalCnt;       
        reStart=n*start;      
        testList = (int *)malloc(reCntSub*sizeof(int ));
        outCnt=0;
        for(i=0; i<cntSub; i++)   {
          x=dataX[i];
          y=dataY[i];
          for(j=0; j<resolY; j++)
            for(k=0; k<resolX; k++) {
              xx=(x-0.5+dx*(0.5+k));
              yy=(y-0.5+dy*(0.5+j));
              if(xx<minX || xx>=maxX || yy<minY || yy>=maxY)  {
                testList[i*n+j*resolX+k]=1;
                outCnt++;
              }
              else   testList[i*n+j*resolX+k]=0;
            }
        }

        //recalculation of saving memory
        reCntSub=cntSub*n-outCnt;
        MPI_Gather(&reCntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        reTotalCnt=0;
        for(i=0; i<nTasks; i++)
          reTotalCnt+=recv[i];
        reStart=0;
        for(i=0; i<myrank; i++)
          reStart+=recv[i];
        MPI_Barrier(MPI_COMM_WORLD);
 
        doubleData = (double *)malloc(reCntSub*sizeof(double ));
//        doubleData2 = (double *)malloc(reCntSub*sizeof(double ));
        intData = (int *)malloc(reCntSub*sizeof(int ));

        index=0;
        for(i=0; i<cntSub; i++)   {
          x=dataX[i];
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++) 
              if(testList[i*n+j*resolX+k]==0)  {
                doubleData[index]=(x-0.5+dx*(0.5+k));
                index++;
              }  else	;
        }
        sprintf(dataName,"%dParticle/x",s);
        saveParticleComp_Double(doubleData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   {
          y=dataY[i];
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  
              if(testList[i*n+j*resolX+k]==0)  {
                doubleData[index]=(y-0.5+dy*(0.5+j));
                index++;
              }  else	;
        } 
        sprintf(dataName,"%dParticle/y",s);
        saveParticleComp_Double(doubleData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  {
              if(testList[i*n+j*resolX+k]==0)  {
//                old2=dataOld2Px[i];
//                old1=dataOld1Px[i];
//                px=dataPx[i];
//                doubleData[index]=calMomentum(old2,old1,px,resolX);
                doubleData[index]=dataPx[i];
                index++;
              }  else	;
            }
        sprintf(dataName,"%dParticle/px",s);
        saveParticleComp_Double(doubleData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  {
              if(testList[i*n+j*resolX+k]==0)  {
//                old2=dataOld2Py[i];
//                old1=dataOld1Py[i];
//                py=dataPy[i];
                doubleData[index]=dataPy[i];
//                doubleData[index]=calMomentum(old2,old1,py,resolX);
                index++;
              }  else	;
            }
        sprintf(dataName,"%dParticle/py",s);
        saveParticleComp_Double(doubleData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  {
              if(testList[i*n+j*resolX+k]==0)  {
//                old2=dataOld2Pz[i];
//                old1=dataOld1Pz[i];
//                py=dataPz[i];
//                doubleData[index]=calMomentum(old2,old1,pz,resolX);
                doubleData[index]=dataPz[i];
                index++;
              }  else	;
            }
        sprintf(dataName,"%dParticle/pz",s);
        saveParticleComp_Double(doubleData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  {
              if(testList[i*n+j*resolX+k]==0)  {
                intData[index]=dataIndex[i];
                index++;
              }  else	;
            }
        sprintf(dataName,"%dParticle/index",s);
        saveParticleComp_Int(intData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        index=0;
        for(i=0; i<cntSub; i++)   
          for(j=0; j<resolY; j++)  
            for(k=0; k<resolX; k++)  {
              if(testList[i*n+j*resolX+k]==0)  {
                intData[index]=dataCores[i];
                index++;
              }  else	;
            }
        sprintf(dataName,"%dParticle/core",s);
        saveParticleComp_Int(intData,outFile,dataName,reTotalCnt,reCntSub,reStart);

        sprintf(dataName,"%dParticle/totalCnt",s);
        if(myrank==0)
          saveIntMeta(outFile,dataName,&reTotalCnt);
        else	;

        //free memory
        free(dataX);
        free(dataY);
        free(dataPx);
        free(dataPy);
        free(dataPz);
        free(dataIndex);
        free(dataCores);
        free(dataOld1Px);
        free(dataOld1Py);
        free(dataOld1Pz);
        free(dataOld2Px);
        free(dataOld2Py);
        free(dataOld2Pz);
        free(testList);
        free(doubleData);
        free(intData);

      }		//End of nSpecies

      if(myrank==0)
        printf("%s is made.\n",outFile);
      break;
    }
    free(recv);

    if(myrank==0) {
      saveIntMeta(outFile,"/minXDomain",&minXDomain);
      saveIntMeta(outFile,"/minYDomain",&minYDomain);
      saveIntMeta(outFile,"/minZDomain",&minZDomain);
    }    else	;
    
    MPI_Finalize();

    return 0;
}

double calMomentum(double old2,double old1,double p,int resolX)
{
  double y1,y2,y3,a,b,c,dx,x,result;

  a=0.5*(old2+p)-old1;
  b=2.0*old1-1.5*old2-0.5*p;
  c=old2;
  dx=1.0/((double)resolX);
  x=2.0-dx*0.5;
  
  result=a*x*x+b*x+c;
  return result;
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

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

