#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "mesh.h"


int main(int argc, char *argv[])
{
   int s,i,j,origin_Offset[3],reOffset[3];
   double ***den1,***den2,***den3,***den4;
   char fileName[100],fileName1[100],outFile[100],dataName[100];
   char testName[100];
   FILE *out;   
   Domain D;
   int myrank, nTasks;
   MPI_Status status; 

   ptclList *p;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   hid_t file_id,group_id;

   if(argc < 9) 
   {  
     printf("mpirun -np N ./resolParticle dimension mode step L M N resolX resolY resolZ\n"); 
     printf("mode 1 -> high resolution change for particle\n"); 
     exit(0); 
   }

   D.dimension=atoi(argv[1]);
   D.mode=atoi(argv[2]);
   D.step=atoi(argv[3]);
   D.L=atoi(argv[4]);
   D.M=atoi(argv[5]);
   D.N=atoi(argv[6]);
   D.resolX=atoi(argv[7]);
   D.resolY=atoi(argv[8]);
   D.resolZ=atoi(argv[9]);
   if(D.L*D.M*D.N!=nTasks)  {
     printf("check nTasks!. Now L=%d,M=%d,N=%d\n",D.L,D.M,D.N);
     exit(0);
   }

   boundary(&D);
   MPI_Barrier(MPI_COMM_WORLD); 

   //save outfile
   sprintf(outFile,"redumpParticle%d.h5",D.step);
   if(myrank==0)     {
       file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
       H5Fclose(file_id);
   }   else ;

   switch (D.mode*3+D.dimension)  {
   
   case (1*3+2) :	//high resolution particle data
     resolHighBoundary(&D);

     den1=memoryAsign(D.nxSub+5,D.nySub+5,1);
     den2=memoryAsign(D.nxSub+5,D.nySub+5,1);
     den3=memoryAsign(D.nxSub+5,D.nySub+5,1);
     den4=memoryAsign(D.nxSub+5,D.nySub+5,1);

     origin_Offset[0]=D.minXSub-D.minXDomain;
     origin_Offset[1]=D.minYSub-D.minYDomain;
     origin_Offset[2]=0;
     reOffset[0]=(D.minXSub-D.minXDomain)*D.resolX+D.biasX;
     reOffset[1]=(D.minYSub-D.minYDomain)*D.resolY+D.biasY;
     reOffset[2]=0;

     sprintf(fileName,"dumpParticle%d.h5",D.step);
     for(s=0; s<D.nSpecies; s++)
     {
       highParticleSave(&D,fileName,s);
       restoreDen(&D,den1,den2,den3,den4,s,"px",D.step,origin_Offset);
       findAveragePx(&D,den1,den2,den3,den4,s);
//       restoreDen(&D,den1,den2,den3,den4,s,"py",D.step,origin_Offset);
//       findAveragePy(&D,den1,den2,den3,den4,s);
//       restoreDen(&D,den1,den2,den3,den4,s,"pz",D.step,origin_Offset);
//       findAveragePz(&D,den1,den2,den3,den4,s);


       sprintf(testName,"denP%d",myrank);
       out=fopen(testName,"w");
       for(i=D.istart; i<D.iend; i++)
       {
         for(j=D.jstart; j<D.jend; j++) 
           fprintf(out,"%d %d %g\n",i,j,den1[i][j][0]);
         fprintf(out,"\n");
       }
       fclose(out);

       sprintf(testName,"particle%d",myrank);
       out=fopen(testName,"w");
       for(i=D.istart; i<D.iend; i++)
         for(j=D.jstart; j<D.jend; j++)  { 
           p=D.particle[i][j][0].head[s]->pt;
           while(p)  {
             p->x=p->x-D.minXDomain;   p->y=p->y-D.minYDomain;
             fprintf(out,"%g %g %g\n",p->x,p->y,p->p1);
             p=p->next;
           }
         }
       fclose(out);

       deleteParticleMemory(&D,s);
     }

     deleteField(den1,D.nxSub+5,D.nySub+5,1);
     deleteField(den2,D.nxSub+5,D.nySub+5,1);
     deleteField(den3,D.nxSub+5,D.nySub+5,1);
     deleteField(den4,D.nxSub+5,D.nySub+5,1);
     cleanSubBoundary(&D);
     break;
   }

/*

    switch (dimension)  {
    case 2 :

      for(s=0; s<nSpecies; s++)
      {

        if(totalCnt>0)
        {
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
*/ 
    cleanBoundary(&D);
 
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
