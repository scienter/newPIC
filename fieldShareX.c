#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void secondFieldShare(Domain D)
{
  int myrank, nTasks;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  void MPI_TransferF_Pukhov_Xminus();
  void MPI_TransferF_Pukhov_Xplus();
  void MPI_TransferF_Pukhov_Yminus();
  void MPI_TransferF_Pukhov_Yplus();

  switch((D.fieldType-1)*3+D.dimension) {
  case (Pukhov-1)*3+2:
    if(D.L>1)  {
      MPI_TransferF_Pukhov_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
      MPI_TransferF_Pukhov_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
    }
    if(D.M>1)  {
      MPI_TransferF_Pukhov_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
      MPI_TransferF_Pukhov_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    }
    break;
  }
}

void firstFieldShare(Domain D)
{
  int myrank, nTasks;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  void MPI_TransferF_Pukhov_Xminus();
  void MPI_TransferF_Pukhov_Xplus();
  void MPI_TransferF_Pukhov_Yminus();
  void MPI_TransferF_Pukhov_Yplus();

  switch((D.fieldType-1)*3+D.dimension) {
  case (Pukhov-1)*3+2:
    if(D.L>1)  {
      MPI_TransferF_Pukhov_Xminus(&D,D.Bx,D.By,D.Bz,D.nySub+5,1,3);
      MPI_TransferF_Pukhov_Xplus(&D,D.Bx,D.By,D.Bz,D.nySub+5,1,3);
    }
    if(D.M>1)  {
//      MPI_TransferF_Pukhov_Yminus(&D,D.Bx,D.By,D.Bz,D.nxSub+5,1,3);
//      MPI_TransferF_Pukhov_Yplus(&D,D.Bx,D.By,D.Bz,D.nxSub+5,1,3);
    }
    break;
  }
}

void MPI_TransferF_Pukhov_Xminus(Domain *D
         ,double ***f1,double ***f2,double ***f3
         ,int ny,int nz,int share)
{
    int i,j,k,numberData,start,end;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 

    MPI_Status status;         

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);
 
    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<jend+3; j++)	
      {
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f1[i+istart][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f2[i+istart][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f3[i+istart][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(D->minusX,D->numMinusX, MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<jend+3; j++)	
        {
          for(k=0; k<nz; k++)
            f1[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusX,D->numMinusX,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<jend+3; j++)	
      {
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f1[i+istart][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f2[i+istart][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f3[i+istart][j][k];
        start+=nz;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(D->minusX,D->numMinusX, MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<jend+3; j++)	
        {
          for(k=0; k<nz; k++)
            f1[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[iend+i][j][k]=D->minusX[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusX,D->numMinusX,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}


void MPI_TransferF_Pukhov_Xplus(Domain *D
           ,double ***f1,double ***f2,double ***f3
           ,int ny,int nz,int share)
{
    int i,j,k,numberData,end,ibegin,ibottom;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 

    MPI_Status status;         
   
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<jend+3; j++)	
      {
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f1[iend-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f2[iend-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f3[iend-i][j][k];
        start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusX,D->numPlusX,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<jend+3; j++)	
        {
          for(k=0; k<nz; k++)
            f1[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(D->plusX,D->numPlusX,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<jend+3; j++)	
      {
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f1[iend-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f2[iend-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusX[start+k]=f3[iend-i][j][k];
        start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(D->plusX,D->numPlusX,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<jend+3; j++)	
        {
          for(k=0; k<nz; k++)
            f1[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[istart-i][j][k]=D->plusX[start+k];
          start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(D->plusX,D->numPlusX,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferF_Xminus(Domain *D,double ***f1,int istart,int iend,int ny,int nz,int share)
{
    int i,j,k,dataNum,start,end;
    int myrank, nTasks,rank; 
    double *shareData;

    MPI_Status status;         

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);
    dataNum=ny*nz*share;
    shareData=(double *)malloc(dataNum*sizeof(double ));

   //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++)
          D->minusX[start+k]=f1[i+istart][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(shareData,dataNum, MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<nz; k++)
            f1[iend+i][j][k]=shareData[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(shareData,dataNum,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++)
          shareData[start+k]=f1[i+istart][j][k];
        start+=nz;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(shareData,dataNum, MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<nz; k++)
            f1[iend+i][j][k]=shareData[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(shareData,dataNum,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
   
    free(shareData);
}

void MPI_TransferF_Xplus(Domain *D,double ***f1,int istart,int iend,int ny,int nz,int share)
{
    int i,j,k,dataNum,start,end;
    int myrank, nTasks,rank; 
    double *shareData;

    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);
    dataNum=ny*nz*(share-1);
    shareData=(double *)malloc(dataNum*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++)
          shareData[start+k]=f1[iend-i][j][k];
        start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(shareData,dataNum,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<nz; k++)
            f1[istart-i][j][k]=shareData[start+k];
          start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(shareData,dataNum,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++)
          shareData[start+k]=f1[iend-i][j][k];
        start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(shareData,dataNum,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	
        {
          for(k=0; k<nz; k++)
            f1[istart-i][j][k]=shareData[start+k];
          start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(shareData,dataNum,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(shareData);
}

void MPI_TransferJ_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share)
{
    int i,j,k,start;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank/(D->M*D->N);

    //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f1[iend+i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f2[iend+i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f3[iend+i][j][k];
        start+=nz;
      }

    if(rank%2==1)
    {
       MPI_Recv(D->XplusJ,D->numPlusXJ,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f2[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f3[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
         }
    }
    else if(rank%2==0 && rank!=D->L-1)
      MPI_Send(D->XplusJ,D->numPlusXJ,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f1[iend+i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f2[iend+i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XplusJ[start+k]=f3[iend+i][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(D->XplusJ,D->numPlusXJ,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f2[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f3[istart+i][j][k]+=D->XplusJ[k+start];
           start+=nz;
         }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(D->XplusJ,D->numPlusXJ,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}
                                                

void MPI_TransferJ_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share)
{
    int i,j,k,start;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank/(D->M*D->N);

    //Transferring even ~ odd cores 
    start=0;
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f1[istart-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f2[istart-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f3[istart-i][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(D->XminusJ,D->numMinusXJ,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=1; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f2[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f3[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
         }
    }
    else if(rank%2==1)
      MPI_Send(D->XminusJ,D->numMinusXJ,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f1[istart-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f2[istart-i][j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->XminusJ[start+k]=f3[istart-i][j][k];
        start+=nz;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(D->XminusJ,D->numMinusXJ,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=1; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f2[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
           for(k=0; k<nz; k++)
             f3[iend-i][j][k]+=D->XminusJ[k+start];
           start+=nz;
         }
    }
    else if(rank%2==0 && rank!=0)
      MPI_Send(D->XminusJ,D->numMinusXJ,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
}
                                                



void MPI_TransferRho_Xplus(Domain *D,double ***f1,int ny,int nz,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank;
    double *Xplus;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank/(D->M*D->N);
    numShare=D->numPlusXJ/3;
    Xplus = (double *)malloc(numShare*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          Xplus[start+k]=f1[iend+i][j][k];
        start+=nz;
      }

    if(rank%2==1)
    {
       MPI_Recv(Xplus,numShare,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[istart+i][j][k]+=Xplus[k+start];
           start+=nz;
         }
    }
    else if(rank%2==0 && rank!=D->L-1)
      MPI_Send(Xplus,numShare,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          Xplus[start+k]=f1[iend+i][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(Xplus,numShare,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[istart+i][j][k]+=Xplus[k+start];
           start+=nz;
         }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(Xplus,numShare,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(Xplus);
}
                                                

void MPI_TransferRho_Xminus(Domain *D,double ***f1,int ny,int nz,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks, rank;
    double *Xminus;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    rank=myrank/(D->M*D->N);
    numShare=D->numMinusXJ/3;
    Xminus = (double *)malloc(numShare*sizeof(double ));


    //Transferring even ~ odd cores 
    start=0;
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          Xminus[start+k]=f1[istart-i][j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(Xminus,numShare,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=1; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[iend-i][j][k]+=Xminus[k+start];
           start+=nz;
         }
    }
    else if(rank%2==1)
      MPI_Send(Xminus,numShare,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)
      {
        for(k=0; k<nz; k++)
          Xminus[start+k]=f1[istart-i][j][k];
        start+=nz;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(Xminus,numShare,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(i=1; i<share; i++)
         for(j=0; j<ny; j++)
         {
           for(k=0; k<nz; k++)
             f1[iend-i][j][k]+=Xminus[k+start];
           start+=nz;
         }
    }
    else if(rank%2==0 && rank!=0)
      MPI_Send(Xminus,numShare,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(Xminus);
}
                                                

