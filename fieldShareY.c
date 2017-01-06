#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_TransferF_Pukhov_Yminus(Domain *D
         ,double ***f1,double ***f2,double ***f3
         ,int nx,int nz,int share)
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

    rank=(myrank%(D->M*D->N))%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f1[i][j+jstart][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f2[i][j+jstart][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f3[i][j+jstart][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,D->numMinusY,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->minusY,D->numMinusY,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)	
      for(j=0; j<share; j++)
      {
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f1[i][jstart+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f2[i][jstart+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->minusY[start+k]=f3[i][jstart+j][k];
        start+=nz;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(D->minusY,D->numMinusY,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jend+j][k]=D->minusY[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->minusY,D->numMinusY,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}


void MPI_TransferF_Pukhov_Yplus(Domain *D
           ,double ***f1,double ***f2,double ***f3
           ,int nx,int nz,int share)
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

    rank=(myrank%(D->M*D->N))%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f1[i][jend-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f2[i][jend-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f3[i][jend-j][k];
        start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->plusY,D->numPlusY,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(D->plusY,D->numPlusY,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f1[i][jend-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f2[i][jend-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->plusY[start+k]=f3[i][jend-j][k];
        start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(D->plusY,D->numPlusY,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jstart-j][k]=D->plusY[start+k];
          start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->M-1)
      MPI_Send(D->plusY,D->numPlusY,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferJ_Yminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share)
{
    int i,j,k,start;
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

    rank=(myrank%(D->M*D->N))%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f1[i][jstart-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f2[i][jstart-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f3[i][jstart-j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(D->YminusJ,D->numMinusYJ,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(D->YminusJ,D->numMinusYJ,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f1[i][jstart-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f2[i][jstart-j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YminusJ[start+k]=f3[i][jstart-j][k];
        start+=nz;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(D->YminusJ,D->numMinusYJ,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jend-j][k]+=D->YminusJ[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(D->YminusJ,D->numMinusYJ,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

}

void MPI_TransferJ_Yplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share)
{
    int i,j,k,start;
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

    rank=(myrank%(D->M*D->N))%D->M;

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f1[i][jend+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f2[i][jend+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f3[i][jend+j][k];
        start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(D->YplusJ,D->numPlusYJ,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(D->YplusJ,D->numPlusYJ,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f1[i][jend+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f2[i][jend+j][k];
        start+=nz;
        for(k=0; k<nz; k++)
          D->YplusJ[start+k]=f3[i][jend+j][k];
        start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(D->YplusJ,D->numPlusYJ,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f2[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
          for(k=0; k<nz; k++)
            f3[i][jstart+j][k]+=D->YplusJ[start+k];
          start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->M-1)
      MPI_Send(D->YplusJ,D->numPlusYJ,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferRho_Yminus(Domain *D,double ***f1,int nx,int nz,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank;
    double *Yminus; 

    MPI_Status status;         

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))%D->M;
    numShare=D->numMinusYJ/3;
    Yminus = (double *)malloc(numShare*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          Yminus[start+k]=f1[i][jstart-j][k];
        start+=nz;
      }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(Yminus,numShare,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend-j][k]+=Yminus[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(Yminus,numShare,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          Yminus[start+k]=f1[i][jstart-j][k];
        start+=nz;
      }
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(Yminus,numShare,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jend-j][k]+=Yminus[start+k];
          start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(Yminus,numShare,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    free(Yminus);
}

void MPI_TransferRho_Yplus(Domain *D,double ***f1,int nx,int nz,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *Yplus;

    MPI_Status status;         
   
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))%D->M;
    numShare=D->numPlusYJ/3;
    Yplus = (double *)malloc(numShare*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          Yplus[start+k]=f1[i][jend+j][k];
        start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(Yplus,numShare,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart+j][k]+=Yplus[start+k];
          start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(Yplus,numShare,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
      {
        for(k=0; k<nz; k++)
          Yplus[start+k]=f1[i][jend+j][k];
        start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(Yplus,numShare,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
        {
          for(k=0; k<nz; k++)
            f1[i][jstart+j][k]+=Yplus[start+k];
          start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->M-1)
      MPI_Send(Yplus,numShare,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    free(Yplus);
}
