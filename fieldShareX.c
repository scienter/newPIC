#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_Transfer2F_Xminus(Domain *D,double ***f1,double ***f2,int ny,int nz,int share)
{
    int i,j,k,num,start,end;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=2*ny*nz*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	 {
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer3F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share)
{
    int i,j,k,num,start,end;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=3*ny*nz*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	
      {
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[i+istart][j][k]; start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	 {
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[i+istart][j][k]; start+=nz;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int nz,int share)
{
    int i,j,k,num,start,end;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=6*ny*nz*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)   {
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f4[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f5[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f6[i+istart][j][k]; start+=nz;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)        {
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f4[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f5[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f6[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f4[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f5[i+istart][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f6[i+istart][j][k]; start+=nz;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)     {
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f4[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f5[iend+i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f6[iend+i][j][k]=data[start+k]; start+=nz;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer2F_Xplus(Domain *D,double ***f1,double ***f2,int ny,int nz,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    kstart=D->kstart;   kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=2*ny*nz*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)      {
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer3F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    kstart=D->kstart;   kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=3*ny*nz*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[iend-i][j][k]; start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)      {
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[iend-i][j][k]; start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int nz,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    kstart=D->kstart;   kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=6*ny*nz*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	{
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f4[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f5[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f6[iend-i][j][k]; start+=nz;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f4[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f5[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f6[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)      {
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f2[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f3[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f4[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f5[iend-i][j][k]; start+=nz;
        for(k=0; k<nz; k++) data[start+k]=f6[iend-i][j][k]; start+=nz;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	{
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f2[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f3[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f4[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f5[istart-i][j][k]=data[start+k]; start+=nz;
          for(k=0; k<nz; k++) f6[istart-i][j][k]=data[start+k]; start+=nz;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
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
                                                

