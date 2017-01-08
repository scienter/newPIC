#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void boundary(Domain *D,External *Ext)
{
     FILE *out;
     int i,j,k,s,rank,rankX,rankY,rankZ,trackStart;
     int remainX,remainY,remainZ,subX,subY,subZ,tmpX,tmpY,tmpZ;
     int nxSub,nySub,numdataUp,numdataBt,numberData;
     int startj,startk,nxSub1D,nySub2D,nzSub3D;
     int minX,maxX,minY,maxY,minZ,maxZ;
     int myrank, nTasks,a;
     double ***memoryAsign(int nx, int ny, int nz);
     double ***memoryAsignJ(int nx, int ny, int nz);
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     D->nxSub=D->nx/D->L;
     subX=D->nxSub;
     remainX=D->nx%D->L;
     minX=maxX=0;

     D->nySub=D->ny/D->M;
     subY=D->nySub;
     remainY=D->ny%D->M;
     minY=maxY=0;

     D->nzSub=D->nz/D->N;
     subZ=D->nzSub;
     remainZ=D->nz%D->N;
     minZ=maxZ=0;

     minX=maxX=D->minXDomain;
     for(rankX=0; rankX<D->L; rankX++)
     {
       if(rankX<remainX)   tmpX=subX+1;
       else                tmpX=subX;
       minX=maxX;
       maxX=minX+tmpX;

       minZ=maxZ=D->minZDomain;
       for(rankZ=0; rankZ<D->N; rankZ++)
       {
         if(rankZ<remainZ)   tmpZ=subZ+1;
         else                tmpZ=subZ;
         minZ=maxZ;
         maxZ=minZ+tmpZ;

         minY=maxY=D->minYDomain;
         for(rankY=0; rankY<D->M; rankY++)
         {
           if(rankY<remainY)   tmpY=subY+1;
           else                tmpY=subY;
           minY=maxY;
           maxY=minY+tmpY;

           rank=rankY+rankZ*D->M+rankX*(D->M*D->N);
           if(myrank==rank)
           {
              D->minXSub=minX;
              D->maxXSub=maxX;
              D->nxSub=tmpX;
              D->minYSub=minY;
              D->maxYSub=maxY;
              D->nySub=tmpY;
              D->minZSub=minZ;
              D->maxZSub=maxZ;
              D->nzSub=tmpZ;
           }
         }
       }
     }

     for(rank=0; rank<nTasks; rank++)
     {
       if(myrank==rank)
       {
         if(D->dimension==1)
         {
           printf("rank=%d, minXSub=%d,maxSub=%d\n",myrank,D->minXSub,D->maxXSub);
         }
         else if(D->dimension==2)
         {
           printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub);
         }
         else if(D->dimension==3)
         {
           printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d,minZSub=%d,maxZSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub,D->minZSub,D->maxZSub);
         }
       }
     }

     //defining next and prev domain
     rankX=myrank/(D->M*D->N);
     rankZ=(myrank%(D->M*D->N))/D->M;
     rankY=(myrank%(D->M*D->N))%D->M;

     D->nextXrank=rankY+rankZ*D->M+(rankX+1)*(D->M*D->N);
     D->prevXrank=rankY+rankZ*D->M+(rankX-1)*(D->M*D->N);
     if(rankX==D->L-1) 
       D->nextXrank=rankY+rankZ*D->M;
     else if(rankX==0) 
       D->prevXrank=rankY+rankZ*D->M+(D->L-1)*(D->M*D->N);

     D->nextYrank=(rankY+1)+rankZ*D->M+rankX*(D->M*D->N);
     D->prevYrank=(rankY-1)+rankZ*D->M+rankX*(D->M*D->N);
     if(rankY==D->M-1) 
       D->nextYrank=rankZ*D->M+rankX*(D->M*D->N);
     else if(rankY==0) 
       D->prevYrank=(D->M-1)+rankZ*D->M+rankX*(D->M*D->N);


     D->nextZrank=rankY+(rankZ+1)*D->M+rankX*(D->M*D->N);
     D->prevZrank=rankY+(rankZ-1)*D->M+rankX*(D->M*D->N);
     if(rankZ==D->N-1) 
       D->nextZrank=rankY+rankX*(D->M*D->N);
     else if(rankZ==0) 
       D->prevZrank=rankY+(D->N-1)*D->M+rankX*(D->M*D->N);


     D->istart=2;
     D->iend=D->nxSub+2;
     D->jstart=0;
     D->jend=1;
     D->kstart=0;
     D->kend=1;
     if(D->dimension>1)  {
        D->jstart=2;
        D->jend=D->nySub+2;
     }
     if(D->dimension>2)  {
        D->kstart=2;
        D->kend=D->nzSub+2;
     }

     // Field setting
     nxSub1D=D->nxSub+5;
     nySub2D=1;
     nzSub3D=1;
     if(D->dimension>1)  
       nySub2D=D->nySub+5;
     if(D->dimension>2)  
       nzSub3D=D->nzSub+5;

     MPI_Barrier(MPI_COMM_WORLD);

    //Share Field
    switch((D->fieldType-1)*3+D->dimension) {
    case (Split-1)*3+1 :
      //first 3 is 3 field variables, 2nd 1 is y, 3rd 1 is z, 4rd 2 or 3 is layer.
      D->numPlusX=3*nySub2D*1*2;
      D->plusX = (double *)malloc(D->numPlusX*sizeof(double ));
      D->numMinusX=3*nySub2D*1*3;
      D->minusX = (double *)malloc(D->numMinusX*sizeof(double ));
      break;
    case (Pukhov-1)*3+2:
      //first 3 is 3 field variables, 2nd 1 is k, 3rd 2 or 3 is layer.
      D->numPlusX=3*nySub2D*1*2;
      D->plusX = (double *)malloc(D->numPlusX*sizeof(double ));
      D->numMinusX=3*nySub2D*1*3;
      D->minusX = (double *)malloc(D->numMinusX*sizeof(double ));
      D->numPlusY=3*nxSub1D*1*2;
      D->plusY = (double *)malloc(D->numPlusY*sizeof(double ));
      D->numMinusY=3*nxSub1D*1*3;
      D->minusY = (double *)malloc(D->numMinusY*sizeof(double ));
      break;
    case (Yee-1)*3+2:
      //first 4 is 4 field variables, 2nd 1 is k, 3rd 2 or 3 is layer.
      D->numPlusX=3*nySub2D*1*2;
      D->plusX = (double *)malloc(D->numPlusX*sizeof(double ));
      D->numMinusX=3*nySub2D*1*3;
      D->minusX = (double *)malloc(D->numMinusX*sizeof(double ));
      D->numPlusY=3*nxSub1D*1*2;
      D->plusY = (double *)malloc(D->numPlusY*sizeof(double ));
      D->numMinusY=3*nxSub1D*1*3;
      D->minusY = (double *)malloc(D->numMinusY*sizeof(double ));
      break;
/*
    case 2 :
      //first 3 is 3 field variables, 2nd 1 is k, 3rd 1 is 1 layer.
      D->numPlusYC=3*(D->nx+5)*1*1;
      D->plusYC = (double *)malloc(D->numPlusYC*sizeof(double ));
      D->numMinusYC=3*(D->nx+5)*1*1;
      D->minusYC = (double *)malloc(D->numMinusYC*sizeof(double ));

      //first 6 is 6 field variables, 2nd 1 is k, 3rd 2 or 3 is layer.
      D->numPlusY=6*(D->nx+5)*1*2;
      D->plusY = (double *)malloc(D->numPlusY*sizeof(double ));
      D->numMinusY=6*(D->nx+5)*1*3;
      D->minusY = (double *)malloc(D->numMinusY*sizeof(double ));

      //density : first 1 is 1 variable.
      D->numPlusDenY=1*(D->nx+5)*1*3;
      D->plusDenY = (double *)malloc(D->numPlusDenY*sizeof(double ));
      D->numMinusDenY=1*(D->nx+5)*1*2;
      D->minusDenY = (double *)malloc(D->numMinusDenY*sizeof(double ));
    break;
    case 3 :
      D->numPlusYC=3*(D->nx+5)*(D->nzSub+5)*1;
      D->plusYC = (double *)malloc(D->numPlusYC*sizeof(double ));
      D->numMinusYC=3*(D->nx+5)*(D->nzSub+5)*1;
      D->minusYC = (double *)malloc(D->numMinusYC*sizeof(double ));
      D->numPlusZC=3*(D->nx+5)*(D->nySub+5)*1;
      D->plusZC = (double *)malloc(D->numPlusZC*sizeof(double ));
      D->numMinusZC=3*(D->nx+5)*(D->nySub+5)*1;
      D->minusZC = (double *)malloc(D->numPlusZC*sizeof(double ));

      D->numPlusY=6*(D->nx+5)*(D->nzSub+5)*2;
      D->plusY = (double *)malloc(D->numPlusY*sizeof(double ));
      D->numMinusY=6*(D->nx+5)*(D->nzSub+5)*3;
      D->minusY = (double *)malloc(D->numMinusY*sizeof(double ));
      D->numPlusZ=6*(D->nx+5)*(D->nySub+5)*2;
      D->plusZ = (double *)malloc(D->numPlusZ*sizeof(double ));
      D->numMinusZ=6*(D->nx+5)*(D->nySub+5)*3;
      D->minusZ = (double *)malloc(D->numMinusZ*sizeof(double ));

      //density : first 1 is 1 variable.
      D->numPlusDenY=1*(D->nx+5)*(D->nzSub+5)*3;
      D->plusDenY = (double *)malloc(D->numPlusDenY*sizeof(double ));
      D->numMinusDenY=1*(D->nx+5)*(D->nzSub+5)*2;
      D->minusDenY = (double *)malloc(D->numMinusDenY*sizeof(double ));
      D->numPlusDenZ=1*(D->nx+5)*(D->nySub+5)*3;
      D->plusDenZ = (double *)malloc(D->numPlusDenZ*sizeof(double ));
      D->numMinusDenZ=1*(D->nx+5)*(D->nySub+5)*2;
      D->minusDenZ = (double *)malloc(D->numMinusDenZ*sizeof(double ));
    break;
*/
    }    

    // current J trasffering boundary
    switch(D->dimension)  {
    case 1 :
      D->numPlusXJ=3*nySub2D*1*3;
      D->XplusJ = (double *)malloc(D->numPlusXJ*sizeof(double ));
      D->numMinusXJ=3*nySub2D*1*2;
      D->XminusJ = (double *)malloc(D->numMinusXJ*sizeof(double ));
      break;
    case 2 :
      D->numPlusXJ=3*nySub2D*1*3;
      D->XplusJ = (double *)malloc(D->numPlusXJ*sizeof(double ));
      D->numMinusXJ=3*nySub2D*1*2;
      D->XminusJ = (double *)malloc(D->numMinusXJ*sizeof(double ));
      D->numPlusYJ=3*nxSub1D*1*3;
      D->YplusJ = (double *)malloc(D->numPlusYJ*sizeof(double ));
      D->numMinusYJ=3*nxSub1D*1*2;
      D->YminusJ = (double *)malloc(D->numMinusYJ*sizeof(double ));
      break;
/*
    case 3 :
      numberData=3*3*(D->nx+5)*(D->nzSub+5);
      D->YplusJ=(double *)malloc(numberData*sizeof(double ));
      numberData=2*3*(D->nx+5)*(D->nzSub+5);
      D->YminusJ=(double *)malloc(numberData*sizeof(double ));
      numberData=3*3*(D->nx+5)*(D->nySub+5);
      D->ZplusJ=(double *)malloc(numberData*sizeof(double ));
      numberData=2*3*(D->nx+5)*(D->nySub+5);
      D->ZminusJ=(double *)malloc(numberData*sizeof(double ));
      break;
*/
    default :
      printf("In boundary's share J, what dimension(%d)?\n",D->dimension);
    }

     //Density memory
     D->Rho=memoryAsign(nxSub1D,nySub2D,nzSub3D);

     switch(D->fieldType)  {
     case Split :
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Pl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sr=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Sl=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ExC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->PlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SrC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->SlC=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jx=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->JxOld=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->JyOld=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->JzOld=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
//       D->JxBoost=memoryAsign(nxSub1D,nySub2D,nzSub3D);
//       D->JyBoost=memoryAsign(nxSub1D,nySub2D,nzSub3D);
//       D->JzBoost=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       break;

     case Yee :
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ey=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ez=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->By=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bz=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ByNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BzNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jx=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       break;
     case Pukhov :
       D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ey=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Ez=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->By=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Bz=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BxNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->ByNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->BzNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
       D->Jx=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jy=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       D->Jz=memoryAsignJ(nxSub1D,nySub2D,nzSub3D);
       break;
     default :
       printf("what field type?\n");
     }

     // Particle setting
     nxSub1D=D->nxSub+3;
     nySub2D=1;
     nzSub3D=1;
     startj=0;
     startk=0;
     if(D->dimension>1) {
       nySub2D=D->nySub+3;
       startj=D->jstart-1;
     }
     if(D->dimension>2) { 
       nzSub3D=D->nzSub+3;
       startk=D->kstart-1;
     }

     D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
     for(i=0; i<nxSub1D; i++) {
       D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
       for(j=0; j<nySub2D; j++) 
         D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
     }

     // setting up particle's pointer
     if(D->dimension==1)
     {
       j=k=0;
       //i starts at 0 because of boost frame
       for(i=0; i<D->iend+1; i++)     {
         D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)      {
           D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j][k].head[s]->pt = NULL;
         }
       }
     }
     else if(D->dimension==2)
     {
       k=0;
       for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
         for(j=D->jstart-1; j<D->jend+1; j++)
         {
           D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
           for(s=0; s<D->nSpecies; s++)
           {
             D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
             D->particle[i][j][k].head[s]->pt = NULL;
           }
         }
     }
     else if(D->dimension==3)
     {
       for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
         for(j=D->jstart-1; j<D->jend+1; j++)
           for(k=D->kstart-1; k<D->kend+1; k++)
           {
             D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
             for(s=0; s<D->nSpecies; s++)
             {
               D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
               D->particle[i][j][k].head[s]->pt = NULL;
             }
           }
     }
   

/*
    //Making track Particle memory
    trackStart=D->dumpStep/D->trackSaveStep;
    if(trackStart==0)  D->trackStart=D->dumpStep;
    else 
      D->trackStart=(trackStart+1)*D->trackSaveStep;

    numberData=(D->maxStep-D->trackStart)/D->trackSaveStep+1;
    D->track = (Track **)malloc(D->idNums*sizeof(Track *));
    for(i=0; i<D->idNums; i++)
      D->track[i] = (Track *)malloc(numberData*sizeof(Track ));
    for(i=0; i<D->idNums; i++)
      for(j=0; j<numberData; j++)
      {
        D->track[i][j].x=0.0;
        D->track[i][j].y=0.0;
        D->track[i][j].z=0.0;
        D->track[i][j].px=0.0;
        D->track[i][j].py=0.0;
        D->track[i][j].pz=0.0;
        D->track[i][j].step=0;
        D->track[i][j].id=0;
        D->track[i][j].core=0;
      }
*/
/*
     D->probe = (Probe **)malloc(D->probeNum*sizeof(Probe *));
     for(i=0; i<D->probeNum; i++)
       D->probe[i] = (Probe *)malloc((D->maxStep+1)*sizeof(Probe ));
     for(i=0; i<D->probeNum; i++)
       for(j=0; j<=D->maxStep; j++)
       {
         D->probe[i][j].Ex=0.0;
         D->probe[i][j].Ey=0.0;
         D->probe[i][j].Ez=0.0;
         D->probe[i][j].Bx=0.0;
         D->probe[i][j].By=0.0;
         D->probe[i][j].Bz=0.0;
       }
*/

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

double ***memoryAsignJ(int nx, int ny, int nz)
{
   int i,j,k;
   double ***current;

   current = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)
   {
     current[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       current[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         current[i][j][k]=0.0;
       }

   return current;
}
