#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

void cleanMemory(Domain *D)
{
    Particle ***particle;
    particle=D->particle;
    LoadList *LL,*tmpLL;
    LaserList *L, *tmpL;
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    int s,nxSub,nySub,nzSub;
    double *plusX;
    void deleteField();
    ptclList *p,*tmp;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    //remove field share
    switch((D->fieldType-1)*3+D->dimension) {
    case 1 :
      break;
    case (Pukhov-1)*3+2:
      free(D->plusX);
      free(D->minusX);
      free(D->plusY);
      free(D->minusY);
      free(D->XplusJ);
      free(D->XminusJ);
      free(D->YplusJ);
      free(D->YminusJ);
      break;
    }

    //remove particles
    if(D->dimension==2)
    {
      k=0;
      for(i=0; i<iend+1; i++)
        for(j=jstart-1; j<jend+1; j++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {	
              tmp=p->next;
              particle[i][j][k].head[s]->pt=tmp; 
              p->next=NULL;
              free(p);
              p=particle[i][j][k].head[s]->pt;
            }
            free(particle[i][j][k].head[s]);
          }
          free(particle[i][j][k].head);
        }

      for(i=0; i<iend+1; i++)
      {
        for(j=jstart-1; j<jend+1; j++)
        {
          free(D->particle[i][j]);
        }
        free(D->particle[i]);
      }
      free(particle);

    }	
	
    else if(D->dimension==3)
    {
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
          for(k=kstart-1; k<=kend; k++)
          {
            for(s=0; s<D->nSpecies; s++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {	
                tmp=p->next;
                particle[i][j][k].head[s]->pt=tmp; 
                p->next=NULL;
                free(p);
                p=particle[i][j][k].head[s]->pt;
              }
              free(particle[i][j][k].head[s]);
            }
            free(particle[i][j][k].head);
          }
    }		
      

      
    LL=D->loadList;
    while(LL->next)
    {
      switch (LL->type)  {
      case Polygon :
        if(LL->xnodes>0)
        {
          free(LL->xpoint);	
          free(LL->xn);	
        }
        if(D->dimension>1)
        {
          if(LL->ynodes>0)
          {
            free(LL->ypoint);
            free(LL->yn);
          }
        }
        if(D->dimension>2)
        {
          if(LL->znodes>0)
          {
            free(LL->zpoint);
            free(LL->zn);
          }
        }
        break;
      case Defined :
        if(LL->numDefined>0)
        {
          if(D->dimension>1)
          {
            free(LL->xPosition);
            free(LL->yPosition);
          }
          if(D->dimension>2)
            free(LL->zPosition);
          if(D->dimension>=2)
            for(i=0; i<LL->numDefined; i++)
              free(LL->define[i]);
        }
        break;
      }

      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;
      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);

    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);

    //remove field
    nxSub=D->nxSub+5;
    nySub=1;
    nzSub=1;
    if(D->dimension>1)  
      nySub=D->nySub+5;
    if(D->dimension>2)  
      nzSub=D->nzSub+5;

    if(D->fieldType==Pukhov)
    {
      deleteField(D->Ex,nxSub,nySub,nzSub);
      deleteField(D->Ey,nxSub,nySub,nzSub);
      deleteField(D->Ez,nxSub,nySub,nzSub);
      deleteField(D->Bx,nxSub,nySub,nzSub);
      deleteField(D->By,nxSub,nySub,nzSub);
      deleteField(D->Bz,nxSub,nySub,nzSub);
      deleteField(D->BxNow,nxSub,nySub,nzSub);
      deleteField(D->ByNow,nxSub,nySub,nzSub);
      deleteField(D->BzNow,nxSub,nySub,nzSub);
      deleteField(D->Jx,nxSub,nySub,nzSub);
      deleteField(D->Jy,nxSub,nySub,nzSub);
      deleteField(D->Jz,nxSub,nySub,nzSub);
    }
/*
    else if(D->fieldType==Split)
    {

      deleteField(D->Ex,nxSub,nySub,nzSub);
      deleteField(D->Pr,nxSub,nySub,nzSub);
      deleteField(D->Pl,nxSub,nySub,nzSub);
      deleteField(D->Bx,nxSub,nySub,nzSub);
      deleteField(D->Sr,nxSub,nySub,nzSub);
      deleteField(D->Sl,nxSub,nySub,nzSub);
      deleteField(D->ExC,nxSub,nySub,nzSub);
      deleteField(D->PrC,nxSub,nySub,nzSub);
      deleteField(D->PlC,nxSub,nySub,nzSub);
      deleteField(D->BxC,nxSub,nySub,nzSub);
      deleteField(D->SrC,nxSub,nySub,nzSub);
      deleteField(D->SlC,nxSub,nySub,nzSub);
      deleteField(D->Jx,nxSub,nySub,nzSub);
      deleteField(D->Jy,nxSub,nySub,nzSub);
      deleteField(D->Jz,nxSub,nySub,nzSub);
      deleteField(D->JxOld,nxSub,nySub,nzSub);
      deleteField(D->JyOld,nxSub,nySub,nzSub);
      deleteField(D->JzOld,nxSub,nySub,nzSub);

    }

    //remove track
    if(D->tracking==ON)
    {
      if(D->idNums>0)
      {
        for(i=0; i<D->idNums; i++)
         free(D->track[i]);
        free(D->track);
        free(D->trackID);
        free(D->trackCore);
        free(D->trackS);
      }
    }

    //remove probe
    if(D->probeNum>0)
    {
      for(i=0; i<D->probeNum; i++)
       free(D->probe[i]);
      free(D->probe);
      free(D->probeX);
      free(D->probeY);
    }


    //remove boost field
    Boost **boost;
    boost=D->boost;

    for(i=0; i<D->nxSub+5; i++)
      free(boost[i]);
    free(boost);
*/
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
