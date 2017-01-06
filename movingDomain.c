#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void movingDomain(Domain *D,int iteration)
{
  int shiftDuration;
  void movingDomain2D_Pukhov();

  switch((D->fieldType-1)*3+D->dimension)  {
  case ((Pukhov-1)*3+1) :
//    movingDomain1D_DSX(D);
    break;
  case ((Yee-1)*3+2) :
    movingDomain2D_Pukhov(D,iteration);
    break;
  case ((Pukhov-1)*3+2) :
    shiftDuration=(int)(1.0/(1.0-D->dtRatio));
    if((iteration-D->nx)%shiftDuration!=0) 
      movingDomain2D_Pukhov(D);
    else	;
    break;
  case ((1-1)*3+3) :
//    movingDomain3D_DSX(D);
    break;
  default:
    printf("what fieldType?(%d) and what dimension?(%d)\n",D->fieldType,D->dimension);
  }
}
/*
void movingDomain1D_DSX(Domain *D)
{
    int i,j,k,istart,iend,s;
    double x;
    ptclList *p,*tmp,*New;      
    Particle ***particle;
    particle=D->particle;
   
    istart=D->istart;
    iend=D->iend;

    j=k=0;
    for(i=istart-1; i<=iend; i++)
      {
        D->ExC[i][j][k]=D->ExC[i+1][j][k];
        D->BxC[i][j][k]=D->BxC[i+1][j][k];
        D->PrC[i][j][k]=D->PrC[i+1][j][k];
        D->PlC[i][j][k]=D->PlC[i+1][j][k];
        D->SrC[i][j][k]=D->SrC[i+1][j][k];
        D->SlC[i][j][k]=D->SlC[i+1][j][k];
        D->Ex[i][j][k]=D->Ex[i+1][j][k];
        D->Bx[i][j][k]=D->Bx[i+1][j][k];
        D->Pr[i][j][k]=D->Pr[i+1][j][k];
        D->Pl[i][j][k]=D->Pl[i+1][j][k];
        D->Sr[i][j][k]=D->Sr[i+1][j][k];
        D->Sl[i][j][k]=D->Sl[i+1][j][k];
        D->JxOld[i][j][k]=D->JxOld[i+1][j][k];
        D->JyOld[i][j][k]=D->JyOld[i+1][j][k];
        D->JzOld[i][j][k]=D->JzOld[i+1][j][k];
        D->Jx[i][j][k]=D->Jx[i+1][j][k];
        D->Jy[i][j][k]=D->Jy[i+1][j][k];
        D->Jz[i][j][k]=D->Jz[i+1][j][k];
      }     
 
      for(i=istart; i<iend; i++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              p->x-=1.0;  
              p->oldX-=1.0;
              p=p->next;
            }
          } 
        }	

//    for(i=istart-1; i<=iend; i++)
//      for(j=0; j<jend+3; j++)
//      {
//        D->JxBoost[i][j][k]=D->JxBoost[i+1][j][k];
//        D->JyBoost[i][j][k]=D->JyBoost[i+1][j][k];
//        D->JzBoost[i][j][k]=D->JzBoost[i+1][j][k];
//      }     
//    i=iend+1;
//    for(j=0; j<jend+3; j++)
//    {
//      D->JxBoost[i][j][k]=0.0;
//      D->JyBoost[i][j][k]=0.0;
//      D->JzBoost[i][j][k]=0.0;
//    }     

    D->minXSub+=1;
    D->maxXSub+=1;
}
*/  
void movingDomain2D_Pukhov(Domain *D)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,shiftDuration;
    double x;
    ptclList *p,*tmp,*New;      
    Particle ***particle;
    particle=D->particle;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    k=0;
    for(i=istart-1; i<iend; i++)
      for(j=0; j<jend+3; j++)
      {
        D->Ex[i][j][k]=D->Ex[i+1][j][k];
        D->Ey[i][j][k]=D->Ey[i+1][j][k];
        D->Ez[i][j][k]=D->Ez[i+1][j][k];
        D->Bx[i][j][k]=D->Bx[i+1][j][k];
        D->By[i][j][k]=D->By[i+1][j][k];
        D->Bz[i][j][k]=D->Bz[i+1][j][k];
        D->Jx[i][j][k]=D->Jx[i+1][j][k];
        D->Jy[i][j][k]=D->Jy[i+1][j][k];
        D->Jz[i][j][k]=D->Jz[i+1][j][k];
      }    

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)  {
            p->x-=1.0;  
            p->oldX-=1.0;
            p=p->next;
          }
        }	

    D->minXSub+=1;
    D->maxXSub+=1;
    D->minXDomain+=1;
}  
/*     
void movingDomain3D_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s;
    double x;
    ptclList *p,*tmp,*New;      
    Particle ***particle;
    particle=D->particle;
   
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    for(i=istart-1; i<iend; i++)
      for(j=0; j<jend+3; j++)
        for(k=0; k<kend+3; k++)
        {
          D->ExC[i][j][k]=D->ExC[i+1][j][k];
          D->BxC[i][j][k]=D->BxC[i+1][j][k];
          D->PrC[i][j][k]=D->PrC[i+1][j][k];
          D->PlC[i][j][k]=D->PlC[i+1][j][k];
          D->SrC[i][j][k]=D->SrC[i+1][j][k];
          D->SlC[i][j][k]=D->SlC[i+1][j][k];
          D->Ex[i][j][k]=D->Ex[i+1][j][k];
          D->Bx[i][j][k]=D->Bx[i+1][j][k];
          D->Pr[i][j][k]=D->Pr[i+1][j][k];
          D->Pl[i][j][k]=D->Pl[i+1][j][k];
          D->Sr[i][j][k]=D->Sr[i+1][j][k];
          D->Sl[i][j][k]=D->Sl[i+1][j][k];
          D->JxOld[i][j][k]=D->JxOld[i+1][j][k];
          D->JyOld[i][j][k]=D->JyOld[i+1][j][k];
          D->JzOld[i][j][k]=D->JzOld[i+1][j][k];
          D->Jx[i][j][k]=D->Jx[i+1][j][k];
          D->Jy[i][j][k]=D->Jy[i+1][j][k];
          D->Jz[i][j][k]=D->Jz[i+1][j][k];
        }
 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            for(s=0; s<D->nSpecies; s++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                p->x-=1.0;  
                p->oldX-=1.0;
                p=p->next;
              }
            } 
          }	
    D->minXSub+=1;
    D->maxXSub+=1;
}  
*/    
