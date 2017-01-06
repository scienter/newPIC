#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

double dampingU(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x<=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

double dampingD(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x>=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

void absorb2D_UD(Domain *D,double *upr,double *btr,double *upd,double *btd,double y,double upL,double bottomL,double LdU,double LdB,double rr,double rd)
{
  double tmp;

  if(y>=upL)  {
    tmp=(y-upL)/LdU;
    *upr=1.0-rr*rr*tmp*tmp;
    *btr=1.0;
    *upd=1.0-rd*rd*tmp*tmp;
    *btd=1.0;
  } else if(y<=bottomL)  {
    tmp=(bottomL-y)/LdB;
    *upr=1.0;
    *btr=1.0-rr*rr*tmp*tmp;
    *upd=1.0;
    *btd=1.0-rd*rd*tmp*tmp;
  } else  {  
    *upr=1.0;
    *btr=1.0;
    *upd=1.0;
    *btd=1.0;
  }
}
/*
void absorb3DC(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,kinit,kfinal,pmlCell;
  double prevSrC,nowSrC,prevPrC,nowPrC,tmp;
  double dx,dt,dy,dz,rr,y,z;
  double damp11,damp12,damp21,damp22,damp31,damp32,damp41,damp42;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  dx=D->dx;
  dy=D->dy;
  dz=D->dz;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
//  rd=D->pmld;

  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        }
      }
    }
    break;
  case UPFRONT :    
    jinit=jend-D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case UPBACK :    
    jinit=jend-D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWNBACK :
    jfinal=jstart+D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWNFRONT :
    jfinal=jstart+D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case FRONT :    
    kinit=kend-D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kinit; k<kend; k++)
    {
      z=k-kinit+1;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp11*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp11*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp21*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp21*D->SlC[i][j][k];
        }
      }
    }
    break;
  case BACK :
    kfinal=kstart+D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kstart; k<kfinal; k++)
    {
      z=kfinal-k;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp11*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp11*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp21*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp21*D->SlC[i][j][k];
        }
      }
    }
    break;
  }
}

void absorb3D(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,kinit,kfinal,pmlCell;
  double prevSr,nowSr,prevPr,nowPr;
  double dx,dt,dy,dz,rr,y,z;
  double damp11,damp12,damp21,damp22,damp31,damp32,damp41,damp42;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  dx=D->dx;
  dy=D->dy;
  dz=D->dz;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
//  rd=D->pmld;

  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case UPFRONT :    
    jinit=jend-D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case UPBACK :    
    jinit=jend-D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Pr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        }
      }
    } 
    break;
  case DOWNBACK :
    jfinal=jstart+D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }
    }
    break;
  case DOWNFRONT :
    jfinal=jstart+D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }
    }
    break;
  case FRONT :    
    kinit=kend-D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kinit; k<kend; k++)
    {
      z=k-kinit+1;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp11*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp11*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp21*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp21*D->Sl[i][j][k];
        }
      }
    }
    break;
  case BACK :
    kfinal=kstart+D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kstart; k<kfinal; k++)
    {
      z=kfinal-k;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp11*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp11*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp21*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp21*D->Sl[i][j][k];
        }
      }
    }
    break;
  }
}

void absorb2DC(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,pmlCell;
  double prevSrC,nowSrC,prevPrC,nowPrC;
  double dx,dt,dy,dz,rr,rd,y,damp11,damp12,damp21,damp22;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  dx=D->dx;
  dy=D->dy;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
  rd=D->pmld;

  
  switch(position) {
  case UP :    
    jinit=jend-pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
        D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
        D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
        D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
        D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
        D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
        D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
        D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
        D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
        D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
      }
    }
    break;
  }
}

void absorb2D(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,pmlCell;
  double prevSr,nowSr,prevPr,nowPr;
  double dx,dt,dy,dz,rr,y,damp11,damp12,damp21,damp22;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  dx=D->dx;
  dy=D->dy;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;

  k=0;
  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
        D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
        D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
        D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
        D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
        D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
        D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
        D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
        D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
        D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
      }
    }
    break;
  }

}
*/
