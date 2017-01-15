#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"

void absorb2D_UD(Domain *D,double *upr,double *btr,double *upd,double *btd,double y,double upL,double bottomL,double LdU,double LdB,double rr,double rd);
void solveR1D_Split(Domain *D);
void solveL1D_Split(Domain *D);
void solve2DC_Split(Domain *D);
void solve2D_Split(Domain *D);
void Bsolve2D_Pukhov(Domain *D);
void Esolve2D_Pukhov(Domain *D);
void Bsolve2D_Yee(Domain *D);
void Esolve2D_Yee(Domain *D);

void fieldSolve(Domain D,double t)
{
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;
  void MPI_Transfer4F_Xminus();
  void MPI_Transfer4F_Xplus();
  void MPI_Transfer3F_Xminus();
  void MPI_Transfer3F_Xplus();
  void MPI_Transfer6F_Xminus();
  void MPI_Transfer6F_Xplus();
  void MPI_Transfer3F_Yminus();
  void MPI_Transfer3F_Yplus();
  void MPI_Transfer6F_Yminus();
  void MPI_Transfer6F_Yplus();

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch((D.fieldType-1)*3+D.dimension) {

  //1D field
  case (Split-1)*3+1:
    if(D.boostOn==OFF)    {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }
    solveR1D_Split(&D);
    solveL1D_Split(&D);
    if(D.L>1) {
      MPI_Transfer6F_Xplus(&D,D.Ex,D.Pr,D.Sr,D.Bx,D.Pl,D.Sl,1,1,3);
      MPI_Transfer6F_Xminus(&D,D.Ex,D.Pr,D.Sr,D.Bx,D.Pl,D.Sl,1,1,3);
    } else ;
    break;

  //split 2D
  case (Split-1)*3+2:
    if(D.boostOn==OFF)    {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }
    solve2DC_Split(&D);
    if(D.L>1)  {
      MPI_Transfer4F_Xminus(&D,D.ExC,D.BxC,D.PlC,D.SlC,D.nySub+5,1,3);
      MPI_Transfer4F_Xplus(&D,D.ExC,D.BxC,D.PrC,D.SrC,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.SrC,D.SlC,D.ExC,D.nxSub+5,1,3);
      MPI_Transfer3F_Yplus(&D,D.PrC,D.PlC,D.BxC,D.nxSub+5,1,3);
    } else	;

    solve2D_Split(&D);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,1,3);
      MPI_Transfer6F_Xplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
      MPI_Transfer6F_Yplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
    } else	;

    break;

  //3D
  case (Split-1)*3+3:
//    solveField3D_DSX(D);
    
//    if(D->pmlOn==ON)
//    {
//      rankM=myrank%D->M;
//      rankN=(int)(myrank/D->M);
//      if(rankM==D->M-1)
//        absorb3D(D,UP);
//      else	;
////      if(rankM==D->M-1 && rankN==D->N-1)
////        absorb3D(D,UPFRONT);
////     if(rankM==D->M-1 && rankN==0)
////        absorb3D(D,UPBACK);
//      if(rankM==0)
//        absorb3D(D,DOWN);
//      else	;
////      if(rankM==0 && rankN==D->N-1)
////        absorb3D(D,DOWNFRONT);
////      if(rankM==0 && rankN==0)
////        absorb3D(D,DOWNBACK);
//      if(rankN==D->N-1)
//        absorb3D(D,FRONT);
//      else	;
//      if(rankN==0)
//        absorb3D(D,BACK);
//      else	;
//    }
//    else	;
//    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //Yee 2D
  case (Yee-1)*3+2:
    Esolve2D_Yee(&D);
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
//        if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }

    Bsolve2D_Yee(&D);
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Bx,D.By,D.Bz,D.nySub+5,1,3);
      MPI_Transfer3F_Xplus(&D,D.Bx,D.By,D.Bz,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Bx,D.By,D.Bz,D.nxSub+5,1,3);
      MPI_Transfer3F_Yplus(&D,D.Bx,D.By,D.Bz,D.nxSub+5,1,3);
    } else	;
    break;

  //Pukhov 2D
  case (Pukhov-1)*3+2:
    Esolve2D_Pukhov(&D);
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)       {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }  else ;

    Bsolve2D_Pukhov(&D);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
      MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;
    break;

  default:
    printf("what fieldType? and what dimension?\n");
  }
}

void solveR1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;
  nxSub=D->nxSub;

  j=k=0;
  for(i=iend-1; i>=istart; i--)
  {
    D->Ex[i][j][k]+=-2.0*pi*dt*D->Jx[i][j][k];
    D->Pr[i][j][k]=D->Pr[i-1][j][k]-pi*dt*D->Jy[i][j][k];
    D->Sr[i][j][k]=D->Sr[i-1][j][k]-pi*dt*D->Jz[i][j][k];
  }  
}

void solveL1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;
  nxSub=D->nxSub;

  j=k=0;
  for(i=istart; i<iend; i++)
  {
    D->Pl[i][j][k]=D->Pl[i+1][j][k]-pi*dt*D->Jy[i+1][j][k];
    D->Sl[i][j][k]=D->Sl[i+1][j][k]-pi*dt*D->Jz[i+1][j][k];
  }  
}

void Bsolve2D_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,oldBx,oldBy,oldBz,ax,ay,bx,by,x1,x2,y1,y2,x,y;
    int minXSub,minYSub;
    double dtOverdx,dtOverdy,LdR,LdL,LdU,LdB,rr,rd;
    double right1r,right2r,left1r,left2r,upr,upd,btr,btd;
    double right1d,right2d,left1d,left2d,up1d,up2d,bt1d,bt2d;
    double rightL,leftL,upL,bottomL,tmp;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
    ay=ax=0.125*dx/dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=dt/dy;
    dtOverdx=dt/dx;

    minXSub=D->minXSub;
    minYSub=D->minYSub;
    rightL=D->nx+D->minXDomain-D->pmlCellRight;
    leftL=D->minXDomain+D->pmlCellLeft;
    upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
    bottomL=(double)(D->minYDomain+D->pmlCellBottom);
    LdR=D->pmlCellRight;
    LdL=D->pmlCellLeft;
    LdU=D->pmlCellUp;
    LdB=D->pmlCellBottom;
    rr=D->pmlr;
    rd=D->pmld;

    right1r=left1r=right2r=left2r=1.0;
    right1d=left1d=right2d=left2d=1.0;
    upr=btr=upd=btd=1.0;

    //Solving B field
    for(i=istart; i<iend; i++)
    {
      x1=(i-istart)+minXSub;
      x2=(i-istart)+minXSub+0.5;
/*
      right1r=1.0; //damping(x1,rightL,LdR,rr);
      right2r=1.0; //damping(x2,rightL,LdR,rr);
      left1r=1.0;  //damping(x1,leftL,LdL,rr);
      left2r=1.0;  //damping(x2,leftL,LdL,rr);
      right1d=1.0;	//damping(x1,rightL,LdR,rd);
      right2d=1.0;	//damping(x2,rightL,LdR,rd);
      left1d=1.0;	//damping(x1,leftL,LdL,rd);
      left2d=1.0;	//damping(x2,leftL,LdL,rd);
*/
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else	;

        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];

        tmp=right1r*left1r*upr*btr*(-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]));
        D->Bx[i][j][k]=right1d*left1d*upd*btd*(oldBx+tmp);
        tmp=right2r*left2r*upr*btr*(dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]));
        D->By[i][j][k]=right2d*left2d*upd*btd*(oldBy+tmp);
        tmp=right2r*left2r*upr*btr*(dtOverdy*(bx*(D->Ex[i][j+1][k]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]))-dtOverdx*(by*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+ay*(D->Ey[i+1][j+1][k]+D->Ey[i+1][j-1][k]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k])));
        D->Bz[i][j][k]=right2d*left2d*upd*btd*(oldBz+tmp);

        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
    }
}



void Esolve2D_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,ax,ay,bx,by,x1,x2,y1,y2,x,y;
    double oldEx,oldEy,oldEz;
    int minXSub,minYSub;
    double dtOverdx,dtOverdy,LdR,LdL,LdU,LdB,rr,rd;
    double right1r,right2r,left1r,left2r,upr,upd,btr,btd;
    double right1d,right2d,left1d,left2d,up1d,up2d,bt1d,bt2d;
    double rightL,leftL,upL,bottomL,tmp;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    ay=ax=0.125*dx/dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=dt/dy;
    dtOverdx=dt/dx;

    minXSub=D->minXSub;
    minYSub=D->minYSub;
    rightL=D->nx+D->minXDomain-D->pmlCellRight;
    leftL=D->minXDomain+D->pmlCellLeft;
    upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
    bottomL=(double)(D->minYDomain+D->pmlCellBottom);
    LdR=D->pmlCellRight;
    LdL=D->pmlCellLeft;
    LdU=D->pmlCellUp;
    LdB=D->pmlCellBottom;
    rr=D->pmlr;
    rd=D->pmld;

    right1r=left1r=right2r=left2r=1.0;
    right1d=left1d=right2d=left2d=1.0;
    upr=btr=upd=btd=1.0;
//    up1d=bt1d=up2d=bt2d=1.0;

    k=0;
    //Solving E field
    for(i=istart; i<iend; i++)
    {
      x1=(i-istart)+minXSub;
      x2=(i-istart)+minXSub+0.5;
/*
      right1r=1.0;	///damping(x1,rightL,LdR,rr);
      right2r=1.0;	//damping(x2,rightL,LdR,rr);
      left1r=1.0;	//damping(x1,leftL,LdL,rr);
      left2r=1.0;	//damping(x2,leftL,LdL,rr);
      right1d=1.0;	//damping(x1,rightL,LdR,rd);
      right2d=1.0;	//damping(x2,rightL,LdR,rd);
      left1d=1.0;	//damping(x1,leftL,LdL,rd);
      left2d=1.0;	//damping(x2,leftL,LdL,rd);
*/
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else	;

        oldEx=D->Ex[i][j][k];
        oldEy=D->Ey[i][j][k];
        oldEz=D->Ez[i][j][k];


        tmp=right2r*left2r*upr*btr*(dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-2.0*pi*dt*D->Jx[i][j][k]);
        D->Ex[i][j][k]=right2d*left2d*upd*btd*(oldEx+tmp);
        tmp=right1r*left1r*upr*btr*(-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-2.0*pi*dt*D->Jy[i][j][k]);
        D->Ey[i][j][k]=right1d*left1d*upd*btd*(oldEy+tmp);
        tmp=right1r*left1r*upr*btr*(dtOverdx*(by*(D->By[i][j][k]-D->By[i-1][j][k])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i-1][j+1][k]-D->By[i-1][j-1][k]))-dtOverdy*(bx*(D->Bx[i][j][k]-D->Bx[i][j-1][k])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i-1][j-1][k]))-2.0*pi*dt*D->Jz[i][j][k]);
        D->Ez[i][j][k]=right1d*left1d*upd*btd*(oldEz+tmp);
      }
    }
}

void Bsolve2D_Yee(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,oldBx,oldBy,oldBz;

    dx=D->dx;    dy=D->dy;    dt=D->dt; 
    nxSub=D->nxSub;    nySub=D->nySub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    k=0;
    //Solving B field
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];
        D->Bx[i][j][k]+=-dt/dy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->By[i][j][k]+=dt/dx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->Bz[i][j][k]+=-dt/dx*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+dt/dy*(D->Ex[i][j+1][k]-D->Ex[i][j][k]);
        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
}

void Esolve2D_Yee(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt;

    dx=D->dx;    dy=D->dy;    dz=D->dz;    dt=D->dt;
    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    k=0;
    //Solving E field
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        D->Ex[i][j][k]+=dt/dy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-2*pi*dt*D->Jx[i][j][k];
        D->Ey[i][j][k]+=-dt/dx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-2*pi*dt*D->Jy[i][j][k];
        D->Ez[i][j][k]+=dt/dx*(D->By[i][j][k]-D->By[i-1][j][k])-dt/dy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-2*pi*dt*D->Jz[i][j][k];
      }
}
//lala

void solve2DC_Split(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,minXSub,minYSub,y;
    double dtOverdx,dtOverdy,LdR,LdL,LdU,LdB,rr,rd;
    double right1r,right2r,left1r,left2r,upr,upd,btr,btd;
    double right1d,right2d,left1d,left2d,up1d,up2d,bt1d,bt2d;
    double rightL,leftL,upL,bottomL,tmp;

    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    dx=D->dx;    dy=D->dy;    dz=D->dz;    dt=D->dt;
    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend; k=0;

    minXSub=D->minXSub; minYSub=D->minYSub;
    rightL=D->nx+D->minXDomain-D->pmlCellRight;
    leftL=D->minXDomain+D->pmlCellLeft;
    upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
    bottomL=(double)(D->minYDomain+D->pmlCellBottom);
    LdR=D->pmlCellRight; LdL=D->pmlCellLeft;
    LdU=D->pmlCellUp;    LdB=D->pmlCellBottom;
    rr=D->pmlr;    rd=D->pmld;

    right1r=left1r=right2r=left2r=1.0;
    right1d=left1d=right2d=left2d=1.0;
    upr=btr=upd=btd=1.0;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else    ;

        tmp=right1r*left1r*upr*btr*(dt/dy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])-pi*dt*(D->JxOld[i][j][k]+D->Jx[i][j][k]));
        D->ExC[i][j][k]=right1d*left1d*upd*btd*(D->ExC[i][j][k]+tmp);
        tmp=right2r*left2r*upr*btr*(-dt/dy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]));
        D->BxC[i][j][k]=right1d*left1d*upd*btd*(D->BxC[i][j][k]+tmp);

        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->Ex[i+1][j+1][k]+D->Ex[i][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i][j][k])-0.5*pi*dt*(D->JyOld[i+1][j][k]+D->Jy[i+1][j][k]));
        D->PlC[i][j][k]=right2r*left2r*upd*btd*(D->PlC[i+1][j][k]+tmp);

        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i][j-1][k])-0.5*pi*dt*(D->JzOld[i+1][j][k]+D->Jz[i+1][j][k]));
        D->SlC[i][j][k]=right2r*left2r*upd*btd*(D->SlC[i+1][j][k]+tmp);
      }	

    // PrC,SrC
    for(i=iend-1; i>=istart; i--)
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else    ;
        tmp=right2r*left2r*upr*btr*(0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]));
        D->PrC[i][j][k]=right2r*left2r*upd*btd*(D->PrC[i-1][j][k]+tmp);
        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]));
        D->SrC[i][j][k]=right2r*left2r*upd*btd*(D->SrC[i-1][j][k]+tmp);
      }	

}

void solve2D_Split(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,minXSub,minYSub,y;
    double dtOverdx,dtOverdy,LdR,LdL,LdU,LdB,rr,rd;
    double right1r,right2r,left1r,left2r,upr,upd,btr,btd;
    double right1d,right2d,left1d,left2d,up1d,up2d,bt1d,bt2d;
    double rightL,leftL,upL,bottomL,tmp;

    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    dx=D->dx;    dy=D->dy;    dz=D->dz;    dt=D->dt;
    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend; k=0;

    minXSub=D->minXSub; minYSub=D->minYSub;
    rightL=D->nx+D->minXDomain-D->pmlCellRight;
    leftL=D->minXDomain+D->pmlCellLeft;
    upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
    bottomL=(double)(D->minYDomain+D->pmlCellBottom);
    LdR=D->pmlCellRight; LdL=D->pmlCellLeft;
    LdU=D->pmlCellUp;    LdB=D->pmlCellBottom;
    rr=D->pmlr;    rd=D->pmld;

    right1r=left1r=right2r=left2r=1.0;
    right1d=left1d=right2d=left2d=1.0;
    upr=btr=upd=btd=1.0;

    // Pl,E1,Sl,B1
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else    ;

        tmp=right1r*left1r*upr*btr*(dt/dy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])-2*pi*dt*D->Jx[i][j][k]);
        D->Ex[i][j][k]=right1d*left1d*upd*btd*(D->Ex[i][j][k]+tmp);
        tmp=right2r*left2r*upr*btr*(-dt/dy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]));
        D->Bx[i][j][k]=right1d*left1d*upd*btd*(D->Bx[i][j][k]+tmp);

        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->ExC[i+1][j+1][k]+D->ExC[i][j+1][k]-D->ExC[i+1][j][k]-D->ExC[i][j][k])-pi*dt*D->Jy[i+1][j][k]);
        D->Pl[i][j][k]=right2r*left2r*upd*btd*(D->Pl[i+1][j][k]+tmp);

        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->BxC[i+1][j][k]+D->BxC[i][j][k]-D->BxC[i+1][j-1][k]-D->BxC[i][j-1][k])-pi*dt*D->Jz[i+1][j][k]);
        D->Sl[i][j][k]=right2r*left2r*upd*btd*(D->Sl[i+1][j][k]+tmp);
      }	

    // Pr,Sr
    for(i=iend-1; i>=istart; i--)
      for(j=jstart; j<jend; j++)
      {
        y=(j-jstart)+minYSub;
        if(D->pmlOn==ON)
          absorb2D_UD(D,&upr,&btr,&upd,&btd,y,upL,bottomL,LdU,LdB,rr,rd);
        else    ;
        tmp=right2r*left2r*upr*btr*(0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jy[i][j][k]);
        D->Pr[i][j][k]=right2r*left2r*upd*btd*(D->Pr[i-1][j][k]+tmp);
        tmp=right2r*left2r*upr*btr*(-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-pi*dt*D->Jz[i][j][k]);
        D->Sr[i][j][k]=right2r*left2r*upd*btd*(D->Sr[i-1][j][k]+tmp);
      }	
}

/*
void solveField3DC_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub; 
    double dx,dy,dz,dt,preB1C;
    double nowPr,nowSr,prevPr,prevSr;
    double nowPrC,nowSrC,prevPrC,prevSrC;
    int myrank,rank,rankM,rankN;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

//    if(D->pmlOn==ON)
//    {
//      rankM=myrank%D->M;
//      if(rankM==D->M-1)
//        jend=jend-D->pmlCell;
//      if(rankM==0)
//        jstart=jstart+D->pmlCell;
//      rankN=(int)(myrank/D->M);
//      if(rankN==D->N-1)
//        kend=kend-D->pmlCell;
//      if(rankN==0)
//        kstart=kstart+D->pmlCell;
//    }

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(k=kstart; k<kend; k++)
      for(j=jstart; j<jend; j++)
      {
        nowPrC=D->PrC[istart-1][j][k];
        nowSrC=D->SrC[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]+=dt/dy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])+dt/dz*(D->Sr[i][j][k]-D->Sr[i][j][k-1]-D->Sl[i][j][k]+D->Sl[i][j][k-1])-pi*dt*(D->JxOld[i][j][k]+D->Jx[i][j][k]);
          D->BxC[i][j][k]+=dt/dz*(D->Pr[i][j][k+1]-D->Pr[i][j][k]+D->Pl[i][j][k+1]-D->Pl[i][j][k])-dt/dy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
          prevPrC=nowPrC;
          nowPrC=D->PrC[i][j][k];
          D->PrC[i][j][k]=prevPrC+0.25*dt/dz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])+0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          D->PlC[i-1][j][k]=D->PlC[i][j][k]+0.25*dt/dz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])-0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
          prevSrC=nowSrC;
          nowSrC=D->SrC[i][j][k];
          D->SrC[i][j][k]=prevSrC-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])+0.25*dt/dz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
          D->SlC[i-1][j][k]=D->SlC[i][j][k]-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-0.25*dt/dz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
        }	//End of i
      }		//End of j,k
}

void solveField3D_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,preB1C;
    double nowPr,nowSr,prevPr,prevSr;
    double nowPrC,nowSrC,prevPrC,prevSrC;
    int myrank,rankM,rankN;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

//    if(D->pmlOn==ON)
//    {
//      rankM=myrank%D->M;
//      if(rankM==D->M-1)
//        jend=jend-D->pmlCell;
//      if(rankM==0)
//        jstart=jstart+D->pmlCell;
//      rankN=(int)(myrank/D->M);
//      if(rankN==D->N-1)
//        kend=kend-D->pmlCell;
//      if(rankN==0)
//        kstart=kstart+D->pmlCell;
//    }

    // Pr,Pl,E1,Sr,Sl,B1
    for(k=kstart; k<kend; k++)
      for(j=jstart; j<jend; j++)
      {
        nowPr=D->Pr[istart-1][j][k];
        nowSr=D->Sr[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]+=dt/dy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])+dt/dz*(D->SrC[i][j][k]-D->SrC[i][j][k-1]-D->SlC[i][j][k]+D->SlC[i][j][k-1])-2*pi*dt*D->Jx[i][j][k];
          D->Bx[i][j][k]+=dt/dz*(D->PrC[i][j][k+1]-D->PrC[i][j][k]+D->PlC[i][j][k+1]-D->PlC[i][j][k])-dt/dy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]);
          prevPr=nowPr;
          nowPr=D->Pr[i][j][k];
          D->Pr[i][j][k]=prevPr+0.25*dt/dz*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j][k-1]-D->BxC[i-1][j][k-1])+0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jy[i][j][k];
          D->Pl[i-1][j][k]=D->Pl[i][j][k]+0.25*dt/dz*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j][k-1]-D->BxC[i-1][j][k-1])-0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jy[i][j][k];
          prevSr=nowSr;
          nowSr=D->Sr[i][j][k];
          D->Sr[i][j][k]=prevSr-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])+0.25*dt/dz*(D->ExC[i][j][k+1]+D->ExC[i-1][j][k+1]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jz[i][j][k];
          D->Sl[i-1][j][k]=D->Sl[i][j][k]-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-0.25*dt/dz*(D->ExC[i][j][k+1]+D->ExC[i-1][j][k+1]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*D->Jz[i][j][k];
        }	//End of i
      }		//End of j,k
}
*/
