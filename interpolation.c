#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

double calBi2D(double ***field,int i,int j,double x,double y);

void interpolation(Domain *D,External *Ext)
{
  void interpolation1D_Split_1st();
  void interpolation2D_Pukhov_1st();
  void interpolation2D_Pukhov_2nd();
//  void interpolation2D_DSX_2nd();
//  void interpolation3D_DSX_1st();
//  void interpolation3D_DSX_2nd();
  void MPI_Transfer3F_Xminus();
  void MPI_Transfer3F_Xplus();
  void MPI_TransferF_Pukhov_Yminus();
  void MPI_TransferF_Pukhov_Yplus();

  switch((D->interpolationType-1)*3+D->dimension)  {
  //1D 1st
  case ((1-1)*3+1) :
    interpolation1D_Split_1st(D,Ext);
    break;
  case ((1-1)*3+2) :
    if(D->L>1)  {
      MPI_Transfer3F_Xminus(D,D->BxNow,D->ByNow,D->BzNow,D->nySub+5,1,3);
      MPI_Transfer3F_Xplus(D,D->BxNow,D->ByNow,D->BzNow,D->nySub+5,1,3);
    } else	;
    if(D->M>1)  {
      MPI_TransferF_Pukhov_Yminus(D,D->BxNow,D->ByNow,D->BzNow,D->nxSub+5,1,3);
      MPI_TransferF_Pukhov_Yplus(D,D->BxNow,D->ByNow,D->BzNow,D->nxSub+5,1,3);
    } else	;

    interpolation2D_Pukhov_1st(D,Ext);
    break;
  case ((2-1)*3+2) :
    if(D->L>1)  {
      MPI_Transfer3F_Xminus(D,D->BxNow,D->ByNow,D->BzNow,D->nySub+5,1,3);
      MPI_Transfer3F_Xplus(D,D->BxNow,D->ByNow,D->BzNow,D->nySub+5,1,3);
    }
    if(D->M>1)  {
      MPI_TransferF_Pukhov_Yminus(D,D->BxNow,D->ByNow,D->BzNow,D->nxSub+5,1,3);
      MPI_TransferF_Pukhov_Yplus(D,D->BxNow,D->ByNow,D->BzNow,D->nxSub+5,1,3);
    }

    interpolation2D_Pukhov_2nd(D,Ext);
    break;
/*
  case ((1-1)*3+2) :
    interpolation2D_DSX_1st(D,Ext);
    break;
  case ((2-1)*3+2) :
    interpolation2D_DSX_2nd(D,Ext);
    break;
  case ((1-1)*3+3) :
    interpolation3D_DSX_1st(D,Ext);
    break;
  case ((2-1)*3+3) :
    interpolation3D_DSX_2nd(D,Ext);
    break;
*/
  default :
    printf("In interpolation, what interpolationType(%d)? and what dimension(%d)?\n",D->interpolationType,D->dimension);
  }
}

void interpolation2D_Pukhov_1st(Domain *D,External *Ext)
{
   int i,j,k,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s,cnt;
   double Bx,By,Bz,Ex,Ey,Ez,extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           i1=(int)(i+x+0.5);
           j1=(int)(j+y+0.5);
           x1=x+0.5-((int)(x+0.5));
           y1=y+0.5-((int)(y+0.5));

           Bx=(1.0-x)*(1.0-y1)*D->BxNow[i][j1-1][k1]
             +     x *(1.0-y1)*D->BxNow[i+1][j1-1][k1]
             +(1.0-x)*     y1 *D->BxNow[i][j1][k1]
             +     x *     y1 *D->BxNow[i+1][j1][k1];
           By=(1.0-x1)*(1.0-y)*D->ByNow[i1-1][j][k1]
             +     x1 *(1.0-y)*D->ByNow[i1][j][k1]
             +(1.0-x1)*     y *D->ByNow[i1-1][j+1][k1]
             +     x1 *     y *D->ByNow[i1][j+1][k1];
           Bz=(1.0-x1)*(1.0-y1)*D->BzNow[i1-1][j1-1][k1]
             +     x1 *(1.0-y1)*D->BzNow[i1][j1-1][k1]
             +(1.0-x1)*     y1 *D->BzNow[i1-1][j1][k1]
             +     x1 *     y1 *D->BzNow[i1][j1][k1];
           Ex=(1.0-x1)*(1.0-y)*D->Ex[i1-1][j][k1]
             +     x1 *(1.0-y)*D->Ex[i1][j][k1]
             +(1.0-x1)*     y *D->Ex[i1-1][j+1][k1]
             +     x1 *     y *D->Ex[i1][j+1][k1];
           Ey=(1.0-x)*(1.0-y1)*D->Ey[i][j1-1][k1]
             +     x *(1.0-y1)*D->Ey[i+1][j1-1][k1]
             +(1.0-x)*     y1 *D->Ey[i][j1][k1]
             +     x *     y1 *D->Ey[i+1][j1][k1];
           Ez=(1.0-x)*(1.0-y)*D->Ez[i][j][k1]
             +     x *(1.0-y)*D->Ez[i+1][j][k1]
             +(1.0-x)*     y *D->Ez[i][j+1][k1]
             +     x *     y *D->Ez[i+1][j+1][k1];

           p->E1=Ex+extE1; p->E2=Ey+extE2; p->E3=Ez+extE3;
           p->B1=Bx+extB1; p->B2=By+extB2; p->B3=Bz+extB3;

           p=p->next;
           cnt++;
         }
       }		//for(s)        
     }		   //for(i,j)
}


void interpolation2D_Pukhov_2nd(Domain *D,External *Ext)
{
   int s,i,j,k,ii,jj,i1,j1,k1,istart,iend,jstart,jend,kstart,kend;
   double E1,E2,E3,B1,B2,B3,extE1,extE2,extE3,extB1,extB2,extB3;
   double x,y,z,xx,yy,zz,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           //edge
           ii=(int)(i+x+0.5);
           jj=(int)(j+y+0.5);
           xx=1.0+x-((int)(0.5+x));
           yy=1.0+y-((int)(0.5+y));
           //side
           i1=i;
           j1=j;
           x1=0.5+x;
           y1=0.5+y;

           E1=calBi2D(D->Ex,i1,jj,x1,yy);
           E2=calBi2D(D->Ey,ii,j1,xx,y1);
           E3=calBi2D(D->Ez,ii,jj,xx,yy);
           B1=calBi2D(D->BxNow,ii,j1,xx,y1);
           B2=calBi2D(D->ByNow,i1,jj,x1,yy);
           B3=calBi2D(D->BzNow,i1,j1,x1,y1);

           p->E1=E1+extE1; p->E2=E2+extE2; p->E3=E3+extE3;
           p->B1=B1+extB1; p->B2=B2+extB2; p->B3=B3+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

double calBi2D(double ***field,int i,int j,double x,double y)
{
  int n;
  double y1,y2,y3,a,b,c,result,data[3];

  for(n=0; n<3; n++)  {
    y1=field[i-1+n][j-1][0];
    y2=field[i-1+n][j][0];
    y3=field[i-1+n][j+1][0];
    a=0.5*(y1+y3)-y2;
    b=2.0*y2-1.5*y1-0.5*y3;
    c=y1;
    data[n]=a*y*y+b*y+c;
  }
  y1=data[0];
  y2=data[1];
  y3=data[2];
  a=0.5*(y1+y3)-y2;
  b=2.0*y2-1.5*y1-0.5*y3;
  c=y1;
  result=a*x*x+b*x+c;

  return result;
}

void interpolation1D_Split_1st(Domain *D,External *Ext)
{
   int i,j,k,i1,istart,iend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3,x,x1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;
   
   j=k=0;
   for(i=istart; i<iend; i++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x; // y=p->y; // z=p->z;
           i1=(int)(i+x+0.5);
           x1=x+0.5-((int)(x+0.5));

           B1=0;
           E1=(1-x1)*D->Ex[i1-1][j][k]
             +    x1*D->Ex[i1][j][k];
           Pr=(1-x1)*D->Pr[i1-1][j][k]
             +    x1*D->Pr[i1][j][k];
           Pl=(1-x1)*D->Pl[i1-1][j][k]
             +    x1*D->Pl[i1][j][k];
           Sr=(1-x1)*D->Sr[i1-1][j][k]
             +    x1*D->Sr[i1][j][k];
           Sl=(1-x1)*D->Sl[i1-1][j][k]
             +    x1*D->Sl[i1][j][k];

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
           cnt++;
         }
       }		//for(s)        
     }		   //for(i,j)
}

/*
void interpolation2D_DSX_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,k,ii,jj,kk,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s;
   double x,y,z,Pr,Pl,Sr,Sl,E1,B1;
   double totalPr,totalPl,totalSr,totalSl,totalE1,totalB1;
   double Wx[3],Wy[3],Wz[3],WxC[3],WyC[3],WzC[3];
   double extE1,extE2,extE3,extB1,extB2,extB3;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;

   k=0; 
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y;  z=p->z;
           WxC[0]=0.5*(1-x)*(1-x);
           WxC[1]=0.75-(0.5-x)*(0.5-x);
           WxC[2]=0.5*x*x;
           WyC[0]=0.5*(1-y)*(1-y);
           WyC[1]=0.75-(0.5-y)*(0.5-y);
           WyC[2]=0.5*y*y;
               
           i1=(int)(i+x+0.5);
           j1=(int)(j+y+0.5);
           x=i+x-i1;
           y=j+y-j1;
           Wx[0]=0.5*(0.5-x)*(0.5-x);
           Wx[1]=0.75-x*x;
           Wx[2]=0.5*(x+0.5)*(x+0.5);
           Wy[0]=0.5*(0.5-y)*(0.5-y);
           Wy[1]=0.75-y*y;
           Wy[2]=0.5*(y+0.5)*(y+0.5);

           Pr=Pl=Sr=Sl=E1=B1=0.0;
           for(jj=0; jj<3; jj++)
             for(ii=0; ii<3; ii++)
             {
               Pr+=D->Pr[i-1+ii][j-1+jj][0]*WxC[ii]*WyC[jj];
               Pl+=D->Pl[i-1+ii][j-1+jj][0]*WxC[ii]*WyC[jj];
               Sr+=D->Sr[i-1+ii][j1-1+jj][0]*WxC[ii]*Wy[jj];
               Sl+=D->Sl[i-1+ii][j1-1+jj][0]*WxC[ii]*Wy[jj];
               E1+=D->Ex[i-1+ii][j1-1+jj][0]*WxC[ii]*Wy[jj];
               B1+=D->Bx[i-1+ii][j-1+jj][0]*WxC[ii]*WyC[jj];
             }

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;
         
           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

}
*/
/*
void interpolation3D_DSX_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,k,ii,jj,kk,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s;
   double x,y,z,Pr,Pl,Sr,Sl,E1,B1;
   double totalPr,totalPl,totalSr,totalSl,totalE1,totalB1;
   double Wx[3],Wy[3],Wz[3],WxC[3],WyC[3],WzC[3];
   double extE1,extE2,extE3,extB1,extB2,extB3;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;
 
   if(D->fieldType==1)
   {   
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
         for(k=kstart; k<kend; k++)
         {
           for(s=0; s<D->nSpecies; s++)
           {
             p=particle[i][j][k].head[s]->pt;
             while(p)
             {
               x=p->x;  y=p->y;  z=p->z;
               WxC[0]=0.5*(1-x)*(1-x);
               WxC[1]=0.75-(0.5-x)*(0.5-x);
               WxC[2]=0.5*x*x;
               WyC[0]=0.5*(1-y)*(1-y);
               WyC[1]=0.75-(0.5-y)*(0.5-y);
               WyC[2]=0.5*y*y;
               WzC[0]=0.5*(1-z)*(1-z);
               WzC[1]=0.75-(0.5-z)*(0.5-z);
               WzC[2]=0.5*z*z;

               
               i1=(int)(i+x+0.5);
               j1=(int)(j+y+0.5);
               k1=(int)(k+z+0.5);
               x=i+x-i1;
               y=j+y-j1;
               z=k+z-k1;
               Wx[0]=0.5*(0.5-x)*(0.5-x);
               Wx[1]=0.75-x*x;
               Wx[2]=0.5*(x+0.5)*(x+0.5);
               Wy[0]=0.5*(0.5-y)*(0.5-y);
               Wy[1]=0.75-y*y;
               Wy[2]=0.5*(y+0.5)*(y+0.5);
               Wz[0]=0.5*(0.5-z)*(0.5-z);
               Wz[1]=0.75-z*z;
               Wz[2]=0.5*(z+0.5)*(z+0.5);

               Pr=Pl=Sr=Sl=E1=B1=0.0;
               for(kk=0; kk<3; kk++)
                 for(jj=0; jj<3; jj++)
                   for(ii=0; ii<3; ii++)
                   {
                     Pr+=D->Pr[i-1+ii][j-1+jj][k1-1+kk]*WxC[ii]*WyC[jj]*Wz[kk];
                     Pl+=D->Pl[i-1+ii][j-1+jj][k1-1+kk]*WxC[ii]*WyC[jj]*Wz[kk];
                     Sr+=D->Sr[i-1+ii][j1-1+jj][k-1+kk]*WxC[ii]*Wy[jj]*WzC[kk];
                     Sl+=D->Sl[i-1+ii][j1-1+jj][k-1+kk]*WxC[ii]*Wy[jj]*WzC[kk];
                     E1+=D->Ex[i-1+ii][j1-1+jj][k1-1+kk]*WxC[ii]*Wy[jj]*Wz[kk];
                     B1+=D->Bx[i-1+ii][j-1+jj][k-1+kk]*WxC[ii]*WyC[jj]*WzC[kk];
                   }


               p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
               p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;
         
               p=p->next;
             }
           }		//for(s)        
         }		   //for(i,j)
   }           //End of fieldType=1

}


void interpolation3D_DSX_1st(Domain *D,External *Ext)
{
   int i,j,k,i1,j1,k1,istart,iend,jstart,jend,kstart,kend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;
   extE2=Ext->E2;
   extE3=Ext->E3;
   extB1=Ext->B1;
   extB2=Ext->B2;
   extB3=Ext->B3;

   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++)
       {
         for(s=0; s<D->nSpecies; s++)
         {
           cnt=0;
           p=particle[i][j][k].head[s]->pt;
           while(p)
           {
             x=p->x;  y=p->y;  z=p->z;
             i1=(int)(i+x+0.5);
             j1=(int)(j+y+0.5);
             k1=(int)(k+z+0.5);
             x1=x+0.5-((int)(x+0.5));
             y1=y+0.5-((int)(y+0.5));
             z1=z+0.5-((int)(z+0.5));

             B1=(1-x1)*(1-y1)*(1-z1)*D->Bx[i1-1][j1-1][k1-1]
               +    x1*(1-y1)*(1-z1)*D->Bx[i1][j1-1][k1-1]
               +(1-x1)*    y1*(1-z1)*D->Bx[i1-1][j1][k1-1]
               +    x1*    y1*(1-z1)*D->Bx[i1][j1][k1-1]
               +(1-x1)*(1-y1)*    z1*D->Bx[i1-1][j1-1][k1]
               +    x1*(1-y1)*    z1*D->Bx[i1][j1-1][k1]
               +(1-x1)*    y1*    z1*D->Bx[i1-1][j1][k1]
               +    x1*    y1*    z1*D->Bx[i1][j1][k1];
             E1=(1-x1)*(1-y)*(1-z)*D->Ex[i1-1][j][k]
               +    x1*(1-y)*(1-z)*D->Ex[i1][j][k]
               +(1-x1)*    y*(1-z)*D->Ex[i1-1][j+1][k]
               +    x1*    y*(1-z)*D->Ex[i1][j+1][k]
               +(1-x1)*(1-y)*    z*D->Ex[i1-1][j][k+1]
               +    x1*(1-y)*    z*D->Ex[i1][j][k+1]
               +(1-x1)*    y*    z*D->Ex[i1-1][j+1][k+1]
               +    x1*    y*    z*D->Ex[i1][j+1][k+1];
             Pr=(1-x1)*(1-y1)*(1-z)*D->Pr[i1-1][j1-1][k]
               +    x1*(1-y1)*(1-z)*D->Pr[i1][j1-1][k]
               +(1-x1)*y1    *(1-z)*D->Pr[i1-1][j1][k]
               +    x1*    y1*(1-z)*D->Pr[i1][j1][k]
               +(1-x1)*(1-y1)*    z*D->Pr[i1-1][j1-1][k+1]
               +    x1*(1-y1)*    z*D->Pr[i1][j1-1][k+1]
               +(1-x1)*    y1*    z*D->Pr[i1-1][j1][k+1]
               +    x1*    y1*    z*D->Pr[i1][j1][k+1];
             Pl=(1-x1)*(1-y1)*(1-z)*D->Pl[i1-1][j1-1][k]
               +    x1*(1-y1)*(1-z)*D->Pl[i1][j1-1][k]
               +(1-x1)*y1    *(1-z)*D->Pl[i1-1][j1][k]
               +    x1*    y1*(1-z)*D->Pl[i1][j1][k]
               +(1-x1)*(1-y1)*    z*D->Pl[i1-1][j1-1][k+1]
               +    x1*(1-y1)*    z*D->Pl[i1][j1-1][k+1]
               +(1-x1)*    y1*    z*D->Pl[i1-1][j1][k+1]
               +    x1*    y1*    z*D->Pl[i1][j1][k+1];
             Sr=(1-x1)*(1-y)*(1-z1)*D->Sr[i1-1][j][k1-1]
               +    x1*(1-y)*(1-z1)*D->Sr[i1][j][k1-1]
               +(1-x1)*    y*(1-z1)*D->Sr[i1-1][j+1][k1-1]
               +    x1*    y*(1-z1)*D->Sr[i1][j+1][k1-1]
               +(1-x1)*(1-y)*    z1*D->Sr[i1-1][j][k1]
               +    x1*(1-y)*    z1*D->Sr[i1][j][k1]
               +(1-x1)*    y*    z1*D->Sr[i1-1][j+1][k1]
               +    x1*    y*    z1*D->Sr[i1][j+1][k1];
             Sl=(1-x1)*(1-y)*(1-z1)*D->Sl[i1-1][j][k1-1]
               +    x1*(1-y)*(1-z1)*D->Sl[i1][j][k1-1]
               +(1-x1)*    y*(1-z1)*D->Sl[i1-1][j+1][k1-1]
               +    x1*    y*(1-z1)*D->Sl[i1][j+1][k1-1]
               +(1-x1)*(1-y)*    z1*D->Sl[i1-1][j][k1]
               +    x1*(1-y)*    z1*D->Sl[i1][j][k1]
               +(1-x1)*    y*    z1*D->Sl[i1-1][j+1][k1]
               +    x1*    y*    z1*D->Sl[i1][j+1][k1];

             p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
             p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)

}

*/
