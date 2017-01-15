#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void loadLaser(Domain *D,LaserList *L,double t)
{
  void loadLaser1D_Split(Domain *D,LaserList *L,double t);
  void loadLaser2D_Split(Domain *D,LaserList *L,double t);
  void loadLaser3D_DSX(Domain *D,LaserList *L,double t);
  void loadLaser2D_Yee_Pukhov(Domain *D,LaserList *L,double t);

  if(D->boostOn==OFF)
  {
    switch((D->fieldType-1)*3+D->dimension)  {
    case (Split-1)*3+1 :
      loadLaser1D_Split(D,L,t);
      break;
    case (Split-1)*3+2 :
      loadLaser2D_Split(D,L,t);
      break;
//    case (Split-1)*3+3 :
//      loadLaser3D_DSX(D,L,t);
//      break;
    case (Yee-1)*3+2 :
      loadLaser2D_Yee_Pukhov(D,L,t);
      break;
    case (Pukhov-1)*3+2 :
      loadLaser2D_Yee_Pukhov(D,L,t);
      break;
    default :
      printf("In loadLaser, what is field_type? and what is dimension?\n");
    }
  }
}

void loadLaser2D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp;
   double x,y,z,r2,w2;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=L->loadPointY;	//y position

   t0=2*rU+L->retard;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*pi*L->omega/D->omega;
   kx=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;

   if(laserOK==1)
   {
     if(L->polarity==2)
     {
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ey[positionX][j][k]=amp;            
       }
     }  
/*         
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
*/
   }     //End of field is OK
}

void loadLaser1D_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat,omega,amp;
   int istart,iend,positionX,j,k,laserOK=0;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   t0=2.0*rU;
   omega=2.0*pi*L->omega/D->omega;

   if(t<2.0*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2.0*rU && t<2.0*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2.0*rU+flat && t<2.0*rU+flat+2.0*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2.0*rU+2.0*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=1;
   else ;
   positionX=L->loadPointX+istart-D->minXSub;

   if(laserOK==1)   {
     j=k=0;
     if(L->polarity==2)     {
       amp=longitudinal*sin(omega*t);
       D->Pr[positionX][j][k]=amp;            
       D->Pl[positionX][j][k]=amp;           
     } else if(L->polarity==3)   {
       amp=longitudinal*sin(omega*t);
       D->Sr[positionX][j][k]=amp;            
       D->Sl[positionX][j][k]=amp;           
     } else ;
   }  else ;   //End of field is OK
}

void loadLaser2D_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp;
   double x,y,z,r2,w2,retard;
   int istart,iend,jstart,jend,kstart,kend,minj,maxj;
   int positionX,rank,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;  k=0;

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;
   retard=L->retard;

   jC=L->loadPointY;	//y position

   t0=2*rU+retard;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*pi*L->omega/D->omega;
   kx=2*pi*D->lambda/L->lambda;

   if(t<2*rU+retard)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU+retard && t<2*rU+flat+retard) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat+retard && t<2*rU+flat+2*rD+retard) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat+retard) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;
   if(laserOK==1)
   {
     switch (L->mode)  {
     case DEFAULT :
       w2=w*w;
       if(L->polarity==2)       {
         for(j=0; j<jend+3; j++)         {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Pr[positionX][j][k]=amp;            
           D->Pl[positionX][j][k]=amp;           
         }
       } 
       else if(L->polarity==3)       {
         for(j=jstart; j<jend; j++)         {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Sr[positionX][j][k]=amp;            
           D->Sl[positionX][j][k]=amp;           
         }
       }
       break;

     case 1 : //adding
       if(L->polarity==2)    {
         for(j=jstart; j<jend; j++)    {
           y=(j-jstart+D->minYSub-jC)*D->dy;    r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Pr[positionX][j][k]+=amp;            
           D->Pl[positionX][j][k]+=amp;           
         }
       }
       else if(L->polarity==3)       {
         for(j=jstart; j<jend; j++)    {
           y=(j-jstart+D->minYSub-jC)*D->dy;   r2=y*y;
           pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
           amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
           D->Sr[positionX][j][k]+=amp;            
           D->Sl[positionX][j][k]+=amp;           
         }
       }
       break;
     default :
       printf("In loadLaser.c, what laser mode?\n");
       break;
     }		//End of switch
   }  else ;   	//End of field is OK
}

/*
void loadLaser3D_DSX(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,w,phi,omega,kx,pphi,amp,Phi;
   double x,y,z,r,r2,w2;
   int istart,iend,jstart,jend,kstart,kend;
   int positionX,rank,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=L->loadPointY;	//y position
   kC=L->loadPointZ;	//z position
//   jC=(int)(D->ny*0.5);	//center position
//   kC=(int)(D->nz*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   omega=2*pi*L->omega/D->omega;
   kx=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;


   if(laserOK==1)
   {
     if(L->polarity==2)
     {
       w2=w*w;
       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         for(k=kstart; k<kend; k++)
         {
           z=(k-kstart+D->minZSub-kC)*D->dz;
           r2=y*y+z*z;
           pphi=x/zR*r2/w2-phi+kx*x-omega*t;
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           D->Pr[positionX][j][k]=amp;            
           D->Pl[positionX][j][k]=amp;           
         } 
       }
     }  
     else if(L->polarity==12)
     {
       for(k=kstart; k<kend; k++)
       {
         z=(k-kstart+D->minZSub-kC)*D->dz;
         if(z>=0) 
         {
           for(j=jstart; j<jend; j++)
           {
             y=(j-jstart+D->minYSub-jC)*D->dy;
//lala
             r=sqrt(y*y+z*z);
             if(r==0)  Phi=0.0;
             else      Phi=acos(y/r);
           
             pphi=Phi-2.0*phi+w0*w0/w/w/zR/zR*0.5*kx*r*r*x+kx*x-omega*t;
             amp=longitudinal*w0/w*sqrt(2.0)*r/w*exp(-r*r/w/w)*cos(pphi);
             D->Pr[positionX][j][k]=amp;            
             D->Pl[positionX][j][k]=amp;           
           }
         } 
         else  
         {
           for(j=jstart; j<jend; j++)
           {
             y=(j-jstart+D->minYSub-jC)*D->dy;
//lala
             r=sqrt(y*y+z*z);
             if(r==0)  Phi=0.0;
             else      Phi=-acos(y/r);           
             pphi=Phi-2.0*phi+w0*w0/w/w/zR/zR*0.5*kx*r*r*x+kx*x-omega*t;
             amp=longitudinal*w0/w*sqrt(2.0)*r/w*exp(-r*r/w/w)*cos(pphi);
             D->Pr[positionX][j][k]=amp;            
             D->Pl[positionX][j][k]=amp;           
           }
         } 
       }
     }  
     else if(L->polarity==13)
     {
       for(k=kstart; k<kend; k++)
       {
         z=(k-kstart+D->minZSub-kC)*D->dz;
         if(z>=0) 
         {
           for(j=jstart; j<jend; j++)
           {
             y=(j-jstart+D->minYSub-jC)*D->dy;
//lala
             r=sqrt(y*y+z*z);
             if(r==0)  Phi=0.0;
             else      Phi=acos(y/r);
           
             pphi=Phi-2.0*phi+w0*w0/w/w/zR/zR*0.5*kx*r*r*x+kx*x-omega*t;
             amp=longitudinal*w0/w*sqrt(2.0)*r/w*exp(-r*r/w/w)*cos(pphi);
             D->Sr[positionX][j][k]=amp;            
             D->Sl[positionX][j][k]=amp;           
           }
         } 
         else  
         {
           for(j=jstart; j<jend; j++)
           {
             y=(j-jstart+D->minYSub-jC)*D->dy;
//lala
             r=sqrt(y*y+z*z);
             if(r==0)  Phi=0.0;
             else      Phi=-acos(y/r);           
             pphi=Phi-2.0*phi+w0*w0/w/w/zR/zR*0.5*kx*r*r*x+kx*x-omega*t;
             amp=longitudinal*w0/w*sqrt(2.0)*r/w*exp(-r*r/w/w)*cos(pphi);
             D->Sr[positionX][j][k]=amp;            
             D->Sl[positionX][j][k]=amp;           
           }
         } 
       }
     }  

   }     //End of laserOK
}
*/

/*
void boostLoadLaser2D(Domain *D,LaserList *L)
{
   double x,x0,rU,rD,beta,gamma,labWaist,f,yy;
   double zR,w0,z,w,phi,omega,k,y,pphi,amp,unitCompen;
   int i,j,jC;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;

   jC=(int)(D->ny*0.5);	//center position
   gamma=D->gamma;
   beta=D->beta;

   labWaist=L->beamWaist*D->lambda;
   zR=L->rayleighLength/gamma;
//   zR=pi/(L->lambda/D->gamma/(1.0+D->beta))*labWaist*labWaist/gamma/D->lambda;
   w0=labWaist/D->lambda;  
   f=-L->focus/gamma;   
   omega=k=2*pi;
   x0=-2*rU;
   unitCompen=gamma*gamma*(1.0-beta*beta);

   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2)
     {
       for(i=2; i<D->nxSub+2; i++)
         for(j=2; j<D->nySub+2; j++)
         {
           x=(i+D->minXSub-2)*D->dx;
           y=(j-2+D->minYSub-jC)*D->dy;
           z=(f+x+x0);          
           w=w0*sqrt(1.0+z*z/zR/zR);
           if(x>=2*x0)
           {
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
//             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
             field[i][j].Pr=amp;            
//             field[i][j].Pl=amp;            
           }
           else      {
             field[i][j].Pr=0.0;            
//             field[i][j].Pl=0.0;            
           }
         }
     }
     else if(L->polarity==3)
     {
       for(i=1; i<=D->nxSub; i++)
         for(j=1; j<=D->nySub; j++)
         {
           x=(i+D->minXSub-1)*D->dx;
           z=f+x+x0;          
           if(x>=2*x0)
           {
             y=(j-1+D->minYSub-jC)*D->dy;
             w=w0*sqrt(1.0+z*z/zR/zR);
             phi=atan(z/zR);
             pphi=z/zR*y*y/w/w-0.5*phi+unitCompen*k*z;
             amp=L->amplitude*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi)*exp(-(x-x0)*(x-x0)/rU/rU);
             field[i][j].Sr=amp;            
             field[i][j].Sl=amp;            
           }
           else      {
             field[i][j].Sr=0.0;            
             field[i][j].Sl=0.0;            
           }
         }
     }
   }     //End of fieldType=1
}
*/
/*
void loadLaser2D(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat;
   double zR,w0,z,w,phi,omega,k,y,pphi,amp;
   int istart,iend,jstart,jend;
   int positionX,positionY,rank,j,jC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=(int)(D->ny*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   z=-L->focus;   
   w=w0*sqrt(1.0+z*z/zR/zR);
   phi=atan(z/zR);
   omega=2*pi*L->omega/D->omega;
   k=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   positionY=L->loadPointY+jstart-D->minYSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
      laserOK=1;


   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2 && laserOK==1)
     {
//         field[positionX][positionY].Pr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Pl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Pr=amp;            
         field[positionX][j].Pl=amp;
       }

     }           
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=z/zR*y*y/w/w-0.5*phi+k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
   }     //End of fieldType=1

}
*/
/*
void loadLaserOpp2D(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat,loadPosition;
   double zR,w0,z,w,phi,omega,k,y,pphi,amp;
   int istart,iend,jstart,jend;
   int positionX,positionY,rank,j,jC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;


   rU=L->rU*D->divisionLambda*D->dt;
   rD=L->rD*D->divisionLambda*D->dt;
   flat=L->flat*D->divisionLambda*D->dt*L->lambda/D->lambda;

   jC=(int)(D->ny*0.5);	//center position

   t0=2*rU;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   loadPosition=L->loadPointX*D->dx;
   z=loadPosition-L->focus;   
   w=w0*sqrt(1.0+z*z/zR/zR);
   phi=atan(z/zR);
   omega=2*pi*L->omega/D->omega;
   k=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub-2;
   positionY=L->loadPointY+jstart-D->minYSub;
   if(positionX>D->minXSub && positionX<=D->maxXSub)
      laserOK=1;


   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     if(L->polarity==2 && laserOK==1)
     {
//         field[positionX][positionY].Pr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Pl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=-z/zR*y*y/w/w-0.5*phi-k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Pr=amp;            
         field[positionX][j].Pl=amp;
       }

     }           
     else if(L->polarity==3 && laserOK==1)
     {
//         field[positionX][positionY].Sr=longitudinal*sin(omega*t);  

//         field[positionX][positionY].Sl=longitudinal*sin(omega*t);

       for(j=jstart; j<jend; j++)
       {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         pphi=-z/zR*y*y/w/w-0.5*phi-k*z-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-y*y/w/w)*sin(pphi);
         field[positionX][j].Sr=amp;            
         field[positionX][j].Sl=amp;
       }

     }
   }     //End of fieldType=1

}
*/
