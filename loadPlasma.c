#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


double maxwellianVelocity(double temperature);
double randomValue(double beta);
void loadPlasma_crystal(Domain *D,LoadList *LL,int s);
void random1D_sobol(double *x,gsl_qrng *q);
void random2D_sobol(double *x,double *y,gsl_qrng *q);
void random3D_sobol(double *x,double *y,double *z);


double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-(x-centerX)*(x-centerX)/gaussCoefX/gaussCoefX);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefX*(x-centerX)*(x-centerX);
    break;
  }
  return result;
}

double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ))/gaussCoefYZ/gaussCoefYZ);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefYZ*((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
    break;
  }
  return result;
}

void loadPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
  void loadPolygonPlasma1D();
  void loadPolygonPlasma2D();
//  void loadPolygonPlasma3D();
  void loadChannelPlasma();
//  void loadPlasma_crystal(Domain *D,LoadList *LL,int s);


  switch((LL->type-1)*3+D->dimension)  {
  //1D
  case ((Polygon-1)*3+1):
    loadPolygonPlasma1D(D,LL,s,iteration); 
//    loadPlasma_crystal(D,LL,s);
    break;

  //2D
  case ((Polygon-1)*3+2):
    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;
//  case ((Polygon-1)*3+3):
//    loadPolygonPlasma3D(D,LL,s,iteration); 
    break;

  default:
    ;
  }
}

void loadPolygonPlasma1D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,intNum,cnt,l,t,modeX;
   double posX,v1,v2,v3,centerX,tmp,weight;
   double ne,randTest,positionX,gaussCoefX,polyCoefX;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   centerX=LL->centerX;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   modeX=LL->modeX;
   weight=1.0/LL->numberInCell;

   srand(myrank-1);

   //position define      
   j=k=0;
   for(i=istart; i<iend; i++)
     {
       gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,1);
       for(l=0; l<LL->xnodes-1; l++)
         {
           posX=(double)(i+D->minXSub-istart);
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
             ne*=tmp;
             ne*=LL->numberInCell;
             intNum=(int)ne;
             randTest=ne-intNum;
             if(randTest>randomValue(1.0)) intNum=intNum+1; else;

             cnt=0;
             while(cnt<intNum)
             {      
//               positionX=randomValue(1.0);
               random1D_sobol(&positionX,q);
               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = 0.0;
               New->oldY=j+0.0;
               New->z = 0.0;
               New->oldZ=k+0.0;

               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=v1;
               New->p2=v2;
               New->p3=v3;
               New->weight=weight;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            

               cnt++; 
             }		//end of while(cnt)

           } else;	
         } 		//end of for(lnodes)  

       gsl_qrng_free(q);
     }			//End of for(i,j)         
}


void loadPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,t;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,tmp,weight;
   double ne,randTest,positionX,positionY,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   k=0;
   centerX=LL->centerX;
   centerY=LL->centerY;
   centerZ=LL->centerZ;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   gaussCoefYZ=LL->gaussCoefYZ;
   polyCoefYZ=LL->polyCoefYZ;
   modeX=LL->modeX;
   modeYZ=LL->modeYZ;
   weight=1.0/LL->numberInCell;

   srand(myrank+1);

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,2);

       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {
           posX=(double)(i+D->minXSub-istart);
           posY=(double)(j+D->minYSub-jstart);
           posZ=(double)(k+D->minZSub-kstart);
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
              posY>=LL->ypoint[t] && posY<LL->ypoint[t+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])*(posY-LL->ypoint[t])+LL->yn[t]);
             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
             ne*=tmp;
             tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
             ne*=tmp;
             ne*=LL->numberInCell;	//it is the double number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             if(randTest>randomValue(1.0)) intNum=intNum+1;
 
             cnt=0;
             while(cnt<intNum)
             {      
//               positionX=randomValue(1.0);
//               positionY=randomValue(1.0);
               random2D_sobol(&positionX,&positionY,q);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k +0;

               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=v1;
               New->p2=v2;
               New->p3=v3;
               New->weight=weight;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            

               cnt++; 
             }		//end of while(cnt)

           }	
         } 		//end of for(lnodes)  
       gsl_qrng_free(q);

     }			//End of for(i,j)
         
}

void loadPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,m,n;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ;
   double ne,randTest,positionX,positionY,positionZ,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   centerX=LL->centerX;
   centerY=LL->centerY;
   centerZ=LL->centerZ;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   gaussCoefYZ=LL->gaussCoefYZ;
   polyCoefYZ=LL->polyCoefYZ;
   modeX=LL->modeX;
   modeYZ=LL->modeYZ;

   srand(myrank+1);

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++)
       {
         for(l=0; l<LL->xnodes-1; l++)
           for(m=0; m<LL->ynodes-1; m++)
             for(n=0; n<LL->znodes-1; n++)
             {
               posX=(double)(i+D->minXSub-istart);
               posY=(double)(j+D->minYSub-jstart);
               posZ=(double)(k+D->minZSub-kstart);
 
               if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
                  posY>=LL->ypoint[m] && posY<LL->ypoint[m+1] &&
                  posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])
               {
                 ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
                 ne*=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
                 ne*=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
                 ne*=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
                 ne*=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
                 ne*=LL->numberInCell;	//it is the double number of superparticles.
                 intNum=(int)ne;
                 randTest=ne-intNum;
             
                 cnt=0;
                 while(cnt<intNum)
                 {               
                   positionX=randomValue(1.0);
                   positionY=randomValue(1.0);
                   positionZ=randomValue(1.0);
  
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j][k].head[s]->pt;
                   particle[i][j][k].head[s]->pt = New;
 
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   New->z = positionZ;
                   New->oldZ=k +positionZ;
  
                   New->E1=New->E2=New->E3=0.0;
                   New->B1=New->B2=New->B3=0.0;
                   v1=maxwellianVelocity(LL->temperature)/velocityC;
                   v2=maxwellianVelocity(LL->temperature)/velocityC;
                   v3=maxwellianVelocity(LL->temperature)/velocityC;
                   New->p1=-D->gamma*D->beta+v1;
                   New->p2=v2;
                   New->p3=v3;
                   New->weight=1.0/LL->numberInCell;
                   LL->index+=1;
                   New->index=LL->index;            
                   New->core=myrank;            
   
                   cnt++; 
                 }		//end of while(cnt)

                 if(randTest>randomValue(1.0))
                 {
                   positionX=randomValue(1.0);
                   positionY=randomValue(1.0);
                   positionZ=randomValue(1.0);

                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j][k].head[s]->pt;
                   particle[i][j][k].head[s]->pt = New;
   
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   New->z = positionZ;
                   New->oldZ=k + positionZ;
                   New->E1=New->E2=New->E3=0.0;
                   New->B1=New->B2=New->B3=0.0;
                   v1=maxwellianVelocity(LL->temperature)/velocityC;
                   v2=maxwellianVelocity(LL->temperature)/velocityC;
                   v3=maxwellianVelocity(LL->temperature)/velocityC;
                   New->p1=-D->gamma*D->beta+v1;
                   New->p2=v2;
                   New->p3=v3;
                   New->weight=1.0/LL->numberInCell;
                   LL->index+=1;
                   New->index=LL->index;            
                   New->core=myrank;            
                 }		//end of if(randTest)
               }		//End of if (l,m,n)
             }			//End of for(l,m,n)

       }		//End of for(i,j,k)         
}

double maxwellianVelocity(double temperature)
{
   double vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=100, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}

void random1D_sobol(double *x,gsl_qrng *q)
{
   double v[1];

   gsl_qrng_get(q,v);
   *x=v[0];
}

void random2D_sobol(double *x,double *y,gsl_qrng *q)
{
   double v[2];

   gsl_qrng_get(q,v);
   *x=v[0];
   *y=v[1];
}

void random3D_sobol(double *x,double *y,double *z)
{
   double v[3];
   gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,3);

   gsl_qrng_get(q,v);
   *x=v[0];
   *y=v[1];
   *z=v[2];

   gsl_qrng_free(q);
}

