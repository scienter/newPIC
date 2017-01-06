#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

double randomValue(double beta);
double maxwellianVelocity(double temperature);
void loadMovingPlasma_crystal(Domain *D,LoadList *LL,int s);
double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX);
double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ);

void loadMovingPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
  int rankX;
  int myrank, nTasks;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  void loadMovingPolygonPlasma2D();
  void loadDefinedPlasma2D();
  void loadMovingPolygonPlasma3D();
  void loadDefinedPlasma3D();
  void loadBoostPlasma2D();

  switch((LL->type-1)*3+D->dimension)
  {
  case 1:
//    loadMovingPlasma_crystal(D,LL,s);
    break;
  case (Polygon-1)*3+2:
    rankX=myrank/(D->M*D->N);
    if(rankX==D->L-1)
      loadMovingPolygonPlasma2D(D,LL,s,iteration); 
    else	;
    break;
  case (Defined-1)*3+2:
    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime)
      loadDefinedPlasma2D(D,LL,s); 
    break;
  case (BoostFrame-1)*3+2:
    loadBoostPlasma2D(D,LL,s,iteration); 
    break;

  //3D
  case (Polygon-1)*3+3:
    loadMovingPolygonPlasma3D(D,LL,s,iteration); 
    break;
  case (Defined-1)*3+3:
    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime)
      loadDefinedPlasma3D(D,LL,s); 
    break;

  default:
    ;
  }
}


void loadBoostPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,l,intNum,cnt;
   int posX,posY,iter,t;
   double tmp,dx,dy;
   double wp,pDt,v1,v2,v3,gamma,beta[2];

   double ne,randTest=0,positionX,positionY;


   ptclList *New,*p;   
   Particle ***particle;
   particle=D->particle;

   dx=D->dx;
   dy=D->dy;

   beta[0]=D->beta;
   beta[1]=1.0;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   srand(iteration);

   //position define      
   iter=0;
   k=0;
   for(i=iend-1; i<=iend; i++)
   {
     for(j=jstart; j<jend; j++)
     {
       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {           
           posX=i+D->minXSub-istart;
           posY=j+D->minYSub-jstart;
           if((posX>=LL->xpoint[l] && posX<LL->xpoint[l+1]) && 
              (posY>=LL->ypoint[t] && posY<LL->ypoint[t+1]))
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])
                *(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])
                *(posY-LL->ypoint[t])+LL->yn[t]);
             ne*=LL->numberInCell*beta[iter];	//it is the double number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;

             while(cnt<intNum)
             {               
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);
  
               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;

               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k+0;

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

               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue(1.0))
             {
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k+0;

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
             } 		//end of randTest
           }		//end of if(x,y nodes)
         }		//end of for(xnodes,ynodes)  

     } 				//End of for(j)
     iter=iter+1;
   }				//End of for(i)
}

void loadDefinedPlasma2D(Domain *D,LoadList *LL,int s)
{
   int ii,i,j,k,istart,iend,jstart,jend,kstart,kend,cnt;
   int l,n,myrank;
   double dx,dy,v1,v2,v3;
   double ne,randTest,positionX,positionY,x,y,xPos,yPos;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   dx=D->dx;
   dy=D->dy;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   cnt=LL->numDefPtcls;
   k=0;

   //position define      
   i=iend-1;
   for(n=0; n<LL->numDefined; n++)
   {
     for(l=0; l<cnt; l++)
     {
       xPos=LL->define[n][l];
       ii=(int)(xPos/dx-D->minXSub+istart);
       if(i==ii)
       {
         yPos=LL->define[n][l+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         if(jstart<=j && j<jend)
         {
           x=xPos/dx-D->minXSub+istart-i;
           y=yPos/dy-D->minYSub+jstart-j;
           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;
  
           New->x = x;
           New->oldX=i+x;
           New->y = y;
           New->oldY=j+y;
           New->z = 0;
           New->oldZ=k +0;
  
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
         }
       }		//if(ii)
     }			//End of for(l)
   }			//End of for(n)
}

void loadDefinedPlasma3D(Domain *D,LoadList *LL,int s)
{
   int ii,i,j,k,istart,iend,jstart,jend,kstart,kend,cnt;
   int l,n,myrank;
   double dx,dy,dz,v1,v2,v3;
   double ne,randTest,positionX,positionY,x,y,z,xPos,yPos,zPos;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   dx=D->dx;
   dy=D->dy;
   dz=D->dz;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
   kstart=D->kstart;
   kend=D->kend;
   cnt=LL->numDefPtcls;

   //position define      
   i=iend-1;
   for(n=0; n<LL->numDefined; n++)
   {
     for(l=0; l<cnt; l++)
     {
       xPos=LL->define[n][l];
       ii=(int)(xPos/dx-D->minXSub+istart);
       if(i==ii)
       {
         yPos=LL->define[n][l+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         zPos=LL->define[n][l+2*cnt];
         k=(int)(zPos/dz-D->minZSub+kstart);
         if(jstart<=j && j<jend && kstart<=k && k<kend)
         {
           x=xPos/dx-D->minXSub+istart-i;
           y=yPos/dy-D->minYSub+jstart-j;
           z=zPos/dz-D->minZSub+kstart-k;
           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;
  
           New->x = x;
           New->oldX=i+x;
           New->y = y;
           New->oldY=j+y;
           New->z = z;
           New->oldZ=k +z;
  
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
         }
       }		//if(ii)
     }			//End of for(l)
   }			//End of for(n)
}

void loadMovingPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,t,myrank;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,tmp;
   double ne,randTest,positionX,positionY,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;
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

   srand(iteration*(myrank+1));

   //position define      
   i=iend-1;
     for(j=jstart; j<jend; j++)
     {
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
             
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue(1.0);
               positionY=randomValue(1.0);

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
               New->p1=-D->gamma*D->beta+v1;
               New->p2=v2;
               New->p3=v3;
               New->weight=1.0/LL->numberInCell;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank; 
             }		//end of if(randTest)
           }	
         } 		//end of for(lnodes)  
     }			//End of for(i,j)
         
}

void loadMovingPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,m,n,myrank;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ;
   double ne,randTest,positionX,positionY,positionZ,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
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

   srand(iteration*(myrank+1));

   //position define      
   i=iend-1;
   for(j=jstart; j<jend; j++)
     for(k=kstart; k<kend; k++)
     {
       for(l=0; l<LL->xnodes-1; l++)
         for(m=0; m<LL->ynodes-1; m++)
           for(n=0; n<LL->znodes-1; n++)
           {
             posX=i+D->minXSub-istart;
             posY=j+D->minYSub-jstart;
             posZ=k+D->minZSub-kstart;
 
             if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
                posY>=LL->ypoint[m] && posY<LL->ypoint[m+1] &&
                posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])
             {
               ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
               ne*=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
               ne*=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
               applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
               applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
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
               }		//end of if(randTest)
             }			//End of if(pos...)
           } 			//end of for(l,m,n)  
     }				//End of for(i,j,k)
         
}

