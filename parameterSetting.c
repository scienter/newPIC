#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>


double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

void parameterSetting(Domain *D,External *Ext, char *input)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *in=NULL;
   int FindParameters();
   int findLoadParameters();
   int findLaserParameters();
   int whatONOFF();
   int whatSaveMode();
   int whatFieldType();
   double minX,maxX,minY,maxY,minZ,maxZ;
   double x,y,z,px,py,pz,gamma;
   double positionX,factor,pMinX,pMaxX,pPosition;
   double normalB,normalE,Ex,Ey,Ez,Bx,By,Bz,dyoverdx,dzoverdx;
   char str[100],name[100],fileName[100];
   int rank,minT,maxT,tmpInt,fail=0,cnt;
   int i,j,k,n,numProbeX,numProbeY,numProbeZ,probeType,id,core,species;
   double lambda,tmpDouble,probeDx,probeDy,probeDz,maxProbeX,minProbeX,maxProbeY,minProbeY,maxProbeZ,minProbeZ,tmp;
   LoadList *LL, *New;
   LaserList *L, *LNew;


   //initially
   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  {
      printf("in [Domain], dimension=?  (1:1D, 2:2D, 3:3D)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"L",input,str)) D->L=atoi(str);
   else  {
      printf("in [Domain], L=?  (Sorry. Fix as L=1)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"M",input,str)) D->M=atoi(str);
   else  {
      printf("in [Domain], M=?  (y directionally dividing number)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"N",input,str)) D->N=atoi(str);
   else  {
      printf("in [Domain], N=?  (z directionally dividing number)\n");
      fail=1;
   }
   if(D->dimension==1)  {
     D->M=1;
     D->N=1;
   }  else if(D->dimension==2) 
     D->N=1;
   if(D->dimension==2)  {
     if(D->M*D->L!=nTasks)  {
       printf("L=%d, M=%d, check the values of L and M.\n",D->L,D->M);
       fail=1;
     }
   }  else if(D->dimension==3)  {
     if(D->M*D->N*D->L!=nTasks)  {
       printf("L=%d, M=%d, N=%d, check the values of L, M, and N.\n",D->L,D->M,D->N);
       fail=1;
     }
   }

   //Field Type
   if(FindParameters("Domain",1,"field_type",input,str)) 
     D->fieldType=whatFieldType(str);
   else  {
      printf("in [Domain], field_type=?  (Split/Yee/Pukhov)\n");
      fail=1;
   }
   //Current Type
   if(FindParameters("Domain",1,"current_order",input,str)) D->currentType=atoi(str);
   else  
      D->currentType=1;
   if(FindParameters("Domain",1,"interpolation_order",input,str)) D->interpolationType=atoi(str);
   else 
      D->interpolationType=1;

   //Boost frame
   if(FindParameters("Domain",1,"boost_gamma",input,str)) D->gamma=atof(str);
   else D->gamma=1;
   if(FindParameters("Domain",1,"boost_ion",input,str)) D->boostIon=whatONOFF(str);
   else D->boostIon=ON;
   if(D->gamma>1)   D->boostOn=ON;
   else             D->boostOn=OFF;
   D->beta=sqrt(1-1.0/D->gamma/D->gamma);

   //save options
   if(FindParameters("Save",1,"field_format",input,str)) 
     D->saveFieldMode=whatSaveMode(str);
   else  
     D->saveFieldMode=TXT;
   if(FindParameters("Save",1,"particle_format",input,str)) 
     D->saveParticleMode=whatSaveMode(str);
   else  
     D->saveParticleMode=TXT;
   if(FindParameters("Save",1,"density_format",input,str)) 
     D->saveDensityMode=whatSaveMode(str);
   else  
     D->saveDensityMode=TXT;
   if(FindParameters("Save",1,"dump_format",input,str)) 
     D->saveDumpMode=whatSaveMode(str);
   else  
     D->saveDumpMode=TXT;
   if(FindParameters("Save",1,"dump_save",input,str)) 
     D->dumpSave=whatONOFF(str);
   else  
     D->dumpSave=OFF;
   if(FindParameters("Save",1,"dump_start",input,str)) 
     D->dumpStart=atoi(str);
   else  
     D->dumpStart=D->saveStart;
   D->dumpStep=0;
   if(FindParameters("Save",1,"field_save",input,str)) 
     D->fieldSave=whatONOFF(str);
   else  
     D->fieldSave=ON;
   if(FindParameters("Save",1,"raman_save",input,str)) 
     D->ramanSave=whatONOFF(str);
   else  
     D->ramanSave=ON;
   if(FindParameters("Save",1,"particle_save",input,str)) 
     D->particleSave=whatONOFF(str);
   else  
     D->particleSave=ON;
   if(FindParameters("Save",1,"density_save",input,str)) 
     D->densitySave=whatONOFF(str);
   else  
     D->densitySave=ON;

   //Domain parameter setting
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], maxStep=? \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_step",input,str)) D->saveStep=atoi(str);
   else  {
      printf("In [Domain], save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  {
      printf("In [Domain], save_start=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"center_save_step",input,str)) D->centerStep=atoi(str);
   else  {
      printf("In [Domain], center_save_step=?\n");
      fail=1;
   }

   if(FindParameters("Save",1,"resolution_change",input,str)) 
      D->resolChange=whatONOFF(str);
   else  D->resolChange=OFF;
   if(FindParameters("Save",1,"resolution_high",input,str)) 
      D->resolHigh=whatONOFF(str);
   else  D->resolHigh=OFF;
   if(FindParameters("Save",1,"resolution_low",input,str)) 
      D->resolLow=whatONOFF(str);
   else  D->resolLow=OFF;
   if(FindParameters("Save",1,"resolution_change_step",input,str)) 
      D->resolStep=atoi(str);
   else  D->resolStep=D->maxTime;
   if(FindParameters("Save",1,"resolution_rate_X",input,str)) 
      D->resolX=atoi(str);
   else  D->resolX=1;
   if(FindParameters("Save",1,"resolution_rate_Y",input,str)) 
      D->resolY=atoi(str);
   else  D->resolY=1;
   if(FindParameters("Save",1,"resolution_rate_Z",input,str)) 
      D->resolZ=atoi(str);
   else  D->resolZ=1;
   if(FindParameters("Domain",1,"minX",input,str)) 
   {
      minX=atof(str);
      minX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("In [Domain], minX=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str)) 
   {
      maxX=atof(str);
      maxX*=D->gamma*(1+D->beta);
   }   else  {
      printf("In [Domain], maxX=? [m].\n");
      fail=1;
   }

   if(D->dimension>1)
   {
     if(FindParameters("Domain",1,"minY",input,str)) 
       minY=atof(str);
     else  {
       printf("In [Domain], minY=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"maxY",input,str)) 
       maxY=atof(str);
     else  {
       printf("In [Domain], maxY=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"dy_over_dx",input,str)) 
       dyoverdx=atof(str);
     else  {
       printf("In [Domain], dy_over_dx=?  [dy/dx]\n");
       fail=1;
     }
   }
   else	;
   if(D->dimension>2)
   {
     if(FindParameters("Domain",1,"minZ",input,str)) 
       minZ=atof(str);
     else  {
       printf("In [Domain], minZ=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"maxZ",input,str)) 
       maxZ=atof(str);
     else  {
       printf("In [Domain], maxZ=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"dz_over_dx",input,str)) 
       dzoverdx=atof(str);
     else  
       dzoverdx=dyoverdx;
   }
   else	;
   if(FindParameters("Domain",1,"moving_domain",input,str)) D->moving=whatONOFF(str);
   else  {
      printf("In [Domain], moving_domain=? [ON/OFF].\n");
      fail=1;
   }   
   if(FindParameters("Domain",1,"lambda",input,str)) 
   {
      D->lambda=atof(str);
      D->lambda*=D->gamma*(1+D->beta);
   }
   else  {
      printf("In [Domain], lambda=? [m]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"division_lambda",input,str)) 
   {
      D->divisionLambda=atof(str);
   }
   else  {
      printf("In [Domain], divisionLambda=? [number of devided wavelength]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dt_ratio",input,str)) 
      D->dtRatio=atof(str);
   else  {
      printf("In [Domain], dt_ratio=? [<1.0]\n");
      fail=1;
   }

   //pml
   if(FindParameters("PML",1,"pml",input,str)) 
     D->pmlOn=whatONOFF(str);
   else  
     D->pmlOn=OFF;
   if(FindParameters("PML",1,"right_pml_cells",input,str))
     D->pmlCellRight=atoi(str);
   else  D->pmlCellRight=1;
   if(FindParameters("PML",1,"left_pml_cells",input,str))
     D->pmlCellLeft=atoi(str);
   else  D->pmlCellLeft=1;
   if(FindParameters("PML",1,"up_pml_cells",input,str))
     D->pmlCellUp=atoi(str);
   else  D->pmlCellUp=1;
   if(FindParameters("PML",1,"down_pml_cells",input,str))
     D->pmlCellBottom=atoi(str);
   else  D->pmlCellBottom=1;
   if(FindParameters("PML",1,"pml_r",input,str))
     D->pmlr=atof(str);
   else  {
     printf("In [PML], pml_r=? (retarding length)\n");
     fail=1;
   }
   if(FindParameters("PML",1,"pml_d",input,str))
     D->pmld=atof(str);
   else  {
     printf("In [PML], pml_d=? (damping length)\n");
     fail=1;
   }


   //External field parameter setting
   if(FindParameters("External",1,"Ex",input,str)) Ex=atof(str);
   else  {
      printf("In [External], Ex=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ey",input,str)) Ey=atof(str);
   else  {
      printf("In [External], Ey=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ez",input,str)) Ez=atof(str);
   else  {
      printf("In [External], Ez=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bx",input,str)) Bx=atof(str);
   else  {
      printf("In [External], Bx=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"By",input,str)) By=atof(str);
   else  {
      printf("In [External], By=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bz",input,str)) Bz=atof(str);
   else  {
      printf("In [External], Bz=? [Tesla]\n");
      fail=1;
   }


   //additional Domain parameters  
   D->dx=1.0/D->divisionLambda;
   D->nx=((int)((maxX-minX)/D->lambda/D->dx));
   D->ny=1;
   D->nz=1;
   tmpDouble=D->dx/(D->gamma*(1+D->beta));
   D->dy=D->dz=1.0;
   D->dt=D->dx*D->dtRatio;
   if(D->dimension>1)
   {
     D->dy=D->dx*dyoverdx;
     D->dy=tmpDouble*dyoverdx;
     if(D->dy<=D->dx*1.5)   {
       printf("dy_over_dx is too low. It must be over than %g.\n", D->gamma*(1+D->beta)*1.5);
       fail=1;
     }
     else	;
     D->ny=((int)((maxY-minY)/D->lambda/D->dy));

     if(D->fieldType==Yee)
       D->dt=D->dtRatio/sqrt(1.0/D->dx/D->dx+1.0/D->dy/D->dy);
     else	;
   }
   else	;
   if(D->dimension>2)
   {
     D->dz=D->dx*dzoverdx;
     D->dz=tmpDouble*dzoverdx;
     if(D->dz<=D->dx*1.5)   {
       printf("dz_over_dx is too low. It must be over than %g.\n", D->gamma*(1+D->beta)*1.5);
       fail=1;
     }
     else	;
     D->nz=((int)((maxZ-minZ)/D->lambda/D->dz));

     if(D->fieldType==Yee)
       D->dt=D->dtRatio/sqrt(1.0/D->dx/D->dx+1.0/D->dy/D->dy+1.0/D->dz/D->dz);
     else	;
   }
   else	;  

   D->minXDomain=D->minYDomain=D->minZDomain=0;
   D->omega=2*pi*velocityC/D->lambda;
   if(D->boostOn==ON)   {
     D->minXSub=-D->nx;
     if(D->dimension>1)
       D->minYDomain=(int)(minY/D->lambda/D->dy);
     else	;
     if(D->dimension>2)
       D->minZDomain=(int)(minZ/D->lambda/D->dz);
     else	;
   }
   else   {
     D->minXSub=0;
     if(D->dimension>1)
       D->minYDomain=(int)(minY/D->lambda/D->dy);
     if(D->dimension>2)
       D->minZDomain=(int)(minZ/D->lambda/D->dz);
   }

   if(myrank==0)
   {
     printf("dx=%g,dt=%g\n",D->dx,D->dt);
     if(D->dimension==2)
     {
       printf("gamma=%g, dx=%g, dy=%g, dt=%g, divisionLambda=%g\n",D->gamma,D->dx,D->dy,D->dt,D->divisionLambda);
     }
     else if(D->dimension==3)
     {
       printf("gamma=%g, dx=%g, dy=%g, dz=%g, dt=%g, divisionLambda=%g\n",D->gamma,D->dx,D->dy,D->dz,D->dt,D->divisionLambda);
     }
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);

   //ID track
   if(FindParameters("Domain",1,"tracking",input,str)) 
     D->tracking=whatONOFF(str);
   else  
     D->tracking=OFF;
   if(D->tracking==ON)
   {
     if(FindParameters("Domain",1,"track_save_step",input,str)) 
       D->trackSaveStep=atoi(str);
     else  
       D->trackSaveStep=1;
     if(myrank==0)   
     {
       in=fopen("trackFile","r");
       if(in!=NULL)  
       {
         if(D->dimension==2)
         {
           cnt=0;
           while(fscanf(in,"%lf %lf %lf %lf %lf %d %d %d"
                   ,&x,&y,&px,&py,&pz,&id,&core,&species)!=EOF)
             cnt++;
           D->idNums=cnt;
         }
         else if(D->dimension==3)   {
           cnt=0;
           while(fscanf(in,"%lf %lf %lf %lf %lf %lf %d %d %d"
                   ,&x,&y,&z,&px,&py,&pz,&id,&core,&species)!=EOF)
             cnt++;
           D->idNums=cnt;
         }
         fclose(in);
       }
       else   {
         printf("'trackFile' is missing.\n");
         fail=1;
       }
     }
     else	;
     MPI_Bcast(&(D->idNums),1,MPI_INT,0,MPI_COMM_WORLD);
   }
   else D->idNums=0;

   if(D->idNums>0)
   {
     D->trackID=(int *)malloc(D->idNums*sizeof(int));
     D->trackCore=(int *)malloc(D->idNums*sizeof(int));
     D->trackS=(int *)malloc(D->idNums*sizeof(int));
    
     if(myrank==0)
     {
       in=fopen("trackFile","r");
       if(D->dimension==2)
       {
         for(i=0; i<D->idNums; i++)
           fscanf(in,"%lf %lf %lf %lf %lf %d %d %d"
                     ,&x,&y,&px,&py,&pz
                     ,&(D->trackID[i]),&(D->trackCore[i]),&(D->trackS[i]));
       }
       else if(D->dimension==3)
       {
         for(i=0; i<D->idNums; i++)
           fscanf(in,"%lf %lf %lf %lf %lf %lf %d %d %d"
                     ,&x,&y,&z,&px,&py,&pz
                     ,&(D->trackID[i]),&(D->trackCore[i]),&(D->trackS[i]));
       }
       fclose(in);
     }	
     else	;
     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Bcast(D->trackID,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->trackCore,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->trackS,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
   }     
/*
   //Probe parameter
   if(FindParameters("Probe",1,"probeType",input,str)) probeType=atoi(str);
   else probeType=0;
   D->probeNum=0;

   if(probeType==0)
   {
     if(FindParameters("Probe",1,"probeNum",input,str)) D->probeNum=atoi(str);
     else  {
       printf("in [Probe], probeNum=?  [ea]\n");
       fail=1;
     }
     if(D->probeNum>0)
     {
       D->probeX=(int *)malloc(D->probeNum*sizeof(int));
       D->probeY=(int *)malloc(D->probeNum*sizeof(int));
       D->probeZ=(int *)malloc(D->probeNum*sizeof(int));
       for(i=0; i<D->probeNum; i++)
       {
         sprintf(name,"probeX%d",i);
         if(FindParameters("Probe",1,name,input,str))   
           D->probeX[i]=((int)(atof(str)/D->lambda/D->dx));      
         else  {
           printf("in [Probe], probeX%d=?\n",i);
           fail=1;
         }
         sprintf(name,"probeY%d",i);
         if(FindParameters("Probe",1,name,input,str))      
           D->probeY[i]=((int)(atof(str)/D->lambda/D->dy));      
         else  {
           printf("in [Probe], probeY%d=?\n",i);
           fail=1;
         }
         sprintf(name,"probeZ%d",i);
         if(FindParameters("Probe",1,name,input,str))      
           D->probeZ[i]=((int)(atof(str)/D->lambda/D->dz));      
         else  {
           printf("in [Probe], probeZ%d=?\n",i);
           fail=1;
         }
       }
     }  
   }
   else if(probeType==1)
   {
     if(FindParameters("Probe",1,"minProbeX",input,str)) minProbeX=atof(str);
     else  {
       printf("in [Probe], minProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeX",input,str)) maxProbeX=atof(str);
     else  {
       printf("in [Probe], maxProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeX",input,str)) numProbeX=atoi(str);
     else 
       numProbeX=1;
     if(FindParameters("Probe",1,"minProbeY",input,str)) minProbeY=atof(str);
     else  {
       printf("in [Probe], minProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeY",input,str)) maxProbeY=atof(str);
     else  {
       printf("in [Probe], maxProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeY",input,str)) numProbeY=atoi(str);
     else 
       numProbeY=1;
     if(FindParameters("Probe",1,"minProbeZ",input,str)) minProbeZ=atof(str);
     else  {
       printf("in [Probe], minProbeZ=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeZ",input,str)) maxProbeZ=atof(str);
     else  {
       printf("in [Probe], maxProbeZ=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeZ",input,str)) numProbeZ=atoi(str);
     else 
       numProbeZ=1;
     if(numProbeX==0 || numProbeY==0 || numProbeZ==0)  {
       printf("in [Probe], it must be that numProbeX or numProbeY or numProbeZ > 0 !!\n");
       fail=1;
     }
           
     probeDx=(maxProbeX-minProbeX)/((double)numProbeX);
     probeDy=(maxProbeY-minProbeY)/((double)numProbeY);
     probeDz=(maxProbeZ-minProbeZ)/((double)numProbeZ);
     D->probeNum=numProbeX*numProbeY*numProbeZ;
     D->probeX=(int *)malloc(D->probeNum*sizeof(int));
     D->probeY=(int *)malloc(D->probeNum*sizeof(int));
     D->probeZ=(int *)malloc(D->probeNum*sizeof(int));

     n=0;
     for(i=0; i<numProbeX; i++)
       for(j=0; j<numProbeY; j++)
         for(k=0; k<numProbeZ; k++)
         {       
           tmpDouble=minProbeX+i*probeDx;
           D->probeX[n]=((int)(tmpDouble/D->lambda/D->dx));      
           tmpDouble=minProbeY+j*probeDy;
           D->probeY[n]=((int)(tmpDouble/D->lambda/D->dy));     
           tmpDouble=minProbeZ+k*probeDz;
           D->probeZ[n]=((int)(tmpDouble/D->lambda/D->dz));     
           n++;
         } 
   }			//End of else if(probeType=2)
*/    
   //additional Boost parameters
   factor=D->gamma*(1+D->beta);
   D->minT=(int)(D->maxStep/factor/factor); 	//boost frame iteration
   D->maxT=(int)((D->maxStep+D->beta*D->nx)/(1+D->beta)-factor*D->gamma*D->minT*D->beta);	//boost frame iteration
   if(myrank==0)
     printf("maxT=%d, nx=%d, ny=%d, nz=%d\n",D->maxT,D->nx,D->ny,D->nz);
   else	;

   //additional external field parameters
   normalB=eMass*D->omega/(-eCharge);
   normalE=normalB*velocityC;
   Ext->E1=Ex/normalE;
   Ext->E2=Ey/normalE;
   Ext->E3=Ez/normalE;
   Ext->B1=Bx/normalB;
   Ext->B2=By/normalB;
   Ext->B3=Bz/normalB;

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input)) 
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

   //Plasma parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
      exit(0);
   else	;

}

int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   double positionX,positionY,positionZ;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"polarity",input,str)) polarity=atoi(str);
   else  polarity=0;

   if(polarity)
   {
     if(FindParameters("Laser",rank,"mode",input,str)) 
        L->mode=atoi(str);
     else  L->mode=0;
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
     {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     }
     else  L->lambda=D->lambda;
  
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"retard",input,str)) L->retard=atof(str);
     else  {
        printf("in [Laser], retard=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  {
        printf("in [Laser], loadPositionX=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionY",input,str)) positionY=atof(str);
     else  positionY=0;
     if(FindParameters("Laser",rank,"loadPositionZ",input,str)) positionZ=atof(str);
     else  positionZ=0;
     if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
     else  {
        printf("in [Laser], beamWaist=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  {
        printf("in [Laser], focus=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  L->flat=0.0;
     if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
     else  L->direction=1;

     //additional laser parameters
     L->polarity=polarity;
     L->omega=2*pi*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dx));   
     L->loadPointY=((int)(positionY/D->lambda/D->dy));   
     L->loadPointZ=((int)(positionZ/D->lambda/D->dz));   
     L->rayleighLength=pi/(L->lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;
     if(fail==1)
        exit(0);
   }
   return polarity;
}

int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   LoadList *New;
   int whatSpecies();
   int whatPlasmaType();
   double whatMass(int species);
   double pointPosition,wp,pDt;
   int whatCharge();
   int whatFunctionMode();
   int whatDefineMode();
   char name[100], str[100];
   int i,n,cnt,species,fail=0;
   double tmp,max,min;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Plasma",rank,"type",input,name)) 
   {
     LL->type = whatPlasmaType(name);
     if(D->boostOn==ON)
       LL->type = BoostFrame;
   }
   else LL->type=0;

   if(LL->type>0)
   {
      if(FindParameters("Plasma",rank,"density",input,str)) 
      {
         LL->density=atof(str);
         LL->density*=D->gamma;
      }
      else  {
         printf("in [Plasma], density=? [m-3]\n");
         fail=1;
      }

/*    
      //testing optimal dx(divisionLambda) size
      wp=sqrt(LL->density*eCharge*eCharge/eMass/eps0);
      pDt=2*pi/wp;   
      pDt=pDt/(2*pi/D->omega)/20.0;
      if(D->dt>pDt)
      {
         printf("dt must be less then %g!\n",pDt);
         printf("So, divisionLambda>%g!\n",1/pDt);
         fail=1;
      }
*/

      if(FindParameters("Plasma",rank,"species",input,name)) 
         species = whatSpecies(name);
      else  species = 0;
      LL->species=species;
      if(FindParameters("Plasma",rank,"numberInCell",input,str)) 
         LL->numberInCell=atoi(str);
      else  {
         printf("in [Plasma], numberInCell=? \n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"startIndex",input,str)) 
         LL->index=atoi(str);
      else  
         LL->index=0;
      if(FindParameters("Plasma",rank,"temperature",input,str))  
         LL->temperature=atof(str);
      else   LL->temperature=0.0;	
      if(FindParameters("Plasma",rank,"given_min_px",input,str)) 
         LL->givenMinPx = atof(str);
      else  LL->givenMinPx = -1e9;
      LL->mass=whatMass(species);
      LL->charge=whatCharge(species);
      LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
//      LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy*D->lambda*D->dz/LL->numberInCell;

      srand(1*(myrank+1));
      switch (LL->type)  {
      case Polygon :
        if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
        else  {
          printf("in [Plasma], Xnodes=?\n");
          printf("Each nodes indicates the point of plasma density changing.\n");
          fail=1;
        }
        if(LL->xnodes>0)
        {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)
          {
            sprintf(name,"X%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xpoint[i] = atof(str)/D->gamma/D->lambda/D->dx;
            else 
            { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xn[i] = atof(str);
            else 
            { printf("Xn%d should be defined.\n",i);  fail=1; } 
          }
        }
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"Ynodes",input,str)) LL->ynodes=atoi(str);
          else  {
            printf("in [Plasma], Ynodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          if(LL->ynodes>0)
          {
            LL->ypoint = (double *)malloc(LL->ynodes*sizeof(double));
            LL->yn = (double *)malloc(LL->ynodes*sizeof(double));   
            for(i=0; i<LL->ynodes; i++)
            {
              sprintf(name,"Y%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) {
                LL->ypoint[i] = atof(str)/D->lambda/D->dy;
              }
              else 
              { printf("Y%d should be defined.\n",i);  fail=1; }
 
              sprintf(name,"Yn%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) 
                LL->yn[i] = atof(str);
              else 
              { printf("Yn%d should be defined.\n",i);  fail=1; } 
            }
          }
        }
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"Znodes",input,str)) LL->znodes=atoi(str);
          else  {
            printf("in [Plasma], Znodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          LL->zpoint = (double *)malloc(LL->znodes*sizeof(double));
          LL->zn = (double *)malloc(LL->znodes*sizeof(double));   
          for(i=0; i<LL->znodes; i++)
          {
            sprintf(name,"Z%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) {
              LL->zpoint[i] = atof(str)/D->lambda/D->dz;
            }
            else 
            { printf("Z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Zn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->zn[i] = atof(str);
            else 
            { printf("Zn%d should be defined.\n",i);  fail=1; } 
          }
        }
        if(FindParameters("Plasma",rank,"centerX",input,str))  
          LL->centerX=atof(str)/D->lambda/D->dx;
        else   LL->centerX=0.0;	
        LL->centerY=0.0;	
        LL->centerZ=0.0;	
        if(FindParameters("Plasma",rank,"gauss_coef_X",input,str))  
          LL->gaussCoefX=atof(str)/D->lambda/D->dx;
        else   LL->gaussCoefX=1.0;
        if(FindParameters("Plasma",rank,"poly_coef_X",input,str))  
          LL->polyCoefX=atof(str)/D->lambda/D->dx;
        else   LL->polyCoefX=0.0;	
        if(FindParameters("Plasma",rank,"function_mode_X",input,str))  
          LL->modeX=whatFunctionMode(str);
        else   LL->modeX=0;	
        if(FindParameters("Plasma",rank,"function_mode_YZ",input,str))  
          LL->modeYZ=whatFunctionMode(str);
        else   LL->modeYZ=0;	
        if(D->dimension>1)
        {	
          if(FindParameters("Plasma",rank,"centerY",input,str))  
            LL->centerY=atof(str)/D->lambda/D->dy;
          else   LL->centerY=0.0;	
          if(FindParameters("Plasma",rank,"gauss_coef_YZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/D->dy;
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"poly_coef_YZ",input,str))  
            LL->polyCoefYZ=atof(str)*D->lambda*D->dy*D->lambda*D->dy;
          else   LL->polyCoefYZ=0.0;	
        }
        else if(D->dimension>2)
        {	
          if(FindParameters("Plasma",rank,"centerZ",input,str))  
            LL->centerZ=atof(str)/D->lambda/D->dz;
          else   LL->centerZ=0.0;	
          if(FindParameters("Plasma",rank,"gauss_coef_YZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/sqrt(D->dy*D->dy+D->dz*D->dz);
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"poly_coef_YZ",input,str))  
            LL->polyCoefYZ=atof(str)*D->lambda*D->lambda*(D->dy*D->dy+D->dz*D->dz);
          else   LL->polyCoefYZ=0.0;	
        }
        else 	;
        break;

      case Defined :
//        srand(time(NULL)*(myrank+1));
        if(FindParameters("Plasma",rank,"define_mode",input,str))  
          LL->defineMode=whatDefineMode(str);
        else   LL->defineMode=byNumber;	
        if(FindParameters("Plasma",rank,"number_defined",input,str))  
          LL->numDefined=atoi(str);
        else   LL->numDefined=0;	
        if(LL->defineMode==byDensity)
        {
          if(FindParameters("Plasma",rank,"minX",input,str))  
            LL->minX=atof(str);
          else   {   printf("minX=?  [m]\n");  exit(0);   }	
          if(FindParameters("Plasma",rank,"maxX",input,str))  
            LL->maxX=atof(str);
          else   {   printf("maxX=?  [m]\n");  exit(0);   }	
          if(D->dimension>1)
          {
            if(FindParameters("Plasma",rank,"minY",input,str))  
              LL->minY=atof(str);
            else   {   printf("minY=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxY",input,str))  
              LL->maxY=atof(str);
            else   {   printf("maxY=?  [m]\n");  exit(0);   }	
          }
          if(D->dimension>2)
          {
            if(FindParameters("Plasma",rank,"minZ",input,str))  
              LL->minZ=atof(str);
            else   {   printf("minZ=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxZ",input,str))  
              LL->maxZ=atof(str);
            else   {   printf("maxZ=?  [m]\n");  exit(0);   }	
          }
        }
        else 	;

        if(FindParameters("Plasma",rank,"xlength_particle",input,str))  
          LL->xLengDef=atof(str)/D->lambda;
        else   LL->xLengDef=0.0;
        if(LL->numDefined>0)	
        {
          LL->xPosition=(double *)malloc(LL->numDefined*sizeof(double));
          shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
          for(i=0; i<LL->numDefined; i++)
          {
            if(LL->defineMode==byNumber)
            {
              sprintf(name,"xPosition%d",i);
              if(FindParameters("Plasma",rank,name,input,str))  
                LL->xPosition[i]=atof(str)/D->lambda;
              else
              { printf("xPosition%d should be defined.\n",i); fail=1;}
            }
            else if(LL->defineMode==byDensity)
            {
              if(myrank==0)
                shareDouble[i]=(LL->minX+randomV()*(LL->maxX-LL->minX))/D->lambda;
              else 	;
            }
            else 	;
          }
          if(LL->defineMode==byDensity)
          {
            MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for(i=0; i<LL->numDefined; i++)
              LL->xPosition[i]=shareDouble[i];
          }
          else	;
          free(shareDouble);
        }
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"ylength_particle",input,str))  
            LL->yLengDef=atof(str)/D->lambda;
          else   LL->yLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->yPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"yPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->yPosition[i]=atof(str)/D->lambda;
                else
                { printf("yPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minY+randomV()*(LL->maxY-LL->minY))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->yPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);
          }
        }		//End of demension>1
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"zlength_particle",input,str))  
            LL->zLengDef=atof(str)/D->lambda;
          else   LL->zLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->zPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"zPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->zPosition[i]=atof(str)/D->lambda;
                else
                { printf("zPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minZ+randomV()*(LL->maxZ-LL->minZ))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->zPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);  
          }
        }
        if(D->dimension==2)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef/D->dx/D->dy*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((cnt*2)*sizeof(double ));
 
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
              }
          }
        }  
        else if(D->dimension==3)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef*LL->zLengDef/D->dx/D->dy/D->dz*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((LL->numDefPtcls*3)*sizeof(double ));
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
                tmp=(double)(randomV());
                LL->define[i][n+2*cnt]=LL->zPosition[i]-0.5*LL->zLengDef+tmp*LL->zLengDef;
              }
          }
/*
        for(i=0; i<LL->numDefined; i++)
          for(n=0; n<cnt; n++)
            printf("n=%d, x=%g,y=%g,z=%g\n",n,LL->define[i][n],LL->define[i][n+cnt],LL->define[i][n+2*cnt]);          
*/
        } 	//End fo dimension=3 : defined Plasma 
        LL->maxLoadTime=(int)((max+1)*D->divisionLambda);
        LL->minLoadTime=(int)((min-1)*D->divisionLambda);
        break;
      }
   
   }	//end of if(species)

   if(fail==1)
      exit(0);

   return LL->type;
}

int whatDefineMode(char *str)
{
   if(strstr(str,"by_number")) 		return byNumber;
   else if(strstr(str,"by_density"))   	return byDensity;
   else 				return byNumber;
}

int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else 				return OFF;
}

int whatSaveMode(char *str)
{
   if(strstr(str,"TXT")) 		return TXT;
   else if(strstr(str,"HDF"))   	return HDF;
   else 				return TXT;
}

int whatFieldType(char *str)
{
   if(strstr(str,"Split")) 		return Split;
   else if(strstr(str,"Yee"))   	return Yee;
   else if(strstr(str,"Pukhov"))   	return Pukhov;
   else 				return 0;
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else if(strstr(str,"AlPlus11"))   	return AlPlus11;
   else return 0;
}

int whatPlasmaType(char *str)
{
   if(strstr(str,"Polygon"))         return Polygon;
   else if(strstr(str,"Defined"))   	return Defined;
   else
   {
     printf("No Plasma type!\n"); 
     exit(0);
   }
   return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)           return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)           return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)           return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)           return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)           return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)           return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)           return (12.0111-6*eMassU)/eMassU;
   else if(species == AlPlus11)         return (26.9815-11*eMassU)/eMassU;
   else {  printf("Species' mass not defined\n");  exit(0);  }
}

int whatCharge(int species)
{
   int fail;

   if(species == Electron) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0)           return 0;
   else if(species == CPlus1)           return 1;
   else if(species == CPlus2)           return 2;
   else if(species == CPlus3)           return 3;
   else if(species == CPlus4)           return 4;
   else if(species == CPlus5)           return 5;
   else if(species == CPlus6)           return 6;
   else if(species == AlPlus11)           return 11;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

int whatFunctionMode(char *str)
{
   if(strstr(str,"Constant")) 		return Constant;
   else if(strstr(str,"Gaussian"))   	return Gaussian;
   else if(strstr(str,"Polynomial"))   	return Polynomial;
   else return 0;
}

