#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    int i,j,k,n,s,iteration=0,filter,boost,filterStep,labSaveStep;
    int rnk,suddenDump=OFF;
    double factor,time_spent;
    clock_t begin,end;
    double t;
    char name[100];
    FILE *out;
    Domain D;  
    LaserList *L;
    LoadList *LL;
    External Ext;
    UPML UPml;
    DPML DPml;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }

    //parameter setting
    parameterSetting(&D,&Ext,argv[1]);
    if(argc >= 3)      D.dumpStep = atoi(argv[2]); else;

    //create mesh
    boundary(&D,&Ext);
    MPI_Barrier(MPI_COMM_WORLD);

    //load plasma or load dump file
    if(argc >= 3)  {   
      iteration=D.dumpStep;
      restoreDump(D,iteration);
      t=D.dt*iteration;
      sprintf(name,"dumpField%d.h5",iteration);
      restoreIntMeta(name,"/minXDomain",&(D.minXDomain),1);
      D.minXSub+=D.minXDomain;
    }  else   {
      LL=D.loadList;
      s=0;
      while(LL->next)      {
        loadPlasma(&D,LL,s,iteration);
        LL=LL->next;
        s++;
      }
      t=0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    //rooping time 
    while(iteration<=D.maxStep+D.resolX)
    {

//      if(D.tracking==ON && iteration%D.trackSaveStep==0)
//      {
//        if(D.dimension==2)
//          trackID(&D,iteration,D.istart,D.iend,D.jstart,D.jend,0,1);
//        else if(D.dimension==3)
//          trackID(&D,iteration,D.istart,D.iend,D.jstart,D.jend,D.kstart,D.kend);
//      }
    
       //calculating running time     
//       end=clock();
//       time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
//       if(myrank==0)
//       {
//         if(time_spent>D.maxTime)
//           suddenDump=ON;
//         else	;       
//       } 
//       else	;
//       if(myrank==0)
//       {
//         for(rnk=1; rnk<nTasks; rnk++)
//           MPI_Send(&suddenDump,1,MPI_INT,rnk,myrank,MPI_COMM_WORLD);
//       }
//       else

 //      MPI_Barrier(MPI_COMM_WORLD);
//       if(suddenDump==ON)
//       {
//         if(D.saveMode==TXT) 	saveDump(D,iteration);
//         else if(D.saveMode==HDF) 	saveDumpHDF(D,iteration);
//         iteration=D.maxStep+1;
//       }
//       else	;

       //resolChange
       if(D.resolChange==ON && D.resolHigh==ON)  {
         if(iteration==D.resolStep)   {
             saveBDump(D,iteration);
             saveEDump(D,iteration);
             saveJDump(D,iteration);
             saveDenParticleHDF(&D,iteration);
         } else	;
         if(iteration+1==D.resolStep)   {
             saveBDump(D,iteration);
             saveJDump(D,iteration);
             saveDenParticleHDF(&D,iteration);
         } else	;
         if(iteration-1==D.resolStep)   {
             saveBDump(D,iteration);
             saveJDump(D,iteration);
             saveDenParticleHDF(&D,iteration);
         } else	;
       }  else if(D.resolChange==ON && D.resolLow==ON)  {
         if(iteration==D.resolStep)   {
             saveEDump(D,iteration);
         } else	;
         if(iteration-D.resolX/2==D.resolStep)   {
             saveBDump(D,iteration);
             saveJDump(D,iteration);
         } else	;
         if(iteration-D.resolX/2+1==D.resolStep)   {
             saveBDump(D,iteration);
             saveJDump(D,iteration);
//             saveDumpParticleResolHDF(&D,iteration);
         } else	;
       }       

       //save File      
       if(iteration%D.saveStep==0 && iteration>=D.saveStart)   {
         saveFile(D,iteration);
         if(D.dumpSave==ON && iteration>=D.dumpStart) 
           saveDump(D,iteration);  else	;
       }  else	;
       //save center field      
       if(iteration%D.centerStep==0)  {
         saveCenterField(&D,iteration);
//         saveCenterDensity(&D,iteration);
       }  else	;
       MPI_Barrier(MPI_COMM_WORLD);

       fieldSolve(D,t);

       interpolation(&D,&Ext);

       particlePush(&D);

       updateCurrent(D);

       if(iteration>=D.nx && D.moving==ON && D.boostOn==OFF 
          && (iteration-D.nx)%D.shiftDuration!=0 )    {
         movingDomain(&D,iteration);
         LL=D.loadList;
         s=0;
         while(LL->next)      {
           loadMovingPlasma(&D,LL,s,iteration);
           LL=LL->next;
           s++;
         }
         rearrangeParticles(&D);
         if(D.L>1)   particleShareX(D);   else	;
         if(D.M>1)   particleShareY(D);   else	;
         removeEdge(D);
       } else       {
          rearrangeParticles(&D);
          if(D.L>1)  particleShareX(D);   else	;
          if(D.M>1)  particleShareY(D);   else	;
          removeEdge(D);
       }

       //time update
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;
       t=D.dt*iteration;  

    }     //end of time roop                  

//    if(D.tracking==ON)
//      saveTracking(&D);

    end=clock();
    time_spent=(end-begin)/CLOCKS_PER_SEC;

    //make 'report' file
    if(myrank==0)
    {
      sprintf(name,"report");
      out = fopen(name,"w");
      fprintf(out,"nx=%d\n",D.nx);
      fprintf(out,"ny=%d\n",D.ny);
      fprintf(out,"nz=%d\n",D.nz);
      fprintf(out,"cores=%d\n",nTasks);
      fprintf(out,"running time=%gm\n",time_spent/60.0);
      fclose(out);
    }
    else	;

    cleanMemory(&D);
    
    MPI_Finalize();

    return 0;
}
