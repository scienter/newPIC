#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreDen(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s,char *dataName,int step,int *offset)
{
   char name[100],fileName[100];
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   sprintf(fileName,"denParticle%d.h5",step);
   sprintf(name,"%d/numberInCell",s);
   if(myrank==0)  {
     restoreIntMeta(fileName,name,&D->numberInCell);
   }    else   ;
   MPI_Bcast(&D->numberInCell,1,MPI_INT,0,MPI_COMM_WORLD);

   switch (D->dimension) {
   case 2 :
     sprintf(name,"%d/%s1",s,dataName);
     restoreFieldComp(den1,fileName,name,D->nx,D->ny,1,D->nxSub,D->nySub,1,D->istart,D->iend,D->jstart,D->jend,0,1,offset);
     sprintf(name,"%d/%s2",s,dataName);
     restoreFieldComp(den2,fileName,name,D->nx,D->ny,1,D->nxSub,D->nySub,1,D->istart,D->iend,D->jstart,D->jend,0,1,offset);
     sprintf(name,"%d/%s3",s,dataName);
     restoreFieldComp(den3,fileName,name,D->nx,D->ny,1,D->nxSub,D->nySub,1,D->istart,D->iend,D->jstart,D->jend,0,1,offset);
     sprintf(name,"%d/%s4",s,dataName);
     restoreFieldComp(den4,fileName,name,D->nx,D->ny,1,D->nxSub,D->nySub,1,D->istart,D->iend,D->jstart,D->jend,0,1,offset);
     break;
   }
}
