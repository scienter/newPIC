#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void MPI_TransferRho_Xplus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Xminus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Yminus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_TransferRho_Yplus(Domain *D,double ***f1,int nx,int nz,int share);
void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void solveDensity1D(Domain *D,int s,double coef);
void solveDensity2D(Domain *D,int s,double coef);
void solveDensity3D(Domain *D,int s,double coef);
void saveCoordHDF(Domain *D,char *fileName);
void density_xdmf(int dimension,char *fileName,int nx,int ny,int nz,int s);
void setZero(double ***den,int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);
void save1D_P_Grid_HDF(Domain *D,int iteration);

void saveP_GridHDF(Domain D,int iteration)
{
   switch (D.dimension) {
   case 1 :
     save1D_P_Grid_HDF(&D,iteration);     
     break;
   default :
     printf("In saveP_GridHDF, what dimension?\n");
   }
}

void save1D_P_Grid_HDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int nxSub,nySub,nzSub,nxSub1D,nySub2D,nzSub3D;
    int rankX,rankY,rankZ,biasX,biasY,biasZ;
    double x,y,z,px,py,pz,xwl,xwr,ywl,ywr,weight;
    double ***den1,***den2;
    int offset[3];
    char dataName[100],fileName[100];
    LoadList *LL;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;
    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    nxSub=D->nxSub;  nySub=D->nySub;  nzSub=D->nzSub;  
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend; j=0; 
    kstart=D->kstart;    kend=D->kend; k=0;
    nxSub1D=D->nxSub+5;    nySub2D=1;    nzSub3D=1;

    LL=D->loadList;
    s=0; while(LL->next) { LL=LL->next; s++; }
    int numberInCell[s];
    LL=D->loadList;
    s=0; while(LL->next) { numberInCell[s]=LL->numberInCell; LL=LL->next; s++;}

    //save denParticle
    sprintf(fileName,"PGrid%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }    else        ;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank==0)  {
        saveIntMeta(fileName,"/nx",&D->nx,1);
        saveIntMeta(fileName,"/ny",&D->ny,1);
        saveIntMeta(fileName,"/nz",&D->nz,1);
        saveIntMeta(fileName,"/minXDomain",&D->minXDomain,1);
        saveIntMeta(fileName,"/minYDomain",&D->minYDomain,1);
        saveIntMeta(fileName,"/minZDomain",&D->minZDomain,1);
        saveIntMeta(fileName,"/nSpecies",&D->nSpecies,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    den1=memoryAsign(nxSub1D,nySub2D,nzSub3D);
    den2=memoryAsign(nxSub1D,nySub2D,nzSub3D);
    nx=D->nx+5; ny=1; nz=1;
    calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);

    offset[0]=D->minXSub-D->minXDomain+biasX;
    offset[1]=0; offset[2]=0;

    for(s=0; s<D->nSpecies; s++)
    {
      if(myrank==0)    {
        file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
        sprintf(dataName,"%d",s);
        group_id=H5Gcreate2(file_id,dataName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Gclose(group_id);
        H5Fclose(file_id);
      }  else         ;
      MPI_Barrier(MPI_COMM_WORLD);

      setZero(den1,nxSub1D,nySub2D,nzSub3D);
      setZero(den2,nxSub1D,nySub2D,nzSub3D);
      for(i=D->istart; i<D->iend; i++)    {
        p=particle[i][j][k].head[s]->pt;
        while(p)            {
          px=p->weight*p->p1;
          xwl=1.0-p->x;  xwr=1.0-xwl;
          den1[i][j][k]+=xwl*px; den2[i][j][k]+=xwr*px;
          p=p->next;
        }
      }
      sprintf(dataName,"%d/px1",s);
      saveFieldComp(den1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      sprintf(dataName,"%d/px2",s);
      saveFieldComp(den2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      setZero(den1,nxSub1D,nySub2D,nzSub3D);
      setZero(den2,nxSub1D,nySub2D,nzSub3D);
      for(i=D->istart; i<D->iend; i++)     {
        p=particle[i][j][k].head[s]->pt;
        while(p)            {
          px=p->weight*p->p2;
          xwl=1.0-p->x;  xwr=1.0-xwl;
          den1[i][j][k]+=xwl*px;
          den2[i][j][k]+=xwr*px;
          p=p->next;
        }
      }
      sprintf(dataName,"%d/py1",s);
      saveFieldComp(den1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      sprintf(dataName,"%d/py2",s);
      saveFieldComp(den2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      setZero(den1,nxSub1D,nySub2D,nzSub3D);
      setZero(den2,nxSub1D,nySub2D,nzSub3D);
      for(i=D->istart; i<D->iend; i++) {
        p=particle[i][j][k].head[s]->pt;
        while(p)            {
          px=p->weight*p->p3;
          xwl=1.0-p->x;  xwr=1.0-xwl;
          den1[i][j][k]+=xwl*px;
          den2[i][j][k]+=xwr*px;
          p=p->next;
        }
      }
      sprintf(dataName,"%d/pz1",s);
      saveFieldComp(den1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      sprintf(dataName,"%d/pz2",s);
      saveFieldComp(den2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

      //save density
      setZero(den1,nxSub1D,nySub2D,nzSub3D);
      setZero(den2,nxSub1D,nySub2D,nzSub3D);
      for(i=D->istart; i<D->iend; i++)    {
        p=particle[i][j][k].head[s]->pt;
        while(p)            {
          weight=p->weight;
          xwl=1.0-p->x;  xwr=1.0-xwl;
          den1[i][j][k]+=xwl*weight;
          den2[i][j][k]+=xwr*weight;
          p=p->next;
        }
      }
      sprintf(dataName,"%d/q1",s);
      saveFieldComp(den1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
      sprintf(dataName,"%d/q2",s);
      saveFieldComp(den2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
    }

    deleteField(den1,nxSub1D,nySub2D,nzSub3D);
    deleteField(den2,nxSub1D,nySub2D,nzSub3D);

    if(myrank==0)  printf("%s\n",fileName);  else;
}

/*
void save2D_P_Grid_HDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int nxSub,nySub,nzSub,nxSub1D,nySub2D,nzSub3D;
    double x,y,z,px,py,pz,xwl,xwr,ywl,ywr,weight;
    double ***den1,***den2,***den3,***den4;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    LoadList *LL;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;      nx=D->nx;
    istart=D->istart;    iend=D->iend;

    nxSub1D=D->nxSub+5;
    nySub2D=1;
    nzSub3D=1;
    if(D->dimension>1) nySub2D=D->nySub+5; else	;
    if(D->dimension>2) nzSub3D=D->nzSub+5; else	;

    LL=D->loadList;
    s=0;
    while(LL->next)      {
      LL=LL->next;
      s++;
    }
    int numberInCell[s];
    LL=D->loadList;
    s=0;
    while(LL->next)      {
      numberInCell[s]=LL->numberInCell;
      LL=LL->next;
      s++;
    }

    switch(D->dimension) {
    //2D
    case 2:
      //save denParticle
      sprintf(name,"denParticle%d.h5",iteration);
      if(myrank==0)     {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
      }    else        ;
      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank==0)  {
        saveIntMeta(name,"/nx",&nx,1);
        saveIntMeta(name,"/ny",&ny,1);
        saveIntMeta(name,"/nz",&nz,1);
      } else      ;
      MPI_Barrier(MPI_COMM_WORLD);

      den1=memoryAsign(nxSub1D,nySub2D,nzSub3D);
      den2=memoryAsign(nxSub1D,nySub2D,nzSub3D);
      den3=memoryAsign(nxSub1D,nySub2D,nzSub3D);
      den4=memoryAsign(nxSub1D,nySub2D,nzSub3D);
      offset[0]=D->minXSub-D->minXDomain;
      offset[1]=D->minYSub-D->minYDomain;
      offset[2]=0;

      for(s=0; s<D->nSpecies; s++)
      {
        if(myrank==0)    {
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(dataName,"%d",s);
          group_id=H5Gcreate2(file_id,dataName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }  else         ;
        MPI_Barrier(MPI_COMM_WORLD);

        setZero(den1,nxSub1D,nySub2D,nzSub3D);
        setZero(den2,nxSub1D,nySub2D,nzSub3D);
        setZero(den3,nxSub1D,nySub2D,nzSub3D);
        setZero(den4,nxSub1D,nySub2D,nzSub3D);
        k=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)            {
              px=p->weight*p->p1;
              xwl=1.0-p->x;  xwr=1.0-xwl;
              ywl=1.0-p->y;  ywr=1.0-ywl;
              den1[i][j][k]+=xwl*ywl*px;
              den2[i][j][k]+=xwr*ywl*px;
              den3[i][j][k]+=xwl*ywr*px;
              den4[i][j][k]+=xwr*ywr*px;
              p=p->next;
            }
          }
        sprintf(dataName,"%d/px1",s);
        saveFieldComp(den1,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/px2",s);
        saveFieldComp(den2,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/px3",s);
        saveFieldComp(den3,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/px4",s);
        saveFieldComp(den4,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

        setZero(den1,nxSub1D,nySub2D,nzSub3D);
        setZero(den2,nxSub1D,nySub2D,nzSub3D);
        setZero(den3,nxSub1D,nySub2D,nzSub3D);
        setZero(den4,nxSub1D,nySub2D,nzSub3D);
        k=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)            {
              px=p->weight*p->p2;
              xwl=1.0-p->x;  xwr=1.0-xwl;
              ywl=1.0-p->y;  ywr=1.0-ywl;
              den1[i][j][k]+=xwl*ywl*px;
              den2[i][j][k]+=xwr*ywl*px;
              den3[i][j][k]+=xwl*ywr*px;
              den4[i][j][k]+=xwr*ywr*px;
              p=p->next;
            }
          }
        sprintf(dataName,"%d/py1",s);
        saveFieldComp(den1,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/py2",s);
        saveFieldComp(den2,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/py3",s);
        saveFieldComp(den3,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/py4",s);
        saveFieldComp(den4,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

        setZero(den1,nxSub1D,nySub2D,nzSub3D);
        setZero(den2,nxSub1D,nySub2D,nzSub3D);
        setZero(den3,nxSub1D,nySub2D,nzSub3D);
        setZero(den4,nxSub1D,nySub2D,nzSub3D);
        k=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)            {
              px=p->weight*p->p3;
              xwl=1.0-p->x;  xwr=1.0-xwl;
              ywl=1.0-p->y;  ywr=1.0-ywl;
              den1[i][j][k]+=xwl*ywl*px;
              den2[i][j][k]+=xwr*ywl*px;
              den3[i][j][k]+=xwl*ywr*px;
              den4[i][j][k]+=xwr*ywr*px;
              p=p->next;
            }
          }
        sprintf(dataName,"%d/pz1",s);
        saveFieldComp(den1,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/pz2",s);
        saveFieldComp(den2,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/pz3",s);
        saveFieldComp(den3,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/pz4",s);
        saveFieldComp(den4,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

        if(myrank==0)  {
          sprintf(dataName,"%d/numberInCell",s);
          saveIntMeta(name,dataName,&numberInCell[s],1);
        } else      ;
        MPI_Barrier(MPI_COMM_WORLD);
      }
      if(myrank==0)  printf("%s\n",name);  else;

      //save density
      sprintf(name,"den%d.h5",iteration);
      if(myrank==0)     {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
      }    else        ;
      if(myrank==0)  {
        saveIntMeta(name,"/nx",&nx,1);
        saveIntMeta(name,"/ny",&ny,1);
        saveIntMeta(name,"/nz",&nz,1);
      } else      ;
      MPI_Barrier(MPI_COMM_WORLD);

      for(s=0; s<D->nSpecies; s++)
      {
        if(myrank==0)    {
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(dataName,"%d",s);
          group_id=H5Gcreate2(file_id,dataName,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }  else         ;
        MPI_Barrier(MPI_COMM_WORLD);

        setZero(den1,nxSub1D,nySub2D,nzSub3D);
        setZero(den2,nxSub1D,nySub2D,nzSub3D);
        setZero(den3,nxSub1D,nySub2D,nzSub3D);
        setZero(den4,nxSub1D,nySub2D,nzSub3D);
        k=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)            {
              weight=p->weight;
              xwl=1.0-p->x;  xwr=1.0-xwl;
              ywl=1.0-p->y;  ywr=1.0-ywl;
              den1[i][j][k]+=xwl*ywl*weight;
              den2[i][j][k]+=xwr*ywl*weight;
              den3[i][j][k]+=xwl*ywr*weight;
              den4[i][j][k]+=xwr*ywr*weight;
              p=p->next;
            }
          }
        sprintf(dataName,"%d/den1",s);
        saveFieldComp(den1,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/den2",s);
        saveFieldComp(den2,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/den3",s);
        saveFieldComp(den3,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%d/den4",s);
        saveFieldComp(den4,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      }

      deleteField(den1,nxSub1D,nySub2D,nzSub3D);
      deleteField(den2,nxSub1D,nySub2D,nzSub3D);
      deleteField(den3,nxSub1D,nySub2D,nzSub3D);
      deleteField(den4,nxSub1D,nySub2D,nzSub3D);

      if(myrank==0)  printf("%s\n",name);  else;
      break;
    }   //End of switch (dimension)
}
*/

//lala
void saveDensityHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,offSetY,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    int offset[3];
    double *rho0;
    double charge;
    char name[100],dataName[100],fileName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=nxSub+2;
    jstart=D->jstart;
    jend=nySub+2;
    kstart=D->kstart;
    kend=nzSub+2;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    s=0;
    rho0 = (double *)malloc((D->nSpecies)*sizeof(double ));
    LL=D->loadList;
    while(LL->next)    {
      if(LL->charge<0)	charge=-1.0*LL->charge;
      else 	        charge=LL->charge;
      rho0[s]=1.0*LL->density;
      LL=LL->next;
      s++;
    }

    nx=D->nx;  ny=D->ny;    nz=D->nz;

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=D->minZSub-D->minZDomain;

    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%ddensity%d.h5",s,iteration);

      if(myrank==0)        {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
        saveIntMeta(name,"/nx",&nx,1);
        saveIntMeta(name,"/ny",&ny,1);
        saveIntMeta(name,"/nz",&nz,1);
      } else      ;
      MPI_Barrier(MPI_COMM_WORLD);

      switch(D->dimension) {
      //1D
      case 1:
        solveDensity1D(D,s,rho0[s]);
        if(D->L>1)  {
          MPI_TransferRho_Xplus(D,D->Rho,1,1,3);
          MPI_TransferRho_Xminus(D,D->Rho,1,1,3);
        }  else     ;

        sprintf(dataName,"%d",s);
        saveFieldComp(D->Rho,name,dataName,nx,1,1,nxSub,1,1,istart,iend,0,1,0,1,offset);

        if(myrank==0)   {
          saveCoordHDF(D,name);
          sprintf(fileName,"%ddensity%d",s,iteration);
          density_xdmf(1,fileName,nx,ny,nz,s);
          printf("%s\n",name);
        }  else  ;
        break;

      //2D
      case 2:
        solveDensity2D(D,s,rho0[s]);
        if(D->L>1)  {
          MPI_TransferRho_Xplus(D,D->Rho,nySub+5,1,3);
          MPI_TransferRho_Xminus(D,D->Rho,nySub+5,1,3);
        }  else     ;
        if(D->M>1)  {
          MPI_TransferRho_Yplus(D,D->Rho,nxSub+5,1,3);
          MPI_TransferRho_Yminus(D,D->Rho,nxSub+5,1,3);
        }  else     ;

        sprintf(dataName,"%d",s);
        saveFieldComp(D->Rho,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

        if(myrank==0)   {
          saveCoordHDF(D,name);
          sprintf(fileName,"%ddensity%d",s,iteration);
          density_xdmf(2,fileName,nx,ny,nz,s);
          printf("%s\n",name);
        }  else  ;

        break;
      }   //End of switch (dimension)

    }		//End of for(s)
    free(rho0);
}

void solveDensity1D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  double Wx[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,weight;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;

  j=k=0;
  for(i=0; i<iend+3; i++)
      D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        weight=p->weight;
        x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        for(ii=0; ii<4; ii++)
            D->Rho[i-1+ii][j][k]+=Wx[ii]*coef*weight;
        p=p->next;
      }
    }
}   

void solveDensity2D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  double Wx[4],Wy[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,weight;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;

  k=0;
  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        weight=p->weight;
        x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        y=p->y; y1=1+y; y2=y; y3=1-y; y4=2-y;
        Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
        Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
        Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
        Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
        for(ii=0; ii<4; ii++)
          for(jj=0; jj<4; jj++)
            D->Rho[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*coef*weight;
        p=p->next;
      }
    }
}   

void solveDensity3D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,kk,istart,iend,jstart,jend,kstart,kend;
  double Wx[4],Wy[4],Wz[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4,weight;;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;

  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      for(k=0; k<kend+3; k++)
        D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {
          weight=p->weight;
          x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
          Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
          Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
          Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
          Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
          y=p->y; y1=1+y; y2=y; y3=1-y; y4=2-y;
          Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
          Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
          Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
          Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
          z=p->z; z1=1+z; z2=z; z3=1-z; z4=2-z;
          Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
          Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
          Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
          Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;
          for(ii=0; ii<4; ii++)
            for(jj=0; jj<4; jj++)
              for(kk=0; kk<4; kk++)
                D->Rho[i-1+ii][j-1+jj][k-1+kk]+=Wx[ii]*Wy[jj]*Wz[kk]*coef*weight;
          p=p->next;
        }
      }
}  

void saveCoordHDF(Domain *D,char *fileName)
{
  int ii,i,nx,ny,nz;
  char name[100];
  double *xtic,*ytic,*ztic;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,tic_id;
  herr_t status;
  hid_t filespace;
  hsize_t dimy[1],dimx[1],dimz[1];
  const char *coorName[] = {"/Y","/X","/Z"};

  nx=D->nx;  ny=D->ny;  nz=D->nz;

  if(myrank==0)
  {
    sprintf(name,"%s",fileName);
    file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);

    switch(D->dimension) {
    //1D
    case 1:
      dimx[0]=nx;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;

      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate2(file_id,coorName[1],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);

      free(xtic);
      break;

    //2D
    case 2:
      dimx[0]=nx;      dimy[0]=ny;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++) xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++) ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      for(ii=0; ii<2; ii++)
      {
        if(ii==0) filespace=H5Screate_simple(1,dimy,NULL);
        else if(ii==1) filespace=H5Screate_simple(1,dimx,NULL);
        dset_id=H5Dcreate2(file_id,coorName[ii],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ii==0 ? ytic : xtic);
        H5Dclose(dset_id);
        H5Sclose(filespace);
      }
      free(xtic);
      free(ytic);
      break;

    //3D
    case 3:
      dimx[0]=nx;      dimy[0]=ny;      dimz[0]=nz;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++) xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++) ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      ztic=(double *)malloc(nz*sizeof(double));
      for(i=0;i<nz;i++) ztic[i]=(i+D->minZDomain)*D->lambda*D->dz;
      filespace=H5Screate_simple(1,dimy,NULL);
      dset_id=H5Dcreate2(file_id,"/Y",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ytic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate2(file_id,"/X",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);

      free(xtic);
      free(ytic);
      free(ztic);
      break;
    }
    H5Fclose(file_id);
  }
  else ;
}

// The number of cells in the X, Y dimensions
void density_xdmf(int dimension,char *fileName,int nx,int ny,int nz,int s)
{
    FILE *xmf = 0;
    char name[100];
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    switch (dimension)  {
    //1D
    case 1 :
      fprintf(xmf, "     <Topology TopologyType=\"1DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VX\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");

      fprintf(xmf, "     <Attribute Name=\"/%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/%d\n",fileName,s);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;

    //2D
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");

      fprintf(xmf, "     <Attribute Name=\"/0%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/%d\n",fileName,s);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
/*
    case 3 :
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",ny,nx,nz);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        %s%d.h5:/Z\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s%d.h5:/X\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s%d.h5:/Y\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(s=0; s<nSpecies; s++)
      {
        fprintf(xmf, "     <Attribute Name=\"%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s%d.h5:/%d\n",name,iteration,s);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
*/
    }
}
                                                          
void setZero(double ***den,int nx, int ny, int nz)
{
   int i,j,k;

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         den[i][j][k]=0.0;
}

 
