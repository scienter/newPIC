#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void saveDump(Domain D,int iteration)
{
  void saveDumpFieldHDF();
  void saveDumpParticleHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((Pukhov-1)*3+2) :
    if(D.saveDumpMode==HDF)  {
      saveDumpFieldHDF(&D,iteration);
      saveDumpParticleHDF(&D,iteration);
      if(myrank==0)  {
        printf("dumpField%d.h5\n",iteration);
        printf("dumpParticle%d.h5\n",iteration);
      }   else	;
    }
    break;
  }

}

void saveJDump(Domain D,int iteration)
{
  void saveJDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((Pukhov-1)*3+2) :
    if(D.saveDumpMode==HDF)  {
      saveJDumpHDF(&D,iteration);
      if(myrank==0) 
        printf("resolJ%d.h5\n",iteration);
    }
    break;
  }
}


void saveBDump(Domain D,int iteration)
{
  void saveBDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((Pukhov-1)*3+2) :
    if(D.saveDumpMode==HDF)  {
      saveBDumpHDF(&D,iteration);
      if(myrank==0) 
        printf("resolB%d.h5\n",iteration);
    }
    break;
  }
}

void saveEDump(Domain D,int iteration)
{
  void saveEDumpHDF();
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  switch ((D.fieldType-1)*3+D.dimension)  {
  case ((Pukhov-1)*3+2) :
    if(D.saveDumpMode==HDF)  {
      saveEDumpHDF(&D,iteration);
      if(myrank==0) 
        printf("resolE%d.h5\n",iteration);
    }
    break;
  }
}

void saveDumpFieldHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,offSetY,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    int offset[3];
    char name[100],name2[100];

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;
    void saveFieldComp();

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=nxSub+2;
    jstart=D->jstart;
    jend=nySub+2;
    kstart=D->kstart;
    kend=nzSub+2;

    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/nSpecies",&D->nSpecies,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
    } else	;
    
//    MPI_Barrier(MPI_COMM_WORLD);
    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }	//End of switch (dimension)
}

void saveDumpParticleResolHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub,cnt,totalCnt,index,start,ii,n,tmp;
    int resolX,resolY,resolZ,minXSub,minYSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    char name[100],name2[100];
    double *saveDouble;
    int *recv,*saveInt,*offSetRank;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    recv = (int *)malloc(nTasks*sizeof(int ));
    offSetRank = (int *)malloc(nTasks*sizeof(int ));

    hid_t file_id,group_id;
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

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      sprintf(name,"resolParticle%d.h5",iteration);
      if(myrank==0)   {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
      }      else ;

      if(myrank==0)  {
        saveIntMeta(name,"/nSpecies",&D->nSpecies,1);
        saveIntMeta(name,"/nx",&D->nx,1);
        saveIntMeta(name,"/ny",&D->ny,1);
        saveIntMeta(name,"/nz",&D->nz,1);
        saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
        saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
        saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
        saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
      } else	;
      MPI_Barrier(MPI_COMM_WORLD);

      istart=D->istart;
      iend=D->iend;
      jstart=D->jstart;
      jend=D->jend;
      minXSub=D->minXSub;
      minYSub=D->minYSub;

      for(s=0; s<D->nSpecies; s++)
      {

        if(myrank==0)    {
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(name2,"%dParticle",s);
          group_id=H5Gcreate2(file_id,name2,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }  else 	;
        MPI_Barrier(MPI_COMM_WORLD);

        k=0;
        cnt=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              cnt++;
              p=p->next;
            }
          }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];
        totalCnt=0;
        for(i=0; i<nTasks; i++)
          totalCnt+=recv[i];
        saveDouble = (double *)malloc(cnt*sizeof(double ));
        saveInt = (int *)malloc(cnt*sizeof(int ));

        for(i=0; i<nTasks; i++)   {
          tmp=0;
          for(ii=0; ii<i; ii++)
            tmp+=recv[ii];
          offSetRank[i]=tmp;
        }

        if(myrank==0)  {
          sprintf(name2,"%dParticle/offSet",s);
          saveIntMeta(name,name2,offSetRank,nTasks);
          sprintf(name2,"%dParticle/totalCnt",s);
          saveIntMeta(name,name2,&totalCnt,1);
        }    else    ;

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->oldX-istart+D->minXSub;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/x",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->oldY-jstart+D->minYSub;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/y",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->p1;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/px",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->p2;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/py",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->p3;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/pz",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveDouble[index]=p->weight;
              index++;
              p=p->next;
            }
          }
        sprintf(name2,"%dParticle/weight",s);
        saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)     {
            p=particle[i][j][k].head[s]->pt;
            while(p)      {
              saveInt[index]=p->index;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/index",s);
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)     {
            p=particle[i][j][k].head[s]->pt;
            while(p)    {
              saveInt[index]=p->core;
              p=p->next;
              index++;
            }
          }
        sprintf(name2,"%dParticle/core",s);
        saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);


        free(saveDouble);
        free(saveInt);

      }	//End of species

    }	//End of switch (dimension)

    if(myrank==0)  printf("%s\n",name);  else	;

    free(recv);
    free(offSetRank);
}

void saveDumpParticleHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
    int nxSub,nySub,nzSub,cnt,totalCnt,index,start,ii,n,tmp;
    int resolX,resolY,resolZ,minXSub,minYSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    char name[100],name2[100];
    double *saveDouble;
    int *recv,*saveInt,*offSetRank;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    recv = (int *)malloc(nTasks*sizeof(int ));
    offSetRank = (int *)malloc(nTasks*sizeof(int ));

    hid_t file_id,group_id;
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

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      sprintf(name,"dumpParticle%d.h5",iteration);
      if(myrank==0)   {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
      }      else ;

      if(myrank==0)  {
        saveIntMeta(name,"/nSpecies",&D->nSpecies,1);
        saveIntMeta(name,"/nx",&D->nx,1);
        saveIntMeta(name,"/ny",&D->ny,1);
        saveIntMeta(name,"/nz",&D->nz,1);
        saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
        saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
        saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      } else	;
      MPI_Barrier(MPI_COMM_WORLD);

      istart=D->istart;
      iend=D->iend;
      jstart=D->jstart;
      jend=D->jend;
      minXSub=D->minXSub;
      minYSub=D->minYSub;

      for(s=0; s<D->nSpecies; s++)
      {

        if(myrank==0)    {
          file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);
          sprintf(name2,"%dParticle",s);
          group_id=H5Gcreate2(file_id,name2,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5Gclose(group_id);
          H5Fclose(file_id);
        }  else 	;
        MPI_Barrier(MPI_COMM_WORLD);

        k=0;
        cnt=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              cnt++;
              p=p->next;
            }
          }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];
        totalCnt=0;
        for(i=0; i<nTasks; i++)
          totalCnt+=recv[i];
        saveDouble = (double *)malloc(cnt*sizeof(double ));
        saveInt = (int *)malloc(cnt*sizeof(int ));

        if(myrank==0)  {
          sprintf(name2,"%dParticle/totalCnt",s);
          saveIntMeta(name,name2,&totalCnt,1);
        }    else    ;

        if(totalCnt>0)
        {
          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->x+i-istart+D->minXSub;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/x",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->y+j-jstart+D->minYSub;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/y",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->p1;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/px",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->p2;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/py",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->p3;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/pz",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)    {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveDouble[index]=p->weight;
                index++;
                p=p->next;
              }
            }
          sprintf(name2,"%dParticle/weight",s);
          saveParticleComp_Double(saveDouble,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)     {
              p=particle[i][j][k].head[s]->pt;
              while(p)      {
                saveInt[index]=p->index;
                p=p->next;
                index++;
              }
            }
          sprintf(name2,"%dParticle/index",s);
          saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);

          index=0;
          for(i=istart; i<iend; i++)
            for(j=jstart; j<jend; j++)     {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                saveInt[index]=p->core;
                p=p->next;
                index++;
              }
            }
          sprintf(name2,"%dParticle/core",s);
          saveParticleComp_Int(saveInt,name,name2,totalCnt,cnt,start);
        }	//End of totalCnt>0
        else	;

        free(saveDouble);
        free(saveInt);

      }	//End of species
    }	//End of switch (dimension)

    if(myrank==0)  printf("%s\n",name);  else	;

    free(recv);
    free(offSetRank);
}

void saveEDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=nxSub+2;
    jstart=D->jstart;
    jend=nySub+2;
    kstart=D->kstart;
    kend=nzSub+2;

    metaDim[0]=1;

    sprintf(name,"resolE%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    } else	;
    MPI_Barrier(MPI_COMM_WORLD);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Ex,name,"/Ex",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ey,name,"/Ey",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Ez,name,"/Ez",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}

void saveJDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();
    void saveIntMeta();

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=nxSub+2;
    jstart=D->jstart;
    jend=nySub+2;
    kstart=D->kstart;
    kend=nzSub+2;

    metaDim[0]=1;

    sprintf(name,"resolJ%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    }	else	;
    MPI_Barrier(MPI_COMM_WORLD);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Jx,name,"/Jx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jy,name,"/Jy",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Jz,name,"/Jz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}


void saveBDumpHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    double x,y,z,y1,y2,y3,a,b,c,dx,dy,dz;
    int offset[3];
    char name[100],name2[100];
    int *saveInt;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    hsize_t metaDim[1];
    herr_t status;
    void saveFieldComp();
    void saveIntMeta();

    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    istart=D->istart;
    iend=nxSub+2;
    jstart=D->jstart;
    jend=nySub+2;
    kstart=D->kstart;
    kend=nzSub+2;

    metaDim[0]=1;

    sprintf(name,"resolB%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D->minZDomain,1);
      saveIntMeta(name,"/nx",&D->nx,1);
      saveIntMeta(name,"/ny",&D->ny,1);
      saveIntMeta(name,"/nz",&D->nz,1);
      saveDoubleMeta(name,"/dtRatio",&D->dtRatio,1);
    }	else	;
    MPI_Barrier(MPI_COMM_WORLD);

    switch((D->fieldType-1)*3+D->dimension) {
    //2D
    case (Pukhov-1)*3+2:
      nx=D->nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
      ny=D->ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

      offset[0]=(D->minXSub-D->minXDomain)+biasX;
      offset[1]=(D->minYSub-D->minYDomain)+biasY;
      offset[2]=0;
      
      saveFieldComp(D->Bx,name,"/Bx",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->By,name,"/By",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      saveFieldComp(D->Bz,name,"/Bz",nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      break;
    }
}

void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
//  if(cnt>0)
    memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

//  if(cnt>0)
    H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet)
{
  int i,j,k;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id,tic_id;
  herr_t status;
  hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
  hsize_t dimsf[1],count[1],offset[1];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);

  dimsf[0]=totalCnt;
  filespace=H5Screate_simple(1,dimsf,NULL);

  count[0]=cnt;
  offset[0]=offSet;
//  if(cnt>0)
    memspace=H5Screate_simple(1,count,NULL);

  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);
 
//  if(cnt>0)
    H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L)
{
  if(L==1) {
    *istart=0;
    *iend+=3;
    *biasX=0;
    *nxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0;
      *nxSub+=2;
      *biasX=0;
    }  else if(rankX==L-1)  {
      *iend+=3;
      *nxSub+=3;
      *biasX=2;
    } else
      *biasX=2;
  }
}

void deleteFieldData(Domain *D,double ***field,int *nxSub,int *nySub,int *nzSub,int rankX,int rankY)
{
   int i,j,k;

   if(rankX==0)  *nxSub-=2;
   else if(rankX==D->L-1)  *nxSub-=3;
   if(rankY==0)  *nySub-=2;
   else if(rankY==D->M-1)  *nySub-=3;

   for(i=0; i<*nxSub+5; i++)
   {
     for(j=0; j<*nySub+5; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
    int ii,i,j,k,start;
    double *field;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=ny;
    dimsf[1]=nx;
    dimsf[2]=nz;
    filespace=H5Screate_simple(3,dimsf,NULL);

    count[0]=nySub;
    count[1]=nxSub;
    count[2]=nzSub;
    offset[0]=offSet[1];
    offset[1]=offSet[0];
    offset[2]=offSet[2];
    memspace=H5Screate_simple(3,count,NULL);

    field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=data[i][j][k];
          start+=nzSub;     
        }  

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void MPI_saveIntArray(int *data,char *fileName,char *dataName,int offSet)
{
    int i;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[1],count[1],offset[1];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=nTasks;
    filespace=H5Screate_simple(1,dimsf,NULL);

    count[0]=1;
    offset[0]=offSet;
    memspace=H5Screate_simple(1,count,NULL);

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      

