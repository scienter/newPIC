#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void removeEdge(Domain D)
{
  void removeEdge1D();
  void removeEdge2D();
  void removeEdge3D();
  
  switch(D.dimension)  {
  case 1 :
    removeEdge1D(&D);
    break;
  case 2 :
    removeEdge2D(&D);
    break;
  case 3 :
    removeEdge3D(&D);
    break;
  default :
    printf("In removeEdge, whate dimension(%d)?\n",D.dimension);
  }
}

/*
void removeEdge2DBoost(Domain *D)
{
    int i,j,istart,iend,jstart,jend,s;
    Particle **particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //remove left boundary
    for(i=0; i<istart; i++)
    {
      for(s=0; s<D->nSpecies; s++)
        for(j=jstart-1; j<=jend; j++)
        { 
          p=particle[i][j].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
        }
    }
    

    //remove right boundary
    i=iend;
    for(s=0; s<D->nSpecies; s++)
      for(j=jstart-1; j<=jend; j++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

    //remove up boundary
    j=jend;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

    //remove bottom boundary
    j=jstart-1;
    for(s=0; s<D->nSpecies; s++)
      for(i=0; i<=iend; i++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

}
*/
void removeEdge1D(Domain *D)
{
    int i,j,k,istart,iend,s;
    Particle ***particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    istart=D->istart;
    iend=D->iend;

    //remove left boundary
    j=k=0;
    for(i=0; i<istart; i++)
      for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            if(p!=NULL)
              free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove right boundary
    i=iend;
    for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {	
          tmp=p->next;
          particle[i][j][k].head[s]->pt=tmp; 
          p->next=NULL;
          free(p);
          p=particle[i][j][k].head[s]->pt;
        }
      }

}


void removeEdge2D(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s,cnt;
    Particle ***particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    //remove left boundary
    k=0;
    for(i=0; i<istart; i++)
      for(j=jstart-1; j<jend+1; j++)
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p!=NULL)
          {
            particle[i][j][k].head[s]->pt=p->next; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove right boundary
    i=iend;
    k=0;
    for(j=jstart-1; j<=jend; j++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {	
          particle[i][j][k].head[s]->pt=p->next; 
          p->next=NULL;
          free(p);
          p=particle[i][j][k].head[s]->pt;
        }
      }

    //remove plusY boundary
    j=jend;
    k=0;
    for(i=0; i<=iend; i++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {	
          particle[i][j][k].head[s]->pt=p->next; 
          p->next=NULL;
          free(p);
          p=particle[i][j][k].head[s]->pt;
        }
      }

    //remove minusY boundary
    j=jstart-1;
    k=0;
    for(i=0; i<=iend; i++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j][k].head[s]->pt; 
        while(p)
        {	
          particle[i][j][k].head[s]->pt=p->next; 
          p->next=NULL;
          free(p);
          p=particle[i][j][k].head[s]->pt;
        }
      }
}

void removeEdge3D(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s;
    Particle ***particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    //remove left boundary
    for(i=0; i<istart; i++)
    {
      for(j=jstart-1; j<=jend; j++)
        for(k=kstart-1; k<=kend; k++)
          for(s=0; s<D->nSpecies; s++)
          { 
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {	
              tmp=p->next;
              particle[i][j][k].head[s]->pt=tmp; 
              p->next=NULL;
              free(p);
              p=particle[i][j][k].head[s]->pt;
            }
          }
    }

    //remove right boundary
    i=iend;
    for(j=jstart-1; j<=jend; j++)
      for(k=kstart-1; k<=kend; k++)
        for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove plusZ boundary
    k=kend;
    for(i=0; i<=iend; i++)
      for(j=jstart-1; j<=jend; j++)
        for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove minusZ boundary
    k=kstart-1;
    for(i=0; i<=iend; i++)
      for(j=jstart-1; j<=jend; j++)
        for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove plusY boundary
    j=jend;
    for(i=0; i<=iend; i++)
      for(k=kstart-1; k<=kend; k++)
        for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

    //remove minusY boundary
    j=jstart-1;
    for(i=0; i<=iend; i++)
      for(k=kstart-1; k<=kend; k++)
        for(s=0; s<D->nSpecies; s++)
        { 
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j][k].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j][k].head[s]->pt;
          }
        }

}
