#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void rearrangeParticles(Domain *D)
{
  void rearrangeParticles1D();
  void rearrangeParticles2D();
  void rearrangeParticles3D();
 
  switch(D->dimension)  {
  case 1 :
    rearrangeParticles1D(D);
    break;
  case 2 :
    rearrangeParticles2D(D);
    break;
  case 3 :
    rearrangeParticles3D(D);
    break;
  default :
    printf("In rearrangeParticles, what dimension(%d)?\n",D->dimension);
  }
}

void rearrangeParticles1D(Domain *D)
{
    Particle ***particle;
    particle=D->particle;

    int i,j,k,s,intX=0,intY=0,intZ=0,cnt,deleteFlag=0;
    int istart,iend;
    double x,y,z;
    ptclList *p,*New,*prev,*tmp;

    istart=D->istart;
    iend=D->iend;

    j=k=0;
    for(i=istart; i<iend; i++)
        for(s=0; s<D->nSpecies; s++)
        {
          cnt=1;
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            if(cnt==1)
              prev=p;
            deleteFlag=0;
              
            x=p->x;
            if(x>=1.0)  {
              intX=(int)x;
              x-=intX;
              deleteFlag=1;
            }
            else if(x<0) {              
              intX=(int)(x-1);
              x-=intX;
              deleteFlag=1;
            } 
            else   intX=0;

            if(deleteFlag==1)
            {
              if(cnt==1)
              {
                p->x=x;    
                particle[i][j][k].head[s]->pt = p->next;
                p->next = particle[i+intX][j][k].head[s]->pt;
                particle[i+intX][j][k].head[s]->pt = p;
                p=particle[i][j][k].head[s]->pt;
                cnt=1;
              }
              else
              {
                prev->next = p->next;
                p->x=x;    
                p->next = particle[i+intX][j][k].head[s]->pt;
                particle[i+intX][j][k].head[s]->pt = p;
                p=prev->next;
              }
            }		//End of if(deleteFlag==1)
            else
            {
              prev=p;
              p=p->next;
              cnt++;
            }              
          }		//End if while(p)
        }		//End of for(s)
}


void rearrangeParticles2D(Domain *D)
{
    Particle ***particle;
    particle=D->particle;

    int i,j,k,s,intX=0,intY=0,intZ=0,cnt,deleteFlag=0;
    int istart,iend,jstart,jend,kstart,kend;
    double x,y,z;
    ptclList *p,*New,*prev,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=intZ=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++)
        {
          cnt=1;
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            if(cnt==1)
              prev=p;
            else	;
            deleteFlag=0;
              
            x=p->x;   y=p->y;   z=p->z;
            if(x>=1.0)  {
              intX=(int)x;
              x-=intX;
              deleteFlag=1;
            }
            else if(x<0.0) {              
              intX=(int)(x-1.0);
              x-=intX;
              deleteFlag=1;
            } 
            else   intX=0;
            if(y>=1.0)  {
              intY=(int)y;
              y-=intY;
              deleteFlag=1;
            }
            else if(y<0.0) {
              intY=(int)(y-1.0);
              y-=intY;
              deleteFlag=1;
            } 
            else   intY=0;       

            if(deleteFlag==1)
            {
              if(cnt==1)
              {
                p->x=x;   p->y=y;   p->z=z;   
                particle[i][j][k].head[s]->pt = p->next;
                p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                p=particle[i][j][k].head[s]->pt;
                cnt=1;
              }
              else
              {
                p->x=x;   p->y=y;   p->z=z;  
                prev->next = p->next;
                p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                p=prev->next;
              }
            }		//End of if(deleteFlag==1)
            else
            {
              prev=p;
              p=p->next;
              cnt++;
            }              
          }		//End if while(p)
        }		//End of for(s)
}

void rearrangeParticles3D(Domain *D)
{
    Particle ***particle;
    particle=D->particle;

    int i,j,k,s,intX=0,intY=0,intZ=0,cnt,deleteFlag=0;
    int istart,iend,jstart,jend,kstart,kend;
    double x,y,z;
    ptclList *p,*New,*prev,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
          for(s=0; s<D->nSpecies; s++)
          {
            cnt=1;
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              if(cnt==1)
                prev=p;
              deleteFlag=0;
            
              x=p->x;   y=p->y;   z=p->z;
              if(x>=1.0)  {
                intX=(int)x;
                x-=intX;
                deleteFlag=1;
              }
              else if(x<0) {              
                intX=(int)(x-1.0);
                x-=intX;
                deleteFlag=1;
              } 
              else   intX=0;
              if(y>=1.0)  {
                 intY=(int)y;
                 y-=intY;
                 deleteFlag=1;
              }
              else if(y<0) {
                intY=(int)(y-1.0);
                y-=intY;
                deleteFlag=1;
              } 
              else   intY=0;       
              if(z>=1.0)  {
                 intZ=(int)z;
                 z-=intZ;
                 deleteFlag=1;
              }
              else if(z<0) {
                intZ=(int)(z-1.0);
                z-=intZ;
                deleteFlag=1;
              } 
              else   intZ=0;       
              if(deleteFlag==1)
              {
                if(cnt==1)
                {
                  particle[i][j][k].head[s]->pt = p->next;
                  p->x=x;   p->y=y;   p->z=z;   
                  p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                  particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                  p=particle[i][j][k].head[s]->pt;
                  cnt=1;
                }
                else
                {
                  prev->next = p->next;
                  p->x=x;   p->y=y;   p->z=z;  
                  p->next = particle[i+intX][j+intY][k+intZ].head[s]->pt;
                  particle[i+intX][j+intY][k+intZ].head[s]->pt = p;
                  p=prev->next;
                }
              }		//End of if(deleteFlag==1)
              else
              {
                prev=p;
                p=p->next;
                cnt++;
              }              
            }	//End if while(p)
          }			//End of for(s)
}

