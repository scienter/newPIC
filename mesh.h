#include "particle.h"
#include "laser.h"
#include "plasma.h"

#define FIRST 	1
#define SECOND 	2
#define THIRD 	3

#define ON	1
#define OFF	0

#define TXT	0
#define HDF	1

#define Split	1
#define Yee	2
#define Pukhov	3

#define UP		1
#define DOWN		2
#define FRONT		3
#define BACK		4
#define UPFRONT		5
#define UPBACK		6
#define DOWNFRONT	7
#define DOWNBACK	8

typedef struct _Domain 
{
   int dimension;

   int fieldType;
   int currentType;
   int interpolationType;

   int maxTime;
   int saveFieldMode;
   int saveParticleMode;
   int saveDensityMode;
   int saveDumpMode;
   int maxStep;
   int saveStep;
   int centerStep;
   int saveStart;
   int dumpStart;
   int dumpSave;
   int dumpStep;
   int fieldSave;
   int ramanSave;
   int particleSave;
   int densitySave;

   int nx;             //Total domain
   int ny;             //Total domain
   int nz;             //Total domain
   int nxSub;          //Each core has sub domain        
   int nySub;          //Each core has sub domain        
   int nzSub;          //Each core has sub domain        
   int istart;
   int iend;
   int jstart;
   int jend;
   int kstart;
   int kend;
   //Each core has start mesh point in total domain
   int minXSub,maxXSub,minYSub,maxYSub,minZSub,maxZSub;     
   int minXDomain,minYDomain,minZDomain;
   int numberInCell;
   int moving;         //Moving domain option. 1:on
   int shiftDuration;

   double lambda;
   double omega;
   double divisionLambda;
   double dt;
   double dtRatio;
   double dx;
   double dy;
   double dz;
   int resolChange;
   int resolHigh;
   int resolLow;
   int resolStep;
   int resolX;
   int resolY;
   int resolZ;
    

   //MPI parameter
   int L;
   int M;
   int N;
   int nextXrank;
   int prevXrank;
   int nextYrank;
   int prevYrank;
   int nextZrank;
   int prevZrank;
   
   //sharing mesh
   double *XplusJ,*XminusJ,*YplusJ,*YminusJ,*ZplusJ,*ZminusJ;
   double *minusDenY,*plusDenY,*minusDenZ,*plusDenZ;
   int numPlusXJ,numMinusXJ,numPlusYJ,numMinusYJ,numPlusZJ,numMinusZJ;
   int numPlusDenY,numMinusDenY,numPlusDenZ,numMinusDenZ;

   double ***Rho;    
   double ***Ex;    
   double ***Ey;    
   double ***Ez;    
   double ***Bx;    
   double ***By;    
   double ***Bz;    
   double ***BxNow;    
   double ***ByNow;    
   double ***BzNow;    
   double ***Jx;    
   double ***Jy;    
   double ***Jz;    

   //split mode
   double ***Pr;    
   double ***Pl;    
   double ***Sr;    
   double ***Sl;
   double ***ExC;
   double ***BxC;      
   double ***PrC;    
   double ***PlC;    
   double ***SrC;    
   double ***SlC;    
   double ***JxOld;    
   double ***JyOld;    
   double ***JzOld;    
   
   struct _Particle ***particle;    
   struct _Boost **boost;    

   //Plasma load
   struct _LoadList *loadList;
   int nSpecies;

   //Plasma load
   struct _LaserList *laserList;
   int nLaser;

   //Boost
   int boostOn;
   int boostIon;
   double gamma;
   double beta;
   int minT;	//boost frame's step
   int maxT;	//boost frame's step
   int boostSaveStep;	//lab frame's step
   
   //Probe
   int probeNum;
   int *probeX;
   int *probeY;
   int *probeZ;
   struct _Probe **probe;
   
   //ID Track
   int tracking;
   int trackSaveStep;
   int trackStart;
   int idNums;
   int *trackID;
   int *trackCore;
   int *trackS;
   struct _Track **track;

   //PML
   int pmlOn;
   int pmlCellRight;   
   int pmlCellLeft;   
   int pmlCellUp;   
   int pmlCellBottom;   
   double pmlr;
   double pmld;

}  Domain; 

typedef struct _Boost
{
   double x;
   double y;
   double E1;
   double B1;
   double Pr;
   double Pl;
   double Sr;
   double Sl;   
}  Boost;

typedef struct _UPML 
{
   double ***PrC;
   double ***PlC;
   double ***SrC;
   double ***SlC;
   double ***ExC;
   double ***BxC;
   double ***Pr;
   double ***Pl;
   double ***Sr;
   double ***Sl;
   double ***Ex;
   double ***Ey;
   double ***Ezx;
   double ***Ezy;
   double ***Bx;
   double ***By;
   double ***Bzx;
   double ***Bzy;
}  UPML;

typedef struct _DPML 
{
   double ***PrC;
   double ***PlC;
   double ***SrC;
   double ***SlC;
   double ***ExC;
   double ***BxC;
   double ***Pr;
   double ***Pl;
   double ***Sr;
   double ***Sl;
   double ***Ex;
   double ***Ey;
   double ***Ezx;
   double ***Ezy;
   double ***Bx;
   double ***By;
   double ***Bzx;
   double ***Bzy;
}  DPML;

typedef struct _Particle 
{
   // Particle List Header
   ptclHead **head;            
}  Particle;

typedef struct _External 
{
   double E1;
   double E2;
   double E3;
   double B1;
   double B2;
   double B3;
}  External;

typedef struct _Probe
{
   double E1;
   double Pr;
   double Pl;
   double B1;
   double Sr;
   double Sl;
}  Probe;

typedef struct _Track
{
   double x;
   double y;
   double z;
   double px;
   double py;
   double pz;
   int step;
   int id;
   int core;
   double wp;
   double kp;
}  Track;

void cleanMemory(Domain *D);
void saveTracking(Domain *D);
void removeEdge(Domain D);
void particleShareZ(Domain *D);
void particleShareY(Domain D);
void particleShareX(Domain D);
void rearrangeParticles(Domain *D);
void movingDomain(Domain *D,int iteration);
void loadMovingPlasma(Domain *D,LoadList *LL,int s,int iteration);
void updateCurrent(Domain D);
void particlePush(Domain *D);
void interpolation(Domain *D,External *Ext);
void fieldSolve(Domain D,double t);
void loadLaser(Domain *D,LaserList *L,double t);
void saveDump(Domain D,int iteration);
void saveBDump(Domain D,int iteration);
void saveEDump(Domain D,int iteration);
void saveJDump(Domain D,int iteration);
void saveDumpParticleHDF(Domain *D,int iteration);
void saveDumpParticleResolHDF(Domain *D,int iteration);
void saveDumpDensityResolHDF(Domain D,int iteration);
void saveDenParticleHDF(Domain *D,int iteration);
void saveFile(Domain D,int iteration);
void firstFieldShare(Domain D);
void secondFieldShare(Domain D);
void trackID(Domain *D,int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend);
void loadPlasma(Domain *D,LoadList *LL,int s,int iteration);
void restoreDump(Domain D,int iteration);
void boundary(Domain *D,External *Ext);
void parameterSetting(Domain *D,External *Ext, char *input);
int FindParameters (char *block, int rank, char *options, char *input, char *ret);
void saveFieldHDF(Domain D,int iteration);
void saveCenterDensity(Domain *D,int iteration);
void saveCenterField(Domain *D,int iteration);
double ***memoryAsign(int nx, int ny, int nz);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

