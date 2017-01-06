typedef struct _Domain
{
  int dimension;
  int mode;
  int step;
  int L,M,N;
  int resolX,resolY,resolZ;
  int rankX,rankY,rankZ;

  int nSpecies;
  int nx,ny,nz;

  int nxSub,nySub,nzSub;
  int minXSub,minYSub,minZSub,maxXSub,maxYSub,maxZSub;
  int reMinXSub,reMinYSub,reMinZSub,reMaxXSub,reMaxYSub,reMaxZSub;
  int minXDomain,minYDomain,minZDomain;
  int maxXDomain,maxYDomain,maxZDomain;

  int istart,jstart,kstart,iend,jend,kend;

  //resize
  int reNx,reNy,reNz;
  int reNxSub,reNySub,reNzSub;
  int reIend,reJend,reKend;
  int minX,maxX,minY,maxY,minZ,maxZ;
  int reMinX,reMaxX,reMinY,reMaxY,reMinZ,reMaxZ;
  int saveIstart,saveIend,saveNxSub,biasX;
  int saveJstart,saveJend,saveNySub,biasY;
  int saveKstart,saveKend,saveNzSub,biasZ;
  int reSaveIstart,reSaveIend,reSaveNxSub;
  int reSaveJstart,reSaveJend,reSaveNySub;
  int reSaveKstart,reSaveKend,reSaveNzSub;

  int *recv,*recvDataCnt,*sharePNum,*coreCnt;
  int *minXSubList,*minYSubList,*minZSubList;
  int *maxXSubList,*maxYSubList,*maxZSubList;
  int *reMinXSubList,*reMinYSubList,*reMinZSubList;
  int *reMaxXSubList,*reMaxYSubList,*reMaxZSubList;
  double **recvData,**sendData;

  //particle momory
  int cntSub,totalCnt,numberInCell;
  int reCntSub,reTotalCnt;
  double *dataX,*dataY,*dataZ; 
  double *dataPx,*dataPy,*dataPz; 
  int *dataIndex,*dataCore,*dataCores;
  double *reDataX,*reDataY,*reDataZ; 
  int *reDataIndex,*reDataCore,*reDataCores;
  struct _Particle ***particle;

} Domain;

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _Particle
{
   ptclHead **head;
}  Particle;

typedef struct _ptclList  {
    double x;
    double y;
    double z;
    double p1;    //momentum  
    double p2;
    double p3;
    int index;
    int core;
    double weight;
    struct _ptclList *next;
} ptclList;


void restoreIntMeta(char *fileName,char *dataName,int *data);
double ***memoryAsign(int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);
void deleteParticleMemory(Domain *D,int s);
void highParticleSave(Domain *D,char *fileName,int s);
void cleanSubBoundary(Domain *D);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void restoreDen(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s,char *dataName,int step,int *offset);
void findAveragePx(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s);
void findAveragePy(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s);
void findAveragePz(Domain *D,double ***den1,double ***den2,double ***den3,double ***den4,int s);
