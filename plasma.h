#define Electron 	1
#define HPlus0 	 	100
#define HPlus1 	 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0          600
#define CPlus1          601
#define CPlus2          602
#define CPlus3          603
#define CPlus4          604
#define CPlus5          605
#define CPlus6          606
#define AlPlus11        1311
#define userDefined   	9999999

#define Polygon    	1
#define Defined    	2
#define Channel    	3
#define BoostFrame	4
#define Circle    	5
#define Exp    		6

#define Constant   	0
#define Gaussian   	1
#define Polynomial   	2

#define byNumber	0
#define byDensity	1

typedef struct _LoadList  {
   int type;
   int species;
   double superP;
   double density;
   double numberInCell;
   double criticalDensity;
   int index;  
   double num;      //exceeded number of particle which is less than 1
   int xnodes;     //longitudinal point number
   double *xn;      //longitudinal density (1 is P->density)
   double *xpoint;    //longitudinal point
   int ynodes;     //transverse point number
   double *yn;      //transverse density (1 is P->density)
   double *ypoint;    //transverse point
   int znodes;     //transverse point number
   double *zn;      //transverse density (1 is P->density)
   double *zpoint;    //transverse point
   double givenMinPx;	//for saveParticle option

   //defined plasma
   int defineMode;
   double *xPosition;
   double *yPosition;
   double *zPosition;
   int numDefined;
   double xLengDef;
   double yLengDef;
   double zLengDef;
   int numDefPtcls;
   double **define;
   int maxLoadTime;
   int minLoadTime;
   double minX;
   double maxX;
   double minY;
   double maxY;
   double minZ;
   double maxZ;

   int pointPosition;   
   double p1;
   double p2;
   double p3;

   double mass;
   int charge;
   
   double temperature;
   
   //applying function
   double centerX;
   double centerY;
   double centerZ;
   double gaussCoefX;
   double polyCoefX;
   double gaussCoefYZ;
   double polyCoefYZ;
   int modeX;
   int modeYZ;

   struct _LoadList *next;
} LoadList;
