

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    double x; 
    double oldX;   
    double y; 
    double oldY;   
    double z; 
    double oldZ;   
    double p1;    //momentum  
    double p2;
    double p3;
    double E1;    
    double E2;    
    double E3;    
    double B1;    
    double B2;    
    double B3;    
    double weight;    
    int index; 
    int core; 
    struct _ptclList *next;
} ptclList;

