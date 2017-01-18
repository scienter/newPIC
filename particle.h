

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
    double p1Old2;    //momentum  
    double p2Old2;
    double p3Old2;
    double p1Old1;    //momentum  
    double p2Old1;
    double p3Old1;
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

