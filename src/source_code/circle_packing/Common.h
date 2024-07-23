#ifndef MaxNV
#define MaxNV 10000
#define MaxNN 100
#endif // the maximum number of varaibles in the objective function

#ifndef MaxVS
#define MaxVS  5  // the maximum number of vacancy sites in the configuration. 
#endif 

#ifndef MaxNC
#define MaxNC  5  // the maximum number of vacancy sites in the configuration. 
#endif 


#ifndef MaxVertice
#define MaxVertice 2000  // the maximum number of vertices in the ploygon container or the ploygon holes. 
#endif  

#ifndef MaxNholes
#define MaxNholes 20
#endif  

#ifndef pi
#define pi 3.1415926535898
#endif // the maximum number of varaibles in the objective function

typedef struct Configuration
{
	double f;
	double *x;
	double CurrentR;
}Solution;

typedef struct FlexibleConfiguration
{
	double L;
	double *x;
} Fsolution;

typedef struct LatticePoint
{
	int Index; 
	double x[2]; 
	double r; 
	double e; 
} Site; 

extern cg_parameter UParm1; 
extern INT n,N,SR,HN;
extern double VacancyRadius; 
extern Solution CS, NS, VS, CCS, OS, TBestS,  CBestS, GlobalS;
extern Fsolution FS;
extern char * Output_File_Name;
extern char * Input_File_Name;
extern double final_time, starting_time, Time_limit;
extern double R, R_max, R_min, alpha, VacancyRadius;
extern double X_min, X_max, Y_min, Y_max; 
extern Site VacancySite[MaxVS]; 
extern int CircleSite[MaxNC];  
extern double pointlist[MaxVertice]; 
extern int numPH[MaxNholes];
extern double holepoint[MaxNholes][MaxVertice];
extern INT numP;
extern int **Adjacent;
extern int *NNeighbors;
extern int **AdjacentDD;
extern int *NNeighborsDD;
extern int **AdjacentB;
extern int *NNeighborsB;
extern int ***AdjacentDDH;
extern int **NNeighborsDDH;
extern int ***AdjacentBH;
extern int **NNeighborsBH;
extern int *TabuTenure;
double myvalue
(
    double   *x,
    INT       n
) ;

double myvalue_sub
(
    double   *x,
    INT       n
) ;

void myvalue_each
(
    double   *x,
    INT       n,
    double   *f
);

void mygrad
(
    double    *g,
    double    *x,
    INT        n
) ;

void mygrad_sub
(
    double    *g,
    double    *x,
    INT        n
) ;

double myvalueNN
(
    double   *x,
    INT       n
) ;

void mygradNN
(
    double    *g,
    double    *x,
    INT        n
) ;

double myvalueRR
(
    double   *x,
    INT       n
);
void mygradRR
(
    double    *g,
    double    *x,
    INT        n
);

void TabuSearch(Solution &S, int SearchDepth, double rr);
void Range();
double Deviation(double arr[],int n);
void Outputing(Solution &S, int n) ;
void Outputing1(Solution &S, char *filename, int n);
void Outputing2(Fsolution &S, int n);
int  Inputing(char *filename) ;
void Out_results(double best , double ave,  double worst, int sr, double AvgTime, double deviation, char *filename, char *inputfile, int N);
void AssignMemery( );
void initial(Solution &S, INT n);
void LocalSearch(Solution &S);
void Shake(Solution &S, int n, double Max_delt);
void Find_Neighbor(double x[], double Dcut);
void Find_Neighbor2(double x[], int num, double y[], double Dcut);
double AdjustRadius(Fsolution &S);
void TSGO(); 
int  judin(double x_i,double y_i,int num,double *plist);
void MBH(Solution &S, double CR);
