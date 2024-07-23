#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

cg_parameter UParm1; 

char * Output_File_Name;
char * Input_File_Name;

INT N;     // number of circles
INT n;     // number of continous variables n=2*N
INT SR;    // success rate of algorithm
INT HN;    // The number of ploygon holes in the instance
double R, R_max, R_min, alpha;  
double pointlist[MaxVertice] = {0};  // the coordinates of vertices of container, (x1,y1,x2,y2, ...., x_numP,y_numP)
int numPH[MaxNholes] = {0};          // The numbers of vertices of ploygon holes in the instance
double holepoint[MaxNholes][MaxVertice] = {0};  // holepoint[i] denotes the coordinates of vertices of i-th ploygon hole, (x1,y1,x2,y2, ...., x_numPH[i],y_numPH[i]) 
INT numP;  // number of vertices in the ploygon contianer  

double	VacancyRadius; 
double final_time, starting_time, Time_limit;

double X_min, X_max, Y_min, Y_max; 

Site VacancySite[MaxVS]; 
int CircleSite[MaxNC];  
int *TabuTenure;

Solution CS;  
Solution OS;
Solution NS;
Solution VS; 
Solution CCS;
Solution CBestS;
Solution GlobalS;
Solution TBestS; 
Fsolution FS;

int **Adjacent;      // The adjacency relations between circles  
int *NNeighbors;     // The numbers of neighboring circles for each circle in the currrent solution 

int **AdjacentDD;    // The neighboring vertices of container for each cirlce in the current solution 
int *NNeighborsDD;   // The number of neighboring vertices for each cirlce 

int **AdjacentB;     // The neighboring edges in the container for each cirlce. For each edge <v1,v2>, the first vertex  v1 is stored, ometting the second vertex. 
int *NNeighborsB;    // Then number of neighboring edges in the container for each cirlce. 

int ***AdjacentDDH;  // The neighboring vertices of ploygon holes for each circle
int **NNeighborsDDH; // The numbers of neighboring vertice in each ploygon hole for each circle
 
int ***AdjacentBH;   // The neighboring edges of ploygon hole for each hole and each cirlce 
int **NNeighborsBH;  // The numbers of neighboring edges for each ploygon hole ans each circle 

