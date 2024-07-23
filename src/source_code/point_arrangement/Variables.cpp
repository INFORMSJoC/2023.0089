#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

cg_parameter UParm1; 

char * Output_File_Name;
char * Input_File_Name;

INT N;     // number of points
INT n;     // number of continous variables n=2*N in the solution
INT SR;    // success rate of algorithm
INT HN;    // The number of ploygon holes in the instance
double R, R_max, R_min, alpha;  
double pointlist[MaxVertice] = {0};  // the coordinates of vertices of container, (x1,y1,x2,y2, ...., x_numP,y_numP)
int numPH[MaxNholes] = {0};          // The numbers of vertices of ploygon holes in the instance
double holepoint[MaxNholes][MaxVertice] = {0};  // holepoint[i] denotes the coordinates of vertices of i-th ploygon hole, (x1,y1,x2,y2, ...., x_numPH[i],y_numPH[i]) 
INT numP;  // number of vertices of the ploygon contianer  

double	VacancyRadius; // The radius of an additional detecting cirlce
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

int **Adjacent;      // The adjacency relations between points  
int *NNeighbors;     // The numbers of neighboring circles for each point in the currrent solution 

int **AdjacentDD;    // The neighboring vertices of container for each point in the current solution 
int *NNeighborsDD;   // The number of neighboring vertices for each point

int **AdjacentB;     // The neighboring edges in the container for each point. For each edge <v1,v2>, the first vertex  v1 is stored, ometting the second vertex. 
int *NNeighborsB;    // Then number of neighboring edges in the container for each point. 

int ***AdjacentDDH;  // The neighboring vertices of ploygon holes for each point
int **NNeighborsDDH; // The numbers of neighboring vertice in each ploygon hole for each point

int ***AdjacentBH;   // The neighboring edges of ploygon hole for each hole and each point
int **NNeighborsBH;  // The numbers of neighboring edges for each ploygon hole ans each point




