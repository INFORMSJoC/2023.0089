#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y

// The function Range() aims to find the minimum rectangle containing the region to be packed. 
void Range()
{
	 X_min  =  10000;
	 X_max  = -10000;
	 Y_min  =  10000;
	 Y_max  = -10000; 
	 
     for(int J = 0; J < numP;J++)
     {
    	int j = 2*J;
    	if(pointlist[j] < X_min) X_min = pointlist[j]; 
    	if(pointlist[j] > X_max) X_max = pointlist[j];  
    	
    	if(pointlist[j+1] < Y_min) Y_min = pointlist[j+1]; 
    	if(pointlist[j]+1 > Y_max) Y_max = pointlist[j+1];  
	 } 
	
}
// The function initial() aims to generate randomly an inital solution S in a rectangle. 
void initial(Solution &S, INT n)   
{

	INT I,i;
	I=0;
	while(I < n/2)
	{
		i = 2*I;
	    double x_i = X_min + (1.0*(rand()%RAND_MAX)/RAND_MAX)*(X_max - X_min);
	    double y_i = Y_min + (1.0*(rand()%RAND_MAX)/RAND_MAX)*(Y_max - Y_min);
		S.x[i]   = x_i;  
		S.x[i+1] = y_i;
		I++;
	}
	S.f = 9999999.0; 
}

