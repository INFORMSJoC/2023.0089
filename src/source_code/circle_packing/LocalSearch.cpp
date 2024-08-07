#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

// The multiple-phase local optimization method for the function E_D(X).                                                                                                      
void LocalSearch(Solution &S)
{

	int	I,i;

    cg_descent (S.x, n, NULL, NULL, 1.e-3, myvalue, mygrad, NULL, NULL) ; // Local minimization considering all neighbors for each circle
   
    double dc5 = 4.0*R; 
    Find_Neighbor(S.x,dc5);  // find the neighbors of each circle, i.e., the neighboring vertices and edges of the container and holes. 
    cg_descent (S.x, n, NULL, NULL, 1.e-6, myvalueNN, mygradNN, NULL, NULL) ; // Local minimization only considering neighbors
    S.f = myvalue(S.x,n) ;  
	
    dc5 = 2.5*R;
    Find_Neighbor(S.x,dc5);  // find the neighbors of each circle
    cg_descent (S.x, n, NULL, NULL, 1.e-13, myvalueNN, mygradNN, NULL, NULL) ; // Local minimization only considering neighbors
    cg_descent (S.x, n, NULL, NULL, 1.e-13, myvalue, mygrad, NULL, NULL) ;
    S.f = myvalue(S.x,n) ; 
    
}

// The procedure of adjusting the common radius R of circles. 
double AdjustRadius(Fsolution &S)
{

	double pho = 5.0;
    alpha = 1.0e1;
    int I,i;

    double dc5 = 5.0*S.x[n];
    double RRR = S.x[n]; 
	cg_descent (S.x, n+1, NULL, NULL, 1.e-12, myvalueRR, mygradRR, NULL, NULL) ;  
	for(int k = 0; k < 10; k++)
	{
	   
       alpha = alpha*pho;
       double dc3 = 4.0*S.x[n];
       cg_descent (S.x, n+1, NULL, NULL, 1.e-13, myvalueRR, mygradRR, NULL, NULL) ;
     //  printf("R=%lf\n",S.x[n]); 
	   if(S.x[n] <= 0 )  S.x[n] = RRR; 
	}
	return S.x[n];
}
