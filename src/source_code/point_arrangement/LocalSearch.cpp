#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"
 
// local optimization method for the function E_D(X).                                                                                                           
void LocalSearch(Solution &S)
{
    cg_descent (S.x, n, NULL, NULL, 1.e-10, myvalue, mygrad, NULL, NULL) ;
    S.f = myvalue(S.x,n) ; 
}

// The procedure of adjusting the minimum distance between points. 
double AdjustRadius(Fsolution &S)
{
    int I,i;
	double pho = 5.0;
    double dc5 = 5.0*S.x[n];
    double RRR = S.x[n]; 
    
    alpha = 1.0e1;
	cg_descent (S.x, n+1, NULL, NULL, 1.e-12, myvalueRR, mygradRR, NULL, NULL) ; 
	for(int k = 0; k < 15; k++)
	{
       alpha = alpha*pho;
       cg_descent (S.x, n+1, NULL, NULL, 1.e-13, myvalueRR, mygradRR, NULL, NULL) ;
	   if(S.x[n] <= 0 )  S.x[n] = RRR; 
	}
	return S.x[n];
}
