#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

// The function MBH() aims to minimize the fucntion E_D(X).  
void MBH(Solution &S, double CR)
{
	
	 int i, step;
     INT MaxNoImprove = 10;// The search depth of MBH. 
     INT NoImprove;
     double delt_f;

     R =  CR;
	 NoImprove = 0 ;
     while(NoImprove < MaxNoImprove)
     {

       for(i=0;i<n;i++) OS.x[i]= S.x[i];
       OS.f = S.f; 
       Shake(OS, n, 0.8*CR); // 0.8*CR is the perturbation strength of the function shake(). 
       LocalSearch(OS);
    
       delt_f = OS.f - S.f ;
       if(delt_f < 0 )
	   {
		  for(i=0;i<n;i++) S.x[i]= OS.x[i];
          S.f = OS.f;
	      NoImprove = 0 ;
	     // printf("************ improved ! f = %12.10e\n", OS.f);
	   }
	   else NoImprove++; 
	   if(S.f < 1.0e-24) break;
     }
     
	 S.CurrentR = R;
	 S.f = myvalue(S.x,n);
}

