#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

// The purpose of the function Area() is to calculate the area of a ploygon whose vertices are stored in a list (denoted by plist). 
double Area(int numpp,double *plist)
{

	double **A; 
	A = new double *[numpp+1];  
	for(int i=0; i< numpp+1; i++) A[i] = new double [2]; 	
	double aarea=0;
	for(int l=0;l<=numpp;l++)
	{
		A[l][0] = plist[2*l];   // x coordinate of the vetex  
		A[l][1] = plist[2*l+1]; // y coordinate of the vetex 
	}
	for(int i=0,j=1; i<numpp; i++,j=(j+1)%numpp) 
	{
	   aarea += A[i][0]*A[j][1]-A[i][1]*A[j][0];
	}
	aarea = fabs(aarea*0.5); 

    for(int i=0;i<=numpp; i++) delete [] A[i] ; 
    delete [] A; 
	return aarea;
}


void TSGO()
{
    int i;
	double ff;
	double area = 0;
    double CR;
	double LargerR;
	double delt_rou = 0.001;
    starting_time = clock();
    
    /* The first stage aims to find an intial radius of circles */
    area =Area(numP,pointlist);
	for(int k=0;k<HN;k++)
	{
		double area_hole = Area(numPH[k],holepoint[k]);
		area -= area_hole;
    }

    double rou = 0.85; // the estimated density of intial configuration 
    CR = sqrt(rou*area/(N*pi) ); // the radius of circles in the initial solution 
    
    while(1)
    {
      initial(NS, n);
      NS.CurrentR = CR; 
      R = CR;
      R_max = CR;
      LocalSearch(NS);

      TabuSearch(NS, 10, R);
      
      for(i=0;i<n;i++) FS.x[i] = TBestS.x[i];
      FS.x[n] = CR;
      AdjustRadius(FS);	
      
      for(int i =0; i<n; i++) CS.x[i] = FS.x[i]; 
      R= FS.x[n]; 
      LocalSearch(CS);
    //  printf("f=%e\n",myvalue(CS.x,n) );
      
      if(FS.x[n] > 0 && myvalue(CS.x,n) < 1.0e-16)
      {
        R_max = FS.x[n];
        for(int i=0; i<n; i++) CBestS.x[i] = CS.x[i];
        CBestS.CurrentR = FS.x[n];
        R = FS.x[n];
        CBestS.f = myvalue(CBestS.x,n);
        final_time = 1.0*(clock()- starting_time)/CLOCKS_PER_SEC ;
        break; 
      }
	}
	
   // printf("intial R =%12.15f\n",R_max) ;
   /*  The second stage aims to find an improving packing configuration */
    while(1.0*(clock()- starting_time)/CLOCKS_PER_SEC < Time_limit)
	{
		LargerR = R_max;
		initial(NS, n);
	    TabuSearch(NS, 50, LargerR); // the search depth of tabu search is set to 50. 
		R = LargerR;
		if(myvalue(TBestS.x,n) < 1.0e-24) 
		{
             for(int i=0;i<n;i++) FS.x[i] = TBestS.x[i];
             FS.x[n] = LargerR;
             AdjustRadius(FS);
             
             for(int i =0; i<n; i++) CS.x[i] = FS.x[i]; 
             R= FS.x[n]; 
             LocalSearch(CS);
           //  printf("f=%e\n",myvalue(CS.x,n) );
             if(FS.x[n] > R_max + 1.0e-10 && myvalue(CS.x,n) < 1.0e-16)
             {
                for(int i=0;i<n;i++) CBestS.x[i] = CS.x[i];
                CBestS.CurrentR = FS.x[n];
                R = CBestS.CurrentR;
                CBestS.f = myvalue(CBestS.x,n);
                R_max =  CBestS.CurrentR;
                final_time = 1.0*(clock()- starting_time)/CLOCKS_PER_SEC ;
               // printf("new R = %12.15f time =%f seconds\n", R_max, final_time) ;
             }
        
	    }

    }

}


