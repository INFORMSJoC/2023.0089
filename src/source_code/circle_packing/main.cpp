#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

int main (int argc, char *argv[])
{
    int i ;
	int runs;
	double Results[100];
	double f_best, f_avg, f_worst, sigma, AvgTime ;
    srand(time(NULL)); 

    N = atoi(argv[1]); 
    Input_File_Name = argv[2]; 
    runs = atoi(argv[3]);
    Time_limit = atof(argv[4]);
    
    Inputing(Input_File_Name); 
    
    Output_File_Name = "ComputationalResults_Cirlces.txt";                                                                        
	n = 2*N;   // number of variables
	AssignMemery(); // allocate space for solutions
	
	f_best = 0.0;
	f_avg  = 0.0;
	f_worst= 9999999;
	AvgTime= 0.0;
	SR = 0;
    //TSGO is performed 'runs' times.
	for(i = 0; i < runs; i ++)
	{
        TSGO(); 
        Results[i] = CBestS.CurrentR;
        f_avg += CBestS.CurrentR;
        if(CBestS.CurrentR > f_best + 1.0e-8) 
		{
		   SR = 1;                                       
		   f_best = CBestS.CurrentR;
		   for(int j=0;j<n;j++)  GlobalS.x[j] = CBestS.x[j];
	       GlobalS.CurrentR = CBestS.CurrentR;
        }
		else if(fabs(f_best - CBestS.CurrentR) < 1.0e-8) 
		{
			SR ++ ;                
		}
        if(CBestS.CurrentR < f_worst) f_worst = CBestS.CurrentR;                                          
        AvgTime += final_time;
	}

	f_avg /= runs;
	AvgTime /= runs;
	sigma = Deviation(Results,runs); 
	
    R= GlobalS.CurrentR;
    cg_descent (GlobalS.x, n, NULL, NULL, 1.e-13, myvalue, mygrad, NULL, NULL) ; 
    GlobalS.f = myvalue(GlobalS.x,n); 
    Outputing1(GlobalS,Input_File_Name, n);
    Out_results(f_best, f_avg,  f_worst, SR, AvgTime, sigma, Output_File_Name,Input_File_Name, N) ;
    return 0;
}



