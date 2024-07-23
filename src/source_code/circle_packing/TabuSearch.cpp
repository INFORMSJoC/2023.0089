#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

void sort(double *a, int length, int* b)
{
    int i,j, t1;
    double t; 
    for(j=0; j<length; j++)
        for(i=0; i<length-1-j; i++)
            if(a[i]<a[i+1])
            {
                t=a[i];
                a[i]=a[i+1];
                a[i+1]=t;
 
                t1=b[i];
                b[i]=b[i+1];
                b[i+1]=t1;
            }
}
// the purpose of function CircleSiteConstuct() is to find the highest-energy circles in the current configration S. 
void CircleSiteConstuct(Solution &S, int CSites[], int MaxN, double rr)
{
	double time_1 = clock();
	double *f, *f1; 
	int *identi;
	f  = new double [N];
	f1 = new double [N];
	identi = new int [N]; 
	double MaxF = 0;
	int worst_k;
	for(int i=0; i<2*N; i++) CCS.x[i] = S.x[i] ; 
	for(int i=0; i<MaxN; i++) CSites[i] = 0;
	CCS.CurrentR = rr; 
	R = rr;  	
    myvalue_each(CCS.x,n,f);
   
	for(int k = 0; k < N; k++)
	{
		identi[k] = k; f1[k] = f[k];
	}	
	sort(f1,N,identi);
	
    for(int k = 0; k< MaxN; k++)
    {
    	CSites[k] = identi[k];
	}
    
	double time_2 = 1.0*(clock()- time_1)/CLOCKS_PER_SEC ;
	delete [] f1;
	delete [] f;
	delete [] identi; 
}

// the purpose of function VacancySiteConstuct() is to find the lowest-energy vacancy sites in the configuration S. 
void VacancySiteConstuct(Solution &S, Site VSites[], int MaxV, double r)
{
	double time_1 = clock();
	int Flag;         // 0-1 varaible to indicate whether a site is very close to one of existing sites
	int NumberOfSite; // current number of the vacancy lattice sites
	int worst_k;      // the worst site in the set of generated vacancy lattice sites
	double MaxE;      // the maximum energy for the generated vacancy lattice sites
	double dist;      // distance between two sites
	Site V;
	
	/* initialization for the sub-local optimization, used in the functions myvalue_sub and mygrad_sub */ 
	VacancyRadius  = r ; 
	for(int i=0; i<2*N; i++) VS.x[i] = S.x[i] ; 
	VS.CurrentR = S.CurrentR;  
	
	cg_default1(&UParm1) ;
	
	/* Start to detect the vacancy lattice sites */ 
	NumberOfSite = 0; 
	int i = 0;
	while(i < 5*N)
	{
		V.x[0] = X_min + (1.0*(rand()%RAND_MAX)/RAND_MAX)*(X_max - X_min);
	    V.x[1] = Y_min + (1.0*(rand()%RAND_MAX)/RAND_MAX)*(Y_max - Y_min);
	    int jud_inhole = 0;
	    for(int hh=0;hh<HN;hh++)
	    {
	    	int jiaodian = judin(V.x[0],V.x[1],numPH[hh],holepoint[hh]); // to check whether the point lies in the hole hh. 
	    	if(jiaodian%2) 
	    	{
	    		jud_inhole = 1;
	    		break;
			}
		}
		if(jud_inhole) continue;
	    int jud_in = judin(V.x[0],V.x[1],numP,pointlist);
	    if(jud_in%2) i++;
	    else continue;

	   
		cg_descent(V.x, 2, NULL, &UParm1, 1.e-12, myvalue_sub, mygrad_sub, NULL, NULL); // Local minimization
		V.e = myvalue_sub(V.x,2);
	
		Flag = 1; 
		MaxE = 0; 
		for(int k = 0; k < NumberOfSite && k < MaxV ; k++)
		{
			dist = (VSites[k].x[0] - V.x[0])*(VSites[k].x[0] - V.x[0]) + (VSites[k].x[1] - V.x[1])*(VSites[k].x[1] - V.x[1]) ; 
			dist = sqrt(dist);
			if(dist < r*0.5) { Flag = 0; break; }
			if(VSites[k].e  > MaxE)  { MaxE = VSites[k].e; worst_k = k; }  
		}
		
		if(Flag == 1) 
		{
	  	   if(NumberOfSite < MaxV)
	  	   {
	  	   	 VSites[NumberOfSite].x[0] = V.x[0] ; 
			 VSites[NumberOfSite].x[1] = V.x[1] ;
			 VSites[NumberOfSite].e    = V.e ;
			 NumberOfSite ++;
		   }
		   else  
		   { 
		      if(V.e < VSites[worst_k].e )
		      {
		     	VSites[worst_k].x[0] = V.x[0] ; 
			    VSites[worst_k].x[1] = V.x[1] ;
			    VSites[worst_k].e    = V.e ; 
			  }
		   	 
		   }
		   
		}
		
	}
		
}


void TabuSearch(Solution &S, int SearchDepth, double rr)
{
 int iter; 
 int NoImprove; 
 int num_tabu_best, num_best ;  // the number of tabu neighbors and non-tabu neighbors
 double tabu_best_E, best_E, EnergyOfNeigboor ; // the energy for one circle or vacancy site
 int best_v[2]; // best non-neighbor solution denoted by (v[0],v[1]), where v[0] denotes the nubmering of circle to be moved and v[1] denotes the numbering of site to be occupyied. 
 int tabu_best_v[2] ; // best tabu-neighbor solution denoted by (v[0],v[1]), where v[0] denotes the nubmering of circle to be moved and v[1] denotes the numbering of site to be occupyied.
 int CircleSelected; // the index of circle to be moved. 
 double TS_time, TS_starting_time;
 
 TS_starting_time = clock();	 
 R = rr; 
 LocalSearch(S);
 for(int i = 0; i< 2*N; i++) TBestS.x[i] = S.x[i] ; 
 TBestS.f = S.f ; 
 TBestS.CurrentR = S.CurrentR ;
 
 for(int i = 0; i< 2*N; i++) CS.x[i] = S.x[i] ; 
 CS.f = S.f ; 
 CS.CurrentR = S.CurrentR ;
 
 for(int i = 0; i< N; i++) TabuTenure[i] = 0; 
 iter = 0;
 NoImprove = 0; 
 // printf("TS \n");
 while(NoImprove < SearchDepth && iter < 200 && 1.0*(clock()- starting_time)/CLOCKS_PER_SEC < Time_limit)
  {
    //A. To construct the insersion neighborhood N1.
 	tabu_best_E = 9999999.0 ; 
    best_E = 9999999.0 ;
    num_tabu_best = 0 ; 
    num_best = 0 ;
   	VacancySiteConstuct(CS, VacancySite, MaxVS, rr); 
   	CircleSiteConstuct(CS, CircleSite, MaxNC, rr);  
   	for(int I= 0; I < MaxNC; I++)
   	{
   		for(int J = 0; J < MaxVS; J++)
   		{
   		    for(int i = 0; i< 2*N; i++) NS.x[i] = CS.x[i] ; 
            NS.f = CS.f ; 
            NS.CurrentR                =  CS.CurrentR ;	
            CircleSelected             =  CircleSite[I]; 
            NS.x[2*CircleSelected]     =  VacancySite[J].x[0]; 
            NS.x[2*CircleSelected + 1] =  VacancySite[J].x[1];  
            LocalSearch(NS);
            EnergyOfNeigboor = myvalue(NS.x,n); 
          //  printf("f=%e\n",EnergyOfNeigboor);
            if( iter >= TabuTenure[ CircleSelected ]  ) //non-tabu 
             {
             	if(EnergyOfNeigboor < best_E) 
				 {
				 	best_E    = EnergyOfNeigboor;
					best_v[0] = CircleSelected; 
				 	best_v[1] = J; 
				 
				 }
				 
				 num_best ++; 
			 }
			 else  //tabu
			 {
			 	if(EnergyOfNeigboor < tabu_best_E) 
				 {
				 	tabu_best_E    = EnergyOfNeigboor; 
					tabu_best_v[0] = CircleSelected; 
				 	tabu_best_v[1] = J; 
				 
				 }
				 num_tabu_best ++; 
			 }
		}
   		
	}
   	
   	 // B. To perform the neighborhood move. 
   	 if( ( num_tabu_best > 0 && tabu_best_E < best_E && tabu_best_E < TBestS.f ) || num_best == 0 )
 	 {
 		CS.x[2*tabu_best_v[0]  ]     =  VacancySite[tabu_best_v[1]].x[0];
 	    CS.x[2*tabu_best_v[0]+1]     =  VacancySite[tabu_best_v[1]].x[1]; 
 	    TabuTenure[ tabu_best_v[0] ] += (iter + 5 + rand()%5 );  
 	  //  printf("\n tabu %d    %d    %d\n", tabu_best_v[0], num_best, num_tabu_best); 
	 }
	 else
	 {
	 	CS.x[2*best_v[0]  ]     =  VacancySite[best_v[1]].x[0];
 	    CS.x[2*best_v[0]+1]     =  VacancySite[best_v[1]].x[1]; 
        TabuTenure[ best_v[0] ] += (iter + 5 + rand()%5 ); 
       // printf("\n non tabu %d   %d   %d\n",best_v[0], num_best,num_tabu_best);  
	 }
	 
	//C. Save the improved sokutions 
	LocalSearch(CS);// printf("current S =%e\n", CS.f);
	MBH(CS, rr); // perform the MBH procedure
//	cg_descent (CS.x, n, NULL, NULL, 1.e-12, myvalue, mygrad, NULL, NULL) ; // Local minimization
//	printf("true f =%e\n",myvalue(CS.x, n)); 
 	if(myvalue(CS.x,n) < TBestS.f)
 	{
 	   for(int i = 0; i< 2*N; i++) TBestS.x[i] = CS.x[i] ; 
       TBestS.f = CS.f ; 
       TBestS.CurrentR = CS.CurrentR ;	
	   NoImprove = 0 ; 
	   TS_time = (double) (1.0*(clock()-TS_starting_time)/CLOCKS_PER_SEC); 	
	  // printf("\n TS iteration :  %8d       %e       %e    %lf \n", iter, CS.f, TBestS.f, TS_time); 
	}
	else NoImprove ++; 
	
	if(TBestS.f < 1.0e-24) break; 
 	iter++; 
 	
  }

}



