#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y

// The purpose of function Find_Neighbor() is find the neighboring vertices and edges on the boundaries of the container and holes for each circle. 
void Find_Neighbor(double x[], double Dcut)
{
	
    int I, J, i, j;
	double x_i, y_i, x_j, y_j, dist,a,b,c;
    for(I = 0; I < N; I ++) 
    {
    	NNeighbors[I] = 0;
    	NNeighborsDD[I] = 0;
		NNeighborsB[I] = 0;
	    for(int hh=0;hh<HN;hh++)
	    {
			NNeighborsDDH[I][hh] = 0;
			NNeighborsBH[I][hh] = 0;
	    }
	}
	for(i=0; i<N; i++)
	{
	   for(j=0; j<MaxNN; j++) 
	   {
	   	    Adjacent[i][j] = -1;
	   	    AdjacentDD[i][j] = -1;
	   	    AdjacentB[i][j] = -1;
		    for(int hh=0;hh<HN;hh++)
		    {
				AdjacentDDH[i][hh][j] = -1;
				AdjacentBH[i][hh][j] = -1;
		    }
	   }
    }   
    
    for(I = 0; I < N; I ++)
    {
	  i = 2*I;
	  x_i = x[i];
	  y_i = x[i+1];

	  for(J = I+1; J < N; J ++)
	  {

		j   = 2*J;
		x_j = x[j];
		y_j = x[j+1];

		dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j));
		if(dist < Dcut)
		{
			Adjacent[I][NNeighbors[I]] = J;
			NNeighbors[I]++; 
		}
	  }
	  for(int k = 0; k < numP; k++)
	  {
	  	int K = 2*k;
	  	dist = sqrt( (x_i-pointlist[K])*(x_i-pointlist[K]) + (y_i - pointlist[K+1])*(y_i - pointlist[K+1]) );
	  	if(dist < Dcut)
		{
			AdjacentDD[I][NNeighborsDD[I]] = k; 
			NNeighborsDD[I]++; 
		}
	  } 
	  for(int P = 0; P < numP; P++)
	  {
		int p=2*P;
		{
			if(pointlist[p+2] == pointlist[p])
			{
				a = 1.0;
				b = 0;
				c = -pointlist[p];
			}
			else
			{
				a = (pointlist[p+3] - pointlist[p+1]) / (pointlist[p+2] - pointlist[p]);
				b = -1.0;
		        c = pointlist[p+1] - a*pointlist[p];
			}
		    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
		    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);

		    double maxx = MAX(pointlist[2+p],pointlist[p]);
		    double minx = MIN(pointlist[2+p],pointlist[p]);
		    if(cdx <= maxx && cdx >= minx)
		    {
			    dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
			    if(dist < Dcut)
				{
					AdjacentB[I][NNeighborsB[I]] = P; 
					NNeighborsB[I]++;
				}
			    else
			    {
			    	double dist1 = sqrt( (x_i - pointlist[p]) * (x_i - pointlist[p]) + (y_i - pointlist[p+1]) * (y_i - pointlist[p+1]) );
			    	double dist2 = sqrt( (x_i - pointlist[p+2]) * (x_i - pointlist[p+2]) + (y_i - pointlist[p+3]) * (y_i - pointlist[p+3]) );
			    	if(dist1 < Dcut || dist2 < Dcut)
					{
						AdjacentB[I][NNeighborsB[I]] = P; 
						NNeighborsB[I]++;
					}
				}
		    } 
		}
	  }	  
      for(int hh=0;hh<HN;hh++)
	  {
	  	for(int k=0;k<numPH[hh];k++)
	    {
		  	int K = 2*k;
		  	dist = sqrt( (x_i-holepoint[hh][K])*(x_i-holepoint[hh][K]) + (y_i - holepoint[hh][K+1])*(y_i - holepoint[hh][K+1]) );
		  	if(dist < Dcut)
			{
				AdjacentDDH[I][hh][NNeighborsDDH[I][hh]] = k; 
				NNeighborsDDH[I][hh]++; 
			}
	    }
	  }
      for(int hh=0;hh<HN;hh++)
	  {
	  	for(int k=0;k<numPH[hh];k++)
	    {
	    	int K = 2*k; 
			{
				if(holepoint[hh][K+2] == holepoint[hh][K])
				{
					a = 1.0;
					b = 0;
					c = -holepoint[hh][K];
				}
				else
				{
					a = (holepoint[hh][K+3] - holepoint[hh][K+1]) / (holepoint[hh][K+2] - holepoint[hh][K]);
					b = -1.0;
			        c = holepoint[hh][K+1] - a*holepoint[hh][K];
				}
			    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
			    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);

			    double maxx = MAX(holepoint[hh][2+K],holepoint[hh][K]);
			    double minx = MIN(holepoint[hh][2+K],holepoint[hh][K]);
			    double maxy = MAX(holepoint[hh][2+K+1],holepoint[hh][K+1]);
	    	    double miny = MIN(holepoint[hh][2+K+1],holepoint[hh][K+1]);
			    if(cdx <= maxx && cdx >= minx)
			    {
				   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
				   if(dist < Dcut)
				   {
						AdjacentBH[I][hh][NNeighborsBH[I][hh]] = k; 
						NNeighborsBH[I][hh]++; 
				   }

			    }
			    else
			    {
			    	double dist1 = sqrt( (x_i - holepoint[hh][K]) * (x_i - holepoint[hh][K]) + (y_i - holepoint[hh][K+1]) * (y_i - holepoint[hh][K+1]) );
			    	double dist2 = sqrt( (x_i - holepoint[hh][K+2]) * (x_i - holepoint[hh][K+2]) + (y_i - holepoint[hh][K+3]) * (y_i - holepoint[hh][K+3]) );
			    	if(dist1 < Dcut || dist2 < Dcut)
					{
						AdjacentBH[I][hh][NNeighborsBH[I][hh]] = k; 
						NNeighborsBH[I][hh]++;
					}
				}
			}
	    }
	  }
	}

}

void Find_Neighbor2(double x[], int num, double y[], double Dcut)
{
	
    int I, J, i, j;
	double x_i, y_i, x_j, y_j, dist,a,b,c;
    for(I = 0; I < num; I ++) 
    {
    	NNeighbors[I] = 0;
    	NNeighborsDD[I] = 0;
		NNeighborsB[I] = 0;
	    for(int hh=0;hh<HN;hh++)
	    {
			NNeighborsDDH[I][hh] = 0;
			NNeighborsBH[I][hh] = 0;
	    }
	}
	for(i=0; i<num; i++)
	{
	   for(j=0; j<MaxNN; j++) 
	   {
	   	    Adjacent[i][j] = -1;
	   	    AdjacentDD[i][j] = -1;
	   	    AdjacentB[i][j] = -1;
		    for(int hh=0;hh<HN;hh++)
		    {
				AdjacentDDH[i][hh][j] = -1;
				AdjacentBH[i][hh][j] = -1;
		    }
	   }
    }   
    
    for(I = 0; I < num; I ++)
    {
	  i = 2*I;
	  x_i = x[i];
	  y_i = x[i+1];

	  for(J = 0; J < N; J ++)
	  {

		j   = 2*J;
		x_j = y[j];
		y_j = y[j+1];

		dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j));
		if(dist < Dcut)
		{
			Adjacent[I][NNeighbors[I]] = J;
			NNeighbors[I]++;
		}
	  }
	  for(int k = 0; k < numP; k++)
	  {
	  	int K = 2*k;
	  	dist = sqrt( (x_i-pointlist[K])*(x_i-pointlist[K]) + (y_i - pointlist[K+1])*(y_i - pointlist[K+1]) );
	  	if(dist < Dcut)
		{
			AdjacentDD[I][NNeighborsDD[I]] = k; 
			NNeighborsDD[I]++; 
		}
	  } 
	  for(int P = 0; P < numP; P++)
	  {
		int p=2*P;
		{
			if(pointlist[p+2] == pointlist[p])
			{
				a = 1.0;
				b = 0;
				c = -pointlist[p];
			}
			else
			{
				a = (pointlist[p+3] - pointlist[p+1]) / (pointlist[p+2] - pointlist[p]);
				b = -1.0;
		        c = pointlist[p+1] - a*pointlist[p];
			}
		    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
		    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);  

		    double maxx = MAX(pointlist[2+p],pointlist[p]);
		    double minx = MIN(pointlist[2+p],pointlist[p]);
		    if(cdx <= maxx && cdx >= minx)
		    {
			    dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
			    if(dist < Dcut)
				{
					AdjacentB[I][NNeighborsB[I]] = P; 
					NNeighborsB[I]++;
				}
			    else
			    {
			    	double dist1 = sqrt( (x_i - pointlist[p]) * (x_i - pointlist[p]) + (y_i - pointlist[p+1]) * (y_i - pointlist[p+1]) );
			    	double dist2 = sqrt( (x_i - pointlist[p+2]) * (x_i - pointlist[p+2]) + (y_i - pointlist[p+3]) * (y_i - pointlist[p+3]) );
			    	if(dist1 < Dcut || dist2 < Dcut)
					{
						AdjacentB[I][NNeighborsB[I]] = P; 
						NNeighborsB[I]++; 
					}
				}
		    } 
		}
	  }	  
      for(int hh=0;hh<HN;hh++)
	  {
	  	for(int k=0;k<numPH[hh];k++)
	    {
		  	int K = 2*k;
		  	dist = sqrt( (x_i-holepoint[hh][K])*(x_i-holepoint[hh][K]) + (y_i - holepoint[hh][K+1])*(y_i - holepoint[hh][K+1]) );
		  	if(dist < Dcut)
			{
				AdjacentDDH[I][hh][NNeighborsDDH[I][hh]] = k;  
				NNeighborsDDH[I][hh]++; 
			}
	    }
	  }
      for(int hh=0;hh<HN;hh++)
	  {
	  	for(int k=0;k<numPH[hh];k++)
	    {
	    	int K = 2*k; 
			{
				if(holepoint[hh][K+2] == holepoint[hh][K])
				{
					a = 1.0;
					b = 0;
					c = -holepoint[hh][K];
				}
				else
				{
					a = (holepoint[hh][K+3] - holepoint[hh][K+1]) / (holepoint[hh][K+2] - holepoint[hh][K]);
					b = -1.0;
			        c = holepoint[hh][K+1] - a*holepoint[hh][K];
				}
			    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
			    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);

			    double maxx = MAX(holepoint[hh][2+K],holepoint[hh][K]);
			    double minx = MIN(holepoint[hh][2+K],holepoint[hh][K]);
			    double maxy = MAX(holepoint[hh][2+K+1],holepoint[hh][K+1]);
	    	    double miny = MIN(holepoint[hh][2+K+1],holepoint[hh][K+1]);
			    if(cdx <= maxx && cdx >= minx)
			    {
				   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
				   if(dist < Dcut)
				   {
						AdjacentBH[I][hh][NNeighborsBH[I][hh]] = k; 
						NNeighborsBH[I][hh]++; 
				   }

			    }
			    else
			    {
			    	double dist1 = sqrt( (x_i - holepoint[hh][K]) * (x_i - holepoint[hh][K]) + (y_i - holepoint[hh][K+1]) * (y_i - holepoint[hh][K+1]) );
			    	double dist2 = sqrt( (x_i - holepoint[hh][K+2]) * (x_i - holepoint[hh][K+2]) + (y_i - holepoint[hh][K+3]) * (y_i - holepoint[hh][K+3]) );
			    	if(dist1 < Dcut || dist2 < Dcut)
					{
						AdjacentBH[I][hh][NNeighborsBH[I][hh]] = k; 
						NNeighborsBH[I][hh]++;
					}
				}
			}
	    }
	  }
	}

}

