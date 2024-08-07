#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

void AssignMemery()
{

	 int i;
	 CS.x      = (double *) malloc (n*sizeof (double)) ; CS.f = 0;
	 NS.x      = (double *) malloc (n*sizeof (double)) ; NS.f = 0;
	 VS.x      = (double *) malloc (n*sizeof (double)) ; VS.f = 0;
	 CCS.x     = (double *) malloc (n*sizeof (double)) ; CCS.f = 0;
     OS.x      = (double *) malloc (n*sizeof (double)) ; OS.f = 0;
     TBestS.x  = (double *) malloc (n*sizeof (double)) ; TBestS.f = 0; 
     CBestS.x  = (double *) malloc (n*sizeof (double)) ; CBestS.f = 0;
     GlobalS.x = (double *) malloc (n*sizeof (double)) ; GlobalS.f = 0;
     FS.x      = (double *) malloc ((n+1)*sizeof (double));   FS.L = 0;
     
     Adjacent   =  new int *[N];
     AdjacentDD =  new int *[N];
     AdjacentB  =  new int *[N];
     AdjacentDDH=  new int **[N];
     AdjacentBH =  new int **[N];
     TabuTenure =  new int [N] ; 
	 for(int i=0; i<N; i++)
	 {
	 	Adjacent[i]    = new int [MaxNN];
	 	AdjacentDD[i]  = new int [MaxNN];
	 	AdjacentB[i]   = new int [MaxNN];
	 	AdjacentDDH[i] = new int *[HN];
	 	AdjacentBH[i]  = new int *[HN]; 
	 } 
	 for(int i=0; i<N; i++)
	 {
	 	for(int j=0; j<HN; j++)
		{
		 	AdjacentDDH[i][j] = new int [MaxNN];
		    AdjacentBH[i][j]  = new int [MaxNN];
	    }
	 }
	 
	 NNeighbors   = new int [N];
	 NNeighborsDD = new int [N];
	 NNeighborsB  = new int [N];
	 NNeighborsDDH= new int *[N];
	 NNeighborsBH = new int *[N];
	 
	 for(int i=0; i<N; i++)
	 {
	 	NNeighborsDDH[i] = new int [MaxNN];
	    NNeighborsBH[i]  = new int [MaxNN];
     }
     
    
 }
 
// The function Inputing() aims to read the instance (i.e., the region to be packed).
int Inputing(char *filename) 
{
	char buff[80];
	FILE *fpRead; 
    sprintf(buff,"%s",filename);
    fpRead = fopen(buff,"r");
    

    if(fpRead==NULL)
    {
    	return 0;
	}
	fscanf(fpRead,"%d", &numP);  

	double pt1,pt2;
	for(int l=0;l<=numP;l++)
	{
		if(!feof(fpRead))
		{
			fscanf(fpRead,"%lf %lf\n",&pt1,&pt2);
			pointlist[2*l]  = pt1;   
			pointlist[2*l+1]= pt2;    
		}
	}
	
	Range();
	

	char hole; 
	int holnum = 0; // the numbering or index of hole
	HN = 0;
	fscanf(fpRead,"%c %d",&hole,&numPH[holnum]); 
//	printf("%c  %d \n", hole, numPH[holnum]);  
	
	if(hole == 'H') 
	{
	
		while(1)
		{
			for(int l=0;l<=numPH[holnum];l++)
			{
				if(!feof(fpRead))
				{
					fscanf(fpRead,"%lf %lf\n",&pt1,&pt2);
					holepoint[holnum][2*l]   = pt1;   //x coordinate of vertex of hole 
					holepoint[holnum][2*l+1] = pt2;   //y coordinate of vertex of hole
				}
	        }
	        
	        holnum++;
	        fscanf(fpRead,"%c %d",&hole,&numPH[holnum]);// printf("%c  %d \n", hole, numPH[holnum]);  
	        if(hole == 'H') ; 
		  	else break; 
		}
		HN = holnum;
	}
	
	else
	{
		//printf("no hole !\n");
	}
		
	fclose(fpRead);
	return 1; 
}

double Deviation(double arr[],int n)
{
    int i;
    double sum = 0,tmp = 0, x_avg;
    for(i = 0; i < n; ++i) sum += arr[i];
    x_avg = sum / n;
    for(i = 0; i < n; ++i)  tmp += (arr[i] - x_avg)*(arr[i] - x_avg);
    return sqrt(tmp/n);
} // standard deviation 

void Outputing(Solution &S, int n)
{
    int i;
    double x,y,z;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%d.txt",n/2);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d    %.10f\n", n/2, S.CurrentR);
    for(i=0;i<n/2;i++)
    {
	   int I = 2*i;
	   fprintf(fp,"%15.15f     %15.15f\n", S.x[I],S.x[I+1]);
	}
	fclose(fp);
}

void Outputing1(Solution &S, char *filename, int n)
{
    int i;
    double x,y,z;
	FILE *fp;
	char buff[80];
    sprintf(buff,"Sol_%d_%s", n/2, filename);
    fp=fopen(buff,"a+");
    fprintf(fp,"%4d     %18.12f\n", n/2, S.CurrentR);
  
	for(i=0;i<N;i++) 
    {
	   int I = 2*i; 
	   fprintf(fp,"%4d", i+1); 
	   fprintf(fp,"     %18.12f",S.x[I]); 
	   fprintf(fp,"     %18.12f\n",S.x[I+1]); 
	  
	} 

	fclose(fp);
}

void Outputing2(Fsolution &S, int n)
{
    int i;
    double x,y,z;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%d.txt",111);
    fp=fopen(buff,"a+");
    fprintf(fp,"%4d     %18.12f\n", n/2, S.x[n]);
 
	for(i=0;i<N;i++) 
    {
	   int I = 2*i; 
	   fprintf(fp,"%4d", i+1); 
	   fprintf(fp,"     %18.12f",S.x[I]); 
	   fprintf(fp,"     %18.12f\n",S.x[I+1]); 
	  
	} 

	fclose(fp);
}

void Out_results(double best , double ave,  double worst, int sr, double AvgTime, double deviation, char *filename, char *inputfile, int N)
{
    int i;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%s",filename);
    fp = fopen(buff,"a+");
    fprintf(fp,"%s    %d   %.10f   %.10f   %.10f  %d  %e   %.2f\n",inputfile, N, best, ave, worst, SR, deviation, AvgTime);
	fclose(fp);
}
 
 
 
