#include <math.h>
#include "cg_user.h"
#include "Common.h"

#define MIN(x,y) x<y?x:y
#define MAX(x,y) x>y?x:y

const int beishu = 1; // this parameter corresponds to the parameter 'alpha' of Eq. (6) in the paper
const int kb = 2; // this parameter corresponds to the parameter 'gama' of Eqs. (8) and (9) in the paper

// These two parameters are used in the minimization of potential energy of additional detecting circle.
const int beishu1 = 3; 
const int kb1 = 3;

// The purpose of function judin() is to identify whether the point (x_i,y_i) locates in the polygon denoted by 'plitst'. 
// The answer is No if the returned value is an even number, and Yes otherwise. 
int judin(double x_i,double y_i,int num,double *plist)
{
	int jds=0;
	for(int j=0;j<num;j++) 
	{
		int J=2*j;
		double maxx = MAX(plist[2+J],plist[J]);
	    double minx = MIN(plist[2+J],plist[J]);
	    double maxy = MAX(plist[2+J+1],plist[J+1]);
	    double miny = MIN(plist[2+J+1],plist[J+1]);
	    if(maxy != miny)   
		{
	    	if( x_i <= minx ) 
	    	{
	    	
	    		if( y_i <= maxy && y_i > miny )  jds++;    
	        }
	        else if(minx < x_i && x_i < maxx)   
	        {
	        	
	        	if( y_i <= maxy && y_i > miny )
	        	{
	        		double kk = (plist[3+J] - plist[1+J]) / (plist[2+J] - plist[J]);   
        	        double xx = plist[J] + (y_i-plist[J+1])/kk;  
        	        if(xx > x_i) jds++;
				}
			}
		}
		
	}
	return jds;
}

// The purpose of function myvalue() is to calculate the value of the objective function E_D(X) 
double myvalue
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, i, j;
    f = 0.0 ;  
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist);   
		if(jds%2==0)
		{
			double xyp[2] = {pointlist[0],pointlist[1]};  
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));

			for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud = (pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 && x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny ) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 )  
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						xyvalue = newvalue;
					}
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + R; 
				f += kb*(overlap*overlap);
			}
		}
		
		if(jds%2!=0 || yzbs==1) 
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
				if(fabs(x_i - pointlist[J]) > R || fabs(y_i - pointlist[J+1]) > R )  
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = R - dist;
					if(overlap > 0.0) f += (overlap*overlap);  
				}
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					    dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					    overlap = R - dist;
					    if(overlap > 0.0) f += (overlap*overlap); 
				    } 
				}
		
			}
        }
        
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			
			if(jiaodian%2!=0) 
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]};  
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny)  
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
						} 
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b); 
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);

					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								xyvalueh = dist;
							}
					    } 
					}
				}
				if(yzbs==0) 
				{
					overlap = xyvalueh + R; 
					f += kb*(overlap*overlap);
				} 
			}
			
			if(jiaodian%2==0 || yzbs==1) 
			{		
				for(j=0;j<numPH[hh];j++)  
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > R || fabs(y_i - holepoint[hh][J+1]) > R ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = R - dist;
						if(overlap > 0.0) f += (overlap*overlap);  
					} 
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = R - dist;
						   if(overlap > 0.0) f += (overlap*overlap); 
					    } 
					}
			
				}								
			} 
				
	    }
	} 
    f*=beishu;
        
	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
		for (J = I+1; J < n/2; J++)
		{
			j = 2*J;
			x_j = x[j];
			y_j = x[j+1];

			if(fabs(x_i - x_j) > 2.0 * R || fabs(y_i - y_j) > 2.0 * R ) continue;   
		
			dist = sqrt( (x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * R - dist ;
			if(overlap > 0.0) f += (overlap*overlap); 
		}

	}
	
    return (f) ;
}

// The purpose of function mygrad() is to calculate the gradient of the fucntion E_D(X)
void mygrad
(
    double   *g,
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c,aa,bb;
    INT I, J, i, j;
    for (i = 0; i < n; i++) g [i] = 0;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist);  
		if(jds%2==0) 
		{
			double xyp[2] = {pointlist[0],pointlist[1]};  
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));
			double cdgx = 0; 
            for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						cdgx = 0;  
						xyvalue = newvalue;
					} 
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b); 
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b); 

				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							aa = a;
							bb = b;
							cdgx = 1;
							xyvalue = dist;
						} 
				    } 
				}
			}
			if(yzbs==0) 
			{
				overlap = xyvalue + R;  
				if(cdgx)  
				{
					g[i]   +=  kb*beishu*2.0*overlap*( (x_i - xyp[0])*aa*aa + (y_i - xyp[1])*aa*bb )/( xyvalue*(aa*aa+bb*bb) ); 
                    g[i+1] +=  kb*beishu*2.0*overlap*( (x_i - xyp[0])*aa*bb + (y_i - xyp[1])*bb*bb )/( xyvalue*(aa*aa+bb*bb) );
				}
				else  
				{
					g[i]   +=  kb*beishu*2.0*overlap*(x_i - xyp[0])/xyvalue; 
                    g[i+1] +=  kb*beishu*2.0*overlap*(y_i - xyp[1])/xyvalue;					
				}

			}
		}
		if(jds%2!=0 || yzbs==1)
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	
				if(fabs(x_i - pointlist[J]) > R || fabs(y_i - pointlist[J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = R - dist;
	                if(overlap > 0.0)
	                {
	                    g[i]   -=  beishu*2.0*overlap*(x_i - pointlist[J])/dist;
	                    g[i+1] -=  beishu*2.0*overlap*(y_i - pointlist[J+1])/dist;				
	
	                } 
 
				}	
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b); 
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b); 
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0)
					   {
                            g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) ); 
                            g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );
					   }

				    } 
				}
		
			}
        }
        
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			if(jiaodian%2!=0) 
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]};  
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
				double cdgx = 0;
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny)  
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
							cdgx = 0;
						} 
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								aa = a;
								bb = b;
								xyvalueh = dist;
								cdgx = 1;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + R; 
					if(cdgx)
					{
						g[i]   +=  kb*beishu*2.0*overlap*( (x_i - xyph[0])*aa*aa + (y_i - xyph[1])*aa*bb )/( xyvalueh*(aa*aa+bb*bb) );
			            g[i+1] +=  kb*beishu*2.0*overlap*( (x_i - xyph[0])*aa*bb + (y_i - xyph[1])*bb*bb )/( xyvalueh*(aa*aa+bb*bb) );
					}
					else
					{
						g[i]   +=  kb*beishu*2.0*overlap*(x_i - xyph[0])/xyvalueh;
			            g[i+1] +=  kb*beishu*2.0*overlap*(y_i - xyph[1])/xyvalueh;
					}


				}
			}
			if(jiaodian%2==0 || yzbs==1)  
			{		
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > R || fabs(y_i - holepoint[hh][J+1]) > R ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = R - dist;
						if(overlap > 0.0) 
                        {
                            g[i]   -=  beishu*2.0*overlap*(x_i - holepoint[hh][J])/dist;
                            g[i+1] -=  beishu*2.0*overlap*(y_i - holepoint[hh][J+1])/dist;
                        } 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = R - dist;
						   if(overlap > 0.0)  
					        {
			          	        g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) );
			                    g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );	
					        } 
					    } 
					}

				}								
			}
				
	    }
	}

	for (I = 0; I < n/2; I++)
	{
 	    i   = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		for (J = I+1; J < n/2; J++)
		{
   	        j   = 2*J;
			x_j = x[j];
			y_j = x[j+1];

            if(fabs(x_i - x_j) > 2.0*R || fabs(y_i - y_j) > 2.0*R ) continue; 

	 	    dist = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
	 	    overlap = 2.0*R - dist;
	 	    if(overlap > 0 )
	 	    {
			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
			}

		}
	}  
	     
    return ;
}


// The purpose of function myvalue_sub() is to calculate the potential energy E(c_i) for an additional detecting circle c_i, based on a given configration X. 
double myvalue_sub
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, j;
    const int i = 0; 
    f = 0.0 ;  
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
  
	x_i = x[i];
	y_i = x[i+1];
	{
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist); 
		if(jds%2==0)
		{
			double xyp[2] = {pointlist[0],pointlist[1]}; 
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));

			for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud = (pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						xyvalue = newvalue;
					} 
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + VacancyRadius;  
				f += kb1*(overlap*overlap);
			}
		}

		if(jds%2!=0 || yzbs==1) 
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	        
				if(fabs(x_i - pointlist[J]) > VacancyRadius || fabs(y_i - pointlist[J+1]) > VacancyRadius ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = VacancyRadius - dist;
					if(overlap > 0.0) f += (overlap*overlap); 
				}
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = VacancyRadius - dist;
					   if(overlap > 0.0) f += (overlap*overlap); 
				    } 
				}
		
			}
        }
        
      
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			
			if(jiaodian%2!=0)
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]}; 
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
						}
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								xyvalueh = dist;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + VacancyRadius; 
					f += kb1*(overlap*overlap);
				}
			}
			
			if(jiaodian%2==0 || yzbs==1)  
			{		
				for(j=0;j<numPH[hh];j++)  
				{
				    J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > VacancyRadius || fabs(y_i - holepoint[hh][J+1]) > VacancyRadius ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = VacancyRadius - dist;
						if(overlap > 0.0) f += (overlap*overlap); 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = VacancyRadius - dist;
						   if(overlap > 0.0) f += (overlap*overlap); 
					    } 
					}
			
				}								
			}
				
	    }
	}
    f*=beishu1;
	{	
		for (int J = 0; J < N; J++)
		{
			j = 2*J;
			x_j = VS.x[j];
			y_j = VS.x[j+1];

			if(fabs(x_i - x_j) > 2.0 * VacancyRadius|| fabs(y_i - y_j) > 2.0 * VacancyRadius ) continue;  
		
			dist = sqrt( (x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * VacancyRadius - dist ;
			if(overlap > 0.0) f += (overlap*overlap); 
		}
	}
    return (f) ;
}

// The purpose of function mygrad_sub() is to calculate the gradient of function E(c_i) for an additional detecting circle c_i, based on a given configration X.
void mygrad_sub
(
    double   *g,
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c,aa,bb;
    INT I, J, j;
    const int i = 0; 
    for(int i1 =0; i1<n;i1++) g[i1] = 0; 

	x_i = x[i];
	y_i = x[i+1];
	{

		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist); 
		if(jds%2==0)
		{
			double xyp[2] = {pointlist[0],pointlist[1]}; 
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));
			double cdgx = 0;
            for(j=0;j<numP;j++) 
			{
				J=2*j;
				double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny)  
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						cdgx = 0;  
						xyvalue = newvalue;
					}
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							aa = a;
							bb = b;
							cdgx = 1;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + VacancyRadius; 
				if(cdgx)  
				{
					g[i]   +=  kb1*beishu1*2.0*overlap*( (x_i - xyp[0])*aa*aa + (y_i - xyp[1])*aa*bb )/( xyvalue*(aa*aa+bb*bb) ); 
                    g[i+1] +=  kb1*beishu1*2.0*overlap*( (x_i - xyp[0])*aa*bb + (y_i - xyp[1])*bb*bb )/( xyvalue*(aa*aa+bb*bb) );
                  
				}
				else 
				{
					g[i]   +=  kb1*beishu1*2.0*overlap*(x_i - xyp[0])/xyvalue; 
                    g[i+1] +=  kb1*beishu1*2.0*overlap*(y_i - xyp[1])/xyvalue;
                    					
				}

			}
		}
		if(jds%2!=0 || yzbs==1)
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	
				if(fabs(x_i - pointlist[J]) > VacancyRadius || fabs(y_i - pointlist[J+1]) > VacancyRadius ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = VacancyRadius - dist;
	                if(overlap > 0.0)
	                {
	                    g[i]   -=  beishu1*2.0*overlap*(x_i - pointlist[J])/dist;
	                    g[i+1] -=  beishu1*2.0*overlap*(y_i - pointlist[J+1])/dist;				
	
	                } 
				}	
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0)
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = VacancyRadius - dist;
					   if(overlap > 0.0)
					   {
                            g[i]   -=  beishu1*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) ); 
                            g[i+1] -=  beishu1*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );
                    
					   }
				    } 
				}
		
			}
        }
        
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			if(jiaodian%2!=0)
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]}; 
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
				double cdgx = 0;
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
							cdgx = 0;
						}
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								aa = a;
								bb = b;
								xyvalueh = dist;
								cdgx = 1;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + VacancyRadius; 
					if(cdgx)
					{
						g[i]   +=  kb1*beishu1*2.0*overlap*( (x_i - xyph[0])*aa*aa + (y_i - xyph[1])*aa*bb )/( xyvalueh*(aa*aa+bb*bb) );
			            g[i+1] +=  kb1*beishu1*2.0*overlap*( (x_i - xyph[0])*aa*bb + (y_i - xyph[1])*bb*bb )/( xyvalueh*(aa*aa+bb*bb) );
                    
					}
					else
					{
						g[i]   +=  kb1*beishu1*2.0*overlap*(x_i - xyph[0])/xyvalueh;
			            g[i+1] +=  kb1*beishu1*2.0*overlap*(y_i - xyph[1])/xyvalueh;
                    
					}


				}
			}
			if(jiaodian%2==0 || yzbs==1)
			{		
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > VacancyRadius || fabs(y_i - holepoint[hh][J+1]) > VacancyRadius ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = VacancyRadius - dist;
						if(overlap > 0.0) 
                        {
                            g[i]   -=  beishu1*2.0*overlap*(x_i - holepoint[hh][J])/dist;
                            g[i+1] -=  beishu1*2.0*overlap*(y_i - holepoint[hh][J+1])/dist;
                        } 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = VacancyRadius - dist;
						   if(overlap > 0.0)  
					        {
			          	        g[i]   -=  beishu1*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) );
			                    g[i+1] -=  beishu1*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );	

					        } 
					    } 
					}

				}								
			}
				
	    }
	}
	{
		for (int J = 0; J < N; J++)
		{
   	        j   = 2*J;
			x_j = VS.x[j];
			y_j = VS.x[j+1];

            if(fabs(x_i - x_j) > 2.0*VacancyRadius || fabs(y_i - y_j) > 2.0*VacancyRadius ) continue; 

	 	    dist = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
	 	    overlap = 2.0*VacancyRadius - dist;
	 	    if(overlap > 0 )
	 	    {
			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;
			}
		}
	}  // Calculate the forces between the circles
  
    return ;
}

// The purpose of function myvalue_each() is to calculate the potential energy for each circle in the solution X. 
void myvalue_each
(
    double   *x,
    INT       n,
    double   *f
)
{
    double dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, i, j;
    for (i = 0; i < N; i++) f[i] = 0;  
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist);
        
		if(jds%2==0)
		{
			double xyp[2] = {pointlist[0],pointlist[1]}; 
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));

			for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud = (pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						xyvalue = newvalue;
					}
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b); 
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + R; 
				f[I] += kb*(overlap*overlap);
			}
		}
		 
		if(jds%2!=0 || yzbs==1) 
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	             
				if(fabs(x_i - pointlist[J]) > R || fabs(y_i - pointlist[J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = R - dist;
					if(overlap > 0.0) f[I] += (overlap*overlap); 
				} 
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0) f[I] += (overlap*overlap); 
				    } 
				}
		
			}
        }
        
         
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]);  
			yzbs = 0;
			 
			if(jiaodian%2!=0) 
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]};  
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 )  
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
						} 
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								xyvalueh = dist;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + R; 
					f[I] += kb*(overlap*overlap);
				} 
			}
			
			if(jiaodian%2==0 || yzbs==1)
			{		
				for(j=0;j<numPH[hh];j++)  
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > R || fabs(y_i - holepoint[hh][J+1]) > R ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = R - dist;
						if(overlap > 0.0) f[I] += (overlap*overlap); 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = R - dist;
						   if(overlap > 0.0) f[I] += (overlap*overlap); 
					    } 
					}
			
				}								
			}
				
	    }
	}
    
	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
		for (J = I+1; J < n/2; J++)
		{
			j = 2*J;
			x_j = x[j];
			y_j = x[j+1];

			if(fabs(x_i - x_j) > 2.0 * R || fabs(y_i - y_j) > 2.0 * R ) continue;   
		
			dist = sqrt( (x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * R - dist ;
			if(overlap > 0.0) 
			{
			    f[I] += (overlap*overlap); 
			    f[J] += (overlap*overlap);
		    }
		}

	}
	
    return ;	
}

// The purpose of function myvlaueNN() is to calculate the value of objective function E_D(X), 
// where only the neighboring vertices and edges of container and holes are used to calculate objective value
// for the sake of acceleration. 
double myvalueNN
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, i, j;
    f = 0.0 ; 
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
        //The points are by dedault considerd to locate in the container. 
		{
			for(j = 0; j < NNeighborsDD[I]; j++) 
			{
				J=2*AdjacentDD[I][j];
	
				if(fabs(x_i - pointlist[J]) > R || fabs(y_i - pointlist[J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i - pointlist[J])*(x_i - pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = R - dist;
					if(overlap > 0.0) f += (overlap*overlap); 
				}	
			} 
			
			for(j = 0; j < NNeighborsB[I]; j++) 
			{				
			    J=2*AdjacentB[I][j];
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0) f += (overlap*overlap); 
				    } 
				}
		
			}
        }
        //The points are by dedault considerd to locate out the holes. 
        for(int hh=0;hh<HN;hh++)
		{
		
			for(j=0;j<NNeighborsDDH[I][hh];j++) 
			{
				J = 2*AdjacentDDH[I][hh][j];
				
				if(fabs(x_i - holepoint[hh][J]) > R || fabs(y_i - holepoint[hh][J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
					overlap = R - dist;
					if(overlap > 0.0) f += (overlap*overlap); 
				}
			}
			for(j=0;j<NNeighborsBH[I][hh];j++) 
			{
				J = 2*AdjacentBH[I][hh][j];
				double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud<=0) 
				{
					if(holepoint[hh][J+2] == holepoint[hh][J])
					{
						a = 1.0;
						b = 0;
						c = -holepoint[hh][J];
					}
					else
					{
						a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
						b = -1.0;
				        c = holepoint[hh][J+1] - a*holepoint[hh][J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0) f += (overlap*overlap); 
				    } 
				}
		
			}							
				
	    }
	}
    f*=beishu;
    
	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
		for (J = 0; J < NNeighbors[I]; J++)
		{
			j = 2*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];
			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * R - dist ;
			if(overlap > 0.0) f += (overlap*overlap); 
		}

	}
	
    return (f) ;
}



// The purpose of function mygradNN() is to calculate the gradient of objective function E_D(X), 
// where only the neighboring vertices and edges of container and holes are used to calculate the gradient
// for the sake of acceleration. 
void mygradNN
(
    double    *g,
    double    *x,
    INT        n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, i, j;
	for (i = 0; i < n; i++) g [i] = 0;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		

         //The points are by dedault considerd to locate in the container. 
		{
			for(j = 0; j < NNeighborsDD[I]; j++) 
			{
				J=2*AdjacentDD[I][j];
	
				if(fabs(x_i - pointlist[J]) > R || fabs(y_i - pointlist[J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = R - dist;
	                if(overlap > 0.0)
	                {
	                    g[i]   -=  beishu*2.0*overlap*(x_i - pointlist[J])/dist;
	                    g[i+1] -=  beishu*2.0*overlap*(y_i - pointlist[J+1])/dist;					
	                } 
				}	
			} 
			
			for(j = 0; j < NNeighborsB[I]; j++) 
			{				
			    J=2*AdjacentB[I][j];
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0)
					   {
                            g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) ); 
                            g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );
					   }
				    } 
				}
		
			}
        }
        
        // The points are by dedault considerd to locate out the holes. 
        for(int hh=0;hh<HN;hh++)
		{
		
			for(j=0;j<NNeighborsDDH[I][hh];j++) 
			{
				J = 2*AdjacentDDH[I][hh][j];
				
				if(fabs(x_i - holepoint[hh][J]) > R || fabs(y_i - holepoint[hh][J+1]) > R ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
					overlap = R - dist;
						if(overlap > 0.0) 
                        {
                            g[i]   -=  beishu*2.0*overlap*(x_i - holepoint[hh][J])/dist;
                            g[i+1] -=  beishu*2.0*overlap*(y_i - holepoint[hh][J+1])/dist;
                        } 
				}
			}
			for(j=0;j<NNeighborsBH[I][hh];j++) 
			{
				J = 2*AdjacentBH[I][hh][j];
				double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud<=0)
				{
					if(holepoint[hh][J+2] == holepoint[hh][J])
					{
						a = 1.0;
						b = 0;
						c = -holepoint[hh][J];
					}
					else
					{
						a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
						b = -1.0;
				        c = holepoint[hh][J+1] - a*holepoint[hh][J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = R - dist;
					   if(overlap > 0.0)  
				       {
		          	        g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) );
		                    g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );	
				       } 
				    } 
				}
		
			}							
				
	    }
	}
    
	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
		for (J = 0; J < NNeighbors[I]; J++)
		{
			j = 2*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];
		
			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * R - dist ;
	 	    if(overlap > 0 )
	 	    {
			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
			}
		}

	}
	
    return  ;
}

// The purpose of this function is to calculat the objective value of the objective function 
// used in the sequential unconstrained minimization technique (SUMT) method 
// for adjusting the minimum distance between circles. 
double myvalueRR
(
    double   *x,
    INT       n
)
{

    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c;
    INT I, J, i, j;
    f = 0.0 ; 
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < (n-1)/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist); 
        
		if(jds%2==0)
		{
			double xyp[2] = {pointlist[0],pointlist[1]}; 
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));

			for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud = (pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						xyvalue = newvalue;
					}
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + x[n-1]; 
				
				f += kb*(overlap*overlap);
			}
		}
		
		if(jds%2!=0 || yzbs==1)
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	            
				if(fabs(x_i - pointlist[J]) > x[n-1] || fabs(y_i - pointlist[J+1]) > x[n-1] ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = x[n-1] - dist;
					if(overlap > 0.0) f += (overlap*overlap); 
				}
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0)
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b); 
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = x[n-1] - dist;
					   if(overlap > 0.0) f += (overlap*overlap); 
				    } 
				}
		
			}
        }
        
        
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			
			if(jiaodian%2!=0)
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]}; 
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
						}
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								xyvalueh = dist;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + x[n-1]; 
					f += kb*(overlap*overlap);
				}
			}
			
			if(jiaodian%2==0 || yzbs==1)
			{		
				for(j=0;j<numPH[hh];j++)  
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > x[n-1] || fabs(y_i - holepoint[hh][J+1]) > x[n-1] ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = x[n-1] - dist;
						if(overlap > 0.0) f += (overlap*overlap); 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0)
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = x[n-1] - dist;
						   if(overlap > 0.0) f += (overlap*overlap); 
					    } 
					}
			
				}								
			}
				
	    }
	}
    f*=beishu;
    
  
	for (I = 0; I < (n-1)/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		
		for (J = I+1; J < (n-1)/2; J++)
		{
			j = 2*J;
			x_j = x[j];
			y_j = x[j+1];

			if(fabs(x_i - x_j) > 2.0 * x[n-1] || fabs(y_i - y_j) > 2.0 * x[n-1] ) continue;   
		
			dist = sqrt( (x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 * x[n-1] - dist ;
			if(overlap > 0.0) f += (overlap*overlap); 
		}

	}
	

	f -= (x[n-1]*x[n-1]/alpha);
    return (f) ;
}


// The purpose of this function is to calculat the gradientof the objective function 
// used in the sequential unconstrained minimization technique (SUMT) method 
// for adjusting the minimum distance between circles. 
void mygradRR
(
    double    *g,
    double    *x,
    INT        n
)
{

    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    double a,b,c,aa,bb;
    INT I, J, i, j;
    for (i = 0; i < n; i++) g [i] = 0;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), R denotes the radius of CIRCLE.
   	for (I = 0; I < (n-1)/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];
		int jds = 0;
		int yzbs = 0;
        jds = judin(x_i,y_i,numP,pointlist); 
		if(jds%2==0) 
		{
			double xyp[2] = {pointlist[0],pointlist[1]}; 
			double xyvalue = sqrt((x_i - xyp[0]) * (x_i - xyp[0]) + (y_i - xyp[1]) * (y_i - xyp[1]));
			double cdgx = 0;
            for(j=0;j<numP;j++) 
			{
				J=2*j;
			    double maxx = MAX(pointlist[2+J],pointlist[J]);
			    double minx = MIN(pointlist[2+J],pointlist[J]);
			    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
	    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]); 
				if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
				{
					yzbs = 1; 
					break; 
				}
				else if(jud < 0 ) 
				{
					double newvalue = sqrt((x_i - pointlist[J]) * (x_i - pointlist[J]) + (y_i - pointlist[J+1]) * (y_i - pointlist[J+1])); 
					if(newvalue < xyvalue)
					{
						xyp[0] = pointlist[J];
						xyp[1] = pointlist[J+1];
						cdgx = 0;  
						xyvalue = newvalue;
					}
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   if(dist < xyvalue)
						{
							xyp[0] = cdx;
							xyp[1] = cdy;
							aa = a;
							bb = b;
							cdgx = 1;
							xyvalue = dist;
						}
				    } 
				}
			}
			if(yzbs==0)
			{
				overlap = xyvalue + x[n-1]; 
				if(cdgx)
				{
					g[i]   +=  kb*beishu*2.0*overlap*( (x_i - xyp[0])*aa*aa + (y_i - xyp[1])*aa*bb )/( xyvalue*(aa*aa+bb*bb) ); 
                    g[i+1] +=  kb*beishu*2.0*overlap*( (x_i - xyp[0])*aa*bb + (y_i - xyp[1])*bb*bb )/( xyvalue*(aa*aa+bb*bb) );
                    g[n-1] +=  kb*beishu*2.0*overlap;
				}
				else 
				{
					g[i]   +=  kb*beishu*2.0*overlap*(x_i - xyp[0])/xyvalue; 
                    g[i+1] +=  kb*beishu*2.0*overlap*(y_i - xyp[1])/xyvalue;
					g[n-1] +=  kb*beishu*2.0*overlap;					
				}

			}
		}
		if(jds%2!=0 || yzbs==1) 
		{
			for(j=0;j<numP;j++) 
			{
				J=2*j;
	
				if(fabs(x_i - pointlist[J]) > x[n-1] || fabs(y_i - pointlist[J+1]) > x[n-1] ) 
				{
				}
				else
				{
					dist = sqrt( (x_i-pointlist[J])*(x_i-pointlist[J]) + (y_i - pointlist[J+1])*(y_i - pointlist[J+1]) );
					overlap = x[n-1] - dist;
	                if(overlap > 0.0)
	                {
	                    g[i]   -=  beishu*2.0*overlap*(x_i - pointlist[J])/dist;
	                    g[i+1] -=  beishu*2.0*overlap*(y_i - pointlist[J+1])/dist;
						g[n-1] +=  beishu*2.0*overlap;				
	                } 
				}	
				
				double jud=(pointlist[2+J] - pointlist[J]) * (y_i - pointlist[J+1]) - (pointlist[2+J+1] - pointlist[J+1]) * (x_i - pointlist[J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
				if(jud>=0) 
				{
					if(pointlist[J+2] == pointlist[J])
					{
						a = 1.0;
						b = 0;
						c = -pointlist[J];
					}
					else
					{
						a = (pointlist[J+3] - pointlist[J+1]) / (pointlist[J+2] - pointlist[J]);
						b = -1.0;
				        c = pointlist[J+1] - a*pointlist[J];
					}
				    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
				    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
	
				    double maxx = MAX(pointlist[2+J],pointlist[J]);
				    double minx = MIN(pointlist[2+J],pointlist[J]);
				    double maxy = MAX(pointlist[2+J+1],pointlist[J+1]);
		    	    double miny = MIN(pointlist[2+J+1],pointlist[J+1]);
				    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
				    {
					   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
					   overlap = x[n-1] - dist;
					   if(overlap > 0.0)
					   {
                            g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) ); 
                            g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );
                            g[n-1] +=  beishu*2.0*overlap;
					   }
				    } 
				}
		
			}
        }
        
        for(int hh=0;hh<HN;hh++)
		{
	        int jiaodian = 0; 
			jiaodian = judin(x_i,y_i,numPH[hh],holepoint[hh]); 
			yzbs = 0;
			if(jiaodian%2!=0)
			{
				double xyph[2] = {holepoint[hh][0],holepoint[hh][1]}; 
				double xyvalueh = sqrt((x_i - xyph[0]) * (x_i - xyph[0]) + (y_i - xyph[1]) * (y_i - xyph[1]));
				double cdgx = 0;
	
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
				    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
				    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
				    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
		    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					double jud = (holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]); 
					if(jud == 0 &&x_i <= maxx && x_i >= minx && y_i <= maxy && y_i >= miny) 
					{
						yzbs = 1; 
						break; 
					}
					else if(jud > 0 ) 
					{
						double newvalueh = sqrt((x_i - holepoint[hh][J]) * (x_i - holepoint[hh][J]) + (y_i - holepoint[hh][J+1]) * (y_i - holepoint[hh][J+1])); 
						if(newvalueh < xyvalueh)
						{
							xyph[0] = holepoint[hh][J];
							xyph[1] = holepoint[hh][J+1];
							xyvalueh = newvalueh;
							cdgx = 0;
						}
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   if(dist < xyvalueh)
							{
								xyph[0] = cdx;
								xyph[1] = cdy;
								aa = a;
								bb = b;
								xyvalueh = dist;
								cdgx = 1;
							}
					    } 
					}
				}
				if(yzbs==0)
				{
					overlap = xyvalueh + x[n-1]; 
					if(cdgx)
					{
						g[i]   +=  kb*beishu*2.0*overlap*( (x_i - xyph[0])*aa*aa + (y_i - xyph[1])*aa*bb )/( xyvalueh*(aa*aa+bb*bb) );
			            g[i+1] +=  kb*beishu*2.0*overlap*( (x_i - xyph[0])*aa*bb + (y_i - xyph[1])*bb*bb )/( xyvalueh*(aa*aa+bb*bb) );
			            g[n-1] +=  kb*beishu*2.0*overlap;
					}
					else
					{
						g[i]   +=  kb*beishu*2.0*overlap*(x_i - xyph[0])/xyvalueh;
			            g[i+1] +=  kb*beishu*2.0*overlap*(y_i - xyph[1])/xyvalueh;
			            g[n-1] +=  kb*beishu*2.0*overlap;
					}


				}
			}
			if(jiaodian%2==0 || yzbs==1)
			{		
				for(j=0;j<numPH[hh];j++) 
				{
					J=2*j;
					
					if(fabs(x_i - holepoint[hh][J]) > x[n-1] || fabs(y_i - holepoint[hh][J+1]) > x[n-1] ) 
					{
					}
					else
					{
						dist = sqrt( (x_i-holepoint[hh][J])*(x_i-holepoint[hh][J]) + (y_i - holepoint[hh][J+1])*(y_i - holepoint[hh][J+1]) );
						overlap = x[n-1] - dist;
						if(overlap > 0.0) 
                        {
                            g[i]   -=  beishu*2.0*overlap*(x_i - holepoint[hh][J])/dist;
                            g[i+1] -=  beishu*2.0*overlap*(y_i - holepoint[hh][J+1])/dist;
                            g[n-1] +=  beishu*2.0*overlap; 
                        } 
					}
										
					double jud=(holepoint[hh][2+J] - holepoint[hh][J]) * (y_i - holepoint[hh][J+1]) - (holepoint[hh][2+J+1] - holepoint[hh][J+1]) * (x_i - holepoint[hh][J]);//≤Ê≥À≈–∂œ «∑Ò‘⁄◊Û ÷±ﬂ 
					if(jud<=0) 
					{
						if(holepoint[hh][J+2] == holepoint[hh][J])
						{
							a = 1.0;
							b = 0;
							c = -holepoint[hh][J];
						}
						else
						{
							a = (holepoint[hh][J+3] - holepoint[hh][J+1]) / (holepoint[hh][J+2] - holepoint[hh][J]);
							b = -1.0;
					        c = holepoint[hh][J+1] - a*holepoint[hh][J];
						}
					    double cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
					    double cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
		
					    double maxx = MAX(holepoint[hh][2+J],holepoint[hh][J]);
					    double minx = MIN(holepoint[hh][2+J],holepoint[hh][J]);
					    double maxy = MAX(holepoint[hh][2+J+1],holepoint[hh][J+1]);
			    	    double miny = MIN(holepoint[hh][2+J+1],holepoint[hh][J+1]);
					    if(cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny)
					    {
						   dist = sqrt( (x_i - cdx) * (x_i - cdx) + (y_i - cdy) * (y_i - cdy) );
						   overlap = x[n-1] - dist;
						   if(overlap > 0.0)  
					        {
			          	        g[i]   -=  beishu*2.0*overlap*( (x_i - cdx)*a*a + (y_i - cdy)*a*b )/( dist*(a*a+b*b) );
			                    g[i+1] -=  beishu*2.0*overlap*( (x_i - cdx)*a*b + (y_i - cdy)*b*b )/( dist*(a*a+b*b) );	
			                    g[n-1] +=  beishu*2.0*overlap;
					        } 
					    } 
					}

				}								
			}
				
	    }
	}

	for (I = 0; I < (n-1)/2; I++)
	{
 	    i   = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		for (J = I+1; J < (n-1)/2; J++)
		{
   	        j   = 2*J;
			x_j = x[j];
			y_j = x[j+1];

            if(fabs(x_i - x_j) > 2.0*x[n-1] || fabs(y_i - y_j) > 2.0*x[n-1] ) continue; 

	 	    dist = sqrt( (x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) );
	 	    overlap = 2.0*x[n-1] - dist;
	 	    if(overlap > 0 )
	 	    {
			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
	         
	         g[n-1] +=  4.0*overlap;
			}
		}
	}  // Calculate the forces between the circles
    
	g[n-1] -=  (2.0*x[n-1]/alpha);
           
    return ;
}
