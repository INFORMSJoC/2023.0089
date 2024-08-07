#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

// The function shake() aims to perturb randomly the input solution S. 
void Shake(Solution &S, int n, double Max_delt)
{

	for(int i=0; i<n; i++)
	{
        S.x[i] += (2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0)*Max_delt;
	}
    
}

