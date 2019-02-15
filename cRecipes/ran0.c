
#include <stdio.h>
#include "cRecipes.h"
#define MAXRAND 32767  
    real ran0(int *idum)
/*
   Random number generator based around rand().  From
   NUMERICAL RECIPES IN C.
   Input idum < 0 to initialize.
   Returns a number between 0 and 1.0.
*/
{
	static float y,maxran,v[98];
	float dum;
	static int iff=0;
	int j;
        int iev;
	unsigned i,k;
 
	if (*idum < 0 || iff == 0) {
		iff=1;
		i=2;
		do {
			k=i;
			i<<=1;
		} while (i);
		maxran=k;
		srand(*idum);
		*idum=1;
		for (j=1;j<=97;j++) dum=rand();
		for (j=1;j<=97;j++) v[j]=rand();
		y=rand();
	}
	j=1+97.0*y/maxran;
	if (j > 97 || j < 1) error("%s\n",
            "ran0 -- Error in random # generator");
	y=v[j];
	v[j]=rand();
      
	return (real) (y/MAXRAND);
}
