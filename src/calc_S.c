/*------------------------------------------------------------------------
 *  Calculate slowness from P-wave velocity model
 *
 *  D. Koehn
 *  Kiel, 12.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_S(float **  Vp, float **  S){

	extern int NX, NY;
	
	/* local variables */
	int i, j;

	/* loop over global grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){

                         S[j][i] = 1.0/Vp[j][i];				

		}
	}

}




