/*------------------------------------------------------------------------
 *  Calculate first arrival traveltimes at receiver positions                           
 *  
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_FA(float ** TT, int ** recpos, int ntr, float * Tmod){

	/* local variables */
	int i;

	for(i=1;i<=ntr;i++){
	    Tmod[i] = TT[recpos[2][i]][recpos[1][i]];
	}
	
}
