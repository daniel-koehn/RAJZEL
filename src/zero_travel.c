/*------------------------------------------------------------------------
 *   Zero FA traveltimes
 *  
 *  
 *   D. Koehn
 *   Kiel, 09/12/2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_travel(float **TT){

	register int i, j, k;
	extern int FW, NX, NY;
        extern float TIME;

	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

                     TT[j][i] = 2.0 * TIME;

		}
	}           
	
}
