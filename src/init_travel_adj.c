/*------------------------------------------------------------------------
 *   Apply boundary conditions for wavefield lambda
 *  
 *  
 *   D. Koehn
 *   Kiel, 11/12/2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_travel_adj(float **lam, int ** recpos, int ntr, float * Tres){

	register int i, j, k, l;
	extern int FW, NX, NY;
        extern float TIME;

        /* set whole matrix to zero */
        init_grad(lam);
	
        /* initialize the interior grid points with non-zero 
           values while preserving the boundary condition for lam */
	for (j=2;j<=NY-1;j++){
		for (i=2;i<=NX-1;i++){

                     lam[j][i] = 2.0 * TIME;

                     /* add data residuals at receiver positions */
                     for (l=1;l<=ntr;l++){
                           
                           if((i==recpos[1][l])&&(j==recpos[2][l])){
			      lam[j][i] = Tres[l];
                           }

                     }

		}
	}           
	
}
