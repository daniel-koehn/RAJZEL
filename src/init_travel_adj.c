/*------------------------------------------------------------------------
 *   Apply boundary conditions for wavefield lambda
 *  
 *  
 *   D. Koehn
 *   Kiel, 11/12/2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_travel_adj(float **lam, int ** recpos, int ntr, float * Tres, float ** TT){

	register int i, j, k, l;
	extern int FW, NX, NY, MYID;
        extern float TIME, DH;

        float ngradTT, eps;
        eps=0.0;

        /* set whole matrix to zero */
        init_grad(lam);

        /* apply boundary condition if receivers are placed at the boundaries */
        /*for (l=1;l<=ntr;l++){

             if(recpos[2][l]==1){
                
                for(i=2;i<=NX-1;i++){
		    if(i==recpos[1][l]){

                      ngradTT = -((TT[2][i]-TT[1][i])/DH);
		      lam[1][i] = Tres[l]/(ngradTT+eps);

		    }
                }
                
             }

        }*/
	
        /* initialize the interior grid points with non-zero 
           values while preserving the boundary condition for lam */
	for (j=2;j<=NY-1;j++){
		for (i=2;i<=NX-1;i++){

                     lam[j][i] = 10.0 * TIME;

                     for (l=1;l<=ntr;l++){
                           
                           if((i==recpos[1][l])&&(j==recpos[2][l])){
			      lam[j][i] = Tres[l];
                           }

                     }

		}
	}          
	
}
