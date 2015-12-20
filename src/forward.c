/*------------------------------------------------------------------------
 *  Solve Eikonal forward problem with Fast sweeping technique
 *
 *  D. Koehn
 *  Kiel, 13.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void forward(float ** S, float ** TT, float * Tmod,  float ** srcpos, int nshots, int ** recpos, int ntr){

	/* declaration of global variables */
        extern int NSHOT1, NSHOT2;
    
        /* declaration of local variables */
        int ishot;

        /* loop over all shots */
	for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

		/* solve forward problem with Eikonal solver */
		eikonal(S,TT,Tmod,srcpos,ishot,nshots,recpos,ntr);

	}
         	    
}
