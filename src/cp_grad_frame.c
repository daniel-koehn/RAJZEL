/*
 * Copy gradient values into boundary frame  
 *
 * Daniel Koehn
 * Kiel, 11/01/2016
 */

#include "fd.h"

void cp_grad_frame(float ** A){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j, jj;

	/* Copy top boundary frame */
	for (i=1;i<=NX;i++){
	    
	       A[2][i] = A[3][i];
               A[1][i] = A[2][i];
		    
	}
	
	/* Copy bottom boundary frame */
	for (i=1;i<=NX;i++){
	    
               A[NY-1][i] = A[NY-2][i];
	       A[NY][i] = A[NY-1][i];
		    
	}
	
	/* Copy left boundary frame */
	for (j=1;j<=NY;j++){
	    
               A[j][2] = A[j][3];
	       A[j][1] = A[j][2];
		    
	}

	
	/* Copy right boundary frame */
	for (j=1;j<=NY;j++){
	    
               A[j][NX-1] = A[j][NX-2];
	       A[j][NX] = A[j][NX-1];
		    
	}

        /* printf("A[NY][NX] = %e \n",A[NY][NX]); */		

}
