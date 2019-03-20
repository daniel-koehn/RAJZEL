/*
 * Assemble gradient from all shots  
 *
 * Daniel Koehn
 * Kiel, 11.12.2015
 *
 */

#include "fd.h"

void ass_grad(float ** grad, float ** grad_shot, float ** S, float ** lam, float **srcpos,  int nshots, int **recpos, int ntr, int ishot){

        /* global variables */
	extern int NX, NY, SWS_TAPER_CIRCULAR_PER_SHOT;

	/* local variables */
	int i, j;
	
	/* assemble gradient from one shots */
	for (i=2;i<=NX-1;i++){
	    for (j=2;j<=NY-1;j++){

               grad_shot[j][i] = lam[j][i]*pow(S[j][i],3.0);
		    
	    }
	}
	
	cp_grad_frame(grad_shot);
	
	/* apply taper at source positions */
	if(SWS_TAPER_CIRCULAR_PER_SHOT==1){
	  taper_grad_shot(grad_shot,srcpos,nshots,recpos,ntr,ishot);
	}

	if(SWS_TAPER_CIRCULAR_PER_SHOT==2){
	  taper_grad_shot1(grad_shot,srcpos,nshots,recpos,ntr,ishot);
	}
	
	/* assemble gradient from all shots */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

	       grad[j][i] += grad_shot[j][i];
		    
	    }
	}
	

}



