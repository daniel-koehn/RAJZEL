/*------------------------------------------------------------------------
 *  Calculate gradient and objective function value
 *
 *  D. Koehn
 *  Kiel, 13.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float grad_obj(float ** grad, float ** S, float ** TT, float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr, int iter){

	/* declaration of global variables */
        extern int SEISMO, NX, NY, NSHOT1, NSHOT2;
    
        /* declaration of local variables */
        int ishot; 
        float L2sum, L2, ** grad_shot;
        grad_shot =  matrix(1,NY,1,NX);

        /* suppress modelled FA traveltime output during FATT */
        SEISMO = 0;

	/* set gradient matrix to zero before next iteration*/
	init_grad(grad);
        L2 = 0.0;

	for (ishot=NSHOT1;ishot<NSHOT2;ishot++){

	     /* set gradient for each shot to zero before next iteration */
	     init_grad(grad_shot); 

	     /* solve forward problem with Eikonal solver */
	     eikonal(S,TT,Tmod,srcpos,ishot,nshots,recpos,ntr);

             /* calculate traveltime residuals at receiver positions */
 	     L2+=calc_FA_res(Tmod,Tobs,Tres,ntr,ishot);

	     /*write_picks(Tres,ntr,ishot);*/
                
             /* calculate adjoint Eikonal wavefield */
             eikonal_adj(TT,Tres,lam,recpos,ntr,ishot); 

             /* assemble gradient for each shot */
	     ass_grad(grad,grad_shot,S,lam,srcpos,nshots,recpos,ntr,ishot);

	} /* end of loop over shots (forward and adjoint) */   

	/* sprintf(filename,"%s_adj.tt",SNAP_FILE);
	writemod(filename,grad,3); */

        /* assemble objective function from all MPI processes */
	L2sum = 0.0;
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        L2 = 0.5 * L2sum;
	
	/*if(MYID==0){
	  printf("L2sum: %e\n", L2sum);
        }*/

        precond(grad,nshots,srcpos,recpos,ntr,iter);

        /* deallocate memory */
        free_matrix(grad_shot,1,NY,1,NX);

return L2;
                	    
}
