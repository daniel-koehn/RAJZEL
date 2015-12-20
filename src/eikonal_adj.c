/*------------------------------------------------------------------------
 *  Eikonal solver for the adjoint traveltime wavefield
 *
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void eikonal_adj(float ** TT, float * Tres, float ** lam, int ** recpos, int ntr, int ishot){

	/* declaration of global variables */
        extern float DH, TIME, TTNORM;
        extern int NX, NY, MYID, INFO;
        extern char SNAP_FILE[STRING_SIZE];
	
        /* declaration of local variables */
        int k;
        float lnorm1, ** lamold;
        char filename[STRING_SIZE];

        lamold =  matrix(1,NY,1,NX);

        if((MYID==0)&&(INFO==1)){
	  printf("\n==============================================\n");
          printf("\n *****  Solve adjoint Eikonal equation  ***** \n");
	  printf("\n==============================================\n\n");				
        }		                                                 

	/* initialize boundary condition for lambda */
	init_travel_adj(lam,recpos,ntr,Tres);

        /* initialize l1 norm of lam - lamold */
        lnorm1 = 10.0 * TIME;

	/* apply fast sweeping method to solve the eikonal equation */
	while(lnorm1>=TTNORM){

                /* save old TT values */
                store_mat(lam,lamold,NX,NY);

                /* sweep with order according to Zhao (2004) */
                sweep_adj(lam,TT,Tres,recpos,ntr,2,NX-1,1,2,NY-1,1);
                sweep_adj(lam,TT,Tres,recpos,ntr,NX-1,2,-1,2,NY-1,1);
                sweep_adj(lam,TT,Tres,recpos,ntr,NX-1,2,-1,NY-1,2,-1);
                sweep_adj(lam,TT,Tres,recpos,ntr,2,NX-1,1,NY-1,2,-1);

		/* calculate l1 norm of lam - lamold */
		lnorm1 = norm1(lam,lamold);

	}

	/* output of First arrival traveltimes*/
        /* sprintf(filename,"%s_adj_shot_%d.tt",SNAP_FILE,ishot);
	writemod(filename,lam,3); */

	free_matrix(lamold,1,NY,1,NX);
    	    
}
