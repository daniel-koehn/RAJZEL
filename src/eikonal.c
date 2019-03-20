/*------------------------------------------------------------------------
 *  Eikonal solver
 *
 *  D. Koehn
 *  Kiel, 10/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void eikonal(float ** S, float ** TT, float * Tmod,  float ** srcpos, int ishot, int nshots, int ** recpos, int ntr){

	/* declaration of global variables */
        extern float DH, TIME, TTNORM;
        extern int NX, NY, SNAP, SEISMO, MYID, INFO;
        extern char SNAP_FILE[STRING_SIZE];
	
        /* declaration of local variables */
        int nxsrc, nysrc, k;
        char filename[STRING_SIZE];
        float lnorm1, ** TTold;

        TTold =  matrix(1,NY,1,NX);

        if((MYID==0)&&(INFO==1)){
	  printf("\n===================================================================================\n");
          printf("\n *****  Solve Eikonal equation for shot %d of %d for MPI process no. %d ********** \n",ishot,nshots, MYID);
	  printf("\n===================================================================================\n\n");				
        }		           
                                      
	/* Eikonal solver */
	nxsrc = iround(srcpos[1][ishot]/DH);
	nysrc = iround(srcpos[2][ishot]/DH); 
	zero_travel(TT);

	/* initialize traveltime field at source positions */
	TT[nysrc][nxsrc] = 0.0;

        /* initialize l1 norm of TT - TTold */
        lnorm1 = 10.0 * TIME;

	/*printf("nxsrc = %d \t nysrc = %d \n",nxsrc,nysrc);*/

	/* apply fast sweeping method to solve the eikonal equation */
	while(lnorm1>=TTNORM){

                /* save old TT values */
                store_mat(TT,TTold,NX,NY);

                /* sweep with order according to Zhao (2004) */
		sweep(S,TT,nxsrc,nysrc,1,NX,1,1,NY,1);
		sweep(S,TT,nxsrc,nysrc,NX,1,-1,1,NY,1);
		sweep(S,TT,nxsrc,nysrc,NX,1,-1,NY,1,-1);
		sweep(S,TT,nxsrc,nysrc,1,NX,1,NY,1,-1);

		/* calculate l1 norm of TT - TTold */
		lnorm1 = norm1(TT,TTold);

	}

	/* output of First arrival traveltimes*/
	if(SNAP==1){
		sprintf(filename,"%s_shot_%d.tt",SNAP_FILE,ishot);
		writemod(filename,TT,3);
	}

        /* calculate TT at receiver positions */
        calc_FA(TT,recpos,ntr,Tmod);

	/* write pick file compatible with DENISE Black-Edition */
	if(SEISMO==1){
		write_picks(Tmod,ntr,ishot);
	}

	free_matrix(TTold,1,NY,1,NX);
    	    
}
