/*------------------------------------------------------------------------
 *  Calculate first arrival traveltime residuals and objective function                           
 *  
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double calc_FA_res(float **TT, float * Tmod, float * Tobs, float * Tres, int ntr, int ishot, int ** recpos){

        /* global variables */
	extern char DATA_DIR[STRING_SIZE];
	extern float DH;
	extern int NX, NY;

	/* local variables */
	int i;
        char pickfile_char[STRING_SIZE];
        float tmp, minres, gradT;
	double l2;

        minres = 1e-6;

        FILE *fp;

        /* open FA pick file */
	sprintf(pickfile_char,"%s_%i.dat",DATA_DIR,ishot);
	fp=fopen(pickfile_char,"r");
	if (fp == NULL) {
		err(" picks_?.dat could not be opened !");
	}

	for(i=1;i<=ntr;i++){

            /* read FA traveltime data */
            fscanf(fp,"%e",&Tobs[i]);  

            /* calculate traveltime residuals */
	    Tres[i] = Tmod[i] - Tobs[i];

            //if(fabs(Tres[i])<=minres){Tres[i] = 0.0;}

            /* calculate objective function */
            l2 += Tres[i] * Tres[i]; 

	    /* scale traveltime residual according to eq. 12 in Taillandier et al. (2009)
	       First-arrival traveltime tomography based on the adjoint-state method
	       lambda = (Tmod - Tobs) / (n * grad(TT)) */
	    if (recpos[2][i] <= NY - 1){

		gradT = (TT[recpos[2][i] + 1][recpos[1][i]] - TT[recpos[2][i]][recpos[1][i]]) / DH;

		if (fabs(gradT) > 0.0)
			Tres[i] = Tres[i] / gradT;
		else
			Tres[i] = 0.0;

	    }else if (NY == recpos[2][i]){

		gradT = (TT[recpos[2][i]][recpos[1][i]] - TT[recpos[2][i] - 1][recpos[1][i]]) / DH;

		if (fabs(gradT) > 0.0)
			Tres[i] = Tres[i] / gradT;
		else
			Tres[i] = 0.0;
	    }

	}

        fclose(fp);

        return l2;
}
