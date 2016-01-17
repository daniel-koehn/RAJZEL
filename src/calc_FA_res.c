/*------------------------------------------------------------------------
 *  Calculate first arrival traveltime residuals and objective function                           
 *  
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float calc_FA_res(float * Tmod, float * Tobs, float * Tres, int ntr, int ishot){

        /* global variables */
	extern char DATA_DIR[STRING_SIZE];

	/* local variables */
	int i;
        char pickfile_char[STRING_SIZE];
        float l2, tmp, minres;

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
	    Tres[i] = Tobs[i] - Tmod[i];

            if(fabs(Tres[i])<=minres){Tres[i] = 0.0;}

            /* calculate objective function */
            l2 += Tres[i] * Tres[i]; 
	}

        fclose(fp);

        return l2;
}
