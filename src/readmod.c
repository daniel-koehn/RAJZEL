/*------------------------------------------------------------------------
 *  Read P-wave velocity model
 *
 *  D. Koehn
 *  Kiel, 7.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readmod(float **  Vp){

	extern int NX, NY, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	int i, j;
	float vp;
	FILE *fp_vs, *fp_vp, *fp_rho;
	char filename[STRING_SIZE];


	fprintf(FP,"\n... reading P-wave velocity models from file...\n");
           
	/* read density and seismic velocities */
	/* ----------------------------------- */
	fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
	sprintf(filename,"%s.vp",MFILE);
	fp_vp=fopen(filename,"r");
	if (fp_vp==NULL) err(" Could not open model file for Vp ! ");           
	   
	/* loop over global grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			Vp[j][i] = vp;				
		}
	}

	fclose(fp_vp);
	
	
	/* each PE writes his model to disk */
	sprintf(filename,"%s.rajzel.vp",MFILE);
	writemod(filename,Vp,3);

}




