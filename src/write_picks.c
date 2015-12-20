/*------------------------------------------------------------------------
 *  Write first arrival traveltime picks                           
 *  
 *  D. Koehn
 *  Kiel, 09/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void write_picks(float * Tmod, int ntr, int ishot){

	/* global variables */
        extern char PICKS_FILE[STRING_SIZE];        
        extern int MYID, INFO;

	/* local variables */
	int i;
        char pickfile_char[STRING_SIZE];
	
        FILE *fp;

        if((MYID==0)&&(INFO==1)){
	    printf("\n Write FA traveltimes at receiver positions ...\n");
        }

        sprintf(pickfile_char,"%s_%i.dat",PICKS_FILE,ishot);

        fp=fopen(pickfile_char,"w");

	for(i=1;i<=ntr;i++){
	    fprintf(fp,"%e \n",Tmod[i]); 
	    /* fwrite(&Tmod[i], sizeof(float), 1, fp); */
	}

        fclose(fp);
	
}
