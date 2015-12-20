/* 
   Reading (distributed) source positions, timeshift, centre frequency 
   and amplitude from SOURCE_FILE.
   
   written by T. Bohlen 
   last update: 12.02.2002
*/

#include "fd.h"

float **sources(int *nsrc){

	/* declaration of global variables */
	extern float DH;
	extern  char SOURCE_FILE[STRING_SIZE];
	extern int MYID, NXG, NYG, QUELLTYP;
	extern FILE *FP;

        /* declaration of local variables */
	float **srcpos;
	int   i, l, isrc=0, current_source=0, nvarin=0,nsrc_dummy;
	float xsrc, ysrc, zsrc, tshift, fc, amp, tan_phi, dz, x;
	char  cline[256];
	FILE *fpsrc;


	if (MYID==0){
			fprintf(FP,"\n Reading source positions, time-shift, centre frequency \n");
			fprintf(FP," and amplitude from file: %s\n",SOURCE_FILE);
			fpsrc=fopen(SOURCE_FILE,"r");
	
			if (fpsrc==NULL) err(" Source file could not be opened !");
			*nsrc=0;

			
			/* read number of source positions */	
			fscanf(fpsrc,"%d",nsrc);
						
			fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,*nsrc);
			
			srcpos=matrix(1,8,1,*nsrc);

			rewind(fpsrc);		
			fgets(cline,255,fpsrc);		/* Dummy fgets for ignoring first line */

			for (l=1;l<=*nsrc;l++){
				fgets(cline,255,fpsrc);

				nvarin=sscanf(cline,"%f%f%f%f%f%f%f%f",&xsrc, &zsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l], &srcpos[8][l]);
				switch(nvarin){
					case 0: xsrc=0.0;
					case 1: zsrc=0.0;
					case 2: ysrc=0.0;
					case 3: if (MYID==0) fprintf(FP," No time shift defined for source %i in %s!\n",l, SOURCE_FILE);
						err("Missing parameter in SOURCE_FILE!");
					case 4: if (MYID==0) fprintf(FP," No frequency defined for source %i in %s!\n",l, SOURCE_FILE);
						err("Missing parameter in SOURCE_FILE!");
					case 5: if (MYID==0) fprintf(FP," No amplitude defined for source %i in %s!\n",l, SOURCE_FILE);
						err("Missing parameter in SOURCE_FILE!");
					case 6: srcpos[7][l]=0.0;
					/*case 7: srcpos[8][l]=QUELLTYP;*/
				}
				if ((srcpos[8][l]!=4) && (nvarin>6)) {
				current_source=(int)srcpos[8][l];
				/*if (MYID==0) fprintf(FP," SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n", l, current_source);*/
				}
				/* fscanf(fpsrc,"%f%f%f%f%f",&xsrc, &ysrc, &tshift, &fc, &amp); */ 
				srcpos[1][l]=xsrc;
				srcpos[2][l]=ysrc;
				srcpos[3][l]=0.0;
				srcpos[4][l]=tshift;
				fc=srcpos[5][l];
			}

/*	fscanf(fpsrc,"%f%f%f%f%f%f",&xsrc, &zsrc, &ysrc, &tshift, &fc, &amp);*/ 	/* old implementation */
			/*	srcpos[1][l]=xsrc;
				srcpos[2][l]=ysrc;
				srcpos[3][l]=zsrc;
				srcpos[4][l]=tshift;
				srcpos[5][l]=fc;
				srcpos[6][l]=amp;
			} */

			fclose(fpsrc);

			/* Compute maximum frequency */
			for (l=1;l<=*nsrc;l++)
				if (srcpos[5][l]>fc) fc=srcpos[5][l];
			fprintf(FP," Maximum frequency defined in %s: %6.2e Hz\n",SOURCE_FILE,fc);

			/* outputs all sources per each subdomain / node*/
		
			if (MYID==0){
				/*fprintf(FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");
	   			for (l=1;l<=*nsrc;l++)
		      			fprintf(FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n\n",
					l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);*/
	      		}
				
			fprintf(FP," Message from function sources (written by PE %d):\n",MYID);			

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
	if (MYID!=0) srcpos=matrix(1,8,1,*nsrc);
	MPI_Bcast(&srcpos[1][1],(*nsrc)*8,MPI_FLOAT,0,MPI_COMM_WORLD);

/*	if (MYID==0){
		fprintf(FP,"\n **Message from function source (written by PE %d):\n",MYID);
		fprintf(FP," Number of global source positions found: %i\n",*nsrc);
		fprintf(FP," x\t\ty\t\tz\t\ttshift\t\tfc\t\tamp\n");
		for (l=1;l<=*nsrc;l++)
			fprintf(FP," %6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\n",
					srcpos[1][l],srcpos[2][l],srcpos[3][l],srcpos[4][l],srcpos[5][l],srcpos[6][l]);
		fprintf(FP,"\n\n");
	}
*/

	return srcpos;
}
