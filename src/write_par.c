/*------------------------------------------------------------------------
 *   Write FD-Parameters to stdout                           
 *   last update 20/02/2001
 *
 *  T. Bohlen
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){

	/* declaration of global variables */
	extern int   NX, NY, NT;
	extern int  SNAP, SNAP_FORMAT, TAPER;
	extern float DH, TIME;
	extern int SEISMO;
	extern int  READMOD, INVTYPE;
	extern float REFREC[4];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char MFILE[STRING_SIZE], JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE];
        extern char PICKS_FILE[STRING_SIZE];
	extern int NP, MYID;
	
	extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, INVMAT;
	extern int GRAD_METHOD;
	extern int FILT_SIZE, MODEL_FILTER;
	extern int FILT_SIZE_GRAD, GRAD_FILTER;
	
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
	extern int SWS_TAPER_FILE;
	extern float SRTRADIUS, EXP_TAPER_GRAD_HOR;
	extern int MIN_ITER;
	extern char INV_MODELFILE[STRING_SIZE];
	extern int nfstart, nf;
	extern int nfstart_jac, nf_jac;
	extern float VPUPPERLIM, VPLOWERLIM;
	
	extern int STEPMAX;
	extern float EPS_SCALE, SCALEFAC, TTNORM, EPS_ADJ;

	extern int NORMALIZE, NLBFGS, N_STREAMER;
        extern float REC_INCR_X, REC_INCR_Y, C1, C2;
	
	extern char MISFIT_LOG_FILE[STRING_SIZE];

        extern float VP0_1, VP0_2, DVP0, GRAD0_1, GRAD0_2, DGRAD0;
        extern char GRIDSEARCH_FILE[STRING_SIZE];

	/* definition of local variables */
	int l;
	
        fprintf(fp," -----------------------  RAJZEL modus  ----------------------\n");
	if (INVMAT==0){
		fprintf(fp," INVMAT=%d: Only solving the Eikonal forward problem.\n",INVMAT);}
	else {
		if(INVMAT==1){
			fprintf(fp," INVMAT=%d: FATT using the adjoint state method.\n",INVMAT);}
		if(INVMAT==2){
			fprintf(fp," INVMAT=%d: Grid Search for 1D linear gradient model.\n",INVMAT);}
	}
        fprintf(fp,"\n");
	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
	fprintf(fp," Grid-spacing (DH): %e meter\n", DH);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");


	fprintf(fp," reading source positions, time delay, centre frequency \n");
	fprintf(fp," and initial amplitude from ASCII-file \n");
	fprintf(fp,"\t%s\n\n",SOURCE_FILE);

	if (SEISMO){
		fprintf(fp," ------------------------- RECEIVER  --------------------------\n");

		fprintf(fp," reading receiver positions from file \n");
		fprintf(fp,"\t%s\n\n",REC_FILE);
		fprintf(fp," reference_point_for_receiver_coordinate_system:\n");
		fprintf(fp," x=%f \ty=%f\t z=%f\n\n",REFREC[1], REFREC[2], REFREC[3]);

                if (N_STREAMER){
                        fprintf(fp," ------------------------- Towed streamer  --------------------------\n");
                        fprintf(fp," Assuming the first N_STREAMER = %d receivers belong to a streamer \n",N_STREAMER);
                        fprintf(fp," Shifting the streamer by ... \n");
                        fprintf(fp," REC_INCR_X=%f in x-direction \n",REC_INCR_X);
                        fprintf(fp," REC_INCR_Y=%f in y-direction \n",REC_INCR_Y);
                }
	}

	if (READMOD){
		fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
		fprintf(fp," names of model-files: \n");
		fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
		fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);
		fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
		fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
		fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);
	}

	if (SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(SNAP){
		case 1:
			fprintf(fp," Write traveltime field");
			break;
		default:
			err(" sorry, incorrect value for SNAP ! \n");
		}

		fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,NX);
		fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,NY);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
		switch (SNAP_FORMAT){
		case 1 :
			err(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			err(" Don't know the format for the Snapshot-data ! \n");
		}
	
		fprintf(fp,"\n\n");
	}
	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  FA travel times  ----------------------\n");
		if (SEISMO){
			fprintf(fp," write FA travel times at receiver positions to ...");
   			fprintf(fp," %s",PICKS_FILE);
		}
				
		fprintf(fp,"\n\n");
	}

		fprintf(fp,"\n");
		fprintf(fp," --------- Eikonal solver convergence criterion ----------------------\n");
		fprintf(fp," Convergence criterion for fast sweeping technique\n");
   		fprintf(fp," TTNORM = %e",TTNORM);		
		fprintf(fp,"\n\n");

        if(INVMAT){
  	        fprintf(fp,"\n");
                fprintf(fp," -----------------------  RAJZEL inversion parameters  ----------------------\n\n");
		fprintf(fp,"\n Maximum number of iterations: %d\n",ITERMAX);
		fprintf(fp,"\n Water level in sweep_adj.c EPS_ADJ: %e\n",EPS_ADJ);
		fprintf(fp,"\n location of measured FA picks : \n ");
		fprintf(fp,"\t%s\n\n",DATA_DIR);
		fprintf(fp,"\n INVTYPE = %d\n\n",INVTYPE);
	
		fprintf(fp," Taper used : \n ");
		fprintf(fp,"\t%d\n",TAPER); 
	
		fprintf(fp," Log file for misfit in each iteration step: \n ");
		fprintf(fp,"\t%s \n\n",MISFIT_LOG_FILE); 
	
		fprintf(fp," Output of inverted models: \n ");
		fprintf(fp,"\t%s (nfstart=%d, nf=%d)\n\n",INV_MODELFILE,nfstart,nf);
	
		fprintf(fp," Output of jacobian: \n ");
		fprintf(fp,"\t%s (nfstart_jac=%d, nf_jac=%d)\n\n\n",JACOBIAN,nfstart_jac,nf_jac);
		
		fprintf(fp," --------------- Gradient tapering -------------------\n");
		if (SWS_TAPER_GRAD_VERT==1){
			fprintf(fp," SWS_TAPER_GRAD_VERT=%d: Vertical taper applied.\n",SWS_TAPER_GRAD_VERT);
			fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d)\n\n",GRADT1,GRADT2,GRADT3,GRADT4);} 
		else	fprintf(fp," SWS_TAPER_GRAD_VERT=%d: No vertical taper applied.\n\n",SWS_TAPER_GRAD_VERT);
	
		if (SWS_TAPER_GRAD_HOR==1){
			fprintf(fp," SWS_TAPER_GRAD_HOR=%d: Horizontal taper applied.\n",SWS_TAPER_GRAD_HOR);
			fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d, EXP_TAPER_GRAD_HOR=%f)\n\n",GRADT1,GRADT2,GRADT3,GRADT4,EXP_TAPER_GRAD_HOR);} 
		else	fprintf(fp," SWS_TAPER_GRAD_HOR=%d: No horizontal taper applied.\n\n",SWS_TAPER_GRAD_HOR);
	
		if (SWS_TAPER_GRAD_SOURCES==1){
			fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: Taper around the sources.\n",SWS_TAPER_GRAD_SOURCES);
			fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);} 
		else	fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: No taper around the sources applied.\n\n",SWS_TAPER_GRAD_SOURCES);
	
		if (SWS_TAPER_CIRCULAR_PER_SHOT==1){
			fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: Taper around the source for each shot.\n",SWS_TAPER_CIRCULAR_PER_SHOT);
			fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);} 
		else	fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: No taper around the sources applied.\n\n",SWS_TAPER_CIRCULAR_PER_SHOT);
		
		fprintf(fp,"\n\n");
		fprintf(fp," --------------- Limits of model parameters -------------------\n");
		fprintf(fp," VPLOWERLIM = %f \t\t VPUPPERLIM = %f \n",VPLOWERLIM,VPUPPERLIM);
	
	
		fprintf(fp,"\n\n");
		fprintf(fp," --------------- Optimization method -------------------\n");
		switch(GRAD_METHOD){
			case 1:
				fprintf(fp," GRAD_METHOD=%d: PCG\n",GRAD_METHOD);
				break;
			case 2:
				fprintf(fp," GRAD_METHOD=%d: LBFGS\n",GRAD_METHOD);
		                fprintf(fp," NLBFGS=%d \n",NLBFGS);
				break;
			case 3:
				fprintf(fp," GRAD_METHOD=%d: \n",GRAD_METHOD);
				break;
			default:
				err(" Sorry, incorrect value for GRAD_METHOD ! \n");
		}
	
	
		fprintf(fp,"\n\n");
		fprintf(fp," --------------- Model smoothing -------------------\n");
		if(MODEL_FILTER==1){
			fprintf(fp," MODEL_FILTER=%d: vp and vs models are filtered after each iteration step.\n",MODEL_FILTER);
			fprintf(fp," (FILT_SIZE=%d)\n",FILT_SIZE);}
		else 	fprintf(fp," MODEL_FILTER=%d: vp and vs models are not filtered after each iteration step.\n",MODEL_FILTER);
	
	
		fprintf(fp,"\n\n");
		fprintf(fp," --------------- Step length estimation -------------------\n");
		fprintf(fp," EPS_SCALE = %f\n",EPS_SCALE);
		fprintf(fp," Wolfe condition parameters C1 = %f, C2 = %f\n",C1, C2);
		fprintf(fp," STEPMAX = %d\n",STEPMAX);
		fprintf(fp," SCALEFAC = %f\n",SCALEFAC);

		fprintf(fp,"\n\n");
		fprintf(fp," --------------- 1D linear gradient grid search parameters -------------------\n");
		fprintf(fp," VP0_1 = %f \t VP0_2 = %f, DVP0 = %f\n",VP0_1,VP0_2,DVP0);
		fprintf(fp," GRAD0_1 = %f \t GRAD0_2 = %f, DGRAD0 = %f\n",GRAD0_1,GRAD0_2,DGRAD0);
		fprintf(fp," writing results of grid search to file \n");
		fprintf(fp,"\t%s\n",GRIDSEARCH_FILE);
			                                                
		fprintf(fp,"\n");
		fprintf(fp," **************************************************************\n");
		fprintf(fp,"\n");

        }
}
