/*------------------------------------------------------------------------
 *  First Arrival traveltime tomography based on the adjoint state method 
 *
 *  D. Koehn
 *  Kiel, 27.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void fatt(float ** Vp, float ** S, float ** TT, float * Tmod, float *Tobs, float *Tres, float ** srcpos, int nshots, int ** recpos, int ntr, char *fileinp1){

	/* declaration of global variables */
        extern int NX, NY, NSHOT1, NSHOT2, GRAD_METHOD, NLBFGS, MYID, ITERMAX, LINESEARCH;
        extern char MISFIT_LOG_FILE[STRING_SIZE];
        extern float PRO;
    
        /* declaration of local variables */
        int ishot, i, j;
        char ext[10];

        /* variables for gradient estimation */
	float  ** lam, ** Vpnp1, ** grad, ** Hgrad, ** gradm;

        /* variables for the L-BFGS method */
	float * rho_LBFGS, * alpha_LBFGS, * beta_LBFGS; 
	float * y_LBFGS, * s_LBFGS, * q_LBFGS, * r_LBFGS;
	int NLBFGS_class, LBFGS_pointer, NLBFGS_vec;

	/*vector for abort criterion*/
	double * L2_hist=NULL;

	/* variables for workflow */
	int nstage, stagemax, iter, iter_true;
        
        /* variables for step-length estimation */
        float eps_scale; 
	double L2, *L2t, diff;
		

        FILE *FP, *FP_stage, *FPL2;

	/* read parameters from workflow-file (stdin) */
	FP=fopen(fileinp1,"r");
	if(FP==NULL) {
		if (MYID == 0){
			printf("\n==================================================================\n");
			printf(" Cannot open RAJZEL workflow input file %s \n",fileinp1);
			printf("\n==================================================================\n\n");
			err(" --- ");
		}
	}

	/* estimate number of lines in FWI-workflow */
	i=0;
	stagemax=0;
	while ((i=fgetc(FP)) != EOF)
	if (i=='\n') ++stagemax;
	rewind(FP);
	stagemax--;
	fclose(FP);


	/* memory allocation for abort criterion*/
	L2_hist = dvector(1,ITERMAX*stagemax);


	/* Variables for the l-BFGS method */
	if(GRAD_METHOD==2){

	  NLBFGS_class = 1;                 /* number of parameter classes */ 
	  NLBFGS_vec = NLBFGS_class*NX*NY;  /* length of one LBFGS-parameter class */
	  LBFGS_pointer = 1;                /* initiate pointer in the cyclic LBFGS-vectors */
	  
	  y_LBFGS  =  vector(1,NLBFGS_vec*NLBFGS);
	  s_LBFGS  =  vector(1,NLBFGS_vec*NLBFGS);

	  q_LBFGS  =  vector(1,NLBFGS_vec);
	  r_LBFGS  =  vector(1,NLBFGS_vec);

	  rho_LBFGS = vector(1,NLBFGS);
	  alpha_LBFGS = vector(1,NLBFGS);
	  beta_LBFGS = vector(1,NLBFGS);
	  
	}

	/* memory of L2 norm */
	L2t = dvector(1,4);

        /* memory for gradient */
   	lam =  matrix(1,NY,1,NX);
   	grad = matrix(1,NY,1,NX);
   	gradm = matrix(1,NY,1,NX);
   	Hgrad = matrix(1,NY,1,NX);
   	Vpnp1 =  matrix(1,NY,1,NX);

	iter_true=1;
	/* Begin of FATT inversion workflow */
	for(nstage=1;nstage<=stagemax;nstage++){

		/* read workflow input file *.inp */
		FP_stage=fopen(fileinp1,"r");
		read_par_inv(FP_stage,nstage,stagemax);

		iter=1;
		/* --------------------------------------
		 * Begin of FATT iteration loop
		 * -------------------------------------- */
		while(iter<=ITERMAX){

			if(GRAD_METHOD==2){
			  
			  /* increase pointer to LBFGS-vector*/
			  if(iter>2){
			    LBFGS_pointer++;
			  }
			  
			  /* if LBFGS-pointer > NLBFGS -> set LBFGS_pointer=1 */ 
			  if(LBFGS_pointer>NLBFGS){LBFGS_pointer=1;}

			}


			if (MYID==0){
			   printf("\n\n\n ------------------------------------------------------------------\n");
			   printf("\n\n\n                   FATT ITERATION %d \t of %d \n",iter,ITERMAX);
			   printf("\n\n\n ------------------------------------------------------------------\n");
			}

			/* Open Log File for L2 norm */
			  if(MYID==0){
			    if(iter_true==1){
			      FPL2=fopen(MISFIT_LOG_FILE,"w");
			    }

			    if(iter_true>1){
			      FPL2=fopen(MISFIT_LOG_FILE,"a");
			    }
			  }

			   /* First-Arrival Traveltime Tomography (FATT) */
			   /* ------------------------------------------ */

			   /* calculate Vp gradient and objective function */
			   L2 = grad_obj(grad,S,TT,lam,Tmod,Tobs,Tres,srcpos,nshots,recpos,ntr,iter);
			   L2t[1] = L2;

			   /* calculate descent directon gradm from gradient grad */
			   descent(grad,gradm);

			   /* estimate search direction waveconv with ... */
			   /* ... non-linear preconditioned conjugate gradient method */
			   if((GRAD_METHOD==1)||(GRAD_METHOD==3)){
                              MPI_Barrier(MPI_COMM_WORLD);
			      PCG(Hgrad,gradm,iter);
			   }

			   /* ... quasi-Newton l-BFGS method */
			   if(GRAD_METHOD==2){
                              MPI_Barrier(MPI_COMM_WORLD);
			      LBFGS(Hgrad,grad,gradm,iter,y_LBFGS,s_LBFGS,rho_LBFGS,alpha_LBFGS,Vp,q_LBFGS,r_LBFGS,beta_LBFGS,LBFGS_pointer,NLBFGS,NLBFGS_vec);
			   }

                           /* ... Descent method */
                           if(GRAD_METHOD==3){
                              MPI_Barrier(MPI_COMM_WORLD);
                              descent(grad,Hgrad);
                           }


			   /* check if search direction is a descent direction, otherwise reset l-BFGS history */
			   check_descent(Hgrad,grad,NLBFGS_vec,y_LBFGS,s_LBFGS,iter);

			   /* Estimate optimum step length ... */
			   MPI_Barrier(MPI_COMM_WORLD);

			   /* ... by line search which satisfies the Wolfe conditions */
                           if(LINESEARCH==1){
			   eps_scale = wolfels(Hgrad,grad,Vp,S,TT,lam,Tmod,Tobs,Tres,srcpos,nshots,recpos,ntr,iter,eps_scale,L2);}
			   
			   /* ... by inexact parabolic line search */
                           if(LINESEARCH==2){
			   eps_scale = parabolicls(Hgrad,grad,Vp,S,TT,lam,Tmod,Tobs,Tres,srcpos,nshots,recpos,ntr,iter,eps_scale,L2);}
			   

			   if(MYID==0){
			      fprintf(FPL2,"%e \t %e \t %d \n",eps_scale,L2t[1],nstage);
			   }

			   /* saving history of final L2*/
			   L2_hist[iter]=L2t[1];

			   /* calculate optimal change in the material parameters */
                           MPI_Barrier(MPI_COMM_WORLD);
			   calc_mat_change_wolfe(Hgrad,Vp,Vpnp1,eps_scale,0);
			   calc_S(Vp,S);

			    if(MYID==0){
			       fclose(FPL2);
			    }

			    /* calculating difference of the actual L2 and before two iterations, dividing with L2_hist[iter-2] provide changing in percent*/
			    diff=1e20;
			    if(iter > 2){
			       diff=fabs((L2_hist[iter-2]-L2_hist[iter])/L2_hist[iter-2]);
			    }

			    /* are convergence criteria satisfied? */	
			    if((diff<=PRO)||(eps_scale<1e-10)){
	
			       /* model output at the end of given workflow stage */
			       model_out(Vp,nstage);
			       iter=0;

			       if(GRAD_METHOD==2){

				  zero_LBFGS(NLBFGS, NLBFGS_vec, y_LBFGS, s_LBFGS, q_LBFGS, r_LBFGS, alpha_LBFGS, beta_LBFGS, rho_LBFGS);
				  LBFGS_pointer = 1;  

			       }

			       if(MYID==0){

				  if(eps_scale<1e-10){

				     printf("\n Steplength estimation failed. Changing to next FATT stage \n");

				  }else{

				    printf("\n Reached the abort criterion of pro=%e and diff=%e \n Changing to next FATT stage \n",PRO,diff);
			
				  }

			       }
		
			       break;
			    }

		iter++;
		iter_true++;

		/* ====================================== */
		} /* end of FATT iteration loop*/
		/* ====================================== */

	} /* End of FATT-workflow loop */
        
        /* memory deallocation */

        /* free memory for abort criterion */
        free_dvector(L2_hist,1,ITERMAX*stagemax);
        free_dvector(L2t,1,4);

        /* free memory for gradient */
        free_matrix(lam,1,NY,1,NX);
        free_matrix(grad,1,NY,1,NX);
        free_matrix(gradm,1,NY,1,NX);
   	free_matrix(Hgrad,1,NY,1,NX);
   	free_matrix(Vpnp1,1,NY,1,NX);

	/* free memory for l-BFGS */
	if(GRAD_METHOD==2){
	  
	  free_vector(y_LBFGS,1,NLBFGS_vec*NLBFGS);
	  free_vector(s_LBFGS,1,NLBFGS_vec*NLBFGS);

	  free_vector(q_LBFGS,1,NLBFGS_vec);
	  free_vector(r_LBFGS,1,NLBFGS_vec);

	  free_vector(rho_LBFGS,1,NLBFGS);
	  free_vector(alpha_LBFGS,1,NLBFGS);
	  free_vector(beta_LBFGS,1,NLBFGS);
	  
	}
        	    
}
