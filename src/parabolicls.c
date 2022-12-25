/*------------------------------------------------------------------------
 * Step length estimation by parabolic line search
 *
 * D. Koehn
 * Kiel, 11.12.2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float parabolicls(float ** Hgrad, float ** grad, float ** Vp, float ** S, float ** TT, 
float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, 
int ** recpos, int ntr, int iter, float alpha, double L2){

extern int MIN_ITER,STEPMAX, NX, NY, MYID;
extern char JACOBIAN[STRING_SIZE];
extern float EPS_SCALE, SCALEFAC;

float opteps_vp, **Vpnp1, **Snp1;
int h, i, j, n, ishot;

/* Variables for step length calculation */
int step1, step2, step3, itests, iteste, stepmax, countstep;
int itest;
float scalefac, tmp, **gt, eps_scale;
float maxgrad, maxvp, *L2t, *epst1;

Vpnp1 =  matrix(1,NY,1,NX);
Snp1 =  matrix(1,NY,1,NX);
gt =  matrix(1,NY,1,NX);
L2t = vector(1,3);
epst1 = vector(1,3);

scalefac = SCALEFAC;  /* scale factor for the step length */
stepmax  = STEPMAX;   /* number of maximum misfit calculations/steplength 2/3*/ 

step1=0;
step2=0;

/* start with first guess for step length alpha */
maxgrad = maximum_m(Hgrad,NX,NY);
  maxvp = maximum_m(Vp,NX,NY);

if(iter==1){alpha = EPS_SCALE * fabs(maxvp/maxgrad);}

countstep=0; /* count number of forward calculations */
L2t[1] = L2;

itests=2;
iteste=2;

while((step2!=1)||(step1!=1)){

      for(itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

          /* test Vp-update */
	  calc_mat_change_wolfe(Hgrad,Vp,Vpnp1,alpha,1);

          /* calculate slowness S from test Vp model*/
          calc_S(Vpnp1,Snp1);
          L2t[itest] = grad_obj(gt,Snp1,TT,lam,Tmod,Tobs,Tres,srcpos,nshots,recpos,ntr,iter);
 
	     
      } /* end of L2 test */

      /* Did not found a step size which reduces the misfit function */
      if((step1==0)&&(L2t[1]<=L2t[2])){

        alpha = alpha/scalefac; 
        countstep++;

      }

      /* Found a step size with L2t[2] < L2t[3] */
      if((step1==1)&&(L2t[2]<L2t[3])){

        epst1[3]=alpha;
        step2=1;

      }

      /* Could not found a step size with L2t[2] < L2t[3]*/
      if((step1==1)&&(L2t[2]>=L2t[3])){

         epst1[3]=alpha;

         /* increase step length to find  a larger misfit function than L2t[2]*/
         alpha = alpha + (alpha/scalefac);
         countstep++;                       
      }         

      /* found a step size which reduces the misfit function */
      if((step1==0)&&(L2t[1]>L2t[2])){

         epst1[2]=alpha; 
         step1=1;
         iteste=3;
         itests=3;
         countstep=0;
	 
         /* find a second step length with a larger misfit function than L2t[2]*/
         alpha = alpha + (alpha/scalefac);

      }

      step3=0;

      if((step1==0)&&(countstep>stepmax)){

         if(MYID==0){printf(" Steplength estimation failed!");} 
         step3=1;
         break;

      }

      if((step1==1)&&(countstep>stepmax)){

         if(MYID==0){
            printf("Could not found a proper 3rd step length which brackets the minimum\n");}
            step1=1;
            step2=1;

      }

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e\n",iteste,itests,step1,step2,eps_scale,countstep,stepmax,scalefac,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

if(step1==1){ /* only find an optimal step length if step1==1 */

   /* calculate optimal step length epsilon for Vp*/
   if(MYID==0){
      printf("===================================================== \n");
      printf("    calculate optimal step length epsilon for Vp      \n");
      printf("===================================================== \n");
   }

   opteps_vp = calc_opt_step(L2t,epst1,1);
   eps_scale = opteps_vp;

}

free_matrix(Vpnp1,1,NY,1,NX);
free_matrix(Snp1,1,NY,1,NX);
free_matrix(gt,1,NY,1,NX);
free_vector(L2t,1,3);
free_vector(epst1,1,3);

return eps_scale;

}

