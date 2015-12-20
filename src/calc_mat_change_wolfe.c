/*------------------------------------------------------------------------
 *   update material parameter for Wolfe line search 
 *   
 *   Daniel Koehn
 *   last update 14.12.2015
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void calc_mat_change_wolfe(float  **  Hgrad, float **  Vp, float **  Vpnp1, float eps_scale, int itest)
{

	/* global variables */
	extern int NX, NY;
	extern char INV_MODELFILE[STRING_SIZE];
	extern float VPUPPERLIM, VPLOWERLIM;

	/* local variables */
	int i, j;
	char modfile[STRING_SIZE];
        float maxgrad, maxvp;

        /* constant step length for debugging */
        /*maxgrad = maximum_m(Hgrad,NX,NY);
        maxvp = maximum_m(Vp,NX,NY);
        eps_scale = 0.01 * maxvp/maxgrad;*/	

	/* loop over local grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){		
		    
		    /* P-wave velocity */
		    Vpnp1[j][i] = Vp[j][i] + eps_scale * Hgrad[j][i]; 	
		  
		    /* apply hard constraints */
	      	    if(Vpnp1[j][i]<VPLOWERLIM){
	               Vpnp1[j][i] = Vp[j][i];
	            }
		      
		    if(Vpnp1[j][i]>VPUPPERLIM){
		       Vpnp1[j][i] = Vp[j][i];
		    }
		      
		      
		    /* P-wave velocity should not be smaller than zero */
		    if(Vpnp1[j][i]<0.0){
		       Vpnp1[j][i] = Vp[j][i];
		    }
  
		    if(itest==0){
		       Vp[j][i] = Vpnp1[j][i]; 
	            } 
		                
		}
	}
	
	if(itest==0){
	   sprintf(modfile,"%s_vp.bin",INV_MODELFILE);
	   writemod(modfile,Vpnp1,3);
	}
                                                
}

