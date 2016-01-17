/*------------------------------------------------------------------------
 *   Apply wavenumber domain filter to gradient
 *   
 *   Daniel Koehn
 *   Kiel, the 2nd of October 2012 
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include <fftw3.h>

void  wavenumber(float ** grad){

	/* declaration of extern variables */
        extern int NX, NY, NXG, NYG, IDX, IDY;
	extern int NPROCX, NPROCY, MYID;
	extern float WD_DAMP;
	extern char JACOBIAN[STRING_SIZE];
	
	/* declaration of local variables */
	int i,j, h, fa, fb, npadx, npady, itr;
	int tracl1, jj, ii, zeropad;
	float gradtmp;
	double npaddx, npaddy;
        float ** gradtmp1;
        
	char jac[STRING_SIZE];
	                    
	/* size of the zeropadding layer */
        zeropad=0;               
                                    
        /*npadx = (int)(pow(2.0, ceil(log((double)(NXG))/log(2.0))+2.0) );
        npady = (int)(pow(2.0, ceil(log((double)(NYG))/log(2.0))+2.0) );*/
        
        npadx = NXG + 2.0*zeropad;
        npady = NYG + 2.0*zeropad;
        
        /* define temporary gradient matrix */
        gradtmp1 = matrix(1,npady,1,npadx);
        
        /* printf("npadx = %d \t npady = %d \n",npadx,npady); */
        
	fftw_complex    *data, *fft_result, *ifft_result;
	fftw_plan       plan_forward, plan_backward;
                                        
	data        = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );
	fft_result  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );
	ifft_result = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * npadx * npady );

        /* damping coefficient */
        /*damp=5e-5;
        damp=6e-5;*/

        if(MYID==0){
	   printf("\n Spatial filter is applied to gradient (written by PE %d)\n",MYID); 
        }	

	/* define real and imaginary part of the gradient */
	for (i=1;i<=npadx;i++){
	   for (j=1;j<=npady;j++){
	        
	        if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>=1+zeropad)&&(j<=NYG+zeropad)){
	            gradtmp1[j][i] = grad[j][i];
		}
		else{
		    gradtmp1[j][i]=0.0;
		}
			
            }
	}
	
	/* Fill padding layer with non-zeros */
	/* top and bottom padding */
        /*for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){
            
                if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j<(1+zeropad))){
                  gradtmp1[j][i]=gradtmp1[1+zeropad][i];
                }
                
                if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>(NYG+zeropad))){
                  gradtmp1[j][i]=gradtmp1[NYG+zeropad][i];
                }
            
            }
        }*/
        
        /* left and right padding */
        /*for (j=1;j<=npady;j++){
            for (i=1;i<=npadx;i++){
                                        
                if(i<(1+zeropad)){
                  gradtmp1[j][i]=gradtmp1[j][1+zeropad];
                }
                                                                                                          
                if(i>(NXG+zeropad)){                
                  gradtmp1[j][i]=gradtmp1[j][NXG+zeropad];         
                }
                                                                                                                                                                         
            }
        }*/
                                                                                                                                                                                             
	
        /* FFT of the  gradient */
        h=0;
        for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){
	                                                                                                               
                data[h][0] = gradtmp1[j][i];
                data[h][1] = 0.0;
	                                                                                                                                                                                                                     
                h++;
	                                                                                                                                                                                                                                                     
            }     
        }
	
	plan_forward  = fftw_plan_dft_2d(npady, npadx, data, fft_result, FFTW_FORWARD, FFTW_ESTIMATE );
        plan_backward = fftw_plan_dft_2d(npady, npadx, fft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE );
	
        /* apply 2D-fft */
	fftw_execute( plan_forward );

	/* Apply Gaussian wavenumber damping */
        h=0;
        for (i=1;i<=npadx;i++){
            for (j=1;j<=npady;j++){ 
                
                if((i<=npadx/2)&&(j<=npady/2)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-1)*(j-1)));
                  fft_result[h][1] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-1)*(j-1)));
                }
                
                if((i>npadx/2)&&(j<=npady/2)){ 
                  fft_result[h][0] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-1)*(j-1)));
                  fft_result[h][1] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-1)*(j-1)));
                }
                
                if((i<=npadx/2)&&(j>npady/2)&&(j<=npady)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-npady)*(j-npady)));          
                  fft_result[h][1] *=exp(-WD_DAMP*((i-1)*(i-1)+(j-npady)*(j-npady))); 
                }
                
                if((i>npadx/2)&&(j>npady/2)&&(j<=npady)){
                  fft_result[h][0] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-npady)*(j-npady)));                            
                  fft_result[h][1] *=exp(-WD_DAMP*((i-npadx)*(i-npadx)+(j-npady)*(j-npady)));                  
                }
                                           
            h++;		   
	    }	   
        }
		
	/* apply 2D-ifft */
	fftw_execute( plan_backward );

	/* output of ifft */
	h=0;
	for (i=1;i<=npadx;i++){
		for (j=1;j<=npady;j++){
			 
			 if((i>=1+zeropad)&&(i<=NXG+zeropad)&&(j>=1+zeropad)&&(j<=NYG+zeropad)){
			    grad[j][i] = ifft_result[h][0];
			 }
			 
			 h++;
		}
	}
	        
	/* free memory */
	fftw_free( data );
	fftw_free( fft_result );
	fftw_free( ifft_result );

        fftw_destroy_plan( plan_forward );
        fftw_destroy_plan( plan_backward );

        free_matrix(gradtmp1,1,npady,1,npadx);
                                

}
