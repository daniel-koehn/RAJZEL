/*------------------------------------------------------------------------
 *  Solve eikonal equation by fast sweeping method according to Zhao (2004)
 *
 *  D. Koehn
 *  Kiel, 09/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void sweep(float ** S, float ** TT, int nxsrc, int nysrc, int nx1, int nx2, int ndx, int ny1, int ny2, int ndy){

	extern int NX, NY;
	extern float DH;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

	
	/* local variables */
	int i, j, h, k;
	float a, b, Tt;

   
	/* sweep over FD-grid */
        h = nx1;
	for (i=1;i<=NX;i++){
                k = ny1;
		for (j=1;j<=NY;j++){	

			/* TT[nysrc][nxsrc] = 0.0; */

			/* model interior */
			if((h>1)&&(h<NX)&&(k>1)&&(k<NY)){
			    a = fminf(TT[k][h-1],TT[k][h+1]);
			    b = fminf(TT[k-1][h],TT[k+1][h]);
			}
		
			/* model borders */
                        /* left */
			if((h==1)&&(k>1)&&(k<NY)){
			    a = fminf(TT[k][h],TT[k][h+1]);
			    b = fminf(TT[k-1][h],TT[k+1][h]);
			}
		
                        /* right */
			if((h==NX)&&(k>1)&&(k<NY)){
			    a = fminf(TT[k][h-1],TT[k][h]);
			    b = fminf(TT[k-1][h],TT[k+1][h]);
			}
		
                        /* top */
			if((k==1)&&(h>1)&&(h<NX)){
			    a = fminf(TT[k][h-1],TT[k][h+1]);
			    b = fminf(TT[k][h],TT[k+1][h]);
			}
		
                        /* bottom */
			if((k==NY)&&(h>1)&&(h<NX)){
			    a = fminf(TT[k][h-1],TT[k][h+1]);
			    b = fminf(TT[k-1][h],TT[k][h]);
			}
		
			/* model corners */
                        /* upper-left */
			if((h==1)&&(k==1)){
			    a = fminf(TT[k][h],TT[k][h+1]);
			    b = fminf(TT[k][h],TT[k+1][h]);
			}
		
                        /* lower-left */
			if((h==1)&&(k==NY)){
			    a = fminf(TT[k][h],TT[k][h+1]);
			    b = fminf(TT[k-1][h],TT[k][h]);
			}

		        /* lower-right */
			if((h==NX)&&(k==NY)){
			    a = fminf(TT[k][h-1],TT[k][h]);
			    b = fminf(TT[k-1][h],TT[k][h]);
			}
		
                        /* upper-right */
			if((h==NX)&&(k==1)){
			    a = fminf(TT[k][h-1],TT[k][h]);
			    b = fminf(TT[k][h],TT[k+1][h]);
			}
		

			/* calculate solution */
			if(fabs(a-b)>=(S[k][h]*DH)){
			    Tt = fminf(a,b) + S[k][h]*DH;
			}

			if(fabs(a-b)<(S[k][h]*DH)){
			    Tt = (a + b + sqrt((2.0*pow(S[k][h],2.0)*pow(DH,2.0)) - pow((a-b),2.0)))/2.0;
			}

			TT[k][h] = fminf(TT[k][h],Tt);

		        k += ndy;
		}
        h += ndx;
	}

}




