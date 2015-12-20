/*------------------------------------------------------------------------
 *  Solve adjoint eikonal equation by fast sweeping method according to Leung & Qian (2006)
 *
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void sweep_adj(float ** lam, float ** TT, float * Tres, int ** recpos, int ntr, int nx1, int nx2, int ndx, int ny1, int ny2, int ndy){

	extern int NX, NY;
	extern float DH, EPS_ADJ;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

	/* local variables */
	int i, j, h, k, l;
	float app, amp, apm, amm;
        float bpp, bmp, bpm, bmm;
        float ap, am , bp, bm;
        float lhs, rhs, lamt;

	/* sweep over FD-grid */
        h = nx1;
	for (i=2;i<=NX-1;i++){
                k = ny1;
		for (j=2;j<=NY-1;j++){

                        /* assemble equation (3.6) in Leung & Qian (2006) */
                        ap = -(TT[k][h+1]-TT[k][h])/DH;
                        am = -(TT[k][h]-TT[k][h-1])/DH;
                        
			bp = -(TT[k+1][h]-TT[k][h])/DH;
                        bm = -(TT[k][h]-TT[k-1][h])/DH;

			app = (ap + fabs(ap))/DH;
			apm = (ap - fabs(ap))/DH;

			amp = (am + fabs(am))/DH;
			amm = (am - fabs(am))/DH;

                        bpp = (bp + fabs(bp))/DH;
			bpm = (bp - fabs(bp))/DH;

			bmp = (bm + fabs(bm))/DH;
			bmm = (bm - fabs(bm))/DH;

                        lhs = (app-amm)/DH + (bpp-bmm)/DH;
                        rhs = (amp*lam[k][h-1]-apm*lam[k][h+1])/DH + (bmp*lam[k-1][h]-bpm*lam[k+1][h])/DH;
                        
                        if(lhs > EPS_ADJ){
                           lamt = rhs/lhs;
                        }else{
                           lamt=0.0;
                        }

                        /* add data residuals at receiver positions */
                        for (l=1;l<=ntr;l++){
                           
                           if((h==recpos[1][l])&&(k==recpos[2][l])){
			      lamt += Tres[l];
                           }

                        }

                        lam[k][h] = fminf(lam[k][h],lamt);

		        k += ndy;
		}
        h += ndx;
	}

}




