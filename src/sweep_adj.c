/*------------------------------------------------------------------------
 *  Solve adjoint eikonal equation by fast sweeping method according to Leung & Qian (2006)
 *
 *  D. Koehn
 *  Kiel, 11/12/2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void sweep_adj(float ** lam, float ** TT, float * Tres, int ** recpos, int ntr, int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int **recflag){

	extern int NX, NY;
	extern float DH, EPS_ADJ;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

	/* local variables */
	int i, j, h, k, l;
	float app, amp, apm, amm;
        float bpp, bmp, bpm, bmm;
        float ap, am , bp, bm;
        float lamt, lhs, rhs;

	/* sweep over FD-grid */
        h = nx1;
	for (i=2;i<=NX-1;i++){
                k = ny1;
		for (j=2;j<=NY-1;j++){

                     if(recflag[k][h]==0){			

                        /* assemble equation (3.6) in Leung & Qian (2006) */
                        ap = -(TT[k][h+1]-TT[k][h])/DH;
                        am = -(TT[k][h]-TT[k][h-1])/DH;
                        
			bp = -(TT[k+1][h]-TT[k][h])/DH;
                        bm = -(TT[k][h]-TT[k-1][h])/DH;

			app = (ap + fabs(ap))/2.0;
			apm = (ap - fabs(ap))/2.0;

			amp = (am + fabs(am))/2.0;
			amm = (am - fabs(am))/2.0;

                        bpp = (bp + fabs(bp))/2.0;
			bpm = (bp - fabs(bp))/2.0;

			bmp = (bm + fabs(bm))/2.0;
			bmm = (bm - fabs(bm))/2.0;

                        
                        /* Leung & Qian (2006) */
                        lhs = (app-amm)/DH + (bpp-bmm)/DH;
                        rhs = (amp*lam[k][h-1]-apm*lam[k][h+1])/DH + (bmp*lam[k-1][h]-bpm*lam[k+1][h])/DH;
                        
                        /* Taillandier et al. (2009) */
                        /*lhs = (apm-amp)/DH + (bpm-bmp)/DH;
                        rhs = (amm*lam[k][h-1]-app*lam[k][h+1])/DH + (bmm*lam[k-1][h]-bpp*lam[k+1][h])/DH;*/
			
			if(lhs>EPS_ADJ){
                            lamt = rhs/lhs;
			}else{
			    lamt = 0.0;
			}			

                        lam[k][h] = fminf(lam[k][h],lamt);
			
			if(TT[k][h]>=INF || TT[k+1][h]>=INF || TT[k-1][h]>=INF || TT[k][h+1]>=INF || TT[k][h-1]>=INF){lam[k][h] = 0.0;}                       

                     }

		        k += ndy;
		}
        h += ndx;
	}

}
