/*------------------------------------------------------------------------
 *   Define receiver flag 
 *  
 *  
 *   D. Koehn
 *   Kiel, 11/01/2016
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_recflag(int **recflag, int ** recpos, int ntr){

	register int i, j, l;
	extern int NX, NY;

        
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

                     recflag[j][i] = 0;

                     /* recflag[j][i] = receiver no. */
                     for (l=1;l<=ntr;l++){
                           
                           if((i==recpos[1][l])&&(j==recpos[2][l])){
			      recflag[j][i] = l;
                           }

                     }

		}
	}           
	
}
