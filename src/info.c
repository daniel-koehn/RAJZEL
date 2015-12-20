/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 13.12.2015
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ***********************************************************\n");
	fprintf(fp," RAJZEL: Parallel 2D First-Arrival Traveltime Modelling and \n");
        fprintf(fp," Tomography Code using an Eikonal Solver and adjoint method \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," Eikonal solver and FATT code                               \n");
	fprintf(fp," written by D. Koehn and D. De Nil                          \n");
	fprintf(fp," Institute of Geosciences, Kiel University, Germany         \n\n");
	fprintf(fp," See COPYING file for copying and redistribution conditions.\n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");

}
