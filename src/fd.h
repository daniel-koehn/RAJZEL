/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD programs          
 *  last update  03/12/2000 
 *
 *  Copyright (c) 1998 T. Bohlen 
 *  See COPYING file for copying and redistribution conditions.
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 96
#define STRING_SIZE 74
#define STRING_SIZE2 256
#define REQUEST_COUNT 4


/* declaration of functions */

void ass_grad(float ** waveconv, float ** waveconv_shot, float ** S, float ** lam, float **srcpos,  int nshots, int **recpos, int ntr, int ishot);

void calc_FA(float ** TT, int ** recpos, int ntr, float * Tmod);

float calc_FA_res(float * Tmod, float * Tobs, float * Tres, int ntr, int ishot);

float calc_mat_change(float  **  waveconv, float **  pi, float **  pinp1, int iter, float eps_scale, int itest, int nfstart);

void calc_mat_change_wolfe(float  **  waveconv, float **  pi, float **  pinp1, float eps_scale, int itest);

float calc_opt_step(float *  L2t, float * epst, int sws);

void calc_S(float **  Vp, float **  S);

void check_descent(float ** waveconv, float ** gradp, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, int iter);

void cp_grad_frame(float ** A);

void descent(float ** grad, float ** gradm);

float dotp(float * vec1, float *vec2, int n1, int n2, int sw);

float dotp_matrix(float ** A, float ** B, int NX, int NY);

void eikonal(float ** S, float ** TT, float * Tmod, float ** srcpos, int ishot, int nshots, int **recpos, int ntr);

void eikonal_adj(float ** TT, float * Tres, float ** lam, int ** recpos, int ntr, int ishot, int **recflag);

void exchange_grad_MPI(float ** grad);

void fatt(float ** Vp, float ** S, float ** TT, float * Tmod, float *Tobs, float *Tres, float ** srcpos, int nshots, int ** recpos, int ntr, char *fileinp1);

void forward(float ** S, float ** TT, float * Tmod,  float ** srcpos, int nshots, int ** recpos, int ntr);

void gauss_filt(float ** waveconv);

float grad_obj(float ** grad, float ** S, float ** TT, float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr, int iter);

void grid_search(float ** Vp, float ** S, float ** TT, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr);

void info_mem(FILE *fp, int NLBFGS_vec, int ntr);

void info(FILE *fp);

void init_grad(float ** A);

void init_MPIshot(int nshots);

void init_recflag(int **recflag, int ** recpos, int ntr);

void init_travel_adj(float **lam, int ** recpos, int ntr, float * Tres, float **TT);

void LBFGS(float ** Hgrad, float ** grad, float ** gradm, int iter, float * y_LBFGS, float * s_LBFGS, float * rho_LBFGS, float * alpha_LBFGS, 
float **Vp, float * q_LBFGS, float * r_LBFGS, float * beta_LBFGS, int LBFGS_pointer, int NLBFGS, int NLBFGS_vec);
           
double LU_decomp(double  **A, double *x, double *b,int n);

float minimum_m(float **mat, int nx, int ny);
float maximum_m(float **mat, int nx, int ny);

float median2d(float **mat, int ny, int nx);

void median_src(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int iter, int sws);

void model(float  ** Vp);

void model_gridsearch(float  ** Vp, float vp0, float grad0);

void model_out(float ** Vp, int iter);		  

float norm(float ** waveconv, int iter, int sws);

float norm1(float ** TT, float ** TTold);

float norm_matrix(float **A, int NX, int NY);

void note(FILE *fp);

float parabolicls(float ** Hgrad, float ** grad, float ** Vp, float ** S, float ** TT, float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr, 
                  int iter, float alpha, float L2);

void PCG(float ** Hgrad, float ** grad, int iter);

void precond(float ** grad, int nsrc, float ** srcpos, int ** recpos, int ntr, int iter);

float readdsk(FILE *fp_in, int format);

void read_par(FILE *fp_in);

void read_par_inv(FILE *fp,int nstage,int stagemax);

void readmod(float  **  Vp);

int **receiver(FILE *fp, int *ntr, int ishot);

void smooth_grad(float ** waveconv);

void  smooth2(float ** grad);

void smooth_model(float ** pinp1, float ** unp1, float ** rho, int iter);

float **sources(int *nsrc);

void solvelin(float  **AA, float *bb, float *x, int e, int method);

void store_mat(float ** A, float ** B, int n, int m);

void sum_grad_MPI(float ** grad);

void sweep(float ** S, float ** TT, int nxsrc, int nysrc, int nx1, int nx2, int ndx, int ny1, int ny2, int ndy);

void sweep_adj(float ** lam, float ** TT, float * Tres, int ** recpos, int ntr, int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int **recflag);

void taper_grad(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int sws);

void taper_grad_shot(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void taper_grad_shot1(float ** waveconv, float **srcpos, int nshots, int **recpos, int ntr, int ishot);

void tripd(float *d, float *e, float *b, int n);

void  wavenumber(float ** grad);

float wolfels(float ** waveconv, float ** gradp, float ** Vp, float ** S, float ** TT, float ** lam, float * Tmod, float * Tobs, float * Tres,  float ** srcpos, int nshots, int ** recpos, int ntr,
              int iter, float alpha, float L2);

void write_gridsearch(float L2, float vp0, float grad0, int count);

void write_par(FILE *fp);

void write_picks(float * Tmod, int ntr, int ishot);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float ** array, int format);

void zero_LBFGS(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS, float * q_LBFGS, float * r_LBFGS, 
                 float * alpha_LBFGS, float * beta_LBFGS, float * rho_LBFGS);

void zero_LBFGS1(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS);

void zero_travel(float **TT);

/* utility functions */
void err(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh); 
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, 
int ndh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);
