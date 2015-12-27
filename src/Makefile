# Makefile for RAJZEL

#--------------------------------------------------------
# edit here:

# source code for model generation

MODEL = 1D_start.c
MODEL_GRID_SEARCH = 1D_start_gridsearch.c
EXEC= ../bin

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4 

# On CRAY T3E
# CC=cc

# On Linux NEC-cluster with Intel-MPI
#CC=mpiicc 
#LFLAGS=-lm -lfftw3 -ftz -lstdc++
#CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
#SFLAGS=-L$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/lib
#IFLAGS=-I$//sfs/fs4/work-sh2/sungw331/fftw-3.3.4/include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
CC=mpicc
LFLAGS=-lm -lfftw3 -lstdc++
CFLAGS=-O3 -fno-stack-protector

# On HLRN system
#CC=mpcc
#LFLAGS=-lm  

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ipo
#LFLAGS=-lmpi -lm -i-static

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $< $(IFLAGS)

FDVEPS_UTIL = \
	absorb.c \
	exchange_par.c \
	info.c \
	initproc.c \
	$(MODEL_EL) \
	note.c \
	outseis.c \
	readdsk.c \
	read_par.c \
	readmod.c \
    	receiver.c \
	sources.c \
	util.c \
	writemod.c \
	write_par.c \
	writedsk.c 

RAJZEL= \
	rajzel.c \
        ass_grad.c \
	calc_FA.c \
	calc_FA_res.c \
        calc_mat_change.c \
        calc_mat_change_wolfe.c \
        calc_opt_step.c \
	calc_opt_step_test.c \
        calc_S.c \
        check_descent.c \
        descent.c \
	dotp.c \
        dotp_matrix.c \
        eikonal.c \
        eikonal_adj.c \
	step_length_est.c \
	median.c \
        fatt.c \
        forward.c \
	gauss_filt.c \
        grad_obj.c \
        grid_search.c \
	info.c \
        info_mem.c \
        init_grad.c \
        init_MPIshot.c \
        init_travel_adj.c \
        LBFGS.c \
	smooth2.c \
	smooth_model.c \
	smooth_grad.c \
	$(MODEL) \
        $(MODEL_GRID_SEARCH) \
	median_src.c \
	model_out.c \
	note.c \
	norm.c \
        norm1.c \
        norm_matrix.c \
	PCG.c \
        precond.c \
	readdsk.c \
	read_par.c \
	read_par_inv.c \
	readmod.c \
    	receiver.c \
	sources.c \
	solvelin.c \
        sweep.c \
        sweep_adj.c \
        store_mat.c \
	taper_grad.c \
	taper_grad_shot.c \
        taper_grad_shot1.c \
	util.c \
	wavenumber.c \
        wolfels.c \
        write_gridsearch.c \
	writemod.c \
	write_par.c \
        write_picks.c \
	writedsk.c \
	zero_LBFGS.c \
        zero_LBFGS1.c \
        zero_travel.c

RAJZEL_OBJ = $(RAJZEL:%.c=%.o)

rajzel:		$(RAJZEL_OBJ)
	$(CC) $(SFLAGS) $(RAJZEL_OBJ) -o $(EXEC)/rajzel $(LFLAGS)

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
