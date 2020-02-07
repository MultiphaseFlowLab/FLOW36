MF	= Makefile

#define shell
SHELL = /bin/sh

#FC	= mpifort
#OLEVEL	= -g #-C to check bounds
#CPPF = -cpp
#FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -J ./../source_code #-std=f2008ts
#FOPTS	= -align none -mcmodel medium -warn all
#FFLAGS	= $(CPPF) $(OLEVEL)

FLIBS    =  -I/usr/include/ -I/usr/local/include -L/usr/local/lib   -lm -ldl #-lfftw3  #
#LIBS    =  -I$(FFTW_INC) -L$(FFTW_LIB) -lm -ldl -lfftw3
FLIB     =   -lfftw3###library at the end: single pass linking(?)this is not a HPC runtime valid solution


##CUDA
CUDA_PATH = /usr/local/cuda
MPI_PATH  = /usr/
#CUDA_FC = /mpi/openmpi-3.1.3/bin/mpifort
#ARC = -Mcuda=cc60
#CUDA_FLAGS = $(OLEVEL) $(ARC)  


VPATH	= ./set_run/sc_compiled/

FLWOBJS =  $(shell cat list_source.list)
#define Fortran objects list
F_OBJS	 = \
		module.o \
		interfaccia.o \
		main.o \
		read_input.o \
		print_start.o \
		define_sizes.o \
		create_plan.o \
		destroy_plan.o \
		initialize.o \
		dump_grid.o \
		write_output.o \
		write_time.o \
		solver.o \
		convective_ns.o \
		wave_numbers.o \
		phys_to_spectral.o \
		spectral_to_phys.o \
		fftx_fwd.o \
		fftx_bwd.o \
		ffty_fwd.o \
		ffty_bwd.o \
		yz2xz.o \
		xz2xy.o \
		xz2yz.o \
		xy2xz.o \
		dctz_fwd.o \
		dctz_bwd.o \
		dz.o \
		time_integration.o \
		hist_term.o \
		helmholtz.o \
		helmholtz_red.o \
		calculate_var.o \
		courant_check.o \
		read_fields.o \
		initialize_phi.o \
		phi_non_linear.o \
		sterm_ch.o \
		statistics.o \
		sim_check.o \
		chop_modes.o \
		initialize_psi.o \
		sterm_surf.o \
		initialize_temp.o \
		sterm_temp.o \
		hist_term_temp.o \
		define_address.o \
		define_dual_grid.o \
		swap_grid.o \
		shrink.o \
		split_comm.o \
		part_fluid_comm.o \
		lagrangian_interpolator.o \
		lagrangian_tracker.o \
		save_flow_comm.o 
		
#define cuda objects list
CUDA_OBJS = cuda_spec_phys.o \
			cuda_tran.o \
			cuda_phys_spec.o \
			init_gpu.o \
			free_gpu.o


#initialize_particle.f90 \

alld:
	make	EXE "PROG = surf_gpu.exe" \
		"GCC = gcc" \
		"GCCFLAGS = -fno-exceptions" \
		"CPP = cpp" \
		"CPPFLAGS = -cpp -D=LINUX" \
		"F90 = mpifort" \
		"F90FLAGS = -O3" \
		"NVCC = nvcc" \
		"NVCCFLAGS = -O3 -use_fast_math -arch=sm_61 --ptxas-options=-v" \
		"INCD = -I"$(CUDA_PATH)"/include " \
		"LIBS = -L"$(CUDA_PATH)"/lib64 -lcudart -lcurand -lcufft -lstdc++ -lfftw3"	
all:
	make	EXE "PROG = surf_gpu.exe" \
		"GCC = gcc" \
		"GCCFLAGS = -fno-exceptions" \
		"CPP = cpp" \
		"CPPFLAGS = -cpp -D=LINUX" \
		"F90 = mpifort" \
		"F90FLAGS = -g" \
		"NVCC = nvcc" \
		"NVCCFLAGS = -g -G -arch=sm_61 --ptxas-options=-v" \
		"INCD = -I"$(CUDA_PATH)"/include -I/usr/include/mpi " \
		"LIBS = -L"$(CUDA_PATH)"/lib64 -lcudart -lcurand -lcufft -lstdc++ -lfftw3"		

.SUFFIXES: .cu .h .c .F90 .f90 .o
.c.o:
	$(GCC) $(GCCFLAGS) -c $< $*.o
.cu.o:
	$(NVCC) $(NVCCFLAGS) $(INCD) -c $< $*.o
.f90.o:
	$(F90) $(CPPFLAGS) $(F90FLAGS) $(FLIBS) -c $<
	
EXE: $(C_OBJS) $(F_OBJS) $(CUDA_OBJS)
	$(F90) -o $(PROG) $(C_OBJS) $(F_OBJS) $(CUDA_OBJS) $(LIBS)

%.o: %.mod

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f *.o *.mod *.exe

#objects rule for iso c binding
#cuda_modules.o:			cuda_modules.f90
#cuda_pass_arrays.o:		cuda_pass_arrays.f90		cuda_modules.o
#fortran objects rule
interfaccia.o:			interfaccia.f90			init_gpu.o		free_gpu.o		cuda_spec_phys.o	cuda_phys_spec.o
#main.o:					main.f90				interfaccia.o
#cuda objects rule
#cuda_surf.o:		cuda_surf.cu	cuda_variables.h	cuda_surf.h		cuda_spec_phys.h	cuda_tran.h	cuda_phys_spec.h
cuda_spec_phys.o:	cuda_spec_phys.cu	cuda_variables.h	cuda_spec_phys.h	cuda_tran.h
cuda_tran.o:		cuda_tran.cu	cuda_variables.h	cuda_tran.h
cuda_phys_spec.o:	cuda_phys_spec.cu	cuda_variables.h	cuda_phys_spec.h	cuda_tran.h
init_gpu.o:				init_gpu.cu				cuda_variables.h		init_gpu.h
free_gpu.o:				free_gpu.cu				cuda_variables.h		free_gpu.h
