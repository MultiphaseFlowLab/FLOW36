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


VPATH	= ./set_run/sc_compiled/


#define Fortran objects list
F_OBJS =  $(shell cat ../../GPU_local/f90_list_gpu.list)

#define CUDA objects list
CUDA_OBJS =  $(shell cat ../../GPU_local/cuda_list_gpu.list)



all:
	make	EXE "PROG = surf36_gpu.exe" \
		"GCC = gcc" \
		"GCCFLAGS = -fno-exceptions " \
		"CPP = cpp" \
		"CPPFLAGS = -cpp -D=LINUX " \
		"F90 = mpifort" \
		"F90FLAGS = -O3" \
		"NVCC = nvcc" \
		"NVCCFLAGS = -O3 -use_fast_math -arch=sm_61 --ptxas-options=-v" \
		"INCD = -I"$(CUDA_PATH)"/include -I/usr/include/mpi " \
		"LIBS = -L"$(CUDA_PATH)"/lib64 -lcudart -lcurand -lcufft -lstdc++ -lfftw3"	
alld:
	make	EXE "PROG = surf36_gpu.exe" \
		"GCC = gcc" \
		"GCCFLAGS = -fno-exceptions " \
		"CPP = cpp" \
		"CPPFLAGS = -cpp -D=LINUX " \
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

#objects rules
#fortran objects rule
interfaccia.o:			interfaccia.f90			init_gpu.o		free_gpu.o		cuda_spec_phys.o	cuda_phys_spec.o
#main.o:					main.f90				interfaccia.o
#cuda objects rule
cuda_spec_phys.o:	cuda_spec_phys.cu	cuda_variables.h	cuda_spec_phys.h	cuda_tran.h
cuda_tran.o:		cuda_tran.cu		cuda_variables.h	cuda_tran.h
cuda_phys_spec.o:	cuda_phys_spec.cu	cuda_variables.h	cuda_phys_spec.h	cuda_tran.h
init_gpu.o:		init_gpu.cu		cuda_variables.h	init_gpu.h
free_gpu.o:		free_gpu.cu		cuda_variables.h	free_gpu.h
