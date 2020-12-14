
#define shell
SHELL = /bin/sh

#set path
MPIFORT_PATH = $(shell which mpifort)
MPI_PATH     = /cineca/prod/opt/compilers/spectrum_mpi/10.3.1/binary/
FFTWLIB      = $(FFTW_LIB)
FFTWINC      = $(FFTW_INC)

FC	= mpifort
RM	= /bin/rm -f
OLEVEL  = -O3
FOPTS   = -cpp -Wall
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    =  -I$(MPI_PATH) -I$(FFTWINC) -L$(FFTWLIB) -lfftw3 -lm -ldl

FLWOBJS =  $(shell cat list_source.list)

MOD     =
OBJS    = $(FLWOBJS)
EXEC    =  ./set_run/sc_compiled/flow36

$(EXEC):        $(OBJS)
				$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
				$(RM) $(OBJS) $(EXEC)

.SUFFIXES: .o

.f90.o:
       	$(FC)  -c $(FFLAGS) $<
