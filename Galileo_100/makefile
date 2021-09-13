
#define shell
SHELL = /bin/sh

#set path
MPIFORT_PATH = $(shell which mpifort)
MPI_PATH     = /cineca/prod/opt/modulefiles/base/compilers/intelmpi/
FFTWLIB      = $(FFTW_LIB)
FFTWINC      = $(FFTW_INC)

FC      = mpiifort
LD      = $(FC)
SWP     =
RM      = /bin/rm -f
MP      =
OLEVEL  = -O3
FOPTS   = -align none -mcmodel medium -cpp -module ./set_run/sc_compiled/ -I./set_run/sc_compiled/
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    = -I$(MPI_PATH)  -I$(FFTW_INC) -L$(FFTW_LIB) -lfftw3 -lm -ldl
LDFLAGS = $(FOPTS) $(OLEVEL)
PROF    =

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
