
#define shell
SHELL = /bin/sh

#set path
MPIFORT_PATH = $(shell which ftn)
MPI_PATH     = $(MPICH_DIR)
FFTWLIB      = /opt/cray/pe/fftw/3.3.8.3/x86_64/lib
FFTWINC      = /opt/cray/pe/fftw/3.3.8.3/x86_64/include


FC      = ftn
RM      = /bin/rm -f
OLEVEL  = -O3
#FOPTS  = -align none -mcmodel medium -cpp -module ./set_run/sc_compiled/ -I./set_run/sc_compiled/
FOPTS   = -cpp
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    =  -I$(MPI_PATH) -I$(FFTWINC) -L$(FFTWLIB) -lfftw3 -lm -ldl


FLWOBJS =  $(shell cat list_source.list)

MOD     =
OBJS    = $(FLWOBJS)
EXEC    =  ./set_run/sc_compiled/surf36_gpu.exe

$(EXEC):        $(OBJS)
        $(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(OBJS) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<

