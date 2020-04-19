#F90	= $(FC)

FC	= mpifort
#FC      = mpiifort
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS   = -march=znver2 -ffree-form -mcmodel=large -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp -J ./set_run/sc_compiled -I./set_run/sc_compiled/  #-std=f2008ts -g
FFLAGS	= $(FOPTS) $(OLEVEL)
LIBS    = -I$(FFTW_ROOT)/include -L$(FFTW_ROOT)/lib -lfftw3 -lm -ldl 
# check for AMD-optimized FFTW, see https://kb.hlrs.de/platforms/index.php/Libraries(Hawk)
#LIBS    = -I$(FFTW_ROOT)/include -L$(FFTW_ROOT)/lib -lamd-fftw -lm -ldl 

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS =  $(shell cat list_source.list)

MOD	=
OBJS	= $(FLWOBJS)
EXEC    =  ./set_run/sc_compiled/flow36

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(OBJS) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
