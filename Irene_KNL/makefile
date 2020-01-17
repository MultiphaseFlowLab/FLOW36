#F90    = $(FC)

FC	= mpifort
LD	= $(FC)
SWP     =
RM	= /bin/rm -f
MP	=
OLEVEL  = -O3
FOPTS   = -align none -mcmodel medium -axMIC-AVX512 -cpp -module ./set_run/sc_compiled/ -I./set_run/sc_compiled/ -I$(FFTW3_INCDIR)
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    = $(FFTW3_FFLAGS) -L$(FFTW3_LIBDIR) $(MKL_FFLAGS) ${MKL_LDFLAGS}
LDFLAGS = $(FOPTS) $(OLEVEL)
PROF    =

FLW_SRC = $(shell cat list_source.list)
FLW_MOD = $(shell cat list_source.list | grep module)

OBJS = $(FLW_SRC:.f90=.o)
MOD = $(FLW_MOD:.f90=.o)

EXEC = ./set_run/sc_compiled/flow36

$(EXEC): $(MOD) $(OBJS)
        $(FC) $(OBJS) $(FFLAGS) $(LIBS) -o $@

%.o: %.f90
        $(FC) $< $(FFLAGS) $(LIBS) -c -o $@

clean:
        $(RM) $(MOD) $(OBJS) $(EXEC)
