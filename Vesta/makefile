#F90    = $(FC)

FC	= mpixlf2003_r
LD	= $(FC)
SWP     =
RM	= /bin/rm -f
MP	=
OLEVEL  = -O3
FOPTS   = -qsuffix=cpp=f90  -qcheck=all -w -qarch=qp -qtune=qp -qmoddir=./set_run/sc_compiled/ -I/set_run/sc_compiled/
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    =  -I/soft/libraries/alcf/current/xl/FFTW3/include -L/soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lm -ldl

LDFLAGS = $(FOPTS) $(OLEVEL)
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
