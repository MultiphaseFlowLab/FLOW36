#FC	= mpipgifort #if using mpi spectrum
FC	= mpif90     #if using hpc-sdk
LD	= $(FC)
RM	= /bin/rm -f
FOPTS   = -fast -Mfree -Mpreprocess  -module=./set_run/sc_compiled
FFLAGS	= $(FOPTS)
LIBS    = -L/opt/software/fftw/3/3.3.10-nvidia-openmpi/lib -I/opt/software/fftw/3/3.3.10-nvidia-openmpi/include -lfftw3 -lm -ldl^M


LDFLAGS	= $(FOPTS)

FLWOBJS =  $(shell cat list_source.list)

OBJS	= $(FLWOBJS)
EXEC    =  ./set_run/sc_compiled/flow36

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(OBJS) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
