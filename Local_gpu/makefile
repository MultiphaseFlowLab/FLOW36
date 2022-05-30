#nvfortran
FC	= mpif90
LD	= $(FC)
RM	= /bin/rm -f
FOPTS	= -fast -acc -gpu=managed -Minfo=accel -Mfree -Mpreprocess -Mconcur -mcmodel=medium -Mlarge_arrays  -module=./set_run/sc_compiled
FFLAGS	= $(FOPTS)
LIBS    = -L/usr/local/cuda/lib64 -lcufft 

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

