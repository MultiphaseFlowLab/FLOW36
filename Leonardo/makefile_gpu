FC    	= pgf90
LD	= $(FC)
RM	= /bin/rm -f
FOPTS   = -fast -acc -gpu=managed -Minfo=accel -Mfree -Mpreprocess  -module=./set_run/sc_compiled -lnvToolsExt
FFLAGS  = $(FOPTS)
#LIBS   = -I$(MPI_PATH) -cudalib=cufft #MPI path if using spectrum
LIBS    = -cudalib=cufft

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
