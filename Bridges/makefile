#F90	= $(FC)

FC	= mpifort
#FC      = mpiifort
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -cpp -module ./set_run/sc_compiled/ -Rset_run/sc_compiled/
#FOPTS   = -align none -mcmodel medium -axMIC-AVX512 -cpp -module ./set_run/sc_compiled/ -I./set_run/sc_compiled/
FFLAGS	= $(FOPTS) $(OLEVEL)
LIBS    =  -lfftw3 -lm -ldl

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
