#F90	= $(FC)

FC	= mpifort
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp -J ./set_run/sc_compiled #-std=f2008ts -g
FFLAGS	= $(FOPTS) $(OLEVEL)
LIBS	=  -I/usr/local/include/ -L/usr/local/lib/ -lfftw3 -lm # -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

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
