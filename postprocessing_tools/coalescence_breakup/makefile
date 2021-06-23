#F90	= $(FC)

FC	= gfortran
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
#FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -J ./modules/ -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -J ./modules/ 
FFLAGS	= $(FOPTS) $(OLEVEL)
LIBS	=  -I/usr/include/ -I/usr/local/include -lfftw3 -lm

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./main.f90 \
./module.f90 \
./read_input.f90 \
./print_start.f90 \
./generate_grid.f90 \
./find_drop.f90 \
./flood_fill.f90 \
./get_center.f90 


MOD	= ./modules/*.mod
OBJS	= $(FLWOBJS)
EXEC    =  ./coalescence_breakup

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC) $(MOD)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
