#F90	= $(FC)

FC	= gfortran
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -J ./modules/ -cpp
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
./wave_numbers.f90 \
./create_plan.f90 \
./destroy_plan.f90 \
./phys_to_spectral.f90 \
./fftx_fwd.f90 \
./ffty_fwd.f90 \
./dctz_fwd.f90 \
./spectral_to_phys.f90 \
./fftx_bwd.f90 \
./ffty_bwd.f90 \
./dctz_bwd.f90 \
./wall_shear.f90 \
./dz.f90


MOD	= ./modules/*.mod
OBJS	= $(FLWOBJS)
EXEC    =  ./wall_shear

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC) $(MOD)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
