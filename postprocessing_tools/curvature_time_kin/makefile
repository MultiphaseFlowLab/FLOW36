#F90	= $(FC)

FC	= mpifort
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp  #-std=f2008ts -g
FFLAGS	= $(FOPTS) $(OLEVEL)
LIBS    =  -I/usr/include/ -I/usr/local/include -L$FFTW_LIB -lfftw3 -lm -ldl

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./main.f90 \
./module.f90 \
./read_input.f90 \
./read_fields.f90 \
./calc_curvature.f90 \
./create_plan.f90 \
./destroy_plan.f90 \
./spectral_to_phys.f90 \
./fftx_bwd.f90 \
./ffty_bwd.f90 \
./dctz_bwd.f90 


MOD	=
OBJS	= $(FLWOBJS)
EXEC    =  ./curvature

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
