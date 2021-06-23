#F90	= $(FC)

FC	= mpifort
LD	= $(FC)
SWP	=
RM	= /bin/rm -f ./get_tke
MP	=
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -fconvert=big-endian -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp
#FOPTS	= -align none -mcmodel medium -warn all
FFLAGS	= $(FOPTS) $(OLEVEL)
#LIBS	=  -I/usr/include/ -I/usr/local/include -L$FFTW_LIB -L/usr/lib -L/usr/lib64 -lfftw3 -lm -ldl
LIBS    =  -I/usr/include/ -I/usr/local/include -L$FFTW_LIB -lfftw3 -lm -ldl
#LIBS	=  -I$//sw/include/ -L$//sw/lib -lfftw3_mpi -lfftw3 -lm
#LIBS	=  -I$//sw/include -I /usr/include/hdf5  -L/usr/lib -L/usr/lib/hdf5 -L/usr/lib64/ -L/usr/lib/hdf5/bin -lhdf5 -lfftw3 -lm
#LIBS	=  -I/usr/include/ -I/usr/include/hdf5  -L$FFTW_LIB -L/usr/lib -L/usr/lib64 -L/usr/lib/hdf5 -L/usr/lib/hdf5/bin -lhdf5 -lfftw3 -lfftw3_mpi -lm -ldl

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./main.f90 \
./module.f90 \
./read_input.f90 \
./create_plan.f90 \
./destroy_plan.f90 \
./read_fields.f90 \
./phys_to_spectral.f90 \
./spectral_to_phys.f90 \
./fftx_fwd.f90 \
./fftx_bwd.f90 \
./ffty_fwd.f90 \
./ffty_bwd.f90 \
./dctz_fwd.f90 \
./dctz_bwd.f90 \
./xy2xz.f90 \
./xz2xy.f90 \
./xz2yz.f90 \
./yz2xz.f90 \
./dz.f90 \
./load_files.f90 


MOD	=
OBJS	= $(FLWOBJS)
EXEC    =  ./get_tke

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM)  $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
