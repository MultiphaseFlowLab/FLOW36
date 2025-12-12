#F90	= $(FC)

FC	= mpif90
LD	= $(FC)
SWP	=
RM	= /bin/rm -f
MP	=
OLEVEL	= -O3
# linux
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp  #-std=f2008ts -g
# mac
#FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000
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
./read_fields.f90 \
./get_interface.f90\
./flood_fill.f90 \


MOD	=
OBJS	= $(FLWOBJS)
EXEC    =  ./mass_center

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC)

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
