# paraview output on Leonardo (GNU compilers + Open MPI + FFTW from the modules), e.g.
#   module load gcc openmpi fftw
# FFTW_LIB / FFTW_INC are set by the fftw module; otherwise pass them: make FFTW_LIB=... FFTW_INC=...

FC	= mpif90
LD	= $(FC)
RM	= /bin/rm -f
OLEVEL	= -O3
FOPTS	= -mcmodel=medium -fconvert=big-endian -ffixed-line-length-140 -fno-align-commons -Wall -cpp
#FOPTS	= -mcmodel=medium -fconvert=big-endian -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp -g
FFLAGS	= $(FOPTS) $(OLEVEL) -I$(FFTW_INC)
LIBS	= -L$(FFTW_LIB) -Wl,-rpath,$(FFTW_LIB) -lfftw3 -lm -ldl

# module.f90 first: the other files use its modules (commondata, fftw3, wavenumber)
FLWOBJS = \
./module.f90 \
./main.f90 \
./read_input.f90 \
./read_fields.f90 \
./generate_output.f90 \
./spectral_to_phys.f90 \
./dctz_bwd.f90 \
./fftx_bwd.f90 \
./ffty_bwd.f90 \
./create_plan.f90 \
./destroy_plan.f90 \
./swap_grid.f90 \
./phys_to_spectral.f90 \
./dctz_fwd.f90 \
./fftx_fwd.f90 \
./ffty_fwd.f90 \
./dz.f90

OBJS	= $(FLWOBJS)
EXEC	= ./read_paraview

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC) *.mod
