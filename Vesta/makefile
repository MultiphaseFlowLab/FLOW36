#F90    = $(FC)

FC	= mpixlf2003_r
LD	= $(FC)
SWP     =
RM	= /bin/rm -f
MP	=
OLEVEL  = -O3
FOPTS   = -qsuffix=cpp=f90  -qcheck=all -w -qarch=qp -qtune=qp -qmoddir=./set_run/sc_compiled/ -I/set_run/sc_compiled/
FFLAGS  = $(FOPTS) $(OLEVEL)
LIBS    =  -I/soft/libraries/alcf/current/xl/FFTW3/include -L/soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3 -lm -ldl

LDFLAGS = $(FOPTS) $(OLEVEL)
PROF	=

FLWOBJS = \
./set_run/sc_compiled/main.f90 \
./set_run/sc_compiled/module.f90 \
./set_run/sc_compiled/read_input.f90 \
./set_run/sc_compiled/print_start.f90 \
./set_run/sc_compiled/define_sizes.f90 \
./set_run/sc_compiled/create_plan.f90 \
./set_run/sc_compiled/destroy_plan.f90 \
./set_run/sc_compiled/initialize.f90 \
./set_run/sc_compiled/dump_grid.f90 \
./set_run/sc_compiled/write_output.f90 \
./set_run/sc_compiled/write_time.f90 \
./set_run/sc_compiled/solver.f90 \
./set_run/sc_compiled/convective_ns.f90 \
./set_run/sc_compiled/wave_numbers.f90 \
./set_run/sc_compiled/phys_to_spectral.f90 \
./set_run/sc_compiled/spectral_to_phys.f90 \
./set_run/sc_compiled/fftx_fwd.f90 \
./set_run/sc_compiled/fftx_bwd.f90 \
./set_run/sc_compiled/ffty_fwd.f90 \
./set_run/sc_compiled/ffty_bwd.f90 \
./set_run/sc_compiled/yz2xz.f90 \
./set_run/sc_compiled/xz2xy.f90 \
./set_run/sc_compiled/xz2yz.f90 \
./set_run/sc_compiled/xy2xz.f90 \
./set_run/sc_compiled/dctz_fwd.f90 \
./set_run/sc_compiled/dctz_bwd.f90 \
./set_run/sc_compiled/dz.f90 \
./set_run/sc_compiled/time_integration.f90 \
./set_run/sc_compiled/hist_term.f90 \
./set_run/sc_compiled/helmholtz.f90 \
./set_run/sc_compiled/helmholtz_red.f90 \
./set_run/sc_compiled/calculate_var.f90 \
./set_run/sc_compiled/courant_check.f90 \
./set_run/sc_compiled/read_fields.f90 \
./set_run/sc_compiled/initialize_phi.f90 \
./set_run/sc_compiled/phi_non_linear.f90 \
./set_run/sc_compiled/sterm_ch.f90 \
./set_run/sc_compiled/statistics.f90 \
./set_run/sc_compiled/sim_check.f90 \
./set_run/sc_compiled/chop_modes.f90 \
./set_run/sc_compiled/initialize_psi.f90 \
./set_run/sc_compiled/sterm_surf.f90 \
./set_run/sc_compiled/initialize_temp.f90 \
./set_run/sc_compiled/sterm_temp.f90 \
./set_run/sc_compiled/hist_term_temp.f90 \
./set_run/sc_compiled/define_address.f90 \
./set_run/sc_compiled/define_dual_grid.f90 \
./set_run/sc_compiled/swap_grid.f90 \
./set_run/sc_compiled/shrink.f90



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
