# define machine
# 0 : local (Mac)
# 1 : local
# 2 : Marconi A1
# 3 : VSC3
# 4 : Vesta (ANL)
# 5 : Marconi A2 KNL
# 6 : Theta (ANL)
# 7 : Bridges (PSC)
machine="7"
echo ""
echo "=============================================================================="
echo "=                                 Running on                                 ="
if [ "$machine" == "0" ]; then
echo "=                             Local machine (Mac)                            ="
cp ./Local_Mac/makefile ./makefile
cp ./Local_Mac/go.sh ./go.sh
# save recovery files in modal space (1) or physical space (0)
savespectral="1"
elif [ "$machine" == "1" ]; then
echo "=                                Local machine                               ="
cp ./Local/makefile ./makefile
cp ./Local/go.sh ./go.sh
savespectral="1"
elif [ "$machine" == "2" ]; then
echo "=                            Marconi A1 Broadwell                            ="
cp ./Marconi/makefile ./makefile
cp ./Marconi/go.sh ./go.sh
module purge
module load profile/global
# load modules
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary
module load fftw/3.3.5--intelmpi--2017--binary
# or
#module load gnu/6.1.0
#module load openmpi/1-10.3--gnu--6.1.0
#module load fftw/3.3.4--openmpi--1-10.3--gnu--6.1.0
savespectral="0"
elif [ "$machine" == "3" ]; then
echo "=                                   VSC-3                                    ="
cp ./VSC-3/makefile ./makefile
cp ./VSC-3/go.sh ./go.sh
module purge
# load modules
module load intel/16 intel-mpi/5.1.3 fftw/3.3.4-DP
# or (but problem with .mod files)
#module load gcc/5.3 intel-mpi/5.1.3 fftw/3.3.4-DP
savespectral="0"
elif [ "$machine" == "4" ]; then
echo "=                                  Vesta                                     ="
cp ./Vesta/makefile ./makefile
cp ./Vesta/go.sh ./go.sh

savespectral="1"
elif [ "$machine" == "5" ]; then
echo "=                              Marconi A2 KNL                                ="
module purge
module load env-knl
module load profile/global
# load modules
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary
module load fftw/3.3.5--intelmpi--2017--binary
cp ./Marconi_KNL/makefile ./makefile
cp ./Marconi_KNL/go.sh ./go.sh

savespectral="0"
elif [ "$machine" == "6" ]; then
echo "=                                  Theta                                     ="
# load modules
#module load gcc
module load fftw
module load craype-hugepages16M

cp ./Theta/makefile ./makefile
cp ./Theta/go.sh ./go.sh

savespectral="0"
elif [ "$machine" == "7" ]; then
echo "=                                  Bridges                                   ="
module purge
# load modules
module load pgi/19.4
module load mpi/pgi_openmpi/19.4
#module load intel
#module load mpi/intel_mpi
module load fftw3/3.3.4
cp ./Bridges/makefile ./makefile
cp ./Bridges/go.sh ./go.sh

savespectral="0"
fi
echo "=============================================================================="
echo ""

################################################################################
# define simulation parameters

# fftw plan craation flag
# 0: FFTW_ESTIMATE, faster plan creation, transforms may be slower
# 1: FFTW_PATIENT, try several algorithms, choose the best one, slower plan creation, transforms may be faster
fftw_flag="0"
# PAY ATTENTION TO VARIABLE TIPE #

# number of grid points (edit only exponent)
ix="1" # integer
iy="8" # integer
iz="8" # integer

# dual grid for surfactant, expansion factors:
exp_x="1" # integer, (2**ix)*exp_x
exp_y="1" # integer, (2**iy)*exp_y
exp_z="1" # integer, (2**iz)*exp_z+1

# parallelization strategy
NYCPU="1" # integer
NZCPU="4" # integer
# running on single shared memory environment (0) or on many (1)
multinode="0" # integer
# number of MPI processes per node
nodesize="68" # integer

################################################################################
# restart flag: 1 restart, 0 new simulation
restart="0" # integer
nt_restart="0" # integer

# initial condition
# 0 : initialize zero velocity
# 1 : laminar Poiseuille flow in x direction
# 2 : laminar Poiseuille flow in y direction
# 3 : read input from file (parallel read)
# 4 : read input from file (serial read)
# 5 : shear flow y direction
# 6 : shear flow x direction
# always keep list of initial conditions updated
incond="0" # integer

# Reynolds number
Re="1.0" # real (double)

# Courant number
Co="0.2" # real (double)

# mean pressure gradient (x and y), defined ad (p_out-p_in)/L
gradpx="0.0" # real (double)
gradpy="0.0" # real (double)

# domain size, divided by pi (z size is always 2, between -1 and 1)
lx="0.001" # real (double)
ly="2.0" # real (double)

# initial time step
nstart="0" # integer

# final time step
nend="5" #integer (up to 8 digits)

# frequency of solution saving in physical space
dump="2000" # integer

# frequency of solution saving in spectral space
sdump="-1" # integer

# frequency of solution saving in spectral space, needed to restart the simulation
#                       from a checkpoint if it crashes (those files are not kept)
failure_dump="500" # integer

# Run time statistics calculation parameters
# frequency of statistics calculation (leave -1 to skip statistics calculation)
st_dump="-1" # integer

# timestep from which begin statistics calculation
start_stats="0" # integer

# flag for mean, rms, skewness, flatness calculation
mean_flag="1" # 0 to skip mean, rms, skewness and flatness calculation, 1 to do it

# flag for budget calculation
budget_flag="1" # 0 to skip budget calculation, 1 to do it

# flag for power spectra calculation
spectra_flag="1" # 0 to skip power spectra calculation, 1 to do it

# dt
dt="1.e-4" # real (exponential)

# 0: no-slip
# 1: free-slip
# 2: y shear flow (+1 at z=1, -1 at z=-1)
# 3: x shear flow (+1 at z=1, -1 at z=-1)
# boundary condition, z=1
bc_upb="0" # integer

# boundary condition, z=-1
bc_lb="0" # integer

################################################################################
# Phase field only
# phase field flag, 0: phase field deactivated, 1: phase field activated
phi_flag="1" # integer

# correction on phi to improve mass conservation
# 0: OFF
# 1: profile-corrected
# 2: flux-corrected
# 3: profile-corrected turned off at the walls
# 4: profile-corrected kill the gradients (filter on gradients lower than threshold 1/(50*Ch)
# 5: flux-corrected kill the gradients (filter on gradients lower than threshold 1/(50*Ch)
phicor_flag="4" # integer

# Value of the parameter lambda used to correct the phi profile (only for phicor_flag 1 or 2)
# Lam=0.3/Ch
lamcorphi="2.5" # real (double)

# matched densities: 1 for matched densities, 0 for rhor < 1, 2 for rhor > 1
matchedrho="1" # integer

# density ratio, phase +1 over phase -1
rhor="1.0" # real (double)

#matched dynamic viscosities: 1 for matched viscosities, 0 for visr < 1 (or non-newtonian), 2 for visr > 1
matchedvis="1" # integer

# dynamic viscosity ratio, phase +1 over phase -1 (not considered when non-newtonian is enabled)
visr="1.0" # real (double)

#non-newtonian phase, 0 deactivaed, 1 phase=+1 is non-newtonian (Carreau model)
non_newtonian="0" # integer

#Ratio between the viscosity at zero and infinity shear rate (Non-newtonian-Carreau model)
muinfmuzero="0.1"

#Exponent for the non-newtonian fluid (phi=+1)-n<1 pseudoplastic
exp_non_new="0.9"

# Weber number
We="1.0" # real (double)

# Cahn number
Ch="0.02" # real (double)

# Peclet number
Pe="100.0" # real (double)

# Froud number
Fr="0.1" # real (double)

# Body force flag, 0: deactivated, 1: activated
body_flag="0" # integer

# Body force coefficient
Bd="4.0" # real (double)

# Body force direction
#  1: positive x direction
# -1: negative x direction
#  2: positive z direction
# -2: negative z direction
#  3: positive y direction
# -3: negative y direction
bodydir="2" # integer

# electric force flag, 0: deactivated, 1: activated
ele_flag="0" # integer

# Stuart number
stuart="1.0" # real (double)

# initial conditions on phi
# 0: only phase -1
# 1: read input from file (parallel read)
# 2: read input from file (serial read)
# 3: 2D drop (radius,height)
# 4: 3D drop (radius,height)
# 5: stratified flow (mean height, sine wave amplitude, sine wave frequency, perturbation amplitude)
# 6: 3D drop array (radius, height, number of drops x direction, number of drops y direction, number of drops z direction),
#                 x, y and z drop centers distance must be at least 2*(radius+5*sqrt(2)*Ch)
#                 otherwise number of drops will be reduced
# 7: Drop attached to the bottom wall z_c=-1 (radius)
# 8: 2x 2D Droplets in kissing mode. (radius, ygap , zgap)
# 9: Layer of phi=+1 (mean height, thickness)
in_condphi="3" # integer
radius="0.5" # real (double)
height="0.0" # real (double)
wave_amp_x="0.0" # real (double)
wave_freq_x="0.0" # real (double)
wave_amp_y="0.0" # real (double)
wave_freq_y="0.0" # real (double)
pert_amp="0.0" # real (double)
num_x="5" # integer
num_y="2" # integer
num_z="3" # integer
ygap="1.0" # real(double)
zgap="0.25" # real (double)
thickness="0.1086" # real (double)

# gravity direction
#  1: positive x direction
# -1: negative x direction
#  2: positive z direction
# -2: negative z direction
#  3: positive y direction
# -3: negative y direction
gravdir="-2" # integer

# buoyancy type
# 0: no buoyancy and gravity
# 1: buoyancy and gravity (rho*g)
# 2: only buoyancy (Delta rho*g)
buoyancy="0" # integer

################################################################################
# Surfactant only
# surfactant flag, 0 : surfactant deactivated, 1 : surfactant activated
psi_flag="1" # integer

# surfactant Peclet number
Pe_psi="100.0" # real (double)

# Ex number
Ex="0.117" # real (double)

# Pi number
PI="1.35" # real (double)

# Elasticity number
El="0.5" # real (double)

# Initial conditions on the surfactant
# 0: initialize constant value (psi_mean)
# 1: read input from file (parallel read)
# 2: initialize equilibrium profile (psi_bulk)
# 3: equilibrium profile multiplied with Y gradient
# 4: equilibrium profile multiplied with Z gradient
# 5: Diffusion Test, angular distribuction
# 6: read input from file (parallel read, fine grid)
in_condpsi="2"
psi_mean="0.01" # real (double)
psi_bulk="0.01" # real (double)
################################################################################
# Temperature only
# temperature flag, 0 : temperature deactivated, 1 : temperature activated
temp_flag="0" # integer

# Rayleigh number
# for Rayleigh-Benard choose Re=sqrt(Ra/Pr)/4
Ra="10000.0" # real (double)

# Prandtl number
Pr="1.0" # real (double)

# boundary conditions
# A*T+B*dT/dz=C  at z=-1
# D*T+E*dT/dz=F  at z=+1
A="1.0" # real (double)
B="0.0" # real (double)
C="1.0" # real (double)
D="1.0" # real (double)
E="0.0" # real (double)
F="-1.0" # real (double)

# initial conditons for the temperature
# 0 : initialize constant temperature (mean_t)
# 1 : read from data file (parallel read)
in_cond_temp="1" # integer
temp_mean="0.0" # real (double)

# 1 activate buoyancy term in N-S, 0 deactivate it (Boussinnesq approximation)
# uses same gravity array as defined in the phase field part
boussinnesq="0" # integer

################################################################################
# Lagrangian Particle Tracking only
part_flag="0" # integer
part_number="1000" # integer
# 1 use tracer particles (implies 1-way coupling), 0 use inertial particles
tracer="0" # integer
# stokes number (in wall units)
stokes="1.0" # real (double)
# drag type, 1 Stokes drag, 0 Schiller-Naumann drag
stokes_drag="1" # integer
# density ratio particle/fluid
dens_part="1.0" # real (double)
# 1 to activate gravity and buoyancy force on particle tracking
part_gravity="1" # integer

# 1 activate two-way coupling, 0 deactivate it
twoway="1" # integer

# initial conditions for the particle position
# 0 : initialize random position
# 1 : read from input file (parallel read, binary file)
# 2 : initialize random position on a x-y plane at height par_plane (part_height)
in_cond_part_pos="0" # integer
part_height="0.0" # real (double) between -1 and +1

# initial conditions for the particle velocity
# 0 : zero velocity
# 1 : fluid velocity at particle position
# 2 : read from input file (parallel read, binary file)
in_cond_part_vel="0" # integer

# end of parameters declaration
################################################################################

echo ""
echo "       FFFFFFF  L        OOO   W           W           333      666"
echo "       F        L       O   O  W     W     W          3   3    6   6"
echo "       F        L       O   O   W   W W   W               3    6"
echo "       FFFFF    L       O   O   W   W W   W              3     6666"
echo "       F        L       O   O    W W   W W                3    6   6"
echo "       F        L       O   O    W W   W W            3   3    6   6"
echo "       F        LLLLLL   OOO      W     W              333      666"
echo ""
echo "                                +hhy/  "
echo "                     ::       -ddddd   "
echo "                     ydy-      :syo-       "
echo "                     +hdho//:/+osyhhdhhyo+:."
echo "                         :oyhddddddddd+//+oshddo"
echo "                               sdddddd-      ."
echo "                               /ddddddo "
echo "                              +hddddhssooshddy      -/osh/"
echo "                             ydds:.        hddo+oyddddh+"
echo "                             hdd: .-:/+osyhddddddddds:"
echo "                            -sddddddddddddddddddho:"
echo "                          oddddddddddddddddys+:."
echo "                          ./+osdddo++//:-."
echo "                             .yho."
echo "                             :."


echo ""
echo "=============================================================================="
echo "=                            START OF NEW RUN                                ="
echo "=                      END OF PARAMETER DECLARATION                          ="
echo "=============================================================================="
echo ""

if [ "$machine" == "0" ]; then
echo ""
echo "==============================OS X Version===================================="
echo ""
echo "                                   ###                                        "
echo "                                 ####                                         "
echo "                                  ###                                         "
echo "                           #######   #######                                  "
echo "                         #####################                                "
echo "                        #####################                                 "
echo "                        ####################                                  "
echo "                        ####################                                  "
echo "                        #####################                                 "
echo "                         ######################                               "
echo "                          ###################                                 "
echo "                            ###############                                   "
echo "                             ####   #####                                     "
echo ""
echo "=============================================================================="
fi

mkdir -p set_run

rm -r ./set_run/go.sh
rm -r ./set_run/sc_compiled
rm -r ./set_run/paraview_output_fg
rm -r ./set_run/stats_calc
rm -r ./set_run/nohup.out
if [ "$restart" == "0" ]; then
 rm -r ./set_run/initial_fields/*
 rm -r ./set_run/results

 mkdir ./set_run/results
 mkdir ./set_run/results/backup
 cp -r ./initial_fields ./set_run/
fi

cp ./restart_copy.sh ./set_run/results/backup/
mkdir ./set_run/sc_compiled

# grid size
NX="$((2**$ix))"
NY="$((2**$iy))"
NZ="$(((2**$iz)+1))"

# set total number of MPI processes requested
if [[ "$multinode" == "1" && "$part_flag" == 1 ]]; then
 NNT="$(($NYCPU*$NZCPU+$nodesize))"
elif [[ "$multinode" == "0"  && "$part_flag" == 1 ]]; then
 NNT="$(($NYCPU*$NZCPU))"
elif [ "$part_flag" == 0 ]; then
 NNT="$(($NYCPU*$NZCPU))"
fi

# copy executable and edit it
cp ./go.sh ./set_run
if [ "$machine" == "0" ]; then
sed -i "" "s/NUMTASKS/$NNT/g" ./set_run/go.sh
else
sed -i "s/NUMTASKS/$NNT/g" ./set_run/go.sh
fi

# copy input file and edit it
cp ./input.f90 ./set_run/sc_compiled
if [ "$machine" == "0" ]; then
sed -i "" "s/restartflag/$restart/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/restart_iteration/$nt_restart/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/initialcondition/$incond/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nxxxxxx/$NX/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nyyyyyy/$NY/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nzzzzzz/$NZ/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Renum/$Re/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/courantnum/$Co/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/gradpx/$gradpx/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/gradpy/$gradpy/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/len_x/$lx/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/len_y/$ly/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nstart/$nstart/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nend/$nend/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nfdump/$dump/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nsdump/$sdump/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/faildump/$failure_dump/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/stats_dump/$st_dump/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/stats_start/$start_stats/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/delta_t/$dt/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/bc_upbound/$bc_upb/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/bc_lowbound/$bc_lb/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/phasephiflag/$phi_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/phaseprofflag/$phicor_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/lamcorphi/$lamcorphi/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/matcheddens/$matchedrho/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/densratio/$rhor/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/matchedvisc/$matchedvis/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/viscratio/$visr/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/nonnewtonian/$non_newtonian/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/muinfmuzero/$muinfmuzero/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/expnonnew/$exp_non_new/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/webernumber/$We/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/cahnnumber/$Ch/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/pecletnumber/$Pe/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/froudnumber/$Fr/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/phinitial_condition/$in_condphi/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/gravitydir/$gravdir/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/gravitytype/$buoyancy/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/bodyforce/$body_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/bodyfcoeff/$Bd/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/bodydirection/$bodydir/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/eleforce/$ele_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/elefcoeff/$stuart/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/surfactantflag/$psi_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/surfpeclet/$Pe_psi/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/exnumber/$Ex/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/pinumber/$PI/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/surfelasticity/$El/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/psinitial_condition/$in_condpsi/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/temperatureflag/$temp_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/rayleighnumb/$Ra/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/prandtlnumb/$Pr/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Aboundary/$A/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Bboundary/$B/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Cboundary/$C/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Dboundary/$D/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Eboundary/$E/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/Fboundary/$F/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/tempinitial_condition/$in_cond_temp/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/particleflag/$part_flag/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/particlenumber/$part_number/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/partstokes/$stokes/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/densityparticle/$dens_part/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/incondpartpos/$in_cond_part_pos/g" ./set_run/sc_compiled/input.f90
sed -i "" "s/incondpartvel/$in_cond_part_vel/g" ./set_run/sc_compiled/input.f90
else
sed -i "s/restartflag/$restart/g" ./set_run/sc_compiled/input.f90
sed -i "s/restart_iteration/$nt_restart/g" ./set_run/sc_compiled/input.f90
sed -i "s/initialcondition/$incond/g" ./set_run/sc_compiled/input.f90
sed -i "s/nxxxxxx/$NX/g" ./set_run/sc_compiled/input.f90
sed -i "s/nyyyyyy/$NY/g" ./set_run/sc_compiled/input.f90
sed -i "s/nzzzzzz/$NZ/g" ./set_run/sc_compiled/input.f90
sed -i "s/Renum/$Re/g" ./set_run/sc_compiled/input.f90
sed -i "s/courantnum/$Co/g" ./set_run/sc_compiled/input.f90
sed -i "s/gradpx/$gradpx/g" ./set_run/sc_compiled/input.f90
sed -i "s/gradpy/$gradpy/g" ./set_run/sc_compiled/input.f90
sed -i "s/len_x/$lx/g" ./set_run/sc_compiled/input.f90
sed -i "s/len_y/$ly/g" ./set_run/sc_compiled/input.f90
sed -i "s/nstart/$nstart/g" ./set_run/sc_compiled/input.f90
sed -i "s/nend/$nend/g" ./set_run/sc_compiled/input.f90
sed -i "s/nfdump/$dump/g" ./set_run/sc_compiled/input.f90
sed -i "s/nsdump/$sdump/g" ./set_run/sc_compiled/input.f90
sed -i "s/faildump/$failure_dump/g" ./set_run/sc_compiled/input.f90
sed -i "s/stats_dump/$st_dump/g" ./set_run/sc_compiled/input.f90
sed -i "s/stats_start/$start_stats/g" ./set_run/sc_compiled/input.f90
sed -i "s/delta_t/$dt/g" ./set_run/sc_compiled/input.f90
sed -i "s/bc_upbound/$bc_upb/g" ./set_run/sc_compiled/input.f90
sed -i "s/bc_lowbound/$bc_lb/g" ./set_run/sc_compiled/input.f90
sed -i "s/phasephiflag/$phi_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/phaseprofflag/$phicor_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/lamcorphi/$lamcorphi/g" ./set_run/sc_compiled/input.f90
sed -i "s/matcheddens/$matchedrho/g" ./set_run/sc_compiled/input.f90
sed -i "s/densratio/$rhor/g" ./set_run/sc_compiled/input.f90
sed -i "s/matchedvisc/$matchedvis/g" ./set_run/sc_compiled/input.f90
sed -i "s/viscratio/$visr/g" ./set_run/sc_compiled/input.f90
sed -i "s/nonnewtonian/$non_newtonian/g" ./set_run/sc_compiled/input.f90
sed -i "s/muinfmuzero/$muinfmuzero/g" ./set_run/sc_compiled/input.f90
sed -i "s/expnonnew/$exp_non_new/g" ./set_run/sc_compiled/input.f90
sed -i "s/webernumber/$We/g" ./set_run/sc_compiled/input.f90
sed -i "s/cahnnumber/$Ch/g" ./set_run/sc_compiled/input.f90
sed -i "s/pecletnumber/$Pe/g" ./set_run/sc_compiled/input.f90
sed -i "s/froudnumber/$Fr/g" ./set_run/sc_compiled/input.f90
sed -i "s/phinitial_condition/$in_condphi/g" ./set_run/sc_compiled/input.f90
sed -i "s/gravitydir/$gravdir/g" ./set_run/sc_compiled/input.f90
sed -i "s/gravitytype/$buoyancy/g" ./set_run/sc_compiled/input.f90
sed -i "s/bodyforce/$body_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/bodyfcoeff/$Bd/g" ./set_run/sc_compiled/input.f90
sed -i "s/bodydirection/$bodydir/g" ./set_run/sc_compiled/input.f90
sed -i "s/eleforce/$ele_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/elefcoeff/$stuart/g" ./set_run/sc_compiled/input.f90
sed -i "s/surfactantflag/$psi_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/surfpeclet/$Pe_psi/g" ./set_run/sc_compiled/input.f90
sed -i "s/exnumber/$Ex/g" ./set_run/sc_compiled/input.f90
sed -i "s/pinumber/$PI/g" ./set_run/sc_compiled/input.f90
sed -i "s/surfelasticity/$El/g" ./set_run/sc_compiled/input.f90
sed -i "s/psinitial_condition/$in_condpsi/g" ./set_run/sc_compiled/input.f90
sed -i "s/temperatureflag/$temp_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/rayleighnumb/$Ra/g" ./set_run/sc_compiled/input.f90
sed -i "s/prandtlnumb/$Pr/g" ./set_run/sc_compiled/input.f90
sed -i "s/Aboundary/$A/g" ./set_run/sc_compiled/input.f90
sed -i "s/Bboundary/$B/g" ./set_run/sc_compiled/input.f90
sed -i "s/Cboundary/$C/g" ./set_run/sc_compiled/input.f90
sed -i "s/Dboundary/$D/g" ./set_run/sc_compiled/input.f90
sed -i "s/Eboundary/$E/g" ./set_run/sc_compiled/input.f90
sed -i "s/Fboundary/$F/g" ./set_run/sc_compiled/input.f90
sed -i "s/tempinitial_condition/$in_cond_temp/g" ./set_run/sc_compiled/input.f90
sed -i "s/particleflag/$part_flag/g" ./set_run/sc_compiled/input.f90
sed -i "s/particlenumber/$part_number/g" ./set_run/sc_compiled/input.f90
sed -i "s/partstokes/$stokes/g" ./set_run/sc_compiled/input.f90
sed -i "s/densityparticle/$dens_part/g" ./set_run/sc_compiled/input.f90
sed -i "s/incondpartpos/$in_cond_part_pos/g" ./set_run/sc_compiled/input.f90
sed -i "s/incondpartvel/$in_cond_part_vel/g" ./set_run/sc_compiled/input.f90
fi
# end of input file editing


# copy source files
cp ./source_code/module.f90 ./set_run/sc_compiled/
cp ./source_code/main.f90 ./set_run/sc_compiled/
cp ./source_code/read_input.f90 ./set_run/sc_compiled/
cp ./source_code/print_start.f90 ./set_run/sc_compiled/
cp ./source_code/define_sizes.f90 ./set_run/sc_compiled/
cp ./source_code/initialize.f90 ./set_run/sc_compiled/
cp ./source_code/dump_grid.f90 ./set_run/sc_compiled/
cp ./source_code/write_time.f90 ./set_run/sc_compiled/
cp ./source_code/solver.f90 ./set_run/sc_compiled/
cp ./source_code/convective_ns.f90 ./set_run/sc_compiled/
cp ./source_code/wave_numbers.f90 ./set_run/sc_compiled/
cp ./source_code/phys_to_spectral.f90 ./set_run/sc_compiled/
cp ./source_code/spectral_to_phys.f90 ./set_run/sc_compiled/
cp ./source_code/write_output.f90 ./set_run/sc_compiled/
cp ./source_code/fftx_fwd.f90 ./set_run/sc_compiled/
cp ./source_code/fftx_bwd.f90 ./set_run/sc_compiled/
cp ./source_code/ffty_fwd.f90 ./set_run/sc_compiled/
cp ./source_code/ffty_bwd.f90 ./set_run/sc_compiled/
cp ./source_code/yz2xz.f90 ./set_run/sc_compiled/
cp ./source_code/xz2xy.f90 ./set_run/sc_compiled/
cp ./source_code/xz2yz.f90 ./set_run/sc_compiled/
cp ./source_code/xy2xz.f90 ./set_run/sc_compiled/
cp ./source_code/dctz_fwd.f90 ./set_run/sc_compiled/
cp ./source_code/dctz_bwd.f90 ./set_run/sc_compiled/
cp ./source_code/dz.f90 ./set_run/sc_compiled/
cp ./source_code/create_plan.f90 ./set_run/sc_compiled/
cp ./source_code/destroy_plan.f90 ./set_run/sc_compiled/
cp ./source_code/time_integration.f90 ./set_run/sc_compiled/
cp ./source_code/hist_term.f90 ./set_run/sc_compiled/
cp ./source_code/helmholtz.f90 ./set_run/sc_compiled/
cp ./source_code/helmholtz_red.f90 ./set_run/sc_compiled/
cp ./source_code/calculate_var.f90 ./set_run/sc_compiled/
cp ./source_code/courant_check.f90 ./set_run/sc_compiled/
cp ./source_code/read_fields.f90 ./set_run/sc_compiled/
cp ./source_code/initialize_phi.f90 ./set_run/sc_compiled/
cp ./source_code/phi_non_linear.f90 ./set_run/sc_compiled/
cp ./source_code/sterm_ch.f90 ./set_run/sc_compiled/
cp ./source_code/statistics.f90 ./set_run/sc_compiled/
cp ./source_code/sim_check.f90 ./set_run/sc_compiled/
cp ./source_code/chop_modes.f90 ./set_run/sc_compiled/
cp ./source_code/initialize_psi.f90 ./set_run/sc_compiled/
cp ./source_code/sterm_surf.f90 ./set_run/sc_compiled/
cp ./source_code/initialize_temp.f90 ./set_run/sc_compiled/
cp ./source_code/sterm_temp.f90 ./set_run/sc_compiled/
cp ./source_code/hist_term_temp.f90 ./set_run/sc_compiled/
cp ./source_code/define_address.f90 ./set_run/sc_compiled/
cp ./source_code/define_dual_grid.f90 ./set_run/sc_compiled/
cp ./source_code/swap_grid.f90 ./set_run/sc_compiled/
cp ./source_code/shrink.f90 ./set_run/sc_compiled/
cp ./source_code/split_comm.f90 ./set_run/sc_compiled/
cp ./source_code/initialize_particle.f90 ./set_run/sc_compiled/
cp ./source_code/part_fluid_comm.f90 ./set_run/sc_compiled/
cp ./source_code/velocity_interpolator.f90 ./set_run/sc_compiled/
cp ./source_code/lagrangian_tracker.f90 ./set_run/sc_compiled/
cp ./source_code/save_flow_comm.f90 ./set_run/sc_compiled/

cp -r ./paraview_output_fg ./set_run
cp -r ./stats_calc ./set_run



if [ "$phi_flag" == "1" ]; then
  if [ "$in_condphi" == "3" ]; then
    echo "$radius                         ! drop radius" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$height                         ! drop height" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "4" ]; then
    echo "$radius                         ! drop radius" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$height                         ! drop height" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "5" ]; then
    echo "$height                         ! mean height" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$wave_amp_x                     ! wave amplitude x" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$wave_freq_x                    ! wave frequency x" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$wave_amp_y                     ! wave amplitude y" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$wave_freq_y                    ! wave frequency y" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$pert_amp                       ! perturbation amplitude" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "6" ]; then
    echo "$radius                         ! drop radius" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$height                         ! drop height" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$num_x                          ! number of drops, x direction" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$num_y                          ! number of drops, y direction" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$num_z                          ! number of drops, z direction" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "7" ]; then
    echo "$radius                         ! drop radius" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$height                         ! drop height" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "8" ]; then
    echo "$radius                         ! drop radius" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$ygap                           ! ygap" >> ./set_run/sc_compiled/input_phase_field.f90
    echo "$zgap                           ! zgap" >> ./set_run/sc_compiled/input_phase_field.f90
  elif [ "$in_condphi" == "9" ]; then
    echo "$thickness                      ! layer thickness" > ./set_run/sc_compiled/input_phase_field.f90
    echo "$height                         ! height" >> ./set_run/sc_compiled/input_phase_field.f90
  fi
fi

if [ "$psi_flag" == "1" ]; then
  if [ "$in_condpsi" == "0" ]; then
    echo "$psi_mean                         ! mean surfactant" > ./set_run/sc_compiled/input_surfactant.f90
  elif [[ "$in_condpsi" -ge "2"  &&  "$in_condpsi" -le "5" ]]; then
    echo "$psi_bulk                         ! bulk surfactant" > ./set_run/sc_compiled/input_surfactant.f90
  fi
fi

if [ "$temp_flag" == "1" ]; then
  if [ "$in_cond_temp" == "0" ]; then
    echo "$temp_mean                         ! mean temperature" > ./set_run/sc_compiled/input_temperature.f90
  fi
fi

if [ "$part_flag" == "1" ]; then
  if [ "$in_cond_part_pos" == "2" ]; then
    echo "$part_height                         ! z height of x-y particle layer" > ./set_run/sc_compiled/input_particle.f90
  fi
fi

endianness=$(echo -n I | od -to2 | head -n1 | cut -f2 -d" " | cut -c6)
echo "nycpu=$NYCPU;" > ./set_run/results/input_param.m
echo "nzcpu=$NZCPU;" >> ./set_run/results/input_param.m
echo "ntask=$NNT;" >> ./set_run/results/input_param.m
echo "nx=$NX;" >> ./set_run/results/input_param.m
echo "ny=$NY;" >> ./set_run/results/input_param.m
echo "nz=$NZ;" >> ./set_run/results/input_param.m
echo "nend=$nend;" >> ./set_run/results/input_param.m
echo "little_endian=$endianness;" >> ./set_run/results/input_param.m

if [ "$machine" == "0" ]; then
sed -i "" "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/nnx/$NX/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/nny/$NY/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/nnz/$NZ/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/phys_to_spectral.f90
sed -i "" "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/phys_to_spectral.f90
sed -i "" "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/spectral_to_phys.f90
sed -i "" "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/spectral_to_phys.f90
sed -i "" "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/convective_ns.f90
sed -i "" "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/wave_numbers.f90
sed -i "" "s/boundary_conditions_up/$bc_upb/g" ./set_run/sc_compiled/calculate_var.f90
sed -i "" "s/boundary_conditions_down/$bc_lb/g" ./set_run/sc_compiled/calculate_var.f90
sed -i "" "s/precisionflag/$fftw_flag/g" ./set_run/sc_compiled/create_plan.f90
sed -i "" "s/physical_dump_frequency/$dump/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/spectral_dump_frequency/$sdump/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/stats_dump_frequency/$st_dump/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/stats_dump_frequency/$st_dump/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/convective_ns.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/courant_check.f90
sed -i "" "s/phicorcompflag/$phicor_flag/g" ./set_run/sc_compiled/sterm_ch.f90
sed -i "" "s/bodycompflag/$body_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/non_newtonian/$non_newtonian/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/courant_check.f90
sed -i "" "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/matched_viscosity/$matchedvis/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/matched_viscosity/$matchedvis/g" ./set_run/sc_compiled/wave_numbers.f90
sed -i "" "s/buoyancytype/$buoyancy/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/initialize_phi.f90
sed -i "" "s/meanflag/$mean_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "" "s/budgetflag/$budget_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "" "s/spectraflag/$spectra_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "" "s/savespectralflag/$savespectral/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "" "s/boussinnesqcompflag/$boussinnesq/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/boussinnesqcompflag/$boussinnesq/g" ./set_run/sc_compiled/print_start.f90
sed -i "" "s/machineflag/$machine/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/machineflag/$machine/g" ./set_run/sc_compiled/courant_check.f90
sed -i "" "s/machineflag/$machine/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/expansionx/$exp_x/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/expansiony/$exp_y/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/expansionz/$exp_z/g" ./set_run/sc_compiled/module.f90
sed -i "" "s/expansionx/$exp_x/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "" "s/expansiony/$exp_y/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "" "s/expansionz/$exp_z/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "" "s/expansionx/$exp_x/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/expansiony/$exp_y/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/expansionz/$exp_z/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/elecompflag/$ele_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "" "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/main.f90
sed -i "" "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/split_comm.f90
sed -i "" "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "" "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/solver.f90
sed -i "" "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/initialize_particle.f90
sed -i "" "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "" "s/machineflag/$machine/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/physical_dump_frequency/$dump/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/spectral_dump_frequency/$sdump/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "" "s/tracerflag/$tracer/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "" "s/tracerflag/$tracer/g" ./set_run/sc_compiled/print_start.f90
sed -i "" "s/stokesflag/$stokes_drag/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "" "s/stokesflag/$stokes_drag/g" ./set_run/sc_compiled/print_start.f90
sed -i "" "s/nnx/$NX/g" ./set_run/sc_compiled/velocity_interpolator.f90
sed -i "" "s/activategravity/$part_gravity/g" ./set_run/sc_compiled/lagrangian_tracker.f90
else
sed -i "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/module.f90
sed -i "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/module.f90
sed -i "s/nnx/$NX/g" ./set_run/sc_compiled/module.f90
sed -i "s/nny/$NY/g" ./set_run/sc_compiled/module.f90
sed -i "s/nnz/$NZ/g" ./set_run/sc_compiled/module.f90
sed -i "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/phys_to_spectral.f90
sed -i "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/phys_to_spectral.f90
sed -i "s/nnycpu/$NYCPU/g" ./set_run/sc_compiled/spectral_to_phys.f90
sed -i "s/nnzcpu/$NZCPU/g" ./set_run/sc_compiled/spectral_to_phys.f90
sed -i "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/convective_ns.f90
sed -i "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/solver.f90
sed -i "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/wave_numbers.f90
sed -i "s/boundary_conditions_up/$bc_upb/g" ./set_run/sc_compiled/calculate_var.f90
sed -i "s/boundary_conditions_down/$bc_lb/g" ./set_run/sc_compiled/calculate_var.f90
sed -i "s/precisionflag/$fftw_flag/g" ./set_run/sc_compiled/create_plan.f90
sed -i "s/physical_dump_frequency/$dump/g" ./set_run/sc_compiled/main.f90
sed -i "s/spectral_dump_frequency/$sdump/g" ./set_run/sc_compiled/main.f90
sed -i "s/stats_dump_frequency/$st_dump/g" ./set_run/sc_compiled/main.f90
sed -i "s/stats_dump_frequency/$st_dump/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/main.f90
sed -i "s/marconi_flag/$marconi_flag/g" ./set_run/sc_compiled/main.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/convective_ns.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/courant_check.f90
sed -i "s/phicorcompflag/$phicor_flag/g" ./set_run/sc_compiled/sterm_ch.f90
sed -i "s/bodycompflag/$body_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/non_newtonian/$non_newtonian/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/main.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/courant_check.f90
sed -i "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/matched_viscosity/$matchedvis/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/matched_viscosity/$matchedvis/g" ./set_run/sc_compiled/wave_numbers.f90
sed -i "s/buoyancytype/$buoyancy/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/matched_density/$matchedrho/g" ./set_run/sc_compiled/initialize_phi.f90
sed -i "s/meanflag/$mean_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "s/budgetflag/$budget_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "s/spectraflag/$spectra_flag/g" ./set_run/sc_compiled/statistics.f90
sed -i "s/savespectralflag/$savespectral/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/main.f90
sed -i "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/sim_check.f90
sed -i "s/boussinnesqcompflag/$boussinnesq/g" ./set_run/sc_compiled/solver.f90
sed -i "s/boussinnesqcompflag/$boussinnesq/g" ./set_run/sc_compiled/print_start.f90
sed -i "s/machineflag/$machine/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/machineflag/$machine/g" ./set_run/sc_compiled/courant_check.f90
sed -i "s/machineflag/$machine/g" ./set_run/sc_compiled/main.f90
sed -i "s/expansionx/$exp_x/g" ./set_run/sc_compiled/module.f90
sed -i "s/expansiony/$exp_y/g" ./set_run/sc_compiled/module.f90
sed -i "s/expansionz/$exp_z/g" ./set_run/sc_compiled/module.f90
sed -i "s/expansionx/$exp_x/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "s/expansiony/$exp_y/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "s/expansionz/$exp_z/g" ./set_run/sc_compiled/swap_grid.f90
sed -i "s/expansionx/$exp_x/g" ./set_run/sc_compiled/solver.f90
sed -i "s/expansiony/$exp_y/g" ./set_run/sc_compiled/solver.f90
sed -i "s/expansionz/$exp_z/g" ./set_run/sc_compiled/solver.f90
sed -i "s/elecompflag/$ele_flag/g" ./set_run/sc_compiled/phi_non_linear.f90
sed -i "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/main.f90
sed -i "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/split_comm.f90
sed -i "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/solver.f90
sed -i "s/particlecompflag/$part_flag/g" ./set_run/sc_compiled/write_output.f90
sed -i "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/solver.f90
sed -i "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/initialize_particle.f90
sed -i "s/twowaycflag/$twoway/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "s/machineflag/$machine/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/phicompflag/$phi_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/psicompflag/$psi_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/tempcompflag/$temp_flag/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/physical_dump_frequency/$dump/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/spectral_dump_frequency/$sdump/g" ./set_run/sc_compiled/save_flow_comm.f90
sed -i "s/tracerflag/$tracer/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "s/tracerflag/$tracer/g" ./set_run/sc_compiled/print_start.f90
sed -i "s/stokesflag/$stokes_drag/g" ./set_run/sc_compiled/lagrangian_tracker.f90
sed -i "s/stokesflag/$stokes_drag/g" ./set_run/sc_compiled/print_start.f90
sed -i "s/nnx/$NX/g" ./set_run/sc_compiled/velocity_interpolator.f90
sed -i "s/activategravity/$part_gravity/g" ./set_run/sc_compiled/lagrangian_tracker.f90
fi

if [ "$machine" == "4" ]; then
sed -i "s/!onlyforvesta/logical	:: mpi_async_protects_nonblocking/g" ./set_run/sc_compiled/xy2xz.f90
sed -i "s/!onlyforvesta/logical	:: mpi_async_protects_nonblocking/g" ./set_run/sc_compiled/xz2xy.f90
sed -i "s/!onlyforvesta/logical	:: mpi_async_protects_nonblocking/g" ./set_run/sc_compiled/xz2yz.f90
sed -i "s/!onlyforvesta/logical	:: mpi_async_protects_nonblocking/g" ./set_run/sc_compiled/yz2xz.f90
fi

if [ "$machine" == "7" ]; then
# OpenMPI requires iadd and number to be integer(KIND=MPI_ADDRESS_KIND)
sed -i "s/integer :: iadd/integer(KIND=MPI_ADDRESS_KIND) :: iadd/g" ./set_run/sc_compiled/xy2xz.f90
sed -i "s/integer :: iadd/integer(KIND=MPI_ADDRESS_KIND) :: iadd/g" ./set_run/sc_compiled/xz2xy.f90
sed -i "s/integer :: iadd/integer(KIND=MPI_ADDRESS_KIND) :: iadd/g" ./set_run/sc_compiled/xz2yz.f90
sed -i "s/integer :: iadd/integer(KIND=MPI_ADDRESS_KIND) :: iadd/g" ./set_run/sc_compiled/yz2xz.f90
sed -i "s/integer :: number/integer(KIND=MPI_ADDRESS_KIND) :: number/g" ./set_run/sc_compiled/initialize_particle.f90
# PGI compiler does not have isnan
sed -i "s/!only for PGI compiler/use, intrinsic :: ieee_arithmetic/g" ./set_run/sc_compiled/courant_check.f90
fi


# only for intel compiler (needed for USE MPI_F08
#source /opt/intel/compilers_and_libraries_2017/linux/bin/compilervars.sh intel64


echo ""
echo "=============================================================================="
echo "=                      BEGINNING OF COMPILATION                              ="
echo "=============================================================================="
echo ""

# double make needed because first one return errors for missing modules, but then creates them,
# second make makes the code executable with the proper module
# modules must be removed to update data inside them when changing simulation parameters like
# nx, ny, nz, nycpu, nzcpu
make

make

echo ""
echo "=============================================================================="
echo "=                            END OF COMPILATION                              ="
echo "=============================================================================="
echo ""

echo ""
echo "=============================================================================="
echo "=                            START RUNNING WITH:                             ="
echo "=============================================================================="
echo "                NYCPU=$NYCPU      NZCPU=$NZCPU    NX=$NX     NY=$NY     NZ=$NZ"
echo ""




if [[  ( "$machine" == "0" ) || ( "$machine" == "1" ) ]]; then
cd ./set_run
./go.sh
cd ..
fi
