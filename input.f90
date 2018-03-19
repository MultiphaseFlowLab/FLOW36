! type of simulation
restartflag                     ! restart flag, if 1 is a restart
restart_iteration               ! number of iteration used as a restart field
initialcondition                ! initial conditions
! grid parameters
nxxxxxx                         ! number of points in x
nyyyyyy                         ! number of points in y
nzzzzzz                         ! number of points in z
! simulation parameters
Renum                           ! Reynolds number
courantnum                      ! Courant number
gradpx                          ! mean pressure gradient x direction: Delta p/ Lx = (p_out-p_in)/Lx
gradpy                          ! mean pressure gradient y direction: Delta p/ Ly
! domain size
len_x                           ! Lx
len_y                           ! Ly
!
nstart                          ! starting timestep
nend                            ! end timestep
nfdump                          ! solution dumping frequency in physical space (if lower than 1 never save solution)
nsdump                          ! solution dumping frequency in spectral space (if lower than 1 never save solution)
faildump                        ! solution dumping frequency in spectral space (used only to save temp files in case the simulation crashes, those files are not kept)
stats_dump                      ! satistics dumping frequency (if lower than 1 never save statistics)
stats_start                     ! time step for starting statistics calculation
delta_t                         ! dt
! boundary conditions 0: no-slip, 1: free-slip
bc_upbound                      ! boundary condition at z=1
bc_lowbound                     ! boundary condition at z=-1
! phase field variables
phasephiflag                    ! 1: phase field avtivated, 0: no phase field
matcheddens                     ! 1: matched densities, 0: otherwise
densratio                       ! density ratio phase +1 over phase -1
matchedvisc                     ! 1: matched viscosities, 0: otherwise
viscratio                       ! dynamic viscosity ratio phase +1 over phase -1
webernumber                     ! Weber number
cahnnumber                      ! Cahn number
pecletnumber                    ! Peclet number
froudnumber                     ! Froud number
phinitial_condition             ! initial conditions on phase variable
gravitydir                      ! direction of gravity
gravitytype                     ! type of buoyancy considered
! surfactant variables
surfactantflag                  ! 1: surfactant activated, 0 : no surfactant
surfpeclet                      ! Peclet number for surfactant
exnumber                        ! Ex number (surfactant)
pinumber                        ! Pi number (surfactant)
surfelasticity                  ! Surface elasticity parameter
psinitial_condition             ! initial condition for the surfactant
! temperature variables
temperatureflag                 ! 1: temperature activated, 0 : temperature deactivated
rayleighnumb                    ! Rayleigh number
prandtlnumb                     ! Prandtl number
Aboundary                       ! boundary condition at z=-1
Bboundary                       ! boundary condition at z=-1
Cboundary                       ! boundary condition at z=-1
Dboundary                       ! boundary condition at z=+1
Eboundary                       ! boundary condition at z=+1
Fboundary                       ! boundary condition at z=+1
tempinitial_condition            ! initial condition on temperature
