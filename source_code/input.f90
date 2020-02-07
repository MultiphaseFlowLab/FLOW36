! type of simulation
0                     ! restart flag, if 1 is a restart
0               ! number of iteration used as a restart field
1                ! initial conditions
! grid parameters
256                         ! number of points in x
256                         ! number of points in y
257                         ! number of points in z
! simulation parameters
10.0                           ! Reynolds number
0.2                      ! Courant number
-1.0                          ! mean pressure gradient x direction: Delta p/ Lx = (p_out-p_in)/Lx
0.0                          ! mean pressure gradient y direction: Delta p/ Ly
0                         ! 1: CPI activated, 0: CPI deactivated
30.6                         ! Power reynolds number
! domain size
4.0                           ! Lx
2.0                           ! Ly
!
0                          ! starting timestep
5                            ! end timestep
500                          ! solution dumping frequency in physical space (if lower than 1 never save solution)
-1                          ! solution dumping frequency in spectral space (if lower than 1 never save solution)
500                        ! solution dumping frequency in spectral space (used only to save temp files in case the simulation crashes, those files are not kept)
-1                      ! satistics dumping frequency (if lower than 1 never save statistics)
0                     ! time step for starting statistics calculation
1.e-3                         ! dt
! boundary conditions 0: no-slip, 1: free-slip
0                      ! boundary condition at z=1
0                     ! boundary condition at z=-1
! phase field variables
0                    ! 1: phase field activated, 0: no phase field
0                   ! 0: Standard model, 1: profile-corrected, 2: flux-corrected
2.5                       ! Value of lambda to correct phi
1                     ! 1: matched densities, 0: otherwise
1.0                       ! density ratio phase +1 over phase -1
1                     ! 1: matched viscosities, 0: otherwise
1.0                       ! dynamic viscosity ratio phase +1 over phase -1
0                   ! non-newtonian fluid phi=+1
0.1                     ! ratio between the viscosity at inf and 0 shear rate (Non-newtonian)
0.9                     ! exponent for the non-newtonian phase (phi=+1)
1.0                     ! Weber number
0.02                      ! Cahn number
100.0                    ! Peclet number
0.1                     ! Froud number
0             ! initial conditions on phase variable
-2                      ! direction of gravity
0                     ! type of buoyancy considered
0                       ! 1: body force activated, 0: body force deactivated
4.0                      ! coefficient of body force
2                   ! direction of body force
0                        ! 1: electric force activated, 0: electric force deactivated
1.0                       ! Stuart number
! surfactant variables
0                  ! 1: surfactant activated, 0 : no surfactant
100.0                      ! Peclet number for surfactant
0.117                        ! Ex number (surfactant)
1.35                        ! Pi number (surfactant)
0.5                  ! Surface elasticity parameter
2             ! initial condition for the surfactant
! temperature variables
0                 ! 1: temperature activated, 0 : temperature deactivated
1000.0                    ! Rayleigh number
1.0                     ! Prandtl number
1.0                       ! boundary condition at z=-1
0.0                       ! boundary condition at z=-1
1.0                       ! boundary condition at z=-1
1.0                       ! boundary condition at z=+1
0.0                       ! boundary condition at z=+1
-1.0                       ! boundary condition at z=+1
0           ! initial condition on temperature
! Lagrangian particle tracking variables
0                    ! 1: particles activated, 0: no particles
0                  ! total number of particles
0                        ! number of sets of particles
3                   ! particle initial condition position
0                   ! particle initial condition velocity
1000                    ! particle saving frequency
10                 ! number of subiteration for particle tracking
