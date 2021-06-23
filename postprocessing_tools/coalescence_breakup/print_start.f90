subroutine print_start

use commondata
use sim_par
use phase_field
use stats
use surfactant
use temperature

write(*,*)
write(*,*) '-----------------------------------------------------------------------'
write(*,'(1x,a)') 'Starting simulation with following parameters:'
if(restart.eq.1)then
 write(*,'(1x,a,i8)') 'Restarting simulation from iteration:',nt_restart
else
 write(*,'(1x,a)') 'New simulation'
endif
if(in_cond.eq.0)then
 write(*,'(1x,a)') 'Zero velocity field'
elseif(in_cond.eq.1)then
 write(*,'(1x,a)') 'Laminar Poiseuille flow in x direction'
elseif(in_cond.eq.2)then
 write(*,'(1x,a)') 'Laminar Poiseuille flow in y direction'
elseif(in_cond.eq.3)then
 write(*,'(1x,a)') 'Read fields from input data'
elseif(in_cond.eq.4)then
 write(*,'(1x,a)') 'Read fields from input data (serial read)'
elseif(in_cond.eq.5)then
 write(*,'(1x,a)') 'Shear flow y direction'
elseif(in_cond.eq.6)then
 write(*,'(1x,a)') 'Shear flow x direction'
else
 write(*,'(1x,a)') 'Dafuq? Check in_cond input value'
endif

write(*,*)
! write(*,'(1x,a40,i0.1,a,i0.1,a,i0.1)') 'Nycpu * Nzcpu = Ntask : ',nycpu,' x ',nzcpu,' = ',ntask
write(*,'(1x,a40,i0.1,a,i0.1,a,i0.1)') 'Nx * Nz * Ny : ',nx,' x ',nz,' x ',ny
write(*,'(1x,a40,f6.3,a,f6.3)') 'Lx * Ly : ',xl,' x ',yl
write(*,*)
write(*,'(1x,a40,i6,a,i6)') 'Timestep from : ',nstart,' to ',nend
write(*,'(1x,a40,i6)') 'Solution saving frequency (physical) : ',ndump
write(*,'(1x,a40,i6)') 'Solution saving frequency (spectral) : ',sdump
! write(*,'(1x,a40,i6)') 'Statistics saving frequency : ',stat_dump
! write(*,'(1x,a40,i6)') 'Statistics calculation first timestep : ',stat_start
write(*,'(1x,a40,es8.2)') 'dt : ',dt
write(*,*)
write(*,'(1x,a40,f8.1)') 'Re : ',Re
! write(*,'(1x,a40,f8.3)') 'Co : ',Co
write(*,'(1x,a40,3(a,f5.1),a)') 'Mean pressure gradient [x,y,z] : ','[',gradpx,',',gradpy,',',0.0d0,']'
if(bc_up.eq.0)then
 write(*,'(1x,a57)') '@ z=1  : no-slip condition'
elseif(bc_up.eq.1)then
 write(*,'(1x,a59)') '@ z=1  : free-slip condition'
elseif(bc_up.eq.2)then
 write(*,'(1x,a62)') '@ z=1  : shear flow y direction'
elseif(bc_up.eq.3)then
 write(*,'(1x,a62)') '@ z=1  : shear flow x direction'
else
 write(*,'(1x,a46)') '@ z=1  : dafuq?'
endif

if(bc_low.eq.0)then
 write(*,'(1x,a57)') '@ z=-1 : no-slip condition'
elseif(bc_low.eq.1)then
 write(*,'(1x,a59)') '@ z=-1 : free-slip condition'
elseif(bc_low.eq.2)then
 write(*,'(1x,a62)') '@ z=-1  : shear flow y direction'
elseif(bc_low.eq.3)then
 write(*,'(1x,a62)') '@ z=-1  : shear flow x direction'
else
 write(*,'(1x,a46)') '@ z=-1 : dafuq?'
endif

if(phi_flag.eq.1)then
 write(*,'(1x,a)') 'Phase field parameters'
 if(in_cond_phi.eq.0)then
  write(*,'(1x,a)') 'Initialize phi=-1 all over the domain'
 elseif(in_cond_phi.eq.1)then
  write(*,'(1x,a)') 'Read phase from input data'
 elseif(in_cond_phi.eq.2)then
  write(*,'(1x,a)') 'Read phase from input data (serial read)'
 elseif(in_cond_phi.eq.3)then
  write(*,'(1x,a)') 'Initialize 2D drop'
 elseif(in_cond_phi.eq.4)then
  write(*,'(1x,a)') 'Initialize 3D drop'
 elseif(in_cond_phi.eq.5)then
  write(*,'(1x,a)') 'Initialize stratified phase field'
 elseif(in_cond_phi.eq.6)then
  write(*,'(1x,a)') 'Initialize 2D array of 3D drops'
elseif(in_cond_phi.eq.6)then
 write(*,'(1x,a)') 'Initialize 2D array of 3D drops'
elseif(in_cond_phi.eq.7)then
 write(*,'(1x,a)') 'Initialize Drop attached to the bottom wall'
elseif(in_cond_phi.eq.8)then
 write(*,'(1x,a)') 'Initialize Drop kissing'
 else
  write(*,'(1x,a)') 'Dafuq? Check in_cond_phi value'
 endif
 write(*,'(1x,a40,f8.5)') 'We : ',we
 write(*,'(1x,a40,f8.5)') 'Ch : ',ch
 write(*,'(1x,a40,f8.1)') 'Pe : ',pe
 write(*,'(1x,a40,f8.5)') 'Fr : ',fr
 write(*,'(1x,a40,f8.3)') 'Rho_r : ',rhor
 write(*,'(1x,a40,f8.3)') 'Vis_r : ',visr
 write(*,'(1x,a40,3(a,f5.1),a)') 'Gravity versor [x,y,z] : ','[',grav(1),',',grav(3),',',grav(2),']'
 if(b_type.eq.0)then
  write(*,'(1x,a63)') 'Buoyancy type  : no gravity and buoyancy'
 elseif(b_type.eq.1)then
  write(*,'(1x,a45)') 'Buoyancy type  : rho*g'
 elseif(b_type.eq.2)then
  write(*,'(1x,a51)') 'Buoyancy type  : Delta rho*g'
 else
  write(*,'(1x,a46)') 'Buoyancy type  : dafuq?'
 endif

 if(psi_flag.eq.1)then
  write(*,'(1x,a)') 'Surfactant parameters'
  if(in_cond_psi.eq.0)then
   write(*,'(1x,a)') 'Initialize constant surfactant field'
  elseif(in_cond_psi.eq.1)then
   write(*,'(1x,a)') 'Read surfactant field from input data'
  elseif(in_cond_psi.eq.2)then
   write(*,'(1x,a)') 'Initialize equilibrium profile'
  elseif(in_cond_psi.eq.3)then
   write(*,'(1x,a)') 'Initialize equilibrium profile + Y Gradient'
  elseif(in_cond_psi.eq.4)then
  write(*,'(1x,a)') 'Initialize equilibrium profile + Z Gradient'
  elseif(in_cond_psi.eq.5)then
  write(*,'(1x,a)') 'Initialize Diffusion Test'
  else
   write(*,'(1x,a)') 'Dafuq? Check in_cond_psi value'
  endif
  write(*,'(1x,a40,f8.1)') 'Pe : ',pe_psi
  write(*,'(1x,a40,f8.5)') 'Ex : ',Ex
  write(*,'(1x,a40,f8.5)') 'Pi : ',P_i
  write(*,'(1x,a40,f8.5)') 'El : ',El
 endif
endif

! if(theta_flag.eq.1)then
!  write(*,'(1x,a)') 'Temperature parameters'
!  if(in_cond_theta.eq.0)then
!   write(*,'(1x,a)') 'Initialize constant temperature'
!  elseif(in_cond_theta.eq.1)then
!   write(*,'(1x,a)') 'Read temperature from input data'
!  else
!   write(*,'(1x,a)') 'Dafuq? Check in_cond_temp value'
!  endif
!  write(*,'(1x,a)') 'Boundary conditions: A*T+B*dT/dz=C'
!  write(*,'(3(f4.1,a))') p_theta(1),'*T+',q_theta(1),'*dT/dz=',r_theta(1)/dble(nx*ny),' at z=-1'
!  write(*,'(3(f4.1,a))') p_theta(2),'*T+',q_theta(2),'*dT/dz=',r_theta(2)/dble(nx*ny),' at z=+1'
!  if(phi_flag.eq.0)then
!   write(*,'(1x,a40,3(a,f5.1),a)') 'Gravity versor [x,y,z] : ','[',grav(1),',',grav(3),',',grav(2),']'
!  endif
! #define boussinnesq boussinnesqcompflag
! #if boussinnesq == 0
!  write(*,'(1x,a)') 'No gravity term in N-S'
! #elif boussinnesq == 1
!  write(*,'(1x,a)') 'Boussinnesq term in N-S'
! #endif
!  write(*,'(1x,a40,es8.1)') 'Ra : ',Ra
!  write(*,'(1x,a40,f8.3)') 'Pr : ',Pr
! endif

write(*,*) '-----------------------------------------------------------------------'


return
end
