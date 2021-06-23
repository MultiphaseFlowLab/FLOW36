subroutine read_input

use commondata
use sim_par
use phase_field
use stats
use surfactant
use temperature

integer :: grav_dir,match_dens,match_visc
double precision :: Lx,Ly
integer :: start_temp, end_temp, delta_temp


open(unit=66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')

read(66,*)
read(66,*) !restart
read(66,*) !nt_restart
read(66,*) !in_cond
read(66,*)
read(66,'(i8)') nx
read(66,'(i8)') ny
read(66,'(i8)') nz
read(66,*)
read(66,'(f16.6)') Re
read(66,*) !Co
read(66,'(f16.6)') gradpx
read(66,'(f16.6)') gradpy
read(66,*)
read(66,'(f16.6)') Lx
read(66,'(f16.6)') Ly
read(66,*)
read(66,'(i8)') nstart
read(66,'(i8)') nend
read(66,'(i8)') ndump
read(66,'(i8)') sdump
read(66,*) !dump_failure
read(66,*) !stat_dump
read(66,*) !stat_start
read(66,'(es16.5)') dt
read(66,*)
read(66,'(i5)') bc_up
read(66,'(i5)') bc_low
read(66,*)
read(66,'(i5)') phi_flag
read(66,*) ! correction flag
read(66,*) ! lambda coefficient
read(66,'(i5)') match_dens
read(66,'(f16.6)') rhor
read(66,'(i5)') match_visc
read(66,'(f16.6)') visr
read(66,*) !non newtonian
read(66,*) !non newtonian
read(66,*) !non newtonian
read(66,'(f16.6)') we
read(66,'(f16.6)') ch
read(66,'(f16.6)') pe
read(66,'(f16.6)') fr
read(66,'(i5)') in_cond_phi
read(66,'(i5)') grav_dir
read(66,'(i5)') b_type
read(66,*) ! body force flag
read(66,*) ! body force coefficient
read(66,*) ! body force direction
read(66,*) ! electric force flag
read(66,*) ! Stuart number
!read(66,*)
!read(66,'(i8)') psi_flag
!read(66,'(f16.6)') Pe_psi
!read(66,'(f16.6)') Ex
!read(66,'(f16.6)') P_i
!read(66,'(f16.6)') El
!read(66,'(i8)') in_cond_psi
!read(66,*)
!read(66,'(i8)') theta_flag
!read(66,'(f16.6)') Ra
!read(66,'(f16.6)') Pr
!read(66,*) !p_theta(1)
!read(66,*) !q_theta(1)
!read(66,*) !r_theta(1)
!read(66,*) !p_theta(2)
!read(66,*) !q_theta(2)
!read(66,*) !r_theta(2)
!read(66,*) !in_cond_theta


xl=Lx*pi
yl=Ly*pi


if(grav_dir.eq.1)then
grav=[1.0d0,0.0d0,0.0d0]
elseif(grav_dir.eq.-1)then
grav=[-1.0d0,0.0d0,0.0d0]
elseif(grav_dir.eq.2)then
grav=[0.0d0,1.0d0,0.0d0]
elseif(grav_dir.eq.-2)then
grav=[0.0d0,-1.0d0,0.0d0]
elseif(grav_dir.eq.3)then
grav=[0.0d0,0.0d0,1.0d0]
elseif(grav_dir.eq.-3)then
grav=[0.0d0,0.0d0,-1.0d0]
endif

!call print_start

close(66)

open(70, file = 'input.f90', form = 'formatted', status= 'old')
read(70, '(i8)') start_temp
read(70, '(i8)') end_temp
read(70, '(i8)') delta_temp
close(70, status = 'keep')

if ((start_temp.ne.-1).and.(start_temp.gt.nstart)) nstart = start_temp
if ((end_temp.ne.-1).and.(end_temp.lt.nend)) nend = end_temp
if ((delta_temp.ne.-1).and.(delta_temp.gt.ndump)) ndump = delta_temp

call print_start

return
end
