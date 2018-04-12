subroutine read_input

use mpi
use commondata
use sim_par
use phase_field
use stats
use surfactant
use temperature

integer :: grav_dir,match_dens,match_visc
double precision :: Lx,Ly

 open(unit=66,file='./sc_compiled/input.f90',form='formatted',status='old',action='read')

 read(66,*)
 read(66,'(i5)') restart
 read(66,'(i8)') nt_restart
 read(66,'(i5)') in_cond
 read(66,*)
 read(66,*) ! skip grid info lines, already in modules
 read(66,*)
 read(66,*)
 read(66,*)
 read(66,'(f16.6)') Re
 read(66,'(f16.6)') Co
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
 read(66,'(i8)') dump_failure
 read(66,'(i8)') stat_dump
 read(66,'(i8)') stat_start
 read(66,'(es16.5)') dt
 read(66,*)
 read(66,'(i5)') bc_up
 read(66,'(i5)') bc_low
 read(66,*)
 read(66,'(i5)') phi_flag
 read(66,'(i5)') phicor_flag
 read(66,'(f16.6)') lamphi
 read(66,'(i5)') match_dens
 read(66,'(f16.6)') rhor
 read(66,'(i5)') match_visc
 read(66,'(f16.6)') visr
 read(66,'(f16.6)') we
 read(66,'(f16.6)') ch
 read(66,'(f16.6)') pe
 read(66,'(f16.6)') fr
 read(66,'(i5)') in_cond_phi
 read(66,'(i5)') grav_dir
 read(66,'(i5)') b_type
 read(66,*)
 read(66,'(i8)') psi_flag
 read(66,'(f16.6)') Pe_psi
 read(66,'(f16.6)') Ex
 read(66,'(f16.6)') P_i
 read(66,'(f16.6)') El
 read(66,'(i8)') in_cond_psi
 read(66,*)
 read(66,'(i8)') theta_flag
 read(66,'(f16.6)') Ra
 read(66,'(f16.6)') Pr
 read(66,'(f16.6)') p_theta(1)
 read(66,'(f16.6)') q_theta(1)
 read(66,'(f16.6)') r_theta(1)
 read(66,'(f16.6)') p_theta(2)
 read(66,'(f16.6)') q_theta(2)
 read(66,'(f16.6)') r_theta(2)
 read(66,'(i8)') in_cond_theta

 r_theta=r_theta*dble(nx*ny)

! if it is a simulation restart, use old flow field data and start from nt_restart
 if(restart.eq.1)then
  nstart=nt_restart
  in_cond=3
  in_cond_phi=1
  in_cond_psi=6
  in_cond_theta=1
 endif

 if((match_dens.eq.1).and.(abs(rhor-1.0d0).lt.0.00000001))then
   if(rank.eq.0) write(*,*) 'Matched densities'
 elseif((match_dens.eq.0).and.(rhor.lt.1.0d0).and.(abs(rhor-1.0d0).gt.0.00000001))then
   if(rank.eq.0) write(*,*) 'Non-matched densities, rhor < 1'
 elseif((match_dens.eq.2).and.(rhor.gt.1.0d0).and.(abs(rhor-1.0d0).gt.0.00000001))then
   if(rank.eq.0) write(*,*) 'Non-matched densities, rhor > 1'
 else
  if(rank.eq.0) write(*,*) 'Error in input parameters: non coherent input for matchedrho and rhor'
  stop
 endif

 ! elseif((match_dens.ne.1).and.(abs(rhor-1.0d0).gt.0.00000001))then
 !   if(rank.eq.0) write(*,*) 'Non-matched densities'
 ! else
 !   if(rank.eq.0) write(*,*) 'Error in input parameters: non coherent input for matchedrho and rhor'
 !   stop
 ! endif


 if((match_visc.eq.1).and.(abs(visr-1.0d0).lt.0.00000001))then
   if(rank.eq.0) write(*,*) 'Matched viscosities'
 elseif((match_visc.ne.1).and.(abs(visr-1.0d0).gt.0.00000001))then
   if(rank.eq.0) write(*,*) 'Non-matched viscosities'
 else
   if(rank.eq.0) write(*,*) 'Error in input parameters: non coherent input for matchedvis and visr'
   stop
 endif

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

 if(rank.eq.0) call print_start

 close(66)


return
end
