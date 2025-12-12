subroutine read_input

use commondata

double precision :: Lx, Ly

integer :: begin,finish,step,ndump,sdump,grav_dir

 open(unit=66,file='./input_curvature.f90',form='formatted',status='old',action='read')

 read(66,'(i8)') begin
 read(66,'(i8)') finish
 read(66,'(i8)') step
 read(66,*)
 read(66,'(i5)') spectral

 close(66)




 open(66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')

 read(66,*)
 read(66,*) !restart
 read(66,*) !nt_restart
 read(66,*) !in_cond
 read(66,*)
 read(66,'(i5)') nx
 read(66,'(i5)') ny
 read(66,'(i5)') nz
 read(66,*)
 read(66,'(f16.6)') Re
 read(66,*) !Co
 read(66,'(f16.6)') gradpx
 read(66,'(f16.6)') gradpy
 read(66,*) !CPI flag
 read(66,*) !power Re
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
 read(66,'(es16.6)') dt
 read(66,*)
 read(66,*) !bc_up
 read(66,*) !bc_low
 read(66,*)
 read(66,'(i5)') phiflag
 read(66,*) ! phicorflag
 read(66,*) ! lambda
 read(66,*) !match_dens
 read(66,'(f16.6)') rhor
 read(66,*) !match_visc
 read(66,'(f16.6)') visr
 read(66,*) !non newtonian
 read(66,*) !non newtonian
 read(66,*) !non newtonian
 read(66,'(f16.6)') we
 read(66,'(f16.6)') ch
 read(66,'(f16.6)') pe
 read(66,'(f16.6)') fr
 read(66,*) !in_cond_phi
 read(66,'(i5)') grav_dir
 read(66,'(i5)') b_type
 read(66,*) !b_type
 read(66,*) !body_flag
 read(66,*) !body_c
 read(66,*) !body_dir
 read(66,*) !ele_flag
 read(66,*) !stuart
 read(66,*)
! read(66,'(i8)') psiflag
! read(66,'(f16.6)') Pe_psi
! read(66,'(f16.6)') Ex
! read(66,'(f16.6)') P_i
! read(66,'(f16.6)') El
! read(66,*) !in_cond_psi
psiflag=0
 close(66,status='keep')


 xl=Lx*pi
 yl=Ly*pi

 if(begin.gt.nstart) nstart=begin

 if((finish.lt.nend).and.(finish.ge.0)) nend=finish

 if(spectral.eq.0)then
  dump=ndump
 else
  dump=sdump
 endif

 if(step.gt.dump) dump=step

 ! gravity unity vector
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

return
end


subroutine print_start

use commondata

write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,'(1x,a)') '                 starting curvature calculation                       '
write(*,'(1x,a,3(i5,a))') 'Nx * Ny * Nz = ',nx,' * ',ny,' * ',nz,' '
write(*,'(1x,a,3(f8.2,a))') 'Lx * Ly * Lz = (',xl,' * ',yl,' * ',2.0d0,')*Re'
write(*,'(1x,a,f8.2)') 'Re = ',re
if(spectral.eq.1)then
 write(*,*) 'Reading fields in modal space'
else
 write(*,*) 'Reading fields in physical space'
endif
write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',dump
write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,*)

return
end
