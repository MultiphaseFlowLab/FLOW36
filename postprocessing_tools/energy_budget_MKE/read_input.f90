subroutine read_input

use mpi
use commondata

double precision :: Lx,Ly
integer :: ndump,sdump,pos,eon,grav_dir
logical :: control
character(len=200) :: rline,fstring

open(unit=66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')
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
 read(66,'(es16.5)') dt
 read(66,*)
 read(66,*) !bc_up
 read(66,*) !bc_low
 read(66,*)
 read(66,'(i5)') phi_flag
 read(66,*) !phicor_flag
 read(66,*) !lamphi
 read(66,'(i5)') matchedrho
 read(66,'(f16.6)') rhor
 read(66,'(i5)') matchedvis
 read(66,'(f16.6)') visr
 read(66,*) !non_newtonian
 read(66,*) !muinfmuzero
 read(66,*) !exp_non_new
 read(66,'(f16.6)') we
 read(66,'(f16.6)') ch
 read(66,*) !pe
 read(66,'(f16.6)') fr
 read(66,*) !in_cond_phi
 read(66,'(i4)') grav_dir
 read(66,'(i4)') b_type
 read(66,*) !body_flag
 read(66,*) !body_c
 read(66,*) !body_dir
 read(66,*) !ele_flag
 read(66,*) !stuart
 read(66,*)
 read(66,'(i8)') psi_flag
 read(66,*) !Pe_psi
 read(66,*) !Ex
 read(66,*) !P_i
 read(66,'(f16.6)') betas
 read(66,*) !in_cond_psi
 read(66,*)
 read(66,*) !theta_flag
 read(66,*) !Ra
 read(66,*) !Pr
 read(66,*) !p_theta(1)
 read(66,*) !q_theta(1)
 read(66,*) !r_theta(1)
 read(66,*) !p_theta(2)
 read(66,*) !q_theta(2)
 read(66,*) !r_theta(2)
 read(66,*) !in_cond_theta
close(66,status='keep')

if(phi_flag.eq.0) psi_flag=0

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

! read grid expansion parameters from module file
open(66,file='../sc_compiled/module.f90',form='formatted',status='old',action='read')
control=.false.
do while(control.eqv..false.)
read(66,'(a200)') rline
  if(rline(1:29)==' integer, parameter :: exp_x=') then
   control=.true.
   ! read exp_x
   pos=scan(rline,'=',.false.)
   eon=scan(rline(pos+1:),',',.false.)
   eon=eon+pos-1
   write(fstring,'(a,i3.3,a)') '(i',eon-pos,')'
   read(rline(pos+1:eon),trim(fstring)) expx
   write(rline(1:eon+1),'(a)') ' '
   ! read exp_y
   pos=scan(rline,'=',.false.)
   eon=scan(rline(pos+1:),',',.false.)
   eon=eon+pos-1
   write(fstring,'(a,i3.3,a)') '(i',eon-pos,')'
   read(rline(pos+1:eon),trim(fstring)) expy
   write(rline(1:eon+1),'(a)') ' '
   ! read exp_z
   pos=scan(rline,'=',.false.)
   eon=len(trim(rline))
   write(fstring,'(a,i3.3,a)') '(i',eon-pos,')'
   read(rline(pos+1:eon),trim(fstring)) expz
  endif
enddo
close(66,status='keep')

if(ndump.eq.-1)then
 if(sdump.eq.-1)then
  write(*,*) 'Error in file saving frequency'
  spectral=-1
 else
  delta=sdump
  spectral=1
 endif
else
 delta=ndump
 spectral=0
endif


xl=Lx*pi
yl=Ly*pi

if(psi_flag.eq.0)then
 expx=1
 expy=1
 expz=1
endif

nxfg=expx*nx
nyfg=expy*ny
nzfg=(nz-1)*expz+1


if(rank.eq.0)then
write(*,*) '---------------------------------------------------------------------------'
write(*,'(6(a12))') 'Re','Ch','We','beta_s','eta_r','rho_r'
write(*,'(6(f12.4))') re,ch,we,betas,visr,rhor
write(*,*)
write(*,'(3(a12))') 'Nx','Ny','Nz'
write(*,'(3(i12))') nx,ny,nz
write(*,*)
write(*,'(3(a12))') 'start','end','delta'
write(*,'(3(i12))') nstart,nend,delta
write(*,*)
write(*,'(2(a12))') 'phi flag','psi flag'
write(*,'(2(i12))') phi_flag,psi_flag
write(*,*)
write(*,'(3x,a)') 'Grid expansion, x,y,z'
write(*,'(6x,3(i6))') expx,expy,expz
write(*,*) '---------------------------------------------------------------------------'
endif



return
end
