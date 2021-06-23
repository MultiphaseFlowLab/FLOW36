subroutine read_input

use commondata
use sim_parameter
use paraview_utils

integer :: ndump,sdump,grav_dir,init,finish,step,ry,rz
double precision :: Lx,Ly

folder='../results'

nycpu=1
nzcpu=ntask

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
read(66,*) !gradpx
read(66,*) !gradpy
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
read(66,*) !dt
read(66,*)
read(66,*) !bc_up
read(66,*) !bc_low
read(66,*)
read(66,'(i5)') phiflag
read(66,*) ! corrected model
read(66,*) ! lambda
read(66,*) !match_dens
read(66,'(f16.6)') rhor
read(66,*) !match_visc
read(66,'(f16.6)') visr
read(66,'(f16.6)') we
read(66,'(f16.6)') ch
read(66,'(f16.6)') pe
read(66,'(f16.6)') fr
read(66,*) !in_cond_phi
read(66,'(i5)') grav_dir
read(66,'(i5)') b_type
read(66,*)
!read(66,'(i8)') psiflag
!read(66,'(f16.6)') Pe_psi
!read(66,'(f16.6)') Ex
!read(66,'(f16.6)') P_i
!read(66,'(f16.6)') El
!read(66,*) !in_cond_psi
psiflag=0

close(66,status='keep')

xl=Lx*pi
yl=Ly*pi


open(66,file='./input.f90',form='formatted',status='old',action='read')

read(66,'(i8)') init
read(66,'(i8)') finish
read(66,'(i8)') step
read(66,*)
read(66,'(i8)') spectral
read(66,*)
read(66,'(i8)') out_spectral
read(66,'(i8)') generate_paraview
read(66,'(i8)') oldstat
read(66,*)
read(66,'(i8)') x_start
read(66,'(i8)') x_end
read(66,'(i8)') dnx
read(66,*)
read(66,'(i8)') y_start
read(66,'(i8)') y_end
read(66,'(i8)') dny
read(66,*)
read(66,'(i8)') z_start
read(66,'(i8)') z_end
read(66,'(i8)') dnz

close(66,status='keep')

begin=nstart
if(init.gt.nstart) nstart=init

if((finish.lt.nend).and.(finish.ge.0)) nend=finish

if(spectral.eq.0)then
 dump=ndump
else
 dump=sdump
endif

if(step.gt.dump) dump=step

if(x_start.lt.1) x_start=1
if(x_end.lt.1) x_end=nx
if(dnx.lt.1) dnx=1

if(y_start.lt.1) y_start=1
if(y_end.lt.1) y_end=ny
if(dny.lt.1) dny=1

if(z_start.lt.1) z_start=1
if(z_end.lt.1) z_end=nz
if(dnz.lt.1) dnz=1

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

! define sizes for parallelization
rz=mod(nz,nzcpu)

fpy=ny

fpz=int((nz-rz)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.rz)then
 fpz=int((nz-rz)/nzcpu)+1
endif


ry=mod(ny,nzcpu)

spx=nx/2+1

spy=int((ny-ry)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.ry)then
 spy=int((ny-ry)/nzcpu)+1
endif


if(rank.eq.0)then
write(*,*)
write(*,'(1x,3(a,i6))') 'Grid size (nx ny nz): ',nx,' x',ny,' x',nz
write(*,'(1x,3(a,f8.2))') 'Domain size (x y z)',xl,' x',yl,' x',2.0d0
write(*,'(1x,a,f8.1)') 'Reynolds: ',Re
write(*,'(1x,3(a,i8))') 'Reading from ',nstart,' to ',nend,' with step ',dump
if(phiflag.eq.1)then
 write(*,*)
 write(*,'(1x,a,f6.4)') 'Density ratio: ',rhor
 write(*,'(1x,a,f6.4)') 'Viscosity ratio: ',visr
 write(*,'(1x,a,f8.4)') 'We: ',We
 write(*,'(1x,a,f8.6)') 'Ch: ',Ch
 write(*,'(1x,a,f8.2)') 'Pe: ',Pe
 write(*,'(1x,a,f8.6)') 'Fr: ',Fr
 if(psiflag.eq.1)then
  write(*,*)
  write(*,'(1x,a,f8.2)') 'Pe_psi: ',Pe_psi
  write(*,'(1x,a,f8.6)') 'Pi: ',P_i
  write(*,'(1x,a,f8.6)') 'E: ',El
  write(*,'(1x,a,f8.6)') 'Ex: ',Ex
 endif
endif
write(*,*)
endif


return
end
