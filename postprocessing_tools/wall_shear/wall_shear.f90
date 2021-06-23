subroutine wall_shear(nstep)

use commondata
use grid
use sim_par
use velocity
use phase_field
use wavenumber

integer :: nstep
integer :: i,j

double precision, dimension(nx,nz,ny) :: tauw
double precision, dimension(nx/2+1,nz,ny,2) :: tauwc
double precision :: tau_u,tau_b,ub,um(nz)

character(len=80) :: path
character(len=8) :: numfile
character(len=100) :: namefile

! write(*,*) 'Reading step ',nstep

path='../results/'
!path='../FLOW36SURF/set_run/results/'

write(numfile,'(i8.8)') nstep
! read u
namefile=trim(path)//'u_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) u
close(666,status='keep')

call phys_to_spectral(u,uc,0)

call dz(uc,tauwc)

call spectral_to_phys(tauwc,tauw,0)

! calculate viscous stress at the wall
! tau_u at z=+1, tau_b at z=-1

tau_u=0.0d0
tau_b=0.0d0

do j=1,ny
  do i=1,nx
    tau_u=tau_u+tauw(i,1,j)
    tau_b=tau_b+tauw(i,nz,j)
  enddo
enddo

! normalize by Re (from non dimensionalization)
tau_u=tau_u/(dble(nx*ny)*Re)
tau_b=tau_b/(dble(nx*ny)*Re)


um=0.0d0
do j=1,ny
  do k=1,nz
    do i=1,nx
      um(k)=um(k)+u(i,k,j)
    enddo
  enddo
enddo
um=um/dble(nx*ny)

ub=0.0d0
do k=2,nz-1
  ub=ub+um(k)*(z(k-1)-z(k+1))/2.0d0
enddo

! divide by channel height
ub=ub/2.0d0

write(*,'(i16,4(es16.5))') nstep,nstep*dt*re,tau_b,tau_u,ub


open(54,file='./output/wall_shear.dat',form='formatted',status='old',position='append')
write(54,'(i16,4(es16.5))') nstep,nstep*dt*re,tau_b,tau_u,ub
close(54,status='keep')



return
end
