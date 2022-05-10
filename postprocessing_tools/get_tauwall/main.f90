program get_interface

use mpi
use commondata
implicit none

integer :: i

call read_input

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(u(nx,nz,ny))

call read_grid

do i=nstart,nend,ndump
  write(*,*) 'Step ',i,' of ',nend
  call compute_tau(i)
enddo


deallocate(x,y,z)
deallocate(u)

return
end program








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_tau(nstep)

use commondata
integer :: nstep,i,j
double precision :: tauwall(nx,ny,2)
double precision :: dz
character(len=40) :: namefile

call read_fields(nstep)

dz=abs(z(2)-z(1))
!print*,"dz",dz

do i=1,nx
  do j=1,ny
    tauwall(i,j,1)=u(i,2,j)/dz
    tauwall(i,j,2)=u(i,nz-1,j)/dz
  enddo
enddo


write(namefile,'(a,i8.8,a)') './output/tauw1_',nstep,'.dat'

open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) tauwall(:,:,1)
close(55)

write(namefile,'(a,i8.8,a)') './output/tauw2_',nstep,'.dat'

open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) tauwall(:,:,2)
close(55)

return
end 
