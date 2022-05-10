program get_interface

use mpi
use commondata
implicit none

integer :: i

call read_input

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(phi(nx,nz,ny))

call read_grid

do i=nstart,nend,ndump
  write(*,*) 'Step ',i,' of ',nend
  call capture_interface(i)
enddo


deallocate(x,y,z)
deallocate(phi)

return
end program








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine capture_interface(nstep)

use commondata
integer :: nstep,i,j,k
double precision :: inte(nx,ny,2)
double precision :: pro,mang,mean(2)
character(len=40) :: namefile

call read_fields(nstep)

do i=1,nx
  do j=1,ny
    !! lower interface
    do k=2,nz/2-1
      pro=phi(i,k-1,j)*phi(i,k,j)
      if (pro .lt. 0.0d0) then
        mang=((phi(i,k,j)-phi(i,k-1,j))/(z(k)-z(k-1)))
        inte(i,j,1)=-phi(i,k-1,j)/mang + z(k-1)
      endif
    enddo
    !! upper interface
    do k=nz/2-1,nz-1
      pro=phi(i,k-1,j)*phi(i,k,j)
      if (pro .lt. 0.0d0) then
        mang=((phi(i,k,j)-phi(i,k-1,j))/(z(k)-z(k-1)))
        inte(i,j,2)=-phi(i,k-1,j)/mang + z(k-1)
      endif
    enddo
  enddo
enddo

! Compute the Mean Value of the Interface
mean=0.0d0
do i=1,nx
  do j=1,ny
    do k=1,2
      mean(k)=mean(k)+inte(i,j,k)/dble(nx*ny)
    enddo
  enddo
enddo

do i=1,nx
  do j=1,ny
    do k=1,2
      inte(i,j,k)=inte(i,j,k)-mean(k)
    enddo
  enddo
enddo


write(namefile,'(a,i8.8,a)') './output/interface1_',nstep,'.dat'

open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) inte(:,:,1)
close(55)

write(namefile,'(a,i8.8,a)') './output/interface2_',nstep,'.dat'

open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) inte(:,:,2)
close(55)

return
end
