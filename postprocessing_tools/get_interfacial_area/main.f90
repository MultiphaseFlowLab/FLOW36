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

dx=6.28318530717959d0/real(nx)
dy=3.14159265358979d0/real(ny)
ch=0.01d0
re=180
dt=0.5e-4

call read_grid



! create output file
open(42,file='./output/area.dat',status='new',form='formatted')
 write(42,'(3(a16,2x))') 'iteration','t^+','area'
close(42,status='keep')


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
double precision :: area
character(len=40) :: namefile

call read_fields(nstep)

do i=1,nx
  do j=1,ny
    do k=2,nz-1
      ! look for points between -0.9 and 0.9
      if (abs(phi(i,k,j)) .le. 0.9d0) then
        area=area + dx*dy*0.5d0*abs(z(k-1)-z(k+1))
      endif
    enddo
  enddo
enddo


area=area/4.1d0/ch

write(*,*) "area", area

open(42,file='./output/area.dat',access='append',form='formatted',status='old')
write(42,'(i16,2x,es16.6,2x,es16.6)') nstep,dble(nstep)*dt*re,area
close(42,status='keep')


return
end
