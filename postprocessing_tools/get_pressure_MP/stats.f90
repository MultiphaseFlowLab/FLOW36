subroutine stats(step)

use commondata
use fields

integer :: step
integer :: i,j,k
character(len=40) :: namefile

pm=0.0d0
prms=0.0d0

do j=1,ny
 do k=1,nz
  do i=1,nx
   pm(k)=pm(k)+press(i,k,j)
  enddo
 enddo
enddo
pm=pm/dble(nx*ny)


do j=1,ny
 do k=1,nz
  do i=1,nx
   prms(k)=prms(k)+(press(i,k,j)-pm(k))**2
  enddo
 enddo
enddo
prms=prms/dble(nx*ny)
prms=(prms)**0.5


write(namefile,'(a,i8.8,a)') './output/pstats_',step,'.dat'
open(666,file=namefile,status='new',form='formatted')
write(666,'(3(a16))') 'z^+','p mean','p rms'
! to get from 0 to 2*Re
do k=nz,1,-1
 write(666,'(3(e16.5))') z(k),pm(k),prms(k)
enddo
close(666,status='keep')

return
end subroutine
