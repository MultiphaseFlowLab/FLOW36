subroutine get_distribution(step)

use commondata
use velocity

integer :: step
integer :: i,j,k

double precision, dimension(nz) :: dz,drop
double precision :: tot

character(len=8) :: nstep

dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do k=2,nz-1
 dz(k)=0.5d0*(z(k-1)-z(k+1))
! write(*,*) dz(k)
enddo

drop=0
do j=1,ny
 do k=1,nz
  do i=1,nx
   if(phi(i,k,j).ge.0.0d0) drop(k)=drop(k)+phi(i,k,j)
  enddo
 enddo
enddo


! normalization
tot=0.0d0
do k=1,nz
 tot=tot+drop(k)*dz(k)
enddo
drop=drop/tot

! write output
write(nstep,'(i8.8)') step
open(456,file='./output/dist_'//nstep//'.dat',status='new',form='formatted')
write(456,'(2(a16))') 'z','PDF phi'
do k=1,nz
 write(456,'(2(es16.6))') z(k),drop(k)
enddo
close(456,status='keep')


return
end

