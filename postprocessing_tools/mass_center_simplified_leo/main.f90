program mass_center

use commondata
use velocity
implicit none

integer :: ierr,i

call read_input

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

allocate(phi(nx,nz,ny))
!allocate(phic(nx/2+1,nz,ny,2))
allocate(phix(nx,nz,ny))
allocate(phiy(nx,nz,ny))
allocate(phiz(nx,nz,ny))
! allocate(div(nx,nz,ny))
! allocate(s1(nx,nz,ny))
! allocate(s2(nx,nz,ny))
! allocate(s3(nx,nz,ny))

call read_grid


! create output file
open(42,file='./output/drop_count.dat',status='new',form='formatted')
 write(42,'(3(a16,2x))') 'iteration','t^+','drop count'
close(42,status='keep')

open(66,file='./output/mass_center.dat',form='formatted',status='new')
 write(66,'(a8,a6,6(a16))') 'nstep','drop','xg','yg','zg','Volume','Area','(wall units)'
close(66,status='keep')

!open(67,file='./output/inertia_tensor.dat',form='formatted',status='new')
! write(67,'(2(a6),6(a16))') 'nstep','drop','a11','a12','a13','a22','a23','a33'
!close(67,status='keep')

!open(68,file='./output/eigen_problem.dat',form='formatted',status='new')
! write(68,'(2(a6),4(a16))') 'nstep','drop','eigenvalue','eigenvector x','eigenvector y','eigenvector z'
!close(68,status='keep')

do i=nstart,nend,dump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8)') 'Step ',i,' out of ',nend
 call read_fields(i)
enddo




deallocate(x)
deallocate(y)
deallocate(z)
deallocate(phi)


end program mass_center
