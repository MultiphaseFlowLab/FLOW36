program get_interface

use mpi
use commondata
implicit none

integer :: ierr,i


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)


allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

allocate(w(nx,nz,ny))
allocate(phi(nx,nz,ny))

call read_grid

open(55,file='./output/int_pos.dat',status='new',form='formatted')
write(55,'(3(a16))')  'time','z_max','z_min'
close(55,status='keep')

do i=nstart,nend,ndump
  write(*,*) 'Step ',i,' of ',nend
  call find_int(i)
enddo


deallocate(x,y,z)
deallocate(w,phi)

return
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_int(nstep)

use commondata

integer :: nstep,i,j,kpos(1)

double precision :: zpos2(nx,ny),zpos(ny)

call read_fields(nstep)


! get interface height
do j=1,ny
  do i=1,nx
    kpos=minloc(dabs(phi(i,:,j)))
    zpos2(i,j)=z(kpos(1))
  enddo
enddo

zpos=(zpos2(1,:)+zpos2(2,:))/2.0d0
! outer units
zpos=zpos/re-1.0d0

!write(*,*) dt,dt*dble(nstep)

open(55,file='./output/int_pos.dat',status='old',form='formatted',access='append')
write(55,'(3(e16.8))') dt*dble(nstep),maxval(zpos),minval(zpos)
close(55,status='keep')

return 
end

