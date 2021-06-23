program get_pressure

use MPI
use mpiIO
use commondata
use vars
use sim_parameter
use wavenumbers
use pdf_calc

integer :: i
!integer :: k(1)
character(len=40) :: namefile
!double precision, allocatable :: samples(:,:)
!double precision :: max_pdf_p,min_pdf_p,max_pdf_t,min_pdf_t,delta_p,delta_t

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)


call read_input

call define_sizes

call wavenumber

call create_plan

! data array initialization
allocate(u(nx,fpz,fpy))
allocate(v(nx,fpz,fpy))
allocate(w(nx,fpz,fpy))
allocate(uc(spx,nz,spy,2))
allocate(vc(spx,nz,spy,2))
allocate(wc(spx,nz,spy,2))
allocate(tke(nx,fpz,fpy))
allocate(tkec(spx,nz,spy,2))
allocate(phi(nx,fpz,fpy))

if(phiflag.eq.1)then
 allocate(phic(spx,nz,spy,2))
endif


allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

write(namefile,'(a)') '../results/x.dat'
open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
read(1) x
close(1,status='keep')

write(namefile,'(a)') '../results/y.dat'
open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
read(2) y
close(2,status='keep')

write(namefile,'(a)') '../results/z.dat'
open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
read(3) z
close(3,status='keep')

!fhandle=313
!open(fhandle,file='./output/list_val.dat',status='new',form='formatted')

total=0
do i=nstart,nend,dump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8)') 'Iteration ',i,' of ',nend
 call load_files(i)
enddo

!close(fhandle,status='keep')

deallocate(u)
deallocate(uc)
deallocate(v)
deallocate(vc)
deallocate(w)
deallocate(wc)
deallocate(tke)
deallocate(tkec)
deallocate(x)
deallocate(y)
deallocate(z)

if(phiflag.eq.1)then
  deallocate(phi)
  deallocate(phic)
endif

deallocate(kx)
deallocate(ky)
deallocate(k2)


!allocate(samples(total,2))
!! read into samples, uper and upar
!open(111,file='./output/list_val.dat',status='old',form='formatted')

!do i=1,total
!  read(111,'(2(e16.6))') samples(i,:)
!enddo

!close(111,status='keep')

!! calculate pdf
!max_pdf_p=maxval(samples(:,1))
!min_pdf_p=minval(samples(:,1))

!max_pdf_t=maxval(samples(:,2))
!min_pdf_t=minval(samples(:,2))

!write(*,*) max_pdf_p,min_pdf_p
!write(*,*) max_pdf_t,min_pdf_t

!delta_p=(max_pdf_p-min_pdf_p)/dble(n_samples-3)
!pdf_p(1,1)=min_pdf_p-delta_p
!do i=2,n_samples
!  pdf_p(i,1)=pdf_p(i-1,1)+delta_p
!enddo

!delta_t=(max_pdf_t-min_pdf_t)/dble(n_samples-3)
!pdf_t(1,1)=min_pdf_t-delta_t
!do i=2,n_samples
!  pdf_t(i,1)=pdf_t(i-1,1)+delta_t
!enddo
!! lowest value of tangential velocity is 0 (PDF of tang velocity modulus)
!if(pdf_t(1,1).lt.0.0d0) pdf_t(1,1)=0.0d0


!! populate pdf
!pdf_p(:,2)=0.0d0
!pdf_t(:,2)=0.0d0
!do i=1,total
!  ! write(*,*) 'Populating pdf, step ',i,' out of ',total
!  k=minloc(abs(samples(i,1)-pdf_p(:,1)))
!  pdf_p(k(1),2)=pdf_p(k(1),2)+1.0d0
!  k=minloc(abs(samples(i,2)-pdf_t(:,1)))
!  pdf_t(k(1),2)=pdf_t(k(1),2)+1.0d0
!enddo

!deallocate(samples)

!! normalize pdf
!pdf_p(:,2)=pdf_p(:,2)/(delta_p*sum(pdf_p(:,2)))
!pdf_t(:,2)=pdf_t(:,2)/(delta_t*sum(pdf_t(:,2)))

!! write output
!open(111,file='./output/pdf_perp.dat',form='formatted',status='new')
!do i=1,n_samples
!  write(111,'(2(e16.6))') pdf_p(i,:)
!enddo
!close(111,status='keep')

!open(112,file='./output/pdf_tang.dat',form='formatted',status='new')
!do i=1,n_samples
!  write(112,'(2(e16.6))') pdf_t(i,:)
!enddo
!close(112,status='keep')



! destroy MPI I/O derived datatype
call mpi_type_free(ftype,ierr)
call mpi_type_free(stype,ierr)

! destroy cartesian communicator
call mpi_comm_free(cart_comm,ierr)

call destroy_plan

call mpi_finalize(ierr)

return
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine define_sizes

use MPI
use mpiIO
use commondata

integer, dimension(3) :: g_size,s_size
integer :: fysize(1),cxsize(1)
integer,allocatable :: cysize(:),fzsize(:)
integer :: dims(2)
integer :: l_comm

logical :: periodic(2),reorder

! creation of cartesian communicator
! create cartesian communicator
dims(1)=nzcpu
dims(2)=nycpu
periodic(1)=.true.
periodic(2)=.true.
reorder=.false.

call mpi_cart_create(mpi_comm_world,2,dims,periodic,reorder,cart_comm,ierr)

allocate(cysize(nzcpu))
allocate(fzsize(nzcpu))
! create derived datatype used in MPI I/O and commit it
call mpi_cart_sub(cart_comm,[.false.,.true.],l_comm,ierr)
call mpi_allgather(fpy,1,mpi_integer,fysize,1,mpi_integer,l_comm,ierr)
call mpi_allgather(spx,1,mpi_integer,cxsize,1,mpi_integer,l_comm,ierr)
call mpi_comm_free(l_comm,ierr)

call mpi_cart_sub(cart_comm,[.true.,.false.],l_comm,ierr)
call mpi_allgather(fpz,1,mpi_integer,fzsize,1,mpi_integer,l_comm,ierr)
call mpi_allgather(spy,1,mpi_integer,cysize,1,mpi_integer,l_comm,ierr)
call mpi_comm_free(l_comm,ierr)

! creation of datatypes for MPI IO
g_size=[nx, nz, ny]
s_size=[nx, fpz, fpy]
fstart=[0, 0, 0]
cstart=[0, 0, 0]
do i=1,mod(rank,nycpu)
  fstart(3)=fstart(3)+fysize(i)
  cstart(1)=cstart(1)+cxsize(i)
enddo
do i=1,floor(real(rank)/real(nycpu))
  fstart(2)=fstart(2)+fzsize(i)
  cstart(3)=cstart(3)+cysize(i)
enddo
deallocate(fzsize)
deallocate(cysize)
! physical space saving
call mpi_type_create_subarray(3,g_size,s_size,fstart,mpi_order_fortran, &
 &     mpi_double_precision,ftype,ierr)

call mpi_type_commit(ftype,ierr)

! spectral space saving
call mpi_type_create_subarray(4,[nx/2+1,nz,ny,2],[spx,nz,spy,2], &
 &     [cstart(1),cstart(2),cstart(3),0],mpi_order_fortran, &
 &     mpi_double_precision,stype,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wavenumber

use commondata
use wavenumbers

allocate(kx(nx/2+1))
allocate(ky(ny))
allocate(k2(nx/2+1,ny))


kx(1)=0.0d0
do i=2,nx/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do i=2,ny/2+1
  ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
  ky(i)=dble(i-1)*2.0d0*pi/yl
enddo

do j=1,ny
  do i=1,nx/2+1
    k2(i,j)=kx(i)*kx(i)+ky(j)*ky(j)
  enddo
enddo


return
end
