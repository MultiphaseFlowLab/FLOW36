module dctz_bwd_module
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dctz_bwd(uin,uout,aliasing)

#define openaccflag openacccompflag
#if openaccflag == 0
use fftw3
#endif
#if openaccflag == 1
use cufft
use openacc
use cufftplans
#endif

implicit none
integer(c_int) :: nsx,nz,npy
integer :: aliasing,i,k,j
real(c_double), dimension(:,:,:,:) :: uin, uout
real(c_double), allocatable :: tin(:),tout(:)
!used only with cuFFT
#if openaccflag == 1
real(c_double), allocatable :: a(:,:,:), b(:,:,:)
complex(c_double_complex), allocatable :: ac(:,:,:), bc(:,:,:)
#endif

! get dimensions
nsx=size(uin,1)
nz=size(uin,2)
npy=size(uin,3)


#if openaccflag == 0
! dealiasing
if(aliasing.eq.1)then
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif

! single dct at a time, try with many dft
allocate(tin(nz),tout(nz))
do i=1,nsx
 do j=1,npy
  do k=1,2
   tin(:)=uin(i,:,j,k)
   tin(nz)=2.0d0*tin(nz)
   call fftw_execute_r2r(plan_z_bwd,tin,tout)
   uout(i,:,j,k)=tout(:)*0.5d0
  enddo
 enddo
enddo
deallocate(tin,tout)
#endif


#if openaccflag == 1
! trick to make DCT faster and in one block (not possible otherwise)
! transpose uin and uout and perform DCT along 1st direction
allocate(a(2*(nz-1),nsx,npy))
allocate(b(2*(nz-1),nsx,npy))
allocate(ac(nz,nsx,npy))
allocate(bc(nz,nsx,npy))
!$acc data copyin(uin) create(a,b,ac,bc) copyout(uout)
!$acc kernels
! dealiasing
if(aliasing.eq.1)then
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
uin(:,nz,:,1)=2.0d0*uin(:,nz,:,1)
uin(:,nz,:,2)=2.0d0*uin(:,nz,:,2)
!$acc end kernels
!$acc kernels
do j=1,npy
 do i=1,nsx
  do k=1,nz-1
   a(k,i,j)=uin(i,k,j,1)
   b(k,i,j)=uin(i,k,j,2)
  enddo
  do k=nz,2*(nz-1)
    a(k,i,j)=uin(i,2*(nz-1)-k+2,j,1)
    b(k,i,j)=uin(i,2*(nz-1)-k+2,j,2)
   enddo
 enddo
enddo
!$acc end kernels
!$acc host_data use_device(a,ac)
 gerr=gerr+cufftExecD2Z(cudaplan_z_bwd,a,ac)
!$acc end host_data
!$acc host_data use_device(b,bc)
gerr=gerr+cufftExecD2Z(cudaplan_z_bwd,b,bc)
!$acc end host_data
!$acc kernels
 do j=1,npy
  do i=1,nsx
   do k=1,nz
    uout(i,k,j,1)=0.5d0*dble(ac(k,i,j))
    uout(i,k,j,2)=0.5d0*dble(bc(k,i,j))
   enddo
  enddo
enddo
!$acc end kernels
!$acc end data
deallocate(a,b,ac,bc)
#endif

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dctz_bwd_fg(uin,uout,aliasing)

#define openaccflag openacccompflag
#if openaccflag == 0
use fftw3
#endif
#if openaccflag == 1
use cufft
use openacc
use cufftplans
#endif

implicit none
integer(c_int) :: nsx,nz,npy
integer :: aliasing,i,k,j
real(c_double), allocatable :: tin(:),tout(:)
real(c_double), dimension(:,:,:,:) :: uin, uout
!used only with cuFFT
#if openaccflag == 1
real(c_double), allocatable :: a(:,:,:), b(:,:,:)
complex(c_double_complex), allocatable :: ac(:,:,:), bc(:,:,:)
#endif

! get dimensions
nsx=size(uin,1)
nz=size(uin,2)
npy=size(uin,3)

#if openaccflag == 0
! dealiasing
if(aliasing.eq.1)then
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif

allocate(tin(nz),tout(nz))
do i=1,nsx
 do j=1,npy
  do k=1,2
   tin(:)=uin(i,:,j,k)
   tin(nz)=2.0d0*tin(nz)
   call fftw_execute_r2r(plan_z_bwd_fg,tin,tout)
   uout(i,:,j,k)=tout(:)*0.5d0
  enddo
 enddo
enddo
deallocate(tin,tout)
#endif


#if openaccflag == 1
! trick to make DCT faster and in one block (see above)
allocate(a(2*(nz-1),nsx,npy))
allocate(b(2*(nz-1),nsx,npy))
allocate(ac(nz,nsx,npy))
allocate(bc(nz,nsx,npy))
!$acc data copyin(uin) create(a,b,ac,bc) copyout(uout)
!$acc kernels
! dealiasing
if(aliasing.eq.1)then
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
uin(:,nz,:,1)=2.0d0*uin(:,nz,:,1)
uin(:,nz,:,2)=2.0d0*uin(:,nz,:,2)
!$acc end kernels
!$acc kernels
do j=1,npy
 do i=1,nsx
  do k=1,nz-1
   a(k,i,j)=uin(i,k,j,1)
   b(k,i,j)=uin(i,k,j,2)
  enddo
  do k=nz,2*(nz-1)
    a(k,i,j)=uin(i,2*(nz-1)-k+2,j,1)
    b(k,i,j)=uin(i,2*(nz-1)-k+2,j,2)
   enddo
 enddo
enddo
!$acc end kernels
!$acc host_data use_device(a,ac)
 gerr=gerr+cufftExecD2Z(cudaplan_z_bwd_fg,a,ac)
!$acc end host_data
!$acc host_data use_device(b,bc)
gerr=gerr+cufftExecD2Z(cudaplan_z_bwd_fg,b,bc)
!$acc end host_data
!$acc kernels
 do j=1,npy
  do i=1,nsx
   do k=1,nz
    uout(i,k,j,1)=0.5d0*dble(ac(k,i,j))
    uout(i,k,j,2)=0.5d0*dble(bc(k,i,j))
   enddo
  enddo
enddo
!$acc end kernels
!$acc end data
deallocate(a,b,ac,bc)
#endif

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dctz_bwd_1d(rin,rout)

#define openaccflag openacccompflag
#if openaccflag == 0
use fftw3
#endif
#if openaccflag == 1
use cufft
use openacc
use cufftplans
#endif

implicit none
integer(c_int) :: nz !sx,nz,npy
real(c_double), dimension(:) :: rin, rout
!use only with cuFFT
#if openaccflag == 1
integer :: k
real(c_double), allocatable :: a(:)
complex(c_double_complex), allocatable :: ac(:)
#endif

! aliasing = 0 by default
! get dimensions
nz=size(rin)

#if openaccflag == 0
! single dct at a time, try with many dft
rin(nz)=2.0d0*rin(nz)
call fftw_execute_r2r(plan_z_bwd,rin,rout)
rout=rout*0.5d0
#endif

#if openaccflag == 1
allocate(a(2*(nz-1)))
allocate(ac(nz))
!$acc kernels
rin(nz)=2.0d0*rin(nz)
do k=1,nz-1
 a(k)=rin(k)
enddo
do k=nz,2*(nz-1)
 a(k)=rin(2*(nz-1)-k+2)
enddo
!$acc end kernels
!$acc host_data use_device(a,ac)
 gerr=gerr+cufftExecD2Z(cudaplan_z_bwd_1d,a,ac)
!$acc end host_data
!$acc kernels
rout=0.5d0*dble(ac)
!$acc end kernels
deallocate(a,ac)
#endif

end subroutine
end module
