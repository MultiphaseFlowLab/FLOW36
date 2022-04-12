module dctz_fwd_module
implicit none
contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dctz_fwd(uin,uout,aliasing)

#define openaccflag openacccompflag
#if openaccflag == 0
use fftw3
#endif
#if openaccflag == 1
use cufft
use openacc
use cufftplans
#endif

integer(c_int) :: nsx,nz,npy
real(c_double), dimension(:,:,:,:) :: uin, uout
real(c_double), allocatable ::  tin(:),tout(:)
integer :: aliasing,i,k,j
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
! single dct at a time, try with many dft
allocate(tin(nz),tout(nz))
do i=1,nsx
 do j=1,npy
  do k=1,2
   tin(:)=uin(i,:,j,k)
   call fftw_execute_r2r(plan_z_bwd,tin,tout)
   uout(i,:,j,k)=tout(:)/dble(nz-1)
  enddo
 enddo
enddo

uout(:,nz,:,:)=0.5d0*uout(:,nz,:,:)

! dealiasing
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
deallocate(tin,tout)
#endif



#if openaccflag == 1
! trick to make DCT faster and in one block (not possible otherwise)
! transpose uin and uout and perform DCT along 1st direction 
allocate(a(2*(nz-1),nsx,npy))
allocate(b(2*(nz-1),nsx,npy))
allocate(ac(nz,nsx,npy))
allocate(bc(nz,nsx,npy))
!trasnpose the z-row and make them even symmetric (no R2R in cuFFT)
!$acc data copyin(uin) create(a,b,ac,bc) copyout(uout)
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
 gerr=gerr+cufftExecD2Z(cudaplan_z_fwd,a,ac)
!$acc end host_data
!$acc host_data use_device(b,bc) 
gerr=gerr+cufftExecD2Z(cudaplan_z_fwd,b,bc)
!$acc end host_data
!$acc kernels
do j=1,npy
  do i=1,nsx
   do k=1,nz
    uout(i,k,j,1)=dble(ac(k,i,j))/dble(nz-1)
    uout(i,k,j,2)=dble(bc(k,i,j))/dble(nz-1)
   enddo
  enddo
enddo
!$acc end kernels
!$acc end data
!$acc kernels
uout(:,nz,:,1)=0.5d0*uout(:,nz,:,1)
uout(:,nz,:,2)=0.5d0*uout(:,nz,:,2)
!$acc end kernels
!dealiasing
!$acc kernels
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
!$acc end kernels
deallocate(a,b,ac,bc)
#endif

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dctz_fwd_fg(uin,uout,aliasing)

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
real(c_double), dimension(:,:,:,:) :: uin, uout
real(c_double), allocatable ::  tin(:),tout(:)
integer :: aliasing,i,k,j
!used only with cuFFT
#if openaccflag == 1
double precision, allocatable :: a(:,:,:), b(:,:,:)
complex(c_double_complex), allocatable :: ac(:,:,:), bc(:,:,:)
#endif

! get dimensions
nsx=size(uin,1)
nz=size(uin,2)
npy=size(uin,3)

#if openaccflag == 0
allocate(tin(nz),tout(nz))
! single dct at a time, try with many dft
do i=1,nsx
 do j=1,npy
  do k=1,2
   tin(:)=uin(i,:,j,k)
   call fftw_execute_r2r(plan_z_bwd_fg,tin,tout)
   uout(i,:,j,k)=tout(:)/dble(nz-1)
  enddo
 enddo
enddo
uout(:,nz,:,:)=0.5d0*uout(:,nz,:,:)
! dealiasing
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
deallocate(tin,tout)
#endif


#if openaccflag == 1
allocate(a(2*(nz-1),nsx,npy))
allocate(b(2*(nz-1),nsx,npy))
allocate(ac(nz,nsx,npy))
allocate(bc(nz,nsx,npy))

!trasnpose the z-row and make them even symmetric (no R2R in cuFFT)
!$acc data copyin(uin) create(a,b,ac,bc) copyout(uout)
!$acc parallel loop collapse(2) async
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
!$acc wait
!$acc host_data use_device(a,ac)
 gerr=gerr+cufftExecD2Z(cudaplan_z_fwd_fg,a,ac)
!$acc end host_data
!$acc host_data use_device(b,bc) 
gerr=gerr+cufftExecD2Z(cudaplan_z_fwd_fg,b,bc)
!$acc end host_data
!$acc kernels
 do j=1,npy
  do i=1,nsx
   do k=1,nz
    uout(i,k,j,1)=dble(ac(k,i,j))/dble(nz-1)
    uout(i,k,j,2)=dble(bc(k,i,j))/dble(nz-1)
   enddo
  enddo
enddo
!$acc end kernels
!$acc end data
!$acc kernels
uout(:,nz,:,1)=0.5d0*uout(:,nz,:,1)
uout(:,nz,:,2)=0.5d0*uout(:,nz,:,2)
!$acc end kernels
!dealiasing
!$acc kernels
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
!$acc end kernels
deallocate(a,b,ac,bc)
#endif

end subroutine
end module
