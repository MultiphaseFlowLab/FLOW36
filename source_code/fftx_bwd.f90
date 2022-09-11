module fftx_bwd_module
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftx_bwd(uc,ur,aliasing)

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
integer :: aliasing,nx,npy,npz
real(c_double), dimension(:,:,:) :: ur
real(c_double), dimension(:,:,:,:) :: uc
complex(c_double_complex), allocatable :: wt(:,:,:)

!! get dimensions
nx=size(ur,1)
npz=size(ur,2)
npy=size(ur,3)
!! temp for storing complex fft results
allocate(wt(nx/2+1,npz,npy))

#if openaccflag == 0
! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))

call fftw_execute_dft_c2r(plan_x_bwd,wt,ur)

ur=ur/dble(nx)
#endif


#if openaccflag == 1
!$acc kernels
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))
!$acc end kernels
!$acc data copyin(wt) copyout(ur)
!$acc host_data use_device(wt,ur)
gerr=gerr+cufftExecZ2D(cudaplan_x_bwd,wt,ur)
!$acc end host_data
!$acc end data
!$acc kernels
ur=ur/dble(nx)
!$acc end kernels
#endif

deallocate(wt)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftx_bwd_fg(uc,ur,aliasing)

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
integer :: aliasing,nx,npy,npz
real(c_double), dimension(:,:,:) :: ur
real(c_double), dimension(:,:,:,:) :: uc
complex(c_double_complex), allocatable :: wt(:,:,:)

!! get dimensions
nx=size(ur,1)
npz=size(ur,2)
npy=size(ur,3)
!! temp for storing complex fft results
allocate(wt(nx/2+1,npz,npy))

#if openaccflag == 0
! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))

call fftw_execute_dft_c2r(plan_x_bwd_fg,wt,ur)

ur=ur/dble(nx)
#endif

#if openaccflag == 1
!$acc kernels
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))
!$acc end kernels
!$acc data copyin(wt) copyout(ur)
!$acc host_data use_device(wt,ur)
gerr=gerr+cufftExecZ2D(cudaplan_x_bwd_fg,wt,ur)
!$acc end host_data
!$acc end data
!$acc kernels
ur=ur/dble(nx)
!$acc end kernels
#endif

deallocate(wt)
end subroutine
end module
