module fftx_fwd_module
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftx_fwd(u,uc,aliasing)

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
integer :: aliasing,i,nx,npz,npy
real(c_double), dimension(:,:,:) :: u
real(c_double), dimension(:,:,:,:) :: uc
complex(c_double_complex), allocatable :: wt(:,:,:)

!! get dimensions
nx=size(u,1)
npz=size(u,2)
npy=size(u,3)
!! temp for storing complex fft results
allocate(wt(nx/2+1,npz,npy))

#if openaccflag == 0
call fftw_execute_dft_r2c(plan_x_fwd,u,wt)

uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))

! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
#endif


#if openaccflag == 1
!$acc data copyin(u)  copyout(wt)
!$acc host_data use_device(u,wt)
gerr=gerr+cufftExecD2Z(cudaplan_x_fwd,u,wt)
!$acc end host_data
!$acc end data
!$acc kernels
uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))
!$acc end kernels
!$acc kernels
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
!$acc end kernels
#endif

deallocate(wt)
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fftx_fwd_fg(u,uc,aliasing)

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
integer :: aliasing,i,nx,npz,npy 
real(c_double), dimension(:,:,:) :: u
real(c_double), dimension(:,:,:,:) :: uc
complex(c_double_complex),  allocatable :: wt(:,:,:)

!! get dimensions (npsix, npsiy, npsiz)
nx=size(u,1)  
npz=size(u,2) 
npy=size(u,3)
!! temp for storing complex fft results
allocate(wt(nx/2+1,npz,npy))

#if openaccflag == 0
call fftw_execute_dft_r2c(plan_x_fwd_fg,u,wt)

uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))

! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
#endif


#if openaccflag == 1
!$acc data copyin(u), copyout(wt)
!$acc host_data use_device(u,wt)
gerr=gerr+cufftExecD2Z(cudaplan_x_fwd_fg,u,wt)
!$acc end host_data
!$acc end data
!$acc kernels
uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))
!$acc end kernels
!$acc kernels
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif
!$acc end kernels
#endif

deallocate(wt)
end subroutine
end module

