module ffty_fwd_module
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ffty_fwd(ui,uo,aliasing)

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
integer :: aliasing,nsx,ny,npz
real(c_double), dimension(:,:,:,:) :: ui, uo
complex(c_double_complex), allocatable :: wt(:,:,:),wot(:,:,:)

!! get dimensions
nsx=size(ui,1)
npz=size(ui,2)
ny=size(ui,3)
!! temp for storing complex fft results
allocate(wt(nsx,npz,ny))
allocate(wot(nsx,npz,ny))

#if openaccflag == 0
wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

call fftw_execute_dft(plan_y_fwd,wt,wot)

uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

! dealiasing
if(aliasing.eq.1)then
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
#endif


#if openaccflag == 1
!$acc kernels
wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))
!$acc end kernels
!$acc data copyin(wt), copyout(wot)
!$acc host_data use_device(wt,wot)
gerr=gerr+cufftExecZ2Z(cudaplan_y_fwd,wt,wot,CUFFT_FORWARD)
!$acc end host_data
!$acc end data
!$acc kernels
uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))
!$acc end kernels
!$acc kernels
if(aliasing.eq.1)then
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
!$acc end kernels
#endif

deallocate(wt,wot)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ffty_fwd_fg(ui,uo,aliasing)

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
integer(c_int) :: nsx,ny,npz
integer :: aliasing
real(c_double), dimension(:,:,:,:) :: ui,uo
complex(c_double_complex), allocatable :: wt(:,:,:),wot(:,:,:)

! get dimensions
nsx=size(ui,1)
npz=size(ui,2)
ny=size(ui,3)
! temp for storing complex fft results
allocate(wt(nsx,npz,ny))
allocate(wot(nsx,npz,ny))

#if openaccflag == 0

wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

call fftw_execute_dft(plan_y_fwd_fg,wt,wot)

uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

! dealiasing
if(aliasing.eq.1)then
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
#endif


#if openaccflag == 1
!$acc kernels
wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))
!$acc end kernels
!$acc data copyin(wt), copyout(wot)
!$acc host_data use_device(wt,wot)
gerr=gerr+cufftExecZ2Z(cudaplan_y_fwd_fg,wt,wot,CUFFT_FORWARD)
!$acc end host_data
!$acc end data
!$acc kernels
uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))
!$acc end kernels
!$acc kernels
if(aliasing.eq.1)then
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 uo(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1)):ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
!$acc end kernels
#endif

deallocate(wt,wot)
end subroutine
end module



