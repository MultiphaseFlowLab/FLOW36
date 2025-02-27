module ffty_bwd_module
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ffty_bwd(ui,uo,aliasing)

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
integer :: aliasing, nsx,ny,npz
real(c_double), dimension(:,:,:,:) :: ui,uo
complex(c_double_complex),allocatable :: wt(:,:,:),wot(:,:,:)

! get dimensions
nsx=size(ui,1)
npz=size(ui,2)
ny=size(ui,3)
! temp for storing complex fft results
allocate(wt(nsx,npz,ny))
allocate(wot(nsx,npz,ny))

#if openaccflag == 0

if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif

wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

call fftw_execute_dft(plan_y_bwd,wt,wot)

uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

uo=uo/dble(ny)
#endif


#if openaccflag == 1
!$acc data copyin(ui) create(wt,wot) copyout(uo)
!$acc kernels
if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))
!$acc end kernels
!$acc host_data use_device(wt,wot)
gerr=gerr+cufftExecZ2Z(cudaplan_y_bwd,wt,wot,CUFFT_INVERSE)
!$acc end host_data
!$acc kernels
uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))
uo=uo/dble(ny)
!$acc end kernels
!$acc end data
#endif

deallocate(wt,wot)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ffty_bwd_fg(ui,uo,aliasing)

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
integer(c_int) :: aliasing,nsx,ny,npz
real(c_double), dimension(:,:,:,:) :: ui, uo
complex(c_double_complex), allocatable :: wt(:,:,:),wot(:,:,:)

! get dimensions
nsx=size(ui,1)
npz=size(ui,2)
ny=size(ui,3)
! temp for storing complex fft results
allocate(wt(nsx,npz,ny))
allocate(wot(nsx,npz,ny))

#if openaccflag == 0
! dealiasing
if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif

wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

call fftw_execute_dft(plan_y_bwd_fg,wt,wot)

uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

uo=uo/dble(ny)
#endif


#if openaccflag == 1
!$acc data copyin(ui) create(wt,wot) copyout(uo)
!$acc kernels
if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif
wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))
!$acc end kernels
!$acc host_data use_device(wt,wot)
gerr=gerr+cufftExecZ2Z(cudaplan_y_bwd_fg,wt,wot,CUFFT_INVERSE)
!$acc end host_data
!$acc kernels
uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))
uo=uo/dble(ny)
!$acc end kernels
!$acc end data
#endif

deallocate(wt,wot)
end subroutine
end module
