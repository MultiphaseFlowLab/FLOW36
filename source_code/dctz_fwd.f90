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
real(c_double) :: tin(nz), tout(nz)
integer :: aliasing,i,k,j
real(c_double), allocatable :: a(:,:,:), b(:,:,:)
complex(c_double_complex), allocatable :: ac(:,:,:), bc(:,:,:)

! get dimensions
nsx=size(uin,1)
nz=size(uin,2)
npy=size(uin,3)

#if openaccflag == 0
! single dct at a time, try with many dft
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
#endif




#if openaccflag == 1
! trick to make DCT faster and in one block (not possible otherwise)
! transpose uin and uout and perform DCT along 1st direction 
allocate(a(2*(nz-1),nsx,npy))
allocate(b(2*(nz-1),nsx,npy))
allocate(ac(nz,nsx,npy))
allocate(bc(nz,nsx,npy))
!trasnpose the z-row and make them even symmetric (no R2R)

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
real(c_double) :: tin(nz),tout(nz)
integer :: aliasing,i,k,j
double precision, allocatable :: a(:,:,:), b(:,:,:)
complex(c_double_complex), allocatable :: ac(:,:,:), bc(:,:,:)
real(c_double), allocatable :: as(:,:), bs(:,:)
complex(c_double_complex), allocatable :: acs(:,:), bcs(:,:)

! get dimensions
nsx=size(uin,1)
nz=size(uin,2)
npy=size(uin,3)

#if openaccflag == 0
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
#endif


#if openaccflag == 1
!input and output for the two FFTs (use FFT to compute DCT)
allocate(a(nsx,2*(nz-1),npy))
allocate(b(nsx,2*(nz-1),npy))
allocate(ac(nsx,nz,npy))
allocate(bc(nsx,nz,npy))
!make the vector even + simmetric
!$acc parallel loop collapse(2)
do i=1,nsx
 do j=1,npy
   do k=1,nz-1
    a(i,k,j)=uin(i,k,j,1)
    b(i,k,j)=uin(i,k,j,2)
   enddo
   do k=nz,2*(nz-1)
    a(i,k,j)=uin(i,2*(nz-1)-k+2,j,1)
    b(i,k,j)=uin(i,2*(nz-1)-k+2,j,2)
   enddo
 enddo
enddo
!$acc end parallel

allocate(as(nsx,2*(nz-1)))
allocate(bs(nsx,2*(nz-1)))
allocate(acs(nsx,nz))
allocate(bcs(nsx,nz))

do j=1,npy
 as=a(:,:,j)
 !$acc data copyin(as)  copyout(acs) 
 !$acc host_data use_device(as,acs) 
  gerr=gerr+cufftExecD2Z(cudaplan_z_fwd,as,acs)
 !$acc end host_data
 !$acc end data
 bs=b(:,:,j)
 !$acc data copyin(bs)  copyout(bcs) 
 !$acc host_data use_device(bs,bcs) 
  gerr=gerr+cufftExecD2Z(cudaplan_z_fwd,bs,bcs)
 !$acc end host_data
 !$acc end data
 !$acc kernels
 ac(:,:,j)=acs
 bc(:,:,j)=bcs
 !$acc end kernels
enddo

!$acc kernels
uout(1:nsx,1:nz,1:npy,1)=dble(ac(1:nsx,1:nz,1:npy))
uout(1:nsx,1:nz,1:npy,2)=dble(bc(1:nsx,1:nz,1:npy))
uout=uout/dble(nz-1)
uout(:,nz,:,1)=0.5d0*uout(:,nz,:,1)
uout(:,nz,:,2)=0.5d0*uout(:,nz,:,2)
!dealiasing
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif
!$acc end kernels
#endif

deallocate(a,b,ac,bc,as,bs,acs,bcs)

end subroutine
end module
