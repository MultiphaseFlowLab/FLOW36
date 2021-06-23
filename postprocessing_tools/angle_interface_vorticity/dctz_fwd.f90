subroutine dctz_fwd(uin,uout,nsx,nz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nsx,nz,npy

real(c_double) :: uin(nsx,nz,npy,2),uout(nsx,nz,npy,2)

real(c_double) :: tin(nz),tout(nz)
integer :: i,j,k
integer :: aliasing


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

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dctz_fwd_red(uin,uout,nsx,nz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nsx,nz,npy

real(c_double) :: uin(nsx,nz,npy),uout(nsx,nz,npy)

real(c_double) :: tin(nz),tout(nz)
integer :: i,j
integer :: aliasing


! single dct at a time, try with many dft

do i=1,nsx
 do j=1,npy
   tin(:)=uin(i,:,j)
   call fftw_execute_r2r(plan_z_bwd,tin,tout)
   uout(i,:,j)=tout(:)/dble(nz-1)
 enddo
enddo

uout(:,nz,:)=0.5d0*uout(:,nz,:)

! dealiasing
if(aliasing.eq.1)then
 uout(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy)=0.0d0
endif

return
end
