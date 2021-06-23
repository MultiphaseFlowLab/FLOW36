subroutine dctz_bwd(uin,uout,nsx,nz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nsx,nz,npy
integer :: aliasing

real(c_double) :: uin(nsx,nz,npy,2),uout(nsx,nz,npy,2)

real(c_double) :: tin(nz),tout(nz)
integer :: i,j,k


! dealiasing
if(aliasing.eq.1)then
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,1)=0.0d0
 uin(1:nsx,floor(2.0*real(nz)/3.0)+1:nz,1:npy,2)=0.0d0
endif


! single dct at a time, try with many dft

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





return
end
