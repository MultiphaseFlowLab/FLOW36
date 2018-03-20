subroutine fftx_bwd(uc,ur,nx,npz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nx,npy,npz
integer :: aliasing

real(c_double) :: ur(nx,npz,npy),uc(nx/2+1,npz,npy,2)
complex(c_double_complex) :: wt(nx/2+1,npz,npy)


! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))


call fftw_execute_dft_c2r(plan_x_bwd,wt,ur)


! code is slower with nested loops in fftw, faster in dct
!plan=fftw_plan_dft_c2r(1,dims,tc,tr,fftw_estimate)
!do i=1,npz
! do j=1,npy
!  tc(:)=wt(:,i,j)
!  call fftw_execute_dft_c2r(plan,tc,tr)
!  ur(:,i,j)=tr(:)
! enddo
!enddo

ur=ur/dble(nx)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fftx_bwd_fg(uc,ur,nx,npz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nx,npy,npz
integer :: aliasing

real(c_double) :: ur(nx,npz,npy),uc(nx/2+1,npz,npy,2)
complex(c_double_complex) :: wt(nx/2+1,npz,npy)


! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

wt(1:nx/2+1,1:npz,1:npy)=dcmplx(uc(1:nx/2+1,1:npz,1:npy,1),uc(1:nx/2+1,1:npz,1:npy,2))


call fftw_execute_dft_c2r(plan_x_bwd_fg,wt,ur)

ur=ur/dble(nx)

return
end
