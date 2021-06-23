subroutine fftx_fwd(u,uc,nx,npz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nx,npy,npz
integer :: aliasing

real(c_double) :: u(nx,npz,npy),uc(nx/2+1,npz,npy,2)
complex(c_double_complex) :: wt(nx/2+1,npz,npy)


!write(*,*) 'verificare funzionamento fftw con funzioni trigonometriche note'

call fftw_execute_dft_r2c(plan_x_fwd,u,wt)



! just for debug
!write(*,'(6(ES16.4E6,1x))') wt(1,1,1),wt(2,1,1),wt(3,1,1)
!write(*,'(6(ES16.4E6,1x))') wt(1,2,1),wt(2,2,1),wt(3,2,1)


uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))


! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fftx_fwd_fg(u,uc,nx,npz,npy,aliasing)

use fftw3
implicit none

integer(c_int) :: nx,npy,npz
integer :: aliasing

real(c_double) :: u(nx,npz,npy),uc(nx/2+1,npz,npy,2)
complex(c_double_complex) :: wt(nx/2+1,npz,npy)

call fftw_execute_dft_r2c(plan_x_fwd_fg,u,wt)


uc(1:nx/2+1,1:npz,1:npy,1)=dble(wt(1:nx/2+1,1:npz,1:npy))
uc(1:nx/2+1,1:npz,1:npy,2)=aimag(wt(1:nx/2+1,1:npz,1:npy))


! dealiasing
if(aliasing.eq.1)then
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,1)=0.0d0
 uc(floor(2.0/3.0*real(nx/2+1))+1:nx/2+1,1:npz,1:npy,2)=0.0d0
endif

return
end
