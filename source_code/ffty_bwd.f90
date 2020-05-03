subroutine ffty_bwd(ui,uo,nsx,npz,ny,aliasing,inloop)

use fftw3
#define GPU_RUN gpucompflag
#if GPU_RUN == 1
use interfaccia
#endif
implicit none

!type(c_ptr) :: plan
integer(c_int) :: nsx,ny,npz
!integer(c_int) :: dims(1)
!integer(c_int) :: inembed(3),onembed(3),istride,ostride,idist,odist
integer :: aliasing, inloop

real(c_double) :: ui(nsx,npz,ny,2),uo(nsx,npz,ny,2)
complex(c_double_complex) :: wt(nsx,npz,ny),wot(nsx,npz,ny)

! dealiasing
if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif


wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

!inembed=[nsx,npz,ny]
!onembed=[nsx,npz,ny]
!istride=nsx*npz
!ostride=nsx*npz
!idist=1
!odist=1

!dims(1)=ny


!! +1=FFT_BACKWARD
!plan=fftw_plan_many_dft(1,dims,nsx*npz,wt,inembed,istride,idist, &
! &    wot,onembed,ostride,odist,+1,FFTW_ESTIMATE)

#if GPU_RUN == 1
  if (inloop == 1) then
    call h_ffty_many_bwd(ui(:,:,:,1),ui(:,:,:,2),uo(:,:,:,1),uo(:,:,:,2),aliasing)
    uo=uo/dble(ny)
  else
	call fftw_execute_dft(plan_y_bwd,wt,wot)
    uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
    uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))
    uo=uo/dble(ny)
  endif
#else
  
call fftw_execute_dft(plan_y_bwd,wt,wot)

!call fftw_destroy_plan(plan)

uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

uo=uo/dble(ny)

#endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ffty_bwd_fg(ui,uo,nsx,npz,ny,aliasing,inloop)

use fftw3
#define GPU_RUN gpucompflag
#if GPU_RUN == 1
use interfaccia
#endif
implicit none

integer(c_int) :: nsx,ny,npz
integer :: aliasing, inloop

real(c_double) :: ui(nsx,npz,ny,2),uo(nsx,npz,ny,2)
complex(c_double_complex) :: wt(nsx,npz,ny),wot(nsx,npz,ny)

! dealiasing
if(aliasing.eq.1)then
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),1)=0.0d0
 ui(1:nsx,1:npz,floor(2.0/3.0*real(ny/2+1))+1:ny-floor(2.0/3.0*real(ny/2)),2)=0.0d0
endif

wt(1:nsx,1:npz,1:ny)=dcmplx(ui(1:nsx,1:npz,1:ny,1),ui(1:nsx,1:npz,1:ny,2))

#if GPU_RUN == 1
  if (inloop == 1) then
    call h_fftymanybwd_fg(ui(:,:,:,1),ui(:,:,:,2),uo(:,:,:,1),uo(:,:,:,2),aliasing)
    uo=uo/dble(ny)
  else
	call fftw_execute_dft(plan_y_bwd_fg,wt,wot)
    uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
    uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

    uo=uo/dble(ny)
  endif
#else
  
  call fftw_execute_dft(plan_y_bwd_fg,wt,wot)

  uo(1:nsx,1:npz,1:ny,1)=dble(wot(1:nsx,1:npz,1:ny))
  uo(1:nsx,1:npz,1:ny,2)=aimag(wot(1:nsx,1:npz,1:ny))

  uo=uo/dble(ny)

#endif


return
end
