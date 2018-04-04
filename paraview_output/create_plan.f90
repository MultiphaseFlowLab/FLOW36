subroutine create_plan

use fftw3
use commondata

integer :: nsx

integer(c_int) :: dims(1)
integer(c_int) :: inembed(3),onembed(3),istride,ostride,idist,odist

real(c_double), allocatable :: b_in(:,:,:)
real(c_double) :: tin(nz),tout(nz)
complex(c_double_complex), allocatable :: b_out(:,:,:), b_t(:,:,:)

! notes on FFTW planner flags:
! EXHAUSTIVE : too slow, do not use it
! PATIENT: high plan creation time, slightly lower fft execution time than FFTW_ESTIMATE
! ESTIMATE: low plan creation time, choose "random" algorithm, may not be the best

! compiler option: 0 FFTW_ESTIMATE, 1 FFTW_PATIENT
write(*,*) 'Creation of transform plans'
write(*,*) 'This may take a while'
#define flag_fftw 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft x
allocate(b_in(nx,nz,ny))
allocate(b_t(nx/2+1,nz,ny))

! fwd
inembed=[nx,nz,ny]
onembed=[nx/2+1,nz,ny]
istride=1
ostride=1
idist=nx
odist=nx/2+1

dims(1)=nx

#if flag_fftw == 0
plan_x_fwd=fftw_plan_many_dft_r2c(1,dims,nz*ny,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_x_fwd=fftw_plan_many_dft_r2c(1,dims,nz*ny,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_PATIENT)
#endif


! bwd
inembed=[nx/2+1,nz,ny]
onembed=[nx,nz,ny]
istride=1
ostride=1
idist=nx/2+1
odist=nx

dims(1)=nx


#if flag_fftw == 0
plan_x_bwd=fftw_plan_many_dft_c2r(1,dims,nz*ny,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw==1
plan_x_bwd=fftw_plan_many_dft_c2r(1,dims,nz*ny,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_PATIENT)
#endif


deallocate(b_in)
deallocate(b_t)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft y
nsx=nx/2+1
allocate(b_t(nsx,nz,ny))
allocate(b_out(nsx,nz,ny))

! fwd
inembed=[nsx,nz,ny]
onembed=[nsx,nz,ny]
istride=nsx*nz
ostride=nsx*nz
idist=1
odist=1

dims(1)=ny


! -1=FFT_FORWARD
#if flag_fftw == 0
plan_y_fwd=fftw_plan_many_dft(1,dims,nsx*nz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_fwd=fftw_plan_many_dft(1,dims,nsx*nz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_PATIENT)
#endif

! bwd
inembed=[nsx,nz,ny]
onembed=[nsx,nz,ny]
istride=nsx*nz
ostride=nsx*nz
idist=1
odist=1

dims(1)=ny


! +1=FFT_BACKWARD
#if flag_fftw == 0
plan_y_bwd=fftw_plan_many_dft(1,dims,nsx*nz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_bwd=fftw_plan_many_dft(1,dims,nsx*nz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_PATIENT)
#endif

deallocate(b_t)
deallocate(b_out)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dct z

! fwd
#if flag_fftw == 0
plan_z_fwd=fftw_plan_r2r_1d(nz,tin,tout, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_fwd=fftw_plan_r2r_1d(nz,tin,tout, FFTW_REDFT00, FFTW_PATIENT)
#endif


! bwd
#if flag_fftw == 0
plan_z_bwd=fftw_plan_r2r_1d(nz,tout,tin, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_bwd=fftw_plan_r2r_1d(nz,tout,tin, FFTW_REDFT00, FFTW_PATIENT)
#endif



return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_plan_fg

use fftw3
use commondata

integer :: nsx,nfx,nfy,nfz

integer(c_int) :: dims(1)
integer(c_int) :: inembed(3),onembed(3),istride,ostride,idist,odist

real(c_double), allocatable :: b_in(:,:,:)
real(c_double) :: tin(exp_z*(nz-1)+1),tout(exp_z*(nz-1)+1)
complex(c_double_complex), allocatable :: b_out(:,:,:), b_t(:,:,:)

! notes on FFTW planner flags:
! EXHAUSTIVE : too slow, do not use it
! PATIENT: high plan creation time, slightly lower fft execution time than FFTW_ESTIMATE
! ESTIMATE: low plan creation time, choose "random" algorithm, may not be the best

! compiler option: 0 FFTW_ESTIMATE, 1 FFTW_PATIENT
#define flag_fftw 0

nfx=nx*exp_x
nfy=ny*exp_y
nfz=(nz-1)*exp_z+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft x
allocate(b_in(nfx,nfz,nfy))
allocate(b_t(nfx/2+1,nfz,nfy))

! fwd
inembed=[nfx,nfz,nfy]
onembed=[nfx/2+1,nfz,nfy]
istride=1
ostride=1
idist=nfx
odist=nfx/2+1

dims(1)=nfx

#if flag_fftw == 0
plan_x_fwd_fg=fftw_plan_many_dft_r2c(1,dims,nfz*nfy,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_x_fwd_fg=fftw_plan_many_dft_r2c(1,dims,nfz*nfy,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_PATIENT)
#endif


! bwd
inembed=[nfx/2+1,nfz,nfy]
onembed=[nfx,nfz,nfy]
istride=1
ostride=1
idist=nfx/2+1
odist=nfx

dims(1)=nfx


#if flag_fftw == 0
plan_x_bwd_fg=fftw_plan_many_dft_c2r(1,dims,nfz*nfy,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw==1
plan_x_bwd_fg=fftw_plan_many_dft_c2r(1,dims,nfz*nfy,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_PATIENT)
#endif


deallocate(b_in)
deallocate(b_t)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft y
nsx=nfx/2+1
allocate(b_t(nsx,nfz,nfy))
allocate(b_out(nsx,nfz,nfy))

! fwd
inembed=[nsx,nfz,nfy]
onembed=[nsx,nfz,nfy]
istride=nsx*nfz
ostride=nsx*nfz
idist=1
odist=1

dims(1)=nfy


! -1=FFT_FORWARD
#if flag_fftw == 0
plan_y_fwd_fg=fftw_plan_many_dft(1,dims,nsx*nfz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_fwd_fg=fftw_plan_many_dft(1,dims,nsx*nfz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_PATIENT)
#endif

! bwd
inembed=[nsx,nfz,nfy]
onembed=[nsx,nfz,nfy]
istride=nsx*nfz
ostride=nsx*nfz
idist=1
odist=1

dims(1)=nfy


! +1=FFT_BACKWARD
#if flag_fftw == 0
plan_y_bwd_fg=fftw_plan_many_dft(1,dims,nsx*nfz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_bwd_fg=fftw_plan_many_dft(1,dims,nsx*nfz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_PATIENT)
#endif

deallocate(b_t)
deallocate(b_out)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dct z

! fwd
#if flag_fftw == 0
plan_z_fwd_fg=fftw_plan_r2r_1d(nfz,tin,tout, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_fwd_fg=fftw_plan_r2r_1d(nfz,tin,tout, FFTW_REDFT00, FFTW_PATIENT)
#endif


! bwd
#if flag_fftw == 0
plan_z_bwd_fg=fftw_plan_r2r_1d(nfz,tout,tin, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_bwd_fg=fftw_plan_r2r_1d(nfz,tout,tin, FFTW_REDFT00, FFTW_PATIENT)
#endif



return
end
