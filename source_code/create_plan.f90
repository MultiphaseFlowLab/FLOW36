subroutine create_plan

use commondata
use par_size
use fftw3 !!until hybrid transfroms, then removed

#define openaccflag openacccompflag
#if openaccflag ==0
use fftw3
#endif
#if openaccflag == 1
use cufft
use cufftplans
#endif

integer :: nsx,npz,npy
integer :: rx,rz,ry
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
#define flag_fftw precisionflag


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft x
allocate(b_in(nx,fpz,fpy))
allocate(b_t(nx/2+1,fpz,fpy))

! fwd
inembed=[nx,fpz,fpy]
onembed=[nx/2+1,fpz,fpy]
istride=1
ostride=1
idist=nx
odist=nx/2+1
dims(1)=nx

#if openaccflag == 0
! create FFTW plan (CPU)
#if flag_fftw == 0
plan_x_fwd=fftw_plan_many_dft_r2c(1,dims,fpz*fpy,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_x_fwd=fftw_plan_many_dft_r2c(1,dims,fpz*fpy,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_x_fwd)
gerr=gerr+cufftPlanMany(cudaplan_x_fwd,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_D2Z,fpz*fpy)
if (gerr.ne.0) write(*,*) "cplanxfwd:", gerr !!to be removed
#endif


! bwd
inembed=[nx/2+1,fpz,fpy]
onembed=[nx,fpz,fpy]
istride=1
ostride=1
idist=nx/2+1
odist=nx

dims(1)=nx


#if openaccflag == 0
#if flag_fftw == 0
! create FFTW plan (CPU)
plan_x_bwd=fftw_plan_many_dft_c2r(1,dims,fpz*fpy,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw==1
plan_x_bwd=fftw_plan_many_dft_c2r(1,dims,fpz*fpy,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_x_bwd)
gerr=gerr+cufftPlanMany(cudaplan_x_bwd,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2D,fpz*fpy)
if (gerr.ne.0) write(*,*) "cplanxbwd:", gerr  !!to be removed
#endif

deallocate(b_in)
deallocate(b_t)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft y
rx=mod(nx/2+1,nycpu)
nsx=int((nx/2+1-rx)/nycpu)
if(mod(rank,nycpu).lt.rx)then
 nsx=int((nx/2+1-rx)/nycpu)+1
endif

rz=mod(nz,nzcpu)
npz=int((nz-rz)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.rz)then
 npz=int((nz-rz)/nzcpu)+1
endif

allocate(b_t(nsx,npz,ny))
allocate(b_out(nsx,npz,ny))

! fwd
inembed=[nsx,npz,ny]
onembed=[nsx,npz,ny]
istride=nsx*npz
ostride=nsx*npz
idist=1
odist=1

dims(1)=ny


#if openaccflag == 0
! create FFTW plan (CPU)
#if flag_fftw == 0
plan_y_fwd=fftw_plan_many_dft(1,dims,nsx*npz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_fwd=fftw_plan_many_dft(1,dims,nsx*npz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_PATIENT)
#endif
#endif


#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_y_fwd)
gerr=gerr+cufftPlanMany(cudaplan_y_fwd,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2Z,nsx*npz)
write(*,*) "cplanyfwd:", gerr !!to be removed
#endif


! bwd
inembed=[nsx,npz,ny]
onembed=[nsx,npz,ny]
istride=nsx*npz
ostride=nsx*npz
idist=1
odist=1

dims(1)=ny

#if openaccflag == 0
! +1=FFT_BACKWARD
#if flag_fftw == 0
plan_y_bwd=fftw_plan_many_dft(1,dims,nsx*npz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_bwd=fftw_plan_many_dft(1,dims,nsx*npz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_y_bwd)
gerr=gerr+cufftPlanMany(cudaplan_y_bwd,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2Z,nsx*npz)
write(*,*) "cplanybwd:", gerr !!to be removed
#endif



deallocate(b_t)
deallocate(b_out)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dct z
ry=mod(ny,nzcpu)
npy=int((ny-ry)/nzcpu)
if(ry.eq.0)then
 ngy=npy
else
 ngy=npy+1
endif

if(floor(real(rank)/real(nycpu)).lt.ry)then
 npy=int((ny-ry)/nzcpu)+1
endif
rz=mod(nz,nzcpu)


!fwd
#if openaccflag == 0
#if flag_fftw == 0
plan_z_fwd=fftw_plan_r2r_1d(nz,tin,tout, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_fwd=fftw_plan_r2r_1d(nz,tin,tout, FFTW_REDFT00, FFTW_PATIENT)
#endif
#endif


#if openaccflag == 1
print*,"fpy",fpy
print*,"npy",npy
print*,"nsx",nsx
print*,"nz",nz
! Done in a single block by transposing the matrix
! Othersiw, it cannot be done as an entire block.
! DCTs by rows and slices are much slower 
inembed=[2*(nz-1),nsx,fpy]
istride=1
idist=2*(nz-1)
onembed=[nz,nsx,fpy]
ostride=1
odist=nz
dims(1)=2*(nz-1)
! Create cuFFT plan (GPU)
gerr=cufftCreate(cudaplan_z_fwd)
gerr=gerr+cufftPlanMany(cudaplan_z_fwd,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_D2Z,nsx*npy)
write(*,*) "cplanzfwd:", gerr !!to be removed
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









subroutine create_plan_fg

use fftw3
use commondata
use par_size
use dual_grid !!until hybrid transfroms, then removed

#define openaccflag openacccompflag
#if openaccflag ==0
use fftw3
#endif
#if openaccflag == 1
use cufft
use cufftplans
#endif

integer :: nsx,npz
integer :: rx,rz
integer(c_int) :: dims(1)
integer(c_int) :: inembed(3),onembed(3),istride,ostride,idist,odist
real(c_double), allocatable :: b_in(:,:,:)
real(c_double) :: tin(npsiz),tout(npsiz)
complex(c_double_complex), allocatable :: b_out(:,:,:), b_t(:,:,:)

! notes on FFTW planner flags:
! EXHAUSTIVE : too slow, do not use it
! PATIENT: high plan creation time, slightly lower fft execution time than FFTW_ESTIMATE
! ESTIMATE: low plan creation time, choose "random" algorithm, may not be the best

! compiler option: 0 FFTW_ESTIMATE, 1 FFTW_PATIENT
#define flag_fftw precisionflag


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft x
allocate(b_in(npsix,fpzpsi,fpypsi))
allocate(b_t(npsix/2+1,fpzpsi,fpypsi))

! fwd
inembed=[npsix,fpzpsi,fpypsi]
onembed=[npsix/2+1,fpzpsi,fpypsi]
istride=1
ostride=1
idist=npsix
odist=npsix/2+1

dims(1)=npsix

#if openaccflag == 0
! create FFTW plan (CPU)
#if flag_fftw == 0
plan_x_fwd_fg=fftw_plan_many_dft_r2c(1,dims,fpzpsi*fpypsi,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_x_fwd_fg=fftw_plan_many_dft_r2c(1,dims,fpzpsi*fpypsi,b_in,inembed,istride,idist,&
 &    b_t,onembed,ostride,odist,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_x_fwd_fg)
gerr=gerr+cufftPlanMany(cudaplan_x_fwd_fg,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_D2Z,fpzpsi*fpypsi)
if (gerr.ne.0) write(*,*) "cplanxfwd_fg:", gerr !!to be removed
#endif


! bwd
inembed=[npsix/2+1,fpzpsi,fpypsi]
onembed=[npsix,fpzpsi,fpypsi]
istride=1
ostride=1
idist=npsix/2+1
odist=npsix

dims(1)=npsix

#if openaccflag == 0
! create FFTW plan (CPU) 
#if flag_fftw == 0
plan_x_bwd_fg=fftw_plan_many_dft_c2r(1,dims,fpzpsi*fpypsi,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_ESTIMATE)
#elif flag_fftw==1
plan_x_bwd_fg=fftw_plan_many_dft_c2r(1,dims,fpzpsi*fpypsi,b_t,inembed,istride,idist, &
 &    b_in,onembed,ostride,odist,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_x_bwd_fg)
gerr=gerr+cufftPlanMany(cudaplan_x_bwd_fg,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2D,fpzpsi*fpypsi)
if (gerr.ne.0) write(*,*) "cplanxbwd_fg:", gerr  !!to be removed
#endif


deallocate(b_in)
deallocate(b_t)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fft y
rx=mod(npsix/2+1,nycpu)
nsx=int((npsix/2+1-rx)/nycpu)
if(mod(rank,nycpu).lt.rx)then
 nsx=int((npsix/2+1-rx)/nycpu)+1
endif

rz=mod(npsiz,nzcpu)
npz=int((npsiz-rz)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.rz)then
 npz=int((npsiz-rz)/nzcpu)+1
endif

allocate(b_t(nsx,npz,npsiy))
allocate(b_out(nsx,npz,npsiy))

! fwd
inembed=[nsx,npz,npsiy]
onembed=[nsx,npz,npsiy]
istride=nsx*npz
ostride=nsx*npz
idist=1
odist=1

dims(1)=npsiy

#if openaccflag == 0
! create FFTW plan (CPU)
#if flag_fftw == 0
plan_y_fwd_fg=fftw_plan_many_dft(1,dims,nsx*npz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_fwd_fg=fftw_plan_many_dft(1,dims,nsx*npz,b_t,inembed,istride,idist, &
 &    b_out,onembed,ostride,odist,-1,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_y_fwd_fg)
gerr=gerr+cufftPlanMany(cudaplan_y_fwd_fg,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2Z,nsx*npz)
write(*,*) "cplanxbwd_fg:", gerr  !!to be removed
#endif

! bwd
inembed=[nsx,npz,npsiy]
onembed=[nsx,npz,npsiy]
istride=nsx*npz
ostride=nsx*npz
idist=1
odist=1

dims(1)=npsiy

#if openaccflag == 0
! create FFTW plan (CPU)
#if flag_fftw == 0
plan_y_bwd_fg=fftw_plan_many_dft(1,dims,nsx*npz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_y_bwd_fg=fftw_plan_many_dft(1,dims,nsx*npz,b_out,inembed,istride,idist, &
 &    b_t,onembed,ostride,odist,+1,FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Create cuFFT plan (GPU)
gerr=gerr+cufftCreate(cudaplan_y_bwd_fg)
gerr=gerr+cufftPlanMany(cudaplan_y_bwd_fg,1,dims,inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_Z2Z,nsx*npz)
write(*,*) "cplanybwd_fg:", gerr  !!to be removed
#endif

deallocate(b_t)
deallocate(b_out)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dct z

! fwd
#if openaccflag == 0
#if flag_fftw == 0
plan_z_fwd_fg=fftw_plan_r2r_1d(npsiz,tin,tout, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_fwd_fg=fftw_plan_r2r_1d(npsiz,tin,tout, FFTW_REDFT00, FFTW_PATIENT)
#endif
#endif

#if openaccflag == 1
! Faster with a large block of transforms
inembed=[2*(npsiz-1),nsx,fpypsi]
istride=1
idist=2*(npsiz-1)
onembed=[npsiz,nsx,fpypsi]
ostride=1
odist=npsiz
! Create cuFFT plan (GPU)
gerr=cufftCreate(cudaplan_z_fwd_fg)
gerr=gerr+cufftPlanMany(cudaplan_z_fwd_fg,1,2*(npsiz-1),inembed,istride,idist,onembed,&
 &    ostride,odist,CUFFT_D2Z,nsx*npy)
write(*,*) "cplanzfwd:", gerr !!to be removed
#endif

! bwd
#if flag_fftw == 0
plan_z_bwd_fg=fftw_plan_r2r_1d(npsiz,tout,tin, FFTW_REDFT00, FFTW_ESTIMATE)
#elif flag_fftw == 1
plan_z_bwd_fg=fftw_plan_r2r_1d(npsiz,tout,tin, FFTW_REDFT00, FFTW_PATIENT)
#endif


return
end
