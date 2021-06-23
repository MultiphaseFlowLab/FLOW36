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
    uin(1:nsx,ceiling(2.0*real(nz)/3.0):nz,1:npy,1)=0.0d0
    uin(1:nsx,ceiling(2.0*real(nz)/3.0):nz,1:npy,2)=0.0d0
  endif


  ! single dct at a time, try with many dft

  do i=1,nsx
    do j=1,npy
      do k=1,2
        tin(:)=uin(i,:,j,k)
        tin(nz)=2.0d0*tin(nz)
        call fftw_execute_r2r(plan_z_bwd,tin,tout)
        uout(i,:,j,k)=tout(:)*0.50d0
      enddo
    enddo
  enddo

  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dctz_bwd_fg(uin,uout,nsx,nz,npy,aliasing)

  use fftw3
  implicit none

  integer(c_int) :: nsx,nz,npy
  integer :: aliasing

  real(c_double) :: uin(nsx,nz,npy,2),uout(nsx,nz,npy,2)

  real(c_double) :: tin(nz),tout(nz)
  integer :: i,j,k


  ! dealiasing
  if(aliasing.eq.1)then
    uin(1:nsx,ceiling(2.0*real(nz)/3.0):nz,1:npy,1)=0.0d0
    uin(1:nsx,ceiling(2.0*real(nz)/3.0):nz,1:npy,2)=0.0d0
  endif


  ! single dct at a time, try with many dft

  do i=1,nsx
    do j=1,npy
      do k=1,2
        tin(:)=uin(i,:,j,k)
        tin(nz)=2.0d0*tin(nz)
        call fftw_execute_r2r(plan_z_bwd_fg,tin,tout)
        uout(i,:,j,k)=tout(:)*0.50d0
      enddo
    enddo
  enddo

  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dctz_bwd_1d(uin,uout,nz,aliasing)

  use fftw3
  implicit none

  integer(c_int) :: nz
  integer :: aliasing
  real(c_double) :: uin(nz),uout(nz)
  real(c_double) :: tin(nz),tout(nz)


  ! dealiasing
  if(aliasing.eq.1)then
    uin(ceiling(2.0*real(nz)/3.0):nz)=0.0d0
  endif

  tin(:)=uin(:)
  tin(nz)=2.0d0*tin(nz)
  call fftw_execute_r2r(plan_z_bwd,tin,tout)
  uout(:)=tout(:)*0.50d0

  return
end subroutine
