subroutine sterm_temp(stheta)

use commondata
use par_size
use velocity
use temperature
use wavenumber

double precision :: stheta(spx,nz,spy,2)
double precision, allocatable, dimension(:,:,:,:) :: a1
double precision, allocatable, dimension(:,:,:) :: a1f

integer :: i,j,k


stheta=0.0d0

! calculate convective term
allocate(a1(spx,nz,spy,2))
allocate(a1f(nx,fpz,fpy))

! u*(d theta /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*thetac(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*thetac(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1f(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo


call phys_to_spectral(a1f,a1,1)

! sum -u*(d theta /dx) to stheta
stheta=-a1

! v*(d theta /dy)

do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*thetac(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*thetac(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1f(i,k,j)=a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(a1f,a1,1)

! sum -v*(d theta /dy) to stheta
stheta=stheta-a1

! w*(d theta /dz)

call dz(thetac,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1f(i,k,j)=a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(a1f,a1,1)

! sum -w*(d theta /dz) to stheta
stheta=stheta-a1

deallocate(a1)
deallocate(a1f)


return
end
