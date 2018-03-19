subroutine sterm_ch(sphi)

use commondata
use par_size
use velocity
use phase_field
use wavenumber

double precision :: sphi(spx,nz,spy,2)
double precision, allocatable, dimension(:,:,:,:) :: a1,a2,a3
double precision, allocatable, dimension(:,:,:) :: a1f

integer :: i,j,k


allocate(a1(spx,nz,spy,2))

! calculate -(1+s)/Pe*nabla^2 phi
call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

sphi=-(1.0d0+s_coeff)/pe*sphi

! calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
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

! sum -u*(d phi /dx) to sphi
sphi=sphi-a1

! v*(d phi /dy)

do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
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

! sum -v*(d phi /dy) to sphi
sphi=sphi-a1

! w*(d phi /dz)

call dz(phic,a1)

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

! sum -w*(d phi /dz) to sphi
sphi=sphi-a1


! calculate 1/Pe nabla^2 phi^3
call spectral_to_phys(phic,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1f(i,k,j)=a1f(i,k,j)**3
    enddo
  enddo
enddo

call phys_to_spectral(a1f,a1,1)

deallocate(a1f)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo


deallocate(a1)
deallocate(a3)

return
end
