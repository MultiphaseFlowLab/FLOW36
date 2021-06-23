subroutine get_div2D(nrx,nry,nrz,ax,ay,az,div2D)

use commondata
use wavenumber

integer :: i,k,j

double precision, dimension(nx,nz,ny) :: nrx,nry,nrz,ax,ay,az,div2D
double precision, allocatable, dimension(:,:,:) :: a1,a2,a3
double precision, allocatable, dimension(:,:,:,:) :: a1c,a2c,a3c,a11c,a22c,a33c


! calculate 2D divergence over the interface

! n x u (x: cross product)
allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))

do j=1,ny
 do k=1,nz
  do i=1,nx
   a1(i,k,j)=nry(i,k,j)*az(i,k,j)-nrz(i,k,j)*ay(i,k,j)
   a2(i,k,j)=nrz(i,k,j)*ax(i,k,j)-nrx(i,k,j)*az(i,k,j)
   a3(i,k,j)=nrx(i,k,j)*ay(i,k,j)-nry(i,k,j)*ax(i,k,j)
  enddo
 enddo
enddo

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))

! calculate curl(n x u) (curl of cross product)
call phys_to_spectral(a1,a1c,0)
call phys_to_spectral(a2,a2c,0)
call phys_to_spectral(a3,a3c,0)

deallocate(a1,a2,a3)

allocate(a11c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-ky(j)*a3c(i,:,j,2)
  a11c(i,:,j,2)=+ky(j)*a3c(i,:,j,1)
  a22c(i,:,j,1)=+kx(j)*a3c(i,:,j,2)
  a22c(i,:,j,2)=-kx(j)*a3c(i,:,j,1)
  a33c(i,:,j,1)=-kx(i)*a2c(i,:,j,2)+ky(j)*a1c(i,:,j,2)
  a33c(i,:,j,2)=+kx(i)*a2c(i,:,j,1)-ky(j)*a1c(i,:,j,1)
 enddo
enddo

call dz(a2c,a3c)
call dz(a1c,a2c)

a11c=a11c-a3c
a22c=a22c+a2c

deallocate(a1c,a2c,a3c)

allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))

call spectral_to_phys(a11c,a1,0)
call spectral_to_phys(a22c,a2,0)
call spectral_to_phys(a33c,a3,0)

deallocate(a11c,a22c,a33c)


! n.curl(n x u) (scalar product)
do j=1,ny
 do k=1,nz
  do i=1,nx
   div2d(i,k,j)=nrx(i,k,j)*a1(i,k,j)+nry(i,k,j)*a2(i,k,j)+nrz(i,k,j)*a3(i,k,j)
  enddo
 enddo
enddo

deallocate(a1,a2,a3)

! ! test on a flat 2D interface: -dw/dz
! allocate(wc(nx/2+1,nz,ny,2))
! allocate(a1c(nx/2+1,nz,ny,2))
! allocate(a2(nx,nz,ny))
!
! call phys_to_spectral(w,wc,0)
! call dz(-wc,a1c)
! call spectral_to_phys(a1c,a2,0)
! write(*,*) maxval(w),maxval(a2)
! deallocate(wc,a1c)
!
! open(1234,file='./output/div2D.dat',status='new',form='unformatted',access='stream',convert='little_endian')
! write(1234) a1(:,(nz+1)/2,:)
! close(1234,status='keep')
!
! open(1234,file='./output/div-w.dat',status='new',form='unformatted',access='stream',convert='little_endian')
! write(1234) a2(:,(nz+1)/2,:)
! close(1234,status='keep')
!
! deallocate(a2)
! ! end of test


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_div2D_FG(nrx,nry,nrz,ax,ay,az,div2D)

use commondata
use wavenumber

integer :: i,k,j

double precision, dimension(nxfg,nzfg,nyfg) :: nrx,nry,nrz,ax,ay,az,div2D
double precision, allocatable, dimension(:,:,:) :: a1,a2,a3
double precision, allocatable, dimension(:,:,:,:) :: a1c,a2c,a3c,a11c,a22c,a33c


! calculate 2D divergence over the interface

! n x u (x: cross product)
allocate(a1(nxfg,nzfg,nyfg))
allocate(a2(nxfg,nzfg,nyfg))
allocate(a3(nxfg,nzfg,nyfg))

do j=1,nyfg
 do k=1,nzfg
  do i=1,nxfg
   a1(i,k,j)=nry(i,k,j)*az(i,k,j)-nrz(i,k,j)*ay(i,k,j)
   a2(i,k,j)=nrz(i,k,j)*ax(i,k,j)-nrx(i,k,j)*az(i,k,j)
   a3(i,k,j)=nrx(i,k,j)*ay(i,k,j)-nry(i,k,j)*ax(i,k,j)
  enddo
 enddo
enddo

allocate(a1c(nxfg/2+1,nzfg,nyfg,2))
allocate(a2c(nxfg/2+1,nzfg,nyfg,2))
allocate(a3c(nxfg/2+1,nzfg,nyfg,2))

! calculate curl(n x u) (curl of cross product)
call phys_to_spectral_fg(a1,a1c,0)
call phys_to_spectral_fg(a2,a2c,0)
call phys_to_spectral_fg(a3,a3c,0)

deallocate(a1,a2,a3)

allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
allocate(a22c(nxfg/2+1,nzfg,nyfg,2))
allocate(a33c(nxfg/2+1,nzfg,nyfg,2))

do j=1,nyfg
 do i=1,nxfg/2+1
  a11c(i,:,j,1)=-kyfg(j)*a3c(i,:,j,2)
  a11c(i,:,j,2)=+kyfg(j)*a3c(i,:,j,1)
  a22c(i,:,j,1)=+kxfg(j)*a3c(i,:,j,2)
  a22c(i,:,j,2)=-kxfg(j)*a3c(i,:,j,1)
  a33c(i,:,j,1)=-kxfg(i)*a2c(i,:,j,2)+kyfg(j)*a1c(i,:,j,2)
  a33c(i,:,j,2)=+kxfg(i)*a2c(i,:,j,1)-kyfg(j)*a1c(i,:,j,1)
 enddo
enddo

call dz_fg(a2c,a3c)
call dz_fg(a1c,a2c)

a11c=a11c-a3c
a22c=a22c+a2c

deallocate(a1c,a2c,a3c)

allocate(a1(nxfg,nzfg,nyfg))
allocate(a2(nxfg,nzfg,nyfg))
allocate(a3(nxfg,nzfg,nyfg))

call spectral_to_phys_fg(a11c,a1,0)
call spectral_to_phys_fg(a22c,a2,0)
call spectral_to_phys_fg(a33c,a3,0)

deallocate(a11c,a22c,a33c)


! n.curl(n x u) (scalar product)
do j=1,nyfg
 do k=1,nzfg
  do i=1,nxfg
   div2d(i,k,j)=nrx(i,k,j)*a1(i,k,j)+nry(i,k,j)*a2(i,k,j)+nrz(i,k,j)*a3(i,k,j)
  enddo
 enddo
enddo

deallocate(a1,a2,a3)

return
end
