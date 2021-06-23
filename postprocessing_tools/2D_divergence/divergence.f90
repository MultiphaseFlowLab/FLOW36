subroutine divergence2D(nstep)

use commondata
use fields
use wavenumber

integer :: nstep
integer :: i,j,k

double precision :: modulo
double precision, dimension(nz) :: um,vm,wm
double precision, dimension(:,:,:), allocatable :: nrx,nry,nrz,a1,a2,a3,phix,phiy,phiz,sigx,sigy,sigz,sigma,phifg
double precision, dimension(:,:,:,:), allocatable :: div2d
double precision, dimension(:,:,:,:), allocatable :: a1c,a2c,a3c,a11c,a22c,a33c
character(len=8) :: step



write(step,'(i8.8)') nstep

write(*,*) 'Step ',nstep,' of ',nend

call read_fields(nstep)

! calculate normal to interface
allocate(nrx(nx,nz,ny))
allocate(nry(nx,nz,ny))
allocate(nrz(nx,nz,ny))
allocate(phic(nx/2+1,nz,ny,2))
allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))

call phys_to_spectral(phi,phic,0)

do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=-kx(i)*phic(i,:,j,2)
  a1c(i,:,j,2)=+kx(i)*phic(i,:,j,1)
  a2c(i,:,j,1)=-ky(j)*phic(i,:,j,2)
  a2c(i,:,j,2)=+ky(j)*phic(i,:,j,1)
 enddo
enddo

call spectral_to_phys(a1c,nrx,0)
call spectral_to_phys(a2c,nry,0)

call dz(phic,a1c)

call spectral_to_phys(a1c,nrz,0)

deallocate(a1c,a2c)

! normalize normal to modulo one
do j=1,ny
 do k=1,nz
  do i=1,nx
   modulo=dsqrt(nrx(i,k,j)**2+nry(i,k,j)**2+nrz(i,k,j)**2)
   nrx(i,k,j)=nrx(i,k,j)/modulo
   nry(i,k,j)=nry(i,k,j)/modulo
   nrz(i,k,j)=nrz(i,k,j)/modulo
  enddo
 enddo
enddo

! end computation of normal to interface

allocate(div2d(expx*nx,expz*(nz-1)+1,expy*ny,4))

allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! u' 2D divergence
! get mean velocity
um=0.0d0
vm=0.0d0
wm=0.0d0
do j=1,ny
 do k=1,nz
  do i=1,nx
   um(k)=um(k)+u(i,k,j)
   vm(k)=vm(k)+v(i,k,j)
   wm(k)=wm(k)+w(i,k,j)
  enddo
 enddo
enddo
um=um/dble(nx*ny)
vm=vm/dble(nx*ny)
wm=wm/dble(nx*ny)

do j=1,ny
 do i=1,nx
  a1(i,:,j)=u(i,:,j)-um
  a2(i,:,j)=v(i,:,j)-vm
  a3(i,:,j)=w(i,:,j)-wm
 enddo
enddo

call get_div2D(nrx,nry,nrz,a1,a2,a3,a1)

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call phys_to_spectral(a1,a1c,0)
call coarse2fine(a1c,a11c)
call spectral_to_phys_fg(a11c,div2d(:,:,:,1),0)
deallocate(a1c,a11c)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! viscous stress 2D divergence

allocate(uc(nx/2+1,nz,ny,2))
allocate(vc(nx/2+1,nz,ny,2))
allocate(wc(nx/2+1,nz,ny,2))

call phys_to_spectral(u,uc,0)
call phys_to_spectral(v,vc,0)
call phys_to_spectral(w,wc,0)

allocate(a3c(nx/2+1,nz,ny,2))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))

call dz(uc,a3c)
call dz(a3c,a11c)

call dz(vc,a3c)
call dz(a3c,a22c)

call dz(wc,a3c)
call dz(a3c,a33c)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=a11c(i,:,j,1)-k2(i,j)*uc(i,:,j,1)
  a11c(i,:,j,2)=a11c(i,:,j,2)-k2(i,j)*uc(i,:,j,2)
  a22c(i,:,j,1)=a22c(i,:,j,1)-k2(i,j)*vc(i,:,j,1)
  a22c(i,:,j,2)=a22c(i,:,j,2)-k2(i,j)*vc(i,:,j,2)
  a33c(i,:,j,1)=a33c(i,:,j,1)-k2(i,j)*wc(i,:,j,1)
  a33c(i,:,j,2)=a33c(i,:,j,2)-k2(i,j)*wc(i,:,j,2)
 enddo
enddo

deallocate(a3c)

call spectral_to_phys(a11c,a1,0)
call spectral_to_phys(a22c,a2,0)
call spectral_to_phys(a33c,a3,0)
deallocate(a11c,a22c,a33c)

a1=a1/re
a2=a2/re
a3=a3/re

call get_div2D(nrx,nry,nrz,a1,a2,a3,a1)

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call phys_to_spectral(a1,a1c,0)
call coarse2fine(a1c,a11c)
call spectral_to_phys_fg(a11c,div2d(:,:,:,2),0)
deallocate(a1c,a11c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! shear stress 2D divergence
! uc,vc,wc already available in spectral space

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))

call dz(vc,a1c)
call dz(uc,a2c)

do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=a1c(i,:,j,1)-ky(j)*wc(i,:,j,2)
  a1c(i,:,j,2)=a1c(i,:,j,2)+ky(j)*wc(i,:,j,1)
  a2c(i,:,j,1)=a2c(i,:,j,1)-kx(i)*wc(i,:,j,2)
  a2c(i,:,j,2)=a2c(i,:,j,2)+kx(i)*wc(i,:,j,1)
  a3c(i,:,j,1)=-ky(j)*uc(i,:,j,2)-kx(i)*vc(i,:,j,2)
  a3c(i,:,j,2)=+ky(j)*uc(i,:,j,1)+kx(i)*vc(i,:,j,1)
 enddo
enddo

deallocate(uc,vc,wc)

call spectral_to_phys(a1c,a1,0)
call spectral_to_phys(a2c,a2,0)
call spectral_to_phys(a3c,a3,0)

deallocate(a1c,a2c,a3c)

call get_div2D(nrx,nry,nrz,a1,a2,a3,a1)

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call phys_to_spectral(a1,a1c,0)
call coarse2fine(a1c,a11c)
call spectral_to_phys_fg(a11c,div2d(:,:,:,3),0)
deallocate(a1c,a11c)

deallocate(a1,a2,a3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Marangoni stress 2D divergence
! gradient of sigma has to be calculated on fine grid


allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))

do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=-kx(i)*phic(i,:,j,2)
  a1c(i,:,j,2)=+kx(i)*phic(i,:,j,1)
  a2c(i,:,j,1)=-ky(j)*phic(i,:,j,2)
  a2c(i,:,j,2)=+ky(j)*phic(i,:,j,1)
 enddo
enddo

call dz(phic,a3c)

allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
allocate(a22c(nxfg/2+1,nzfg,nyfg,2))
allocate(a33c(nxfg/2+1,nzfg,nyfg,2))
call coarse2fine(a1c,a11c)
call coarse2fine(a2c,a22c)
call coarse2fine(a3c,a33c)
deallocate(a1c,a2c,a3c)

allocate(phix(nxfg,nzfg,nyfg))
allocate(phiy(nxfg,nzfg,nyfg))
allocate(phiz(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a11c,phix,0)
call spectral_to_phys_fg(a22c,phiy,0)
call spectral_to_phys_fg(a33c,phiz,0)
deallocate(a11c,a22c,a33c)


allocate(a3c(nxfg/2+1,nzfg,nyfg,2))
allocate(sigma(nxfg,nzfg,nyfg))
sigma=max(1.0d0+betas*log(1.0d0-psi),0.5d0)

call phys_to_spectral_fg(sigma,a3c,0)
deallocate(sigma)

allocate(a1c(nxfg/2+1,nzfg,nyfg,2))
allocate(a2c(nxfg/2+1,nzfg,nyfg,2))
do j=1,nyfg
 do i=1,nxfg/2+1
  a1c(i,:,j,1)=-kxfg(i)*a3c(i,:,j,2)
  a1c(i,:,j,2)=+kxfg(i)*a3c(i,:,j,1)
  a2c(i,:,j,1)=-kyfg(j)*a3c(i,:,j,2)
  a2c(i,:,j,2)=+kyfg(j)*a3c(i,:,j,1)
 enddo
enddo

allocate(sigx(nxfg,nzfg,nyfg))
allocate(sigy(nxfg,nzfg,nyfg))
allocate(sigz(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a1c,sigx,0)
call spectral_to_phys_fg(a2c,sigy,0)

call dz_fg(a3c,a1c)
call spectral_to_phys_fg(a1c,sigz,0)

deallocate(a1c,a2c,a3c)

! assemble Marangoni stresses
allocate(a1(nxfg,nzfg,nyfg))
allocate(a2(nxfg,nzfg,nyfg))
allocate(a3(nxfg,nzfg,nyfg))

do j=1,nyfg
 do k=1,nzfg
  do i=1,nxfg
   a1(i,k,j)=+sigx(i,k,j)*(phiy(i,k,j)**2+phiz(i,k,j)**2)-sigy(i,k,j)*phix(i,k,j)*phiy(i,k,j)-sigz(i,k,j)*phix(i,k,j)*phiz(i,k,j)
   a2(i,k,j)=-sigx(i,k,j)*phix(i,k,j)*phiy(i,k,j)+sigy(i,k,j)*(phix(i,k,j)**2+phiz(i,k,j)**2)-sigz(i,k,j)*phiy(i,k,j)*phiz(i,k,j)
   a3(i,k,j)=-sigx(i,k,j)*phix(i,k,j)*phiz(i,k,j)-sigy(i,k,j)*phiy(i,k,j)*phiz(i,k,j)+sigz(i,k,j)*(phix(i,k,j)**2+phiy(i,k,j)**2)
  enddo
 enddo
enddo

a1=3.0d0/dsqrt(8.0d0)*Ch/We*a1
a2=3.0d0/dsqrt(8.0d0)*Ch/We*a2
a3=3.0d0/dsqrt(8.0d0)*Ch/We*a3

deallocate(phix,phiy,phiz,sigx,sigy,sigz)


! take normals to fine grid
allocate(a1c(nx/2+1,nz,ny,2))
call phys_to_spectral(nrx,a1c,0)
deallocate(nrx)
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call coarse2fine(a1c,a11c)
deallocate(a1c)
allocate(nrx(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a11c,nrx,0)
deallocate(a11c)

allocate(a1c(nx/2+1,nz,ny,2))
call phys_to_spectral(nry,a1c,0)
deallocate(nry)
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call coarse2fine(a1c,a11c)
deallocate(a1c)
allocate(nry(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a11c,nry,0)
deallocate(a11c)

allocate(a1c(nx/2+1,nz,ny,2))
call phys_to_spectral(nrz,a1c,0)
deallocate(nrz)
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
call coarse2fine(a1c,a11c)
deallocate(a1c)
allocate(nrz(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a11c,nrz,0)
deallocate(a11c)

call get_div2D_FG(nrx,nry,nrz,a1,a2,a3,div2d(:,:,:,4))


deallocate(nrx,nry,nrz,a1,a2,a3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! correlation div_2D and surfactant concentration at the interface
! div2d contains 2D divergence of u', viscous stress, shear stress, Marangoni stress in FG


! take phi in fine grid (same grid as surfactant)
allocate(a1c(nxfg/2+1,nzfg,nyfg,2))
call coarse2fine(phic,a1c)
deallocate(phic)
allocate(phifg(nxfg,nzfg,nyfg))
call spectral_to_phys_fg(a1c,phifg,0)
deallocate(a1c)


open(256,file='./output/correlation_'//step//'.dat',form='formatted',status='new')
write(256,'(5(a20))') 'psi','div 2D up','div 2D visc stress','div 2D shear stress','div 2D Mar stress'

do j=1,expy*ny
 do k=1,expz*(nz-1)+1
  do i=1,expx*nx
   if(dabs(phifg(i,k,j)).lt.0.2d0) write(256,'(5(f20.8))') psi(i,k,j),div2d(i,k,j,1),div2d(i,k,j,2),div2d(i,k,j,3),div2d(i,k,j,4)
  enddo
 enddo
enddo

close(256,status='keep')


deallocate(div2d,phifg)

return
end
