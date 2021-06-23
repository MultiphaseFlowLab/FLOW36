subroutine calc_pressure(step)

use commondata
use fields
use wavenumber

integer :: step
integer :: i,j

double precision, dimension(nx,nz,ny) :: rho,eta
double precision, dimension(nx/2+1,nz,ny,2) :: stermc
double precision, dimension(nxfg,nzfg,nyfg) :: sigma
double precision, dimension(nxfg/2+1,nzfg,nyfg,2) :: phicfg
double precision, dimension(nx,ny) :: bcp
double precision, dimension(nx/2+1,ny,2) :: bcc_p,bcc_m
double precision, dimension(nx/2+1,ny,2,2) :: bcc
double precision, dimension(2) :: p,q,zp
double precision, allocatable, dimension(:,:,:) :: a11,a12,a13,a22,a33,a23, a1,a2,a3
double precision, allocatable, dimension(:,:,:,:) :: a11c,a12c,a13c,a22c,a33c,a23c, a1c,a2c,a3c


write(*,'(1x,a,i8,a,i8)') 'Step ',step,' of ',nend

call read_fields(step)

! u,v,w,phi and psi available in physical and spectral space
if(phi_flag.eq.1)then
  rho=1.0d0+(rhor-1.0d0)/2.0d0*(phi+1.0d0)
  eta=1.0d0+(visr-1.0d0)/2.0d0*(phi+1.0d0)
else
  rho=1.0d0
  eta=1.0d0
endif

! assemble S term for calculation of pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! viscous term
allocate(a11(nx,nz,ny))
allocate(a12(nx,nz,ny))
allocate(a13(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a33(nx,nz,ny))
allocate(a23(nx,nz,ny))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))
allocate(a13c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))
allocate(a23c(nx/2+1,nz,ny,2))

call dz(uc,a13c)
call dz(vc,a23c)
call dz(wc,a33c)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*uc(i,:,j,2)
  a11c(i,:,j,2)=+kx(i)*uc(i,:,j,1)
  a22c(i,:,j,1)=-ky(j)*vc(i,:,j,2)
  a22c(i,:,j,2)=+ky(j)*vc(i,:,j,1)
  a12c(i,:,j,1)=-ky(j)*uc(i,:,j,2)-kx(i)*vc(i,:,j,2)
  a12c(i,:,j,2)=+ky(j)*uc(i,:,j,1)+kx(i)*vc(i,:,j,1)
  a13c(i,:,j,1)=a13c(i,:,j,1)-kx(i)*wc(i,:,j,2)
  a13c(i,:,j,2)=a13c(i,:,j,2)+kx(i)*wc(i,:,j,1)
  a23c(i,:,j,1)=a23c(i,:,j,1)-ky(j)*wc(i,:,j,2)
  a23c(i,:,j,2)=a23c(i,:,j,2)+ky(j)*wc(i,:,j,1)
 enddo
enddo
a11c=2.0d0*a11c
a22c=2.0d0*a22c
a33c=2.0d0*a33c

if((phi_flag.eq.1).and.(matchedvis.ne.1))then
  call spectral_to_phys(a11c,a11,0)
  call spectral_to_phys(a12c,a12,0)
  call spectral_to_phys(a13c,a13,0)
  call spectral_to_phys(a22c,a22,0)
  call spectral_to_phys(a23c,a23,0)
  call spectral_to_phys(a33c,a33,0)

  a11=eta*a11
  a12=eta*a12
  a13=eta*a13
  a22=eta*a22
  a23=eta*a23
  a33=eta*a33

  call phys_to_spectral(a11,a11c,0)
  call phys_to_spectral(a12,a12c,0)
  call phys_to_spectral(a13,a13c,0)
  call phys_to_spectral(a22,a22c,0)
  call phys_to_spectral(a23,a23c,0)
  call phys_to_spectral(a33,a33c,0)
endif

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))

! divergence of stress tensor
call dz(a13c,a1c)
call dz(a23c,a2c)
call dz(a33c,a3c)

do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=a1c(i,:,j,1)-kx(i)*a11c(i,:,j,2)-ky(j)*a12c(i,:,j,2)
  a1c(i,:,j,2)=a1c(i,:,j,2)+kx(i)*a11c(i,:,j,1)+ky(j)*a12c(i,:,j,1)
  a2c(i,:,j,1)=a2c(i,:,j,1)-kx(i)*a12c(i,:,j,2)-ky(j)*a22c(i,:,j,2)
  a2c(i,:,j,2)=a2c(i,:,j,2)+kx(i)*a12c(i,:,j,1)+ky(j)*a22c(i,:,j,1)
  a3c(i,:,j,1)=a3c(i,:,j,1)-kx(i)*a13c(i,:,j,2)-ky(j)*a23c(i,:,j,2)
  a3c(i,:,j,2)=a3c(i,:,j,2)+kx(i)*a13c(i,:,j,1)+ky(j)*a23c(i,:,j,1)
 enddo
enddo

! create S term (scalar S term, divergence of vector S term)
! divergence of stress vector
call dz(a3c,stermc)

do j=1,ny
 do i=1,nx/2+1
  stermc(i,:,j,1)=stermc(i,:,j,1)-kx(i)*a1c(i,:,j,2)-ky(j)*a2c(i,:,j,2)
  stermc(i,:,j,2)=stermc(i,:,j,2)+kx(i)*a1c(i,:,j,1)+ky(j)*a2c(i,:,j,1)
 enddo
enddo

stermc=stermc/re

deallocate(a11,a12,a13,a22,a33,a23)
deallocate(a11c,a12c,a13c,a22c,a33c,a23c)
deallocate(a1c,a2c,a3c)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convective term
allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))
allocate(a11c(nx/2+1,nz,ny,2)) ! dudx
allocate(a12c(nx/2+1,nz,ny,2)) ! dudy
allocate(a13c(nx/2+1,nz,ny,2)) ! dudz
allocate(a22c(nx/2+1,nz,ny,2)) ! dvdx
allocate(a33c(nx/2+1,nz,ny,2)) ! dvdy
allocate(a23c(nx/2+1,nz,ny,2)) ! dvdz
allocate(a1c(nx/2+1,nz,ny,2))  ! dwdx
allocate(a2c(nx/2+1,nz,ny,2))  ! dwdy
allocate(a3c(nx/2+1,nz,ny,2))  ! dwdz

call dz(uc,a13c)
call dz(vc,a23c)
call dz(wc,a3c)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*uc(i,:,j,2)
  a11c(i,:,j,2)=+kx(i)*uc(i,:,j,1)
  a22c(i,:,j,1)=-kx(i)*vc(i,:,j,2)
  a22c(i,:,j,2)=+kx(i)*vc(i,:,j,1)
  a1c(i,:,j,1)=-kx(i)*wc(i,:,j,2)
  a1c(i,:,j,2)=+kx(i)*wc(i,:,j,1)
  a12c(i,:,j,1)=-ky(j)*uc(i,:,j,2)
  a12c(i,:,j,2)=+ky(j)*uc(i,:,j,1)
  a33c(i,:,j,1)=-ky(j)*vc(i,:,j,2)
  a33c(i,:,j,2)=+ky(j)*vc(i,:,j,1)
  a2c(i,:,j,1)=-ky(j)*wc(i,:,j,2)
  a2c(i,:,j,2)=+ky(j)*wc(i,:,j,1)
 enddo
enddo

call spectral_to_phys(a11c,a1,0)
call spectral_to_phys(a33c,a2,0)
call spectral_to_phys(a3c,a3,0)

a1=a1*a1+a2*a2+a3*a3

call spectral_to_phys(a12c,a2,0)
call spectral_to_phys(a22c,a3,0)
a1=a1+2.0d0*a2*a3

call spectral_to_phys(a13c,a2,0)
call spectral_to_phys(a1c,a3,0)
a1=a1+2.0d0*a2*a3

call spectral_to_phys(a23c,a2,0)
call spectral_to_phys(a2c,a3,0)
a1=a1+2.0d0*a2*a3

if((phi_flag.eq.1).and.(matchedrho.ne.1)) a1=rho*a1

call phys_to_spectral(a1,a1c,0)
! sum to S term
stermc=stermc-a1c

deallocate(a1,a2,a3)
deallocate(a11c,a12c,a13c,a22c,a23c,a33c,a1c,a2c,a3c)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interface term (calculated in FG, then ported to CG before adding it to S term)
if(phi_flag.eq.1)then
 allocate(a1(nxfg,nzfg,nyfg))
 allocate(a2(nxfg,nzfg,nyfg))
 allocate(a3(nxfg,nzfg,nyfg))
 allocate(a11(nxfg,nzfg,nyfg))
 allocate(a1c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a2c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a3c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a12c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a13c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a22c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a33c(nxfg/2+1,nzfg,nyfg,2))
 allocate(a23c(nxfg/2+1,nzfg,nyfg,2))

 ! port phi to phifg (FG)
 call coarse2fine(phic,phicfg)

 do j=1,nyfg
  do i=1,nxfg/2+1
   a1c(i,:,j,1)=-kxfg(i)*phicfg(i,:,j,2)
   a1c(i,:,j,2)=+kxfg(i)*phicfg(i,:,j,1)
   a2c(i,:,j,1)=-kyfg(j)*phicfg(i,:,j,2)
   a2c(i,:,j,2)=+kyfg(j)*phicfg(i,:,j,1)
  enddo
 enddo
 call dz_fg(phicfg,a3c)

 call spectral_to_phys_fg(a1c,a1,0)
 call spectral_to_phys_fg(a2c,a2,0)
 call spectral_to_phys_fg(a3c,a3,0)

 ! sigma
 if(psi_flag.eq.1)then
  sigma=max((1.0d0+betas*log(1.0d0-psi)),0.5d0)
 else
  sigma=1.0d0
 endif

 ! assemble Korteweg tensor and multiply by sigma
 a11=(a2*a2+a3*a3)*sigma
 call phys_to_spectral_fg(a11,a11c,0)
 a11=-a1*a2*sigma
 call phys_to_spectral_fg(a11,a12c,0)
 a11=-a1*a3*sigma
 call phys_to_spectral_fg(a11,a13c,0)
 a11=(a1*a1+a3*a3)*sigma
 call phys_to_spectral_fg(a11,a22c,0)
 a11=-a2*a3*sigma
 call phys_to_spectral_fg(a11,a23c,0)
 a11=(a1*a1+a2*a2)*sigma
 call phys_to_spectral_fg(a11,a33c,0)

 deallocate(a1,a2,a3,a11)

 ! divergenge of K*sigma tensor
 call dz_fg(a13c,a1c)
 call dz_fg(a23c,a2c)
 call dz_fg(a33c,a3c)

 do j=1,nyfg
  do i=1,nxfg/2+1
   a1c(i,:,j,1)=a1c(i,:,j,1)-kxfg(i)*a11c(i,:,j,2)-kyfg(j)*a12c(i,:,j,2)
   a1c(i,:,j,2)=a1c(i,:,j,2)+kxfg(i)*a11c(i,:,j,1)+kyfg(j)*a12c(i,:,j,1)
   a2c(i,:,j,1)=a2c(i,:,j,1)-kxfg(i)*a12c(i,:,j,2)-kyfg(j)*a22c(i,:,j,2)
   a2c(i,:,j,2)=a2c(i,:,j,2)+kxfg(i)*a12c(i,:,j,1)+kyfg(j)*a22c(i,:,j,1)
   a3c(i,:,j,1)=a3c(i,:,j,1)-kxfg(i)*a13c(i,:,j,2)-kyfg(j)*a23c(i,:,j,2)
   a3c(i,:,j,2)=a3c(i,:,j,2)+kxfg(i)*a13c(i,:,j,1)+kyfg(j)*a23c(i,:,j,1)
  enddo
 enddo

 deallocate(a12c,a13c,a22c,a23c,a33c)

 ! divergence of surface forces (vector)
 call dz_fg(a3c,a11c)
 do j=1,nyfg
  do i=1,nxfg/2+1
   a11c(i,:,j,1)=a11c(i,:,j,1)-kxfg(i)*a1c(i,:,j,2)-kyfg(j)*a2c(i,:,j,2)
   a11c(i,:,j,2)=a11c(i,:,j,2)+kxfg(i)*a1c(i,:,j,1)+kyfg(j)*a2c(i,:,j,1)
  enddo
 enddo

 ! port from FG to CG
 allocate(a12c(nx/2+1,nz,ny,2))
 call fine2coarse(a11c,a12c)

 ! sum to S term (scalar)
 stermc=stermc+3.0d0/dsqrt(8.0d0)*Ch/We*a12c

 deallocate(a11c,a12c,a1c,a2c,a3c)
endif
! end of S term calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! assemble BC for Helmholtz problem
allocate(a1c(nx/2+1,nz,ny,2))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))
allocate(a13c(nx/2+1,nz,ny,2))
allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))

call phys_to_spectral(eta,a11c,0)

call dz(uc,a12c)
do j=1,ny
 do i=1,nx/2+1
  a13c(i,:,j,1)=-kx(i)*a11c(i,:,j,2)
  a13c(i,:,j,2)=+kx(i)*a11c(i,:,j,1)
 enddo
enddo
call spectral_to_phys(a12c,a2,0)
call spectral_to_phys(a13c,a3,0)
a1=a2*a3

call dz(vc,a12c)
do j=1,ny
 do i=1,nx/2+1
  a13c(i,:,j,1)=-ky(j)*a11c(i,:,j,2)
  a13c(i,:,j,2)=+ky(j)*a11c(i,:,j,1)
 enddo
enddo
call spectral_to_phys(a12c,a2,0)
call spectral_to_phys(a13c,a3,0)
a1=a1+a2*a3

call dz(wc,a12c)
call dz(a12c,a13c)
call spectral_to_phys(a13c,a2,0)
! BC in physical space (entire domain, apply to z=+1, z=-1)
a1=a1+eta*a2

! BC z=+1 --> k=1   --> bcc_p (complex space)
! BC z=-1 --> k=nz  --> bcc_m (complex space)
bcp=a1(:,1,:)
call phys_to_spectral_2D(bcp,bcc_p,0)
bcp=a1(:,nz,:)
call phys_to_spectral_2D(bcp,bcc_m,0)


deallocate(a1c,a11c,a12c,a13c)
deallocate(a1,a2,a3)

! solve Helmholtz problem
! stermc, with BC in bcc_p, bcc_m
bcc(:,:,:,1)=bcc_m
bcc(:,:,:,2)=bcc_p

bcc=bcc/re

! BC are:
! p*var+q*d(var)/dz=r
p(1)=0.0d0
p(2)=0.0d0
q(1)=1.0d0
q(2)=1.0d0
zp(1)=-1.0d0
zp(2)=+1.0d0
call helmholtz(stermc,k2,p,q,bcc,zp)

! mode zero, kx=ky=k2=0
allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
allocate(a22c(nxfg/2+1,nzfg,nyfg,2))
allocate(a11(nxfg,nzfg,nyfg))
allocate(a22(nxfg,nzfg,nyfg))
allocate(a12c(nx/2+1,nz,ny,2))
do j=1,nyfg
 do i=1,nxfg/2+1
  a11c(i,:,j,1)=-kxfg(i)*phicfg(i,:,j,2)
  a11c(i,:,j,2)=+kxfg(i)*phicfg(i,:,j,1)
  a22c(i,:,j,1)=-kyfg(j)*phicfg(i,:,j,2)
  a22c(i,:,j,2)=+kyfg(j)*phicfg(i,:,j,1)
 enddo
enddo

call spectral_to_phys_fg(a11c,a11,0)
call spectral_to_phys_fg(a22c,a22,0)

a11=(a11*a11+a22*a22)*sigma

call phys_to_spectral_fg(a11,a11c,0)
call fine2coarse(a11c,a12c)

allocate(a1(nx,nz,ny))
allocate(a1c(nx/2+1,nz,ny,2))
a1=rho*w*w
call phys_to_spectral(a1,a1c,0)

stermc(1,:,1,1)=3.0d0/dsqrt(8.0d0)*Ch/We*a12c(1,:,1,1)-a1c(1,:,1,1)
stermc(1,:,1,2)=0.0d0
deallocate(a1c,a11c,a22c,a12c)
deallocate(a1,a11,a22)

call spectral_to_phys(stermc,press,0)

! calculate pressure statistics
if(statflag.eq.1)then
  call stats(step)
endif

! generate Paraview file with pressure field
if(paraflag.eq.1)then
  call generate_output(step)
endif

! calculate TKE energy budget
if(budflag.eq.1)then
  call calc_budget(step)
endif



return
end
