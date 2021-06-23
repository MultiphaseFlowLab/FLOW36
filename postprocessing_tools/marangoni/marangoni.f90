subroutine marangoni(nstep)

use commondata
use fields
use wavenumber

integer :: nstep
integer :: i,j,k

double precision :: modulo,thu,thup,thsh,thvis,cx,cy,cz
double precision, dimension(nz) :: um,vm,wm
double precision, dimension(:,:,:), allocatable :: phix,phiy,phiz,a11,a22,a33,a12,a13,a23
double precision, dimension(:,:,:), allocatable :: sfx,sfy,sfz,nrx,nry,nrz,shx,shy,shz,visx,visy,visz
double precision, dimension(:,:,:), allocatable :: sigma
double precision, dimension(:,:,:,:), allocatable :: a11c,a22c,a33c,a1c,a2c,a3c
character(len=8) :: step



write(step,'(i8.8)') nstep

write(*,*) 'Step ',nstep,' of ',nend

call read_fields(nstep)


! calculate Korteweg stress tensor (coarse grid)
! phi gradients
allocate(phic(nx/2+1,nz,ny,2))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))

call phys_to_spectral(phi,phic,0)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*phic(i,:,j,2)
  a11c(i,:,j,2)=+kx(i)*phic(i,:,j,1)
  a22c(i,:,j,1)=-ky(j)*phic(i,:,j,2)
  a22c(i,:,j,2)=+ky(j)*phic(i,:,j,1)
 enddo
enddo
call dz(phic,a33c)

allocate(phix(nx,nz,ny))
allocate(phiy(nx,nz,ny))
allocate(phiz(nx,nz,ny))

call spectral_to_phys(a11c,phix,0)
call spectral_to_phys(a22c,phiy,0)
call spectral_to_phys(a33c,phiz,0)

deallocate(phic,a11c,a22c,a33c)

! assemble Korteweg stress tensor
allocate(a11(nx,nz,ny))
allocate(a12(nx,nz,ny))
allocate(a13(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a23(nx,nz,ny))
allocate(a33(nx,nz,ny))

! allocate normal to interface, needed later
allocate(nrx(nx,nz,ny))
allocate(nry(nx,nz,ny))
allocate(nrz(nx,nz,ny))

do j=1,ny
 do k=1,nz
  do i=1,nx
   modulo=dsqrt(phix(i,k,j)**2+phiy(i,k,j)**2+phiz(i,k,j)**2)
   nrx(i,k,j)=phix(i,k,j)/modulo
   nry(i,k,j)=phiy(i,k,j)/modulo
   nrz(i,k,j)=phiz(i,k,j)/modulo
   a11(i,k,j)=phiy(i,k,j)**2+phiz(i,k,j)**2
   a12(i,k,j)=-phix(i,k,j)*phiy(i,k,j)
   a13(i,k,j)=-phix(i,k,j)*phiz(i,k,j)
   a22(i,k,j)=phix(i,k,j)**2+phiz(i,k,j)**2
   a23(i,k,j)=-phiy(i,k,j)*phiz(i,k,j)
   a33(i,k,j)=phix(i,k,j)**2+phiy(i,k,j)**2
  enddo
 enddo
enddo

deallocate(phix,phiy,phiz)


! calculate surface tension gradient (fine grid)
allocate(sigma(nx*expx,(nz-1)*expz+1,ny*expy))

if(psi_flag.eq.1)then
 sigma=1.0d0+betas*log(1.0d0-psi)
else
 sigma=1.0d0
endif

allocate(phic(nx*expx/2+1,(nz-1)*expz+1,ny*expy,2))
call phys_to_spectral_fg(sigma,phic,0)
deallocate(sigma)

allocate(a11c(nx*expx/2+1,(nz-1)*expz+1,ny*expy,2))
allocate(a22c(nx*expx/2+1,(nz-1)*expz+1,ny*expy,2))
allocate(a33c(nx*expx/2+1,(nz-1)*expz+1,ny*expy,2))

do j=1,expy*ny
 do i=1,expx*nx/2+1
  a11c(i,:,j,1)=-kxfg(i)*phic(i,:,j,2)
  a11c(i,:,j,2)=+kxfg(i)*phic(i,:,j,1)
  a22c(i,:,j,1)=-kyfg(j)*phic(i,:,j,2)
  a22c(i,:,j,2)=+kyfg(j)*phic(i,:,j,1)
 enddo
enddo
call dz_fg(phic,a33c)

deallocate(phic)

! a11c,a22c,a33c: gradients in spectral space (fine grid)
allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))
call fine2coarse(a11c,a1c)
call fine2coarse(a22c,a2c)
call fine2coarse(a33c,a3c)
deallocate(a11c,a22c,a33c)

allocate(phix(nx,nz,ny))
allocate(phiy(nx,nz,ny))
allocate(phiz(nx,nz,ny))
call spectral_to_phys(a1c,phix,0)
call spectral_to_phys(a2c,phiy,0)
call spectral_to_phys(a3c,phiz,0)
deallocate(a1c,a2c,a3c)

! Korteweg stress tensor: a11,a12,a13,a22,a23,a33 (CG)
! surface tension gradient: phix,phiy,phiz (CG)

! calculate surface forces
allocate(sfx(nx,nz,ny))
allocate(sfy(nx,nz,ny))
allocate(sfz(nx,nz,ny))

do j=1,ny
 do k=1,nz
  do i=1,nx
   sfx(i,k,j)=a11(i,k,j)*phix(i,k,j)+a12(i,k,j)*phiy(i,k,j)+a13(i,k,j)*phiz(i,k,j)
   sfy(i,k,j)=a12(i,k,j)*phix(i,k,j)+a22(i,k,j)*phiy(i,k,j)+a23(i,k,j)*phiz(i,k,j)
   sfz(i,k,j)=a13(i,k,j)*phix(i,k,j)+a23(i,k,j)*phiy(i,k,j)+a33(i,k,j)*phiz(i,k,j)
   ! normalize to modulo 1
   modulo=dsqrt(sfx(i,k,j)**2+sfy(i,k,j)**2+sfz(i,k,j)**2)
   sfx(i,k,j)=sfx(i,k,j)/modulo
   sfy(i,k,j)=sfy(i,k,j)/modulo
   sfz(i,k,j)=sfz(i,k,j)/modulo
  enddo
 enddo
enddo

deallocate(phix,phiy,phiz)
deallocate(a11,a22,a33,a12,a13,a23)

! surface forces: sfx,sfy,sfz (CG)


! calculate correlations with Marangoni stresses unit-length vector

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


! scalar product Marangoni and:
! 1) tangential velocity (total and prime)
! 2) tangential viscous stresses (total and prime?)


! assemble shear stress vector
allocate(uc(nx/2+1,nz,ny,2))
allocate(vc(nx/2+1,nz,ny,2))
allocate(wc(nx/2+1,nz,ny,2))

call phys_to_spectral(u,uc,0)
call phys_to_spectral(v,vc,0)
call phys_to_spectral(w,wc,0)

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

allocate(shx(nx,nz,ny))
allocate(shy(nx,nz,ny))
allocate(shz(nx,nz,ny))

call spectral_to_phys(a1c,shx,0)
call spectral_to_phys(a2c,shy,0)
call spectral_to_phys(a3c,shz,0)

deallocate(a1c,a2c,a3c)


! assemble viscous stress vector

allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))
allocate(a11c(nx/2+1,nz,ny,2))

call dz(uc,a11c)
call dz(a11c,a1c)

call dz(vc,a11c)
call dz(a11c,a2c)

call dz(wc,a11c)
call dz(a11c,a3c)

deallocate(a11c)

do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=a1c(i,:,j,1)-k2(i,j)*uc(i,:,j,1)
  a1c(i,:,j,2)=a1c(i,:,j,2)-k2(i,j)*uc(i,:,j,2)
  a2c(i,:,j,1)=a2c(i,:,j,1)-k2(i,j)*vc(i,:,j,1)
  a2c(i,:,j,2)=a2c(i,:,j,2)-k2(i,j)*vc(i,:,j,2)
  a3c(i,:,j,1)=a3c(i,:,j,1)-k2(i,j)*wc(i,:,j,1)
  a3c(i,:,j,2)=a3c(i,:,j,2)-k2(i,j)*wc(i,:,j,2)
 enddo
enddo

deallocate(uc,vc,wc)

allocate(visx(nx,nz,ny))
allocate(visy(nx,nz,ny))
allocate(visz(nx,nz,ny))

call spectral_to_phys(a1c,visx,0)
call spectral_to_phys(a2c,visy,0)
call spectral_to_phys(a3c,visz,0)

deallocate(a1c,a2c,a3c)



! Marangoni force: sfx,sfy,sfz
! shear stresses: shx,shy,shz
! normal to interface: nrx,nry,nrz

! component tangential to interface:  V-(n.V)V

open(632,file='./output/correlation_'//step//'.dat',form='formatted',status='new')
write(632,'(4(a16),a)') 'Mar * u_t','Mar * up_t','Mar * sh_t','Mar * vis_t', &
 &  '   (if no values, Marangoni stresses are zero everywhere)'

do j=1,ny
 do k=1,nz
  do i=1,nx
   ! correlation only at interface
   if((dabs(phi(i,k,j)).le.0.2d0).and.(isnan(sfz(i,k,j)).eqv..false.))then
    ! get components tangential to interface
    ! u (total)
    cx=u(i,k,j)-(u(i,k,j)*nrx(i,k,j)+v(i,k,j)*nry(i,k,j)+w(i,k,j)*nrz(i,k,j))*nrx(i,k,j)
    cy=v(i,k,j)-(u(i,k,j)*nrx(i,k,j)+v(i,k,j)*nry(i,k,j)+w(i,k,j)*nrz(i,k,j))*nry(i,k,j)
    cz=w(i,k,j)-(u(i,k,j)*nrx(i,k,j)+v(i,k,j)*nry(i,k,j)+w(i,k,j)*nrz(i,k,j))*nrz(i,k,j)
    modulo=dsqrt(cx**2+cy**2+cz**2)
    thu=(sfx(i,k,j)*cx+sfy(i,k,j)*cy+sfz(i,k,j)*cz)/modulo
    ! u' (fluctuations)
    cx=(u(i,k,j)-um(k))-((u(i,k,j)-um(k))*nrx(i,k,j)+(v(i,k,j)-vm(k))*nry(i,k,j)+(w(i,k,j)-wm(k))*nrz(i,k,j))*nrx(i,k,j)
    cy=(v(i,k,j)-vm(k))-((u(i,k,j)-um(k))*nrx(i,k,j)+(v(i,k,j)-vm(k))*nry(i,k,j)+(w(i,k,j)-wm(k))*nrz(i,k,j))*nry(i,k,j)
    cz=(w(i,k,j)-wm(k))-((u(i,k,j)-um(k))*nrx(i,k,j)+(v(i,k,j)-vm(k))*nry(i,k,j)+(w(i,k,j)-wm(k))*nrz(i,k,j))*nrz(i,k,j)
    modulo=dsqrt(cx**2+cy**2+cz**2)
    thup=(sfx(i,k,j)*cx+sfy(i,k,j)*cy+sfz(i,k,j)*cz)/modulo
    ! shear
    cx=shx(i,k,j)-(shx(i,k,j)*nrx(i,k,j)+shy(i,k,j)*nry(i,k,j)+shz(i,k,j)*nrz(i,k,j))*nrx(i,k,j)
    cy=shy(i,k,j)-(shx(i,k,j)*nrx(i,k,j)+shy(i,k,j)*nry(i,k,j)+shz(i,k,j)*nrz(i,k,j))*nry(i,k,j)
    cz=shz(i,k,j)-(shx(i,k,j)*nrx(i,k,j)+shy(i,k,j)*nry(i,k,j)+shz(i,k,j)*nrz(i,k,j))*nrz(i,k,j)
    modulo=dsqrt(cx**2+cy**2+cz**2)
    thsh=(sfx(i,k,j)*cx+sfy(i,k,j)*cy+sfz(i,k,j)*cz)/modulo
    ! viscous
    cx=visx(i,k,j)-(visx(i,k,j)*nrx(i,k,j)+visy(i,k,j)*nry(i,k,j)+visz(i,k,j)*nrz(i,k,j))*nrx(i,k,j)
    cy=visy(i,k,j)-(visx(i,k,j)*nrx(i,k,j)+visy(i,k,j)*nry(i,k,j)+visz(i,k,j)*nrz(i,k,j))*nry(i,k,j)
    cz=visz(i,k,j)-(visx(i,k,j)*nrx(i,k,j)+visy(i,k,j)*nry(i,k,j)+visz(i,k,j)*nrz(i,k,j))*nrz(i,k,j)
    modulo=dsqrt(cx**2+cy**2+cz**2)
    thvis=(sfx(i,k,j)*cx+sfy(i,k,j)*cy+sfz(i,k,j)*cz)/modulo
    ! write to file
    write(632,'(4(f16.8))') thu,thup,thsh,thvis
   endif
  enddo
 enddo
enddo

close(632,status='keep')

deallocate(shx,shy,shz)
deallocate(sfx,sfy,sfz)
deallocate(visx,visy,visz)
deallocate(nrx,nry,nrz)

return
end
