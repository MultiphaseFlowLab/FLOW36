subroutine calc_curvature(nstep)

use commondata
use wavenumbers
use pdf_calc

! double precision, dimension(nx/2+1,nz,ny,2) :: sphix,sphiy,sphiz,tmp1
double precision :: gradm
double precision, dimension(nx,nz,ny) :: phix,phiy,phiz

integer :: i,j,k

character(len=40) :: namefile
character(len=8) :: time


! calculate gradient with spectral method
! call spectral_to_phys(phic,phi,1)
!
!
! call dz(phic,sphiz)
!
! call spectral_to_phys(sphiz,phiz,1)
!
! ! calculate d phi/dx and d phi/dy
! do j=1,ny
!   do k=1,nz
!     do i=1,nx/2+1
!       sphix(i,k,j,1)=-kx(i)*phic(i,k,j,2)
!       sphix(i,k,j,2)=kx(i)*phic(i,k,j,1)
!       sphiy(i,k,j,1)=-ky(j)*phic(i,k,j,2)
!       sphiy(i,k,j,2)=ky(j)*phic(i,k,j,1)
!     enddo
!   enddo
! enddo
!
! call spectral_to_phys(sphix,phix,1)
! call spectral_to_phys(sphiy,phiy,1)


!  calculate gradient with central finite differences
phix=0.0d0
phiy=0.0d0
phiz=0.0d0
do j=2,ny-1
  do k=2,nz-1
    do i=2,nx-1
      phix(i,k,j)=(phi(i+1,k,j)-phi(i-1,k,j))/(x(i+1)-x(i-1))
      phiy(i,k,j)=(phi(i,k,j+1)-phi(i,k,j-1))/(y(j+1)-y(j-1))
      phiz(i,k,j)=(phi(i,k+1,j)-phi(i,k-1,j))/(z(k-1)-z(k+1))
    enddo
  enddo
enddo





do j=1,ny
  do k=1,nz
    do i=1,nx
      gradm=((phix(i,k,j))**2+(phiy(i,k,j))**2+(phiz(i,k,j))**2)**0.5
      ! gradm contains the modulus of the gradient of phi
      ! normalized gradient of phi
      phix(i,k,j)=phix(i,k,j)/gradm
      phiy(i,k,j)=phiy(i,k,j)/gradm
      phiz(i,k,j)=phiz(i,k,j)/gradm
      if(abs(phi(i,k,j)).gt.0.8)then
        phix(i,k,j)=0.0d0
        phiy(i,k,j)=0.0d0
        phiz(i,k,j)=0.0d0
      endif
    enddo
  enddo
enddo



! with central finite differences
kv=0.0d0
do j=2,ny-1
  do k=2,nz-1
    do i=2,nx-1
      kv(i,k,j)=(phix(i+1,k,j)-phix(i-1,k,j))/(x(i+1)-x(i-1)) &
 &             +(phiy(i,k,j+1)-phiy(i,k,j-1))/(y(j+1)-y(j-1)) &
 &             +(phiz(i,k+1,j)-phiz(i,k-1,j))/(z(k-1)-z(k+1))
    enddo
  enddo
enddo


! with spectral method
! call phys_to_spectral(phix,sphix,1)
! call phys_to_spectral(phiy,sphiy,1)
! call phys_to_spectral(phiz,sphiz,1)
!
!
! call dz(sphiz,tmp1)
!
!
! ! calculate curvature
! do j=1,ny
!   do k=1,nz
!     do i=1,nx/2+1
!       tmp1(i,k,j,1)=tmp1(i,k,j,1)-kx(i)*sphix(i,k,j,2)-ky(j)*sphiy(i,k,j,2)
!       tmp1(i,k,j,2)=tmp1(i,k,j,2)+kx(i)*sphix(i,k,j,1)+ky(j)*sphiy(i,k,j,1)
!     enddo
!   enddo
! enddo
!
! call spectral_to_phys(tmp1,kv,1)



write(time,'(i8.8)') nstep
namefile='./output/k_'//time//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='new',convert='little_endian')
write(666) kv
close(666,status='keep')



! get max and min for curvature calculation
do j=1,ny
  do k=1,nz
    do i=1,nx
! phi rescaled by a factor 50 to reduce oscillation in gradient calculation
! rescaling does not affect normal calculation (normal is normalized to modulo 1)
      if(abs(phi(i,k,j)).le.threshold/50.0d0)then
       mink=min(mink,kv(i,k,j))
       maxk=max(maxk,kv(i,k,j))
      endif
    enddo
  enddo
enddo


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_pdf(nstep)

use commondata
use pdf_calc

double precision :: bin

integer :: nstep
integer :: i,j,k,q


! allocate array if first step
if(nstep.eq.nstart)then
  allocate(pdf(nset-1))
  allocate(axis(nset))
  pdf=0
  bin=(maxk-mink)/dble(nset-1)
  axis(1)=mink
  do i=2,nset
   axis(i)=axis(i-1)+bin
  enddo
endif


do j=1,ny
 do k=1,nz
  do i=1,nx
   if(abs(phi(i,k,j)).lt.threshold)then
    do q=1,nset-1
     if((kv(i,k,j).ge.axis(q)).and.(kv(i,k,j).lt.axis(q+1))) pdf(q)=pdf(q)+1
    enddo
   endif
  enddo
 enddo
enddo



if(nstep.eq.nend)then
 open(13,file='./output/pdf.dat',form='formatted',status='new')
 do i=1,nset-1
  write(13,'(2(es16.8))') axis(i)/re,dble(pdf(i))
 enddo
 write(13,'(2(es16.8))') axis(nset)/re,0.0d0
 close(13,status='keep')
  deallocate(pdf)
  deallocate(axis)
endif

return
end
