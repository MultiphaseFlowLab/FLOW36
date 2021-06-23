subroutine topology_par(nstep)

use commondata
use wavenumber

double precision, dimension(nxf,nzf,nyf) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,def,rot
double precision, dimension(nxf/2+1,nzf,nyf,2) :: tmpc,tmpc2
double precision :: dv(nzf),norm(3)
double precision :: dx,dy,delta,int_thick
double precision, allocatable :: pdf(:,:),pdf_norm(:,:),axis(:)

integer :: nstep
integer :: i,j,k,l
integer :: nbins,found

character(len=8) :: step

! tanh(3/sqrt(2))=0.97
!int_thick=tanh(3.0d0/dsqrt(2.0d0))
int_thick=0.5d0

! calculate flow topology parameter
! Q=-1 purely rotational flow, Q=0 pure shear flow, Q=+1 purely elongational flow
! u derivatives
do j=1,nyf
  do i=1,nxf/2+1
    tmpc(i,:,j,1)=-kx(i)*uc(i,:,j,2)
    tmpc(i,:,j,2)=+kx(i)*uc(i,:,j,1)
    tmpc2(i,:,j,1)=-ky(j)*uc(i,:,j,2)
    tmpc2(i,:,j,2)=+ky(j)*uc(i,:,j,1)
  enddo
enddo

call spectral_to_phys_fg(tmpc,a11,0)
call spectral_to_phys_fg(tmpc2,a12,0)
call dz(uc,tmpc)
call spectral_to_phys_fg(tmpc,a13,0)

! v derivatives
do j=1,nyf
  do i=1,nxf/2+1
    tmpc(i,:,j,1)=-kx(i)*vc(i,:,j,2)
    tmpc(i,:,j,2)=+kx(i)*vc(i,:,j,1)
    tmpc2(i,:,j,1)=-ky(j)*vc(i,:,j,2)
    tmpc2(i,:,j,2)=+ky(j)*vc(i,:,j,1)
  enddo
enddo

call spectral_to_phys_fg(tmpc,a21,0)
call spectral_to_phys_fg(tmpc2,a22,0)
call dz(vc,tmpc)
call spectral_to_phys_fg(tmpc,a23,0)

! w derivatives
do j=1,nyf
  do i=1,nxf/2+1
    tmpc(i,:,j,1)=-kx(i)*wc(i,:,j,2)
    tmpc(i,:,j,2)=+kx(i)*wc(i,:,j,1)
    tmpc2(i,:,j,1)=-ky(j)*wc(i,:,j,2)
    tmpc2(i,:,j,2)=+ky(j)*wc(i,:,j,1)
  enddo
enddo

call spectral_to_phys_fg(tmpc,a31,0)
call spectral_to_phys_fg(tmpc2,a32,0)
call dz(wc,tmpc)
call spectral_to_phys_fg(tmpc,a33,0)

def=(a11)**2+(a22)**2+(a33)**2+0.5d0*((a13+a31)**2+(a12+a21)**2+(a23+a32)**2)
rot=0.5d0*((a32-a23)**2+(a13-a31)**2+(a21-a12)**2)

Qtop=(def-rot)/(def+rot)
! write(*,*) maxval(Qtop),minval(Qtop)

! get PDF of flow topology parameter

! find volume of each cell for PDF normalization
dx=xfg(nxf)/dble(nxf-1)
dy=yfg(nyf)/dble(nyf-1)

dv(1)=dx*dy*(zfg(1)-zfg(2))*0.5d0
dv(nzf)=dx*dy*(zfg(nzf-1)-zfg(nzf))*0.5d0
do k=2,nzf-1
  dv(k)=dx*dy*(zfg(k-1)-zfg(k+1))*0.5d0
enddo

! compute PDF, split interface, droplet and carrier contributions
nbins=201
! axis, droplet, carrier, interface
allocate(pdf(nbins,4))
allocate(pdf_norm(nbins,3))
allocate(axis(nbins+1))

delta=(1.0d0-(-1.0d0))/dble(nbins-1)
pdf=0.0d0
pdf(1,1)=-1.0d0
do k=2,nbins
  pdf(k,1)=pdf(k-1,1)+delta
enddo

axis(1)=-1.0d0
axis(nbins+1)=1.0d0
do k=2,nbins
  axis(k)=(pdf(k-1,1)+pdf(k,1))*0.5d0
enddo

do j=1,nyf
 do k=1,nzf
  do i=1,nxf
   found=0
   l=1
   do while((found.eq.0).and.(l.le.nbins))
    if((Qtop(i,k,j).ge.axis(l)).and.(Qtop(i,k,j).lt.axis(l+1)))then
     found=1
    else
     l=l+1
    endif
   enddo
    ! control on z to remve effects of spotting at wall
   if(found.eq.1)then
    if(phi(i,k,j).ge.0.0d0.and.(z(k).gt.10.0d0).and.(z(k).lt.re-10.0d0))then
     pdf(l,2)=pdf(l,2)+dv(k)
    elseif(phi(i,k,j).lt.0.0d0)then
     pdf(l,3)=pdf(l,3)+dv(k)
    endif
    if((dabs(phi(i,k,j)).le.int_thick).and.(z(k).gt.10.0d0).and.(z(k).lt.re-10.0d0)) pdf(l,4)=pdf(l,4)+dv(k)
   endif
  enddo
 enddo
enddo

! PDF normalization
norm=0.d0
do l=1,nbins
 norm(1)=norm(1)+(axis(l+1)-axis(l))*pdf(l,2)
 norm(2)=norm(2)+(axis(l+1)-axis(l))*pdf(l,3)
 norm(3)=norm(3)+(axis(l+1)-axis(l))*pdf(l,4)
enddo

pdf_norm(:,1)=pdf(:,2)/norm(1)
pdf_norm(:,2)=pdf(:,3)/norm(2)
pdf_norm(:,3)=pdf(:,4)/norm(3)

! write PDF to file
write(step,'(i8.8)') nstep
open(1234,file='./output/PDF_'//step//'.dat',form='formatted',status='new')
write(1234,'(7(a16,2x))') 'axis','Q drops','Q carrier','Q interface','Q drops norm','Q carrier norm','Q interface norm'
do l=1,nbins
 write(1234,'(7(es16.5,2x))') pdf(l,:),pdf_norm(l,:)
enddo
close(1234,status='keep')


deallocate(pdf,pdf_norm,axis)

return
end
