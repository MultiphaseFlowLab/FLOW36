subroutine calc_curvature(nstep)

use commondata
! use wavenumbers
use pdf_calc

double precision :: gradm,mean,rms,curv
double precision, dimension(nx,nz,ny) :: phix,phiy,phiz

integer :: i,j,k
integer :: count

character(len=40) :: namefile
character(len=8) :: time


!  calculate gradient with central finite differences, better results than with spectral methods
phix=0.0d0
phiy=0.0d0
phiz=0.0d0
! x derivative
phix(1,:,:)=(phi(2,:,:)-phi(1,:,:))/(x(2)-x(1))
phix(nx,:,:)=(phi(nx,:,:)-phi(nx-1,:,:))/(x(nx)-x(nx-1))
do i=2,nx-1
  phix(i,:,:)=(phi(i+1,:,:)-phi(i-1,:,:))/(x(i+1)-x(i-1))
enddo
! y derivative
phiy(:,:,1)=(phi(:,:,2)-phi(:,:,1))/(y(2)-y(1))
phiy(:,:,ny)=(phi(:,:,ny)-phi(:,:,ny-1))/(y(ny)-z(ny-1))
do j=2,ny-1
  phiy(:,:,j)=(phi(:,:,j+1)-phi(:,:,j-1))/(y(j+1)-y(j-1))
enddo
! z derivative
phiz(:,1,:)=(phi(:,2,:)-phi(:,1,:))/(z(1)-z(2))
phiz(:,nz,:)=(phi(:,nz,:)-phi(:,nz-1,:))/(z(nz-1)-z(nz))
do k=2,nz-1
  phiz(:,k,:)=(phi(:,k+1,:)-phi(:,k-1,:))/(z(k-1)-z(k+1))
enddo

! normalize ||gradient(phi)||=1
do j=1,ny
  do k=1,nz
    do i=1,nx
      gradm=((phix(i,k,j))**2+(phiy(i,k,j))**2+(phiz(i,k,j))**2)**0.5
      ! gradm contains the modulus of the gradient of phi
      ! normalized gradient of phi
      phix(i,k,j)=phix(i,k,j)/gradm
      phiy(i,k,j)=phiy(i,k,j)/gradm
      phiz(i,k,j)=phiz(i,k,j)/gradm
      ! avoid overshoots when grad(phi)~0
      ! if(abs(phi(i,k,j)).gt.0.8)then
      !   phix(i,k,j)=0.0d0
      !   phiy(i,k,j)=0.0d0
      !   phiz(i,k,j)=0.0d0
      ! endif
    enddo
  enddo
enddo


! calculate divergence with central finite differences
! add x direction derivatives
kv=0.0d0
kv(1,:,:)=kv(1,:,:)+(phix(2,:,:)-phix(1,:,:))/(x(2)-x(1))
kv(nx,:,:)=kv(nx,:,:)+(phix(nx,:,:)-phix(nx-1,:,:))/(x(nx)-x(nx-1))
do i=2,nx-1
  kv(i,:,:)=kv(i,:,:)+(phix(i+1,:,:)-phix(i-1,:,:))/(x(i+1)-x(i-1))
enddo
! add y direction derivatives
kv(:,:,1)=kv(:,:,1)+(phiy(:,:,2)-phiy(:,:,1))/(y(2)-y(1))
kv(:,:,ny)=kv(:,:,ny)+(phiy(:,:,ny)-phiy(:,:,ny-1))/(y(ny)-y(ny-1))
do j=2,ny-1
  kv(:,:,j)=kv(:,:,j)+(phiy(:,:,j+1)-phiy(:,:,j-1))/(y(j+1)-y(j-1))
enddo
! add z direction derivatives
kv(:,1,:)=kv(:,1,:)+(phiz(:,2,:)-phiz(:,1,:))/(z(1)-z(2))
kv(:,nz,:)=kv(:,nz,:)+(phiz(:,nz,:)-phiz(:,nz-1,:))/(z(nz-1)-z(nz))
do k=2,nz-1
  kv(:,k,:)=kv(:,k,:)+(phiz(:,k+1,:)-phiz(:,k-1,:))/(z(k-1)-z(k+1))
enddo


! write one file per timestep with all curvature values at the interface
write(time,'(i8.8)') nstep
namefile='./output/k_'//time//'.dat'
open(75,file=trim(namefile),form='formatted',status='new')
count=0
mean=0.0d0
do j=1,ny
  do k=1,nz
    do i=1,nx
      if(abs(phi(i,k,j)).le.threshold)then
      ! if phi was rescaled in read_fields for improved stability
      !if(abs(phi(i,k,j)).le.threshold/50.0d0)then
        count=count+1
        write(75,'(es16.8)') kv(i,k,j)
        ! calculate mean over time
        mean=mean+kv(i,k,j)
      endif
    enddo
  enddo
enddo
close(75,status='keep')

mean=mean/dble(count)

! calculate rms over time
rms=0.0d0
open(75,file=trim(namefile),form='formatted',status='old')
do i=1,count
  read(75,'(es16.8)') curv
  rms=rms+(curv-mean)**2
enddo
close(75,status='keep')
! get rms
rms=(rms/dble(count))**0.5

write(fhandle,'(i16,2x,3(es16.8,2x))') nstep,dble(nstep)*Re*dt,mean,rms
call flush(fhandle)

write(*,*) 'Number of samples: ',count

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
