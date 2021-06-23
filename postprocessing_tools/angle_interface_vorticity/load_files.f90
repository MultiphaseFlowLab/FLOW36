subroutine load_files(nstep)

use mpi
use commondata
use vars
use sim_parameter
use wavenumbers

double precision :: md
double precision, dimension(nx,nz,ny) :: phix,phiy,phiz,cosine
double precision, dimension(nx/2+1,nz,ny,2) :: tmp_c
double precision, dimension(nx,nz,ny,3) :: vor
double precision :: xp,yp,zp,cosp

integer :: nstep
integer :: i,j,k,c_int
integer :: ip,jp,kp

character(len=8) :: time

logical :: check


write(time,'(i8.8)') nstep

! read fields in physical space
inquire(file=trim(folder)//'/u_'//time//'.dat',exist=check)
if(check.eqv..false.)then
  if(rank.eq.0)write(*,*) 'Missing u input file n=',time,' , stopping program'
  call mpi_barrier(mpi_comm_world,ierr)
  call exit(0)
endif
call read_fields(u,nstep,'u    ')
inquire(file=trim(folder)//'/v_'//time//'.dat',exist=check)
if(check.eqv..false.)then
  if(rank.eq.0)write(*,*) 'Missing v input file n=',time,' , stopping program'
  call mpi_barrier(mpi_comm_world,ierr)
  call exit(0)
endif
call read_fields(v,nstep,'v    ')
inquire(file=trim(folder)//'/w_'//time//'.dat',exist=check)
if(check.eqv..false.)then
  if(rank.eq.0)write(*,*) 'Missing w input file n=',time,' , stopping program'
  call mpi_barrier(mpi_comm_world,ierr)
  call exit(0)
endif
call read_fields(w,nstep,'w    ')
if(phiflag.eq.1)then
  inquire(file=trim(folder)//'/phi_'//time//'.dat',exist=check)
  if(check.eqv..false.)then
    if(rank.eq.0)write(*,*) 'Missing phi input file n=',time,' , stopping program'
    call mpi_barrier(mpi_comm_world,ierr)
    call exit(0)
  endif
  call read_fields(phi,nstep,'phi  ')
else
  phi=0.0d0
  rhor=1.0d0
endif
! write(*,*) maxval(phi),minval(phi)


!call cm_velocity(phi,u,ucm)
!call cm_velocity(phi,v,vcm)
!call cm_velocity(phi,w,wcm)
!write(*,*) ucm,vcm,wcm

call phys_to_spectral(u,uc,0)
call phys_to_spectral(v,vc,0)
call phys_to_spectral(w,wc,0)

call phys_to_spectral(phi,phic,0)

! calculation of gradient in modal space
! calculate phi gradient on x-y plane
do j=1,ny
  do i=1,nx/2+1
    tmp_c(i,:,j,1)=-kx(i)*phic(i,:,j,2)
    tmp_c(i,:,j,2)=+kx(i)*phic(i,:,j,1)
  enddo
enddo
call spectral_to_phys(tmp_c,phix,0)

do j=1,ny
  do i=1,nx/2+1
    tmp_c(i,:,j,1)=-ky(j)*phic(i,:,j,2)
    tmp_c(i,:,j,2)=+ky(j)*phic(i,:,j,1)
  enddo
enddo
call spectral_to_phys(tmp_c,phiy,0)

call dz(phic,tmp_c)
call spectral_to_phys(tmp_c,phiz,0)



! normalize gradient of phi to get outward pointing unitary normal
do j=1,ny
  do k=1,nz
    do i=1,nx
      md=sqrt(phix(i,k,j)**2+phiy(i,k,j)**2+phiz(i,k,j)**2)
      phix(i,k,j)=-phix(i,k,j)/md
      phiy(i,k,j)=-phiy(i,k,j)/md
      phiz(i,k,j)=-phiz(i,k,j)/md
    enddo
  enddo
enddo


! calculate vorticity
! x component
call dz(vc,tmp_c)

tmp_c=-tmp_c
do j=1,ny
  do i=1,nx/2+1
    tmp_c(i,:,j,1)=tmp_c(i,:,j,1)-ky(j)*wc(i,:,j,2)
    tmp_c(i,:,j,2)=tmp_c(i,:,j,2)+ky(j)*wc(i,:,j,1)
  enddo
enddo

call spectral_to_phys(tmp_c,vor(:,:,:,1),0)

! y component
call dz(uc,tmp_c)

do j=1,ny
  do i=1,nx/2+1
    tmp_c(i,:,j,1)=tmp_c(i,:,j,1)+kx(i)*wc(i,:,j,2)
    tmp_c(i,:,j,2)=tmp_c(i,:,j,2)-kx(i)*wc(i,:,j,1)
  enddo
enddo

call spectral_to_phys(tmp_c,vor(:,:,:,2),0)

! z component
do j=1,ny
  do i=1,nx/2+1
    tmp_c(i,:,j,1)=-kx(i)*vc(i,:,j,2)+ky(j)*uc(i,:,j,2)
    tmp_c(i,:,j,2)=+kx(i)*vc(i,:,j,1)-ky(j)*uc(i,:,j,1)
  enddo
enddo

call spectral_to_phys(tmp_c,vor(:,:,:,3),0)


! normalize vorticity to get unitary modulus
do j=1,ny
  do k=1,nz
    do i=1,nx
      md=((vor(i,k,j,1))**2+(vor(i,k,j,2))**2+(vor(i,k,j,3))**2)**0.5
      vor(i,k,j,:)=vor(i,k,j,:)/md
    enddo
  enddo
enddo


! calculate cosine of angle between interface normal and vorticity
do j=1,ny
  do k=1,nz
    do i=1,nx
      cosine(i,k,j)=phix(i,k,j)*vor(i,k,j,1)+phiy(i,k,j)*vor(i,k,j,2)+phiz(i,k,j)*vor(i,k,j,3)
    enddo
  enddo
enddo


c_int=0

open(444,file='./output/samples_'//time//'.dat',form='formatted',status='new',access='stream',position='rewind')

do j=1,ny
  do k=1,nz
    do i=1,nx
      ip=mod(i,nx)+1
      jp=mod(j,ny)+1
      kp=mod(k,nz)+1
      if(phi(i,k,j)*phi(ip,k,j).le.0.0d0) then   ! x direction intersection
        c_int=c_int+1
        call get_interface(x,nx,xl,i,ip,xp,phi(:,k,j))
        call get_cosine(x,nx,xl,i,ip,xp,cosine(:,k,j),cosp)
        ! write output
        write(444,'(1(e16.6))') cosp
      endif
      if(phi(i,k,j)*phi(i,k,jp).le.0.0d0) then   ! y direction intersection
        c_int=c_int+1
        call get_interface(y,ny,yl,j,jp,yp,phi(i,k,:))
        call get_cosine(y,ny,yl,j,jp,yp,cosine(i,k,:),cosp)
        ! write output
        write(444,'(1(e16.6))') cosp
      endif
      if(phi(i,k,j)*phi(i,kp,j).le.0.0d0) then   ! z direction intersection
        c_int=c_int+1
        call get_interface(z,nz,2.0d0,k,kp,zp,phi(i,:,j))
        call get_cosine(z,nz,2.0d0,k,kp,zp,cosine(i,:,j),cosp)
        ! write output
        write(444,'(1(e16.6))') cosp
      endif
    enddo
  enddo
enddo

close(444,status='keep')

open(444,file='./output/n_samples_'//time//'.dat',form='formatted',status='new')
write(444,'(a20,i20)') 'number of samples:',c_int
close(444,status='keep')



total=total+c_int

write(*,*) 'Interface points: ',c_int

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_gradient(ax,dim,len,ip,ipp,xp,phix,phiy,phiz,nrmx,nrmy,nrmz)

use commondata

integer :: dim
integer :: ip,ipp

double precision, dimension(dim) :: ax
double precision, dimension(dim) :: phix,phiy,phiz
double precision :: nrmx,nrmy,nrmz,len
double precision :: m,xp,mod


if(ipp.lt.ip) then
! apply periodicity
 if(xp.lt.ax(ip)) then
  m=(xp+len-ax(ip))/(ax(ipp)+len-ax(ip))
 else
  m=(xp-ax(ip))/(ax(ipp)+len-ax(ip))
 endif
else
 m=(xp-ax(ip))/(ax(ipp)-ax(ip))
endif

nrmx=phix(ip)+m*(phix(ipp)-phix(ip))
nrmy=phiy(ip)+m*(phiy(ipp)-phiy(ip))
nrmz=phiz(ip)+m*(phiz(ipp)-phiz(ip))

mod=(nrmx**2+nrmy**2+nrmz**2)**0.5
nrmx=nrmx/mod
nrmy=nrmy/mod
nrmz=nrmz/mod

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_cosine(ax,dim,len,ip,ipp,xp,var,valp)

use commondata

integer :: dim
integer :: ip,ipp

double precision, dimension(dim) :: ax
double precision, dimension(dim) :: var
double precision :: valp,len
double precision :: m,xp


if(ipp.lt.ip) then
! apply periodicity
 if(xp.lt.ax(ip)) then
  m=(xp+len-ax(ip))/(ax(ipp)+len-ax(ip))
 else
  m=(xp-ax(ip))/(ax(ipp)+len-ax(ip))
 endif
else
 m=(xp-ax(ip))/(ax(ipp)-ax(ip))
endif

valp=var(ip)+m*(var(ipp)-var(ip))

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cm_velocity(phi,vel,velcm)

use commondata

double precision, dimension(nx,nz,ny) :: phi,vel
double precision :: velcm
double precision :: num,den,dx,dy,dz(nz)

integer :: i,j,k

dx=x(2)-x(1)
dy=y(2)-y(1)
dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do k=2,nz-1
 dz(k)=0.5d0*(z(k-1)-z(k+1))
enddo

num=0.0d0
den=0.0d0

do j=1,ny
 do k=1,nz
  do i=1,nx
   if(phi(i,k,j).gt.0.0d0)then
    num=num+vel(i,k,j)*dx*dy*dz(k)
    den=den+dx*dy*dz(k)
   endif
  enddo
 enddo
enddo

velcm=num/den

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_interface(ax,dim,len,ip,ipp,xpos,var)

use commondata

integer :: dim,ip,ipp

double precision, dimension(dim) :: ax,var
double precision :: xpos,len

! write(*,*) maxval(ax)

if(ipp.lt.ip) then
! apply periodicity
  xpos=ax(ip)-var(ip)/(var(ipp)-var(ip))*(ax(ipp)+len-ax(ip))
else
  xpos=ax(ip)-var(ip)/(var(ipp)-var(ip))*(ax(ipp)-ax(ip))
endif

if(xpos.gt.len) xpos=xpos-len

return
end
