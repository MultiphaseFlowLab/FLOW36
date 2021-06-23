subroutine get_interface(nstep)

use commondata

integer :: nstep
integer :: i,j,k,id,jd,kd
integer :: top(nx,nz,ny),s_drop(nx,nz,ny)
integer :: drop_count

character(len=8) :: time

write(time,'(i8.8)') nstep
open(14,file='./output/psiV_'//time//'.dat',status='new',form='formatted')
write(14,'(3(a16))') 'Drop','V (+ units)','psi'
close(14,status='keep')


do j=1,ny
  do k=1,nz
    do i=1,nx
      if(phi(i,k,j).ge.0.0d0)then
        top(i,k,j)=1
      else
        top(i,k,j)=0
      endif
    enddo
  enddo
enddo

drop_count=0

! flood fill algorithm
do jd=1,ny
 do kd=1,nz
  do id=1,nx
   if(top(id,kd,jd).gt.0)then
    drop_count=drop_count+1
    write(*,'(2x,a,i3,a)') 'New drop, ',drop_count,' drops'
    ! single drop part
    s_drop=0
    s_drop(id,kd,jd)=1
    ! flood fill algorithm
    call flood_fill(top,s_drop,id,jd,kd)
    ! remove drops already done from top
    top=top-s_drop
    ! calculate drop volume and center of mass
    call get_center(s_drop,nstep,drop_count)
    ! new drop calculation
   endif
  enddo
 enddo
enddo


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_center(s_drop,nstep,drop_count)

use commondata

double precision :: sum_xy(nz)
double precision :: int_1,dx,dy,dz(nz),Vol,psival

integer :: i,j,k
integer :: s_drop(nx,nz,ny)
integer :: nstep,drop_count

logical :: answer

character(len=8) :: time

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)
! dx, dy are uniform in space, dz=dz(z), integral of phi=+1 on V can be reduced to
! a sum of integrals over z
sum_xy=sum(sum(dble(s_drop),3),1)
! sum_xy=0.0d0
! do j=1,ny
!   do k=1,nz
!     do i=1,nx
!       if(s_drop(i,k,j).eq.1) sum_xy(k)=sum_xy(k)+phi(i,k,j)
!     enddo
!   enddo
! enddo

dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do k=2,nz-1
  dz(k)=(z(k-1)-z(k+1))*0.5d0
enddo

int_1=0.0d0
do k=1,nz
  int_1=int_1+sum_xy(k)*dz(k)
enddo
int_1=int_1*dx*dy

Vol=int_1*Re**3 ! volume in + units

write(time,'(i8.8)') nstep
open(14,file='./output/psiV_'//time//'.dat',status='old',form='formatted',access='append')
! do loop to find interfacial values of surfactant concentration
do j=1,ny
  do k=1,nz
    do i=1,nx
      if(s_drop(i,k,j).eq.1)then
        ! check if it is interface
        call check_interface(answer,i,j,k,s_drop)
        if(answer.eqv..true.)then
        call get_psi(i,k,j,psival)
        write(14,'(i16,2(es16.6))') drop_count,Vol,psival
        endif
      endif
    enddo
  enddo
enddo

close(14,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_interface(answer,i,j,k,s_drop)

use commondata

logical :: answer,boundary

integer :: s_drop(nx,nz,ny)

integer :: i,j,k
integer :: ii,jj

boundary=.false.

! could be done with elseif
if(s_drop(i,k,j).eq.1)then
 ! check orthogonal direction
 !xdir
 ii=mod(i,nx)+1
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
 ii=mod(i+nx-2,nx)+1
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
! ydir
 jj=mod(j,ny)+1
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
 jj=mod(j+ny-2,ny)+1
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
! zdir
 if(s_drop(i,k+1,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k-1,j).eq.0)then
  boundary=.true.
 endif
! ! check on edges (all 3x3x3 cube)
!  ! i+1,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ii=mod(i,nx)+1
!  jj=mod(j,ny)+1
!  if(s_drop(ii,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i+1,j-1
!  jj=mod(j+ny-2,ny)+1
!  if(s_drop(ii,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i+1,j
!  if(s_drop(ii,k-1,j).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,j).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,j).eq.0)then
!   boundary=.true.
!  endif
!  ! i,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  jj=mod(j,ny)+1
!  if(s_drop(i,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(i,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(i,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i,j-1
!  jj=mod(j+ny-2,ny)+1
!  if(s_drop(i,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(i,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(i,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i,j
!  if(s_drop(i,k-1,j).eq.0)then
!   boundary=.true.
!  endif
!  ! if(s_drop(i,k,j).eq.0)then
!  !  boundary=.true.
!  ! endif
!  if(s_drop(i,k+1,j).eq.0)then
!   boundary=.true.
!  endif
!  ! i-1,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ii=mod(i+nx-2,nx)+1
!  jj=mod(j,ny)+1
!  if(s_drop(ii,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i-1,j-1
!  jj=mod(j+ny-2,ny)+1
!  if(s_drop(ii,k-1,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,jj).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,jj).eq.0)then
!   boundary=.true.
!  endif
!  ! i-1,j
!  if(s_drop(ii,k-1,j).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k,j).eq.0)then
!   boundary=.true.
!  endif
!  if(s_drop(ii,k+1,j).eq.0)then
!   boundary=.true.
!  endif

endif


! if at least one neighbour is equal to zero and the node value is 1 it is an interface node
if((s_drop(i,k,j).eq.1).and.(boundary.eqv..true.))then
 answer=.true.
else
 answer=.false.
endif


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_psi(i,k,j,psival)

use commondata

integer :: i,j,k, ii,jj,kk
double precision :: psival

ii=minloc(abs(xfg-x(i)),1)
jj=minloc(abs(yfg-y(j)),1)
kk=minloc(abs(zfg-z(k)),1)

psival=psi(ii,kk,jj)

return
end
