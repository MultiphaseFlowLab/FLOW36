subroutine dz(u,duz)
! be careful, do not use same input and output in the call,
! e.g. call dz(u,u) <---- NO!

use commondata
use par_size

double precision :: u(spx,nz,spy,2),duz(spx,nz,spy,2)

integer :: i,j,k,m

! one kernel: each (i,j,m) runs the recursion over k (same operations and order as the
! array-syntax version, which ran the k loop on the host with one small kernel per k)
!$acc parallel loop collapse(3)
do m=1,2
  do j=1,spy
    do i=1,spx
      duz(i,nz,j,m)=0.0d0
      duz(i,nz-1,j,m)= 2.0d0*(nz-1)*u(i,nz,j,m)
      do k=nz-2,1,-1
        duz(i,k,j,m)=duz(i,k+2,j,m)+2.0d0*dble(k)*u(i,k+1,j,m)
      enddo
    enddo
  enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dz_red(u,duz)

use commondata
use par_size

double precision :: u(spx,nz,spy),duz(spx,nz,spy)

integer :: i,j,k

! one kernel: each (i,j) runs the recursion over k (see dz)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    duz(i,nz,j)=0.0d0
    duz(i,nz-1,j)= 2.0d0*(nz-1)*u(i,nz,j)
    do k=nz-2,1,-1
      duz(i,k,j)=duz(i,k+2,j)+2.0d0*dble(k)*u(i,k+1,j)
    enddo
  enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dz_fg(u,duz)
! be careful, do not use same input and output in the call,
! e.g. call dz(u,u) <---- NO!

use dual_grid

double precision :: u(spxpsi,npsiz,spypsi,2),duz(spxpsi,npsiz,spypsi,2)

integer :: i,j,k,m

! one kernel: each (i,j,m) runs the recursion over k (see dz)
!$acc parallel loop collapse(3)
do m=1,2
  do j=1,spypsi
    do i=1,spxpsi
      duz(i,npsiz,j,m)=0.0d0
      duz(i,npsiz-1,j,m)= 2.0d0*(npsiz-1)*u(i,npsiz,j,m)
      do k=npsiz-2,1,-1
        duz(i,k,j,m)=duz(i,k+2,j,m)+2.0d0*dble(k)*u(i,k+1,j,m)
      enddo
    enddo
  enddo
enddo

return
end
