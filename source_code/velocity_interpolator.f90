subroutine lagran4(ppos,pvel)

use commondata
use grid
use particle
use sim_par

integer, dimension(1) :: pos_ind
integer, dimension(5) :: iarr,jarr,karr
integer :: i0,i1,i2,i3,i4,j0,j1,j2,j3,j4,k0,k1,k2,k3,k4
integer :: i,j,k

double precision, dimension(3) :: ppos,pvel
double precision, dimension(3) :: x0,x1,x2,x3,x4
double precision, dimension(5,3) :: pol

#define NX nnx

! particles in + units, flow in - units
! particles: z in [0,2*Re], flow in [-1,+1]
pos_ind=minloc(abs(x*re-ppos(1)))
i2=pos_ind(1)
pos_ind=minloc(abs(y*re-ppos(2)))
j2=pos_ind(1)
pos_ind=minloc(abs((z+1.0d0)*re-ppos(3)))
k2=pos_ind(1)

! initialize other indexes, check also for periodicity
! x indexes (periodic)
i0=i2-2
i1=i2-1
i3=i2+1
i4=i2+2
do while(i0.lt.1)
 i0=i0+nx
enddo
do while(i1.lt.1)
 i1=i1+nx
enddo
do while(i3.gt.nx)
 i3=i3-nx
enddo
do while(i4.gt.nx)
 i4=i4-nx
enddo

! y indexes (periodic)
j0=j2-2
j1=j2-1
j3=j2+1
j4=j2+2
do while(j0.lt.1)
 j0=j0+ny
enddo
do while(j1.lt.1)
 j1=j1+ny
enddo
do while(j3.gt.ny)
 j3=j3-ny
enddo
do while(j4.gt.ny)
 j4=j4-ny
enddo

! z indexes
k0=k2-2
k1=k2-1
k3=k2+1
k4=k2+2
do while(k0.lt.1)
 k0=k0+1
 k1=k1+1
 k2=k2+1
 k3=k3+1
 k4=k4+1
enddo
do while(k4.gt.nz)
 k0=k0-1
 k1=k1-1
 k2=k2-1
 k3=k3-1
 k4=k4-1
enddo

iarr=[i0,i1,i2,i3,i4]
jarr=[j0,j1,j2,j3,j4]
karr=[k0,k1,k2,k3,k4]

! points for 4th order Lagrangian interpolation (wall units)
x0=[x(i0),y(j0),(z(k0)+1.0d0)]*re
x1=[x(i1),y(j1),(z(k1)+1.0d0)]*re
x2=[x(i2),y(j2),(z(k2)+1.0d0)]*re
x3=[x(i3),y(j3),(z(k3)+1.0d0)]*re
x4=[x(i4),y(j4),(z(k4)+1.0d0)]*re

! Lagrangian polynomials
! 0th polynomial
pol(1,1)=(ppos(1)-x1(1))/(x0(1)-x1(1))* &
 &       (ppos(1)-x2(1))/(x0(1)-x2(1))* &
 &       (ppos(1)-x3(1))/(x0(1)-x3(1))* &
 &       (ppos(1)-x4(1))/(x0(1)-x4(1))
pol(1,2)=(ppos(2)-x1(2))/(x0(2)-x1(2))* &
 &       (ppos(2)-x2(2))/(x0(2)-x2(2))* &
 &       (ppos(2)-x3(2))/(x0(2)-x3(2))* &
 &       (ppos(2)-x4(2))/(x0(2)-x4(2))
pol(1,3)=(ppos(3)-x1(3))/(x0(3)-x1(3))* &
 &       (ppos(3)-x2(3))/(x0(3)-x2(3))* &
 &       (ppos(3)-x3(3))/(x0(3)-x3(3))* &
 &       (ppos(3)-x4(3))/(x0(3)-x4(3))
! 1st polynomial
pol(2,1)=(ppos(1)-x0(1))/(x1(1)-x0(1))* &
 &       (ppos(1)-x2(1))/(x1(1)-x2(1))* &
 &       (ppos(1)-x3(1))/(x1(1)-x3(1))* &
 &       (ppos(1)-x4(1))/(x1(1)-x4(1))
pol(2,2)=(ppos(2)-x0(2))/(x1(2)-x0(2))* &
 &       (ppos(2)-x2(2))/(x1(2)-x2(2))* &
 &       (ppos(2)-x3(2))/(x1(2)-x3(2))* &
 &       (ppos(2)-x4(2))/(x1(2)-x4(2))
pol(2,3)=(ppos(3)-x0(3))/(x1(3)-x0(3))* &
 &       (ppos(3)-x2(3))/(x1(3)-x2(3))* &
 &       (ppos(3)-x3(3))/(x1(3)-x3(3))* &
 &       (ppos(3)-x4(3))/(x1(3)-x4(3))
! 2nd polynomial
pol(3,1)=(ppos(1)-x0(1))/(x2(1)-x0(1))* &
 &       (ppos(1)-x1(1))/(x2(1)-x1(1))* &
 &       (ppos(1)-x3(1))/(x2(1)-x3(1))* &
 &       (ppos(1)-x4(1))/(x2(1)-x4(1))
pol(3,2)=(ppos(2)-x0(2))/(x2(2)-x0(2))* &
 &       (ppos(2)-x1(2))/(x2(2)-x1(2))* &
 &       (ppos(2)-x3(2))/(x2(2)-x3(2))* &
 &       (ppos(2)-x4(2))/(x2(2)-x4(2))
pol(3,3)=(ppos(3)-x0(3))/(x2(3)-x0(3))* &
 &       (ppos(3)-x1(3))/(x2(3)-x1(3))* &
 &       (ppos(3)-x3(3))/(x2(3)-x3(3))* &
 &       (ppos(3)-x4(3))/(x2(3)-x4(3))
! 3rd polynomial
pol(4,1)=(ppos(1)-x0(1))/(x3(1)-x0(1))* &
 &       (ppos(1)-x1(1))/(x3(1)-x1(1))* &
 &       (ppos(1)-x2(1))/(x3(1)-x2(1))* &
 &       (ppos(1)-x4(1))/(x3(1)-x4(1))
pol(4,2)=(ppos(2)-x0(2))/(x3(2)-x0(2))* &
 &       (ppos(2)-x1(2))/(x3(2)-x1(2))* &
 &       (ppos(2)-x2(2))/(x3(2)-x2(2))* &
 &       (ppos(2)-x4(2))/(x3(2)-x4(2))
pol(4,3)=(ppos(3)-x0(3))/(x3(3)-x0(3))* &
 &       (ppos(3)-x1(3))/(x3(3)-x1(3))* &
 &       (ppos(3)-x2(3))/(x3(3)-x2(3))* &
 &       (ppos(3)-x4(3))/(x3(3)-x4(3))
! 4th polynomial
pol(5,1)=(ppos(1)-x0(1))/(x4(1)-x0(1))* &
 &       (ppos(1)-x1(1))/(x4(1)-x1(1))* &
 &       (ppos(1)-x2(1))/(x4(1)-x2(1))* &
 &       (ppos(1)-x3(1))/(x4(1)-x3(1))
pol(5,2)=(ppos(2)-x0(2))/(x4(2)-x0(2))* &
 &       (ppos(2)-x1(2))/(x4(2)-x1(2))* &
 &       (ppos(2)-x2(2))/(x4(2)-x2(2))* &
 &       (ppos(2)-x3(2))/(x4(2)-x3(2))
pol(5,3)=(ppos(3)-x0(3))/(x4(3)-x0(3))* &
 &       (ppos(3)-x1(3))/(x4(3)-x1(3))* &
 &       (ppos(3)-x2(3))/(x4(3)-x2(3))* &
 &       (ppos(3)-x3(3))/(x4(3)-x3(3))

pvel=0.0d0

#if NX > 2
do j=1,5
 do k=1,5
  do i=1,5
   pvel(1)=pvel(1)+uf(iarr(i),karr(k),jarr(j))*pol(i,1)*pol(j,2)*pol(k,3)
   pvel(2)=pvel(2)+vf(iarr(i),karr(k),jarr(j))*pol(i,1)*pol(j,2)*pol(k,3)
   pvel(3)=pvel(3)+wf(iarr(i),karr(k),jarr(j))*pol(i,1)*pol(j,2)*pol(k,3)
  enddo
 enddo
enddo
#else
! for 2D simulation no interpolation along x
do j=1,5
 do k=1,5
  pvel(1)=pvel(1)+uf(i2,karr(k),jarr(j))*pol(j,2)*pol(k,3)
  pvel(2)=pvel(2)+vf(i2,karr(k),jarr(j))*pol(j,2)*pol(k,3)
  pvel(3)=pvel(3)+wf(i2,karr(k),jarr(j))*pol(j,2)*pol(k,3)
 enddo
enddo
#endif

return
end subroutine
