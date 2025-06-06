subroutine define_dual_grid(cg_cxsize,cg_cysize)

use commondata
use dual_grid
use par_size
use mpi
use mpiIO

integer :: cg_cysize(nzcpu),cg_cxsize(nycpu)
integer :: fg_cysize(nzcpu),fg_cxsize(nycpu)
integer :: fg_pysize(nycpu),fg_pzsize(nzcpu)
integer :: rx,ry,rz
integer :: i,j

npsix=exp_x*nx
npsiy=exp_y*ny
npsiz=exp_z*(nz-1)+1

if((rank.eq.0).and.((exp_x.ne.1).or.(exp_y.ne.1).or.(exp_z.ne.1))) then
 write(*,'(18x,3(3x,a2,2x))') 'nx','ny','nz'
 write(*,'(1x,a,4x,3(i5,2x))') 'Coarse grid: ',nx,ny,nz
 write(*,'(1x,a,4x,3(i5,2x))') 'Fine grid:   ',npsix,npsiy,npsiz
 write(*,*) '-----------------------------------------------------------------------'
endif


! fine grid size in physical space
ry=mod(npsiy,nycpu)
rz=mod(npsiz,nzcpu)

fpypsi=int((npsiy-ry)/nycpu)
if(mod(rank,nycpu).lt.ry)then
 fpypsi=int((npsiy-ry)/nycpu)+1
endif

fpzpsi=int((npsiz-rz)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.rz)then
 fpzpsi=int((npsiz-rz)/nzcpu)+1
endif


fg_pysize=int((npsiy-ry)/nycpu)
do i=1,nycpu
  if(mod(i-1,nycpu).lt.ry) fg_pysize(i)=fg_pysize(i)+1
enddo

fg_pzsize=int((npsiz-rz)/nzcpu)
do j=1,nzcpu
  if(floor(real(nycpu*(j-1))/real(nycpu)).lt.rz) fg_pzsize(j)=fg_pzsize(j)+1
enddo

fstartpsi=[0,0,0]

do i=1,mod(rank,nycpu)
  fstartpsi(3)=fstartpsi(3)+fg_pysize(i)
enddo
do i=1,floor(real(rank)/real(nycpu))
  fstartpsi(2)=fstartpsi(2)+fg_pzsize(i)
enddo


! coarse grid: rank, cxsize, cysize
do j=1,nzcpu
  do i=1,nycpu
    cg_size(i,j,1)=nycpu*(j-1)+(i-1)
    cg_size(i,j,2)=cg_cxsize(i)
    cg_size(i,j,3)=cg_cysize(j)
  enddo
enddo


rx=mod(npsix/2+1,nycpu)
fg_cxsize=int((npsix/2+1-rx)/nycpu)
do i=1,nycpu
  if(mod(i-1,nycpu).lt.rx) fg_cxsize(i)=fg_cxsize(i)+1
enddo


ry=mod(npsiy,nzcpu)
fg_cysize=int((npsiy-ry)/nzcpu)
do j=1,nzcpu
  if(floor(real(nycpu*(j-1))/real(nycpu)).lt.ry) fg_cysize(j)=fg_cysize(j)+1
enddo

! fine grid: rank, cxsize, cysize
do j=1,nzcpu
  do i=1,nycpu
    fg_size(i,j,1)=nycpu*(j-1)+(i-1)
    fg_size(i,j,2)=fg_cxsize(i)
    fg_size(i,j,3)=fg_cysize(j)
  enddo
enddo

cg_size(1,:,4)=1
cg_size(:,1,5)=1
fg_size(1,:,4)=1
fg_size(:,1,5)=1
do i=2,nycpu
  cg_size(i,:,4)=cg_size(i-1,:,4)+cg_size(i-1,:,2)
  fg_size(i,:,4)=fg_size(i-1,:,4)+fg_size(i-1,:,2)
enddo
do j=2,nzcpu
  cg_size(:,j,5)=cg_size(:,j-1,5)+cg_size(:,j-1,3)
  fg_size(:,j,5)=fg_size(:,j-1,5)+fg_size(:,j-1,3)
enddo

! cg_size: (rank,cxsize,cysize,cxstart,cystart) for each MPI task on coarse grid
! fg_size: (rank,cxsize,cysize,cxstart,cystart) for each MPI task on fine grid

! if(rank.eq.0) then
!   write(*,*)
!   write(*,'(a6,2x,a12,4x,a12)') 'Rank','Coarse grid','Fine grid'
!   do j=1,nzcpu
!     do i=1,nycpu
!       ! dimension (local per rank)
!       ! write(*,'(i6,2x,2(i6),4x,2(i6))') cg_size(i,j,1),cg_size(i,j,2),cg_size(i,j,3),fg_size(i,j,2),fg_size(i,j,3)
!       ! starting indexes (local per rank)
!       write(*,'(i6,2x,2(i6),4x,2(i6))') cg_size(i,j,1),cg_size(i,j,4),cg_size(i,j,5),fg_size(i,j,4),fg_size(i,j,5)
!     enddo
!   enddo
! endif


! allocate correct sizes for each rank
! find indexes corresponding to rank
i=mod(rank,nycpu)+1
j=floor(real(rank)/real(nycpu))+1
spxpsi=fg_size(i,j,2)
spypsi=fg_size(i,j,3)
cstartpsi(1)=fg_size(i,j,4)-1
cstartpsi(2)=0
cstartpsi(3)=fg_size(i,j,5)-1
! write(*,*) rank,cg_size(i,j,1),i,j
! write(*,*) rank,cstartpsi(1),spxpsi,cstartpsi(3),spypsi


! create datatype for MPI I/O operations on fine grid
! physical space
call mpi_type_create_subarray(3,[npsix,npsiz,npsiy],[npsix,fpzpsi,fpypsi], &
&     fstartpsi,mpi_order_fortran, &
&     mpi_double_precision,ftype_fg,ierr)

call mpi_type_commit(ftype_fg,ierr)

! spectral space, save only dealiased grid
! mx=floor(2.0d0/3.0d0*dble(npsix/2+1))
! my=floor(2.0d0/3.0d0*dble(npsiy/2+1))+floor(2.0d0/3.0d0*dble(npsiy/2))
!
! call mpi_type_create_subarray(4,[npsix/2+1,npsiz,npsiy,2],[spxpsi,npsiz,spypsi,2], &
!  &     [cstartpsi(1),cstartpsi(2),cstartpsi(3),0],mpi_order_fortran, &
!  &     mpi_double_precision,stype_fg,ierr)
!
! call mpi_type_commit(stype_fg,ierr)

return
end
