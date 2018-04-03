subroutine shrink_domain(full,dealiased)

use commondata
use par_size
use shrink_grid

double precision :: full(spx,nz,spy,2), dealiased(dimc(1),dimc(2),dimc(3),2)

integer :: my,mz

my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
mz=floor(2.0d0/3.0d0*dble(nz))

if((sx1.ne.-1).and.(sy1.ne.-1))then
  if(up.eq.0) then
    dealiased(sx1:sx2,1:mz,sy1:sy2,:)=full(sx1:sx2,1:mz,sy1:sy2,:)
  elseif(up.eq.1)then
    dealiased(sx1:sx2,1:mz,1:sy2-sy1+1,:)=full(sx1:sx2,1:mz,sy1:sy2,:)
  endif
  if(sy3.ne.-1)then
    dealiased(sx1:sx2,1:mz,sy2+1:my,:)=full(sx1:sx2,1:mz,sy3:sy4,:)
  endif
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shrink_domain_fg

use commondata
use dual_grid
use shrink_grid



return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shrink_mapping

use commondata
use par_size
use shrink_grid
use mpi
use mpiIO

integer :: mx,my,myl,myu
integer :: starty
integer :: color

mx=floor(2.0d0/3.0d0*dble(nx/2+1))
my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
myl=floor(2.0d0/3.0d0*dble(ny/2+1))
myu=ny-floor(2.0d0/3.0d0*dble(ny/2))+1

sx1=-1
sx2=-1
sy1=-1
sy2=-1
sy3=-1
sy4=-1

if(cstart(1)+1.le.mx) then
  sx1=1
  if(cstart(1)+spx.le.mx) then
    sx2=spx
  else
    sx2=mx-cstart(1)
  endif
endif

if(cstart(3)+1.le.myl)then
  starty=cstart(3)
  sy1=1
  up=0
  if(cstart(3)+spy.le.myl) then
    sy2=spy
  else
    sy2=myl-cstart(3)
    if(cstart(3)+spy.ge.myu)then
      sy3=myu-cstart(3)
      sy4=spy
    endif
  endif
elseif((cstart(3)+1.gt.myl).and.(cstart(3)+1.lt.myu)) then
  if(cstart(3)+spy.ge.myu) then
    starty=myl
    sy1=myu-cstart(3)
    sy2=spy
    up=1
  endif
elseif(cstart(3)+1.ge.myu) then
  starty=myl+(cstart(3)-myu)
  sy1=1
  sy2=spy
  up=1
endif

dimc=0

if(sx1.ne.-1) then
  dimc(1)=sx2-sx1+1
endif

dimc(2)=floor(2.0d0/3.0d0*dble(nz))

if(sy1.ne.-1) then
  dimc(3)=sy2-sy1+1
  if(sy3.ne.-1) then
    dimc(3)=dimc(3)+sy4-sy3+1
  endif
endif

write(*,*) rank,dimc,cstart(1),0,starty,mx,dimc(2),my

color=MPI_UNDEFINED
if((sx1.ne.-1).and.(sy1.ne.-1)) then
  call mpi_type_create_subarray(4,[mx,dimc(2),my,2],[dimc(1),dimc(2),dimc(3),2], &
   &     [cstart(1),0,starty,0],mpi_order_fortran, &
   &     mpi_double_precision,stype,ierr)

  call mpi_type_commit(stype,ierr)

  color=1
endif

call mpi_comm_split(mpi_comm_world,color,0,sp_save_comm,ierr)


! old method, saves variable with aliased modes, now fiel size reduced to 30%
! call mpi_type_create_subarray(4,[nx/2+1,nz,ny,2],[spx,nz,spy,2], &
!  &     [cstart(1),cstart(2),cstart(3),0],mpi_order_fortran, &
!  &     mpi_double_precision,stype,ierr)

! call mpi_type_commit(stype,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shrink_mapping_fg

use commondata
use dual_grid
use shrink_grid
use mpi
use mpiIO


return
end
