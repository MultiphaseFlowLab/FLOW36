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

subroutine shrink_domain_fg(full,dealiased)

use commondata
use dual_grid
use shrink_grid

double precision :: full(spxpsi,npsiz,spypsi,2), dealiased(dimc_fg(1),dimc_fg(2),dimc_fg(3),2)

integer :: my,mz

my=floor(2.0d0/3.0d0*dble(npsiy/2+1))+floor(2.0d0/3.0d0*dble(npsiy/2))
mz=floor(2.0d0/3.0d0*dble(npsiz))

if((sx1_fg.ne.-1).and.(sy1_fg.ne.-1))then
  if(up_fg.eq.0) then
    dealiased(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)=full(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)
  elseif(up_fg.eq.1)then
    dealiased(sx1_fg:sx2_fg,1:mz,1:sy2_fg-sy1_fg+1,:)=full(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)
  endif
  if(sy3_fg.ne.-1)then
    dealiased(sx1_fg:sx2_fg,1:mz,sy2_fg+1:my,:)=full(sx1_fg:sx2_fg,1:mz,sy3_fg:sy4_fg,:)
  endif
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine expand_domain(dealiased,full)

use commondata
use par_size
use shrink_grid

double precision :: full(spx,nz,spy,2), dealiased(dimc(1),dimc(2),dimc(3),2)

integer :: my,mz

my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
mz=floor(2.0d0/3.0d0*dble(nz))

full=0.0d0
if((sx1.ne.-1).and.(sy1.ne.-1))then
  if(up.eq.0) then
    full(sx1:sx2,1:mz,sy1:sy2,:)=dealiased(sx1:sx2,1:mz,sy1:sy2,:)
  elseif(up.eq.1)then
    full(sx1:sx2,1:mz,sy1:sy2,:)=dealiased(sx1:sx2,1:mz,1:sy2-sy1+1,:)
  endif
  if(sy3.ne.-1)then
    full(sx1:sx2,1:mz,sy3:sy4,:)=dealiased(sx1:sx2,1:mz,sy2+1:my,:)
  endif
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine expand_domain_fg(dealiased,full)

use commondata
use dual_grid
use shrink_grid

double precision :: full(spxpsi,npsiz,spypsi,2), dealiased(dimc_fg(1),dimc_fg(2),dimc_fg(3),2)

integer :: my,mz

my=floor(2.0d0/3.0d0*dble(npsiy/2+1))+floor(2.0d0/3.0d0*dble(npsiy/2))
mz=floor(2.0d0/3.0d0*dble(npsiz))

full=0.0d0
if((sx1_fg.ne.-1).and.(sy1_fg.ne.-1))then
  if(up_fg.eq.0) then
    full(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)=dealiased(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)
  elseif(up_fg.eq.1)then
    full(sx1_fg:sx2_fg,1:mz,sy1_fg:sy2_fg,:)=dealiased(sx1_fg:sx2_fg,1:mz,1:sy2_fg-sy1_fg+1,:)
  endif
  if(sy3_fg.ne.-1)then
    full(sx1_fg:sx2_fg,1:mz,sy3_fg:sy4_fg,:)=dealiased(sx1_fg:sx2_fg,1:mz,sy2_fg+1:my,:)
  endif
endif

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

integer :: mx,my,myl,myu
integer :: starty
integer :: color

mx=floor(2.0d0/3.0d0*dble(npsix/2+1))
my=floor(2.0d0/3.0d0*dble(npsiy/2+1))+floor(2.0d0/3.0d0*dble(npsiy/2))
myl=floor(2.0d0/3.0d0*dble(npsiy/2+1))
myu=npsiy-floor(2.0d0/3.0d0*dble(npsiy/2))+1

sx1_fg=-1
sx2_fg=-1
sy1_fg=-1
sy2_fg=-1
sy3_fg=-1
sy4_fg=-1

if(cstartpsi(1)+1.le.mx) then
  sx1_fg=1
  if(cstartpsi(1)+spxpsi.le.mx) then
    sx2_fg=spxpsi
  else
    sx2_fg=mx-cstartpsi(1)
  endif
endif

if(cstartpsi(3)+1.le.myl)then
  starty=cstartpsi(3)
  sy1_fg=1
  up_fg=0
  if(cstartpsi(3)+spypsi.le.myl) then
    sy2_fg=spypsi
  else
    sy2_fg=myl-cstartpsi(3)
    if(cstartpsi(3)+spypsi.ge.myu)then
      sy3_fg=myu-cstartpsi(3)
      sy4_fg=spypsi
    endif
  endif
elseif((cstartpsi(3)+1.gt.myl).and.(cstartpsi(3)+1.lt.myu)) then
  if(cstartpsi(3)+spypsi.ge.myu) then
    starty=myl
    sy1_fg=myu-cstartpsi(3)
    sy2_fg=spypsi
    up_fg=1
  endif
elseif(cstartpsi(3)+1.ge.myu) then
  starty=myl+(cstartpsi(3)-myu)
  sy1_fg=1
  sy2_fg=spypsi
  up_fg=1
endif

dimc_fg=0

if(sx1_fg.ne.-1) then
  dimc_fg(1)=sx2_fg-sx1_fg+1
endif

dimc_fg(2)=floor(2.0d0/3.0d0*dble(npsiz))

if(sy1_fg.ne.-1) then
  dimc_fg(3)=sy2_fg-sy1_fg+1
  if(sy3_fg.ne.-1) then
    dimc_fg(3)=dimc_fg(3)+sy4_fg-sy3_fg+1
  endif
endif


color=MPI_UNDEFINED
if((sx1_fg.ne.-1).and.(sy1_fg.ne.-1)) then
  call mpi_type_create_subarray(4,[mx,dimc_fg(2),my,2],[dimc_fg(1),dimc_fg(2),dimc_fg(3),2], &
   &     [cstartpsi(1),0,starty,0],mpi_order_fortran, &
   &     mpi_double_precision,stype_fg,ierr)

  call mpi_type_commit(stype_fg,ierr)

  color=1
endif

call mpi_comm_split(mpi_comm_world,color,0,sp_save_comm_fg,ierr)


return
end
