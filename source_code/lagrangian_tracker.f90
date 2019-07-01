subroutine lagrangian_tracker

use mpi
use commondata
use particle
use sim_par

#define twowayc twowaycflag

integer :: i,j


! loop over particles
do i=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,2)
  call lagran4(xp(i,:),up(i,:))
  ! tracers only
  xp(i,1)=xp(i,1)+dt*re*up(i,1)
  xp(i,2)=xp(i,2)+dt*re*up(i,2)
  xp(i,3)=xp(i,3)+dt*re*up(i,3)
  ! tracer velocity is fluid velocity at particle position
  call lagran4(xp(i,:),up(i,:))

  ! check periodicity
  if(xp(i,1).lt.0.0d0) then
   xp(i,1)=xp(i,1)+xl*re
  elseif(xp(i,1).gt.xl*re)then
   xp(i,1)=xp(i,1)-xl*re
  endif
  if(xp(i,2).lt.0.0d0) then
   xp(i,2)=xp(i,2)+yl*re
  elseif(xp(i,2).gt.yl*re)then
   xp(i,2)=xp(i,2)-yl*re
  endif
  ! check rebound (elastic)
  if(xp(i,3).lt.0)then
   xp(i,3)=-xp(i,3)
   up(i,3)=-up(i,3)
  elseif(xp(i,3).gt.2.0d0*re)then
   xp(i,3)=4.0d0*Re-xp(i,3)
   up(i,3)=-up(i,3)
  endif

enddo

call mpi_win_fence(0,window_xp,ierr)
call mpi_win_fence(0,window_up,ierr)


#if twowayc == 1
! calculate feedback forces on fluid (two-way coupling only)
fb_x=0.0d0
fb_y=0.0d0
fb_z=0.0d0
do i=0,ntask_sh-1
 if(rank_loc.eq.i)then
  do j=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,2)
   ! add force contribution, to be filled in
   ! fb_x=fb_x+something(j)
   ! fb_y=fb_y+something(j)
   ! fb_z=fb_z+something(j)
  enddo
 endif
 ! synchronize windows
 call mpi_win_fence(0,window_fx,ierr)
 call mpi_win_fence(0,window_fy,ierr)
 call mpi_win_fence(0,window_fz,ierr)
enddo

#endif

return
end subroutine
