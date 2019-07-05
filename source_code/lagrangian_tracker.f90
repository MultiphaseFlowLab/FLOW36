subroutine lagrangian_tracker

use mpi
use commondata
use particle
use sim_par

#define twowayc twowaycflag
#define tracer tracerflag

double precision, dimension(3) :: for
integer :: i,j

! loop over particles
do i=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,1)+part_index(rank_loc+1,2)
#if tracer == 1
  ! tracers only
  call lagran4(xp(i,:),up(i,:))
  xp(i,1)=xp(i,1)+dt*re*up(i,1)
  xp(i,2)=xp(i,2)+dt*re*up(i,2)
  xp(i,3)=xp(i,3)+dt*re*up(i,3)
  ! tracer velocity is fluid velocity at particle position
  call lagran4(xp(i,:),up(i,:))
#elif tracer == 0
  ! intertial particles
  call calculate_forces(xp(i,:),up(i,:),for)

! write(*,*) i,xp(i,:),up(i,:),for

  ! time integration
  xp(i,:)=xp(i,:)+dt*re*up(i,:)
  up(i,:)=up(i,:)+dt*re*for



#endif

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


#if twowayc == 1 && tracer == 0
! calculate feedback forces on fluid (two-way coupling only)
fb_x=0.0d0
fb_y=0.0d0
fb_z=0.0d0
do i=0,ntask_sh-1
 if(rank_loc.eq.i)then
  do j=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,1)+part_index(rank_loc+1,2)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_forces(pos,vel,for)

use particle
use sim_par
use phase_field

double precision, dimension(3) :: pos,vel
double precision, dimension(3) :: velf,for
double precision :: re_p

#define stokes_drag stokesflag
#define part_grav activategravity

! get fluid velocity at particle position
call lagran4(pos,velf)

! particle Reynolds number
re_p=abs(dsqrt((velf(1)-vel(1))**2+(velf(2)-vel(2))**2+(velf(3)-vel(3))**2))*d_par

! forces calculation
for=0.0d0

#if stokes_drag == 1
! drag force (Stokes drag)
for=(velf-vel)/stokes
#elif stokes_drag == 0
! drag force (with Schiller-Naumann correction)
for=(velf-vel)/stokes*(1.0d0+0.15d0*re_p**0.687)
#endif

#if part_grav == 1
! buoyancy and gravity, in grav ordering is x,z,y
for(1)=for(1)+1.0d0/(re*Fr**2)*(dens_part-1.0d0)/dens_part*grav(1)
for(2)=for(2)+1.0d0/(re*Fr**2)*(dens_part-1.0d0)/dens_part*grav(3)
for(3)=for(3)+1.0d0/(re*Fr**2)*(dens_part-1.0d0)/dens_part*grav(2)
#endif



return
end
