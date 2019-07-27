subroutine lagrangian_tracker

use mpi
use commondata
use particle
use sim_par

#define twowayc twowaycflag
#define tracer tracerflag

double precision, dimension(3) :: for
integer :: i,j

! loop over nset
do j=1,nset
  ! loop over particles
  do i=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,1)+part_index(rank_loc+1,2)
#if tracer == 1
    ! tracers only
    call lagran4(xp(i,:,j),up(i,:,j))
    xp(i,1,j)=xp(i,1,j)+dt*re*up(i,1,j)
    xp(i,2,j)=xp(i,2,j)+dt*re*up(i,2,j)
    xp(i,3,j)=xp(i,3,j)+dt*re*up(i,3,j)
    ! tracer velocity is fluid velocity at particle position
    call lagran4(xp(i,:,j),up(i,:,j))
#elif tracer == 0
    ! inertial particles
    if(abs(stokes(j)).lt.1.0e-15)then
     ! for multiple sets of particles: when using tracer St=0, thus for=NaN
     for=0.0d0
    else
     call calculate_forces(xp(i,:,j),up(i,:,j),for,j)
    endif

    ! time integration
    xp(i,:,j)=xp(i,:,j)+dt*re*up(i,:,j)
    up(i,:,j)=up(i,:,j)+dt*re*for
#endif

    ! check periodicity
    if(xp(i,1,j).lt.0.0d0) then
     xp(i,1,j)=xp(i,1,j)+xl*re
    elseif(xp(i,1,j).gt.xl*re)then
     xp(i,1,j)=xp(i,1,j)-xl*re
    endif
    if(xp(i,2,j).lt.0.0d0) then
     xp(i,2,j)=xp(i,2,j)+yl*re
    elseif(xp(i,2,j).gt.yl*re)then
     xp(i,2,j)=xp(i,2,j)-yl*re
    endif
    ! check rebound (elastic), rebound at particle surface
    if(xp(i,3,j).lt.d_par(j)/2.0d0)then
     xp(i,3,j)=d_par(j)-xp(i,3,j)
     up(i,3,j)=-up(i,3,j)
    elseif(xp(i,3,j).gt.2.0d0*re-d_par(j)/2.0d0)then
     xp(i,3,j)=4.0d0*Re-d_par(j)-xp(i,3,j)
     up(i,3,j)=-up(i,3,j)
    endif
  enddo
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

subroutine calculate_forces(pos,vel,for,case)

use particle
use sim_par
use phase_field

double precision, dimension(3) :: pos,vel
double precision, dimension(3) :: velf,for
double precision :: re_p
integer :: case

#define stokes_drag stokesflag
#define part_grav activategravity

! get fluid velocity at particle position
call lagran4(pos,velf)

! particle Reynolds number
re_p=abs(dsqrt((velf(1)-vel(1))**2+(velf(2)-vel(2))**2+(velf(3)-vel(3))**2))*d_par(case)

! forces calculation
for=0.0d0

#if stokes_drag == 1
! drag force (Stokes drag)
for=(velf-vel)/stokes(case)
#elif stokes_drag == 0
! drag force (with Schiller-Naumann correction)
for=(velf-vel)/stokes(case)*(1.0d0+0.15d0*re_p**0.687)
#endif

#if part_grav == 1
! buoyancy and gravity, in grav ordering is x,z,y
for(1)=for(1)+1.0d0/(re*Fr**2)*(dens_part(case)-1.0d0)/dens_part(case)*grav(1)
for(2)=for(2)+1.0d0/(re*Fr**2)*(dens_part(case)-1.0d0)/dens_part(case)*grav(3)
for(3)=for(3)+1.0d0/(re*Fr**2)*(dens_part(case)-1.0d0)/dens_part(case)*grav(2)
#endif



return
end
