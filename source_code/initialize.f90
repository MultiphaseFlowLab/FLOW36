subroutine initialize

use commondata
use sim_par
use par_size
use grid
use velocity
use sterms
use wavenumber
use dual_grid
use particle

double precision, dimension(nx,fpz,fpy) :: fgradpx,fgradpy
double precision :: beta2(spx,spy),k2l(spx,spy)

integer :: j,k

logical :: checkf,checks

character(len=8) :: time


allocate(u(nx,fpz,fpy))
allocate(v(nx,fpz,fpy))
allocate(w(nx,fpz,fpy))

allocate(uc(spx,nz,spy,2))
allocate(vc(spx,nz,spy,2))
allocate(wc(spx,nz,spy,2))

allocate(u_fg(npsix,fpzpsi,fpypsi))
allocate(v_fg(npsix,fpzpsi,fpypsi))
allocate(w_fg(npsix,fpzpsi,fpypsi))

allocate(uc_fg(spxpsi,npsiz,spypsi,2))
allocate(vc_fg(spxpsi,npsiz,spypsi,2))
allocate(wc_fg(spxpsi,npsiz,spypsi,2))

allocate(s1_o(spx,nz,spy,2))
allocate(s2_o(spx,nz,spy,2))
allocate(s3_o(spx,nz,spy,2))

if(part_flag.eq.0)then
 allocate(forx(nx,fpz,fpy))
 allocate(fory(nx,fpz,fpy))
 allocate(forz(nx,fpz,fpy))
endif

! transform mean presure gradient in spectral space
allocate(sgradpx(spx,nz,spy,2))
allocate(sgradpy(spx,nz,spy,2))

fgradpx=gradpx
fgradpy=gradpy

call phys_to_spectral(fgradpx,sgradpx,0)
call phys_to_spectral(fgradpy,sgradpy,0)


! declarations of differents initial conditions for velocity
if(in_cond.eq.0)then
  if(rank.eq.0) write(*,*) 'Initializing zero flow field'
  u=0.0d0
  v=0.0d0
  w=0.0d0
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
elseif(in_cond.eq.1)then
  if(rank.eq.0) write(*,*) 'Initializing laminar Poiseuille flow in x direction'
  do k=1,fpz
   u(:,k,:)=re/2.0d0*(1.0d0-z(fstart(2)+k)*z(fstart(2)+k))
  enddo
  v=0.0d0
  w=0.0d0
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
elseif(in_cond.eq.2)then
  if(rank.eq.0) write(*,*) 'Initializing laminar Poiseuille flow in y direction'
  do k=1,fpz
   v(:,k,:)=re/2.0d0*(1.0d0-z(fstart(2)+k)*z(fstart(2)+k))
  enddo
  u=0.0d0
  w=0.0d0
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
elseif(in_cond.eq.3)then
  if(rank.eq.0) write(*,*) 'Initializing velocity fields from data file'
  write(time,'(I8.8)') nt_restart
  if(restart.eq.1)then
    inquire(file=trim(folder)//'/u_'//time//'.dat',exist=checkf)
    inquire(file=trim(folder)//'/uc_'//time//'.dat',exist=checks)
  else
    checkf=.true.
  endif
  if(checkf.eqv..true.)then
    call read_fields(u,nt_restart,'u    ',restart)
    call read_fields(v,nt_restart,'v    ',restart)
    call read_fields(w,nt_restart,'w    ',restart)
    ! transform physical variable to spectral space
    call phys_to_spectral(u,uc,0)
    call phys_to_spectral(v,vc,0)
    call phys_to_spectral(w,wc,0)
  elseif(checks.eqv..true.)then
    call read_fields_s(uc,nt_restart,'uc   ',restart)
    call read_fields_s(vc,nt_restart,'vc   ',restart)
    call read_fields_s(wc,nt_restart,'wc   ',restart)
    ! transform to physical space
    call spectral_to_phys(uc,u,0)
    call spectral_to_phys(vc,v,0)
    call spectral_to_phys(wc,w,0)
  else
    if(rank.eq.0) write(*,'(1x,a,a,a)') 'Missing flow input files ',time,' , stopping simulation'
    call exit(0)
  endif
elseif(in_cond.eq.4)then
  if(rank.eq.0) write(*,*) 'Initializing velocity fields from data file (serial read)'
  call read_fields_serial(u,nt_restart,'u    ',restart)
  call read_fields_serial(v,nt_restart,'v    ',restart)
  call read_fields_serial(w,nt_restart,'w    ',restart)
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
elseif(in_cond.eq.5)then
  if(rank.eq.0) write(*,*) 'Initializing shear flow y direction'
  u=0.0d0
  w=0.0d0
  do k=1,fpz
   v(:,k,:)=z(fstart(2)+k)
  enddo
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
elseif(in_cond.eq.6)then
  if(rank.eq.0) write(*,*) 'Initializing shear flow x direction'
  v=0.0d0
  w=0.0d0
  do k=1,fpz
   u(:,k,:)=z(fstart(2)+k)
  enddo
  ! transform physical variable to spectral space
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)
else
  if(rank.eq.0) write(*,*) 'Dafuq? Check in_cond input value'
  stop
endif


! definition of boundary conditions
! index 1: z=-1 ; index 2: z=1
zp(1)=-1.0d0
zp(2)=1.0d0
if(bc_up.eq.0)then
 p_u(2)=1.0d0
 q_u(2)=0.0d0
 r_u(2)=0.0d0
 p_v(2)=1.0d0
 q_v(2)=0.0d0
 r_v(2)=0.0d0
 p_w(2)=1.0d0
 q_w(2)=0.0d0
 r_w(2)=0.0d0
 p_o(2)=1.0d0
 q_o(2)=0.0d0
 r_o(2)=0.0d0
elseif(bc_up.eq.1)then
 p_u(2)=0.0d0
 q_u(2)=1.0d0
 r_u(2)=0.0d0
 p_v(2)=0.0d0
 q_v(2)=1.0d0
 r_v(2)=0.0d0
 p_w(2)=1.0d0
 q_w(2)=0.0d0
 r_w(2)=0.0d0
 p_o(2)=0.0d0
 q_o(2)=1.0d0
 r_o(2)=0.0d0
elseif(bc_up.eq.2)then
  p_u(2)=1.0d0
  q_u(2)=0.0d0
  r_u(2)=0.0d0
  p_v(2)=1.0d0
  q_v(2)=0.0d0
  r_v(2)=1.0d0*dble(nx*ny)
  p_w(2)=1.0d0
  q_w(2)=0.0d0
  r_w(2)=0.0d0
  p_o(2)=1.0d0
  q_o(2)=0.0d0
  r_o(2)=0.0d0
elseif(bc_up.eq.3)then
  p_u(2)=1.0d0
  q_u(2)=0.0d0
  r_u(2)=1.0d0*dble(nx*ny)
  p_v(2)=1.0d0
  q_v(2)=0.0d0
  r_v(2)=0.0d0
  p_w(2)=1.0d0
  q_w(2)=0.0d0
  r_w(2)=0.0d0
  p_o(2)=1.0d0
  q_o(2)=0.0d0
  r_o(2)=0.0d0
else
 write(*,*) 'Check your boundary conditions at z=+1'
 stop
endif

if(bc_low.eq.0)then
 p_u(1)=1.0d0
 q_u(1)=0.0d0
 r_u(1)=0.0d0
 p_v(1)=1.0d0
 q_v(1)=0.0d0
 r_v(1)=0.0d0
 p_w(1)=1.0d0
 q_w(1)=0.0d0
 r_w(1)=0.0d0
 p_o(1)=1.0d0
 q_o(1)=0.0d0
 r_o(1)=0.0d0
elseif(bc_low.eq.1)then
 p_u(1)=0.0d0
 q_u(1)=1.0d0
 r_u(1)=0.0d0
 p_v(1)=0.0d0
 q_v(1)=1.0d0
 r_v(1)=0.0d0
 p_w(1)=1.0d0
 q_w(1)=0.0d0
 r_w(1)=0.0d0
 p_o(1)=0.0d0
 q_o(1)=1.0d0
 r_o(1)=0.0d0
elseif(bc_low.eq.2)then
  p_u(1)=1.0d0
  q_u(1)=0.0d0
  r_u(1)=0.0d0
  p_v(1)=1.0d0
  q_v(1)=0.0d0
  r_v(1)=-1.0d0*dble(nx*ny)
  p_w(1)=1.0d0
  q_w(1)=0.0d0
  r_w(1)=0.0d0
  p_o(1)=1.0d0
  q_o(1)=0.0d0
  r_o(1)=0.0d0
elseif(bc_low.eq.3)then
  p_u(1)=1.0d0
  q_u(1)=0.0d0
  r_u(1)=-1.0d0*dble(nx*ny)
  p_v(1)=1.0d0
  q_v(1)=0.0d0
  r_v(1)=0.0d0
  p_w(1)=1.0d0
  q_w(1)=0.0d0
  r_w(1)=0.0d0
  p_o(1)=1.0d0
  q_o(1)=0.0d0
  r_o(1)=0.0d0
else
 write(*,*) 'Check your boundary conditions at z=-1'
 stop
endif


! solve auxiliary Helmholtz problems (influence matrix method, needed to solve w Helmholtz equation)
allocate(wa2(spx,nz,spy))
allocate(wa3(spx,nz,spy))

wa2=0.0d0
wa3=0.0d0

do j=1,spy
  do i=1,spx
    beta2(i,j)=(1.0d0+gamma*k2(i+cstart(1),j+cstart(3)))/gamma
    k2l(i,j)=k2(i+cstart(1),j+cstart(3))
  enddo
enddo


! solve auxiliary problem 2 (see notes)
! psi_2
call helmholtz_red(wa2,beta2,[1.0d0,1.0d0],[0.0d0,0.0d0],[1.0d0,0.0d0],zp)
! w_2
call helmholtz_red(wa2,k2l,[1.0d0,1.0d0],[0.0d0,0.0d0],[0.0d0,0.0d0],zp)

! solve auxiliary problem 3 (see notes)
! psi_3
call helmholtz_red(wa3,beta2,[1.0d0,1.0d0],[0.0d0,0.0d0],[0.0d0,1.0d0],zp)
! w_3
call helmholtz_red(wa3,k2l,[1.0d0,1.0d0],[0.0d0,0.0d0],[0.0d0,0.0d0],zp)

! wa2: solution of Helmholtz problem 2 (w_2)
! wa3: solution of Helmholtz problem 3 (w_3)
! wiil be used in calculate_w to solve the Helmholtz problem with the influence matrix method


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine destroy

use velocity
use particle
use sterms

deallocate(u)
deallocate(v)
deallocate(w)

deallocate(uc)
deallocate(vc)
deallocate(wc)

deallocate(u_fg)
deallocate(v_fg)
deallocate(w_fg)

deallocate(uc_fg)
deallocate(vc_fg)
deallocate(wc_fg)

deallocate(s1_o)
deallocate(s2_o)
deallocate(s3_o)

deallocate(sgradpx)
deallocate(sgradpy)

deallocate(wa2)
deallocate(wa3)

if(part_flag.eq.0)then
 deallocate(forx,fory,forz)
endif

return
end
