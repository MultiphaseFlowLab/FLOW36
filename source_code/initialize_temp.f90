subroutine initialize_temp

use commondata
use par_size
use grid
use temperature
use sterms
use sim_par

double precision :: theta_k
character(len=8) :: time
logical :: checkf,checks

allocate(theta(nx,fpz,fpy))
allocate(thetac(spx,nz,spy,2))
allocate(stheta_o(spx,nz,spy,2))


if(in_cond_theta.eq.0)then
  if(rank.eq.0)write(*,*) 'Initializing constant temperature'
  open(66,file='./sc_compiled/input_temperature.f90',form='formatted',status='old')
  read(66,'(f16.8)') theta_k
  close(66,status='keep')
  theta=theta_k
  call phys_to_spectral(theta,thetac,0)
elseif(in_cond_theta.eq.1)then
  if(rank.eq.0)write(*,*) 'Reading from data file (parallel read)'
  write(time,'(I8.8)') nt_restart
  if(restart.eq.1)then
    inquire(file=trim(folder)//'/T_'//time//'.dat',exist=checkf)
    inquire(file=trim(folder)//'/Tc_'//time//'.dat',exist=checks)
  else
    checkf=.true.
  endif
  if(checkf.eqv..true.)then
    call read_fields(theta,nt_restart,'T    ',restart)
    ! transform physical variable to spectral space
    call phys_to_spectral(theta,thetac,0)
  elseif(checks.eqv..true.)then
    call read_fields_s(thetac,nt_restart,'Tc   ',restart)
    ! transform to physical space
    call spectral_to_phys(thetac,theta,0)
  else
    if(rank.eq.0) write(*,'(1x,a,a,a)') 'Missing temperature input file ',time,' , stopping simulation'
    call exit(0)
  endif
else
  if(rank.eq.0)write(*,*) 'Check initial condition value on theta'
  stop
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine destroy_theta

use temperature
use sterms

deallocate(theta)
deallocate(thetac)
deallocate(stheta_o)

return
end
