subroutine save_flow_comm(i)

use commondata
use sim_par
use velocity
use phase_field
use stats
use surfactant
use temperature

#define machine machineflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define fdump physical_dump_frequency
#define sdump spectral_dump_frequency

integer :: i
character(len=5) :: namevar


#if fdump > 0
 ! save fields at end of timestep according to ndump value
 if(mod(i,ndump).eq.0) then
  if(rank.eq.0) write(*,*) 'saving solution in physical space'
  call spectral_to_phys(uc,u,0,0)
  namevar='u'
  call write_output(u,i,namevar)
  call spectral_to_phys(vc,v,0,0)
  namevar='v'
  call write_output(v,i,namevar)
  call spectral_to_phys(wc,w,0,0)
  namevar='w'
  call write_output(w,i,namevar)
#if phiflag == 1
   call spectral_to_phys(phic,phi,0,0)
   namevar='phi'
   call write_output(phi,i,namevar)
#if psiflag == 1
    call spectral_to_phys_fg(psic_fg,psi_fg,0)
    namevar='psi'
    call write_output_fg(psi_fg,i,namevar)
#endif
#endif
#if tempflag == 1
  call spectral_to_phys(thetac,theta,0,0)
  namevar='T'
  call write_output(theta,i,namevar)
#endif
 endif
#endif

#if sdump > 0
#if machine != 2 && machine != 5
 ! save fields at end of timestep according to ndump value
 if(mod(i,sdump).eq.0) then
  if(rank.eq.0) write(*,*) 'saving solution in spectral space'
  namevar='uc'
  call write_output_spectral(uc,i,namevar)
  namevar='vc'
  call write_output_spectral(vc,i,namevar)
  namevar='wc'
  call write_output_spectral(wc,i,namevar)
#if phiflag == 1
   namevar='phic'
   call write_output_spectral(phic,i,namevar)
#if psiflag == 1
    namevar='psic'
    call write_output_spectral_fg(psic_fg,i,namevar)
#endif
#endif
#if tempflag == 1
  namevar='Tc'
  call write_output_spectral(thetac,i,namevar)
#endif
 endif
#endif
#endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_flow_comm_final(i)

use commondata
use sim_par
use velocity
use phase_field
use stats
use surfactant
use temperature

#define machine machineflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define fdump physical_dump_frequency
#define sdump spectral_dump_frequency

integer :: i
character(len=5) :: namevar


#if fdump > 0
 ! save fields at end of timestep according to ndump value
 if(mod(i,ndump).ne.0) then
  if(rank.eq.0) write(*,*) 'Writing final fields in physical space'
  call spectral_to_phys(uc,u,0,0)
  namevar='u'
  call write_output(u,i,namevar)
  call spectral_to_phys(vc,v,0,0)
  namevar='v'
  call write_output(v,i,namevar)
  call spectral_to_phys(wc,w,0,0)
  namevar='w'
  call write_output(w,i,namevar)
#if phiflag == 1
   call spectral_to_phys(phic,phi,0,0)
   namevar='phi'
   call write_output(phi,i,namevar)
#if psiflag == 1
    call spectral_to_phys_fg(psic_fg,psi_fg,0)
    namevar='psi'
    call write_output_fg(psi_fg,i,namevar)
#endif
#endif
#if tempflag == 1
  call spectral_to_phys(thetac,theta,0,0)
  namevar='T'
  call write_output(theta,i,namevar)
#endif
 endif
#endif

#if sdump > 0
#if machine != 2 && machine != 5
 ! save fields at end of timestep according to ndump value
 if(mod(i,sdump).ne.0) then
  if(rank.eq.0) write(*,*) 'Writing final fields in spectral space'
  namevar='uc'
  call write_output_spectral(uc,i,namevar)
  namevar='vc'
  call write_output_spectral(vc,i,namevar)
  namevar='wc'
  call write_output_spectral(wc,i,namevar)
#if phiflag == 1
   namevar='phic'
   call write_output_spectral(phic,i,namevar)
#if psiflag == 1
    namevar='psic'
    call write_output_spectral_fg(psic_fg,i,namevar)
#endif
#endif
#if tempflag == 1
  namevar='Tc'
  call write_output_spectral(thetac,i,namevar)
#endif
 endif
#endif
#endif

return
end
