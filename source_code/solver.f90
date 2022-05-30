subroutine solver(ntime)

use commondata
use par_size
use sim_par
use velocity
use sterms
use grid
use phase_field
use surfactant
use temperature
use particle
use dual_grid

double precision, allocatable, dimension(:,:,:,:) :: s1,s2,s3,sphi,hphi,spsi,hpsi,stheta,htheta
double precision, dimension(spx,nz,spy,2) :: h1,h2,h3,h,omega

integer :: ntime

#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define partflag particlecompflag
#define twowayc twowaycflag
#define boussinnesq boussinnesqcompflag
#define match_dens matched_density
#define expx expansionx
#define expy expansiony
#define expz expansionz
#define cpiflag cpicompflag


! flow part
if(rank.lt.flow_comm_lim)then
 allocate(s1(spx,nz,spy,2))
 allocate(s2(spx,nz,spy,2))
 allocate(s3(spx,nz,spy,2))

 !$acc kernels
 s1=0.0d0
 s2=0.0d0
 s3=0.0d0
 !$acc end kernels

 ! calculate convective terms of N-S equation and store them in s1,s2,s3
 call convective_ns(s1,s2,s3)

! add mean pressure gradient to S term
#if cpiflag == 0
#if match_dens == 2
  ! rescale NS equation if rhor > 1 for improved stability
  !$acc kernels
  s1=s1-sgradpx/rhor
  s2=s2-sgradpy/rhor
  !$acc end kernels
#else
  !$acc kernels
  s1=s1-sgradpx
  s2=s2-sgradpy
  !$acc end kernels
#endif
#elif cpiflag == 1
#if match_dens == 2
  ! rescale NS equation if rhor > 1 for improved stability
  !$acc kernels
  s1=s1-sgradpx*dabs(gradpx)/rhor
  s2=s2-sgradpy*dabs(gradpy)/rhor
  !$acc end kernels
#else
  !$acc kernels
  s1=s1-sgradpx*dabs(gradpx)
  s2=s2-sgradpy*dabs(gradpy)
  !$acc end kernels
#endif
#endif


! add non-linear part of phase field to N-S non-linear terms
#if phiflag == 1
 call phi_non_linear(s1,s2,s3)
#endif

! add temperature contribution to N-S non-linear terms (Boussinesq)
#if tempflag == 1
#if boussinnesq == 1
 ! write(*,*) grav*Ra/(16.0d0*Pr*Re**2)
 s1=s1-grav(1)*Ra/(16.0d0*Pr*Re**2)*thetac
 s2=s2-grav(3)*Ra/(16.0d0*Pr*Re**2)*thetac
 s3=s3-grav(2)*Ra/(16.0d0*Pr*Re**2)*thetac
#endif
#endif

 ! time integration of the non-linear terms
 ! for first time step the code uses an explicit Euler algorithm,
 ! while from the second time step on an Adams-Bashforth algorithm
 if (ntime.eq.nstart+1)then
   call euler(s1,s2,s3,h1,h2,h3)
 else
   call adams_bashforth(s1,s2,s3,h1,h2,h3)
 endif

 !GPU kernels lead to performance derating (GPU version)
 s1_o=s1
 s2_o=s2
 s3_o=s3

 deallocate(s1)
 deallocate(s2)
 deallocate(s3)

 ! add linear contribution to history term
 call hist_term(h1,h2,h3,h)

 ! solve Helmholtz equation for w
 call calculate_w(h)

 ! solve Helmholtz equation for omega_z and store in h1
 call calculate_omega(h1,h2,omega)

 ! calculate u,v from continuity and vorticity definition
 call calculate_uv(omega,h1,h2)

 call chop_modes(uc)
 call chop_modes(vc)
 call chop_modes(wc)

 ! only if surfactant calculated on finer grid than phase field
#if expx != 1 || expy != 1 || expz != 1
 ! for phase field calculation
 call spectral_to_phys(uc,u,1)
 call spectral_to_phys(vc,v,1)
 call spectral_to_phys(wc,w,1)
 ! for surfactant calculation
 call coarse2fine(uc,uc_fg)
 call spectral_to_phys_fg(uc_fg,u_fg,1)
 call coarse2fine(vc,vc_fg)
 call spectral_to_phys_fg(vc_fg,v_fg,1)
 call coarse2fine(wc,wc_fg)
 call spectral_to_phys_fg(wc_fg,w_fg,1)
#else
 call spectral_to_phys(uc,u,1)
 call spectral_to_phys(vc,v,1)
 call spectral_to_phys(wc,w,1)
 uc_fg=uc
 vc_fg=vc
 wc_fg=wc
 u_fg=u
 v_fg=v
 w_fg=w
#endif

 ! Cahn-Hilliard equation solution
#if phiflag == 1
 allocate(sphi(spx,nz,spy,2))
 allocate(hphi(spx,nz,spy,2))
 ! calculate non-linear terms of Cahn-Hilliard equation
 call sterm_ch(sphi)

 if (ntime.eq.nstart+1)then
   call euler_phi(sphi,hphi)
 else
   call adams_bashforth_phi(sphi,hphi)
 endif

 sphi_o=sphi

 deallocate(sphi)

 ! history term
 hphi=hphi+phic

 call calculate_phi(hphi)

 deallocate(hphi)

 call chop_modes(phic)


 ! Cahn-Hilliard equation for surfactant
#if psiflag == 1
 allocate(spsi(spxpsi,npsiz,spypsi,2))
 allocate(hpsi(spxpsi,npsiz,spypsi,2))
 ! calculate non-linear terms of surfactant Cahn-Hilliard equation
 call sterm_surf(spsi)

 if (ntime.eq.nstart+1)then
   call euler_psi(spsi,hpsi)
 else
   call adams_bashforth_psi(spsi,hpsi)
 endif

 spsi_o=spsi

 deallocate(spsi)

 hpsi=hpsi+psic_fg

 call calculate_psi(hpsi)

 deallocate(hpsi)

 call chop_modes_fg(psic_fg)
#endif

! surfactant part executed only iff the phase field is activated
#endif


! Temperature transport equation
#if tempflag == 1
 allocate(stheta(spx,nz,spy,2))
 allocate(htheta(spx,nz,spy,2))

 ! calculate non-linear terms of temperature equation
 call sterm_temp(stheta)

 if (ntime.eq.nstart+1)then
   call euler_theta(stheta,htheta)
 else
   call adams_bashforth_theta(stheta,htheta)
 endif

 stheta_o=stheta

 deallocate(stheta)

 ! assemble history term
 call hist_term_temp(htheta)

 call calculate_theta(htheta)

 deallocate(htheta)

 call chop_modes(thetac)
#endif
endif


! particle part
#if partflag == 1
if(rank.ge.leader)then
 call lagrangian_tracker
endif

call get_velocity

#if twowayc == 1
call get_2WCforces
#endif
#endif


return
end
