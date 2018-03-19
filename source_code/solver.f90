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

double precision, allocatable, dimension(:,:,:,:) :: s1, s2, s3, sphi,hphi,spsi,hpsi,stheta,htheta
double precision, dimension(spx,nz,spy,2) :: h1, h2, h3, h, omega

integer :: ntime

#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define boussinnesq boussinnesqcompflag
#define match_dens matched_density

allocate(s1(spx,nz,spy,2))
allocate(s2(spx,nz,spy,2))
allocate(s3(spx,nz,spy,2))


s1=0.0d0
s2=0.0d0
s3=0.0d0


! calculate convective terms of N-S equation and store them in s1,s2,s3
call convective_ns(s1,s2,s3)


! add mean pressure gradient to S term
#if match_dens == 2
! rescale NS equation if rhor > 1 for improved stability
s1=s1-sgradpx/rhor
s2=s2-sgradpy/rhor
#else
s1=s1-sgradpx
s2=s2-sgradpy
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



! Cahn-Hilliard equation solution
#if phiflag == 1
allocate(sphi(spx,nz,spy,2))
allocate(hphi(spx,nz,spy,2))
!sphi=0.0d0
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



! Cahn-Hilliard equation for surfactant
#if psiflag == 1
allocate(spsi(spx,nz,spy,2))
allocate(hpsi(spx,nz,spy,2))
! calculate non-linear terms of surfactant Cahn-Hilliard equation
call sterm_surf(spsi)

if (ntime.eq.nstart+1)then
  call euler_psi(spsi,hpsi)
else
  call adams_bashforth_psi(spsi,hpsi)
endif

spsi_o=spsi

deallocate(spsi)

hpsi=hpsi+psic

call calculate_psi(hpsi)

if (rank.eq.0) then
  print*,'Average Concentration Surfactant (Not Normalized):',psic(1,1,1,1)
endif

deallocate(hpsi)

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

#endif



!call chop_modes(uc)
!call chop_modes(vc)
!call chop_modes(wc)
!#if phiflag == 1
!call chop_modes(phic)
!#endif

return
end
