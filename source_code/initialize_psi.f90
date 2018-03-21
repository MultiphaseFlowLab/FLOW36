subroutine initialize_psi

use commondata
use par_size
use grid
use surfactant
use sim_par
use sterms
use phase_field
use dual_grid

double precision :: psi_k

logical :: checkf,checks

character(len=8) :: time

allocate(psi(nx,fpz,fpy))
allocate(psic(spx,nz,spy,2))
allocate(psi_fg(npsix,fpzpsi,fpypsi))
allocate(psic_fg(spxpsi,npsiz,spypsi,2))

allocate(spsi_o(spxpsi,npsiz,spypsi,2))


if(in_cond_psi.eq.0)then
  if(rank.eq.0) write(*,*) 'Initializing constant surfactant field'
  open(66,file='./sc_compiled/input_surfactant.f90',status='old',form='formatted')
  read(66,'(f16.6)') psi_k
  close(66,status='keep')
  psi_fg=psi_k
  call phys_to_spectral_fg(psi_fg,psic_fg,0)
  call fine2coarse(psic_fg,psic)
  call spectral_to_phys(psic,psi,0)
elseif(in_cond_psi.eq.1)then
  if(rank.eq.0) write(*,*) 'Initializing surfactant from data file (parallel read)'
  write(time,'(I8.8)') nt_restart
  if(restart.eq.1)then
    inquire(file=trim(folder)//'/psi_'//time//'.dat',exist=checkf)
    inquire(file=trim(folder)//'/psic_'//time//'.dat',exist=checks)
  else
    checkf=.true.
  endif
  if(checkf.eqv..true.)then
    call read_fields(psi,nt_restart,'psi  ',restart)
    ! transform physical variable to spectral space
    call phys_to_spectral(psi,psic,0)
  elseif(checks.eqv..true.)then
    call read_fields_s(psic,nt_restart,'psic ',restart)
    ! transform to physical space
    call spectral_to_phys(psic,psi,0)
  else
    if(rank.eq.0) write(*,'(1x,a,a,a)') 'Missing surfactant input file ',time,' , stopping simulation'
    call exit(0)
  endif
  call coarse2fine(phic,phic_fg)
  call spectral_to_phys_fg(phic_fg,phi_fg,0)
elseif(in_cond_psi.eq.2)then
  if(rank.eq.0) write(*,*) 'Initializing equilibrium profile'
  call psi_eq
  call phys_to_spectral_fg(psi_fg,psic_fg,0)
  call fine2coarse(psic_fg,psic)
  call spectral_to_phys(psic,psi,1)
elseif(in_cond_psi.eq.3)then
  if(rank.eq.0) write(*,*) 'Initializing equilibrium profile multiplied with Y gradient'
  ! linear profile, calculated on coarse grid
  call psi_grady
  call phys_to_spectral(psi,psic,0)
  call coarse2fine(psic,psic_fg)
  call spectral_to_phys_fg(psic_fg,psi_fg,0)
elseif(in_cond_psi.eq.4)then
  if(rank.eq.0) write(*,*) 'Initializing equilibrium profile multiplied with Z gradient'
  ! linear profile, calculated on coarse grid
  call psi_gradz
  call phys_to_spectral(psi,psic,0)
  call coarse2fine(psic,psic_fg)
  call spectral_to_phys_fg(psic_fg,psi_fg,0)
elseif(in_cond_psi.eq.5)then
  if(rank.eq.0) write(*,*) 'Initializing Diffusion Test 2D (NON ANCORA FATTO !!!)'
  call psi_diff_test
  call phys_to_spectral(psi,psic,0)
  call coarse2fine(psic,psic_fg)
  call spectral_to_phys_fg(psic_fg,psi_fg,0)
else
  if(rank.eq.0)write(*,*) 'Check initial condition value on psi'
  stop
endif

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psi_eq

use commondata
use par_size
use phase_field
use grid
use phase_field
use surfactant
use dual_grid

double precision psi_k, psi_c
integer :: i,j,k

open(66,file='./sc_compiled/input_surfactant.f90',status='old',form='formatted')
read(66,'(f16.6)') psi_k
close(66,status='keep')

call coarse2fine(phic,phic_fg)
call spectral_to_phys_fg(phic_fg,phi_fg,0)

do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
       psi_c=exp(-(1.0d0-(phi_fg(i,k,j))**2)**2/(4.0d0*P_i)-(1.0d0-(phi_fg(i,k,j))**2)/(2.0d0*Ex*P_i))
       psi_fg(i,k,j)=psi_k/(psi_k+psi_c*(1.0d0-psi_k))
    enddo
  enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psi_grady

use commondata
use par_size
use phase_field
use grid
use phase_field
use surfactant

double precision psi_k, psi_c
integer :: i,j,k,kg,jg

open(66,file='./sc_compiled/input_surfactant.f90',status='old',form='formatted')
read(66,'(f16.6)') psi_k
close(66,status='keep')

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      jg=fstart(3)+j
      kg=fstart(2)+k
      psi_c=exp(-(1.0d0-(phi(i,k,j))**2)**2/(4.0d0*P_i)-(1.0d0-(phi(i,k,j))**2)/(2.0d0*Ex*P_i))
      psi(i,k,j)=psi_k/(psi_k+psi_c*(1.0d0-psi_k))-psi_k
      psi(i,k,j)=psi(i,k,j)*(-dtanh(y(jg)-0.48d0*yl)+1.0d0)
      psi(i,k,j)=psi_k+psi(i,k,j)
    enddo
  enddo
enddo

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psi_gradz

use commondata
use par_size
use phase_field
use grid
use phase_field
use surfactant

double precision psi_k, psi_c
integer :: i,j,k,kg,jg

open(66,file='./sc_compiled/input_surfactant.f90',status='old',form='formatted')
read(66,'(f16.6)') psi_k
close(66,status='keep')

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      jg=fstart(3)+j
      kg=fstart(2)+k
      psi_c=exp(-(1.0d0-(phi(i,k,j))**2)**2/(4.0d0*P_i)-(1.0d0-(phi(i,k,j))**2)/(2.0d0*Ex*P_i))
      psi(i,k,j)=psi_k/(psi_k+psi_c*(1.0d0-psi_k))-psi_k
      psi(i,k,j)=-0.5d0*(-1.0d0+z(kg))*psi(i,k,j)
      psi(i,k,j)=psi_k+psi(i,k,j)
    enddo
  enddo
enddo

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psi_diff_test

use commondata
use par_size
use phase_field
use grid
use phase_field
use surfactant
use sim_par

double precision psi_k, psi_c, theta
integer :: i,j,k,kg,jg

open(66,file='./sc_compiled/input_surfactant.f90',status='old',form='formatted')
read(66,'(f16.6)') psi_k
close(66,status='keep')

do j=1,fpy
  do k=1,fpz
    do i=1,nx
       jg=fstart(3)+j
       kg=fstart(2)+k
       psi_c=exp(-(1.0d0-(phi(i,k,j))**2)**2/(4.0d0*P_i)-(1.0d0-(phi(i,k,j))**2)/(2.0d0*Ex*P_i))
       psi(i,k,j)=psi_k/(psi_k+psi_c*(1.0d0-psi_k))-psi_k
       theta=datan(z(kg)/(y(jg)-0.5d0*yl))
       if ((y(jg)-0.5d0*yl) < 0.0d0) theta=theta+pi
       psi(i,k,j)=psi(i,k,j)*0.5d0*(1-dcos(theta))
       psi(i,k,j)=psi_k+psi(i,k,j)
    enddo
  enddo
enddo

return
end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine destroy_psi

use surfactant
use sterms

deallocate(psi)
deallocate(psic)
deallocate(psi_fg)
deallocate(psic_fg)

deallocate(spsi_o)

return
end
