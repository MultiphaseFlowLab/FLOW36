subroutine initialize_phi

use commondata
use par_size
use grid
use phase_field
use velocity_old
use velocity
use sterms
use sim_par
use dual_grid


logical :: checkf,checks

character(len=8) :: time

allocate(phi(nx,fpz,fpy))
allocate(phic(spx,nz,spy,2))
allocate(phi_fg(npsix,fpzpsi,fpypsi))
allocate(phic_fg(spxpsi,npsiz,spypsi,2))

allocate(sphi_o(spx,nz,spy,2))

#define match_dens matched_density
#if match_dens != 1
allocate(ucp(spx,nz,spy,2))
allocate(vcp(spx,nz,spy,2))
allocate(wcp(spx,nz,spy,2))
! initialize previous velocity to current value (first time step initialization)
ucp=uc
vcp=vc
wcp=wc
#endif

s_coeff=((4.0d0*pe*ch**2)/dt)**0.5d0


! initial conditions + transform in physical space
if(in_cond_phi.eq.0)then
  if(rank.eq.0)write(*,*) 'Initializing phi=-1'
  phi=-1.0d0
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.1)then
  if(rank.eq.0)write(*,*) 'Initializing phase fields from data file (parallel read)'
  write(time,'(I8.8)') nt_restart
  if(restart.eq.1)then
    inquire(file=trim(folder)//'/phi_'//time//'.dat',exist=checkf)
    inquire(file=trim(folder)//'/phic_'//time//'.dat',exist=checks)
  else
    checkf=.true.
  endif
  if(checkf.eqv..true.)then
    call read_fields(phi,nt_restart,'phi  ',restart)
    ! transform physical variable to spectral space
    call phys_to_spectral(phi,phic,0)
  elseif(checks.eqv..true.)then
    call read_fields_s(phic,nt_restart,'phic ',restart)
    ! transform to physical space
    call spectral_to_phys(phic,phi,0,0)
  else
    if(rank.eq.0) write(*,'(1x,a,a,a)') 'Missing phase field input file ',time,' , stopping simulation'
    call exit(0)
  endif
elseif(in_cond_phi.eq.2)then
  if(rank.eq.0) write(*,*) 'Initializing phase fields from data file (serial read)'
  call read_fields_serial(phi,nt_restart,'phi  ',restart)
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.3)then
  if(rank.eq.0) write(*,*) 'Initializing 2D drop'
  call drop_2d
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.4)then
  if(rank.eq.0) write(*,*) 'Initializing 3D drop'
  call drop_3d
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.5)then
  if(rank.eq.0) write(*,*) 'Initializing stratified phase field'
  call stratified
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.6)then
  if(rank.eq.0) write(*,*) 'Initializing array of drops in x-y plane'
  call drop_array
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.7)then
  if(rank.eq.0) write(*,*) 'Initializing 2D droplet attached to the bottom wall'
  call drop_wall
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.8)then
  if(rank.eq.0) write(*,*) 'Initializing 2D droplets closed to Kiss <3'
  call drop_kiss
  call phys_to_spectral(phi,phic,0)
elseif(in_cond_phi.eq.9)then
  if(rank.eq.0) write(*,*) 'Initializing Layer of phi=+1'
  call layer
  call phys_to_spectral(phi,phic,0)
else
  if(rank.eq.0)write(*,*) 'Check initial condition value on phi'
  stop
endif


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drop_2d

use commondata
use par_size
use phase_field
use grid

double precision :: radius, height
double precision :: y_c,z_c

integer :: jm,jp,km,kp
integer :: j,k,kg,jg


 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') radius
 read(66,'(f16.8)') height

 close(66,status='keep')

! initialize 2D drop at x=Lx/2, z=height and radius given

jm=minloc(dabs(y-yl/2+radius),1)
jp=minloc(dabs(y-yl/2-radius),1)

km=minloc(dabs(z-height+radius),1)
kp=minloc(dabs(z-height-radius),1)

y_c=(y(jm)+y(jp))/2.0d0
z_c=(z(km)+z(kp))/2.0d0


do j=1,fpy
  do k=1,fpz
    jg=fstart(3)+j
    kg=fstart(2)+k
    phi(:,k,j)=-dtanh((dsqrt((y(jg)-y_c)**2+(z(kg)-z_c)**2)-radius)/(ch*dsqrt(2.0d0)))
  enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine drop_3d

use commondata
use par_size
use phase_field
use grid

double precision :: radius, height
double precision :: x_c,y_c,z_c

integer :: im,ip,jm,jp,km,kp
integer :: jg,kg,i,j,k

 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') radius
 read(66,'(f16.8)') height

 close(66,status='keep')

! initialize 3D drop at x=Lx/2, y=Ly/2, z=height and radius given

im=minloc(dabs(x-xl/2+radius),1)
ip=minloc(dabs(x-xl/2-radius),1)

jm=minloc(dabs(y-yl/2+radius),1)
jp=minloc(dabs(y-yl/2-radius),1)

km=minloc(dabs(z-height+radius),1)
kp=minloc(dabs(z-height-radius),1)

x_c=(x(im)+x(ip))/2.0d0
y_c=(y(jm)+y(jp))/2.0d0
z_c=(z(km)+z(kp))/2.0d0


do j=1,fpy
  do k=1,fpz
    do i=1,nx
      jg=fstart(3)+j
      kg=fstart(2)+k
      phi(i,k,j)=-dtanh((dsqrt((x(i)-x_c)**2+(y(jg)-y_c)**2+(z(kg)-z_c)**2)-radius)/(ch*dsqrt(2.0d0)))
    enddo
  enddo
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stratified

use commondata
use par_size
use phase_field
use grid

double precision :: height,wave_amp_x,wave_freq_x,wave_amp_y,wave_freq_y,pert_amp
double precision :: k_c, rd

integer :: i,j,k,jg,kg

 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') height
 read(66,'(f16.8)') wave_amp_x
 read(66,'(f16.8)') wave_freq_x
 read(66,'(f16.8)') wave_amp_y
 read(66,'(f16.8)') wave_freq_y
 read(66,'(f16.8)') pert_amp

 close(66,status='keep')

! maximum amplitude: surface in [h-amp,h+amp]

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      kg=k+fstart(2)
      jg=j+fstart(3)
      call random_number(rd)
      k_c=height+wave_amp_x*dsin(wave_freq_x*x(i))+wave_amp_y*dsin(wave_freq_y*y(jg))+pert_amp*(2.0d0*rd-1.0d0)
      phi(i,k,j)=-dtanh((k_c-z(kg))/(ch*dsqrt(2.0d0)))
    enddo
  enddo
enddo


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine drop_array

use commondata
use par_size
use phase_field
use grid

double precision :: radius, height
double precision :: xdist,ydist,zdist, dist, x_c, y_c, z_c

integer :: num_x, num_y, num_z
integer :: i,j,k, id,jd,kd, jg,kg

 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') radius
 read(66,'(f16.8)') height
 read(66,'(i8)') num_x
 read(66,'(i8)') num_y
 read(66,'(i8)') num_z

 close(66,status='keep')

! minimum distance between drop centers
! tanh(3)=0.9951~1
dist=2.0d0*(radius+2.0d0*2.0d0**0.5*Ch)


do while(xl.le.dble(num_x)*dist)
  num_x=num_x-1
enddo

do while(yl.le.dble(num_y)*dist)
  num_y=num_y-1
enddo

do while(2.0d0.le.dble(num_z)*dist)
  num_z=num_z-1
enddo

if(rank.eq.0) write(*,'(1x,a,i8,a,i8,a,i8,a,i8,a)')  &
 &        'Start simulation with',num_x,' x',num_y,' x',num_z,' =',num_x*num_y*num_z,' drops'

phi=0.0d0


if(num_x*num_y*num_z.gt.0)then
  xdist=xl/dble(num_x)
  ydist=yl/dble(num_y)
  zdist=2.0d0/dble(num_z)

  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        jg=fstart(3)+j
        kg=fstart(2)+k
        do kd=1,num_z
          if(num_z.eq.1)then
            z_c=height
          else
            z_c=(dble(kd-1)+0.5d0)*zdist-1.0d0
          endif
          do jd=1,num_y
            do id=1,num_x
              x_c=(dble(id-1)+0.5d0)*xdist
              y_c=(dble(jd-1)+0.5d0)*ydist
              phi(i,k,j)=phi(i,k,j)+dtanh((radius-dsqrt((x(i)-x_c)**2+(y(jg)-y_c)**2+(z(kg)-z_c)**2))/(ch*dsqrt(2.0d0)))+1.0d0
            enddo
          enddo
        enddo
        if(phi(i,k,j).ge.1.99999d0)then
          phi(i,k,j)=2.0d0
        endif
        if(phi(i,k,j).le.0.00001d0)then
          phi(i,k,j)=0.0d0
        endif
      enddo
    enddo
  enddo



  phi=phi-1.0d0

endif

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine drop_wall

use commondata
use par_size
use phase_field
use grid

double precision :: radius, height
double precision :: y_c

integer :: jm,jp
integer :: j,k,kg,jg


 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') radius
 read(66,'(f16.8)') height
 close(66,status='keep')

jm=minloc(dabs(y-yl/2+radius),1)
jp=minloc(dabs(y-yl/2-radius),1)

y_c=(y(jm)+y(jp))/2.0d0

do j=1,fpy
  do k=1,fpz
    jg=fstart(3)+j
    kg=fstart(2)+k
    phi(:,k,j)=-dtanh((dsqrt((y(jg)-y_c)**2+(z(kg)+1.0d0)**2)-radius)/(ch*dsqrt(2.0d0)))
  enddo
enddo

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine drop_kiss

use commondata
use par_size
use phase_field
use grid

double precision :: radius
double precision :: y_c,z_c,ygap,zgap

integer :: j,k,kg,jg


 open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

 read(66,'(f16.8)') radius
 read(66,'(f16.8)') ygap
 read(66,'(f16.8)') zgap
 close(66,status='keep')


  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        jg=fstart(3)+j
        kg=fstart(2)+k
        y_c=(yl+ygap)*0.5d0
        z_c=-0.5d0*zgap
        phi(i,k,j)=phi(i,k,j)+dtanh((radius-dsqrt((y(jg)-y_c)**2+(z(kg)-z_c)**2))/(ch*dsqrt(2.0d0)))+1.0d0
!        if(phi(i,k,j).ge.1.99999d0) phi(i,k,j)=2.0d0
!        if(phi(i,k,j).le.0.00001d0) phi(i,k,j)=0.0d0
      enddo
    enddo
  enddo

  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        jg=fstart(3)+j
        kg=fstart(2)+k
        y_c=(yl-ygap)*0.5d0
        z_c= 0.5d0*zgap
        phi(i,k,j)=phi(i,k,j)+dtanh((radius-dsqrt((y(jg)-y_c)**2+(z(kg)-z_c)**2))/(ch*dsqrt(2.0d0)))+1.0d0
!        if(phi(i,k,j).ge.1.99999d0) phi(i,k,j)=2.0d0
!        if(phi(i,k,j).le.0.00001d0) phi(i,k,j)=0.0d0
      enddo
    enddo
  enddo



  phi=phi-1.0d0

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine layer

use commondata
use par_size
use phase_field
use grid

double precision :: thickness,height
double precision :: zlow,ztop

integer :: indc,k,k_g

open(66,file='./sc_compiled/input_phase_field.f90',status='old',form='formatted')

read(66,'(f16.8)') thickness
read(66,'(f16.8)') height
close(66,status='keep')

zlow=height-thickness/2.0d0
ztop=height+thickness/2.0d0

indc=minloc(dabs(z-height),1)

do k=1,fpz
  k_g=fstart(2)+k
  if(k_g.lt.indc)then
    phi(:,k,:)=-dtanh((z(k_g)-ztop)/(ch*dsqrt(2.0d0)))
  else
    phi(:,k,:)=+dtanh((z(k_g)-zlow)/(ch*dsqrt(2.0d0)))
  endif
enddo


return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine destroy_phi

use phase_field
use velocity_old
use sterms

deallocate(phi)
deallocate(phic)
deallocate(phi_fg)
deallocate(phic_fg)

deallocate(sphi_o)

#define match_dens 1
#if match_dens != 1
deallocate(ucp)
deallocate(vcp)
deallocate(wcp)
#endif

return
end
