subroutine read_fields(nstep)

use commondata

integer :: nstep
integer :: mx,my,mz,myl,myu
integer :: mx_fg,my_fg,mz_fg,myl_fg,myu_fg

double precision, allocatable, dimension(:,:,:,:) :: fl,fl_fg,sh,sh_fg

character(len=40) :: namedir,namefile
character(len=8) :: numfile

logical :: check

namedir='../results/'
write(numfile,'(i8.8)') nstep

mx=floor(2.0d0/3.0d0*dble(nx/2+1))
my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
mz=floor(2.0d0/3.0d0*dble(nz))

myl=floor(2.0d0/3.0d0*dble(ny/2+1))
myu=ny-floor(2.0d0/3.0d0*dble(ny/2))+1

mx_fg=floor(2.0d0/3.0d0*dble((exp_x*nx)/2+1))
my_fg=floor(2.0d0/3.0d0*dble(exp_y*ny/2+1))+floor(2.0d0/3.0d0*dble(exp_y*ny/2))
mz_fg=floor(2.0d0/3.0d0*dble(exp_z*(nz-1)+1))

myl_fg=floor(2.0d0/3.0d0*dble(exp_y*ny/2+1))
myu_fg=exp_y*ny-floor(2.0d0/3.0d0*dble(exp_y*ny/2))+1


allocate(fl(nx/2+1,nz,ny,2))
allocate(fl_fg((exp_x*nx)/2+1,exp_z*(nz-1)+1,exp_y*ny,2))

allocate(sh(mx,mz,my,2))
allocate(sh_fg(mx_fg,mz_fg,my_fg,2))



if(spectral.eq.0)then
 namefile=trim(namedir)//'u_'//numfile//'.dat'
 ! check if u file exists; if u exists we assume that also v and w exist
 inquire(file=trim(namefile),exist=check)

 if(check.eqv..true.)then
  ! reading u
  write(*,*) 'Reading step ',nstep,' out of ',nend
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) u
  close(666,status='keep')
  ! port to fine grid
  call phys_to_spectral(u,fl,0)
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,u_fg,0)
  ! reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  ! port to fine grid
  call phys_to_spectral(v,fl,0)
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,v_fg,0)
  ! reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')
  ! port to fine grid
  call phys_to_spectral(w,fl,0)
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,w_fg,0)
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phi_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) phi
    close(669,status='keep')
    ! port to fine grid
    call phys_to_spectral(phi,fl,0)
    call coarse2fine(fl,fl_fg)
    call spectral_to_phys_fg(fl_fg,phi_fg,0)
  endif
  if(psiflag.eq.1)then
    ! reading psi fine grid
    namefile=trim(namedir)//'psi_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) psi_fg
    close(670,status='keep')
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'T_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) theta
    close(671,status='keep')
    ! port to fine grid
    call phys_to_spectral(theta,fl,0)
    call coarse2fine(fl,fl_fg)
    call spectral_to_phys_fg(fl_fg,theta_fg,0)
  endif
  ! generate paraview output file
  call generate_output(nstep)
 endif
else
 namefile=trim(namedir)//'uc_'//numfile//'.dat'
 ! check if u file exists; if u exists we assume that also v and w exist
 inquire(file=trim(namefile),exist=check)

 if(check.eqv..true.)then
  ! reading u
  write(*,*) 'Reading step ',nstep,' out of ',nend
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) sh
  close(666,status='keep')
  fl=0.0d0
  fl(1:mx,1:mz,1:myl,:)=sh(:,:,1:myl,:)
  fl(1:mx,1:mz,myu:ny,:)=sh(:,:,myl+1:my,:)
  ! port to fine grid
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,u_fg,0)
  ! reading v
  namefile=trim(namedir)//'vc_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) sh
  close(667,status='keep')
  fl=0.0d0
  fl(1:mx,1:mz,1:myl,:)=sh(:,:,1:myl,:)
  fl(1:mx,1:mz,myu:ny,:)=sh(:,:,myl+1:my,:)
  ! port to fine grid
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,v_fg,0)
  ! reading w
  namefile=trim(namedir)//'wc_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) sh
  close(668,status='keep')
  fl=0.0d0
  fl(1:mx,1:mz,1:myl,:)=sh(:,:,1:myl,:)
  fl(1:mx,1:mz,myu:ny,:)=sh(:,:,myl+1:my,:)
  ! port to fine grid
  call coarse2fine(fl,fl_fg)
  call spectral_to_phys_fg(fl_fg,w_fg,0)
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phic_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) sh
    close(669,status='keep')
    fl=0.0d0
    fl(1:mx,1:mz,1:myl,:)=sh(:,:,1:myl,:)
    fl(1:mx,1:mz,myu:ny,:)=sh(:,:,myl+1:my,:)
    ! port to fine grid
    call coarse2fine(fl,fl_fg)
    call spectral_to_phys_fg(fl_fg,phi_fg,0)
  endif
  if(psiflag.eq.1)then
    ! reading psi fg
    namefile=trim(namedir)//'psic_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) sh_fg
    close(670,status='keep')
    fl_fg=0.0d0
    fl_fg(1:mx_fg,1:mz_fg,1:myl_fg,:)=sh_fg(:,:,1:myl_fg,:)
    fl_fg(1:mx_fg,1:mz_fg,myu_fg:ny*exp_y,:)=sh_fg(:,:,myl_fg+1:my_fg,:)
    call spectral_to_phys_fg(fl_fg,psi_fg,0)
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'Tc_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) sh
    close(671,status='keep')
    fl=0.0d0
    fl(1:mx,1:mz,1:myl,:)=sh(:,:,1:myl,:)
    fl(1:mx,1:mz,myu:ny,:)=sh(:,:,myl+1:my,:)
    ! port to fine grid
    call coarse2fine(fl,fl_fg)
    call spectral_to_phys_fg(fl_fg,theta_fg,0)
  endif
  ! generate paraview output file
  call generate_output(nstep)
 endif
endif


deallocate(fl,fl_fg,sh,sh_fg)

return
end




subroutine read_grid

use commondata

double precision :: dx,dy
integer :: i

character(len=40) :: namefile

 write(namefile,'(a)') '../results/x.dat'
 open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
! open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(1) x
 close(1,status='keep')


 write(namefile,'(a)') '../results/y.dat'
 open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
! open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(2) y
 close(2,status='keep')


 write(namefile,'(a)') '../results/z.dat'
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
! open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(3) z
 close(3,status='keep')

 x=x*re
 y=y*re
 z=(z+1.0d0)*re


dx=maxval(x)/dble(nx*exp_x-1)
x_fg(1)=0.0d0
do i=2,nx*exp_x
  x_fg(i)=x_fg(i-1)+dx
enddo

dy=maxval(y)/dble(ny*exp_y-1)
y_fg(1)=0.0d0
do i=2,ny*exp_y
  y_fg(i)=y_fg(i-1)+dy
enddo

do i=1,(nz-1)*exp_z+1
  z_fg(i)=dcos(((i-1)*pi)/((nz-1)*exp_z+1-1))
enddo
z_fg=(z_fg+1.0d0)*re

return
end
