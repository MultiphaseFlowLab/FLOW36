subroutine read_fields(nstep)

use commondata

integer :: nstep
integer :: mx,my,mz,myl,myu
integer :: mx_fg,my_fg,mz_fg,myl_fg,myu_fg

double precision, allocatable, dimension(:,:,:,:) :: tmp,tmp_fg
double precision :: psi_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny)
double precision :: psic_fg((exp_x*nx)/2+1,exp_z*(nz-1)+1,exp_y*ny,2)

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

allocate(tmp(mx,mz,my,2))
allocate(tmp_fg(mx_fg,mz_fg,my_fg,2))

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
  ! reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  ! reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phi_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) phi
    close(669,status='keep')
  endif
  if(psiflag.eq.1)then
    ! ! reading psi
    ! namefile=trim(namedir)//'psi_'//numfile//'.dat'
    ! open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    ! read(670) psi
    ! close(670,status='keep')
    ! reading psi fine grid
    namefile=trim(namedir)//'psi_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) psi_fg
    close(670,status='keep')
    call phys_to_spectral_fg(psi_fg,psic_fg,0)
    call fine2coarse(psic_fg,psic)
    call spectral_to_phys(psic,psi,0)
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'T_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) theta
    close(671,status='keep')
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
  read(666) tmp
  close(666,status='keep')
  uc=0.0d0
  uc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  uc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  ! reading v
  namefile=trim(namedir)//'vc_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) tmp
  close(667,status='keep')
  vc=0.0d0
  vc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  vc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  ! reading w
  namefile=trim(namedir)//'wc_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) tmp
  close(668,status='keep')
  wc=0.0d0
  wc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  wc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phic_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) tmp
    close(669,status='keep')
    phic=0.0d0
    phic(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
    phic(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
    call spectral_to_phys(phic,phi,0)
  endif
  if(psiflag.eq.1)then
    ! ! reading psi
    ! namefile=trim(namedir)//'psic_'//numfile//'.dat'
    ! open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    ! read(670) tmp
    ! close(670,status='keep')
    ! psic=0.0d0
    ! psic(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
    ! psic(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
    ! call spectral_to_phys(psic,psi,0)
    ! reading psi fg
    namefile=trim(namedir)//'psic_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) tmp_fg
    close(670,status='keep')
    psic_fg=0.0d0
    psic_fg(1:mx_fg,1:mz_fg,1:myl_fg,:)=tmp_fg(:,:,1:myl_fg,:)
    psic_fg(1:mx_fg,1:mz_fg,myu_fg:ny*exp_y,:)=tmp_fg(:,:,myl_fg+1:my_fg,:)
    call fine2coarse(psic_fg,psic)
    call spectral_to_phys(psic,psi,0)
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'thetac_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) tmp
    close(671,status='keep')
    thetac=0.0d0
    thetac(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
    thetac(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
    call spectral_to_phys(thetac,theta,0)
  endif
  ! transform variables to physical space
  call spectral_to_phys(uc,u,0)
  call spectral_to_phys(vc,v,0)
  call spectral_to_phys(wc,w,0)
  ! generate paraview output file
  call generate_output(nstep)
 endif
endif

deallocate(tmp)

return
end




subroutine read_grid

use commondata

character(len=40) :: namefile

 write(namefile,'(a)') '../results/x.dat'
 open(1,file=trim(namefile),form='unformatted',access='stream',status='old')
 read(1) x
 close(1,status='keep')


 write(namefile,'(a)') '../results/y.dat'
 open(2,file=trim(namefile),form='unformatted',access='stream',status='old')
 read(2) y
 close(2,status='keep')


 write(namefile,'(a)') '../results/z.dat'
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old')
 read(3) z
 close(3,status='keep')

 x=x*re
 y=y*re
 z=z*re

return
end
