subroutine read_fields(nstep)

use commondata

integer :: nstep
integer :: mx,my,mz,myl,myu
integer :: mx_fg,my_fg,mz_fg,myl_fg,myu_fg

double precision, allocatable, dimension(:,:,:,:) :: tmp,tmp_fg
double precision, allocatable :: inp(:,:,:),inpc(:,:,:,:)
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

allocate(inp(nx,nz,ny))
allocate(inpc(nx/2+1,nz,ny,2))


if(spectral.eq.0)then
 namefile=trim(namedir)//'u_'//numfile//'.dat'
 ! check if u file exists; if u exists we assume that also v and w exist
 inquire(file=trim(namefile),exist=check)

 if(check.eqv..true.)then
  ! reading u
  write(*,*) 'Reading step ',nstep,' out of ',nend
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) inp
  close(666,status='keep')
  call phys_to_spectral(inp,inpc,0)
  call coarse2fine(inpc,uc)
  call spectral_to_phys_fg(uc,u,0)
  ! reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) inp
  close(667,status='keep')
  call phys_to_spectral(inp,inpc,0)
  call coarse2fine(inpc,vc)
  call spectral_to_phys_fg(vc,v,0)
  ! reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) inp
  close(668,status='keep')
  call phys_to_spectral(inp,inpc,0)
  call coarse2fine(inpc,wc)
  call spectral_to_phys_fg(wc,w,0)
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phi_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) inp
    close(669,status='keep')
    call phys_to_spectral(inp,inpc,0)
    call coarse2fine(inpc,phic)
    call spectral_to_phys_fg(phic,phi,0)
  endif
  if(psiflag.eq.1)then
    namefile=trim(namedir)//'psi_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) psi
    close(670,status='keep')
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'T_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) inp
    close(671,status='keep')
    call phys_to_spectral(inp,inpc,0)
    call coarse2fine(inpc,thetac)
    call spectral_to_phys_fg(thetac,theta,0)
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
  inpc=0.0d0
  inpc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  inpc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call coarse2fine(inpc,uc)
  call spectral_to_phys_fg(uc,u,0)
  ! reading v
  namefile=trim(namedir)//'vc_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) tmp
  close(667,status='keep')
  inpc=0.0d0
  inpc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  inpc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call coarse2fine(inpc,vc)
  call spectral_to_phys_fg(vc,v,0)
  ! reading w
  namefile=trim(namedir)//'wc_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) tmp
  close(668,status='keep')
  inpc=0.0d0
  inpc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  inpc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call coarse2fine(inpc,wc)
  call spectral_to_phys_fg(wc,w,0)
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phic_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) tmp
    close(669,status='keep')
    inpc=0.0d0
    inpc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
    inpc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
    call coarse2fine(inpc,phic)
    call spectral_to_phys_fg(phic,phi,0)
  endif
  if(psiflag.eq.1)then
    namefile=trim(namedir)//'psic_fg_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) tmp_fg
    close(670,status='keep')
    psic=0.0d0
    psic(1:mx_fg,1:mz_fg,1:myl_fg,:)=tmp_fg(:,:,1:myl_fg,:)
    psic(1:mx_fg,1:mz_fg,myu_fg:ny*exp_y,:)=tmp_fg(:,:,myl_fg+1:my_fg,:)
    call spectral_to_phys_fg(psic,psi,0)
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'Tc_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) tmp
    close(671,status='keep')
    inpc=0.0d0
    inpc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
    inpc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
    call coarse2fine(inpc,thetac)
    call spectral_to_phys_fg(thetac,theta,0)
  endif
  ! generate paraview output file
  call generate_output(nstep)
 endif
endif

deallocate(tmp)
deallocate(tmp_fg)
deallocate(inp)
deallocate(inpc)

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


! fine grid
dx=xl/dble(exp_x*nx-1)
dy=yl/dble(exp_y*ny-1)

do i=1,exp_x*nx
  xfg(i)=dble(i-1)*dx
enddo

do i=1,exp_y*ny
  yfg(i)=dble(i-1)*dy
enddo

do i=1,exp_z*(nz-1)+1
  zfg(i)=dcos(((i-1)*pi)/((exp_z*(nz-1)+1)-1))
enddo

return
end
