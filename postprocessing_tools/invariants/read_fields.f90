subroutine read_fields(nstep)

use commondata

integer :: nstep
integer :: mx,my,mz,myl,myu
integer :: mx_fg,my_fg,mz_fg,myl_fg,myu_fg

double precision, allocatable, dimension(:,:,:,:) :: tmp,tmp_fg

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
  write(*,*) 'Reading step ',nstep,' out of ',nend
  ! reading u
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(666) u
  close(666,status='keep')
  call phys_to_spectral(u,uc,0)

  ! reading v
  namefile=trim(namedir)//'v_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) v
  close(667,status='keep')
  call phys_to_spectral(v,vc,0)

  ! reading w
  namefile=trim(namedir)//'w_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) w
  close(668,status='keep')
  call phys_to_spectral(w,wc,0)

  ! reading phi
  namefile=trim(namedir)//'phi_'//numfile//'.dat'
  open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(669) phi
  close(669,status='keep')

  ! generate paraview output file
  call get_invariants(nstep)
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
  call spectral_to_phys(uc,u,0)

  ! reading v
  namefile=trim(namedir)//'vc_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) tmp
  close(667,status='keep')
  vc=0.0d0
  vc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  vc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call spectral_to_phys(vc,v,0)

  ! reading w
  namefile=trim(namedir)//'wc_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) tmp
  close(668,status='keep')
  wc=0.0d0
  wc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  wc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call spectral_to_phys(wc,w,0)

  ! reading phi
  namefile=trim(namedir)//'phic_'//numfile//'.dat'
  open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(669) tmp
  close(669,status='keep')
  phic=0.0d0
  phic(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
  phic(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
  call spectral_to_phys(phic,phi,0)

  ! generate paraview output file
  call get_invariants(nstep)
 endif
endif

deallocate(tmp)
deallocate(tmp_fg)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_grid

use commondata

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

return
end
