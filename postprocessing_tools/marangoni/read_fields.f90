subroutine read_fields(nstep)

use commondata
use fields

integer :: nstep
integer :: mx,my,mz,myl,myu
integer :: mx_fg,my_fg,mz_fg,myl_fg,myu_fg

double precision, allocatable, dimension(:,:,:,:) :: tmp,tmp_fg
! double precision, allocatable :: inp(:,:,:),inpc(:,:,:,:)

character(len=40) :: namedir,namefile
character(len=8) :: numfile


namedir='../results/'
write(numfile,'(i8.8)') nstep

if(spectral.eq.1)then
 allocate(uc(nx/2,nz,ny,2))
 allocate(vc(nx/2,nz,ny,2))
 allocate(wc(nx/2,nz,ny,2))
 ! read from spectral space
 mx=floor(2.0d0/3.0d0*dble(nx/2+1))
 my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
 mz=floor(2.0d0/3.0d0*dble(nz))
 myl=floor(2.0d0/3.0d0*dble(ny/2+1))
 myu=ny-floor(2.0d0/3.0d0*dble(ny/2))+1
 mx_fg=floor(2.0d0/3.0d0*dble((expx*nx)/2+1))
 my_fg=floor(2.0d0/3.0d0*dble(expy*ny/2+1))+floor(2.0d0/3.0d0*dble(expy*ny/2))
 mz_fg=floor(2.0d0/3.0d0*dble(expz*(nz-1)+1))
 myl_fg=floor(2.0d0/3.0d0*dble(expy*ny/2+1))
 myu_fg=expy*ny-floor(2.0d0/3.0d0*dble(expy*ny/2))+1
 allocate(tmp(mx,mz,my,2))
 allocate(tmp_fg(mx_fg,mz_fg,my_fg,2))
 ! u
 namefile=trim(namedir)//'uc_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) tmp
 close(666,status='keep')
 uc=0.0d0
 uc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
 uc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
 call spectral_to_phys(uc,u,0)
 ! v
 namefile=trim(namedir)//'vc_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) tmp
 close(666,status='keep')
 vc=0.0d0
 vc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
 vc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
 call spectral_to_phys(vc,v,0)
 ! w
 namefile=trim(namedir)//'wc_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) tmp
 close(666,status='keep')
 wc=0.0d0
 wc(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
 wc(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
 call spectral_to_phys(wc,w,0)
 deallocate(uc,vc,wc)
 if(phi_flag.eq.1)then
   allocate(phic(nx/2,nz,ny,2))
   namefile=trim(namedir)//'phic_'//numfile//'.dat'
   open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(666) tmp
   close(666,status='keep')
   phic=0.0d0
   phic(1:mx,1:mz,1:myl,:)=tmp(:,:,1:myl,:)
   phic(1:mx,1:mz,myu:ny,:)=tmp(:,:,myl+1:my,:)
   call spectral_to_phys(phic,phi,0)
   deallocate(phic)
 endif
 if(psi_flag.eq.1)then
   allocate(psic(nxfg/2,nzfg,nyfg,2))
   namefile=trim(namedir)//'psic_fg_'//numfile//'.dat'
   open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(666) tmp
   close(666,status='keep')
   psic=0.0d0
   psic(1:mx_fg,1:mz_fg,1:myl_fg,:)=tmp_fg(:,:,1:myl_fg,:)
   psic(1:mx_fg,1:mz_fg,myu_fg:nyfg,:)=tmp_fg(:,:,myl_fg+1:my_fg,:)
   call spectral_to_phys_fg(psic,psi,0)
   deallocate(psic)
 endif
elseif(spectral.eq.0)then
 ! read from physical space
 ! u
 namefile=trim(namedir)//'u_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) u
 close(666,status='keep')
 ! v
 namefile=trim(namedir)//'v_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) v
 close(666,status='keep')
 ! w
 namefile=trim(namedir)//'w_'//numfile//'.dat'
 open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 read(666) w
 close(666,status='keep')
 if(phi_flag.eq.1)then
   namefile=trim(namedir)//'phi_'//numfile//'.dat'
   open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(666) phi
   close(666,status='keep')
 endif
 if(psi_flag.eq.1)then
   namefile=trim(namedir)//'psi_fg_'//numfile//'.dat'
   open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
   read(666) psi
   close(666,status='keep')
 endif
elseif(spectral.eq.-1)then
 ! error in parameters
 write(*,*) 'Error in saving frequency'
endif


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
