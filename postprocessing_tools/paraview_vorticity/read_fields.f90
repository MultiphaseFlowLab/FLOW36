subroutine read_fields(nstep)

use commondata

integer :: nstep

character(len=40) :: namedir,namefile
character(len=8) :: numfile

logical :: check

namedir='../results/'
write(numfile,'(i8.8)') nstep


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
    ! reading psi
    namefile=trim(namedir)//'psi_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) psi
    close(670,status='keep')
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
  read(666) uc
  close(666,status='keep')
  ! reading v
  namefile=trim(namedir)//'vc_'//numfile//'.dat'
  open(667,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(667) vc
  close(667,status='keep')
  ! reading w
  namefile=trim(namedir)//'wc_'//numfile//'.dat'
  open(668,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  read(668) wc
  close(668,status='keep')
  if(phiflag.eq.1)then
    ! reading phi
    namefile=trim(namedir)//'phic_'//numfile//'.dat'
    open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(669) phic
    close(669,status='keep')
    call spectral_to_phys(phic,phi,0)
  endif
  if(psiflag.eq.1)then
    ! reading psi
    namefile=trim(namedir)//'psic_'//numfile//'.dat'
    open(670,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(670) psic
    close(670,status='keep')
    call spectral_to_phys(psic,psi,0)
  endif
  if(tempflag.eq.1)then
    ! reading theta
    namefile=trim(namedir)//'thetac_'//numfile//'.dat'
    open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
    read(671) thetac
    close(671,status='keep')
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

return
end




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
