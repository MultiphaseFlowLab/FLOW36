subroutine read_fields(nstep)

use commondata

integer :: nstep

character(len=40) :: namedir
character(len=8) :: numfile

namedir='../results/'
write(numfile,'(i8.8)') nstep


if(spectral.eq.1)then
!  open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! ! open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
!  read(66) phic
!  close(66,status='keep')
!  call spectral_to_phys(phic,phi,0)
!if(psiflag.eq.1)then
!!  open(66,file=trim(namedir)//'psic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
!  open(66,file=trim(namedir)//'psic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
!  read(66) psic
!  close(66,status='keep')
!endif
else
! read in physical space
 open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
 read(66) phi
 close(66,status='keep')

 open(66,file=trim(namedir)//'psi_fg_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'psi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
 read(66) psi
 close(66,status='keep')
endif


write(*,*) maxval(phi),minval(phi) !, maxval(psi),minval(psi)


call get_interface(nstep)

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

write(*,'(a,2(f9.4),a)') 'x range: ',x(1),x(nx),' (outer units)'
write(*,'(a,2(f9.4),a)') 'y range: ',y(1),y(ny),' (outer units)'
write(*,'(a,2(f9.4),a)') 'z range: ',z(1),z(nz),' (outer units)'


return
end
