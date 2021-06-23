subroutine read_fields(nstep)

use commondata

integer :: nstep

character(len=40) :: namedir
character(len=8) :: numfile

namedir='../results/'
write(numfile,'(i8.8)') nstep

if(spectral.eq.1)then
! read in spectral space
 open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(66) phic
 close(66,status='keep')
 ! outward pointing normal: change sign of phic
 phic=-phic
 call spectral_to_phys(phic,phi,0)
else
! read in physical space
 open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(66) phi
 close(66,status='keep')
 ! outward pointing normal: change sign of phi
 phi=-phi
 ! call phys_to_spectral(phi,phic,0)
endif

open(676,file=trim(namedir)//'psi_fg_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
read(676) psi
close(676,status='keep')

! phi rescaled by a factor 50 to reduce oscillation in gradient calculation
! rescaling does not affect normal calculation (normal is normalized to modulo 1)
! phi=phi/50.0d0
! phic=phic/50.0d0

call calc_curvature(nstep)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_grid

use commondata

character(len=40) :: namefile

 write(namefile,'(a)') '../results/x.dat'
 open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(1) x
 close(1,status='keep')


 write(namefile,'(a)') '../results/y.dat'
 open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(2) y
 close(2,status='keep')


 write(namefile,'(a)') '../results/z.dat'
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(3) z
 close(3,status='keep')

return
end
