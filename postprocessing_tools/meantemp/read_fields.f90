subroutine read_fields(nstep)

use commondata

integer :: nstep
character(len=40) :: namedir,namefile
character(len=8) :: numfile


namedir='../results/'
write(numfile,'(i8.8)') nstep



!write(*,*) 'Reading step ',nstep,' out of ',nend,' , flow and particles'
! reading phi
namefile=trim(namedir)//'phi_fg_'//numfile//'.dat'
open(669,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(669) phi
close(669,status='keep')

! reading theta 
namefile=trim(namedir)//'T_'//numfile//'.dat'
open(671,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(671) theta
close(671,status='keep')


return
end


