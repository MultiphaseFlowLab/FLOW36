subroutine read_fields(nstep)

use commondata

integer :: nstep

character(len=40) :: namedir,namefile
character(len=8) :: numfile


namedir='../results/'
write(numfile,'(i8.8)') nstep

namefile=trim(namedir)//'w_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) w
close(666,status='keep')

namefile=trim(namedir)//'phi_'//numfile//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) phi
close(666,status='keep')

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
