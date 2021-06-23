subroutine write_output(stats)

use commondata

double precision :: stats(nz,13)

integer :: i

open(66,status='replace',file='./output/statistics.dat',form='formatted')

write(66,'(a,i8,a,3(i5),a)') 'Statistics gathered on ',counter,' flow fields, on a ',nx,ny,nz,' grid (nx,ny,nz)'
write(66,'(13(a12,2x))') 'z','u mean','v mean','w mean','u rms','v rms','w rms','u skw','v skw','w skw','u flt','v flt','w flt'
write(66,*)

do i=1,nz
 write(66,'(f12.5,2x,12(es12.5,2x))') stats(i,1:13)
enddo

 close(66,status='keep')

return
end

