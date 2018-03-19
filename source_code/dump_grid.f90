subroutine dump_grid

use commondata
use grid

 open(unit=1,file=trim(folder)//'/x.dat',status='replace',form='unformatted',access='stream')
 open(unit=2,file=trim(folder)//'/y.dat',status='replace',form='unformatted',access='stream')
 open(unit=3,file=trim(folder)//'/z.dat',status='replace',form='unformatted',access='stream')

 write(1) x

 write(2) y

 write(3) z

 close(1,status='keep')
 close(2,status='keep')
 close(3,status='keep')

return
end
