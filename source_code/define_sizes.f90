subroutine define_sizes

use commondata
use par_size

integer :: rx,ry,rz


! define proper size for each rank in physical space
ry=mod(ny,nycpu)
rz=mod(nz,nzcpu)

fpy=int((ny-ry)/nycpu)
if(mod(rank,nycpu).lt.ry)then
 fpy=int((ny-ry)/nycpu)+1
endif

fpz=int((nz-rz)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.rz)then
 fpz=int((nz-rz)/nzcpu)+1
endif


! define proper size for each rank in spectral space
rx=mod(nx/2+1,nycpu)
ry=mod(ny,nzcpu)

spx=int((nx/2+1-rx)/nycpu)
if(mod(rank,nycpu).lt.rx)then
 spx=int((nx/2+1-rx)/nycpu)+1
endif

spy=int((ny-ry)/nzcpu)
if(floor(real(rank)/real(nycpu)).lt.ry)then
 spy=int((ny-ry)/nzcpu)+1
endif


!write(*,'(i4,1x,2(3(i4),4x))') rank,nx,fpz,fpy,spx,nz,spy

return
end
