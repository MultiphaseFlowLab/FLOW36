subroutine read_input

use commondata

double precision :: Lx, Ly

 open(unit=66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')

 read(66,*)
 read(66,*)
 read(66,*)
 read(66,*)
 read(66,*)
 read(66,'(i5)') nx
 read(66,'(i5)') ny
 read(66,'(i5)') nz
 read(66,*)
 read(66,*)
 close(66)

 open(unit=66,file='input_int.f90',form='formatted',status='old',action='read')

 read(66,'(i8)') nstart
 read(66,'(i8)') ndump
 read(66,'(i8)') sdump
 read(66,'(i8)') nend

 close(66)


 call print_start

return
end










subroutine print_start

use commondata

write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,'(1x,a)') '                   starting post-processing                           '
write(*,'(1x,a,3(i5,a))') 'Nx * Ny * Nz = ',nx,' * ',ny,' * ',nz,' '
write(*,'(1x,a,3(f8.2,a))') 'Lx * Ly * Lz = (',xl,' * ',yl,' * ',2.0d0,')*Re'
write(*,'(1x,a,f8.2)') 'Re = ',re
if(spectral.eq.1)then
 write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',sdump
else
 write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',ndump
endif
write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,*)

return
end
