subroutine read_input

use commondata
use pdf_calc

double precision :: Lx, Ly

 open(unit=66,file='./input_curvature.f90',form='formatted',status='old',action='read')

 read(66,'(i5)') nx
 read(66,'(i5)') ny
 read(66,'(i5)') nz
 read(66,*)
 read(66,'(f16.8)') Lx
 read(66,'(f16.8)') Ly
 read(66,*)
 read(66,'(i8)') nstart
 read(66,'(i8)') nend
 read(66,'(i8)') ndump
 read(66,*)
 read(66,'(f16.8)') Re
 read(66,'(f16.8)') dt
 read(66,*)
 read(66,'(i5)') spectral
 read(66,*)
 read(66,'(i8)') nset
 read(66,'(f16.8)') threshold
 read(66,*)
 read(66,'(i5)') exp_x
 read(66,'(i5)') exp_y
 read(66,'(i5)') exp_z

 xl=Lx*pi
 yl=Ly*pi

 close(66)


 if(rank.eq.0) call print_start

return
end


subroutine print_start

use commondata

write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,'(1x,a)') '                 starting curvature calculation                       '
write(*,'(1x,a,3(i5,a))') 'Nx * Ny * Nz = ',nx,' * ',ny,' * ',nz,' '
write(*,'(1x,a,3(f8.2,a))') 'Lx * Ly * Lz = (',xl,' * ',yl,' * ',2.0d0,')*Re'
write(*,'(1x,a,f8.2)') 'Re = ',re
if(spectral.eq.1)then
 write(*,*) 'Reading fields in modal space'
else
 write(*,*) 'Reading fields in physical space'
endif
write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',ndump
write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,*)

return
end
