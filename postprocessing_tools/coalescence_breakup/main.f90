program main
! modules
use sim_par
use grid
use phase_field
use commondata
use velocity
! declaration of variables
integer :: i
!integer :: etape, nb,  timestep


! code

call read_input

call generate_grid

allocate(phi(nx,nz,ny))
allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))
! time advancement


open(70,file='./output/mass_center_global.dat',form='formatted',status='new',access='append')
write(70,'(a8,a8,8x,a8,8x,a8,8x,a8,8x,a8)') 'nstep','ndrop','xg','yg','zg','volume'
close(70,status='keep')

!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!Fichier global qui contient le nombre de droplets à chaque timestep
open(42, file='./output/drop_count.dat',form='formatted',status='new', access='append')
!!!!!!!!!!!!!!!!!
write(42,'(a8,a8,a6)') 'nstep','nRe','ndrop'
!write(42,'(a8,2x,a8,2x,a8)' 'nstep','nRe','ndrop' ! avec les espaces : syntaxe à revoir
close(42, status = 'keep')

!!!!!!!!!!!!!!!!!




do i=nstart,nend,ndump
!coeur du code
  call get_interface(i)
!  write(*,'(1x,a,i8.8)') 'Do something at iteration: ',i
! call some subroutine

enddo


!!! Creer nouveau fichier avec juste une seule fois le timestep et le nombre de droplets
!!! à chque à timestep
!! replace write with a read

!open(67,file='./output/nb_droplet.dat',form='formatted',status='new', access='append')
!open(unit = 66,file='./output/mass_center.dat',form='formatted',status='old', action='read')

!read(66,'(a8,a6)') nstep,ndrop
!write(67,'(a8,a6)') nstep, ndrop


!close(67, status = 'keep')
!close(66, status = 'keep')



deallocate(phi)
! free all allocated memory an terminate code
deallocate(x)
deallocate(y)
deallocate(z)


return
end program
