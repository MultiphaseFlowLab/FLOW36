program main

! modules
use commondata
use sim_par
use grid
use wavenumber
use velocity
use phase_field

! declaration of variables
integer :: i

! code

call read_input

call generate_grid

call wave_numbers

call create_plan

allocate(u(nx,nz,ny))
allocate(uc(nx/2+1,nz,ny,2))

write(*,'(5(a16))') 'step','t^+','tau_w at z=-1','tau_w at z=+1','u bulk'


open(54,file='./output/wall_shear.dat',form='formatted',status='new',position='append')
write(54,'(5(a16))') 'step','t^+','tau_w at z=-1','tau_w at z=+1','u bulk'
close(54,status='keep')

! time advancement
do i=nstart,nend,ndump
  call wall_shear(i)
enddo




! free all allocated memory and terminate code
deallocate(x)
deallocate(y)
deallocate(z)

deallocate(u)
deallocate(uc)


deallocate(kx)
deallocate(ky)

call destroy_plan


return
end program
