program grid_check


integer :: nx,ny,nz
double precision, parameter :: pi=3.14159265358979
double precision :: lx,ly,Re,Ch_in

integer :: k, npoints
double precision :: dx,dy,dz,dxp,dyp,dzp,Ch,n_x,n_y,n_z
double precision, allocatable :: z(:), deltaz(:)


open(66,file='./input_par.f90',status='old',form='formatted')
read(66,'(i8)') nx
read(66,'(i8)') ny
read(66,'(i8)') nz
read(66,*)
read(66,'(f16.8)') Lx
read(66,'(f16.8)') Ly
read(66,*)
read(66,'(f16.8)') Re
read(66,*)
read(66,'(f16.8)') Ch_in
close(66,status='keep')

write(*,'(a)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
write(*,*)
write(*,'(1x,a,i6,a,i6,a,i6)') 'Grid: ',nx,' x',ny,' x',nz
write(*,'(1x,a,f6.3,a,f6.3,a,f6.3)') 'Domain size: ',lx,' pi x',ly,' pi x',2.0d0
write(*,'(1x,a,f10.3)') 'Re=',Re
write(*,'(1x,a,f8.5)') 'Ch=',Ch_in
write(*,*)
write(*,'(a)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'


dx=Lx*pi/real(nx-1)
dy=Ly*pi/real(ny-1)
allocate(z(nz))
allocate(deltaz(nz-1))
do k = 1, nz
  z(k)=dcos(((k-1)*pi)/(nz-1))
enddo
deltaz=z(1:nz-1)-z(2:nz)
dz=maxval(deltaz)
deallocate(z)
deallocate(deltaz)

dxp=dx*Re
dyp=dy*Re
dzp=dz*Re

write(*,*)
write(*,'(1x,a)') 'Closed channel, max grid spacing: 15 w.u.'
write(*,'(2x,3(a,f10.3,7x))') 'dx+=',dxp,' dy+=',dyp,' dz+=',dzp
write(*,*)

! minimum number of points across the interface
npoints=4
Ch=dble(npoints-1)*dx/4.1d0
Ch=max(Ch,dble(npoints-1)*dy/4.1d0)
Ch=max(Ch,dble(npoints-1)*dz/4.1d0)

n_x=4.1d0*Ch_in/dx+1.0d0
n_y=4.1d0*Ch_in/dy+1.0d0
n_z=4.1d0*Ch_in/dz+1.0d0

write(*,'(1x,a)') 'Minimum Cahn number (4 points across the interface)'
write(*,'(2x,a,f8.4)') 'Ch=',Ch
write(*,'(1x,a,f8.5)') 'Number of points with Ch=',Ch_in
write(*,'(2x,3(a,f5.1,7x))') 'n_x=',n_x,'n_y=',n_y,'n_z=',n_z
write(*,*)
write(*,'(a)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'


return
end program
