program meant

use commondata
implicit none
integer :: i,nstep

write(*,*) "Read input files"

call read_input

allocate(phi(nx,nz,ny))
allocate(theta(nx,nz,ny))

nstep=nstart

! read fields (phi and T)
do i=nstart,nend,ndump

 call read_fields(nstep)
 
 call meantemp(nstep)

 nstep=nstep + ndump
enddo



!clean up everything
deallocate(phi,theta)

end program meant


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine meantemp(step)

use commondata
integer :: i,k,j,step,ncar,ndrop,bufz
double precision :: tcar,tdrop,time

bufz=30
ndrop=0
ncar=0
tdrop=0.0d0
tcar=0.0d0

time=dt*step*re

! conditional average
do i=1,nx
 do j=1,ny
  do k=bufz,nz-bufz
   if (phi(i,k,j) .gt. 0.0d0) then 
    tdrop= tdrop + theta(i,k,j) 
    ndrop= ndrop +1
   endif
   if (phi(i,k,j) .lt. 0.0d0) then 
    tcar= tcar + theta(i,k,j) 
    ncar= ncar +1
    endif
  enddo
 enddo
enddo

!normalize
tdrop=tdrop/ndrop
tcar=tcar/ncar

write(*,*) step, time, tcar, tdrop

!write output
!write(namefile,'(a,i8.8,a)') './output/temp_mean',nstep,'.dat'
!open(10,file=namefile,status='new',form='formatted')
!write(10,'(2(es16.8))') ax(m), mean(m)
!close(10)


return 
end


