subroutine get_interface(nstep)

use commondata

integer :: nstep
integer :: i,j,k,id,jd,kd
integer :: top(nx,nz,ny),s_drop(nx,nz,ny)
integer :: drop_count


do j=1,ny
  do k=1,nz
    do i=1,nx
      if(phi(i,k,j).ge.0.0d0)then
        top(i,k,j)=1
      else
        top(i,k,j)=0
      endif
    enddo
  enddo
enddo

drop_count=0

! flood fill algorithm
do jd=1,ny
 do kd=1,nz
  do id=1,nx
   if(top(id,kd,jd).gt.0)then
    drop_count=drop_count+1
    write(*,'(2x,a,i3,a)') 'New drop, ',drop_count,' drops'
    ! single drop part
    s_drop=0
    s_drop(id,kd,jd)=1
    ! flood fill algorithm
    call flood_fill(top,s_drop,id,jd,kd)
    ! remove drops already done from top
    top=top-s_drop
    ! new drop calculation
   endif
  enddo
 enddo
enddo

write(*,'(2x,a,i4)') 'Number of drops: ',drop_count
write(*,*)

open(42,file='./output/drop_count.dat',access='append',form='formatted',status='old')
 write(42,'(i16,2x,es16.6,2x,i16)') nstep,dble(nstep)*dt*re,drop_count
close(42,status='keep')

return
end

