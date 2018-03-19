subroutine calc_mean

use commondata

integer :: i,j,k

double precision :: mean_t(nz)

mean_t=0.0d0

do k=1,nz
  do j=1,ny
    do i=1,nx
      mean_t(k)=mean_t(k)+u(i,k,j)
    enddo
  enddo
enddo

mean_t=mean_t/dble(nx*ny)

counter=counter+1

do k=1,nz
  do j=1,ny
    do i=1,nx
      mean_u(k)=mean_u(k)+u(i,k,j)
      rms_u(k)=rms_u(k)+(u(i,k,j)-mean_t(k))**2
      skw_u(k)=skw_u(k)+(u(i,k,j)-mean_t(k))**3
      flt_u(k)=flt_u(k)+(u(i,k,j)-mean_t(k))**4
    enddo
  enddo
enddo




return
end
