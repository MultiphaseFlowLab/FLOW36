subroutine coarse2fine(varc,varf)

use commondata
use par_size
use dual_grid
use mpi

double precision :: varc(spx,nz,spy,2), varf(spxpsi,npsiz,spypsi,2)
double precision, allocatable, asynchronous, dimension(:,:,:,:) :: bufs,bufr
double precision :: rx,ry

integer :: i,il,iu,jl,ju
integer :: npx,npy,numel
integer :: tag
integer :: ii,jj

#define expx expansionx
#define expy expansiony
#define expz expansionz


! if(rank.eq.0) write(*,*)
! if(rank.eq.0) write(*,*) 'Coarse to fine'
!$acc kernels
varf=0.0d0
!$acc end kernels

#if expx == 1 && expy == 1 && expz == 1
  ! if no grid expansion
  !$acc kernels
  varf=varc
  !$acc end kernels
#elif expx == 1 && expy == 1
  ! if grid expanded in z only no need for MPI communications
  varf(:,1:nz,:,:)=varc(:,1:nz,:,:)
#else
  rx=mod(nx/2+1,nycpu)
  npx=int((nx/2+1-rx)/nycpu)
  if(rx.ne.0) npx=npx+1
  ry=mod(ny,nzcpu)
  npy=int((ny-ry)/nzcpu)
  if(ry.ne.0) npy=npy+1
  numel=npx*nz*npy*2
  allocate(bufs(npx,nz,npy,2))
  allocate(bufr(npx,nz,npy,2))

  do i=1,ntask
    if(c2fadd(i,1).ne.-1)then
      if((rank.eq.c2fadd(i,1)).and.(rank.eq.c2fadd(i,2)))then
        ! rank sender=rank receiver, no MPI communications
        ! write(*,*) 'rank ',rank,c2fadd(i,1),' int copy to ',c2fadd(i,2)
        il=max(cstart(1)+1,cstartpsi(1)+1)
        iu=min(cstart(1)+1+spx-1,cstartpsi(1)+1+spxpsi-1)
        ! copy 1:ny/2+1
        jl=max(cstart(3)+1,cstartpsi(3)+1)
        ju=min(cstart(3)+1+spy-1,cstartpsi(3)+1+spypsi-1,ny/2+1)
        !$acc kernels
        if(jl.le.ju) varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,jl-cstartpsi(3):ju-cstartpsi(3),:)= &
          & varc(il-cstart(1):iu-cstart(1),1:nz,jl-cstart(3):ju-cstart(3),:)
        !$acc end kernels 
        ! copy ny/2+2:ny
        ! minimum distance from top, maximum distance from top
        jl=min(ny-cstart(3),npsiy-cstartpsi(3),ny/2-1)
        ju=max(ny-(cstart(3)+1+spy-1),npsiy-(cstartpsi(3)+1+spypsi-1))
        !$acc kernels
        if(jl.ge.ju) varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,npsiy-jl-cstartpsi(3)+1:npsiy-ju-cstartpsi(3),:)= &
          & varc(il-cstart(1):iu-cstart(1),1:nz,ny-jl-cstart(3)+1:ny-ju-cstart(3),:)
        !$acc end kernels
      else
        ! write(*,*) 'rank ',rank,c2fadd(i,1),' MPI send to ',c2fadd(i,2)
        !$acc kernels
        bufs=0.0d0
        bufs(1:spx,1:nz,1:spy,1:2)=varc
        !$acc end kernels

        ! tag=1000000*c2fadd(i,1)+c2fadd(i,2)
        tag=13
        if(rank.eq.c2fadd(i,1)) then
          ! send
          !$acc host_data use_device(bufs)
          call mpi_ssend(bufs,numel,MPI_double_precision,c2fadd(i,2),tag,flow_comm,ierr)
          !$acc end host_data
        elseif(rank.eq.c2fadd(i,2)) then
          ! receive, sender=receiver already accounted for in internal copy
          !$acc host_data use_device(bufr)
          call mpi_recv(bufr,numel,MPI_double_precision,c2fadd(i,1),tag,flow_comm,mpi_status_ignore,ierr)
          !$acc end host_data
          ! update fine grid variable
          ii=mod(c2fadd(i,1),nycpu)+1
          jj=int((c2fadd(i,1)-(ii-1))/nycpu)+1
          il=max(cg_size(ii,jj,4),cstartpsi(1)+1)
          iu=min(cg_size(ii,jj,4)+cg_size(ii,jj,2)-1,cstartpsi(1)+1+spxpsi-1)
          ! copy 1:ny/2+1
          jl=max(cg_size(ii,jj,5),cstartpsi(3)+1)
          ju=min(cg_size(ii,jj,5)+cg_size(ii,jj,3)-1,cstartpsi(3)+1+spypsi-1,ny/2+1)
          !$acc kernels
          if(jl.le.ju) varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,jl-cstartpsi(3):ju-cstartpsi(3),:)= &
            & bufr(il-cg_size(ii,jj,4)+1:iu-cg_size(ii,jj,4)+1,1:nz,jl-cg_size(ii,jj,5)+1:ju-cg_size(ii,jj,5)+1,:)
          !$acc end kernels
          ! copy ny/2+2:ny
          ! minimum distance from top, maximum distance from top
          jl=min(ny-cg_size(ii,jj,5)+1,npsiy-cstartpsi(3),ny/2-1)
          ju=max(ny-(cg_size(ii,jj,5)+cg_size(ii,jj,3)-1),npsiy-(cstartpsi(3)+1+spypsi-1))
          !$acc kernels
          if(jl.ge.ju) varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,npsiy-jl-cstartpsi(3)+1:npsiy-ju-cstartpsi(3),:)= &
            & bufr(il-cg_size(ii,jj,4)+1:iu-cg_size(ii,jj,4)+1,1:nz,ny-jl+1-cg_size(ii,jj,5)+1:ny-ju+1-cg_size(ii,jj,5),:)
          !$acc end kernels
        endif
      endif
    endif
  enddo
  deallocate(bufs)
  deallocate(bufr)
#endif


! renormalize
!$acc kernels
varf=varf*dble(exp_x*exp_y)
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fine2coarse(varf,varc)

use commondata
use par_size
use dual_grid
use mpi

double precision :: varc(spx,nz,spy,2), varf(spxpsi,npsiz,spypsi,2)
double precision, allocatable, asynchronous, dimension(:,:,:,:) :: bufs,bufr
double precision :: rx,ry

integer :: i,il,iu,jl,ju
integer :: npx,npy,numel
integer :: tag
integer :: ii,jj

#define expx expansionx
#define expy expansiony
#define expz expansionz


! if(rank.eq.0) write(*,*)
! if(rank.eq.0) write(*,*) 'Fine to coarse'

!$acc kernels
varc=0.0d0
!$acc end kernels

#if expx == 1 && expy == 1 && expz == 1
  !$acc kernels
  varc=varf
  !$acc end kernels
#elif expx == 1 && expy == 1
  ! if grid expanded in z only no need for MPI communications
  varc(:,1:nz,:,:)=varf(:,1:nz,:,:)
#else
  rx=mod(npsix/2+1,nycpu)
  npx=int((npsix/2+1-rx)/nycpu)
  if(rx.ne.0) npx=npx+1
  ry=mod(npsiy,nzcpu)
  npy=int((npsiy-ry)/nzcpu)
  if(ry.ne.0) npy=npy+1
  numel=npx*npsiz*npy*2
  allocate(bufs(npx,npsiz,npy,2))
  allocate(bufr(npx,npsiz,npy,2))
  do i=1,ntask
    if(f2cadd(i,1).ne.-1)then
      if((rank.eq.f2cadd(i,1)).and.(rank.eq.f2cadd(i,2)))then
        ! rank sender=rank receiver, no MPI communications
        ! write(*,*) 'rank ',rank,f2cadd(i,1),' int copy to ',f2cadd(i,2)
        il=max(cstart(1)+1,cstartpsi(1)+1)
        iu=min(cstart(1)+1+spx-1,cstartpsi(1)+1+spxpsi-1)
        ! copy 1:ny/2+1
        jl=max(cstart(3)+1,cstartpsi(3)+1)
        ju=min(cstart(3)+1+spy-1,cstartpsi(3)+1+spypsi-1,ny/2+1)
        !$acc kernels
        if(jl.le.ju) varc(il-cstart(1):iu-cstart(1),1:nz,jl-cstart(3):ju-cstart(3),:)= &
          & varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,jl-cstartpsi(3):ju-cstartpsi(3),:)
        !$acc end kernels
        ! copy ny/2+2:ny
        ! minimum distance from top, maximum distance from top
        jl=min(ny-cstart(3),npsiy-cstartpsi(3),ny/2-1)
        ju=max(ny-(cstart(3)+1+spy-1),npsiy-(cstartpsi(3)+1+spypsi-1))
        !$acc kernels
        if(jl.ge.ju) varc(il-cstart(1):iu-cstart(1),1:nz,ny-jl-cstart(3)+1:ny-ju-cstart(3),:)= &
          & varf(il-cstartpsi(1):iu-cstartpsi(1),1:nz,npsiy-jl-cstartpsi(3)+1:npsiy-ju-cstartpsi(3),:)
        !$acc end kernels
      else
        !$acc kernels
        bufs=0.0d0
        bufs(1:spxpsi,1:npsiz,1:spypsi,1:2)=varf
        !$acc end kernels
        ! tag=1000000*f2cadd(i,1)+f2cadd(i,2)
        tag=15
        if(rank.eq.f2cadd(i,1)) then
          ! send
          !$acc host_data use_device(bufs)
          call mpi_ssend(bufs,numel,MPI_double_precision,f2cadd(i,2),tag,flow_comm,ierr)
          !$acc end host_data
        elseif(rank.eq.f2cadd(i,2)) then
          ! receive, sender=receiver already accounted for in internal copy
          !$acc host_data use_device(bufr)
          call mpi_recv(bufr,numel,MPI_double_precision,f2cadd(i,1),tag,flow_comm,mpi_status_ignore,ierr)
          !$acc end host_data
          ! update fine grid variable
          ii=mod(f2cadd(i,1),nycpu)+1
          jj=int((f2cadd(i,1)-(ii-1))/nycpu)+1
          il=max(cstart(1)+1,fg_size(ii,jj,4))
          iu=min(cstart(1)+1+spx-1,fg_size(ii,jj,4)+fg_size(ii,jj,2)-1)
          ! copy 1:ny/2+1
          jl=max(cstart(3)+1,fg_size(ii,jj,5))
          ju=min(cstart(3)+1+spy-1,fg_size(ii,jj,5)+fg_size(ii,jj,3)-1,ny/2+1)
          !$acc kernels
          if(jl.le.ju) varc(il-cstart(1):iu-cstart(1),1:nz,jl-cstart(3):ju-cstart(3),:)= &
            & bufr(il-fg_size(ii,jj,4)+1:iu-fg_size(ii,jj,4)+1,1:nz,jl-fg_size(ii,jj,5)+1:ju-fg_size(ii,jj,5)+1,:)
          !$acc end kernels
          ! copy ny/2+2:ny
          ! minimum distance from top, maximum distance from top
          jl=min(ny-cstart(3),npsiy-fg_size(ii,jj,5)+1,ny/2-1)
          ju=max(ny-(cstart(3)+1+spy-1),npsiy-(fg_size(ii,jj,5)+fg_size(ii,jj,3)-1))
          !$acc kernels
          if(jl.ge.ju) varc(il-cstart(1):iu-cstart(1),1:nz,ny-jl-cstart(3)+1:ny-ju-cstart(3),:)= &
            & bufr(il-fg_size(ii,jj,4)+1:iu-fg_size(ii,jj,4)+1,1:nz,npsiy-jl+1-fg_size(ii,jj,5)+1:npsiy-ju+1-fg_size(ii,jj,5),:)
          !$acc end kernels
        endif
      endif
    endif
  enddo
  deallocate(bufs)
  deallocate(bufr)
#endif


! renormalize
!$acc kernels
varc=varc/dble(exp_x*exp_y)
!$acc end kernels

return
end
