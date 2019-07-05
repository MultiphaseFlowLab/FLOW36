subroutine initialize_stats

use commondata
use stats
use mpi
use par_size
use grid
use sim_par
use wavenumber

double precision :: st_glob(nz,12),pm(nz),prms(nz),budget_data(nz,5)
! u,v and w power spectra at three different z^+ locations
double precision :: stream_ps(nx/2+1,9), span_ps(ny/2+1,9)

integer :: i

flowiter=0

#define mean meanflag
#define budget budgetflag
#define spectra spectraflag

call mpi_cart_sub(cart_comm,[.false.,.true.],plane_comm,ierr)
call mpi_cart_sub(cart_comm,[.true.,.false.],col_comm,ierr)

if(stat_start.eq.nstart)then
#if mean == 1
 call mean_calc(st_glob)
#endif
#if budget == 1
 call budget_calc(pm,prms,budget_data)
#endif
#if spectra == 1
 call power_spectra(stream_ps,span_ps)
#endif
else
! when starting statistic calculation from stat_start skip the initial step: not included in the statistics calculation
#if mean == 1
 st_glob=0.0d0
#endif
#if budget == 1
 pm=0.0d0
 prms=0.0d0
 budget_data=0.0d0
#endif
#if spectra == 1
 stream_ps=0.0d0
 span_ps=0.0d0
#endif
 flowiter=-1
endif

if((rank.eq.0).and.(restart.eq.0))then

#if mean == 1
 open(66,file=trim(folder)//'/stats.dat',form='formatted',status='new')
 write(66,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(66,'(13(a12,2x))') 'z','u mean','v mean','w mean','rms u','rms v','rms w', &
   &   'skw u','skw v','skw w','flt u','flt v','flt w'
 write(66,*)
 do i=1,nz
  write(66,'(13(ES12.5,2x))') (z(i)+1.0d0)*re,st_glob(i,:)
 enddo
 close(66,status='keep')
#endif
#if budget == 1
 open(67,file=trim(folder)//'/budget.dat',form='formatted',status='new')
 write(67,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(67,'(8(a12,2x))') 'z','p_m','p_rms','turb prod','turb transp','pressure b','visc diss','pseudo diss'
 write(67,*)
 do i=1,nz
  write(67,'(8(ES12.5,2x))') (z(i)+1.0d0)*re,pm(i),prms(i),budget_data(i,:)
 enddo
 close(67,status='keep')
#endif
#if spectra == 1
 open(68,file=trim(folder)//'/power_xspectra.dat',form='formatted',status='new')
 write(68,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(68,'(a16,2x,9(a10,f6.0,2x))') 'kx','stream u',5.0,'stream v',5.0,'stream w',5.0, &
 &                                      'stream u',15.0,'stream v',15.0,'stream w',15.0, &
 &                                      'stream u',re,'stream v',re,'stream w',re
 write(68,*)
 do i=1,nx/2+1
  write(68,'(10(ES16.5,2x))') kx(i),stream_ps(i,:)
 enddo

 close(68,status='keep')

 open(69,file=trim(folder)//'/power_yspectra.dat',form='formatted',status='new')
 write(69,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(69,'(a16,2x,9(a10,f6.0,2x))') 'ky','span u',5.0,'span v',5.0,'span w',5.0, &
 &                                      'span u',15.0,'span v',15.0,'span w',15.0, &
 &                                      'span u',re,'span v',re,'span w',re
 write(69,*)
 do i=1,ny/2+1
  write(69,'(10(ES16.5,2x))') ky(i),span_ps(i,:)
 enddo

 close(69,status='keep')
#endif

elseif((rank.eq.0).and.(restart.eq.1)) then

#if mean == 1
 open(66,file=trim(folder)//'/stats.dat',form='formatted',status='old')
 read(66,'(22x,i8)') flowiter
 close(66,status='keep')
#endif
#if budget == 1
 open(67,file=trim(folder)//'/budget.dat',form='formatted',status='old')
 read(67,'(22x,i8)') flowiter
 close(67,status='keep')
#endif
#if spectra == 1
 open(68,file=trim(folder)//'/power_xspectra.dat',form='formatted',status='old')
 read(68,'(22x,i8)') flowiter
 close(68,status='keep')
#endif
 flowiter=flowiter-1
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine del_old_stats

use commondata
use stats
use mpi

! #define mean meanflag
! #define budget budgetflag
! #define spectra spectraflag

! ! delete old stat file
! if(rank.eq.0)then
! #if mean == 1
!  open(66,file=trim(folder)//'/stats_old.dat',form='formatted',status='old')
!  close(66,status='delete')
! #endif
! #if budget == 1
!  open(67,file=trim(folder)//'/budget_old.dat',form='formatted',status='old')
!  close(67,status='delete')
! #endif
! #if spectra == 1
!  open(68,file=trim(folder)//'/power_xspectra_old.dat',form='formatted',status='old')
!  close(68,status='delete')
!
!  open(69,file=trim(folder)//'/power_yspectra_old.dat',form='formatted',status='old')
!  close(69,status='delete')
! #endif
! endif

call mpi_comm_free(plane_comm,ierr)
call mpi_comm_free(col_comm,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine statistics

use commondata
use stats
use velocity
use grid
use sim_par
use wavenumber

double precision :: st_glob(nz,12),st_glob_old(nz,12),tmp,pm(nz),p_old(nz),prms(nz),prms_old(nz)
double precision :: budget_data(nz,5),budget_old(nz,5)
! u,v and w power spectra at three different z^+ locations
double precision :: stream_ps(nx/2+1,9), span_ps(ny/2+1,9),stream_old(nx/2+1,9), span_old(ny/2+1,9)

integer :: i

flowiter=flowiter+1

#define mean meanflag
#define budget budgetflag
#define spectra spectraflag

call spectral_to_phys(uc,u,0)
call spectral_to_phys(vc,v,0)
call spectral_to_phys(wc,w,0)


#if mean == 1
 call mean_calc(st_glob)
#endif
#if budget == 1
 call budget_calc(pm,prms,budget_data)
#endif
#if spectra == 1
 call power_spectra(stream_ps,span_ps)
#endif



if(rank.eq.0)then
! read old stats files
#if mean == 1
 open(66,file=trim(folder)//'/stats.dat',form='formatted',status='old')
 !read from old stats file
 read(66,*)
 read(66,*)
 read(66,*)
 do i=1,nz
  read(66,'(13(ES12.5,2x))') tmp,st_glob_old(i,:)
 enddo
 close(66,status='delete')
#endif
#if budget == 1
 open(76,file=trim(folder)//'/budget.dat',form='formatted',status='old')
 !read from old stats file
 read(76,*)
 read(76,*)
 read(76,*)
 do i=1,nz
  read(76,'(8(ES12.5,2x))') tmp,p_old(i),prms_old(i),budget_old(i,:)
 enddo
 close(76,status='delete')
#endif
#if spectra == 1
 open(86,file=trim(folder)//'/power_xspectra.dat',form='formatted',status='old')
 !read from old stats file
 read(86,*)
 read(86,*)
 read(86,*)
 do i=1,nx/2+1
  read(86,'(10(ES16.5,2x))') tmp,stream_old(i,:)
 enddo
 close(86,status='delete')

 open(87,file=trim(folder)//'/power_yspectra.dat',form='formatted',status='old')
 !read from old stats file
 read(87,*)
 read(87,*)
 read(87,*)
 do i=1,ny/2+1
  read(87,'(10(ES16.5,2x))') tmp,span_old(i,:)
 enddo
 close(87,status='delete')
#endif

! mean over time
#if mean == 1
 st_glob=(dble(flowiter)*st_glob_old+st_glob)/dble(flowiter+1)
 open(67,file=trim(folder)//'/stats.dat',form='formatted',status='new')
 write(67,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(67,'(13(a12,2x))') 'z','u mean','v mean','w mean','rms u','rms v','rms w', &
   &   'skw u','skw v','skw w','flt u','flt v','flt w'
 write(67,*)
 do i=1,nz
  write(67,'(13(ES12.5,2x))') (z(i)+1.0d0)*re,st_glob(i,:)
 enddo
 close(67,status='keep')
#endif
#if budget == 1
 pm=(dble(flowiter)*p_old+pm)/dble(flowiter+1)
 prms=(dble(flowiter)*prms_old+prms)/dble(flowiter+1)
 budget_data=(dble(flowiter)*budget_old+budget_data)/dble(flowiter+1)
 open(77,file=trim(folder)//'/budget.dat',form='formatted',status='new')
 write(77,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(77,'(8(a12,2x))') 'z','p_m','p_rms','turb prod','turb transp','pressure b','visc diss','pseudo diss'
 write(77,*)
 do i=1,nz
  write(77,'(8(ES12.5,2x))') (z(i)+1.0d0)*re,pm(i),prms(i),budget_data(i,:)
 enddo
 close(77,status='keep')
#endif
#if spectra == 1
 stream_ps=(dble(flowiter)*stream_old+stream_ps)/dble(flowiter+1)
 span_ps=(dble(flowiter)*span_old+span_ps)/dble(flowiter+1)
 open(87,file=trim(folder)//'/power_xspectra.dat',form='formatted',status='new')
 write(87,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(87,'(a16,2x,9(a10,f6.0,2x))') 'kx','stream u',5.0,'stream v',5.0,'stream w',5.0, &
 &                                      'stream u',15.0,'stream v',15.0,'stream w',15.0, &
 &                                      'stream u',re,'stream v',re,'stream w',re
 write(87,*)
 do i=1,nx/2+1
  write(87,'(10(ES16.5,2x))') kx(i),stream_ps(i,:)
 enddo

 close(87,status='keep')

 open(88,file=trim(folder)//'/power_yspectra.dat',form='formatted',status='new')
 write(88,'(a,i8,a,3(i6),a)') 'Statistics gathered on',flowiter+1,' flow fields, on a',nx,ny,nz,' grid (nx,ny,nz)'
 write(88,'(a16,2x,9(a10,f6.0,2x))') 'ky','span u',5.0,'span v',5.0,'span w',5.0, &
 &                                      'span u',15.0,'span v',15.0,'span w',15.0, &
 &                                      'span u',re,'span v',re,'span w',re
 write(88,*)
 do i=1,ny/2+1
  write(88,'(10(ES16.5,2x))') ky(i),span_ps(i,:)
 enddo

 close(88,status='keep')
#endif
endif



return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mean_calc(st_glob)

use commondata
use par_size
use sim_par
use velocity
use mpi
use stats

! st_loc[t](:,1:3) : u_mean, v_mean, w_mean
! st_loc[t](:,4:6) : rms_u, rms_v, rms_w
! st_loc[t](:,7,9) : skw_u, skw_v, skw_w
! st_loc[t](:,10:12) : flt_u, flt_v, flt_w

! st_glob[t](:,1:3) : u_mean, v_mean, w_mean
! st_glob[t](:,4:6) : rms_u, rms_v, rms_w
! st_glob[t](:,7,9) : skw_u, skw_v, skw_w
! st_glob[t](:,10:12) : flt_u, flt_v, flt_w

double precision :: st_loc(fpz,12), st_loct(fpz,12), st_glob(nz,12), st_globt(nz,12)

integer :: i,j

st_loc=0.0d0
st_globt=0.0d0

do j=1,fpy
  do i=1,nx
    st_loc(:,1)=st_loc(:,1)+u(i,:,j)
    st_loc(:,2)=st_loc(:,2)+v(i,:,j)
    st_loc(:,3)=st_loc(:,3)+w(i,:,j)
  enddo
enddo

! exchange data among ranks in same x-y plane
call mpi_allreduce(st_loc,st_loct,fpz*12,mpi_double_precision,mpi_sum,plane_comm,ierr)

! each rank has the instantaneous mean velocity in its x-y plane
st_loct(:,1:3)=st_loct(:,1:3)/dble(nx*ny)

! calculate rms, skewness, flatness
do j=1,fpy
  do i=1,nx
    st_loct(:,4)=st_loct(:,4)+(u(i,:,j)-st_loct(:,1))**2
    st_loct(:,5)=st_loct(:,5)+(v(i,:,j)-st_loct(:,2))**2
    st_loct(:,6)=st_loct(:,6)+(w(i,:,j)-st_loct(:,3))**2
    st_loct(:,7)=st_loct(:,7)+(u(i,:,j)-st_loct(:,1))**3
    st_loct(:,8)=st_loct(:,8)+(v(i,:,j)-st_loct(:,2))**3
    st_loct(:,9)=st_loct(:,9)+(w(i,:,j)-st_loct(:,3))**3
    st_loct(:,10)=st_loct(:,10)+(u(i,:,j)-st_loct(:,1))**4
    st_loct(:,11)=st_loct(:,11)+(v(i,:,j)-st_loct(:,2))**4
    st_loct(:,12)=st_loct(:,12)+(w(i,:,j)-st_loct(:,3))**4
  enddo
enddo

call mpi_allreduce(st_loct,st_loc,fpz*12,mpi_double_precision,mpi_sum,plane_comm,ierr)


st_loc=st_loc/dble(nx*ny)
st_loc(:,1:3)=st_loct(:,1:3)
st_globt(fstart(2)+1:fstart(2)+fpz,:)=st_loc



if(mod(rank,nycpu).eq.0)then
 call mpi_reduce(st_globt,st_glob,nz*12,mpi_double_precision,mpi_sum,0,col_comm,ierr)
 do i=1,nz
   st_glob(i,4)=(st_glob(i,4))**0.5
   st_glob(i,5)=(st_glob(i,5))**0.5
   st_glob(i,6)=(st_glob(i,6))**0.5
   st_glob(i,7)=st_glob(i,7)/(st_glob(i,4))**3
   st_glob(i,8)=st_glob(i,8)/(st_glob(i,5))**3
   st_glob(i,9)=st_glob(i,9)/(st_glob(i,6))**3
   st_glob(i,10)=st_glob(i,10)/(st_glob(i,4))**4
   st_glob(i,11)=st_glob(i,11)/(st_glob(i,5))**4
   st_glob(i,12)=st_glob(i,12)/(st_glob(i,6))**4
 enddo
 st_glob(1,7:12)=0.0d0
 st_glob(nz,7:12)=0.0d0
endif


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine budget_calc(pm,prms,budget_data)

use commondata
use par_size
use sim_par
use velocity
use mpi
use stats
use wavenumber

double precision :: f(spx,nz,spy,2),press(nx,fpz,fpy),ptm(nz),pm(nz),ptrms(nz),prms(nz)
double precision, dimension(spx,nz,spy,2) :: s1,s2,s3
double precision :: ftmp(nx,fpz,fpy),f2tmp(nx,fpz,fpy),stmp(spx,nz,spy,2),s2tmp(spx,nz,spy,2)
double precision :: rhs(2)
double precision, dimension(nz,5) :: budget_data
double precision, dimension(nz) :: um,vm,wm,utm,vtm,wtm,tmp1,tmp2

integer :: i,j,k

! pressure solver
call sterm_pressure(s1,s2,s3)
call dz(s3,f)

do j=1,spy
  do k=1,nz
    do i=1,spx
      f(i,k,j,1)=f(i,k,j,1)-kx(cstart(1)+i)*s1(i,k,j,2)-ky(cstart(3)+j)*s2(i,k,j,2)
      f(i,k,j,2)=f(i,k,j,2)+kx(cstart(1)+i)*s1(i,k,j,1)+ky(cstart(3)+j)*s2(i,k,j,1)
    enddo
  enddo
enddo


! set boundary conditions
call dz(wc,s2tmp)
call dz(s2tmp,stmp)

! solve Helmholtz equation for pressure
do j=1,spy
  do i=1,spx
    rhs(1)=0.5d0*stmp(i,1,j,1)
    rhs(2)=0.5d0*stmp(i,1,j,1)
    do k=2,nz
      rhs(1)=rhs(1)+(zp(1))**(k-1)*stmp(i,k,j,1)
      rhs(2)=rhs(2)+stmp(i,k,j,1)
    enddo
    rhs=rhs/re
    pm=f(i,:,j,1)
    call helmholtz_rred(pm,k2(i+cstart(1),j+cstart(3)),[0.0d0,0.0d0],[1.0d0,1.0d0],rhs,zp)
    f(i,:,j,1)=pm

    rhs(1)=0.5d0*stmp(i,1,j,2)
    rhs(2)=0.5d0*stmp(i,1,j,2)
    do k=2,nz
      rhs(1)=rhs(1)+(zp(1))**(k-1)*stmp(i,k,j,2)
      rhs(2)=rhs(2)+stmp(i,k,j,2)
    enddo
    rhs=rhs/re
    pm=f(i,:,j,2)
    call helmholtz_rred(pm,k2(i+cstart(1),j+cstart(3)),[0.0d0,0.0d0],[1.0d0,1.0d0],rhs,zp)
    f(i,:,j,2)=pm
  enddo
enddo

! solve for k2=0
call spectral_to_phys(wc,w,1)
pm=0.0d0
do k=1,fpz
  do j=1,fpy
    do i=1,nx
      pm(fstart(2)+k)=pm(fstart(2)+k)+w(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

call mpi_reduce(pm,ptm,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)

if(rank.eq.0)then
  ptm=ptm/dble(nx*ny)
  call dctz_fwd(ptm,ptm,1,nz,1,1)
  f(1,:,1,1)=-dble(nx*ny)*ptm
  f(1,:,1,2)=0.0d0
endif

call spectral_to_phys(f,press,0)

! take mean in x,y
pm=0.0d0
do j=1,fpy
  do i=1,nx
    pm(fstart(2)+1:fstart(2)+fpz)=pm(fstart(2)+1:fstart(2)+fpz)+press(i,:,j)
  enddo
enddo

call mpi_allreduce(pm,ptm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)

pm=ptm/dble(nx*ny)

! calculate pressure rms
ptrms=0.0d0
do k=1,fpz
  do j=1,fpy
    do i=1,nx
      ptrms(fstart(2)+k)=ptrms(fstart(2)+k)+(press(i,k,j)-pm(fstart(2)+k))**2
    enddo
  enddo
enddo

call mpi_reduce(ptrms,prms,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)

prms=prms/dble(nx*ny)

do k=1,nz
  prms(k)=(prms(k))**0.5
enddo

! energy budget calculation
budget_data=0.0d0

! turbulent production
call spectral_to_phys(uc,u,1)
call spectral_to_phys(vc,v,1)
call spectral_to_phys(wc,w,1)

call dz(uc,stmp)
call spectral_to_phys(stmp,ftmp,1)

utm=0.0d0
vtm=0.0d0
wtm=0.0d0
ptrms=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      utm(fstart(2)+k)=utm(fstart(2)+k)+u(i,k,j)
      vtm(fstart(2)+k)=vtm(fstart(2)+k)+v(i,k,j)
      wtm(fstart(2)+k)=wtm(fstart(2)+k)+w(i,k,j)
      ptrms(fstart(2)+k)=ptrms(fstart(2)+k)+ftmp(i,k,j)
    enddo
  enddo
enddo
utm=utm/dble(nx*ny)
vtm=vtm/dble(nx*ny)
wtm=wtm/dble(nx*ny)
ptrms=ptrms/dble(nx*ny)

call mpi_allreduce(utm,um,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(vtm,vm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(wtm,wm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(ptrms,ptm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)

tmp1=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+ptm(fstart(2)+k)*(u(i,k,j)-um(fstart(2)+k))*(w(i,k,j)-wm(fstart(2)+k))
    enddo
  enddo
enddo

call mpi_reduce(tmp1,tmp2,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)
tmp2=tmp2/dble(nx*ny)

budget_data(1:nz,1)=tmp2

! turbulent transport
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=(u(i,k,j)-um(fstart(2)+k))**2*(w(i,k,j)-wm(fstart(2)+k))+ &
 &                (v(i,k,j)-vm(fstart(2)+k))**2*(w(i,k,j)-wm(fstart(2)+k))+ &
 &                (w(i,k,j)-wm(fstart(2)+k))**3
    enddo
  enddo
enddo
ftmp=ftmp*0.5d0

call phys_to_spectral(ftmp,stmp,1)
call dz(stmp,s2tmp)
call spectral_to_phys(s2tmp,ftmp,1)

tmp1=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+ftmp(i,k,j)
    enddo
  enddo
enddo

call mpi_reduce(tmp1,tmp2,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)
tmp2=tmp2/dble(nx*ny)

budget_data(1:nz,2)=tmp2

! pressure
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=-(press(i,k,j)-pm(fstart(2)+k))*(w(i,k,j)-wm(fstart(2)+k))
    enddo
  enddo
enddo

call phys_to_spectral(ftmp,stmp,1)
call dz(stmp,s2tmp)
call spectral_to_phys(s2tmp,ftmp,1)

tmp1=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+ftmp(i,k,j)
    enddo
  enddo
enddo

call mpi_reduce(tmp1,tmp2,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)
tmp2=tmp2/dble(nx*ny)

budget_data(1:nz,3)=tmp2

! kinetic energy dissipation
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=0.5d0*((u(i,k,j)-um(fstart(2)+k))**2+(v(i,k,j)-vm(fstart(2)+k))**2+ &
 &                (w(i,k,j)-wm(fstart(2)+k))**2)
    enddo
  enddo
enddo

call phys_to_spectral(ftmp,stmp,1)
call dz(stmp,s2tmp)
call dz(s2tmp,stmp)
call spectral_to_phys(stmp,ftmp,1)

tmp1=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+ftmp(i,k,j)/re
    enddo
  enddo
enddo

call mpi_reduce(tmp1,tmp2,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)
tmp2=tmp2/dble(nx*ny)

budget_data(1:nz,4)=tmp2

! pseudo-dissipation
! u prime
! dx
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=u(i,k,j)-um(fstart(2)+k)
    enddo
  enddo
enddo

call phys_to_spectral(ftmp,stmp,1)

do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-kx(cstart(1)+i)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=kx(cstart(1)+i)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

tmp1=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dy
do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-ky(cstart(3)+j)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=ky(cstart(3)+j)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dz
call dz(stmp,s2tmp)

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! v prime
! dx
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=v(i,k,j)-vm(fstart(2)+k)
    enddo
  enddo
enddo

call phys_to_spectral(ftmp,stmp,1)

do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-kx(cstart(1)+i)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=kx(cstart(1)+i)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dy
do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-ky(cstart(3)+j)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=ky(cstart(3)+j)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dz
call dz(stmp,s2tmp)

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! w prime
! dx
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ftmp(i,k,j)=w(i,k,j)-wm(fstart(2)+k)
    enddo
  enddo
enddo

call phys_to_spectral(ftmp,stmp,1)

do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-kx(cstart(1)+i)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=kx(cstart(1)+i)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dy
do j=1,spy
  do k=1,nz
    do i=1,spx
      s2tmp(i,k,j,1)=-ky(cstart(3)+j)*stmp(i,k,j,2)
      s2tmp(i,k,j,2)=ky(cstart(3)+j)*stmp(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! dz
call dz(stmp,s2tmp)

call spectral_to_phys(s2tmp,f2tmp,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      tmp1(fstart(2)+k)=tmp1(fstart(2)+k)+f2tmp(i,k,j)*f2tmp(i,k,j)
    enddo
  enddo
enddo

! gather terms
call mpi_reduce(tmp1,tmp2,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)
tmp2=-tmp2/(re*dble(nx*ny))

budget_data(1:nz,5)=tmp2

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sterm_pressure(s1,s2,s3)

use commondata
use velocity
use par_size
use wavenumber
use phase_field

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision, allocatable, dimension(:,:,:) :: uu,uv,uw,vv,vw,ww
double precision, allocatable, dimension(:,:,:,:) :: uuc,uvc,uwc,vvc,vwc,wwc
double precision, allocatable, dimension(:,:,:,:) :: uucx,uvcy,uwcz, uvcx,vvcy,vwcz, uwcx,vwcy,wwcz


integer :: indx,indy
integer :: i,j,k

! transform variables back to physical space and perform dealiasing
call spectral_to_phys(uc,u,1)
call spectral_to_phys(vc,v,1)
call spectral_to_phys(wc,w,1)

allocate(uu(nx,fpz,fpy))
allocate(uv(nx,fpz,fpy))
allocate(uw(nx,fpz,fpy))
allocate(vv(nx,fpz,fpy))
allocate(vw(nx,fpz,fpy))
allocate(ww(nx,fpz,fpy))

allocate(uuc(spx,nz,spy,2))
allocate(uvc(spx,nz,spy,2))
allocate(uwc(spx,nz,spy,2))
allocate(vvc(spx,nz,spy,2))
allocate(vwc(spx,nz,spy,2))
allocate(wwc(spx,nz,spy,2))


! form products uu,uv,uw, ...
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      uu(i,k,j)=u(i,k,j)*u(i,k,j)
      uv(i,k,j)=u(i,k,j)*v(i,k,j)
      uw(i,k,j)=u(i,k,j)*w(i,k,j)
      vv(i,k,j)=v(i,k,j)*v(i,k,j)
      vw(i,k,j)=v(i,k,j)*w(i,k,j)
      ww(i,k,j)=w(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo


! transform products uu,uv, ...  in spectral space
call phys_to_spectral(uu,uuc,1)
call phys_to_spectral(uv,uvc,1)
call phys_to_spectral(uw,uwc,1)
call phys_to_spectral(vv,vvc,1)
call phys_to_spectral(vw,vwc,1)
call phys_to_spectral(ww,wwc,1)

deallocate(uu)
deallocate(uv)
deallocate(uw)
deallocate(vv)
deallocate(vw)
deallocate(ww)


! take derivatives
allocate(uucx(spx,nz,spy,2))
allocate(uvcy(spx,nz,spy,2))
allocate(uwcz(spx,nz,spy,2))
allocate(uvcx(spx,nz,spy,2))
allocate(vvcy(spx,nz,spy,2))
allocate(vwcz(spx,nz,spy,2))
allocate(uwcx(spx,nz,spy,2))
allocate(vwcy(spx,nz,spy,2))
allocate(wwcz(spx,nz,spy,2))

! x derivatives
indx=cstart(1)
do i=1,spx
  uucx(i,:,:,1)=-uuc(i,:,:,2)*kx(indx+i)
  uucx(i,:,:,2)=uuc(i,:,:,1)*kx(indx+i)
  uvcx(i,:,:,1)=-uvc(i,:,:,2)*kx(indx+i)
  uvcx(i,:,:,2)=uvc(i,:,:,1)*kx(indx+i)
  uwcx(i,:,:,1)=-uwc(i,:,:,2)*kx(indx+i)
  uwcx(i,:,:,2)=uwc(i,:,:,1)*kx(indx+i)
enddo


! y derivatives
indy=cstart(3)
do j=1,spy
  uvcy(:,:,j,1)=-uvc(:,:,j,2)*ky(indy+j)
  uvcy(:,:,j,2)=uvc(:,:,j,1)*ky(indy+j)
  vvcy(:,:,j,1)=-vvc(:,:,j,2)*ky(indy+j)
  vvcy(:,:,j,2)=vvc(:,:,j,1)*ky(indy+j)
  vwcy(:,:,j,1)=-vwc(:,:,j,2)*ky(indy+j)
  vwcy(:,:,j,2)=vwc(:,:,j,1)*ky(indy+j)
enddo


! z derivatives
call dz(uwc,uwcz)
call dz(vwc,vwcz)
call dz(wwc,wwcz)



deallocate(uuc)
deallocate(uvc)
deallocate(uwc)
deallocate(vvc)
deallocate(vwc)
deallocate(wwc)


! form non linear terms
s1=-(uucx+uvcy+uwcz)
s2=-(uvcx+vvcy+vwcz)
s3=-(uwcx+vwcy+wwcz)


deallocate(uucx)
deallocate(uvcy)
deallocate(uwcz)
deallocate(uvcx)
deallocate(vvcy)
deallocate(vwcz)
deallocate(uwcx)
deallocate(vwcy)
deallocate(wwcz)


!!! conditional compilation, add phase variable contribution for non-matched densities
!!#define match_dens matched_density
!!#define phiflag phicompflag
!!#if phiflag == 1
!!#if match_dens == 0
!!! case for non-matched densities

!!allocate(uu(nx,fpz,fpy))
!!allocate(vv(nx,fpz,fpy))
!!allocate(ww(nx,fpz,fpy))

!!call spectral_to_phys(s1,uu,1)
!!call spectral_to_phys(s2,vv,1)
!!call spectral_to_phys(s3,ww,1)

!!call spectral_to_phys(phic,phi,1)

!!!
!!do j=1,fpy
!!  do k=1,fpz
!!    do i=1,nx
!!      uu(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phi(i,k,j)+1.0d0))*uu(i,k,j)
!!      vv(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phi(i,k,j)+1.0d0))*vv(i,k,j)
!!      ww(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phi(i,k,j)+1.0d0))*ww(i,k,j)
!!    enddo
!!  enddo
!!enddo

!!call phys_to_spectral(uu,s1,1)
!!call phys_to_spectral(vv,s2,1)
!!call phys_to_spectral(ww,s3,1)

!!deallocate(uu)
!!deallocate(vv)
!!deallocate(ww)

!!#endif
!!#endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine power_spectra(stream_ps,span_ps)

use commondata
use par_size
use sim_par
use velocity
use mpi
use stats
use grid
use fftw3

double precision, dimension(nz) :: um,vm,wm,utm,vtm,wtm
double precision, dimension(nx,ny) :: up,vp,wp,utp,vtp,wtp
double precision :: zplus
! u,v and w power spectra at three different z^+ locations
double precision :: stream_ps(nx/2+1,9), span_ps(ny/2+1,9)

complex(c_double_complex), dimension(nx/2+1,ny) :: upcx,vpcx,wpcx
complex(c_double_complex), dimension(nx,ny/2+1) :: upcy,vpcy,wpcy

integer :: i,j,k
integer :: zind, flag
integer :: pl_comm,red_comm

type(c_ptr) :: plan_x,plan_y

call spectral_to_phys(uc,u,1)
call spectral_to_phys(vc,v,1)
call spectral_to_phys(wc,w,1)

utm=0.0d0
vtm=0.0d0
wtm=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      utm(fstart(2)+k)=utm(fstart(2)+k)+u(i,k,j)
      vtm(fstart(2)+k)=vtm(fstart(2)+k)+v(i,k,j)
      wtm(fstart(2)+k)=wtm(fstart(2)+k)+w(i,k,j)
    enddo
  enddo
enddo

call mpi_allreduce(utm,um,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(vtm,vm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(wtm,wm,nz,mpi_double_precision,mpi_sum,flow_comm,ierr)

um=um/dble(nx*ny)
vm=vm/dble(nx*ny)
wm=wm/dble(nx*ny)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! power spectra at z^+=5
zplus=5.0d0
! find which ranks contain same z^+
zind=minloc(dabs((z+1.0d0)*re-zplus),1)

if((fstart(2)+1.le.zind).and.(fstart(2)+fpz.ge.zind))then
  flag=1
else
 flag=0
endif

! create subcommunicator among ranks that contains z^+
call mpi_comm_split(flow_comm,flag,0,pl_comm,ierr)

utp=0.0d0
vtp=0.0d0
wtp=0.0d0
if(flag.eq.1)then
  do j=1,fpy
    do i=1,nx
        utp(i,fstart(3)+j)=u(i,zind-fstart(2),j)-um(zind)
        vtp(i,fstart(3)+j)=v(i,zind-fstart(2),j)-vm(zind)
        wtp(i,fstart(3)+j)=w(i,zind-fstart(2),j)-wm(zind)
    enddo
  enddo

! reduce to rank floor(rank/nycpu)*nycpu
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)

! set flag to zero, only reduce rank keeps flag=1
  if(rank.ne.floor(real(rank)/real(nycpu))*nycpu)then
    flag=0
  else
    utp=up
    vtp=vp
    wtp=wp
  endif
endif

call mpi_comm_free(pl_comm,ierr)

! if rank 0 is not in the new communicator it must be included
if(rank.eq.0) flag=1
call mpi_comm_split(flow_comm,flag,0,red_comm,ierr)

! reduce velocity fluctuations to rank 0
if(flag.eq.1)then
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
endif

call mpi_comm_free(red_comm,ierr)

! calculate power spectra, done only from rank 0
if(rank.eq.0)then
! create plan for fft x and fft y
  plan_x=fftw_plan_many_dft_r2c(1,[nx],ny,up,[nx,ny],1,nx,upcx,[nx/2+1,ny],1,nx/2+1,FFTW_ESTIMATE)
  plan_y=fftw_plan_many_dft_r2c(1,[ny],nx,up,[nx,ny],nx,1,upcy,[nx,ny/2+1],nx,1,FFTW_ESTIMATE)

! perform Fourier transform
  call fftw_execute_dft_r2c(plan_x,up,upcx)
  call fftw_execute_dft_r2c(plan_y,up,upcy)
  call fftw_execute_dft_r2c(plan_x,vp,vpcx)
  call fftw_execute_dft_r2c(plan_y,vp,vpcy)
  call fftw_execute_dft_r2c(plan_x,wp,wpcx)
  call fftw_execute_dft_r2c(plan_y,wp,wpcy)

! calculate power spectra
  stream_ps=0.0d0
  span_ps=0.0d0

  do j=1,ny
    do i=1,nx/2+1
      stream_ps(i,1)=stream_ps(i,1)+(dble(upcx(i,j)))**2+(aimag(upcx(i,j)))**2
      stream_ps(i,2)=stream_ps(i,2)+(dble(vpcx(i,j)))**2+(aimag(vpcx(i,j)))**2
      stream_ps(i,3)=stream_ps(i,3)+(dble(wpcx(i,j)))**2+(aimag(wpcx(i,j)))**2
    enddo
  enddo

  do j=1,ny/2+1
    do i=1,nx
      span_ps(j,1)=span_ps(j,1)+(dble(upcy(i,j)))**2+(aimag(upcy(i,j)))**2
      span_ps(j,2)=span_ps(j,2)+(dble(vpcy(i,j)))**2+(aimag(vpcy(i,j)))**2
      span_ps(j,3)=span_ps(j,3)+(dble(wpcy(i,j)))**2+(aimag(wpcy(i,j)))**2
    enddo
  enddo

  call fftw_destroy_plan(plan_x)
  call fftw_destroy_plan(plan_y)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! power spectra at z^+=15
zplus=15.0d0
! find which ranks contain same z^+
zind=minloc(dabs((z+1.0d0)*re-zplus),1)

if((fstart(2)+1.le.zind).and.(fstart(2)+fpz.ge.zind))then
  flag=1
else
 flag=0
endif

! create subcommunicator among ranks that contains z^+
call mpi_comm_split(flow_comm,flag,0,pl_comm,ierr)

utp=0.0d0
vtp=0.0d0
wtp=0.0d0
if(flag.eq.1)then
  do j=1,fpy
    do i=1,nx
        utp(i,fstart(3)+j)=u(i,zind-fstart(2),j)-um(zind)
        vtp(i,fstart(3)+j)=v(i,zind-fstart(2),j)-vm(zind)
        wtp(i,fstart(3)+j)=w(i,zind-fstart(2),j)-wm(zind)
    enddo
  enddo

! reduce to rank floor(rank/nycpu)*nycpu
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)

! set flag to zero, only reduce rank keeps flag=1
  if(rank.ne.floor(real(rank)/real(nycpu))*nycpu)then
    flag=0
  else
    utp=up
    vtp=vp
    wtp=wp
  endif
endif

call mpi_comm_free(pl_comm,ierr)

! if rank 0 is not in the new communicator it must be included
if(rank.eq.0) flag=1
call mpi_comm_split(flow_comm,flag,0,red_comm,ierr)

! reduce velocity fluctuations to rank 0
if(flag.eq.1)then
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
endif

call mpi_comm_free(red_comm,ierr)

! calculate power spectra, done only from rank 0
if(rank.eq.0)then
! create plan for fft x and fft y
  plan_x=fftw_plan_many_dft_r2c(1,[nx],ny,up,[nx,ny],1,nx,upcx,[nx/2+1,ny],1,nx/2+1,FFTW_ESTIMATE)
  plan_y=fftw_plan_many_dft_r2c(1,[ny],nx,up,[nx,ny],nx,1,upcy,[nx,ny/2+1],nx,1,FFTW_ESTIMATE)

! perform Fourier transform
  call fftw_execute_dft_r2c(plan_x,up,upcx)
  call fftw_execute_dft_r2c(plan_y,up,upcy)
  call fftw_execute_dft_r2c(plan_x,vp,vpcx)
  call fftw_execute_dft_r2c(plan_y,vp,vpcy)
  call fftw_execute_dft_r2c(plan_x,wp,wpcx)
  call fftw_execute_dft_r2c(plan_y,wp,wpcy)

! calculate power spectra
  do j=1,ny
    do i=1,nx/2+1
      stream_ps(i,4)=stream_ps(i,4)+(dble(upcx(i,j)))**2+(aimag(upcx(i,j)))**2
      stream_ps(i,5)=stream_ps(i,5)+(dble(vpcx(i,j)))**2+(aimag(vpcx(i,j)))**2
      stream_ps(i,6)=stream_ps(i,6)+(dble(wpcx(i,j)))**2+(aimag(wpcx(i,j)))**2
    enddo
  enddo

  do j=1,ny/2+1
    do i=1,nx
      span_ps(j,4)=span_ps(j,4)+(dble(upcy(i,j)))**2+(aimag(upcy(i,j)))**2
      span_ps(j,5)=span_ps(j,5)+(dble(vpcy(i,j)))**2+(aimag(vpcy(i,j)))**2
      span_ps(j,6)=span_ps(j,6)+(dble(wpcy(i,j)))**2+(aimag(wpcy(i,j)))**2
    enddo
  enddo

  call fftw_destroy_plan(plan_x)
  call fftw_destroy_plan(plan_y)

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! power spectra at z^+=Re
zplus=re
! find which ranks contain same z^+
zind=minloc(dabs((z+1.0d0)*re-zplus),1)

if((fstart(2)+1.le.zind).and.(fstart(2)+fpz.ge.zind))then
  flag=1
else
 flag=0
endif

! create subcommunicator among ranks that contains z^+
call mpi_comm_split(flow_comm,flag,0,pl_comm,ierr)

utp=0.0d0
vtp=0.0d0
wtp=0.0d0
if(flag.eq.1)then
  do j=1,fpy
    do i=1,nx
        utp(i,fstart(3)+j)=u(i,zind-fstart(2),j)-um(zind)
        vtp(i,fstart(3)+j)=v(i,zind-fstart(2),j)-vm(zind)
        wtp(i,fstart(3)+j)=w(i,zind-fstart(2),j)-wm(zind)
    enddo
  enddo

! reduce to rank floor(rank/nycpu)*nycpu
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,pl_comm,ierr)

! set flag to zero, only reduce rank keeps flag=1
  if(rank.ne.floor(real(rank)/real(nycpu))*nycpu)then
    flag=0
  else
    utp=up
    vtp=vp
    wtp=wp
  endif

endif

call mpi_comm_free(pl_comm,ierr)

! if rank 0 is not in the new communicator it must be included
if(rank.eq.0) flag=1
call mpi_comm_split(flow_comm,flag,0,red_comm,ierr)

! reduce velocity fluctuations to rank 0
if(flag.eq.1)then
  call mpi_reduce(utp,up,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(vtp,vp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
  call mpi_reduce(wtp,wp,nx*ny,mpi_double_precision,mpi_sum,0,red_comm,ierr)
endif

call mpi_comm_free(red_comm,ierr)

! calculate power spectra, done only from rank 0
if(rank.eq.0)then
! create plan for fft x and fft y
  plan_x=fftw_plan_many_dft_r2c(1,[nx],ny,up,[nx,ny],1,nx,upcx,[nx/2+1,ny],1,nx/2+1,FFTW_ESTIMATE)
  plan_y=fftw_plan_many_dft_r2c(1,[ny],nx,up,[nx,ny],nx,1,upcy,[nx,ny/2+1],nx,1,FFTW_ESTIMATE)

! perform Fourier transform
  call fftw_execute_dft_r2c(plan_x,up,upcx)
  call fftw_execute_dft_r2c(plan_y,up,upcy)
  call fftw_execute_dft_r2c(plan_x,vp,vpcx)
  call fftw_execute_dft_r2c(plan_y,vp,vpcy)
  call fftw_execute_dft_r2c(plan_x,wp,wpcx)
  call fftw_execute_dft_r2c(plan_y,wp,wpcy)

! calculate power spectra
  do j=1,ny
    do i=1,nx/2+1
      stream_ps(i,7)=stream_ps(i,7)+(dble(upcx(i,j)))**2+(aimag(upcx(i,j)))**2
      stream_ps(i,8)=stream_ps(i,8)+(dble(vpcx(i,j)))**2+(aimag(vpcx(i,j)))**2
      stream_ps(i,9)=stream_ps(i,9)+(dble(wpcx(i,j)))**2+(aimag(wpcx(i,j)))**2
    enddo
  enddo

  do j=1,ny/2+1
    do i=1,nx
      span_ps(j,7)=span_ps(j,7)+(dble(upcy(i,j)))**2+(aimag(upcy(i,j)))**2
      span_ps(j,8)=span_ps(j,8)+(dble(vpcy(i,j)))**2+(aimag(vpcy(i,j)))**2
      span_ps(j,9)=span_ps(j,9)+(dble(wpcy(i,j)))**2+(aimag(wpcy(i,j)))**2
    enddo
  enddo

  call fftw_destroy_plan(plan_x)
  call fftw_destroy_plan(plan_y)
endif

! mean over all samples and normalize transform
! x: mean over ny samples, normalize transform by nx^2 (square of the modulus)
! y: mean over nx samples, normalize transform by ny^2 (square of the modulus)
stream_ps=stream_ps/dble(nx**2*ny)
span_ps=span_ps/dble(nx*ny**2)


return
end
