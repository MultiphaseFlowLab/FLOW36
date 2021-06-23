program main

integer, parameter :: nbins=201
integer :: start,delta,last
integer :: nlines,length
integer :: i,j,index,io

double precision, allocatable :: data(:,:)
double precision :: tmpd(4),pdf(nbins,9),bins(nbins-1),dx

character(len=200) :: path

call read_input(start,last,delta)

length=0
do i=start,last,delta
 write(*,'(a30,i8,a,i8)') 'Reading number of lines, step ',i,' of ',last
 write(path,'(a,i8.8,a)') '../marangoni/output/correlation_',i,'.dat'
 call count_lines(nlines,path,1)
 length=length+nlines
 ! write(*,*) length,nlines
enddo

! write(*,*) length

if(length.ne.0)then
 allocate(data(length,4))

 index=1
 do i=start,last,delta
  write(*,'(a30,i8,a,i8)') 'Generating PDF, step ',i,' of ',last
  write(path,'(a,i8.8,a)') '../marangoni/output/correlation_',i,'.dat'
  open(455,file=path,status='old',form='formatted')
  read(455,*)
  do
   read(455,'(4(f16.8))',iostat=io) tmpd
   if(io.ne.0)then
    exit
   else
    data(index,:)=tmpd
    ! write(*,*) data(index,:)
    index=index+1
   endif
  enddo
  close(455,status='keep')
 enddo

 ! read in all data, generate PDF
 pdf=0.0d0
 ! samples between -1 and 1, generate pdf axis
 dx=2.0d0/dble(nbins-1)
 pdf(1,1)=-1.0d0
 do i=2,nbins
  pdf(i,1)=pdf(i-1,1)+dx
  ! write(*,*) pdf(i,1)
  bins(i-1)=(pdf(i,1)+pdf(i-1,1))/2.0d0
 enddo


 do i=1,length
  ! PDF u
  if(data(i,1).lt.bins(1))then
   ! write(*,*) '.lt. ',bins(1)
   pdf(1,2)=pdf(1,2)+1
  elseif(data(i,1).ge.bins(nbins-1))then
   ! write(*,*) '.ge. ',bins(nbins-1)
   pdf(nbins,2)=pdf(nbins,2)+1
  else
   do j=1,nbins-2
    ! write(*,*) bins(j),' .lt. x .le. ',bins(j+1)
    if(data(i,1).ge.bins(j).and.data(i,1).lt.bins(j+1))then
     pdf(j+1,2)=pdf(j+1,2)+1
     exit
    endif
   enddo
  endif

  ! PDF u'
  if(data(i,2).lt.bins(1))then
   ! write(*,*) '.lt. ',bins(1)
   pdf(1,3)=pdf(1,3)+1
  elseif(data(i,2).ge.bins(nbins-1))then
   ! write(*,*) '.ge. ',bins(nbins-1)
   pdf(nbins,3)=pdf(nbins,3)+1
  else
   do j=1,nbins-2
    ! write(*,*) bins(j),' .lt. x .le. ',bins(j+1)
    if(data(i,2).ge.bins(j).and.data(i,2).lt.bins(j+1))then
     pdf(j+1,3)=pdf(j+1,3)+1
     exit
    endif
   enddo
  endif

  ! PDF shear stress
  if(data(i,3).lt.bins(1))then
   ! write(*,*) '.lt. ',bins(1)
   pdf(1,4)=pdf(1,4)+1
  elseif(data(i,3).ge.bins(nbins-1))then
   ! write(*,*) '.ge. ',bins(nbins-1)
   pdf(nbins,4)=pdf(nbins,4)+1
  else
   do j=1,nbins-2
    ! write(*,*) bins(j),' .lt. x .le. ',bins(j+1)
    if(data(i,3).ge.bins(j).and.data(i,3).lt.bins(j+1))then
     pdf(j+1,4)=pdf(j+1,4)+1
     exit
    endif
   enddo
  endif

  ! PDF viscous stress
  if(data(i,4).lt.bins(1))then
   ! write(*,*) '.lt. ',bins(1)
   pdf(1,5)=pdf(1,5)+1
  elseif(data(i,4).ge.bins(nbins-1))then
   ! write(*,*) '.ge. ',bins(nbins-1)
   pdf(nbins,5)=pdf(nbins,5)+1
  else
   do j=1,nbins-2
    ! write(*,*) bins(j),' .lt. x .le. ',bins(j+1)
    if(data(i,4).ge.bins(j).and.data(i,4).lt.bins(j+1))then
     pdf(j+1,5)=pdf(j+1,5)+1
     exit
    endif
   enddo
  endif
 enddo

 deallocate(data)


 ! normalize PDF
 tmpd(1)=(bins(1)+1.0d0)*pdf(1,2)
 tmpd(2)=(bins(1)+1.0d0)*pdf(1,3)
 tmpd(3)=(bins(1)+1.0d0)*pdf(1,4)
 tmpd(4)=(bins(1)+1.0d0)*pdf(1,5)
 do j=1,nbins-2
  tmpd(1)=tmpd(1)+(bins(j+1)-bins(j))*pdf(j+1,2)
  tmpd(2)=tmpd(2)+(bins(j+1)-bins(j))*pdf(j+1,3)
  tmpd(3)=tmpd(3)+(bins(j+1)-bins(j))*pdf(j+1,4)
  tmpd(4)=tmpd(4)+(bins(j+1)-bins(j))*pdf(j+1,5)
 enddo
 tmpd(1)=tmpd(1)+(1.0d0-bins(nbins-1))*pdf(nbins,2)
 tmpd(2)=tmpd(2)+(1.0d0-bins(nbins-1))*pdf(nbins,3)
 tmpd(3)=tmpd(3)+(1.0d0-bins(nbins-1))*pdf(nbins,4)
 tmpd(4)=tmpd(4)+(1.0d0-bins(nbins-1))*pdf(nbins,5)

 pdf(:,6)=pdf(:,2)/tmpd(1)
 pdf(:,7)=pdf(:,3)/tmpd(2)
 pdf(:,8)=pdf(:,4)/tmpd(3)
 pdf(:,9)=pdf(:,5)/tmpd(4)

 ! write PDF to file (non-normalized and normalized)
 open(666,file='./output/pdf.dat',status='new',form='formatted')
 write(666,'(9(a24))') 'axis','pdf u (non-norm)','pdf up (non-norm)','pdf shear (non-norm)','pdf visc (non-norm)', &
 &                            'pdf u','pdf up','pdf shear','pdf visc'
 do i=1,nbins
  write(666,'(9(es24.6))') pdf(i,:)
 enddo
 close(666,status='keep')

endif

return
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine count_lines(nlines,path,header)

integer :: nlines,io,header

double precision :: rd(4)

character(len=200) :: path

nlines = 0
! write(*,*) trim(path)

open(56,file=path,form='formatted',status='old')
! remove header
if(header.eq.1) read(56,*)

do
 read(56,'(4(f16.8))',iostat=io) rd
 if(io.ne.0)then
  exit
 else
  if((isnan(rd(1)).eqv..false.).or.(isnan(rd(2)).eqv..false.).or.(isnan(rd(3)).eqv..false.).or.(isnan(rd(4)).eqv..false.))then
   nlines=nlines+1
  endif
 endif
enddo

close(56,status='keep')



return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_input(start,last,delta)

integer :: start,last,delta
integer :: tmpst,tmpls,tmpdlt

open(456,file='../sc_compiled/input.f90',form='formatted',status='old')
read(456,*) ! type of simulation
read(456,*) !0                     ! restart flag, if 1 is a restart
read(456,*) !0               ! number of iteration used as a restart field
read(456,*) !0                ! initial conditions
read(456,*) ! grid parameters
read(456,*) !64                         ! number of points in x
read(456,*) !32                         ! number of points in y
read(456,*) !129                         ! number of points in z
read(456,*) ! simulation parameters
read(456,*) !1.0                           ! Reynolds number
read(456,*) !0.2                      ! Courant number
read(456,*) !0.0                          ! mean pressure gradient x direction: Delta p/ Lx = (p_out-p_in)/Lx
read(456,*) !0.0                          ! mean pressure gradient y direction: Delta p/ Ly
read(456,*) !0                         ! 1: CPI activated, 0: CPI deactivated
read(456,*) !30.6                         ! Power reynolds number
read(456,*) ! domain size
read(456,*) !4.0                           ! Lx
read(456,*) !2.0                           ! Ly
read(456,*) !
read(456,'(i8)') start
read(456,'(i8)') last
read(456,'(i8)') delta
read(456,*) !-1                          ! solution dumping frequency in spectral space (if lower than 1 never save solution)
close(456,status='keep')

open(543,file='./input.f90',form='formatted',status='old')
read(543,'(i8)') tmpst
read(543,'(i8)') tmpls
read(543,'(i8)') tmpdlt
close(543,status='keep')

if(tmpst.ne.-1.and.tmpst.gt.start) start=tmpst
if(tmpls.ne.-1.and.tmpls.lt.last) last=tmpls
if(tmpdlt.ne.-1.and.tmpdlt.gt.delta) delta=tmpdlt

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
