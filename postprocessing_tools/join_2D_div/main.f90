program main

integer :: start,delta,last
integer :: nlines,length
integer :: i,index,io

double precision, allocatable :: data(:,:)
double precision :: tmpd(5)

character(len=200) :: path

call read_input(start,last,delta)

length=0
do i=start,last,delta
 write(*,'(a30,i8,a,i8)') 'Reading number of lines, step ',i,' of ',last
 write(path,'(a,i8.8,a)') '../2D_divergence/output/correlation_',i,'.dat'
 call count_lines(nlines,path,1)
 length=length+nlines
 ! write(*,*) length,nlines
enddo

! write(*,*) length

if(length.ne.0)then
 allocate(data(length,5))

 index=1
 do i=start,last,delta
  write(*,'(a30,i8,a,i8)') 'Reading in data, step ',i,' of ',last
  write(path,'(a,i8.8,a)') '../2D_divergence/output/correlation_',i,'.dat'
  open(455,file=path,status='old',form='formatted')
  read(455,*)
  do
   read(455,'(5(f20.8))',iostat=io) tmpd
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

 open(65,file='./output/lines.dat',form='formatted',status='new')
 write(65,'(i16)') length
 close(65,status='keep')

 open(66,file='./output/psi.dat',form='unformatted',status='new',access='stream')
 write(66) data(:,1)
 close(66,status='keep')

 open(66,file='./output/2D_div_up.dat',form='unformatted',status='new',access='stream')
 write(66) data(:,2)
 close(66,status='keep')

 open(66,file='./output/2D_div_visc.dat',form='unformatted',status='new',access='stream')
 write(66) data(:,3)
 close(66,status='keep')

 open(66,file='./output/2D_div_shear.dat',form='unformatted',status='new',access='stream')
 write(66) data(:,4)
 close(66,status='keep')

 open(66,file='./output/2D_div_mar.dat',form='unformatted',status='new',access='stream')
 write(66) data(:,5)
 close(66,status='keep')

 deallocate(data)

endif

return
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine count_lines(nlines,path,header)

integer :: nlines,io,header

double precision :: rd(5)

character(len=200) :: path

nlines = 0
! write(*,*) trim(path)

open(56,file=path,form='formatted',status='old')
! remove header
if(header.eq.1) read(56,*)

do
 read(56,'(5(f20.8))',iostat=io) rd
 if(io.ne.0)then
  exit
 else
  if((isnan(rd(1)).eqv..false.).or.(isnan(rd(2)).eqv..false.).or.(isnan(rd(3)).eqv..false.) &
 &                             .or.(isnan(rd(2)).eqv..false.).or.(isnan(rd(3)).eqv..false.))then
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
