subroutine join_samples

use commondata

integer :: index,dump
double precision, allocatable, dimension(:,:) :: data
double precision :: tmpd(5)
character(len=200) :: path

if(spectral.eq.1)then
 dump=sdump
else
 dump=ndump
endif

! samples inside droplet
if(samples_in.gt.0)then
 allocate(data(samples_in,5))
 index=1
 do i=nstart,nend,dump
  write(*,'(a30,i8,a,i8)') 'Reading in data (in), step ',i,' of ',nend
  write(path,'(a,i8.8,a)') './output/inv_in_',i,'.dat'
  open(455,file=trim(path),status='old',form='formatted')
  read(455,*) ! header
  do
   read(455,'(5(es18.6))',iostat=io) tmpd
   if(io.ne.0)then
    exit
   else
    data(index,:)=tmpd
    index=index+1
   endif
  enddo
  close(455,status='delete')
 enddo

 ! write output files
 ! data contains (in order): Q_v,R_v,Q_s,R_s,Q_omega
 open(55,file='./output/nlines_in.dat',status='new',form='formatted')
 write(55,'(i20)') samples_in
 close(55,status='keep')

 open(556,file='./output/Q_v_in.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,1)
 close(556,status='keep')

 open(556,file='./output/R_v_in.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,2)
 close(556,status='keep')

 open(556,file='./output/Q_s_in.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,3)
 close(556,status='keep')

 open(556,file='./output/R_s_in.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,4)
 close(556,status='keep')

 open(556,file='./output/Q_omega_in.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,5)
 close(556,status='keep')

 deallocate(data)
endif


! samples outside droplet
if(samples_out.gt.0)then
 allocate(data(samples_out,5))
 index=1
 do i=nstart,nend,dump
  write(*,'(a30,i8,a,i8)') 'Reading in data (out), step ',i,' of ',nend
  write(path,'(a,i8.8,a)') './output/inv_out_',i,'.dat'
  open(455,file=trim(path),status='old',form='formatted')
  read(455,*) ! header
  do
   read(455,'(5(es18.6))',iostat=io) tmpd
   if(io.ne.0)then
    exit
   else
    data(index,:)=tmpd
    index=index+1
   endif
  enddo
  close(455,status='delete')
 enddo

 ! write output files
 ! data contains (in order): Q_v,R_v,Q_s,R_s,Q_omega
 open(55,file='./output/nlines_out.dat',status='new',form='formatted')
 write(55,'(i20)') samples_out
 close(55,status='keep')

 open(556,file='./output/Q_v_out.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,1)
 close(556,status='keep')

 open(556,file='./output/R_v_out.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,2)
 close(556,status='keep')

 open(556,file='./output/Q_s_out.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,3)
 close(556,status='keep')

 open(556,file='./output/R_s_out.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,4)
 close(556,status='keep')

 open(556,file='./output/Q_omega_out.dat',status='new',form='unformatted',access='stream')
 write(556) data(:,5)
 close(556,status='keep')

 deallocate(data)
endif


return
end
