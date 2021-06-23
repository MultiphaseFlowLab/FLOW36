subroutine initialize_jpdf

use jpdf

integer :: i

double precision, dimension(2) :: Qvb,Rvb,Qsb,Rsb,Qob
double precision :: dqv,drv,dqs,drs,dqo

nbins=500

allocate(QvRv(nbins,nbins,2))
allocate(QsRs(nbins,nbins,2))
allocate(QsQo(nbins,nbins,2))

QvRv=0.0d0
QsRs=0.0d0
QsQo=0.0d0

allocate(Qva(nbins+1))
allocate(Rva(nbins+1))
allocate(Qsa(nbins+1))
allocate(Rsa(nbins+1))
allocate(Qoa(nbins+1))

Qvb=[-1.0d0*10.0d0**(5),1.0d0*10.0d0**(5)]
Rvb=[-1.0d0*10.0d0**(7),1.0d0*10.0d0**(7)]
Qsb=[-1.0d0*10.0d0**(6),0.0d0]
Rsb=[-1.0d0*10.0d0**(8),1.0d0*10.0d0**(8)]
Qob=[-1.0d0*10.0d0**(6),1.0d0*10.0d0**(6)]

write(*,'(a20,2(es10.2))') 'Qv bounds ',Qvb(:)
write(*,'(a20,2(es10.2))') 'Rv bounds ',Rvb(:)
write(*,'(a20,2(es10.2))') 'Qs bounds ',Qsb(:)
write(*,'(a20,2(es10.2))') 'Rs bounds ',Rsb(:)
write(*,'(a20,2(es10.2))') 'Qo bounds ',Qob(:)

! define axis
dqv=(Qvb(2)-Qvb(1))/dble(nbins)
drv=(Rvb(2)-Rvb(1))/dble(nbins)
dqs=(Qsb(2)-Qsb(1))/dble(nbins)
drs=(Rsb(2)-Rsb(1))/dble(nbins)
dqo=(Qob(2)-Qob(1))/dble(nbins)

Qva(1)=Qvb(1)
Rva(1)=Rvb(1)
Qsa(1)=Qsb(1)
Rsa(1)=Rsb(1)
Qoa(1)=Qob(1)

do i=2,nbins+1
 Qva(i)=Qva(i-1)+dqv
 Rva(i)=Rva(i-1)+drv
 Qsa(i)=Qsa(i-1)+dqs
 Rsa(i)=Rsa(i-1)+drs
 Qoa(i)=Qoa(i-1)+dqo
enddo

! write(*,*) Qva(1),Qva(nbins+1)
! write(*,*) Rva(1),Rva(nbins+1)
! write(*,*) Qsa(1),Qsa(nbins+1)
! write(*,*) Rsa(1),Rsa(nbins+1)
! write(*,*) Qoa(1),Qoa(nbins+1)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_jpdf(Qv,Rv,Qs,Rs,Qo,phi)

use jpdf

integer :: i,j

double precision :: Qv,Rv,Qs,Rs,Qo,phi

logical :: ch_qvrv,ch_qsrs,ch_qsqo

ch_qvrv=.false.
ch_qsrs=.false.
ch_qsqo=.false.

do j=1,nbins
 do i=1,nbins
  ! QvRv
  if((Qv.ge.Qva(j)).and.(Qv.lt.Qva(j+1)).and.(Rv.ge.Rva(i)).and.(Rv.lt.Rva(i+1)))then
   if(phi.gt.0.0d0)then
    QvRv(i,j,1)=QvRv(i,j,1)+1
   else
    QvRv(i,j,2)=QvRv(i,j,2)+1
   endif
   ch_qvrv=.true.
   exit
  endif
 enddo
enddo

do j=1,nbins
 do i=1,nbins
  ! QsRs
  if((Qs.ge.Qsa(j)).and.(Qs.lt.Qsa(j+1)).and.(Rs.ge.Rsa(i)).and.(Rs.lt.Rsa(i+1)))then
   if(phi.gt.0.0d0)then
    QsRs(i,j,1)=QsRs(i,j,1)+1
   else
    QsRs(i,j,2)=QsRs(i,j,2)+1
   endif
   ch_qsrs=.true.
   exit
  endif
enddo
enddo

do j=1,nbins
 do i=1,nbins
  ! -QsQo
  if((Qo.ge.Qoa(j)).and.(Qo.lt.Qoa(j+1)).and.(Qs.ge.Qsa(i)).and.(Qs.lt.Qsa(i+1)))then
   if(phi.gt.0.0d0)then
    QsQo(i,j,1)=QsQo(i,j,1)+1
   else
    QsQo(i,j,2)=QsQo(i,j,2)+1
   endif
   ch_qsqo=.true.
   exit
  endif
 enddo
enddo

! write if  value out of bounds!
if(ch_qvrv.eqv..false.) write(*,*) 'QvRv out of bounds',Qv,Rv
if(ch_qsrs.eqv..false.) write(*,*) 'QsRs out of bounds',Qs,Rs
if(ch_qsqo.eqv..false.) write(*,*) 'QoQs out of bounds',Qo,Qs



return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_jpdf

use jpdf

integer :: i

character(len=200) :: fstring


open(555,file='./output/jpdf_QvRv_in.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Rva(1:nbins)+Rva(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qva(i)+Qva(i+1)),QvRv(i,:,1)
enddo
close(555,status='keep')

open(555,file='./output/jpdf_QvRv_out.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Rva(1:nbins)+Rva(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qva(i)+Qva(i+1)),QvRv(i,:,2)
enddo
close(555,status='keep')


open(555,file='./output/jpdf_QsRs_in.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Rsa(1:nbins)+Rsa(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qsa(i)+Qsa(i+1)),QsRs(i,:,1)
enddo
close(555,status='keep')

open(555,file='./output/jpdf_QsRs_out.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Rsa(1:nbins)+Rsa(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qsa(i)+Qsa(i+1)),QsRs(i,:,2)
enddo
close(555,status='keep')


open(555,file='./output/jpdf_QsQo_in.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Qoa(1:nbins)+Qoa(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qsa(i)+Qsa(i+1)),QsQo(i,:,1)
enddo
close(555,status='keep')

open(555,file='./output/jpdf_QsQo_out.dat',status='new',form='formatted')
write(fstring,'(a,i8,a)') '(a16,',nbins,'(es16.5))'
write(555,fstring) '',0.5d0*(Qoa(1:nbins)+Qoa(2:nbins+1))
write(fstring,'(a,i8,a)') '(',nbins+1,'(es16.5))'
do i=1,nbins
 write(555,fstring) 0.5d0*(Qsa(i)+Qsa(i+1)),QsQo(i,:,2)
enddo
close(555,status='keep')


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine clear_jpdf

use jpdf

deallocate(QvRv,QsRs,QsQo)
deallocate(Qva,Rva,Qsa,Rsa,Qoa)

return
end
