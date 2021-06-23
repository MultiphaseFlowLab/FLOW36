program create_pdf


double precision, allocatable :: samples(:),pdf(:,:)
double precision :: val,maxv,minv,delta

integer :: n_samples,i,j,count,pos(1)
integer :: start,finish,step,nbins

character(len=8) :: time


open(44,file='pdf_input.dat',form='formatted',status='old')
read(44,'(i8)') start
read(44,'(i8)') finish
read(44,'(i8)') step
read(44,*)
read(44,'(i8)') nbins
close(44,status='keep')


open(25,file='./output/total_s.dat',status='replace',form='formatted')

count=0

do i=start,finish,step
  write(time,'(i8.8)') i
  write(*,'(a,i8,a,i8)') 'Iteration ',i,' of ',finish
  open(23,file='./output/n_samples_'//time//'.dat',status='old',form='formatted')
  read(23,'(20x,i20)') n_samples
  close(23,status='keep')
  count=count+n_samples

  open(24,file='./output/samples_'//time//'.dat',status='old',form='formatted')
  allocate(samples(n_samples))
  do j=1,n_samples
    read(24,'(e16.6)') samples(j)
    write(25,'(e16.6)') samples(j)
  enddo
  close(24,status='keep')

  ! generate pdf
  minv=minval(samples)
  maxv=maxval(samples)
  delta=(maxv-minv)/dble(nbins-3)

  allocate(pdf(nbins,2))
  pdf(:,2)=0.0d0
  pdf(1,1)=minv-delta
  do j=2,nbins
   pdf(j,1)=pdf(j-1,1)+delta
  enddo

  ! populate PDF
  do j=1,n_samples
    pos=minloc(abs(samples(j)-pdf(:,1)))
    pdf(pos(1),2)=pdf(pos(1),2)+1.0d0
  enddo

  deallocate(samples)

  ! normalize PDF
  pdf(:,2)=pdf(:,2)/(delta*sum(pdf(:,2)))
  ! write output PDF
  open(434,file='./output/pdf_'//time//'.dat',form='formatted',status='replace')
  do j=1,nbins
    write(434,'(2(e16.6))') pdf(j,1),pdf(j,2)
  enddo
  close(434,status='keep')

  deallocate(pdf)
enddo

close(25,status='keep')

write(*,'(a,i12)') 'Total number of samples: ',count

! read samples
allocate(samples(count))

open(25,file='./output/total_s.dat',status='old',form='formatted')
do i=1,count
  read(25,'(e16.6)') samples(i)
enddo
close(25,status='delete')

! generate pdf
minv=minval(samples)
maxv=maxval(samples)

delta=(maxv-minv)/dble(nbins-3)

allocate(pdf(nbins,2))
pdf(:,2)=0.0d0
pdf(1,1)=minv-delta
do i=2,nbins
 pdf(i,1)=pdf(i-1,1)+delta
enddo


! populate PDF
do i=1,count
  pos=minloc(abs(samples(i)-pdf(:,1)))
  pdf(pos(1),2)=pdf(pos(1),2)+1.0d0
enddo

deallocate(samples)

! normalize PDF
pdf(:,2)=pdf(:,2)/(delta*sum(pdf(:,2)))


! write output PDF
open(434,file='./output/pdf.dat',form='formatted',status='replace')
do i=1,nbins
 write(434,'(2(e16.6))') pdf(i,1),pdf(i,2)
enddo

close(434,status='keep')


deallocate(pdf)

return
end

