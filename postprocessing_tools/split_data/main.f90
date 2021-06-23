program split_data

logical :: ex
integer :: step,count,ios,step_old
character(len=8) :: nstep
double precision :: xg,yg,zg,V,A

open(66,file='../mass_center/output/mass_center.dat',status='old',form='formatted')
read(66,*,iostat=ios)

step_old=-7
do
 read(66,'(i8,i6,5(es16.6))',iostat=ios) step,count,xg,yg,zg,V,A
 if(ios.lt.0) then
   write(*,*) 'EOF reached, stopping'
   exit
 endif

 if(step.ne.step_old)  write(*,*) 'Step: ',step 
 step_old=step

 write(nstep,'(i8.8)') step
 inquire(file='./output/mass_center_'//nstep//'.dat',exist=ex)
 if(ex.eqv..false.) then
  open(68,file='./output/mass_center_'//nstep//'.dat',status='new',form='formatted',access='append')
  write(68,'(5(a16))') 'xg','yg','zg','V','A'
 else
  open(68,file='./output/mass_center_'//nstep//'.dat',status='old',form='formatted',access='append')
 endif
 write(68,'(5(es16.6))') xg,yg,zg,V,A 
 close(68,status='keep')

enddo


return
end program
