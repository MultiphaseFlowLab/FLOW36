program change_grid

use commondata

character(len=5) :: namevar

write(*,*)
write(*,*) 'Starting program change grid'

call read_input


! porting u
write(*,*)
namevar='u'
call port_grid(namevar)


! porting v
write(*,*)
namevar='v'
call port_grid(namevar)

! porting w
write(*,*)
namevar='w'
call port_grid(namevar)

! porting phi
if(phiflag.eq.1)then
  write(*,*)
  namevar='phi'
!  namevar='T'
  call port_grid(namevar)
  if(psiflag.eq.1)then
    write(*,*)
    namevar='psi'
    call port_grid(namevar)
  endif
endif

end program



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_input

use input
use output
use commondata

 open(66,file='./input.f90',form='formatted',status='old')
! input grid
 read(66,*)
 read(66,'(i8)') inx
 read(66,'(i8)') inz
 read(66,'(i8)') iny
! output grid
 read(66,*)
 read(66,'(i8)') onx
 read(66,'(i8)') onz
 read(66,'(i8)') ony
! read phi flag
 read(66,*)
 read(66,'(i8)') phiflag
 read(66,*)
 read(66,'(i8)') psiflag

 close(66,status='keep')


write(*,'(1x,3(a,i6))') 'Change grid from: ',inx,' x',iny,' x',inz
write(*,'(15x,3(a,i6))') 'to: ',onx,' x',ony,' x',onz

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine port_grid(namevar)

use input
use output
use commondata

character(len=5) :: namevar
character(len=40) :: fname

double precision, allocatable :: iu(:,:,:,:), ou(:,:,:,:)
double precision, allocatable :: iru(:,:,:), oru(:,:,:)

integer :: tinx,tiny,tinz

nx=inx
ny=iny
nz=inz

allocate(iu(inx/2+1,inz,iny,2))
allocate(iru(nx,nz,ny))

! reading input variable
fname='./input/'//trim(namevar)//'.dat'
write(*,*) 'Reading  ',trim(fname)

open(unit=66,file=trim(fname),form='unformatted',status='old',access='stream')!,convert='little_endian')
!open(unit=66,file=trim(fname),form='unformatted',status='old',access='stream',convert='big_endian')
 read(66) iru
 close(66,status='keep')

! debug only
! check for correct normalization


call create_plan
call phys_to_spectral(iru,iu,0)
write(*,'(a16,2(f20.16))') 'Input max, min',maxval(iru),minval(iru)
deallocate(iru)

call destroy_plan

! end debug only


! porting to output grid
allocate(ou(onx/2+1,onz,ony,2))
tinx=inx
tiny=iny
tinz=inz

if(inx.gt.onx) tinx=onx
if(iny.gt.ony) tiny=ony
if(inz.gt.onz) tinz=onz


ou(1:tinx/2+1,1:tinz,1:tiny/2+1,1:2)=iu(1:tinx/2+1,1:tinz,1:tiny/2+1,1:2)
ou(1:tinx/2+1,1:tinz,ony-(tiny/2-1)+1:ony,1:2)=iu(1:tinx/2+1,1:tinz,tiny/2+2:tiny,1:2)


! renormalize
ou=ou*onx/inx*ony/iny/1.d0;

! transform to physical space
nx=onx
ny=ony
nz=onz

call create_plan

allocate(oru(nx,nz,ny))

call spectral_to_phys(ou,oru,0)

write(*,'(a16,2(f20.16))') 'Output max, min ',maxval(oru),minval(oru)

call destroy_plan

deallocate(iu)
deallocate(ou)

! output in physical space
fname='./output/'//trim(namevar)//'.dat'
write(*,*) 'Writing  ',trim(fname)

open(unit=67,file=trim(fname),form='unformatted',status='new',access='stream',convert='little_endian')
!open(unit=67,file=trim(fname),form='unformatted',status='new',access='stream',convert='big_endian')
 write(67) oru
 close(67,status='keep')

deallocate(oru)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
