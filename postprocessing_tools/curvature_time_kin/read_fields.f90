subroutine read_fields(nstep)

use commondata

double precision,dimension(nz) :: um,vm,wm

integer :: nstep,i,j,k

character(len=40) :: namedir
character(len=8) :: numfile

namedir='../results/'
write(numfile,'(i8.8)') nstep

if(spectral.eq.1)then
! read in spectral space
 open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(66) phic
 close(66,status='keep')
 ! outward pointing normal: change sign of phic
 phic=-phic
 call spectral_to_phys(phic,phi,0)
else
! read in physical space
 open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(66) phi
 close(66,status='keep')
 ! outward pointing normal: change sign of phi
 phi=-phi
 ! call phys_to_spectral(phi,phic,0)
endif

open(899,file=trim(namedir)//'u_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
read(899) u
close(899,status='keep')

open(898,file=trim(namedir)//'v_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
read(898) v
close(898,status='keep')

open(897,file=trim(namedir)//'w_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
read(897) w
close(897,status='keep')


um=0.0d0
vm=0.0d0
wm=0.0d0
do j=1,ny
  do k=1,nz
    do i=1,nx
      um(k)=um(k)+u(i,k,j)
      vm(k)=vm(k)+v(i,k,j)
      wm(k)=wm(k)+w(i,k,j)
    enddo
  enddo
enddo
um=um/dble(nx*ny)
vm=vm/dble(nx*ny)
wm=wm/dble(nx*ny)

do j=1,ny
  do k=1,nz
    do i=1,nx
      kin(i,k,j)=0.5d0*((u(i,k,j)-um(k))**2+(v(i,k,j)-vm(k))**2+(w(i,k,j)-wm(k))**2)
    enddo
  enddo
enddo

! phi rescaled by a factor 50 to reduce oscillation in gradient calculation
! rescaling does not affect normal calculation (normal is normalized to modulo 1)
! phi=phi/50.0d0
! phic=phic/50.0d0

call calc_curvature(nstep)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_grid

use commondata

character(len=40) :: namefile

 write(namefile,'(a)') '../results/x.dat'
 open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(1) x
 close(1,status='keep')


 write(namefile,'(a)') '../results/y.dat'
 open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(2) y
 close(2,status='keep')


 write(namefile,'(a)') '../results/z.dat'
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
!  open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(3) z
 close(3,status='keep')

return
end
