subroutine read_fields(nstep)

use commondata
use velocity

integer :: nstep

character(len=40) :: namedir
character(len=8) :: numfile

namedir='../results/'
write(numfile,'(i8.8)') nstep


!allocate(u(nx,nz,ny))
!allocate(uc(nx/2+1,nz,ny,2))
!allocate(v(nx,nz,ny))
!allocate(vc(nx/2+1,nz,ny,2))
!allocate(w(nx,nz,ny))
!allocate(wc(nx/2+1,nz,ny,2))
!allocate(press(nx,nz,ny))
!allocate(pressc(nx/2+1,nz,ny,2))
!if(psiflag.eq.1)then
! allocate(psi(nx,nz,ny))
! allocate(psic(nx/2+1,nz,ny,2))
!endif

if(spectral.eq.1)then
! read in spectral space
!! open(66,file=trim(namedir)//'uc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'uc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) uc
! close(66,status='keep')
!! open(66,file=trim(namedir)//'vc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'vc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) vc
! close(66,status='keep')
!! open(66,file=trim(namedir)//'wc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'wc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) wc
! close(66,status='keep')
! open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 open(66,file=trim(namedir)//'phic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
 read(66) phic
 close(66,status='keep')
 call spectral_to_phys(phic,phi,0)
! open(66,file='../get_pressure/output/pc_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! read(66) pressc
! close(66,status='keep')
!if(psiflag.eq.1)then
!!  open(66,file=trim(namedir)//'psic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
!  open(66,file=trim(namedir)//'psic_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
!  read(66) psic
!  close(66,status='keep')
!endif
else
! read in physical space
!! open(66,file=trim(namedir)//'u_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'u_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) u
! close(66,status='keep')
! call phys_to_spectral(u,uc,0)
!! open(66,file=trim(namedir)//'v_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'v_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) v
! close(66,status='keep')
! call phys_to_spectral(v,vc,0)
!! open(66,file=trim(namedir)//'w_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! open(66,file=trim(namedir)//'w_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
! read(66) w
! close(66,status='keep')
! call phys_to_spectral(w,wc,0)
! open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
 read(66) phi
 close(66,status='keep')
 call phys_to_spectral(phi,phic,0)
! open(66,file='../get_pressure/output/p_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
! read(66) press
! close(66,status='keep')
! call phys_to_spectral(press,pressc,0)
!if(psiflag.eq.1)then
!!  open(66,file=trim(namedir)//'psi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
!  open(66,file=trim(namedir)//'psi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
!  read(66) psi
!  close(66,status='keep')
!  call phys_to_spectral(psi,psic,0)
!endif
endif


! calculate N-S rhs
!call get_sterm


!deallocate(u)
!deallocate(v)
!deallocate(w)
!deallocate(uc)
!deallocate(vc)
!deallocate(wc)
!deallocate(press)
!deallocate(pressc)
!if(psiflag.eq.1)then
! deallocate(psi)
! deallocate(psic)
!endif

call get_interface(nstep)

return
end




subroutine read_grid

use commondata

character(len=40) :: namefile

 write(namefile,'(a)') '../results/x.dat'
! open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 open(1,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(1) x
 close(1,status='keep')


 write(namefile,'(a)') '../results/y.dat'
! open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 open(2,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(2) y
 close(2,status='keep')


 write(namefile,'(a)') '../results/z.dat'
! open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
 open(3,file=trim(namefile),form='unformatted',access='stream',status='old',convert='big_endian')
 read(3) z
 close(3,status='keep')

write(*,'(a,2(f9.4),a)') 'x range: ',x(1),x(nx),' (outer units)'
write(*,'(a,2(f9.4),a)') 'y range: ',y(1),y(ny),' (outer units)'
write(*,'(a,2(f9.4),a)') 'z range: ',z(1),z(nz),' (outer units)'


return
end
