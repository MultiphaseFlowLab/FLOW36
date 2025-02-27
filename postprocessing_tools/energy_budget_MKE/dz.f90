subroutine dz(varin,duz)
  ! be careful, do not use same input and output in the call,
  ! e.g. call dz(u,u) <---- NO!

  use commondata

  double precision :: varin(nx/2+1,nz,ny,2),duz(nx/2+1,nz,ny,2)

  integer :: k,m

  do m=1,2
    duz(:,nz,:,m)=0.0d0
    duz(:,nz-1,:,m)= 2.0d0*(nz-1)*varin(:,nz,:,m)
    do k=nz-2,1,-1
      duz(:,k,:,m)=duz(:,k+2,:,m)+2.0d0*dble(k)*varin(:,k+1,:,m)
    enddo
  enddo

  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dz_fg(varin,duz)
  ! be careful, do not use same input and output in the call,
  ! e.g. call dz(u,u) <---- NO!

  use commondata

  double precision :: varin(nxfg/2+1,nzfg,nyfg,2),duz(nxfg/2+1,nzfg,nyfg,2)

  integer :: k

  do m=1,2
    duz(:,nzfg,:,m)=0.0d0
    duz(:,nzfg-1,:,m)= 2.0d0*(nzfg-1)*varin(:,nzfg,:,m)
    do k=nzfg-2,1,-1
      duz(:,k,:,m)=duz(:,k+2,:,m)+2.0d0*dble(k)*varin(:,k+1,:,m)
    enddo
  enddo

  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dz_1d(varin,duz)
  ! be careful, do not use same input and output in the call,
  ! e.g. call dz(u,u) <---- NO!

  use commondata

  double precision :: varin(nz),duz(nz)

  integer :: k

  duz(nz) = 0.0d0
  duz(nz-1) = 2.0d0*(nz-1)*varin(nz)
  do k=nz-2,1,-1
    duz(k) = duz(k+2) + 2.0d0*dble(k)*varin(k+1)
  enddo

  return
end subroutine
