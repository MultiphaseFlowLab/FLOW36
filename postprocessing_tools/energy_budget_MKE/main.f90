program main

  use mpi
  use commondata
  use fields
  use wavenumber

  integer :: i,j

  call mpi_init(ierr)

  call mpi_comm_size(mpi_comm_world,ntask,ierr)
  call mpi_comm_rank(mpi_comm_world,rank,ierr)

  call read_input

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))
  allocate(pm(nz))
  allocate(prms(nz))
  allocate(press(nx,nz,ny))
  allocate(u(nx,nz,ny))
  allocate(v(nx,nz,ny))
  allocate(w(nx,nz,ny))
  allocate(uc(nx/2+1,nz,ny,2))
  allocate(vc(nx/2+1,nz,ny,2))
  allocate(wc(nx/2+1,nz,ny,2))
  allocate(phi(nx,nz,ny))
  allocate(phic(nx/2+1,nz,ny,2))
  if(phi_flag.eq.0) phi=-1.0d0
  if(psi_flag.eq.1) then
    allocate(psi(expx*nx,expz*(nz-1)+1,expy*ny))
    allocate(psic((expx*nx)/2+1,expz*(nz-1)+1,expy*ny,2))
  endif

  allocate(kx(nx/2+1))
  allocate(ky(ny))
  allocate(k2(nx/2+1,ny))
  allocate(kxfg(nxfg/2+1))
  allocate(kyfg(nyfg))
  allocate(k2fg(nxfg/2+1,nyfg))

  ! wavenumbers
  kx(1)=0.0d0
  do i=2,nx/2+1
    kx(i)=dble(i-1)*2.0d0*pi/xl
  enddo

  ky(1)=0.0d0
  do i=2,ny/2+1
    ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
    ky(i)=dble(i-1)*2.0d0*pi/yl
  enddo

  do j=1,ny
    do i=1,nx/2+1
      k2(i,j)=kx(i)*kx(i)+ky(j)*ky(j)
    enddo
  enddo

  kxfg(1)=0.0d0
  do i=2,nxfg/2+1
    kxfg(i)=dble(i-1)*2.0d0*pi/xl
  enddo

  kyfg(1)=0.0d0
  do i=2,nyfg/2+1
    kyfg(nyfg-i+2)=-dble(i-1)*2.0d0*pi/yl
    kyfg(i)=dble(i-1)*2.0d0*pi/yl
  enddo

  do j=1,nyfg
    do i=1,nxfg/2+1
      k2fg(i,j)=kxfg(i)*kxfg(i)+kyfg(j)*kyfg(j)
    enddo
  enddo

  ! plan creation
  call create_plan
  call create_plan_fg
  call create_plan_2D

  call read_grid

  ! call phys_to_spectral(u,uc,0)
  ! call spectral_to_phys(uc,u,0)
  ! call phys_to_spectral_fg(psi,psic,0)
  ! call spectral_to_phys_fg(psic,psi,0)

  open(1456,file='./output/integral_mean.dat',status='new',form='formatted')
  write(1456,'(a10,14(a15))') 'step','t^+','press_work c','visc_diff_c','visc_diss_c','prod_c','turb_diff_c','gravity_c', &
    &                                 'press_diff_d','visc_diff_d','visc_diss_d','prod_d','turb_diff_d','gravity_d','interf'
  close(1456,status='keep')


  do i=nstart,nend,delta
    call calc_budget_mean(i)
  enddo
  ! run last step if not already run
  if(mod(nend,delta).ne.0) call calc_budget_mean(nend)



  call destroy_plan

  deallocate(x,y,z)
  deallocate(pm,prms)
  deallocate(press)
  deallocate(u,v,w,uc,vc,wc)
  deallocate(phi,phic)
  if(psi_flag.eq.1) deallocate(psi,psic)
  deallocate(kx,ky,k2,kxfg,kyfg,k2fg)


  call mpi_finalize(ierr)

  return
end program main
