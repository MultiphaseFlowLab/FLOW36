! subroutine explore(dist)
!
! use commondata
! use mpi
!
! double precision, dimension(nx,nz,ny) :: sgn,dist
! integer :: i,j,k
! integer :: il,iu,jl,ju,kl,ku,delta,spx,spy,spz,ii,jj,kk
! ! logical, dimension(nx,nz,ny) :: mask
!
! ! double precision :: seed,table(nx,nz,ny),eudist(nx,nz,ny),posit
! ! double precision :: seed,eudist(nx,nz,ny),posit,tmp(nx)
! double precision :: seed,posit,tmp(nx),len
! ! double precision, allocatable, dimension(:,:,:) :: eudist
! double precision :: ts,te
!
! dist=0.0d0
! call sign_f(phi,sgn)
!
! delta=int(max(nx,ny,nz)/4)
! ! spx=min(nx,2*delta+1)
! ! spy=min(ny,2*delta+1)
! ! spz=min(nz,2*delta+1)
! ! allocate(eudist(spx,spz,spy))
! len=dsqrt((x(nx)-x(1))**2+(y(ny)-y(1))**2+(2*Re)**2)
!
!
! ts=mpi_wtime()
! do j=1,ny
!  do k=1,nz
!   do i=1,nx
!    ! with expanding to spherical shell
!    ! reinitialize posit
!    posit=len
!
!
!
!
!    ! with min distance but missing periodicity
!    ! seed=sgn(i,k,j)
!    !   ! look for interface + missing periodicity
!    !   ! eudist=(xyz(:,:,:,1)-x(i))**2+(xyz(:,:,:,2)-z(k))**2+(xyz(:,:,:,3)-y(j))**2
!    !   ! find min within mask
!    !   posit=minval(eudist,seed*sgn.le.0)
!    !   ! dist(i,k,j)=-seed*dsqrt(posit)
!    !   ! tmp fits in cache, small speed-up
!    !   tmp(i)=-seed*dsqrt(posit)
!    enddo
!    dist(:,k,j)=tmp
!  enddo
! enddo
!
! te=mpi_wtime()
! write(*,*) te-ts,' seconds'
!
!
! ! ! max distance lower than diagonal of domain
! ! ! not really needed without inner loop
! ! posit=dsqrt((x(nx)-x(1))**2+(y(ny)-y(1))**2+(2*Re)**2)
! !
! ! if(seed.eq.0.0d0)then
! !  ! found interface
! !  dist(ic,kc,jc)=0.0d0
! ! else
! !  ! look for interface
! !  eudist=dsqrt((xyz(:,:,:,1)-x(ic))**2+(xyz(:,:,:,2)-z(kc))**2+(xyz(:,:,:,3)-y(jc))**2)
! !  ! find min within mask
! !  mask=seed*sgn.le.0
! !  posit=minval(eudist,mask)
! !  dist(ic,kc,jc)=-seed*posit
! !
! !
! !  ! table=seed*sgn
! !  ! do j=1,ny
! !  !  do k=1,nz
! !  !   do i=1,nx
! !  !    if(table(i,k,j).le.0.0d0) posit=min(posit,dsqrt((x(i)-x(ic))**2+(y(j)-y(jc))**2+(z(k)-z(kc))**2))
! !  !   enddo
! !  !  enddo
! !  ! enddo
! !  ! ! distance from interface, positive outside droplet, negative inside droplet
! !  ! dist(ic,kc,jc)=-seed*posit
! ! endif
!
!
! return
! end



! subroutine explore(dist)
!
! use commondata
! use mpi
!
! double precision, dimension(nx,nz,ny) :: dist,sgn
! double precision, dimension(nx) :: tmp
! double precision, allocatable, dimension(:,:,:) :: esgn,eudist
! double precision, allocatable, dimension(:) :: xd,yd,zd
! integer :: i,j,k  , ii,jj,kk
! integer :: il,iu,jl,ju,kl,ku
! double precision :: ts,te,seed
!
!
! dist=0.0d0
!
! call sign_f(phi,sgn)
! ! expand domain in x and y to consider periodicity
! allocate(esgn(nx+2*deltax,nz,ny+2*deltay))
!
! esgn(1+deltax:nx+deltax,:,1+deltay:ny+deltay)=sgn
! ! periodic directions x
! esgn(1:deltax,:,1+deltay:ny+deltay)=sgn(nx-deltax+1:nx,:,:)
! esgn(nx+deltax+1:nx+2*deltax,:,1+deltay:ny+deltay)=sgn(1:deltax,:,:)
! ! periodic directions y
! esgn(1+deltax:nx+deltax,:,1:deltay)=sgn(:,:,ny-deltay+1:ny)
! esgn(1+deltax:nx+deltax,:,ny+deltay+1:ny+2*deltay)=sgn(:,:,1:deltay)
!
!
! allocate(eudist(2*deltax+1,2*deltaz+1,2*deltay+1))
! allocate(xd(2*deltax+1))
! allocate(yd(2*deltay+1))
! allocate(zd(2*deltaz+1))
!
! ts=mpi_wtime()
! do j=1+deltay,ny+deltay
!  do k=1,nz
!   do i=1+deltax,nx+deltax
! ! do j=1,ny
! !  do k=1,nz
! !   do i=1,nx
!    ! look for closest point in box 2*deltax+1,2*deltaz+1,2*deltay+1
!    ! il=max(1,i-deltax)
!    ! iu=min(nx,i+deltax)
!    ! jl=max(1,j-deltay)
!    ! ju=min(ny,j+deltay)
!    il=i-deltax
!    iu=i+deltax
!    jl=j-deltay
!    ju=j+deltay
!    kl=max(1,k-deltaz)
!    ku=min(nz,k+deltaz)
!
!
!    do ii=1,iu-il+1
!      xd(ii)=(xx(il+ii-1)-xx(i))**2
!    end do
!
!    do jj=1,ju-jl+1
!      yd(jj)=(yy(jj+jl-1)-yy(j))**2
!    end do
!
!    do kk=1,ku-kl+1
!      zd(kk)=(z(kk+kl-1)-z(k))**2
!    end do
!
!    do jj=1,ju-jl+1
!      do kk=1,ku-kl+1
!        do ii=1,iu-il+1
!          eudist(ii,kk,jj)=xd(ii)+yd(jj)+zd(kk)
!        enddo
!      enddo
!    enddo
!
!
!
!  !   eudist(:,1:ku-kl+1,:)=(xyz(il:iu,kl:ku,jl:ju,1)-xx(i))**2 &
!  ! &                      +(xyz(il:iu,kl:ku,jl:ju,2)-z(k))**2 &
!  ! &                      +(xyz(il:iu,kl:ku,jl:ju,3)-yy(j))**2
!
!    seed=esgn(i,k,j)
!    tmp(i-deltax)=minval(eudist(:,1:ku-kl+1,:),seed*esgn(il:iu,kl:ku,jl:ju).le.0)
!    tmp(i-deltax)=-seed*dsqrt(tmp(i-deltax))
!
!   enddo
!   dist(:,k,j-deltay)=tmp
!  enddo
! enddo
!
! te=mpi_wtime()
! write(*,*) te-ts,' seconds'
!
!
! deallocate(eudist)
! deallocate(esgn)
! deallocate(xd,yd,zd)
!
!
! return
! end
