subroutine initializeGPU()

  use commondata
  use velocity
  use par_size
  use CUDA_modu
  use grid
  use sim_par

  !implicit none

  integer :: i,j,k,ind,index


!!velocity components
  allocate(ur_t(spx*nz*spy))
  allocate(uc_t(spx*nz*spy))

  allocate(vr_t(spx*nz*spy))
  allocate(vc_t(spx*nz*spy))
  
  allocate(wr_t(spx*nz*spy))
  allocate(wc_t(spx*nz*spy))
  
  !grid parameters
  allocate(z_t(nz))
  allocate(fstart_t(3))
  
  
  !pass integers
  nx_f   = nx
  spx_f  = spx
  spy_f  = spy
  nz_f   = nz
  xl_f   = xl
  yl_f   = yl
  dt_f   = dt
  
  
  write(*,*)"passing..."
  ind=1
  do k=1,spy
    do j=1,nz
      do i=1,spx
        ur_t(ind) = uc(i,j,k,1)        
        uc_t(ind) = uc(i,j,k,2)

        vr_t(ind) = vc(i,j,k,1)        
        vc_t(ind) = vc(i,j,k,2)
      
        wr_t(ind) = wc(i,j,k,1)        
        wc_t(ind) = wc(i,j,k,2)
        
        ind = ind + 1
      end do
    end do
  end do
	  
	  
  do i=1,nz
    z_t(i) = z(i)
  enddo
  
  do i=1,3
    fstart_t(i) = fstart(i)
  enddo
  
	  
  !!pass the pointers
  !double arrays
  ur_f = c_loc(ur_t(1))
  uc_f = c_loc(uc_t(1))

  vr_f = c_loc(vr_t(1))
  vc_f = c_loc(vc_t(1))
  
  wr_f = c_loc(wr_t(1))
  wc_f = c_loc(wc_t(1))
  
  !integer arrays
  z_f      = c_loc(z_t(1))
  fstart_f = c_loc(fstart_t(1))
  
  

end subroutine initializeGPU

subroutine passback()

  use commondata
  use velocity
  use par_size
  use CUDA_modu
  
  integer :: i,j,k,ind,index
  
  !allocate(uur_t(spx*nz*spy))
  !allocate(uvr_t(spx*nz*spy))
  !allocate(uwr_t(spx*nz*spy))
  !allocate(vvr_t(spx*nz*spy))
  !allocate(vwr_t(spx*nz*spy))
  !allocate(wwr_t(spx*nz*spy))
  
  !allocate(uuc_t(spx*nz*spy))
  !allocate(uvc_t(spx*nz*spy))
  !allocate(uwc_t(spx*nz*spy))
  !allocate(vvc_t(spx*nz*spy))
  !allocate(vwc_t(spx*nz*spy))
  !allocate(wwc_t(spx*nz*spy))

  !call c_f_pointer(uur_f, uur_t,[nz*spy*spx])
  !call c_f_pointer(uvr_f, uvr_t,[nz*spy*spx])
  !call c_f_pointer(uwr_f, uwr_t,[nz*spy*spx])
  !call c_f_pointer(vvr_f, vvr_t,[nz*spy*spx])
  !call c_f_pointer(vwr_f, vwr_t,[nz*spy*spx])
  !call c_f_pointer(wwr_f, wwr_t,[nz*spy*spx])
  
  !call c_f_pointer(uuc_f, uuc_t,[nz*spy*spx])
  !call c_f_pointer(uvc_f, uvc_t,[nz*spy*spx])
  !call c_f_pointer(uwc_f, uwc_t,[nz*spy*spx])
  !call c_f_pointer(vvc_f, vvc_t,[nz*spy*spx])
  !call c_f_pointer(vwc_f, vwc_t,[nz*spy*spx])
  !call c_f_pointer(wwc_f, wwc_t,[nz*spy*spx])
    
  ind=1

  do j=1,spy!fpy
    do i=1,spx!nx
      do k=1,nz!fpz
        !uuc(i,k,j,1) = uur_t(ind)
        !uuc(i,k,j,2) = uuc_t(ind)
        
        ind = ind + 1
      enddo
    enddo
  enddo
  
end subroutine passback

subroutine deallocateGPUarrays
  
  use CUDA_modu
  !implicit none
  deallocate(ur_t)
  deallocate(uc_t)
  !deallocate(uout_t)

  deallocate(vr_t)
  deallocate(vc_t)
  !deallocate(vout_t)
  
  deallocate(wr_t)
  deallocate(wc_t)
  !deallocate(wout_t)
  
  
  deallocate(z_t)
  deallocate(fstart_t)
  
end subroutine deallocateGPUarrays