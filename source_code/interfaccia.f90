module interfaccia
!*********************************************************************
!
! Interface to call Cuda C kernels in Fortran to perform backwards 
! and forward transformations;
! Fortran is column major so the memory is stored from left to right
! in the order of the dimensions;
!
! copyright  Multiphase Flow Laboratory, University of Udine 
! author     D. Di Giusto jan 2020
! 
!*********************************************************************
  interface
     subroutine h_initialize_gpu(spx_f,spy_f,nz_f,nx_f,f_handle) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int),value :: spx_f
       integer(c_int),value :: spy_f
       integer(c_int),value :: nz_f
       integer(c_int),value :: nx_f
       integer(c_int) :: f_handle
     end subroutine h_initialize_gpu
     
     subroutine h_free_gpu() bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
     end subroutine h_free_gpu
     
     subroutine h_chebyshev_back(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_chebyshev_back

     subroutine h_ffty_back(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_ffty_back

     subroutine h_fftx_back(in_r,in_c,out_r,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       integer(c_int), value :: aliasing
     end subroutine h_fftx_back
     
     subroutine h_fftx_fwd(in_r,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_fftx_fwd
     
     subroutine h_ffty_fwd(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_ffty_fwd     
     
     
     subroutine h_chebyshev_fwd(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_chebyshev_fwd
     
     subroutine h_ffty_many_bwd(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_ffty_many_bwd
     
     subroutine h_ffty_many_fwd(in_r,in_c,out_r,out_c,aliasing) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), dimension(*) :: in_r
       real(c_double), dimension(*) :: in_c
       real(c_double), dimension(*) :: out_r
       real(c_double), dimension(*) :: out_c
       integer(c_int), value :: aliasing
     end subroutine h_ffty_many_fwd
     
  end interface     
end module interfaccia