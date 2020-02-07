module CUDA_modu

  use, intrinsic :: iso_c_binding

  integer, parameter :: dp = selected_real_kind(15, 307)

  !!integer controls
  integer(c_int), bind(c) :: spx_f
  integer(c_int), bind(c) :: spy_f
  integer(c_int), bind(c) :: nz_f
  integer(c_int), bind(c) :: nx_f
  
  !!doubles
  real (c_double), bind(c) :: xl_f
  real (c_double), bind(c) :: yl_f
  real (c_double), bind(c) :: dt_f
  
  !!grid size
  type(c_ptr), bind(c) :: z_f
  type(c_ptr), bind(c) :: fstart_f
  
  real(dp), dimension(:), pointer :: z_t
  integer, dimension(:), pointer :: fstart_t
  
  
  !!velocity components
  type(c_ptr), bind(c) :: ur_f
  type(c_ptr), bind(c) :: uc_f
  type(c_ptr), bind(c) :: uout_f
  
  real(dp), dimension(:), pointer :: ur_t
  real(dp), dimension(:), pointer :: uc_t
  real(dp), dimension(:), pointer :: uout_t
  
  
  type(c_ptr), bind(c) :: vr_f
  type(c_ptr), bind(c) :: vc_f
  type(c_ptr), bind(c) :: vout_f

  real(dp), dimension(:), pointer :: vr_t
  real(dp), dimension(:), pointer :: vc_t
  real(dp), dimension(:), pointer :: vout_t
  

  type(c_ptr), bind(c) :: wr_f
  type(c_ptr), bind(c) :: wc_f
  type(c_ptr), bind(c) :: wout_f

  real(dp), dimension(:), pointer :: wr_t
  real(dp), dimension(:), pointer :: wc_t
  real(dp), dimension(:), pointer :: wout_t

  
  
end module CUDA_modu
