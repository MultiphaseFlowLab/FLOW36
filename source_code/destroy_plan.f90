subroutine destroy_plan

#define openaccflag openacccompflag
#if openaccflag == 0
use fftw3
#endif
#if openaccflag == 1
use cufft
use openacc
use cufftplans
#endif

implicit none

#if openaccflag == 0
call fftw_destroy_plan(plan_x_fwd)
call fftw_destroy_plan(plan_y_fwd)
call fftw_destroy_plan(plan_z_fwd)

call fftw_destroy_plan(plan_x_bwd)
call fftw_destroy_plan(plan_y_bwd)
call fftw_destroy_plan(plan_z_bwd)

call fftw_destroy_plan(plan_x_fwd_fg)
call fftw_destroy_plan(plan_y_fwd_fg)
call fftw_destroy_plan(plan_z_fwd_fg)

call fftw_destroy_plan(plan_x_bwd_fg)
call fftw_destroy_plan(plan_y_bwd_fg)
call fftw_destroy_plan(plan_z_bwd_fg)
#endif


#if openaccflag == 1
gerr=gerr+cufftDestroy(cudaplan_x_fwd)
gerr=gerr+cufftDestroy(cudaplan_y_fwd)
gerr=gerr+cufftDestroy(cudaplan_z_fwd)

gerr=gerr+cufftDestroy(cudaplan_x_bwd)
gerr=gerr+cufftDestroy(cudaplan_y_bwd)
gerr=gerr+cufftDestroy(cudaplan_z_bwd)

gerr=gerr+cufftDestroy(cudaplan_x_fwd_fg)
gerr=gerr+cufftDestroy(cudaplan_y_fwd_fg)
gerr=gerr+cufftDestroy(cudaplan_z_fwd_fg)

gerr=gerr+cufftDestroy(cudaplan_x_bwd_fg)
gerr=gerr+cufftDestroy(cudaplan_y_bwd_fg)
gerr=gerr+cufftDestroy(cudaplan_z_bwd_fg)

gerr=gerr+cufftDestroy(cudaplan_z_fwd_1d)
gerr=gerr+cufftDestroy(cudaplan_z_bwd_1d)
if (gerr.ne.0) write(*,*) "Error in cufftDestroy"
#endif

return
end



