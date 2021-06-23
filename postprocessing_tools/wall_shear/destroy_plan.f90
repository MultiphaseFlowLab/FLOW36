subroutine destroy_plan

use fftw3
implicit none


call fftw_destroy_plan(plan_x_fwd)
call fftw_destroy_plan(plan_y_fwd)
call fftw_destroy_plan(plan_z_fwd)

call fftw_destroy_plan(plan_x_bwd)
call fftw_destroy_plan(plan_y_bwd)
call fftw_destroy_plan(plan_z_bwd)


return
end
