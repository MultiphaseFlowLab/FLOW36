subroutine destroy_plan

use fftw3
implicit none


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

call fftw_destroy_plan(plan_x_fwd_2D)
call fftw_destroy_plan(plan_y_fwd_2D)

return
end
