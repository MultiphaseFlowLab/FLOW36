subroutine get_velocity

use mpi
use commondata
use particle


! get velocity, only ranks in comm_comm are involved
if(leader.eq.0)then


else


endif




return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_2wforces




return
end
