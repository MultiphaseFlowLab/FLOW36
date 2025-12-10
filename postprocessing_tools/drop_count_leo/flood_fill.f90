subroutine flood_fill(top, s_drop, id, jd, kd)
  use commondata
  implicit none

  ! Arrays now dimensioned (nx, nz, ny)
  integer :: top(nx,nz,ny), s_drop(nx,nz,ny)
  integer :: id, jd, kd            ! (x, z, y)
  integer :: stack_i(nx*nz*ny)
  integer :: stack_k(nx*nz*ny)
  integer :: stack_j(nx*nz*ny)
  integer :: topi, topk, topj
  integer :: ni, nk, nj
  integer :: sp
  integer :: ddx(18), ddz(18), ddy(18)
  integer :: n

  ! 18-connectivity offsets (unchanged)
  data ddx / &
      1,-1, 0, 0, 0, 0, &       ! faces
      1, 1,-1,-1, 1, 1,-1,-1, & ! edges in xy
      1,-1, 0, 0 /              ! edges xz/yz
  data ddy / &
      0, 0, 1,-1, 0, 0, &       ! faces
      1,-1, 1,-1, 0, 0, 0, 0, & ! edges
      0, 0, 1,-1 /
  data ddz / &
      0, 0, 0, 0, 1,-1, &       ! faces
      0, 0, 0, 0, 1,-1, 1,-1, & ! edges
      1,-1, 1,-1 /

  ! Starting voxel must be pore
  if (top(id, jd, kd) /= 1) return

  ! Initialize DFS stack
  sp = 1
  stack_i(sp) = id        ! x
  stack_k(sp) = jd        ! z (middle index)
  stack_j(sp) = kd        ! y (third index)

  s_drop(id, jd, kd) = 1

  ! Depth-first flood fill
  do while (sp > 0)

    topi = stack_i(sp)
    topk = stack_k(sp)
    topj = stack_j(sp)
    sp = sp - 1

    do n = 1, 18

      ! --- PERIODIC in X (first index) ---
      ni = 1 + mod(topi - 1 + ddx(n) + nx, nx)

      ! --- NON-PERIODIC in Z (second index) ---
      nk = topk + ddz(n)
      if (nk < 1 .or. nk > nz) cycle

      ! --- PERIODIC in Y (third index) ---
      nj = 1 + mod(topj - 1 + ddy(n) + ny, ny)

      ! Check pore & unvisited
      if (top(ni, nk, nj) == 1 .and. s_drop(ni, nk, nj) == 0) then
        s_drop(ni, nk, nj) = 1
        sp = sp + 1
        stack_i(sp) = ni
        stack_k(sp) = nk
        stack_j(sp) = nj
      end if

    end do

  end do

end subroutine flood_fill

