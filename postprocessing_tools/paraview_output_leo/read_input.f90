subroutine read_input

use commondata

integer :: temp_s,temp_e
double precision :: Lx,Ly

 open(unit=66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')

 read(66,*)
 read(66,*)
 read(66,*)
 read(66,*)
 read(66,*)
 read(66,'(i5)') nx
 read(66,'(i5)') ny
 read(66,'(i5)') nz
 read(66,*)
 read(66,'(f16.6)') Re
 read(66,*) ! Co
 read(66,*) ! gradpx
 read(66,*) ! gradpy
 read(66,*) ! CPI flag
 read(66,*) ! repow
 read(66,*)
 read(66,'(f16.8)') Lx
 read(66,'(f16.8)') Ly
 read(66,*)
 read(66,'(i8)') nstart
 read(66,'(i8)') nend
 read(66,'(i8)') ndump
 read(66,'(i8)') sdump
 read(66,*) !dump_failure
 read(66,*) !stat_dump
 read(66,*) !stat_start
 read(66,'(es16.5)') dt
 read(66,*)
 read(66,*) !bc_up
 read(66,*) !bc_low
 read(66,*)
 read(66,*) !phi_flag
 read(66,*) !phicor_flag
 read(66,*) !lamphi
 read(66,*) !match_dens
 read(66,*) !rhor
 read(66,*) !match_visc
 read(66,*) !visr
 read(66,*) !non_newtonian
 read(66,*) !muinfmuzero
 read(66,*) !exp_non_new
 read(66,'(f16.6)') we
 read(66,'(f16.6)') ch
 read(66,*) !pe
 read(66,*) !fr
 read(66,*) !in_cond_phi
 read(66,*) !grav_dir
 read(66,*) !b_type
 read(66,*) !body_flag
 read(66,*) !body_c
 read(66,*) !body_dir
 read(66,*) !sgradp_flag
 read(66,*) !sgradp_dir
 read(66,*) !rep_flag
 read(66,*) !rep_c
 read(66,*)
 read(66,*) !psi_flag
 read(66,*) !Pe_psi
 read(66,*) !Ex
 read(66,*) !P_i
 read(66,'(f16.6)') betas
 read(66,*) !in_cond_psi
 read(66,*)
 read(66,*) !theta_flag
 read(66,*) !Ra
 read(66,*) !Pr
 read(66,*) !p_theta(1)
 read(66,*) !q_theta(1)
 read(66,*) !r_theta(1)
 read(66,*) !p_theta(2)
 read(66,*) !q_theta(2)
 read(66,*) !r_theta(2)
 read(66,*) !in_cond_theta
 read(66,*)
 read(66,*) !part_flag
 read(66,'(i12)') part_number
 read(66,'(i5)') nset
 read(66,*) !in_cond_part_pos
 read(66,*) !in_cond_part_vel
 read(66,'(i8)') part_dump
 read(66,*) !subiteration
 close(66)

 xl=Lx*pi
 yl=Ly*pi

 open(unit=68,file='./input_paraview.f90',form='formatted',status='old',action='read')

 read(68,'(i5)') spectral
 read(68,*)
 read(68,'(i5)') x_start
 read(68,'(i5)') x_end
 read(68,'(i5)') dnx
 read(68,*)
 read(68,'(i5)') y_start
 read(68,'(i5)') y_end
 read(68,'(i5)') dny
 read(68,*)
 read(68,'(i5)') z_start
 read(68,'(i5)') z_end
 read(68,'(i5)') dnz
 read(68,*)
 read(68,'(i8)') temp_s
 read(68,'(i8)') temp_e
 read(68,*)
 read(68,'(i5)') uflag
 read(68,'(i5)') vflag
 read(68,'(i5)') wflag
 read(68,'(i5)') phiflag
 read(68,'(i5)') psiflag
 read(68,'(i5)') tempflag
 read(68,'(i5)') upflag
 read(68,'(i5)') vorflag
 read(68,'(i5)') strflag
 read(68,'(i5)') topflag
 read(68,'(i5)') marflag
 read(68,'(i5)') div2dflag
 read(68,'(i5)') partposflag
 read(68,'(i5)') partvelflag
 read(68,*)
 read(68,'(i5)') exp_x
 read(68,'(i5)') exp_y
 read(68,'(i5)') exp_z

nxf=exp_x*nx
nyf=exp_y*ny
nzf=exp_z*(nz-1)+1

 x_start=max(x_start,1)
 if(x_end.eq.0) x_end=nxf

 y_start=max(y_start,1)
 if(y_end.eq.0) y_end=nyf

 z_start=max(z_start,1)
 if(z_end.eq.0) z_end=nzf

 if(temp_s.ne.-1)then
  nstart=temp_s
 endif

 if(temp_e.ne.-1)then
  nend=temp_e
 endif

 close(68,status='keep')

 if(rank.eq.0) call print_start

return
end










subroutine print_start

use commondata

write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,'(1x,a)') '                   starting post-processing                           '
write(*,'(1x,a,3(i5,a))') 'Nx * Ny * Nz = ',nx,' * ',ny,' * ',nz,' '
write(*,'(1x,a,3(f8.2,a))') 'Lx * Ly * Lz = (',xl,' * ',yl,' * ',2.0d0,')*Re'
write(*,'(1x,a,f8.2)') 'Re = ',re
if(spectral.eq.1)then
 write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',sdump
else
 write(*,'(1x,a,i8,a,i8,a,i5)') 'from ',nstart,' to ',nend,' with step ',ndump
endif
write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,'(1x,a)') '                            fields                                    '
write(*,'(3(a19,i2,3x))') 'u',uflag,'v',vflag,'w',wflag
write(*,'(3(a19,i2,3x))') 'phi',phiflag,'psi',psiflag,'temperature',tempflag
write(*,'(3(a19,i2,3x))') 'up (vec)',upflag,'vor (vec)',vorflag,'strain rate (vec)',strflag
write(*,'(3(a19,i2,3x))') 'topology parameter',topflag,'Marangoni (vec)',marflag,'2D divergence',div2dflag
write(*,'(2(a19,i2,3x))') 'part pos',partposflag,'part vel',partvelflag

write(*,'(1x,a)') '----------------------------------------------------------------------'
write(*,*)

return
end
