subroutine get_interface(nstep)

use commondata
use phase_field
use sim_par
use velocity

integer :: nstep
integer :: i,j,k,id,jd,kd
integer :: top(nx,nz,ny),s_drop(nx,nz,ny)
integer :: drop_count

character(len=40) :: namedir
character(len=8) :: numfile

logical :: check

namedir='../results/'
write(numfile,'(i8.8)') nstep


inquire(file=trim(namedir)//'phi_'//numfile//'.dat',exist=check)

if(check.eqv..true.)then
 open(75,file='./output/center_timestep_'//numfile//'.dat',form='formatted',status='new',access='append')
 write(75,'(8(a8,6x))') 'drop','xG','yG','zG','uG','vG','wG','volume'
 close(75,status='keep')

 !On accède aux données de position obtenues dans le dossier Results
 !open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='big_endian')
 open(66,file=trim(namedir)//'phi_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(66) phi
 close(66,status='keep')
! write(*,*)maxval(phi)

 !On accède aux données de vitesse obtenues dans le dossier results
 open(80,file=trim(namedir)//'u_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(80) u
 close(80,status='keep')
! write(*,*)maxval(u)

 open(81,file=trim(namedir)//'v_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(81) v
 close(81,status='keep')
! write(*,*)maxval(v)

 open(81,file=trim(namedir)//'w_'//numfile//'.dat',access='stream',status='old',form='unformatted',convert='little_endian')
 read(81) w
 close(81,status='keep')
! write(*,*)maxval(w)

 ! Cette partie change la matrice 3D de nombres réels en une matrice 3D d'entiers
 ! Dans la matrice d'etiers : si =1 : est dans la droplet , =0 : en dehors de la droplet

 do j=1,ny
   do k=1,nz
     do i=1,nx
       if(phi(i,k,j).ge.0.0d0)then
         top(i,k,j)=1
       else
         top(i,k,j)=0
       endif
     enddo
   enddo
 enddo

 ! Fin de conversion de la matrice de réels en integers
 ! top(i,k,j)=1 veut dire que le point (i,k,j) appartient à la droplet
 ! top(i,k,j)=0 veut dire que le point (i,k,j) n'appartient pas à la droplet
 !!! top est la matrice 3D qui permet de dire si un élément est dans une droplet ou non

 drop_count=0

 ! flood fill algorithm
 do jd=1,ny
  do kd=1,nz
   do id=1,nx
    if(top(id,kd,jd).gt.0)then
     drop_count=drop_count+1
     write(*,'(2x,a,i10,a,i3,a)') 'Step ',nstep,', New drop, ',drop_count,' drops'
     ! single drop part
     s_drop=0 ! on réinitialise la matrice single droplet avec des zéros
     s_drop(id,kd,jd)=1 ! 1er indice dans la droplet, on actualise la matrice single droplet
     ! flood fill algorithm
     ! grâce à l'appel récursif de flood_fill, on remplit entièrement la matrice single droplet
     ! en explorant tous les voisins à partir du premier indice appartenant à la droplet
     ! on obtient alors une matrice avec des 0 partout sauf au niveau de la dite droplet ( 1 )
     ! La différence avec top, c'est que single droplet ne contient qu'une droplet !
     call flood_fill(top,s_drop,id,jd,kd)
     ! remove drops already done from top
     ! on enlève la droplet simple qu'on vient de traiter à top
     ! top contient donc au début de l'algo toutes les droplets, puis de moins en moins
     ! au fur et à mesure que celles-ci sont trouvées et comptabilisées
     top=top-s_drop
     !! ici qu'il faut code
     ! calculate drop volume and center of mass
     call get_center(s_drop,nstep,drop_count)
     ! new drop calculation
    endif
   enddo
  enddo
 enddo

 ! On ouvre le fichier créé dans le main
 ! ce fichier contient le nb total de droplet à chaque timestep
 open(42,file='./output/drop_count.dat',access='append',form='formatted',status='old')
 ! On vient de calculer le nb total de droplet pour le timestep donné
 ! on actualise donc le fichier
 write(42,'(i16,2x,es16.6,2x,i16)') nstep,dble(nstep)*dt*Re,drop_count
 ! on ferme le fichier
 close(42,status='keep')

 !close(70,status='keep')
 !write(*,'(2x,a,i4)') 'Number of drops: ',drop_count
 !write(*,*)
endif


return
end
