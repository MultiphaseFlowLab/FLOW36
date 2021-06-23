recursive subroutine flood_fill(top,s_drop,id,jd,kd)

use commondata

integer :: top(nx,nz,ny),s_drop(nx,nz,ny)

integer :: id,jd,kd
integer :: inew,jnew

!On commence par actualiser la valeur de s_drop (single drop)
! vaut 1 car contenu dans une droplet
s_drop(id,kd,jd)=1



! search only in x, y and z direction (6 directions)
if(id+1.le.nx)then ! Si l'indice suivant est toujours dans le domaine, alors ...
  ! Si il n'y est plus : else => on va regarder le bords

!Si l'indice suivant appartient à une droplet et n'a toujours pas été actualisé Dans
! s_drop, alors on fait appel au code récursivement
! va actualiser la valeur de s_drop à 1,et regarder pour l'indice d'encore après
! et ainsi de suite
  if((top(id+1,kd,jd).eq.1).and.(s_drop(id+1,kd,jd).eq.0)) call flood_fill(top,s_drop,id+1,jd,kd)

else ! Recherche aux bords
  if((top(1,kd,jd).eq.1).and.(s_drop(1,kd,jd).eq.0)) call flood_fill(top,s_drop,1,jd,kd)
endif

if(id-1.ge.1)then
  if((top(id-1,kd,jd).eq.1).and.(s_drop(id-1,kd,jd).eq.0)) call flood_fill(top,s_drop,id-1,jd,kd)
else
  if((top(nx,kd,jd).eq.1).and.(s_drop(nx,kd,jd).eq.0)) call flood_fill(top,s_drop,nx,jd,kd)
endif

if(jd+1.le.ny)then
  if((top(id,kd,jd+1).eq.1).and.(s_drop(id,kd,jd+1).eq.0)) call flood_fill(top,s_drop,id,jd+1,kd)
else
  if((top(id,kd,1).eq.1).and.(s_drop(id,kd,1).eq.0)) call flood_fill(top,s_drop,id,1,kd)
endif

if(jd-1.ge.1)then
  if((top(id,kd,jd-1).eq.1).and.(s_drop(id,kd,jd-1).eq.0)) call flood_fill(top,s_drop,id,jd-1,kd)
else
  if((top(id,kd,ny).eq.1).and.(s_drop(id,kd,ny).eq.0)) call flood_fill(top,s_drop,id,ny,kd)
endif

if(kd+1.le.nz)then
  if((top(id,kd+1,jd).eq.1).and.(s_drop(id,kd+1,jd).eq.0)) call flood_fill(top,s_drop,id,jd,kd+1)
! no periodicity in z
! else
!   if((top(id,1,jd).eq.1).and.(s_drop(id,1,jd).eq.0)) call flood_fill(top,s_drop,id,jd,1)
endif

if(kd-1.ge.1)then
  if((top(id,kd-1,jd).eq.1).and.(s_drop(id,kd-1,jd).eq.0)) call flood_fill(top,s_drop,id,jd,kd-1)
! no periodicity in z
! else
!   if((top(id,nz,jd).eq.1).and.(s_drop(id,nz,jd).eq.0)) call flood_fill(top,s_drop,id,jd,nz)
endif


! search on diagonals (20 directions)
! search in x+1,y+1,[z-1,z,z+1]  --> 3 cases
inew=mod(id+1,nx)+1
jnew=mod(jd+1,ny)+1
if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
if(kd+1.le.nz)then
 if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
endif
if(kd-1.ge.1)then
 if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
endif

! search in x+1,y-1,[z-1,z,z+1]  --> 3 cases and in x+1,y,[z-1,z+1]  --> 2 cases
jnew=mod(ny+jd-1,ny)+1
if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
if(kd+1.le.nz)then
 if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
 if((top(inew,kd+1,jd).eq.1).and.(s_drop(inew,kd+1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd+1)
endif
if(kd-1.ge.1)then
 if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
 if((top(inew,kd-1,jd).eq.1).and.(s_drop(inew,kd-1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd-1)
endif

! search in x-1,y+1,[z-1,z,z+1]  --> 3 cases
inew=mod(nx+id-1,nx)+1
jnew=mod(jd+1,ny)+1
if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
if(kd+1.le.nz)then
 if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
endif
if(kd-1.ge.1)then
 if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
endif

! search in x-1,y-1,[z-1,z,z+1]  --> 3 cases and in x-1,y,[z-1,z+1]  --> 2 cases
jnew=mod(ny+jd-1,ny)+1
if((top(inew,kd,jnew).eq.1).and.(s_drop(inew,kd,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd)
if(kd+1.le.nz)then
 if((top(inew,kd+1,jnew).eq.1).and.(s_drop(inew,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd+1)
 if((top(inew,kd+1,jd).eq.1).and.(s_drop(inew,kd+1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd+1)
endif
if(kd-1.ge.1)then
 if((top(inew,kd-1,jnew).eq.1).and.(s_drop(inew,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,inew,jnew,kd-1)
 if((top(inew,kd-1,jd).eq.1).and.(s_drop(inew,kd-1,jd).eq.0)) call flood_fill(top,s_drop,inew,jd,kd-1)
endif

! search in x,[y-1,y+1],[z-1,z+1]  --> 4 cases
jnew=mod(jd+1,ny)+1
if(kd+1.le.nz)then
 if((top(id,kd+1,jnew).eq.1).and.(s_drop(id,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd+1)
endif
if(kd-1.ge.1)then
 if((top(id,kd-1,jnew).eq.1).and.(s_drop(id,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd-1)
endif

jnew=mod(ny+jd-1,ny)+1
if(kd+1.le.nz)then
 if((top(id,kd+1,jnew).eq.1).and.(s_drop(id,kd+1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd+1)
endif
if(kd-1.ge.1)then
 if((top(id,kd-1,jnew).eq.1).and.(s_drop(id,kd-1,jnew).eq.0)) call flood_fill(top,s_drop,id,jnew,kd-1)
endif



return
end
