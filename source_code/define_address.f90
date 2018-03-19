subroutine define_address

use commondata
use dual_grid

integer :: i,j,ii,jj,cnt,k
integer :: snd,rcv,db
character(len=40) :: fmt


! coarse to fine addresses
c2fadd=-1
cnt=1

do j=1,nzcpu
  do i=1,nycpu
    ii=mod(rank,nycpu)+1
    jj=floor(real(rank)/real(nycpu))+1
    ! rank is sender
    if((cg_size(ii,jj,4).le.fg_size(i,j,4)+fg_size(i,j,2)-1).and. &
     & (cg_size(ii,jj,4)+cg_size(ii,jj,2)-1.ge.fg_size(i,j,4))) then
      if((cg_size(ii,jj,5)+cg_size(ii,jj,3).le.ny/2+1).and. &
       & (cg_size(ii,jj,5).le.fg_size(i,j,5)+fg_size(i,j,3)-1).and. &
       & (cg_size(ii,jj,5)+cg_size(ii,jj,3)-1.ge.fg_size(i,j,5))) then
        c2fadd(cnt,1)=rank
        c2fadd(cnt,2)=nycpu*(j-1)+(i-1)
        cnt=cnt+1
        ! write(*,*) 'cond 1  ',rank,' to ',nycpu*(j-1)+(i-1)
      elseif((cg_size(ii,jj,5).gt.ny/2+1).and. &
       & (npsiy-ny+cg_size(ii,jj,5).le.fg_size(i,j,5)+fg_size(i,j,3)-1).and. &
       & (npsiy-ny+cg_size(ii,jj,5)+cg_size(ii,jj,3)-1.ge.fg_size(i,j,5))) then
        c2fadd(cnt,1)=rank
        c2fadd(cnt,2)=nycpu*(j-1)+(i-1)
        cnt=cnt+1
        ! write(*,*) 'cond 2  ',rank,' to ',nycpu*(j-1)+(i-1)
      elseif(((cg_size(ii,jj,5).le.ny/2+1).and.(cg_size(ii,jj,5)+cg_size(ii,jj,3).gt.ny/2+1)))then
        if((cg_size(ii,jj,5).le.fg_size(i,j,5)+fg_size(i,j,3)-1).and. &
          & (ny/2+1.ge.fg_size(i,j,5))) then
          c2fadd(cnt,1)=rank
          c2fadd(cnt,2)=nycpu*(j-1)+(i-1)
          cnt=cnt+1
          ! write(*,*) 'cond 3.1',rank,' to ',nycpu*(j-1)+(i-1)
        endif
        if((npsiy-ny/2+2.le.fg_size(i,j,5)+fg_size(i,j,3)-1).and. &
          & (npsiy-ny+cg_size(ii,jj,5)+cg_size(ii,jj,3)-1.ge.fg_size(i,j,5))) then
          ! check for no double
          snd=rank
          rcv=nycpu*(j-1)+(i-1)
          db=0
          do k=1,cnt
            if((c2fadd(k,1).eq.snd).and.(c2fadd(k,2).eq.rcv)) then
              db=1
              exit
            endif
          enddo
          if(db.eq.0)then
            c2fadd(cnt,1)=rank
            c2fadd(cnt,2)=nycpu*(j-1)+(i-1)
            cnt=cnt+1
            ! write(*,*) 'cond 3.2',rank,' to ',nycpu*(j-1)+(i-1)
          endif
        endif
      endif
    endif
    ! rank is receiver, check to avoid including addresses already included
    if((cg_size(i,j,4).le.fg_size(ii,jj,4)+fg_size(ii,jj,2)-1).and. &
     & (cg_size(i,j,4)+cg_size(i,j,2)-1.ge.fg_size(ii,jj,4))) then
      if((cg_size(i,j,5)+cg_size(i,j,3).le.ny/2+1).and. &
       & (cg_size(i,j,5).le.fg_size(ii,jj,5)+fg_size(ii,jj,3)-1).and. &
       & (cg_size(i,j,5)+cg_size(i,j,3)-1.ge.fg_size(ii,jj,5))) then
        snd=nycpu*(j-1)+(i-1)
        rcv=rank
        db=0
        do k=1,cnt
          if((c2fadd(k,1).eq.snd).and.(c2fadd(k,2).eq.rcv)) then
            db=1
            exit
          endif
        enddo
        if(db.eq.0)then
          c2fadd(cnt,1)=nycpu*(j-1)+(i-1)
          c2fadd(cnt,2)=rank
          cnt=cnt+1
          ! write(*,*) 'cond 1  ',nycpu*(j-1)+(i-1),' to ',rank
        endif
      elseif((cg_size(i,j,5).gt.ny/2+1).and. &
       & (npsiy-ny+cg_size(i,j,5).le.fg_size(ii,jj,5)+fg_size(ii,jj,3)-1).and. &
       & (npsiy-ny+cg_size(i,j,5)+cg_size(i,j,3)-1.ge.fg_size(ii,jj,5))) then
        snd=nycpu*(j-1)+(i-1)
        rcv=rank
        db=0
        do k=1,cnt
          if((c2fadd(k,1).eq.snd).and.(c2fadd(k,2).eq.rcv)) then
            db=1
            exit
          endif
        enddo
        if(db.eq.0)then
          c2fadd(cnt,1)=nycpu*(j-1)+(i-1)
          c2fadd(cnt,2)=rank
          cnt=cnt+1
          ! write(*,*) 'cond 2  ',nycpu*(j-1)+(i-1),' to ',rank
        endif
      elseif(((cg_size(i,j,5).le.ny/2+1).and.(cg_size(i,j,5)+cg_size(i,j,3).gt.ny/2+1)))then
        if((cg_size(i,j,5).le.fg_size(ii,jj,5)+fg_size(ii,jj,3)-1).and. &
          & (ny/2+1.ge.fg_size(ii,jj,5))) then
          snd=nycpu*(j-1)+(i-1)
          rcv=rank
          db=0
          do k=1,cnt
            if((c2fadd(k,1).eq.snd).and.(c2fadd(k,2).eq.rcv)) then
              db=1
              exit
            endif
          enddo
          if(db.eq.0)then
            c2fadd(cnt,1)=nycpu*(j-1)+(i-1)
            c2fadd(cnt,2)=rank
            cnt=cnt+1
            ! write(*,*) 'cond 3.1  ',nycpu*(j-1)+(i-1),' to ',rank
          endif
        endif
        if((npsiy-ny/2+2.le.fg_size(ii,jj,5)+fg_size(ii,jj,3)-1).and. &
          & (npsiy-ny+cg_size(i,j,5)+cg_size(i,j,3)-1.ge.fg_size(ii,jj,5))) then
          snd=nycpu*(j-1)+(i-1)
          rcv=rank
          db=0
          do k=1,cnt
            if((c2fadd(k,1).eq.snd).and.(c2fadd(k,2).eq.rcv)) then
              db=1
              exit
            endif
          enddo
          if(db.eq.0)then
            c2fadd(cnt,1)=nycpu*(j-1)+(i-1)
            c2fadd(cnt,2)=rank
            cnt=cnt+1
            ! write(*,*) 'cond 3.2  ',nycpu*(j-1)+(i-1),' to ',rank
          endif
        endif
      endif
    endif
  enddo
enddo

! format string for address printing
write(fmt,'(a,i0,a,i0,a)') '(i4,3x,',ntask,'(i4),4x,',ntask,'(i4))'


! write 1st line sender, 2nd line receiver
! write(*,fmt) rank,c2fadd


! fine to coarse addresses
! inverse mapping fine to coarse
f2cadd(:,1)=c2fadd(:,2)
f2cadd(:,2)=c2fadd(:,1)

! write 1st line, 2nd line
! write(*,fmt) rank,f2cadd

return
end
