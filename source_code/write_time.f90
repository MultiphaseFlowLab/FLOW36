subroutine write_time(string,t)

implicit none

double precision :: t

character(len=50) :: string

if(t.lt.1.0d0)then
 write(*,'(1x,A,F16.4,A)') trim(string),1000.0d0*t,' milliseconds'
elseif(t.lt.60.0d0)then
 write(*,'(1x,A,F16.4,A)') trim(string),t,' seconds'
elseif(t.lt.3600.0d0)then
 write(*,'(1x,A,F16.4,A)') trim(string),t/60.0d0,' minutes' 
elseif(t.lt.86400.0d0)then
 write(*,'(1x,A,F16.4,A)') trim(string),t/3600.0d0,' hours'
else
 write(*,'(1x,A,F16.4,A)') trim(string),t/86400.0d0,' days'
endif


return
end
