
SUBROUTINE JULIAN_DATE(year,month,day,hour,dt)

! The subroutine JULIAN_DATE updates the julian date
! after timestep dt

implicit none

real(8), parameter:: hour_dur = 3600.d0 ! The duration of hour, s 

!Input variables
real(8)   , intent(in)   :: dt  ! Timestep, s

!Input/output variables
integer(4), intent(inout):: year
integer(4), intent(inout):: month
integer(4), intent(inout):: day

real(8)   , intent(inout):: hour

!Local variables
integer(4) ndaym(12)
data       ndaym /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

SAVE ndaym

if (mod(year,4)==0) then 
  ndaym(2) = 29 !Checking for leap-year
else
  ndaym(2) = 28
endif

hour = hour + dt/hour_dur
if (hour > 24.d0) then
  hour = hour - 24.d0
  day  = day  + 1
  if (day > ndaym(month)) then
    day   = 1
    month = month + 1
    if (month > 12) then
      month = 1
      year  = year + 1
    endif 
  endif
endif

END SUBROUTINE JULIAN_DATE



INTEGER(4) FUNCTION JULIAN_DAY(year,month,day)

implicit none

integer(4), intent(in):: year
integer(4), intent(in):: month
integer(4), intent(in):: day

integer(4) i

integer(4) ndaym(12)
data       ndaym /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

SAVE ndaym

if (mod(year,4)==0) then ! Checking for leap-year
  ndaym(2) = 29
else
  ndaym(2) = 28
endif

if (month == 1) then
  JULIAN_DAY = day
else
  JULIAN_DAY = 0 
  do i = 1, month-1
    JULIAN_DAY = JULIAN_DAY + ndaym(i)
  enddo
  JULIAN_DAY = JULIAN_DAY + day
endif

END FUNCTION JULIAN_DAY



SUBROUTINE DATEMINUS(n,year,month,day,hour,year1,month1,day1,hour1)
use driving_params
implicit none
real(8) hour
integer(4) year,month,day
integer(4) n,ndaym(12)
character year1*4,month1*2,day1*2,hour1*2
data ndaym /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      
SAVE
      
if (n == 1) then
! Minus 1 month      
  if (month == 1) then
    write (month1,'(i2)') 12
    write (year1, '(i4)') year - 1
  else
    write (month1,'(i2)') month - 1
    write (year1, '(i4)') year
  endif 
elseif (n == 2) then
! Minus 1 day      
  if (day == 1) then
    if (month == 1) then
      write (month1,'(i2)') 12
      write (year1, '(i4)') year - 1
      write (day1,  '(i2)') 31
    elseif (mod(year,4) == 0.and.month == 3) then !Checking for leap-year
      write (month1,'(i2)') month - 1
      write (year1, '(i4)') year
      write (day1,  '(i2)') ndaym(month-1) + 1 ! 29 days 
    else
      write (month1,'(i2)') month - 1
      write (year1, '(i4)') year
      write (day1,  '(i2)') ndaym(month-1) 
    endif
  else
    write (month1,'(i2)') month
    write (year1, '(i4)') year
    write (day1,  '(i2)') day - 1 
  endif
elseif (n == 3) then
! Minus 1 hour      
  if (int(hour) == 0) then
    if (day == 1) then
      if (month == 1) then
        write (month1,'(i2)') 12
        write (year1, '(i4)') year-1
        write (day1,  '(i2)') 31
        write (hour1, '(i2)') 23  !24-int(interval)
      elseif (mod(year,4) == 0.and.month == 3) then !Checking for leap-year
        write (month1,'(i2)') month-1
        write (year1, '(i4)') year
        write (day1,  '(i2)') ndaym(month-1) + 1 ! 29 days
        write (hour1, '(i2)') 23 !24-int(interval)
	  else
        write (month1,'(i2)') month-1
        write (year1, '(i4)') year
        write (day1,  '(i2)') ndaym(month-1) 
        write (hour1, '(i2)') 23 !24-int(interval)
      endif
    else
      write (month1,'(i2)') month
      write (year1, '(i4)') year
      write (day1,  '(i2)') day-1
      write (hour1, '(i2)') 23 !24-int(interval)
    endif
  else
    write (month1,'(i2)') month
    write (year1, '(i4)') year
    write (day1,  '(i2)') day
    write (hour1, '(i2)') int(hour)-1 !int(interval)
  endif
else
  print*, 'Invalid meaning of control variable N in &
&  subroutine DATEMINUS: STOP'
  STOP
endif
      
END SUBROUTINE DATEMINUS
      
                  
SUBROUTINE TIMESTR(n,y,m,d,h,timestring)
implicit none
integer n,i
character y*4,m*2,d*2,h*2
character(len=n) timestring 
      
if (n==10) then
  write (timestring, '(a4,3a2)') y,m,d,h 
elseif (n==8) then
  write (timestring, '(a4,2a2)') y,m,d 
elseif (n==6) then
  write (timestring, '(a4,a2)') y,m 
endif 
      
do i = 1, n
  if (timestring(i:i)==' '.or.timestring(i:i)==char(0)) &
&  timestring(i:i)='0'
enddo
      
END SUBROUTINE TIMESTR
