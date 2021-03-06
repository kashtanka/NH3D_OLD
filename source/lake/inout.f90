
function fext2(fname,ext)
character*80 fname,fext2
character*3 ext
parameter(n=80)

fext2=fname

do 10 i=30,1,-1
  if(fext2(i:i).eq.'.') then
    fext2(i:n)=char(0)
    go to 11
 endif
10 continue
11 continue

do 20 i=30,1,-1
  if(fext2(i:i).ne.' ' .and. fext2(i:i).ne.char(0)) then
    fext2(i+1:i+4)='.'//ext
    go to 21
  endif
20 continue
   write(*,*) 'erro em fext'

21 continue
   return
end
                 
subroutine readgrd_lake(nunit,var,i0,i1,j0,j1)
!
! reads 2d array from a grd file
!
implicit real*8(a-h,o-z)
dimension var(i0:i1,j0:j1)
character*4 dsaa

read(nunit,'(a4)') dsaa
read(nunit,*) nx,ny
read(nunit,*) xmin,xmax
read(nunit,*) ymin,ymax
read(nunit,*) zmin,zmax
      
do j=j0,j1
  read(nunit,*) (var(i,j),i=i0,i1)
enddo
      
return
end
      
subroutine wrigrd_lake(nunit,z,dx,dy,i0,i1,j0,j1)

implicit none
integer i0,i1,j0,j1,nunit,i,j
real(8) z(i0:i1,j0:j1)
real(8) xmin,xmax,ymin,ymax,zmin,zmax,dx,dy
      
zmin=z(i0,j0)
zmax=z(i0,j0)

do j=j0,j1
  do i=i0,i1
    zmin=min(zmin,z(i,j))
    zmax=max(zmax,z(i,j))
  enddo
enddo

xmin=float(i0)*dx
xmax=float(i1)*dx
ymin=float(j0)*dy
ymax=float(j1)*dy
      
write(nunit,'(a4)') 'DSAA'
write(nunit,*) (i1-i0+1),(j1-j0+1)
write(nunit,*) xmin,xmax
write(nunit,*) ymin,ymax
write(nunit,*) zmin,zmax

do j=j0,j1
  write(nunit,*) (z(i,j),i=i0,i1)
enddo
     
return
end


SUBROUTINE MON_OUT(ix,iy,nx,ny,year,month,day,hour)
use ARRAYS
implicit none

! Input variables
! Reals
real(8), intent(in) :: hour

! Integers
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny

! Local variables
! Reals
real(8), allocatable :: accum_var (:,:,:,:)
real(8), allocatable :: var       (:,:,:,:)
real(8), allocatable :: tsteps    (:,:)

real(8), external:: DZETA

! Integers
integer(4), allocatable :: month_old(:,:)
integer(4) :: out_unit
integer(4) :: i ! Loop index

! Characters
character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: outpath*60
character :: timestring*6
character :: coords_point*6

! Logicals
logical :: firstcall

common /out/ outpath
data firstcall /.true./
data out_unit /600/

SAVE

if (firstcall) then
  allocate (month_old(1:nx, 1:ny))
  allocate (tsteps    (1:nx, 1:ny))  
  allocate (var      (1:4, 1:nx, 1:ny, 1:M+1) )
  allocate (accum_var(1:4, 1:nx, 1:ny, 1:M+1) )
  month_old(:,:) = month
  tsteps   (:,:) = 0.d0  
  accum_var(:,:,:,:) = 0.d0
  var      (:,:,:,:) = 0.d0
endif

var(1,ix,iy,1:M+1) = Tw1 (1:M+1)
var(2,ix,iy,1:M+1) = Sal1(1:M+1)
var(3,ix,iy,1:M)   = E1  (1:M)
var(4,ix,iy,1:M)   = eps1(1:M)

accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
tsteps(ix,iy) = tsteps(ix,iy) + 1.d0

if (month_old(ix,iy)/=month) then
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:)/tsteps(ix,iy)
  call DATEMINUS (1,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR(6,year1,month1,day1,hour1,timestring)
  call CHECK_UNIT(out_unit)
  write (coords_point,'(2i3)') ix, iy
  open(out_unit,file=outpath(1:len_trim(outpath))//'monthly/'// &
 & 'Profiles'//coords_point//timestring//'.dat', status='unknown')
  write (out_unit,*) '1 - depth, m' 
  write (out_unit,*) '2 - temperature, C'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - turbulent kinetic energy, m**2/s**2'
  write (out_unit,*) '5 - disspation rate, m**2/s**3'
  do i=1,M+1 
    write (out_unit,'(f10.6,f9.3,3e12.4)') &
    & -dzeta(dble(i))*h1, accum_var(1:4,ix,iy,i)
  enddo
  close(out_unit)
  month_old(ix,iy) = month
  tsteps(ix,iy) = 0.d0
  accum_var(:,ix,iy,:) = 0.d0
endif

if (firstcall) firstcall=.false.
END SUBROUTINE MON_OUT


SUBROUTINE DAY_OUT(ix,iy,nx,ny,year,month,day,hour)
use ARRAYS
implicit none

! Input variables
! Reals
real(8), intent(in) :: hour

! Integers
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny

! Local variables
! Reals
real(8), allocatable :: tsteps   (:,:)
real(8), allocatable :: accum_var(:,:,:,:)
real(8), allocatable :: var      (:,:,:,:)

real(8), external:: DZETA

! Integers
integer(4), allocatable :: day_old(:,:)
integer(4) :: i ! Loop index
integer(4) :: out_unit


character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: outpath*60
character :: timestring*8
character :: coords_point*6

logical :: firstcall

common /out/ outpath
data firstcall /.true./
data out_unit /600/
      
SAVE

if (firstcall) then
  allocate (day_old   (1:nx, 1:ny))
  allocate (tsteps    (1:nx, 1:ny))  
  allocate (var       (1:4, 1:nx, 1:ny, 1:M+1) )
  allocate (accum_var (1:4, 1:nx, 1:ny, 1:M+1) )
  day_old  (:,:) = day
  tsteps   (:,:) = 0.d0  
  accum_var(:,:,:,:) = 0.d0
  var      (:,:,:,:) = 0.d0
endif

var(1,ix,iy,1:M+1) = Tw1 (1:M+1)
var(2,ix,iy,1:M+1) = Sal1(1:M+1)
var(3,ix,iy,1:M)   = E1  (1:M)
var(4,ix,iy,1:M)   = eps1(1:M)

accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
tsteps(ix,iy) = tsteps(ix,iy) + 1.d0

if (day_old(ix,iy)/=day) then
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:)/tsteps(ix,iy)
  call DATEMINUS (2,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR (8,year1,month1,day1,hour1,timestring)
  call CHECK_UNIT(out_unit)
  write (coords_point, '(2i3)') ix, iy
  open(out_unit, file=outpath(1:len_trim(outpath))//'daily/'// &
  & 'Profiles'//coords_point//timestring//'.dat',status='unknown')
  write (out_unit,*) '1 - depth, m' 
  write (out_unit,*) '2 - temperature, C'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - turbulent kinetic energy, m**2/s**2'
  write (out_unit,*) '5 - disspation rate, m**2/s**3'
  do i=1,M+1 
    write (out_unit,'(f10.6,f9.3,3e12.4)') &
    &  -dzeta(dble(i))*h1, accum_var(1:4,ix,iy,i)
  enddo
  close(out_unit)
  day_old(ix,iy) = day
  accum_var(:,ix,iy,:) = 0.d0
  tsteps(ix,iy) = 0.d0
endif

if (firstcall) firstcall=.false.
END SUBROUTINE DAY_OUT


SUBROUTINE HOUR_OUT(ix,iy,nx,ny,year,month,day,hour)

use ARRAYS, only:  &
 & Tw1, Sal1,      &
 & E1 , eps1,      &
 & h1 ,            &
 & PEMF,           &
 & PT_DOWN, PSAL_DOWN, &
 & PDENS_DOWN,     &
 & row,            &
 & k_turb_T_flux,  &
 & T_massflux,     & 
 & H_mixed_layer,  & 
 & w_conv_scale,   &
 & T_conv_scale,   &
 & Buoyancy0,      &
 & S,              &
 & Gen,            &
 & TKE_turb_trans, &
 & k3_mid,         &
 & u1,             &
 & v1

use DRIVING_PARAMS

implicit none

! Input variables
! Reals
real(8), intent(in) :: hour

! Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day

! Local variables
! Reals
real(8), allocatable :: accum_var (:,:,:,:)
real(8), allocatable :: var       (:,:,:,:)

real(8), external :: DZETA

! Integers
integer(4), parameter :: n_var = 21
integer(4), allocatable :: hour_old(:,:)
integer(4), allocatable :: tsteps(:,:)
integer(4) :: i
integer(4) :: out_unit

character :: month1*2
character :: year1*4
character :: day1*2
character :: hour1*2
character :: outpath*60
character :: timestring*10
character :: coords_point*6

logical :: firstcall

common /out/ outpath
data firstcall /.true./
data out_unit /600/
      
SAVE

if (firstcall) then
  allocate (accum_var (n_var,1:nx,1:ny,1:M+1) )
  allocate (var       (n_var,1:nx,1:ny,1:M+1) )
  allocate (tsteps    (1:nx,1:ny) )
  allocate (hour_old  (1:nx,1:ny) )
  accum_var(:,:,:,:) = 0.d0
  var(:,:,:,:) = 0.d0
  tsteps(:,:) = 0
  hour_old(:,:) = int(hour)
endif

var(1,ix,iy,1:M)   = PEMF          (1:M)
var(2,ix,iy,1:M+1) = PDENS_DOWN    (1:M+1)
var(3,ix,iy,1:M+1) = PT_DOWN       (1:M+1)
var(4,ix,iy,1:M+1) = PSAL_DOWN     (1:M+1)
var(5,ix,iy,1:M)   = k_turb_T_flux (1:M)
var(6,ix,iy,1:M)   = T_massflux    (1:M)
var(7,ix,iy,1:M+1) = row           (1:M+1)
var(8,ix,iy,1:M+1) = Tw1           (1:M+1)
var(9,ix,iy,1:M+1) = Sal1          (1:M+1)
var(10,ix,iy,1:M)  = E1            (1:M)
var(11,ix,iy,1:M)  = eps1          (1:M)

if (scale_output == 1) then
  var(12,ix,iy,1)    = H_mixed_layer ! Scale
  var(13,ix,iy,1)    = Buoyancy0     ! Scale
  var(14,ix,iy,1)    = w_conv_scale  ! Scale
  var(15,ix,iy,1)    = T_conv_scale  ! Scale
elseif (scale_output == 0) then
  var(12,ix,iy,1)    = 1.d0 ! No scaling
  var(13,ix,iy,1)    = 1.d0 ! No scaling
  var(14,ix,iy,1)    = 1.d0 ! No scaling
  var(15,ix,iy,1)    = 1.d0 ! No scaling
else
  print*, 'Scale_output must be 0 or 1: STOP'  
  STOP
endif
      
var(16,ix,iy,1:M)  = S              (1:M)
var(17,ix,iy,1:M)  = Gen            (1:M)
var(18,ix,iy,1:M)  = TKE_turb_trans (1:M)
      
var(19,ix,iy,2:M)  = k3_mid         (2:M)

var(20,ix,iy,1:M+1)  = u1(1:M+1)
var(21,ix,iy,1:M+1)  = v1(1:M+1)

accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) + var(:,ix,iy,:)
tsteps(ix,iy) = tsteps(ix,iy) + 1

if (hour_old(ix,iy)/=int(hour)) then
  accum_var(:,ix,iy,:) = accum_var(:,ix,iy,:) / real(tsteps(ix,iy))
  call DATEMINUS(3,year,month,day,hour,year1,month1,day1,hour1)
  call TIMESTR(10,year1,month1,day1,hour1,timestring)
  call CHECK_UNIT(out_unit)
  write (coords_point, '(2i3)') ix, iy
  
! Writing to the file Profiles<yyyymmddhh>.dat  
  open(out_unit, file=outpath(1:len_trim(outpath))// &
  & 'hourly/'//'Profiles'//coords_point//timestring//'.dat', &
  &  status='unknown')
  write (out_unit,*) '1 - depth, normalized' 
  write (out_unit,*) '2 - temperature, normalized'
  write (out_unit,*) '3 - salinity, kg/kg' 
  write (out_unit,*) '4 - water density, kg/m**3'  
  write (out_unit,*) '5 - turbulent kinetic energy, normalized'
  write (out_unit,*) '6 - disspation rate, normalized'
  write (out_unit,*) '7 - eddy diffusivity (TKE**2/dissipation), m**2/s'
  write (out_unit,*) '8  - mass flux, m/s'
  write (out_unit,*) '9  - downdraft temperature, C' 
  write (out_unit,*) '10 - downdraft salinity, kg/kg' 
  write (out_unit,*) '11 - downdraft density, kg/m**3'
  write (out_unit,*) '12 - x-component of speed, m/s'
  write (out_unit,*) '13 - y-component of speed, m/s'
  write (out_unit,'(11(i6,4x))') 1,2,3,4,5,6,7,8,9,10,11 
  do i=1,M 
    write (out_unit,'(f14.6,f12.5,e10.3,f10.3,3e10.3,4e14.7,2f11.5)')      &
    & -DZETA(dble(i))*h1                         /accum_var(12,ix,iy,1),  &
    & (accum_var(8,ix,iy,i) - scale_output*maxval(accum_var(8,ix,iy,:)) ) &
    &  /accum_var(15,ix,iy,1), &
    & accum_var(9,ix,iy,i),                                               &
    & accum_var(7,ix,iy,i),                                               &
    & accum_var(10,ix,iy,i)    / (accum_var(13,ix,iy,1)*accum_var(12,ix,iy,1) /       &
    &  accum_var(14,ix,iy,1) ),                                           &
    & accum_var(11,ix,iy,i)    / accum_var(13,ix,iy,1),                   &
    & accum_var(19,ix,iy,i),                                              & !accum_var(10,i)**2 / accum_var(11,i),
    & accum_var(1,ix,iy,i),                                               &
    & accum_var(3,ix,iy,i),                                               &
    & accum_var(4,ix,iy,i),                                               &
    & accum_var(2,ix,iy,i),                                               &
    & accum_var(20,ix,iy,i),                                              &
    & accum_var(21,ix,iy,i)
  enddo 
  close (out_unit)
  
! Writing to the file EDMF_profiles<yyyymmddhh>.dat
  open(out_unit, file=outpath(1:len_trim(outpath))//  &
  & 'hourly/'//'EDMF_profiles'//coords_point//timestring//'.dat',   &
  & status='unknown')
  write (out_unit,*) '1 - depth, normalized by mixed layer depth' 
  write (out_unit,*) '2 - "k - flux" of temperature, normalized'
  write (out_unit,*) '3 - mass flux of temperature,  normalized'
  write (out_unit,*) '4 - total flux of temperature, normalized'
  write (out_unit,'( 4(i10,6x) )') 1,2,3,4
  do i=1, M
    write (out_unit,'(f18.6,3e16.7)')                   &
    & -DZETA(dble(i+0.5))*h1/ accum_var(12,ix,iy,1),          &
    & accum_var(5,ix,iy,i)/(accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1)), &
    & accum_var(6,ix,iy,i)/(accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1)), &
    & (accum_var(5,ix,iy,i) + accum_var(6,ix,iy,i) ) /              &
    & (accum_var(14,ix,iy,1)*accum_var(15,ix,iy,1))
  enddo
  close (out_unit)
  
! Writing to the file TKE_budget<yyyymmddhh>.dat
  open(out_unit, file=outpath(1:len_trim(outpath))//    &
  & 'hourly/'//'TKE_budget'//coords_point//timestring//'.dat',        &
  &  status='unknown')
  write (out_unit,*) '1 - depth, normalised with mixed layer depth'
  write (out_unit,*) '2 - shear production, normalised with surface buoyancy flux'
  write (out_unit,*) '3 - buoyancy source,  normalised with surface buoyancy flux'
  write (out_unit,*) '4 - dissipation rate, normalised with surface buoyancy flux'
  write (out_unit,*) '5 - turbulent transport, normalised with surface buoyancy flux'
  write (out_unit,'( 5(i10,6x) )') 1,2,3,4,5
  do i = 1, M
    write (out_unit,'(f18.6,4e16.7)')          &
    & -dzeta(dble(i+0.5))*h1/ accum_var(12,ix,iy,1), &
    & accum_var(17,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(16,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(11,ix,iy,i) / accum_var(13,ix,iy,1), &
    & accum_var(18,ix,iy,i) / accum_var(13,ix,iy,1)
  enddo
  close (out_unit)
  hour_old(ix,iy) = int(hour)
  accum_var(:,ix,iy,:) = 0.d0
  tsteps(ix,iy)  = 0
endif

if (firstcall) firstcall=.false.
END SUBROUTINE HOUR_OUT          
      


SUBROUTINE EVERYSTEP_OUT(ix,iy,nx,ny)
use ARRAYS
implicit none

! Input variables
! Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny

! Local variables
! Reals
real(8), external :: DZETA
real(8), parameter :: ACC = 1.d-20

! Integers
integer(4), parameter :: nmax_files = 10 ! The maximal number of files
integer(4) :: out_unit(1:nmax_files)
integer(4), allocatable :: n_unit(:,:)
integer(4) :: count_unit
integer(4) :: i ! Loop index

! Characters
character :: outpath*60
character :: coords_point*6

! Logicals
logical, allocatable :: firstcall(:,:)

common /out/ outpath
data out_unit(1) /605/
data count_unit /0/
      
SAVE

if (.not.allocated(firstcall)) then
  allocate (firstcall(1:nx, 1:ny) )
  allocate (n_unit   (1:nx, 1:ny) )
  firstcall = .true.
  n_unit = 0
endif  
      
if (firstcall(ix,iy)) then
  count_unit = count_unit + 1
  if (count_unit > nmax_files) then
    print*, 'Too much output files opened in Lake model: STOP'
    STOP
  endif
  if (count_unit > 1) then
    out_unit(count_unit) = out_unit(count_unit-1) + 1
  endif
  call CHECK_UNIT(out_unit(count_unit))
  write (coords_point, '(2i3)') ix, iy
  open (out_unit(count_unit), file=outpath(1:len_trim(outpath))//'everystep/'// &
  & 'Profiles'//coords_point//'.dat',  status='unknown')
  write (out_unit(count_unit),*) '1 - depth, m' 
  write (out_unit(count_unit),*) '2 - temperature, C'
  write (out_unit(count_unit),*) '3 - salinity, kg/kg' 
  write (out_unit(count_unit),*) '4 - turbulent kinetic energy, m**2/s**2'
  write (out_unit(count_unit),*) '5 - disspation rate, m**2/s**3'
  write (out_unit(count_unit),*) '6 - eddy diffusivity (TKE**2/dissipation), m**2/s'
  write (out_unit(count_unit),*) '7 - k-flux of temperature, K*m/s'
  n_unit(ix,iy) = out_unit(count_unit)
endif
      
write (n_unit(ix,iy),*) 'nstep = ', nstep
do i=1,M 
  write (n_unit(ix,iy), '(f10.6,f9.3,5e12.4)') &
  & -dzeta(dble(i))*h1,Tw1(i),Sal1(i),E1(i),eps1(i),  &
  & E1(i)**2/(eps1(i)+ACC), k_turb_T_flux(i)
enddo
      
if (firstcall(ix,iy)) firstcall(ix,iy)=.false.
END SUBROUTINE EVERYSTEP_OUT
      
      
SUBROUTINE SERIES_OUT(ix,iy,nx,ny,year,month,day,hour,tsw)
use DRIVING_PARAMS
use ARRAYS, only: &
& row,            &
& Tw1,Tskin,      &
& zinv,deltaskin, &
& time,           &
& h1,l1,hs1,ls1,  &
& H_mixed_layer 
use ATMOS,  only:    &
& eflux0_kinem,      &
& turb_density_flux, &
& hw, xlew
implicit none

! Input variables
! Reals
real(8), intent(in) :: tsw
real(8), intent(in) :: hour

! Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day

! Local variables
! Reals
real(8), parameter :: hour_sec = 60.*60.
real(8) :: T_mean

real(8), external :: VARMEAN

! Integers
integer(4), parameter :: nmax_files = 30 ! The maximal number of files
integer(4) :: out_unit(1:nmax_files)
integer(4) :: count_unit
integer(4), allocatable :: count_out(:,:)
integer(4), allocatable :: n_unit(:,:,:)

! Characters
character :: outpath*60
character :: coords_point*6

! Logicals
logical, allocatable :: firstcall(:,:)

common /out/ outpath
data count_unit /0/
data out_unit(1) /650/
      
SAVE

if (.not.allocated(firstcall)) then
  allocate (firstcall(1:nx, 1:ny) )
  allocate (count_out(1:nx, 1:ny) )
  allocate (n_unit   (1:3, 1:nx, 1:ny) )
  firstcall(:,:) = .true.
!  n_unit = 0
endif  
            
if (firstcall(ix,iy)) then

  write (coords_point, '(2i3)') ix, iy
  
  count_unit = count_unit + 1
  if (count_unit > 1) then
    out_unit(count_unit) = out_unit(count_unit-1) + 1
  endif
  call CHECK_UNIT(out_unit(count_unit))
    
  open (out_unit(count_unit),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'layers'//coords_point//'.dat',  status='unknown')
  write (out_unit(count_unit),*)'Col. 1 - year'
  write (out_unit(count_unit),*)'Col. 2 - month'
  write (out_unit(count_unit),*)'Col. 3 - day'
  write (out_unit(count_unit),*)'Col. 4 - hour'
  write (out_unit(count_unit),*)'Col. 5 - the time from the start of integration, hours'
  write (out_unit(count_unit),*)'Col. 6 - water layer thickness, m'
  write (out_unit(count_unit),*)'Col. 7 - ice layer thickness,   m'
  write (out_unit(count_unit),*)'Col. 8 - snow layer thickness,  m'
  write (out_unit(count_unit),*)'Col. 9 - bottom ice thickness,  m'
  n_unit(1,ix,iy) = out_unit(count_unit)

  count_unit = count_unit + 1
  out_unit(count_unit) = out_unit(count_unit-1) + 1
  call CHECK_UNIT(out_unit(count_unit))
  
  open (out_unit(count_unit),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'T_fluxes'//coords_point//'.dat',status='unknown')
  write (out_unit(count_unit),*)'Col. 1 - year'
  write (out_unit(count_unit),*)'Col. 2 - month'
  write (out_unit(count_unit),*)'Col. 3 - day'
  write (out_unit(count_unit),*)'Col. 4 - hour'
  write (out_unit(count_unit),*)'Col. 5 - the time from the start of integration, hours'
  write (out_unit(count_unit),*)'Col. 6 - surface temperature, C'
  write (out_unit(count_unit),*)'Col. 7 - skin temperature, C'
  write (out_unit(count_unit),*)'Col. 8 - skin thickness, m'
  write (out_unit(count_unit),*)'Col. 9 - bottom temperature of water, C'
  write (out_unit(count_unit),*)'Col. 10 - mean temperature of water coloumn, C'
  write (out_unit(count_unit),*)'Col. 11 - sensible heat flux,    W/m**2'
  write (out_unit(count_unit),*)'Col. 12 - latent heat flux,     W/m**2'
  n_unit(2,ix,iy) = out_unit(count_unit)

  count_unit = count_unit + 1
  out_unit(count_unit) = out_unit(count_unit-1) + 1
  call CHECK_UNIT(out_unit(count_unit))

  open (out_unit(count_unit),file=outpath(1:len_trim(outpath))//'time_series/'// &
  & 'conv_series'//coords_point//'.dat',status='unknown')
  write (out_unit(count_unit),*)'Col. 1 - year'
  write (out_unit(count_unit),*)'Col. 2 - month'
  write (out_unit(count_unit),*)'Col. 3 - day'
  write (out_unit(count_unit),*)'Col. 4 - hour'
  write (out_unit(count_unit),*)'Col. 5 - the time from the start of integration, hours'
  write (out_unit(count_unit),*)'Col. 6 - surface temperature, C'
  write (out_unit(count_unit),*)'Col. 7 - surface density, kg/m**3'
  write (out_unit(count_unit),*)'Col. 8 - heat flux downwards at the surface, K*m/s'
  write (out_unit(count_unit),*)'Col. 9 - turbulent density flux at the surface, kg/m**2/s'
  write (out_unit(count_unit),*)'Col. 10 - inversion depth (mass flux diagnostics),    m'
  write (out_unit(count_unit),*)'Col. 11 - mixed layer depth,                          m'
  n_unit(3,ix,iy) = out_unit(count_unit)
  
  if (count_unit > nmax_files) then
    print*, 'Too much output files is opened in SERIES_OUT of Lake model: STOP'
    STOP
  endif
  
  count_out(ix,iy) = int(time/(dt_out*hour_sec))
endif
      
if (int(time/(dt_out*hour_sec))>=count_out(ix,iy)) then
  T_mean = varmean(Tw1(1:M+1),1)
  write (n_unit(1,ix,iy),'(3i6,f7.2,f10.2,4f9.4)') year,month,day ,hour, &
  &                                      time/hour_sec ,            &
  &                                      h1  ,l1   ,hs1 ,ls1        
  write (n_unit(2,ix,iy),'(3i6,f7.2,f10.2,2f9.4,f9.6,4f9.4)')   year,month,day,hour, &
  &                                      time/hour_sec ,                &
  &                                      tsw-273.15,Tskin(1),deltaskin,Tw1(M+1),  &
  &                                      T_mean,hw,xlew
  write (n_unit(3,ix,iy),'(3i6,f7.2,f10.2,2f9.4,2e15.5,2f9.4)')      &
  &                                      year,month,day,hour,        &
  &                                      time/hour_sec,              &
  &                                      Tw1(1),row(1),eflux0_kinem, &
  &                                      turb_density_flux,zinv,H_mixed_layer
  count_out(ix,iy) = count_out(ix,iy) + 1 
endif 
      
if (firstcall(ix,iy)) firstcall(ix,iy)=.false.
END SUBROUTINE SERIES_OUT

      

SUBROUTINE CHECK_UNIT(nunit)
implicit none

! Input variables
integer(4), intent(inout) :: nunit

! Local variables
logical :: unit_opened

inquire (unit=nunit,opened=unit_opened)

if (unit_opened) then
  do while (unit_opened)
    print*, 'The unit ', nunit, 'is attempted &
    & to be connected to a file, while already connected: incrementing unit'
    print*, 'New unit number is ', nunit + 1
    nunit = nunit + 1
    inquire (unit=nunit,opened=unit_opened)
  enddo
!  STOP
endif

END SUBROUTINE CHECK_UNIT


FUNCTION GETVARVAL(n1,n2,line,name)
implicit none

real(8) :: GETVARVAL

! Input variables
integer,       intent(in):: n2,n1
character*200, intent(in):: line
character*25,  intent(in):: name

! Local variables
real(8) work

read (line((n2+1):100),*) work
print*, name//' = ', work

GETVARVAL = work      

RETURN
END FUNCTION GETVARVAL


FUNCTION IGETVARVAL(n1,n2,line,name)
implicit none

integer(4) :: IGETVARVAL

! Input variables
integer(4),    intent(in):: n2,n1
character*200, intent(in):: line
character*25,  intent(in):: name

! Local variables
integer(4) iwork

read (line((n2+1):100),*) iwork
print*, name//' = ', iwork

IGETVARVAL = iwork      

RETURN
END FUNCTION IGETVARVAL
