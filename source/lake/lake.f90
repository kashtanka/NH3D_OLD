
SUBROUTINE INIT_LAKE(nx,ny,fnmap1,dt, &
 & h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,&
     &         init_T_1	  )


!The subroutine INIT_LAKE implements
!initialisation of the model (but not initial conditions)

use PHYS_PARAMETERS
use PHYS_CONSTANTS 
use NUMERIC_PARAMS
use PHYS_CONSTANTS2
use DRIVING_PARAMS 
use ATMOS
use ARRAYS
implicit none
real(8),intent(out)::h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1 
	 integer(8),intent(out):: init_T_1

!Input variables
!Reals
real(8), intent(in) :: dt

!Integers
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny

!Characters
character(len=*), intent(in) :: fnmap1

!Output variables

!Local variables
!Reals
real(8), parameter :: hour_sec = 60.*60.
real(8) :: x1

!Integers
integer(4) :: n_unit
integer(4) :: i
integer(4) :: j

!Logicals
logical :: uniform_depth

!Characters
character(len=80) :: path2
character(len=80) :: fnmap

data uniform_depth /.true./
data n_unit /999/

!External procedures
character(len=80), external:: fext2

fnmap = fnmap1

call DEFINE_PHYS_PARAMETERS
call DEFINE_PHYS_CONSTANTS
call DEFINE_PHYS_CONSTANTS2
call DEFINE_DRIVING_PARAMS( &
 & h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,&
     &         init_T_1	  )
call ALLOCARRAYS(nx,ny)
call PBLDAT
call COMSOILFORLAKE

if (error_cov==1) then
  if (runmode==2) then
    print*, 'The error covariance estimation algorithm works &
    &only in standalone runs: STOP'
    STOP
  endif
! if (spinup_day/=0) then
!  print*, 'No spinup is allowed when calculating 
!& error covariances: STOP'
!  STOP
! endif 
  if (assim == 1) then
    print*, 'assim = 1 and error_cov = 1: these regimes could not &
    &be turned on simultaneously: STOP'
    STOP
  endif
  if (nx>1.or.ny>1) then
    print*, 'Error covariance calculation is currently adjusted &
    &only for one-point stand-alone runs of the lake model: STOP'
    STOP
  endif
endif
if (assim==1.and.(nx>1.or.ny>1)) then
  print*, 'Data assimilation is currently adjusted &
  &only for one-point stand-alone lake model: STOP'
  STOP
endif

if (dt_out<dt/hour_sec) then
  print*, 'dt_out must be larger or equal to timestep: STOP'
  STOP
endif
   
if (runmode == 2.and.(.not.uniform_depth)) then
  path2 = fext2(fnmap,'dep')
  call CHECK_UNIT(n_unit)
  open (n_unit, file=path(1:len_trim(path))// &
  & path2(1:len_trim(path2)), status='old')
  call READGRD_LAKE(n_unit,dep_2d,0,nx-1,0,ny-1)
  close (n_unit)
  dep_av = 0.
  x1 = 0.
  do i = 1, nx-1
    do j = 1, ny-1
      if (dep_2d(i,j) > 0.) then
        dep_av = dep_av + dep_2d(i,j)
        x1 = x1 + 1.
      endif
    enddo
  enddo
  dep_av = dep_av / x1
endif

print*, 'Lake model is initialized'

END SUBROUTINE INIT_LAKE


SUBROUTINE LAKE &
& (tempair1,humair1,pressure1,uwind1,vwind1, &
& longwave1,shortwave1,precip1,Sflux01,zref1,dtl, &
& h10, l10, ls10, hs10, Ts0, Tb0, Tbb0, h_ML0, extwat, extice, &
& kor_in, trib_inflow, Sals0, Salb0, fetch, phi, lam, us0, vs0, &
& Tm, alphax, alphay, a_veg, c_veg, h_veg, area_lake, &
& tsw,hw1,xlew1,cdmw1, &
& ix,iy,nx,ny,year,month,day,hour,init_T,flag_assim,flag_print, &
& outpath1)

!MAIN SUBROUTINE OF THE MODEL 
!Calculates all the parameters of state of a lake (temperature, currents,
!eddy diffusivity coefficients, snow, ice, soil characteristics, etc.) 
!at the next timestep (i.e. implements transition of variables 
!from i-th time step to (i+1)-th )

!In atmospheric model must be called once each timestep,
!or once each N timesteps, where N = dt_lake/dt_atmos 
!(dt_lake - timestep in the lake model, dt_atmos - timestep the atmospheric model) 

!INPUT VARIABLES:
!tempair1--- air temperature, K;
!humair1 --- specific humidity of air, kg/kg;
!pressure1--- air pressure, Pa;
!uwind1  --- x-component of wind, m/s;
!vwind1  --- y-component of wind, m/s;
!longwave1--- longwave downward radiation, W/m**2;
!shortwave1--- net solar radiation, W/m**2;
!precip1    --- precipitation intensity, m/s;
!xref1   --- the height (m) of level in atmosphere, where temperature, 
!humidity and wind are measured (calculated by atmospheric model);
!trib_inflow1  --- tributary inflow rate, m/s;
!dtl--- timestep (for stability must be not larger than 10-15 min), s;
!ix--- current number of grid point in X-direction;
!iy--- current number of grid point in Y-direction;
!nx--- total number of grid points in X-direction of atmospheric model;
!ny--- total number of grid points in Y-direction of atmospheric model;

!OUTPUT VARIABLES:   
!tsw--- surface temperature of a lake, K;
!hw1--- sensible heat flux from lake, W/m**2; 
!xlew1   --- latent heat flux from lake, W/m**2; 
!cdmw    --- exchange coefficient, m/s ( == wind_speed*exchange_nondimensional_coeficient)

use NUMERIC_PARAMS
use DRIVING_PARAMS
use ATMOS
use PHYS_CONSTANTS2
use ARRAYS
use ALLOC, only: &
& watice
use ASSIM_VAR
use MODI_MASSFLUX_CONVECTION
use CONVECTPAR
use WATER_DENSITY, only: &
& WATER_DENS_TS
use PHYS_FUNC, only: &
& TURB_SCALES, &
& MELTPNT, &
& WATER_FREEZE_MELT, &
& TURB_DENS_FLUX, &
& SINH0, &
& WATER_ALBEDO, &
& EXTINCT_SNOW, &
& UNFRWAT, &
& WI_MAX, &
& WL_MAX

implicit none


!------------------------ MAIN VARIABLES ------------------------
!  arrays:    Tw1 and Tw2 - temperature of water, C
! Ti1 and Ti2 - temperature of ice, C
! T - temperature of snow, C
!  functions of time:h1 and h2 - thickness of water, m
! l1 and l2 - thickness of ice, m
! hs1 - thickness of snow, m  
! flag_ice_surf - shows if ice is melting on the top surface, n/d  
!  constants: cw - specific heat of water, J/(kg*K) 
! ci - specific heat of ice, J/(kg*K)
! row0 - density of water, kg/m**3
! roi - density of ice, kg/m**3
! lamw - thermal conductivity of water, J*sec/(m**3*K)
! lami - thermal conductivity of ice, J*sec/(m**3*K)
! L - specific heat of melting, J/kg
! Lwv - specific heat of evaporation, J/kg
! Liv - specific heat of sublimation, J/kg    
! ddz - size of grid element (space), m
! dt - size of grid element (time), sec
! kstratwater - coefficient for linear initial conditions in water
! kstratice - coefficient for linear initial conditions in ice 
!  boundary conditions:
! eFlux - total heat flux on the upper surface,
!   shortwave(1-A)-T**4+longwave-H-LE, J/(m**2*sec)
! Elatent - latent heat flux, J/(m**2*sec) 
! Erad - shortwave radiation, penetrated below a surface, J/(m**2*sec)
! tempair - air temperature (2 m), C
! precip - precipitation, m/sec
!(M) number of layers in water and snow   

!Parameters

integer(4), parameter :: vector_length = 350

!Input variables
!Reals
real(8), intent(in) :: tempair1
real(8), intent(in) :: humair1
real(8), intent(in) :: pressure1
real(8), intent(in) :: uwind1
real(8), intent(in) :: vwind1
real(8), intent(in) :: longwave1
real(8), intent(in) :: zref1
real(8), intent(in) :: shortwave1
real(8), intent(in) :: precip1
real(8), intent(in) :: Sflux01
real(8), intent(in) :: dtl
real(8), intent(in) :: hour

real(8), intent(in) :: h10
real(8), intent(in) :: l10
real(8), intent(in) :: ls10
real(8), intent(in) :: hs10
real(8), intent(in) :: Ts0
real(8), intent(in) :: Tb0
real(8), intent(in) :: Tbb0
real(8), intent(in) :: h_ML0
real(8), intent(in) :: extwat
real(8), intent(in) :: extice
real(8), intent(in) :: kor_in
real(8), intent(in) :: trib_inflow
real(8), intent(in) :: Sals0
real(8), intent(in) :: Salb0
real(8), intent(in) :: fetch
real(8), intent(in) :: phi
real(8), intent(in) :: lam
real(8), intent(in) :: us0
real(8), intent(in) :: vs0
real(8), intent(in) :: Tm
real(8), intent(in) :: alphax
real(8), intent(in) :: alphay
real(8), intent(in) :: a_veg
real(8), intent(in) :: c_veg
real(8), intent(in) :: h_veg
real(8), intent(in) :: area_lake

!Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day
integer(4), intent(in) :: init_T

!Characters
character(len=60), intent(in) :: outpath1

!Logicals
logical, intent(in) :: flag_assim
logical, intent(in) :: flag_print

!Output variables
!Reals
real(8), intent(out) :: hw1
real(8), intent(out) :: xlew1
real(8), intent(out) :: cdmw1
real(8), intent(out) :: tsw

!Local variables
!Reals
real(8) :: l2
real(8) :: h2
real(8) :: dhwhigh
real(8) :: dhwlow
real(8) :: dhwp
real(8) :: dhwe
real(8) :: dhihigh
real(8) :: dhilow
real(8) :: dhip
real(8) :: totalevap
real(8) :: dhis
real(8) :: totalprecip
real(8) :: totalmelt
real(8) :: totalpen
real(8) :: totalwat
real(8) :: dhif
real(8) :: ls2
real(8) :: snowmass
real(8) :: snowmass_init
real(8) :: totalerad
real(8) :: totalhflux
real(8) :: totalpen1
real(8) :: Tf_old1
real(8) :: pi
real(8) :: year_sec
real(8) :: h_ML0zv
real(8) :: dhwls
real(8) :: x
real(8) :: kor
real(8) :: Ti10 ! Initial temperature of the ice surface (if l10 > 0)

real(8), dimension(1:vector_length) :: a
real(8), dimension(1:vector_length) :: b
real(8), dimension(1:vector_length) :: c
real(8), dimension(1:vector_length) :: d 

real(8) :: lams(ms) 
real(8) :: dz
real(8) :: ST
real(8) :: PGR
real(8) :: TGROLD
real(8) :: QGROLD
real(8) :: RADIAT
real(8) :: WSOLD
real(8) :: SNOLD
real(8) :: ZS
real(8) :: thsoil
real(8) :: whsoil
real(8) :: bOLD
real(8) :: RF1
real(8) :: RF2
real(8) :: SNMELT
real(8) :: HSNOW
real(8) :: HS
real(8) :: ES
real(8) :: TGRNEW
real(8) :: QGRNEW
real(8) :: WSNEW
real(8) :: SNNEW
real(8) :: RUNOFF
real(8) :: ElatOld
real(8) :: HFold
real(8) :: PRSold
real(8) :: extinct(ms)

real(8) :: AL
real(8) :: DLT
real(8) :: DVT
real(8) :: ALLL
real(8) :: DL
real(8) :: ALV
real(8) :: DV
real(8) :: Z
real(8) :: T
real(8) :: WL
real(8) :: WV
real(8) :: WI
real(8) :: dens

real(8) :: dt
real(8) :: q(110)
real(8) :: roughness
real(8) :: emissivity
real(8) :: albedo
real(8) :: aM
real(8) :: bM
real(8) :: relhums
real(8) :: b0
real(8) :: tau_air
real(8) :: tau_i
real(8) :: tau_gr

real(8), dimension(1:vector_length) :: Temp

real(8) :: flux1
real(8) :: flux2

real(8), external :: DZETA
real(8), external :: DZETAI


!Integers    
integer(4) :: i
integer(4) :: j
integer(4) :: i_ML
integer(4) :: flag_snow
integer(4) :: itop
integer(4) :: mon 
integer(4) :: nmonth_old
integer(4) :: nhour_old
integer(4) :: iter
integer(4) :: xday
integer(4) :: xhour
integer(4) :: nstep_meas
integer(4) :: n_1cm
integer(4) :: n_5cm
integer(4) :: layer_case

!Logicals
logical :: soil_ts_ext_iter
logical :: uniform_depth
logical :: flag

!Characters
character(len=60) :: outpath

common /watericesnowarr/ lams,q
common /SOILDAT/ dz(ms),itop
common /bL/ ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD, &
& ZS,thsoil,whsoil,bOLD,RF1,RF2,SNMELT,HSNOW, &
& HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF, &
& ElatOld,HFold,PRSold,extinct
common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML), &
& ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML), &
& dens(ms)
common /ts_ext_iter/ Tf_old1,soil_ts_ext_iter 
common /surface/ roughness,emissivity,albedo,aM,bM,relhums
common /out/ outpath


data uniform_depth /.true./
data nstep_meas /0/
 
SAVE
  
!BODY OF PROGRAM

rosoil = 1200.  
pi= 4.d0*atan(1.)

year_sec  = 60*60*24*365
    
tempair= tempair1-273.15
humair= humair1
pressure    = pressure1
uwind = uwind1
vwind = vwind1
zref  = zref1
precip= precip1
Sflux0= Sflux01*(roa0/row0)
outpath= outpath1

if (kor_in == -999.d0) then
  kor = 2.d0*omega*dsin(phi*pi/180.d0)
else
  kor = kor_in
endif

if(ifrad == 1) then
  longwave  = longwave1
  shortwave = shortwave1
elseif (ifrad == 0) then
  longwave  = 0.d0
  shortwave = 0.d0
endif

if (uwind==0) uwind=0.1
if (vwind==0) vwind=0.1
wind = dsqrt(uwind**2+vwind**2)

!Control of the input atmospheric forcing
if (tempair>60.or.tempair<-90) then
  print*, 'The air temperature ', tempair, 'is illegal: STOP'
  STOP
elseif (dabs(humair)>1.) then
  print*, 'The air humidity ', humair, 'is illegal: STOP'
  STOP
elseif (pressure > 110000 .or. pressure < 90000.) then
  print*, 'The air pressure ', pressure, 'is illegal: STOP'
  STOP
elseif (dabs(uwind)>200.) then
  print*, 'The x-component of wind ', uwind, 'is illegal: STOP'
  STOP
elseif (dabs(vwind)>200.) then
  print*, 'The y-component of wind ', vwind, 'is illegal: STOP'
  STOP
elseif (dabs(longwave)>1000.) then
  print*, 'The longwave radiation ',longwave,'is illegal: STOP'
  STOP
elseif (dabs(shortwave)>1400.) then
  print*, 'The shortwave radiation ', shortwave, 'is illegal: STOP'
  STOP
elseif (dabs(precip)>1.) then
  print*, 'The atmospheric precipitation ',precip, &
  & 'is illegal: STOP'
  STOP
endif

!Erad = shortwave*(1-albedo)
dt = dtl
dt_keps = dt/real(nstep_keps)

if (init(ix,iy)==0) then
! Specification of initial profiles   
  h1 = h10
  l1 = l10
  hs1 = hs10
  ls1 = ls10
  do i=1,M 
    E1(i) = 10.d0**(-5.5) !+ float(i)/float(M)*(1.d-10  -  1.d-8)
    eps1(i)= 1.d-9 !+ float(i)/float(M)*(1.d-18  -  1.d-14)
  enddo
  do i=1,ns
    if (i<0.7*ns) then 
      Tsoil1(i) = Tb0 + float(i-1)/(0.7*float(ns))*(Tbb0-Tb0)
    else
      Tsoil1(i) = Tbb0
    endif 
    Sals1(i) = Salb0 ! assuming, that salinity in ground is the same 
                     ! as one in near bottom layer of water
    if (Tsoil1(i)>Meltpnt(Sals1(i))) then
      wi1(i) = 0.
      wl1(i) = wl_max(por(i),rosoil(i),0.d0)-0.01
    else
      wl1(i) = unfrwat(Tsoil1(i),i)
      wi1(i) = wi_max(por(i),rosoil(i))-wl1(i)-0.01
    endif
  enddo

  i_ML = int(M*h_ML0/h1)
  Sal1(1:max(i_ML,1)) = Sals0
  u1(1:max(i_ML,1)) = us0
  v1(1:max(i_ML,1)) = vs0
  do i=max(i_ML+1,2),M+1
    Sal1(i) = Sals0 + (Salb0-Sals0)*float(i-i_ML)/float(M+1-i_ML)
    u1(i) = us0 + (0.d0 - us0)*float(i-i_ML)/float(M+1-i_ML)
    v1(i) = vs0 + (0.d0 - vs0)*float(i-i_ML)/float(M+1-i_ML)
  enddo
 
  if (init_T == 2) then
    h_ML0zv = (Tm-0.5*(Tb0+Ts0))/(Ts0/h1-0.5*(Ts0+Tb0)/h1)
    i_ML = int(M*h_ML0zv/h1)  
  endif 
  if (init_T == 1 .or. init_T == 2) then
    Tw1(1:max(i_ML,1)) = Ts0
    do i=max(i_ML+1,2), M+1 
      Tw1(i) = Ts0 + (Tb0-Ts0)*float(i-i_ML)/float(M+1-i_ML)
    enddo
! The initial temperature profile is given from the input file
  elseif (init_T == 3) then
    do i = 1, M+1
      z_full  (i) = DZETA(dble(i))*h1
    enddo
    call LININTERPOL (zTinitprof,Tinitprof,Tinitlength, &
    & z_full,Tw1,M+1,flag)
    if (.not.flag) then
      print*, 'The error while interpolating the initial &
      &temperature profile: terminating program'
      STOP
    endif
  endif
  
  if (skin == 0) then
    Tskin(1:2) = 0.d0
  else
    Tskin(1) = Tw1(1)
  endif  
    
 
  Ti1 = 0.
  Tis1 = 0. 
  Ti10 = dmin1(tempair,-1.d-1) ! Initial ice surface temperature, Celsius
  if (l1 /= 0) then
    do i = 1, Mice+1
      Ti1(i)= Ti10 - (Ti10 - Meltpnt(Sals0))* &
      & float(i-1)/float(Mice) ! Linear profile, if dzetai-grid is regular,
    enddo   ! water layer underneath is assumed to exist
  endif
  
  if (hs1==0.) then   
    flag_snow = 0
    flag_snow_init = 1
  else
    flag_snow = 1
    flag_snow_init = 0
    hs1 = dmax1(2.d-2, hs1)
    hs1 = int(hs1/0.01)*0.01
    n_5cm = int(hs1/0.05)
    n_1cm = int((hs1-n_5cm*0.05)/0.01)
    if (n_1cm==0) then
      n_1cm=1; hs1=hs1+0.01
    endif
    itop = ms - (n_1cm+n_5cm)
    do i=itop,itop+n_1cm-1
      dz(i) = 0.01
      T(i) = -5.
      wl(i) = 0
      dens(i) = 150.
    enddo
    do i=itop+n_1cm,ms-1
      dz(i) = 0.05
      T(i) = -5.
      wl(i) = 0
      dens(i) = 150.
    enddo
    dz(ms) = 0.5d0*dz(ms-1)
    T(ms) = -5.d0
    wl(ms) = 0.d0
    dens(ms) = 150.d0
  endif
  
  snmelt = 0.d0 
  cdmw2 = 1.d-15
  velfrict_prev = 1.d-2 
  eflux0_kinem = 0.d0
  Elatent = 0.

  totalevap = 0 
  totalmelt = 0 
  totalprecip = 0. 
  totalwat = 0
  totalpen = 0
  totalpen1 = 0
  time = 0
  nstep = 0
  dhwfsoil = 0.
  dhw = 0.
  dhw0 = 0.
  dhi = 0.
  dhi0 = 0.
  dls0 = 0.

  lamw(:) = 1.d+3

  par = 1
  mon = 1  
  iter = 0 
  nmonth_old = 1
  nhour_old = 21
  xday=1
  xhour=0


  if (runmode==2.and.(uniform_depth .eqv. .false.)) then
    if (ix<=nx-1.and.iy<=ny-1) then
      h1 = dep_2d(ix,iy)
    else
      h1 = dep_av
    endif  
    if (h1<=0) then
      print*, 'Negative or zero initial depth at ix=',ix, &
      & 'iy=',iy,'h1=',h1,':terminating program'
      STOP 
    endif
  endif

! if (assim == 1.or.error_cov == 1) then
!   call measured_data_reader
! endif

else
!ASSIGNING VALUES FROM PREVIOUS TIME STEP (time) IN CURRENT POINT (ix,iy)
  l1 = l1_2d(ix,iy)
  watice(ix,iy)=l1_2d(ix,iy) 
  h1 = h1_2d(ix,iy)
  ls1 = ls1_2d(ix,iy) 
  hs1 = hs1_2d(ix,iy)

  do i=1,M+1
    u1(i)=u_2d(i,ix,iy)
    v1(i)=v_2d(i,ix,iy)
  enddo

  do i=1,M+1
    E1(i)=E_2d(i,ix,iy)
    eps1(i)=eps_2d(i,ix,iy)
  enddo
    
  do i=1,ns 
    Tsoil1(i) = Tsoil1_2d(i,ix,iy)
    Sals1(i) = Sals1_2d(i,ix,iy)
    wi1(i)=wi1_2d(i,ix,iy)
    wl1(i)=wl1_2d(i,ix,iy)
  enddo
    
  do i=1,M+1 
    Tw1(i) = Tw1_2d(i,ix,iy)
    Sal1(i) = Sal1_2d(i,ix,iy)
  enddo
  
  Tskin(1) = Tskin_2d(ix,iy)

 
  do i = 1, Mice + 1 
    Ti1(i) = Ti1_2d(i,ix,iy)
    Tis1(i) = Tis1_2d(i,ix,iy)
  enddo

  flag_snow = fl_sn_2d(ix,iy)
  flag_snow_init = fl_sn_init_2d(ix,iy)

  itop = itop_2d(ix,iy)

  do i=max(1,itop),ms
    dz(i) = dz_2d(i,ix,iy)
    T(i) = T_2d(i,ix,iy)
    wl(i) = wl_2d(i,ix,iy)
    dens(i) = dens_2d(i,ix,iy)
  enddo

  snmelt = snmelt_2d(ix,iy)
  
  cdmw2 = cdm2(ix,iy)

  totalevap = 0 
  totalmelt = 0 
  totalprecip = 0. 
  totalwat = 0
  totalpen = 0
  totalpen1 = 0
  
  time = time_2d(ix,iy)
  nstep = nstep_2d(ix,iy)
  dhwfsoil = dhwfsoil_2d(ix,iy)
  Elatent = Elatent_2d(ix,iy)
  elatent_prev=Elatent_2d(ix,iy)
  hflux_prev=hflux_2d(ix,iy)
  Radbal_prev=Radbal_2d(ix,iy)
  dhw = dhw_2d(ix,iy)
  dhw0 = dhw0_2d(ix,iy)
  dhi = dhi_2d(ix,iy)
  dhi0 = dhi0_2d(ix,iy)
  dls0 = dls0_2d(ix,iy)
  velfrict_prev = velfrict_2d(ix,iy)
  eflux0_kinem = eflux0_kinem_2d(ix,iy)

  lamw(1:M) = lamw_2d(ix,iy,1:M)

endif

!if (flag_assim==.true.) then
!Correcting model fields at the current timestep 
!using least-square data assimilation   
! do i = 1, M+1
!  j = (i-1)*n_modvar
!  x_mod(1,j+1)= Tw1(i)
!  x_mod(1,j+2)= E1(i)
! enddo
! call assim_lst_sq 
!endif
   
time = time + dt
nstep = nstep + 1

layer_case = 1

if (h1 > 0  .and. l1 > 0) layer_case = 2 
if (h1 ==0  .and. l1 > 0) layer_case = 3
if (h1 ==0  .and. l1 ==0) layer_case = 4

if (layer_case == 1 .and. &
& WATER_FREEZE_MELT(Tw1(1), 0.5*ddz(1)*h1, Meltpnt(Sal1(1)), +1) .and. &
& h1 - min_ice_thick * roi/row0 > min_water_thick) then
  layer_case = 2
endif 

if (layer_case == 3 .and. &
& WATER_FREEZE_MELT(Ti1(Mice+1), 0.5*ddzi(Mice)*l1, Meltpnt(0.d0), -1) .and. &
& l1 - min_water_thick * row0/roi > min_ice_thick) then
  layer_case = 2
endif  


if (layer_case == 3 .and. &
& WATER_FREEZE_MELT(Ti1(1), 0.5*ddzi(1)*l1, Meltpnt(0.d0), -1) .and. &
& l1 - min_water_thick * row0/roi > min_ice_thick .and. flag_snow == 0) then
  h1 = min_water_thick
  Tw1 = Meltpnt(0.d0) + T_phase_threshold 
  Sal1 = 1.d-5
  ls1 = l1 - min_water_thick * row0/roi
  Tis1 = Ti1
  l2 = 0
  Ti2 = 0
  if (ls1<0.009999) print*, 'Thin layer! - ls'
  layer_case = 1 
endif

  
if1: IF (layer_case == 1) THEN 

!---------------------------- CASE 1: LIQUID WATER LAYER AT THE TOP ("Summer")-----------------------

!Check the bottom layer water temperature
if (WATER_FREEZE_MELT(Tw1(M+1), 0.5*ddz(M)*h1, Meltpnt(Sal1(M+1)), +1) .and. &
& ls1 == 0 .and. h1 - min_ice_thick * roi/row0 > min_water_thick) then
  ls1 = min_ice_thick
  Tis1 = Meltpnt(Sal1(M+1)) - T_phase_threshold
  h1 = h1 - min_ice_thick*roi/row0
endif   

!Check the liquid temperature to be above melting point
do i=2,M
  if (Tw1(i)<Meltpnt(Sal1(i))) Tw1(i)=Meltpnt(Sal1(i))
enddo

!Calculation of water density profile at the current timestep,
!taking into account only temperature
do i = 1, M+1
  row(i) = WATER_DENS_TS(Tw1(i),Sal1(i)) !0.d0
enddo

!The salinity and salinity flux at the surface are passed to
!TURB_DENS_FLUX as zeros since they are currently
!______
!not taken into account in calculation of density flux w'row'	
 
turb_density_flux = TURB_DENS_FLUX(eflux0_kinem,0.d0,Tw1(1),0.d0)

if (massflux == 1) then
!Updating density profile for massflux convection
!Updating z-grid properties 
  do i = 1, M+1
    z_full  (i) = DZETA(dble(i))*h1
  enddo
  do i = 1, M
    dz_full (i) = ddz(i)  *h1
    z_half  (i) = DZETA(dble(i)+0.5d0)*h1
  enddo
  call MASSFLUX_CONVECTION( &
  & nstep,dt,z_half,z_full,dz_full, &
  & turb_density_flux,eflux0_kinem, &
! The latent heat flux and salinity flux at the surface are passed as zeros,
! as far as currently they are not taken 
! into account in massflux conection parameterization
  & 0.d0,0.d0, &
! The surface wind stress components are passed as zeros,
! as far as currently they are not taken 
! into account in massflux conection parameterization
  & 0.d0,0.d0, &
  & u1,v1,w1,E1(1:M), &
  & u1,v1,w1,row,Tw1,Sal1, &
  & PEMF,PDENS_DOWN,PT_DOWN,PSAL_DOWN,zinv,1,M)
! Debugging purposes only
! PEMF = 0.d0
! call TEMP_MASSFLUX(PEMF,PT_DOWN,ddz,h1,dt,Tw2,T_massflux)
  do i = 1, M
    pt_down_f(i) = 0.5d0*(PT_DOWN(i) + PT_DOWN(i+1))
  enddo
endif

! The solar radiation, absorbed by the water surface
if (varalb == 0) then
  Erad = shortwave*(1-albedoofwater)
elseif (varalb == 1) then
  Erad = shortwave* &
  & (1-WATER_ALBEDO( SINH0(year,month,day,hour,phi) ) ) ! Note: the time here must be a local
endif  ! time, not UTC

!Radiation fluxes in layers
do i=1,M
  SR(i)=Erad*(1-sabs)*dexp(-extwat*DZETA(dble(i+0.5))*h1)
enddo
if (ls1 /= 0.d0) then
  do i = 1, Mice
    SRdi(i) = Erad*(1-sabs)*dexp(-extwat*h1)* &
    & (1-albedoofice)*dexp(-extice*DZETAI(dble(i+0.5))*ls1)
  enddo
endif

!Calculation of the layer's thickness increments
if (ls1==0) then
  ice=0;snow=0;water=1;deepice=0
  dhwp = precip*dt
  dhwe = - Elatent/(Lwv*row0)*dt
  dhw = dhwp + dhwe + dhwfsoil
  dhw0 = dhwp + dhwe
  dls = 0.
else
  ice=0;snow=0;water=1;deepice=1
  flux1 = -lamw(M)*(Tw1(M+1)-Tw1(M))/(ddz(M)*h1) + SR(M) + &
  & cw*row0*(dhw-dhw0)*(Tw1(M+1)-Tw1(M))/(2.*dt) + &
  & cw*row0*PEMF(M)*(pt_down_f(M)-0.5d0*(Tw1(M+1)+Tw1(M)) )
  flux2 = -lami*(Tis1(2)-Tis1(1))/(ddzi(1)*ls1) + SRdi(1) + &
  & ci*roi*dls0*(Tis1(2)-Tis1(1))/(2.*dt)
  dhwls = dt*(flux1 - flux2)/(row0*Lwi)
!  dhwls = -lamw(M)*dt/(ddz(M)*h1*row0*Lwi)*(Tw1(M+1) - Tw1(M))+
!&  lami*dt/(ddz(1)*ls1*row0*Lwi)*(Tis1(2) - Tis1(1))+
!&  (1-albedoofice)*Erad*(1-sabs)*dexp(-extwat*h1)/(row0*Lwi)*dt 
  dls = -dhwls*row0/roi
  dls0 = dls
  dhwp = precip*dt
  dhwe = - Elatent/(Lwv*row0)*dt
  dhw = dhwp + dhwe + dhwfsoil + dhwls
  dhw0 = dhwp + dhwe
endif

!Calculation of water current's velocities and turbulent characteristics
call MOMENT_SOLVER(ix, iy, nx, ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, &
& alphax, alphay, dt, b0, tau_air, tau_i, tau_gr)
do i = 1, nstep_keps
  call K_EPSILON(ix,iy,nx,ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, dt_keps, &
  & b0, tau_air, tau_i, tau_gr)
enddo

  
do i=1,M
  lamw(i)=lamw0+KT(i) !+KC(i) !+KLengT(i)
enddo
lamsal = lamw/(cw*row0) * alsal

!Calculation of the whole temperature profile
call SOIL_COND_HEAT_COEF()
CALL S_DIFF(dt)
CALL T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
& extwat, extice, fetch, dt)
call SOILFORLAKE &
& (ix,iy,nx,ny,year,month,day,hour,phi, &
& extwat,extice,fetch,dt,a,b,c,d,Temp)

!Diagnostic calculations
do i = 1, M
  k_turb_T_flux(i) = - lamw(i)/(cw*row0)*( Tw2(i+1) - Tw2(i) )/ &
  & (ddz(i)*h1)
  T_massflux(i) = PEMF(i)*(pt_down_f(i)-0.5d0*(Tw2(i+1)+Tw2(i)))
enddo

!Adding the heat input from tributaries
if (tribheat==1) call TRIBTEMP(dt,area_lake,Tw2)

!The scales of turbulence
if (scale_output == 1) then
  call TURB_SCALES(k_turb_T_flux, T_massflux, eflux0_kinem, h1, &
  & turb_density_flux, Buoyancy0, H_mixed_layer, w_conv_scale, &
  & T_conv_scale)
endif
  
h2 = h1 + dhw
ls2 = ls1 + dls
l2 = 0.

ENDIF if1   

if2: IF (layer_case==2) THEN

!---------------------------CASE 2: WATER AND ICE ("Winter")------------------------------


! Creation of the initial thin ice layer
if (l1 == 0.) then
  l1 = min_ice_thick
  h1 = h1 - min_ice_thick * roi/row0
  if (tempair<0.) then
    do i=1, Mice+1
      Ti1(i) = tempair + float(i-1)/float(Mice) * &
      & (Meltpnt(Sal1(1))-tempair)
    enddo 
  else 
    Ti1 = Meltpnt(Sal1(1)) - T_phase_threshold
  endif
end if 

!Check for water temperature to be higher than melting point
do i=2,M
  if (Tw1(i)<Meltpnt(Sal1(i))) Tw1(i)=Meltpnt(Sal1(i))
enddo

!Creation of the initial thin water layer
if (h1 == 0.) then
  h1 = min_water_thick
  Sal1 = 1.d-5
  Tw1 = Meltpnt(0.d0) + T_phase_threshold
  l1 = l1 - min_water_thick * row0/roi
end if

!Creation of the initial thin bottom ice layer
if (WATER_FREEZE_MELT(Tw1(M+1), 0.5*ddz(M)*h1, Meltpnt(Sal1(M+1)), +1) .and. &
& ls1 == 0 .and. h1 - min_ice_thick * roi/row0 > min_water_thick) then
  ls1 = min_ice_thick
  Tis1 = Meltpnt(Sal1(M+1)) - T_phase_threshold
  h1 = h1 - min_ice_thick * roi/row0
endif

!Calculation of water density profile at the current timestep,
!taking into account only temperature
do i = 1, M+1
  row(i) = WATER_DENS_TS(Tw1(i),Sal1(i)) !0.d0
enddo

!Currently the EDMF parameterization of convection is not implemented in the code
!during 'winter' (if the ice layer is present over
!water layer). This seems to be reasonable, since the temperature profile
!is stable in this case - until spring convection under thin ice
!develops.

PEMF = 0.d0
pt_down_f = 0.d0
  
!CASE 2.1: WATER, ICE AND SNOW 

if (flag_snow == 1) then

!Radiation fluxes in layers
  SR_botsnow = shortwave*(1-albedoofsnow)
  if (flag_snow_init /= 1) then
    do i = itop, ms-1
      SR_botsnow = SR_botsnow*EXTINCT_SNOW(dens(i))**dz(i)
    enddo
  endif
  do i = 1, Mice
    SRi(i) = SR_botsnow*dexp(-extice*DZETAI(dble(i+0.5))*l1)
  enddo
  do i = 1, M
    SR(i) = SR_botsnow*dexp(-extice*l1)* &
    & dexp(-extwat*DZETA(dble(i+0.5))*h1)
  enddo
  if (ls1 /= 0.d0) then
    do i = 1, Mice
      SRdi(i) = SR_botsnow*dexp(-extice*l1)*dexp(-extwat*h1)* &
      & (1-albedoofice)*dexp(-extice*DZETAI(dble(i+0.5))*ls1)
    enddo
  endif
  

!Calculation of the layer's thickness increments
  flux1=-lami*(Ti1(Mice+1)-Ti1(Mice))/(ddzi(Mice)*l1)+SRi(Mice) + &
  & ci*roi*(dhi-dhi0)*(Ti1(Mice+1)-Ti1(Mice))/(2.d0*dt)
  flux2 = -lamw(1)*(Tw1(2)-Tw1(1))/(ddz(1)*h1) + SR(1) + &
  & cw*row0*dhw0*(Tw1(2)-Tw1(1))/(2.*dt) + &
  & cw*row0*PEMF(1)*(pt_down_f(1)-0.5d0*(Tw1(2)+Tw1(1)) )
  dhwlow = dt*(flux1 - flux2)/(row0*Lwi)
! dhwlow = lamw(1)*dt/(ddz(1)*h1*row0*Lwi)*(Tw1(2) - Tw1(1))   VS,06.2007
!&  -lami*dt/(ddz(M)*l1*row0*Lwi)*(Ti1(M+1) - Ti1(M))VS,06.2007 
  dhilow = -dhwlow*row0/roi
  if (Ti1(1) > Meltpnt(0.d0) + T_phase_threshold) then
    dhwhigh = (Ti1(1) - Meltpnt(0.d0) - T_phase_threshold) * &
    & ci*roi*l1*ddzi(1)*0.5d0/(Lwi*row0)
    dhihigh = -dhwhigh*row0/roi
  else
    dhwhigh = 0.d0
    dhihigh = 0.d0
  endif
  if (ls1/=0) then
    ice=1;snow=1;water=1;deepice=1
    flux1 = -lamw(M)*(Tw1(M+1)-Tw1(M))/(ddz(M)*h1) + SR(M) + &
    & cw*row0*(dhw-dhw0)*(Tw1(M+1)-Tw1(M))/(2.d0*dt) + &
    & cw*row0*PEMF(M)*(pt_down_f(M)-0.5d0*(Tw1(M+1)+Tw1(M)) )
    flux2 = -lami*(Tis1(2)-Tis1(1))/(ddzi(1)*ls1) + SRdi(1) + &
    & ci*roi*dls0*(Tis1(2)-Tis1(1))/(2.d0*dt)
    dhwls = dt*(flux1 - flux2)/(row0*Lwi)
!  dhwls=-lamw(M)*dt/(ddz(M)*h1*row0*Lwi)*(Tw1(M+1) - Tw1(M))  VS,06.2007
!&   +lami*dt/(ddz(1)*ls1*row0*Lwi)*(Tis1(2) - Tis1(1))    VS,06.2007
    dls = -dhwls*row0/roi
    dhw = dhwhigh + dhwlow + snmelt*dt + dhwfsoil + dhwls
    dhw0 = dhwlow + dhwhigh + snmelt*dt
    dls0 = dls
  else
    ice=1;snow=1;water=1;deepice=0
    dls = 0.
    dhwls = 0.
    dhw = dhwhigh + dhwlow + snmelt*dt + dhwfsoil
    dhw0 = dhwlow + dhwhigh + snmelt*dt
  endif
  dhi = dhilow + dhihigh 
  dhi0 = dhihigh

!Calculation of water current's velocities and turbulent characteristics
  if (Turbpar/=1) then
    call MOMENT_SOLVER(ix, iy, nx, ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, &
    & alphax, alphay, dt, b0, tau_air, tau_i, tau_gr)  
    do i = 1, nstep_keps
      call K_EPSILON(ix,iy,nx,ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, dt_keps, &
      & b0, tau_air, tau_i, tau_gr)
    enddo
    do i=1,M
      lamw(i)=lamw0+KT(i) !+KC(i)
    enddo  
  else
    lamw=lamw0*3.d0
  endif
  lamsal = lamw/(cw*row0) * alsal
  

!Calculation of the whole temperature profile 
  call SNOW_COND_HEAT_COEF()
  call SOIL_COND_HEAT_COEF()
  call S_DIFF(dt)
  call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat, extice, fetch, dt)
  call SOILFORLAKE( &
  & ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat,extice,fetch,dt,a,b,c,d,Temp)
  call SNOWTEMP(ix,iy,nx,ny,year,month,day,hour,snowmass, &
  & snowmass_init,a,b,c,d,Temp,phi,extwat,extice,fetch,dt)

  h2 = h1 + dhw
  l2 = l1 + dhi 
  ls2 = ls1 + dls
  
else

!CASE 2.2: WATER AND ICE WITHOUT SNOW

!Radiation fluxes in layers
  Erad = shortwave*(1-albedoofice)
  do i = 1, Mice
    SRi(i) = Erad*dexp(-extice*DZETAI(dble(i+0.5))*l1)
  enddo
  do i = 1, M
    SR(i) = Erad*dexp(-extice*l1)* &
    & dexp(-extwat*DZETA(dble(i+0.5))*h1)
  enddo
  if (ls1 /= 0.d0) then
    do i = 1, Mice
      SRdi(i) = Erad*dexp(-extice*l1)*dexp(-extwat*h1)* &
      & (1-albedoofice)*dexp(-extice*DZETAI(dble(i+0.5))*ls1)
   enddo
  endif

!Calculation of the layer's thickness increments
  if (Ti1(1) > Meltpnt(0.d0) + T_phase_threshold) then
    dhwhigh = (Ti1(1) - Meltpnt(0.d0) - T_phase_threshold) * &
    & ci*roi*ddzi(1)*l1*0.5d0/(Lwi*row0)
    dhis = 0. 
    dhihigh = -dhwhigh*row0/roi
  endif
! Creation of the initial thin snow layer - to be considered at the next timestep
  if (precip>0.and.tempair<0.and.l1>0.05) then
    flag_snow = 1
    hs1 = 0.02
    itop = ms-2
    dhip = 0
  else
    dhip = precip*row0/roi*dt
  end if
! write(*,*) 'I am here',Mice,ubound(Ti1),ubound(ddzi),
!&  ubound(SRi)
! STOP
  flux1=-lami*(Ti1(Mice+1)-Ti1(Mice))/(ddzi(Mice)*l1)+SRi(Mice) + &
  &  ci*roi*(dhi-dhi0)*(Ti1(Mice+1)-Ti1(Mice))/(2.d0*dt)
  flux2 = -lamw(1)*(Tw1(2)-Tw1(1))/(ddz(1)*h1) + SR(1) + &
  &  cw*row0*dhw0*(Tw1(2)-Tw1(1))/(2.*dt) + &
  &  cw*row0*PEMF(1)*(pt_down_f(1)-0.5d0*(Tw1(2)+Tw1(1)) )
  dhwlow = dt*(flux1 - flux2)/(row0*Lwi)
! dhwlow = lamw(1)*dt/(ddz(1)*h1*row0*Lwi)*(Tw1(2) - Tw1(1))    VS,06.2007
!&  -lami*dt/(ddz(M)*l1*row0*Lwi)*(Ti1(M+1) - Ti1(M))VS,06.2007
  dhilow = -dhwlow*row0/roi

  if (Ti1(1) < Meltpnt(0.d0) + T_phase_threshold) then
    dhis = - Elatent/(roi*Liv)*dt
    dhwhigh = 0
    dhihigh = 0
  endif
  dhi = dhihigh + dhilow + dhis + dhip
  dhi0 = dhis + dhihigh + dhip
  if (ls1/=0) then
    ice=1;snow=0;water=1;deepice=1 
    flux1 = -lamw(M)*(Tw1(M+1)-Tw1(M))/(ddz(M)*h1) + SR(M) + &
    & cw*row0*(dhw-dhw0)*(Tw1(M+1)-Tw1(M))/(2.d0*dt) + &
    & cw*row0*PEMF(M)*(pt_down_f(M)-0.5d0*(Tw1(M+1)+Tw1(M)) )
    flux2 = -lami*(Tis1(2)-Tis1(1))/(ddzi(1)*ls1) + SRdi(1) + &
    & ci*roi*dls0*(Tis1(2)-Tis1(1))/(2.d0*dt)
    dhwls = dt*(flux1 - flux2)/(row0*Lwi)
!  dhwls =-lamw(M)*dt/(ddz(M)*h1*row0*Lwi)*(Tw1(M+1) - Tw1(M))  VS,06.2007
!&   +lami*dt/(ddz(1)*ls1*row0*Lwi)*(Tis1(2) - Tis1(1))    VS,06.2007
    dls = -dhwls*row0/roi
    dls0 = dls
    dhw = dhwhigh + dhwfsoil + dhwls + dhwlow
    dhw0 = dhwhigh+dhwlow
  else
    ice=1;snow=0;water=1;deepice=0
    dhw = dhwhigh + dhwfsoil + dhwlow
    dhw0 = dhwhigh + dhwlow
    dls = 0.
  endif

!Calculation of water current's velocities and turbulent characteristics
  if (Turbpar/=1) then
    call MOMENT_SOLVER(ix, iy, nx, ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, &
    & alphax, alphay, dt, b0, tau_air, tau_i, tau_gr)  
    do i = 1, nstep_keps
      call K_EPSILON(ix,iy,nx,ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, dt_keps, &
      & b0, tau_air, tau_i, tau_gr)
    enddo
    do i=1,M
      lamw(i)=lamw0+KT(i) !+KC(i)
    enddo  
  else
    lamw=lamw0*3.d0
  endif
  lamsal = lamw/(cw*row0) * alsal

!Calculation of the whole temperature profile
  call SOIL_COND_HEAT_COEF()
  call S_DIFF(dt)
  call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat, extice, fetch, dt)
  call SOILFORLAKE( &
  & ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat,extice,fetch,dt,a,b,c,d,Temp)

  h2 = h1 + dhw
  l2 = l1 + dhi
  ls2 = ls1 + dls
 
endif

if (tribheat==1) call TRIBTEMP(dt,area_lake,Tw2)
 
ENDIF if2

!CASE 3: ICE WITHOUT WATER

if3: IF (layer_case==3) THEN
    
!CASE 3.1: ICE WITH SNOW

if (flag_snow == 1) then

!Radiation fluxes
  SR_botsnow = shortwave*(1-albedoofsnow)
  if (flag_snow_init /= 1) then
    do i = itop, ms-1
      SR_botsnow = SR_botsnow*EXTINCT_SNOW(dens(i))**dz(i)
    enddo
  endif
  do i = 1, Mice
    SRi(i) = SR_botsnow*dexp(-extice*DZETAI(dble(i+0.5))*l1)
  enddo

  if (Ti1(1) > Meltpnt(0.d0) + T_phase_threshold) then
    dhwhigh = (Ti1(1) - Meltpnt(0.d0) - T_phase_threshold) * &
    & ci*roi*l1*ddzi(1)*0.5d0/(Lwi*row0)
    dhif = -dhwhigh*row0/roi
  else
    dhwhigh = 0
    dhif = 0
  endif
  dhw = snmelt*dt + dhwhigh
  dhi = dhif
  dhi0 = dhi
  ice=1;snow=1;water=0;deepice=0

  call SNOW_COND_HEAT_COEF()
  call SOIL_COND_HEAT_COEF()
  call S_DIFF(dt)
  call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat, extice, fetch, dt)
  call SOILFORLAKE( &
  & ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat,extice,fetch,dt,a,b,c,d,Temp)
  call SNOWTEMP(ix,iy,nx,ny,year,month,day,hour,snowmass, &
  & snowmass_init,a,b,c,d,Temp,phi,extwat,extice,fetch,dt)

  l2 = l1 + dhi 
  h2 = h1 + dhw
  ls2 = 0.

else

!CASE 3.2: ICE WITHOUT SNOW !

!Radiation fluxes
  Erad = shortwave*(1-albedoofice)
  do i = 1, Mice
    SRi(i) = Erad*dexp(-extice*DZETAI(dble(i+0.5))*l1)
  enddo

  if (precip>0.and. tempair<0) then
    flag_snow = 1
    hs1 = 0.02
    dhip = 0
  else
    dhip = precip*row0/roi*dt
  end if
  if (Ti1(1) > Meltpnt(0.d0) + T_phase_threshold) then
    dhw = (Ti1(1) - Meltpnt(0.d0) - T_phase_threshold) * &
    & ci*roi*l1*ddzi(1)*0.5d0/(Lwi*row0)
    dhihigh = - dhw*row0/roi
    dhis = 0
  else
    dhis = - Elatent/(roi*Liv)*dt
    dhw = 0.
    dhihigh = 0
  endif
  dhi = dhis + dhihigh + dhip
  dhi0 = dhi
  ice=1;snow=0;water=0;deepice=0 

  call SOIL_COND_HEAT_COEF()
  call S_DIFF(dt)
  call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat, extice, fetch, dt)
  call SOILFORLAKE( &
  & ix,iy,nx,ny,year,month,day,hour,phi, &
  & extwat,extice,fetch,dt,a,b,c,d,Temp)
    
  l2 = l1 + dhi
  h2 = h1 + dhw
  ls2 = 0.
 
endif

ENDIF if3 

! CASE 4:SNOW ANS SOIL WITHOUT ICE AND WATER
    
if4: IF (layer_case==4) THEN
  print*, 'The non-operational case: &
  &there is no water and ice layers: STOP'
  STOP
  if (flag_snow==1) then
    ice=0;snow=1;water=0;deepice=0
    call SNOW_COND_HEAT_COEF()
    call SOIL_COND_HEAT_COEF()
    call S_DIFF(dt)
    call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
    & extwat, extice, fetch, dt)
    call SOILFORLAKE( &
    & ix,iy,nx,ny,year,month,day,hour,phi, &
    & extwat,extice,fetch,dt,a,b,c,d,Temp)
    call SNOWTEMP(ix,iy,nx,ny,year,month,day,hour,snowmass, &
    & snowmass_init,a,b,c,d,Temp,phi,extwat,extice,fetch,dt)
  else
    ice=0;snow=0;water=0;deepice=0
    call SOIL_COND_HEAT_COEF()
    call S_DIFF(dt)
    call T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
    & extwat, extice, fetch, dt)
    call SOILFORLAKE( &
    & ix,iy,nx,ny,year,month,day,hour,phi, &
    & extwat,extice,fetch,dt,a,b,c,d,Temp)
  endif 
ENDIF if4


!TRIBUTARY INFLOW to a lake

h2 = h2 + trib_inflow/year_sec*dt 

if (h2 < 0.005.and.ls2/=0.and.l2==0) then
  ls2 = ls2 + h2*row0/roi
  if (ls2 < 0.005) then 
    ls2 = 0; Tis2 = 0
  endif
  h2 = 0.
  Tw2 = 0.
  Sal2 = 0.
endif

if ((h2<0.005.and.h2>0).and.l2 /= 0) then
  l2 = l2 + h2*row0/roi
  if (l2 < 0.005) then 
    l2 = 0; Ti2 = 0
  endif
  h2 = 0.
  Tw2 = 0.
  Sal2 = 0. 
end if
    
if (h2 < 0.005.and.l2/=0.and.ls2/=0) then
  Ti2 = Ti2*l2/(l2+ls2) + Ti2*ls2/(l2+ls2)
  l2 = l2 + h2*row0/roi + ls2 
  ls2 = 0.
  Tis2 = 0.
  h2 = 0 
  Tw2 = 0.
  Sal2 = 0.
endif
   
if (l2 < 0.005.and.h2 /= 0) then
  h2 = h2 + l2*roi/row0
  if (h2 < 0.005) then 
    h2 = 0; Tw2 = 0; Sal2 = 0.
  endif
  l2 = 0.
  Ti2 = 0.
end if

if (h2 < 0.005.and. l2 == 0) then
  h2 = 0.
  Tw2 = 0.
  Sal2 = 0. 
end if

if (l2 < 0.005.and. h2 == 0) then
  l2 = 0.
  Ti2 = 0.
end if
   
if ((hs1 < 0.005.and.hs1>0).or.(hs1>0.and.l2==0.and.h2/=0)) then
! h2 = h2 + (totalprecips-totalevaps-totalmelts)
  h2 = h2 + snowmass/row0
  if (h2 < 0.005) then
    l2 = l2 + row0*h2/roi
    h2 = 0; Tw2 = 0; Sal2 = 0.
    if (l2 < 0.005) then 
      l2 = 0; Ti2 = 0
    endif
  endif
  hs1 = 0.
  flag_snow = 0. 
  flag_snow_init = 1
end if
 
if (ls2 < 0.005.and.h2/=0) then
  h2 = h2 + ls2*roi/row0
  if (h2 < 0.005) then 
    h2 = 0; Tw2 = 0; Sal2 = 0.
  endif
  ls2 = 0
  Tis2 = 0.
endif 

if (ls2>0.and.h2 == 0) then
  l2 = ls2
  Ti2 = Tis2
  if (l2 < 0.005) then 
    l2 = 0; Ti2 = 0
  endif
  ls2 = 0
  Tis2 = 0.
endif
    
h1    = h2
l1    = l2
ls1   = ls2
Tw1   = Tw2
Ti1   = Ti2
Tis1  = Tis2
Sal1  = Sal2
Sals1 = Sals2
u1 = u2
v1 = v2
Tskin(1) = Tskin(2)
    
! CALCULATION OF SUMMARY FLUXES !
totalpen = totalpen + dhwfsoil 
totalpen1 = totalpen1 + dhwfsoil
totalevap = totalevap + Elatent/(row0*Lwv)*dt
totalprecip = totalprecip + precip*dt
totalhflux = totalhflux + hflux*dt
totalerad = totalerad + erad*dt

if (l1==0) tsw = Tw1(1) + 273.15
if (l1/=0.and.flag_snow==0) tsw = Ti1(1) + 273.15
if (flag_snow==1) tsw = T(itop) + 273.15

if (init(ix,iy)==0) then
  init(ix,iy)=1  
endif

!VALUES AT NEXT TIME STEP (time + dt) IN CURRENT POINT (ix,iy)

l1_2d(ix,iy)  = l1
watice(ix,iy)=l1_2d(ix,iy) 
h1_2d(ix,iy)  = h1
ls1_2d(ix,iy)  = ls1
hs1_2d(ix,iy) = hs1

do i=1,M+1
  u_2d(i,ix,iy)=u1(i)
  v_2d(i,ix,iy)=v1(i)
enddo

do i=1,M+1
  E_2d(i,ix,iy)=E1(i)
  eps_2d(i,ix,iy)=eps1(i)
enddo
    
do i=1,ns 
  Tsoil1_2d(i,ix,iy)=Tsoil1(i)
  Sals1_2d(i,ix,iy) = Sals1(i) 
  wi1_2d(i,ix,iy)=wi1(i)
  wl1_2d(i,ix,iy)=wl1(i)
enddo
    
do i=1,M+1 
  Tw1_2d(i,ix,iy) = Tw1(i)
  Sal1_2d(i,ix,iy) = Sal1(i) 
enddo

Tskin_2d(ix,iy) = Tskin(1)

do i=1,Mice+1
  Ti1_2d(i,ix,iy) = Ti1(i)
  Tis1_2d(i,ix,iy) = Tis1(i)
enddo
 
fl_sn_2d(ix,iy) = flag_snow
fl_sn_init_2d(ix,iy) = flag_snow_init

itop_2d(ix,iy) = itop

do i=max(1,itop),ms
  dz_2d(i,ix,iy) = dz(i)
  T_2d(i,ix,iy) = T(i)
  wl_2d(i,ix,iy) = wl(i)
  dens_2d(i,ix,iy) = dens(i)
enddo

snmelt_2d(ix,iy) = snmelt

cdm2(ix,iy) = cdmw
time_2d(ix,iy) = time
nstep_2d(ix,iy) = nstep
dhwfsoil_2d(ix,iy) = dhwfsoil
Elatent_2d(ix,iy) = Elatent
Radbal_2d(nx,ny)=Radbal
hflux_2d(nx,ny)=hflux
dhw_2d(ix,iy) = dhw
dhw0_2d(ix,iy) = dhw0
dhi_2d(ix,iy) = dhi
dhi0_2d(ix,iy) = dhi0
dls0_2d(ix,iy) = dls0
velfrict_2d(ix,iy) = velfrict
eflux0_kinem_2d(ix,iy) = eflux0_kinem

lamw_2d(ix,iy,1:M) = lamw(1:M)

hw1 = hw
xlew1 = xlew
cdmw1 = cdmw

!if (error_cov==1) then
! if (int(nstep*dt/(interval*60*60))>=nstep_meas) then
!  nstep_meas = nstep_meas + 1
!  if (nstep_meas<=ntim) then
!   do i = 1, M+1
!    j = (i-1)*n_modvar
!    x_mod(nstep_meas,j+1)= Tw2(i)
!    x_mod(nstep_meas,j+2)= E2(i)
!   enddo 
!  endif
! endif
!endif 
!if (error_cov==1.and.nstep_meas==ntim) then
! call HL_method
! error_cov = 0
!endif

!Diagnoctics

!call fluxes

if (ix==10.and.iy==10) print*, 'Lake', &
& 'T1=',Tw1(1),'T2=',Tw1(2),'h1=',h1,'S=',shortwave, &
& 'hh=',hflux,'LE=',Elatent, &
& 'eflux=',eflux, 'Longwave=', longwave, &
& 'Bal=',(shortwave*(1-albedoofwater)*(sabs)+longwave- &
&  surfrad-hflux-Elatent-eflux)/(cw*row0*ddz(1)*h1/2)*dt

if (flag_print) then ! The output in ASCII files 
  if (monthly_out ==1) call MON_OUT(ix,iy,nx,ny,year,month,day,hour)
  if (daily_out   ==1) call DAY_OUT(ix,iy,nx,ny,year,month,day,hour)
  if (hourly_out  ==1) call HOUR_OUT(ix,iy,nx,ny,year,month,day,hour)
  if (everystep   ==1) call EVERYSTEP_OUT(ix,iy,nx,ny)
  if (time_series ==1) call SERIES_OUT(ix,iy,nx,ny,year,month,day,hour,tsw)
endif
 
if (runmode == 1)  then
  if (mod(nstep,nscreen)==0) then
    write (*,'(3a25)') 'timestep', 'N of point', &
    & 'Surface temperature'
    write (*,'(2i25,f25.1)') nstep, ix, tsw
    write (*,'(4a8)') 'Year', 'Month', 'Day', 'Hour'
    i = year 
    if (year<1000) i = year+1000
    write (*,'(3i8,f8.2)') i, month, &
    & day   , hour
  endif
endif
  
!FORMATS!
7   format (f7.3, 4i5,37f7.3) 
60  format (f5.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
! 61  format (i5, f6.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
61  format (f6.2, f5.2, 2e16.5, 2f8.3, e15.5, f11.5, f7.1)
62  format (f6.2, f5.2, 2f12.3, 2f8.3, f11.5, f7.1,3f12.4,2f7.2,3f7.1)
80  format (f7.3, 10f9.2)
90  format (f7.3, 12f9.2)   
100 format (a5, 2a16, 2a8, a15, a11, a7) 
!FINISHING PROGRAM! 
  
END SUBROUTINE LAKE

!-----------------------------------------------------------------------------------
!TEMPERATURE EQUATION SOLVER
!-----------------------------------------------------------------------------------

SUBROUTINE T_SOLVER(ix,iy,nx,ny,year,month,day,hour,phi, &
& extwat, extice, fetch, dt)

!T_SOLVER implements iterations to find water surface temperature 
!at the next time step    

use DRIVING_PARAMS
use ARRAYS
implicit none

!Input variables
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny

integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day

real(8), intent(in) :: hour
real(8), intent(in) :: phi
real(8), intent(in) :: extwat
real(8), intent(in) :: extice
real(8), intent(in) :: fetch
real(8), intent(in) :: dt

!Local variables
real(8), allocatable :: Tsurf_prev(:,:)
integer(4), allocatable :: surftyp_prev(:,:)
integer(4), allocatable :: surftyp(:,:)

real(8) :: dt_scan
real(8) :: Tsurf1
real(8) :: Tsurf2
real(8) :: Tsurf3
real(8) :: Bal1
real(8) :: Bal2
real(8) :: Bal3

integer(4) :: iter
integer(4) :: maxiter
integer(4) :: n_unit

logical firstcall
data firstcall /.true./
data n_unit /666/

!External functions
real(8), external :: HEATBALANCE

SAVE

if (firstcall) then
  call CHECK_UNIT(n_unit)
  open (n_unit,file=path(1:len_trim(path))// &
  & 'results/debug/iter_T.dat')
  maxiter = 10
  dt_scan = 5.

  allocate (Tsurf_prev(nx,ny))
  allocate (surftyp_prev(nx,ny))
  allocate (surftyp(nx,ny))
endif

if (snow==1) then
  surftyp(ix,iy) = 3
elseif (ice==1) then
  surftyp(ix,iy) = 2
elseif (water==1) then
  surftyp(ix,iy) = 4
else
  surftyp(ix,iy) = 1
endif
call SURF_CHAR(surftyp(ix,iy),year,month,day,hour,phi)
 
if (nstep > 1.and.surftyp(ix,iy) == surftyp_prev(ix,iy)) then
  Tsurf2 = Tsurf_prev(ix,iy) - 10.
elseif (surftyp(ix,iy) == 4) then
  Tsurf2 = -10.
else
  Tsurf2 = -90.
endif

!print*, 'In T_SOLVER: ', Tsurf2, surftyp(ix,iy)

!SCANNING INTERVAL OF SURFACE TEMPERATURE [-90 C, ...]

Bal1=1.;Bal2=1.
do while (Bal1*Bal2>0)
  Tsurf1=Tsurf2
  Tsurf2=Tsurf2+dt_scan
  call T_diff(1,Tsurf1,dt)
  Bal1 = HEATBALANCE(Tsurf1,surftyp(ix,iy),extwat,extice,fetch,dt) 
  call T_diff(1,Tsurf2,dt)
  Bal2 = HEATBALANCE(Tsurf2,surftyp(ix,iy),extwat,extice,fetch,dt) 
  if (Tsurf2>80.) then
    print*, 'Severe: iterations for the surface temperature &
    &do not converge: at the point', ix, iy, 'STOP', 'nstep = ', nstep, lamw
!    print*, 'Temperature limit in scan process is exceeded: STOP'
    STOP
  endif
enddo

iter=0
Bal3=1.

!CHORDE METHOD TO FIND SURFACE TEMPERATURE
if (Bal1 /= 0.d0.and.Bal2 /= 0.d0) then
  cy1:do while (dabs(Bal3)>0.1)
    iter=iter+1
    call T_diff(1,Tsurf1,dt)
    Bal1 = HEATBALANCE(Tsurf1,surftyp(ix,iy),extwat,extice,fetch,dt) 
    call T_diff(1,Tsurf2,dt)
    Bal2 = HEATBALANCE(Tsurf2,surftyp(ix,iy),extwat,extice,fetch,dt) 
    Tsurf3=(Tsurf1*Bal2-Tsurf2*Bal1)/(Bal2-Bal1)
    call T_diff(1,Tsurf3,dt)
    Bal3 = HEATBALANCE(Tsurf3,surftyp(ix,iy),extwat,extice,fetch,dt)
    if (Bal1*Bal3<0.) then
      Tsurf2=Tsurf3
    elseif (Bal2*Bal3<0.) then
      Tsurf1=Tsurf3
    endif
    if (iter>maxiter) then
      write (666,*) Bal3
      exit cy1
    endif
  enddo cy1
elseif (Bal1 == 0.d0) then
  Tsurf3 = Tsurf1
elseif (Bal2 == 0.d0) then
  Tsurf3 = Tsurf2
endif

call T_diff(0,Tsurf3,dt)

Tsurf_prev(ix,iy) = Tsurf3
surftyp_prev(ix,iy) = surftyp(ix,iy)

if (firstcall) firstcall = .false.
END SUBROUTINE T_SOLVER


SUBROUTINE T_DIFF(surf,Tsurf,dt)

!T_DIFF calculates temperature profile in water, soil, ice and snow,
!with known temparature at the surface 

use NUMERIC_PARAMS
use PHYS_CONSTANTS2
use DRIVING_PARAMS
use ARRAYS
use PHYS_FUNC, only: &
& MELTPNT 
implicit none

!Parameters
integer(4), parameter :: soil_indic = 1
integer(4), parameter :: ice_indic = 2
integer(4), parameter :: snow_indic = 3
integer(4), parameter :: water_indic = 4
integer(4), parameter :: deepice_indic = 5

integer(4), parameter :: vector_length = 350

!Variables
real(8) dt,Tsurf,Hflow
real(8), dimension(1:ms) :: Tsn,lams,q,cs
real(8), dimension(1:vector_length) :: a,b,c,d,Temp
real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
real(8) dz
integer(4) i,j,itop,surf
   

common /snow_char/ Tsn,cs
common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML), &
& ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
common /SOILDAT/ dz(ms),itop
common /watericesnowarr/ lams,q
!data num /1/

SAVE

Hflow=0.
!Meltpnt - Melting point temperature, C degrees 
!Meltpnt=0.
!ddz=1./float(M)VS,06.2007

if (surf==1) then
  if (snow==1) then
    c(itop)=1.
    b(itop)=0.
    d(itop)=Tsurf
    if (itop<=ms-2) then
      call DIFF_COEF(a,b,c,d,itop+1,ms-1,itop+1,ms-1,snow_indic,dt) 
    endif
    if (ice==1) then
!----------------SNOW-ICE INTERFACE
      a(ms)=-(lams(ms-1)+lams(ms))/(2.*dz(ms-1))
      b(ms)=-lami/(ddzi(1)*l1)+ci*roi*dhi0/(2.*dt)
      c(ms)=a(ms)+b(ms)-(cs(ms)*dens(ms)*dz(ms)/(2.*dt) + &
      & ci*roi*ddzi(1)*l1/(2.*dt))
      d(ms)=-Ti1(1)*(cs(ms)*dens(ms)*dz(ms)/(2.*dt) + &
      & ci*roi*ddzi(1)*l1/(2.*dt))-SR_botsnow + SRi(1)
!-----------------------------------
      call DIFF_COEF(a,b,c,d,2,Mice,ms+1,ms+Mice-1,ice_indic,dt) 
      if (water==1) then
        a(ms+Mice)=0.
        c(ms+Mice)=1.
        d(ms+Mice)=Meltpnt(Sal2(1))
        call PROGONKA (a,b,c,d,Temp,itop,ms+Mice)
        do i=itop, ms
          Tsn(i)=Temp(i)
        enddo
        do i=ms, ms+Mice
          Ti2(i-ms+1)=Temp(i)
        enddo 
      else
!----------------ICE-SOIL INTERFACE-------------------------
        a(ms+Mice)=-lami/(ddzi(Mice)*l1)+ci*roi*(dhi-dhi0)/(2.*dt)
        b(ms+Mice)=-(lamsoil(1)+lamsoil(2))/(2.*dzs(1))
        c(ms+Mice)=a(ms+Mice)+b(ms+Mice) - &
        & (csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & ci*roi*ddzi(Mice)*l1/(2.*dt))
        d(ms+Mice)=-Ti1(Mice+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & ci*roi*ddzi(Mice)*l1/(2.*dt))-SRi(Mice)
!----------------------------------------------------------- 
        call DIFF_COEF(a,b,c,d,2,ns-1,ms+Mice+1,ms+Mice+ns-2,soil_indic,dt)
        c(ms+Mice+ns-1)=1.
        a(ms+Mice+ns-1)=1.
        d(ms+Mice+ns-1)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
        call PROGONKA (a,b,c,d,Temp,itop,ms+Mice+ns-1)
        do i=itop, ms
          Tsn(i)=Temp(i)
        enddo
        do i=ms, ms+Mice
          Ti2(i-ms+1)=Temp(i)
        enddo 
        do i=ms+Mice, ms+Mice+ns-1
          Tsoil2(i-ms-Mice+1)=Temp(i)
        enddo
      endif
    else
!-------------------SNOW-SOIL INTERFACE-------------------------
      a(ms)=-(lams(ms-1)+lams(ms))/(2.*dz(ms-1))
      b(ms)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
      c(ms)=a(ms)+b(ms)-(cs(ms)*dens(ms)*dz(ms)/(2.*dt) + &
      & csoil(1)*rosoil(1)*dzs(1)/(2.*dt))
      d(ms)=-Tsoil1(1)*(cs(ms)*dens(ms)*dz(ms)/(2.*dt) + &
      & csoil(1)*rosoil(1)*dzs(1)/(2.*dt))-SR_botsnow
!---------------------------------------------------------------
      call DIFF_COEF(a,b,c,d,2,ns-1,ms+1,ms+ns-2,soil_indic,dt) 
      c(ms+ns-1)=1.
      a(ms+ns-1)=1.
      d(ms+ns-1)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
      call PROGONKA (a,b,c,d,Temp,itop,ms+ns-1)
      do i=itop, ms
        Tsn(i)=Temp(i)
      enddo
      do i=ms, ms+ns-1
        Tsoil2(i-ms+1)=Temp(i)
      enddo
    endif
  elseif (ice==1) then
    b(1)=0.
    c(1)=1.
    d(1)=Tsurf
    call DIFF_COEF(a,b,c,d,2,Mice,2,Mice,ice_indic,dt)  
    if (water==1) then
      a(Mice+1)=0.
      c(Mice+1)=1.
      d(Mice+1)=Meltpnt(Sal2(1))
      call PROGONKA (a,b,c,d,Temp,1,Mice+1)
      do i=1,Mice+1
        Ti2(i)=Temp(i)
      enddo
    else
!----------------ICE-SOIL INTERFACE-------------------------
      a(Mice+1)=-lami/(ddzi(Mice)*l1)+ci*roi*(dhi-dhi0)/(2.*dt)
      b(Mice+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
      c(Mice+1)=a(Mice+1) + b(Mice+1) - &
      & (csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & ci*roi*ddzi(Mice)*l1/(2.*dt))
      d(Mice+1)=-Ti1(Mice+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & ci*roi*ddzi(Mice)*l1/(2.*dt))-SRi(Mice)
!-----------------------------------------------------------
      call DIFF_COEF(a,b,c,d,2,ns-1,Mice+2,Mice+ns-1,soil_indic,dt) 
      c(Mice+ns)=1.
      a(Mice+ns)=1.
      d(Mice+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
      call PROGONKA (a,b,c,d,Temp,1,Mice+ns)
      do i=1,Mice+1
        Ti2(i)=Temp(i)
      enddo
      do i=Mice+1,Mice+ns
        Tsoil2(i-Mice)=Temp(i)
      enddo
    endif
  elseif (water==1) then  
    b(1)=0.
    c(1)=1.
    d(1)=Tsurf
    call DIFF_COEF(a,b,c,d,2,M,2,M,water_indic,dt) 
    if (deepice==1) then
      a(M+1)=0.
      c(M+1)=1.
      d(M+1)=Meltpnt(Sal2(M+1))
      call PROGONKA(a,b,c,d,Temp,1,M+1) 
      do i=1,M+1
        Tw2(i)=Temp(i)
      enddo
    else
!---------------------WATER-SOIL INTERFACE------------------
      a(M+1)=-lamw(M)/(ddz(M)*h1)+cw*row0*(dhw-dhw0)/(2.*dt) + &
      & 0.5d0*cw*row0*PEMF(M)
      b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
      c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & cw*row0*ddz(M)*h1/(2.*dt))-cw*row0*PEMF(M)
      d(M+1)=-Tw1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & cw*row0*ddz(M)*h1/(2.*dt))-SR(M) - &
      & cw*row0*PEMF(M)*pt_down_f(M)
!-----------------------------------------------------------
      call DIFF_COEF(a,b,c,d,2,ns-1,M+2,M+ns-1,soil_indic,dt)
      c(M+ns)=1.
      a(M+ns)=1.
      d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
      call PROGONKA (a,b,c,d,Temp,1,M+ns)
      do i=1,M+1
        Tw2(i)=Temp(i)
      enddo
      do i=M+1,M+ns
        Tsoil2(i-M)=Temp(i)
      enddo
    endif
  else
    b(1)=0.
    c(1)=1.
    d(1)=Tsurf  
    call DIFF_COEF(a,b,c,d,2,ns-1,2,ns-1,soil_indic,dt) 
    c(ns)=1.
    a(ns)=1.
    d(ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
    call PROGONKA (a,b,c,d,Temp,1,ns)
    do i=1,ns
      Tsoil2(i)=Temp(i)
    enddo
  endif
else   
  if (water==0) then
    RETURN
  else
    if (ice==1) then
      b(1)=0.
      c(1)=1.
      d(1)=Meltpnt(Sal2(1))
      call DIFF_COEF(a,b,c,d,2,M,2,M,water_indic,dt)  
      if (deepice==1) then
        a(M+1)=0.
        c(M+1)=1.
        d(M+1)=Meltpnt(Sal2(M+1))
        call PROGONKA (a,b,c,d,Temp,1,M+1)
        do i=1,M+1
          Tw2(i)=Temp(i)
        enddo
        b(1)=0.
        c(1)=1.
        d(1)=Meltpnt(Sal2(M+1))
        call DIFF_COEF(a,b,c,d,2,Mice,2,Mice,deepice_indic,dt)
!----------------DEEPICE-SOIL INTERFACE-------------------------
        a(Mice+1)=-lami/(ddzi(Mice)*ls1)+ci*roi*(dls-dls0)/(2.*dt)
        b(Mice+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
        c(Mice+1)=a(Mice+1)+b(Mice+1) - &
        & (csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & ci*roi*ddzi(Mice)*ls1/(2.*dt))
        d(Mice+1)=-Tis1(Mice+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & ci*roi*ddzi(Mice)*ls1/(2.*dt))-SRdi(Mice)
!-----------------------------------------------------------
        call DIFF_COEF(a,b,c,d,2,ns-1,Mice+2,Mice+ns-1,soil_indic,dt)  
        c(Mice+ns)=1.
        a(Mice+ns)=1.
        d(Mice+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
        call PROGONKA (a,b,c,d,Temp,1,Mice+ns)
        do i=1,Mice+1
          Tis2(i)=Temp(i)
        enddo
        do i=Mice+1,Mice+ns
          Tsoil2(i-Mice)=Temp(i)
        enddo
      else
!---------------------WATER-SOIL INTERFACE------------------
        a(M+1)=-lamw(M)/(ddz(M)*h1)+cw*row0*(dhw-dhw0)/(2.*dt) + &
        & 0.5d0*cw*row0*PEMF(M)
        b(M+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
        c(M+1)=a(M+1)+b(M+1)-(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & cw*row0*ddz(M)*h1/(2.*dt))-cw*row0*PEMF(M)
        d(M+1)=-Tw1(M+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
        & cw*row0*ddz(M)*h1/(2.*dt))-SR(M) - &
        & cw*row0*PEMF(M)*pt_down_f(M)
!-----------------------------------------------------------
        call DIFF_COEF(a,b,c,d,2,ns-1,M+2,M+ns-1,soil_indic,dt)
        c(M+ns)=1.
        a(M+ns)=1.
        d(M+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
        call PROGONKA (a,b,c,d,Temp,1,M+ns)
        do i=1,M+1
          Tw2(i)=Temp(i)
        enddo
        do i=M+1,M+ns
          Tsoil2(i-M)=Temp(i)
        enddo
      endif
    elseif (deepice==1) then    
      b(1)=0.
      c(1)=1.
      d(1)=Meltpnt(Sal2(M+1))
      call DIFF_COEF(a,b,c,d,2,Mice,2,Mice,deepice_indic,dt)
!----------------DEEPICE-SOIL INTERFACE-------------------------
      a(Mice+1)=-lami/(ddzi(Mice)*ls1)+ci*roi*(dls-dls0)/(2.*dt)
      b(Mice+1)=-(lamsoil(1)+lamsoil(2))/(2*dzs(1))
      c(Mice+1)=a(Mice+1)+b(Mice+1) - &
      & (csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & ci*roi*ddzi(Mice)*ls1/(2.*dt))
      d(Mice+1)=-Tis1(Mice+1)*(csoil(1)*rosoil(1)*dzs(1)/(2.*dt) + &
      & ci*roi*ddzi(Mice)*ls1/(2.*dt))-SRdi(Mice)
!-----------------------------------------------------------
      call DIFF_COEF(a,b,c,d,2,ns-1,Mice+2,Mice+ns-1,soil_indic,dt)  
      c(Mice+ns)=1.
      a(Mice+ns)=1.
      d(Mice+ns)=Hflow*dzs(ns-1)/(lamsoil(ns-1)+lamsoil(ns))
      call PROGONKA (a,b,c,d,Temp,1,Mice+ns)
      do i=1,Mice+1
        Tis2(i)=Temp(i)
      enddo
      do i=Mice+1,Mice+ns
        Tsoil2(i-Mice)=Temp(i)
      enddo
    endif
  endif
endif

END SUBROUTINE T_DIFF
    



!---------------------------------------------------------------------------------


SUBROUTINE DIFF_COEF(a,b,c,d,n0,n1,m0,m1,substr,dt)

!-------------------DEFINES COEFFICIENTS FOR SOLVING A THERMAL DIFFUSIVITY-----------------
!-------------------EQUATION BY FACTORIXATION METHOD----------------------------------------

use NUMERIC_PARAMS
use PHYS_CONSTANTS2  
use DRIVING_PARAMS
use ARRAYS
implicit none

!Parameters
integer(4), parameter :: soil_indic = 1
integer(4), parameter :: ice_indic = 2
integer(4), parameter :: snow_indic = 3
integer(4), parameter :: water_indic = 4
integer(4), parameter :: deepice_indic = 5
integer(4), parameter :: water_salinity_indic = 6
integer(4), parameter :: soil_salinity_indic = 7

integer(4), parameter :: vector_length = 350

!Input variables
real(8), intent(in) :: dt

integer(4), intent(in) :: n0
integer(4), intent(in) :: n1
integer(4), intent(in) :: m0
integer(4), intent(in) :: m1
integer(4), intent(in) :: substr

real(8), intent(out), dimension(1:vector_length) :: a,b,c,d

!Common variables
integer(4) :: itop
real(8), dimension(1:ms) :: Tsn, cs
real(8), dimension(1:ms) :: lams, q
real(8) :: dz
real(8) :: AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens

common /snow_char/ Tsn,cs
common /watericesnowarr/ lams,q
common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML), &
& ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
common /SOILDAT/ dz(ms),itop

!External functions
real(8), external :: DZETA
real(8), external :: DZETAI

!Local variables
real(8) :: lam1
real(8) :: lam2
real(8) :: dzmean

integer(4) :: i
integer(4) :: j

SAVE

if (m1-m0/=n1-n0) then
  print*, 'Error: diff_coef'
  STOP
endif

SELECT CASE (substr)
!CASE 1: SOIL
  case (soil_indic)
    do i = m0, m1
      j = i - m0 + n0
      lam2 = 0.5d0*(lamsoil(j)+lamsoil(j+1))
      lam1 = 0.5d0*(lamsoil(j-1)+lamsoil(j))
      dzmean = 0.5d0*(dzs(j-1)+dzs(j))
      a(i)=-lam1/dzs(j-1)
      b(i)=-lam2/dzs(j)
      c(i)=-(lam1/dzs(j-1)+lam2/dzs(j)) - &
      & 1.d0/dt*(csoil(j)*rosoil(j)*dzmean)
      d(i)=-Tsoil1(j)/dt*(csoil(j)*rosoil(j)*dzmean)
    enddo 
!CASE 2: ICE
  case (ice_indic)
    do i = m0, m1
      j = i - m0 + n0
      c(i)=-(lami/ddzi(j-1)+lami/ddzi(j))/l1 - &
      & 1.d0/dt*0.5d0*(ddzi(j)+ddzi(j-1))*l1*ci*roi
      a(i)=-lami/(ddzi(j-1)*l1) + &
      & ci*roi*dhi*DZETAI(dble(j))/(2.d0*dt) - &
      & ci*roi*dhi0/(2.d0*dt)
      b(i)=-lami/(ddzi(j)*l1) - &
      & ci*roi*dhi*DZETAI(dble(j))/(2.d0*dt) + &
      & ci*roi*dhi0/(2.d0*dt)  
      d(i)=-Ti1(j)/dt*0.5d0*(ddzi(j)+ddzi(j-1))*l1*ci*roi + &
      & SRi(j) - SRi(j-1) 
    enddo
!CASE   3: SNOW  
  case (snow_indic)
    do i = m0, m1
      j = i - m0 + n0
      lam2 = 0.5d0*(lams(j)+lams(j+1))
      lam1 = 0.5d0*(lams(j-1)+lams(j))
      dzmean = 0.5d0*(dz(j-1)+dz(j))
      a(i)=-lam1/dz(j-1)
      b(i)=-lam2/dz(j)
      c(i)=-(lam1/dz(j-1)+lam2/dz(j)) - &
      & 1.d0/dt*(cs(j)*dens(j)*dzmean)
      d(i)=-T(j)/dt*(cs(j)*dens(j)*dzmean)
    enddo 
!CASE 4: WATER
  case (water_indic)
    do i = m0, m1
      j = i - m0 + n0
      c(i)=-(lamw(j)/ddz(j)+lamw(j-1)/ddz(j-1))/h1 - &
      & 1.d0/dt*0.5d0*(ddz(j-1)+ddz(j))*h1*cw*row0 - &
      & cw*row0*0.5d0*PEMF(j-1) + cw*row0*0.5d0*PEMF(j)
      a(i)=-lamw(j-1)/(ddz(j-1)*h1)+cw*row0*dhw*DZETA(dble(j))/ &
      & (2.d0*dt)-cw*row0*dhw0/(2.d0*dt)+cw*row0*0.5d0*PEMF(j-1)
      b(i)=-lamw(j)/(ddz(j)*h1)-cw*row0*dhw*DZETA(dble(j))/ &
      & (2.d0*dt)+cw*row0*dhw0/(2.d0*dt)-cw*row0*0.5d0*PEMF(j)
      d(i)=-Tw1(j)/dt*0.5d0*(ddz(j)+ddz(j-1))*h1*cw*row0 + &
      & SR(j) - SR(j-1) - cw*row0*PEMF(j-1)*pt_down_f(j-1) + &
      & cw*row0*PEMF(j)*pt_down_f(j)
    enddo 
!CASE 5: DEEP ICE
  case (deepice_indic)
    do i = m0, m1
      j = i - m0 + n0
      c(i)=-(lami/ddzi(j)+lami/ddzi(j-1))/ls1 - &
      & 1.d0/dt*0.5d0*(ddzi(j)+ddzi(j-1))*ls1*ci*roi
      a(i)=-lami/(ddzi(j-1)*ls1)+ci*roi*dls*DZETAI(dble(j))/ &
      & (2.d0*dt)-ci*roi*dls0/(2.d0*dt)
      b(i)=-lami/(ddzi(j)*ls1)-ci*roi*dls*DZETAI(dble(j))/ &
      & (2.d0*dt)+ci*roi*dls0/(2.d0*dt)
      d(i)=-Tis1(j)/dt*0.5d0*(ddzi(j)+ddzi(j-1))*ls1*ci*roi + &
      & SRdi(j) - SRdi(j-1)
    enddo
!CASE 6: SALINITY IN WATER
 case (water_salinity_indic)
   do i = m0, m1
     j = i - m0 + n0
     c(i)=-(lamsal(j)/ddz(j)+lamsal(j-1)/ddz(j-1))/h1 - &
     & 1.d0/dt*0.5d0*(ddz(j)+ddz(j-1))*h1  
     a(i)=-lamsal(j-1)/(ddz(j-1)*h1)+dhw*DZETA(dble(j))/(2.d0*dt) - &
     & dhw0/(2.d0*dt)
     b(i)=-lamsal(j)/(ddz(j)*h1)-dhw*DZETA(dble(j))/(2.d0*dt) + &
     & dhw0/(2.d0*dt) 
     d(i)=-Sal1(j)/dt*0.5d0*(ddz(j)+ddz(j-1))*h1   
   enddo 
!CASE 7: SALINITY IN THE SOIL  
 case (soil_salinity_indic)
   do i = m0, m1
     j = i - m0 + n0
     a(i) = -wsoil(j-1)*dt/(dzs(j-1)+dzs(j)) 
     b(i) = wsoil(j)*dt/(dzs(j-1)+dzs(j))   
     c(i) = -1.d0-wsoil(j)*dt/(dzs(j-1)+dzs(j)) + &
     & wsoil(j-1)*dt/(dzs(j-1)+dzs(j))   
     d(i) = -Sals1(j)  
   enddo 
END SELECT
  
END SUBROUTINE DIFF_COEF




SUBROUTINE SURF_CHAR(surftyp,year,month,day,hour,phi)
use DRIVING_PARAMS
use PHYS_CONSTANTS2
use ARRAYS
use ATMOS, only: &
& WIND, &
& VELFRICT_PREV
use SFCFLX, only: &
& SFCFLX_ROUGHNESS    
use PHYS_FUNC, only: &
& SINH0, &
& WATER_ALBEDO 

!----------------DEFINES SURFACE CHARACTERISTICS:-------------------
!----------------roughness,emissivity,albedo,relative humidity,-----
!----------------coefficients in Magnus formula---------------------


implicit none

integer(4), intent(in) :: surftyp
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day
real(8)   , intent(in) :: hour
real(8)   , intent(in) :: phi

real(8) roughness,emissivity,albedo,aM,bM,niu,relhums
real(8) x1,x2,x3,x4,x5

common /surface/ roughness,emissivity,albedo,aM,bM,relhums
 
SAVE 
    
!niu is air molecular viscosity
niu = 1.5d-5 !1.007E-6
select case (surftyp)
  case (1)
    albedo = albedoofsoil
    roughness = 0.05
    aM = aMagw 
    bM = bMagw
    emissivity = emissivityofsoil
    relhums = 0.5
  case (2) 
    albedo = albedoofice
    emissivity = emissivityofice
    aM = aMagi
    bM = bMagi 
    roughness = 0.00001
    relhums = 0.7
  case (3)
    albedo = albedoofsnow
    emissivity = emissivityofsnow
    aM = aMagi
    bM = bMagi 
    roughness = 0.001
    relhums = 0.7
  case (4)
    velfrict_prev = dmax1(velfrict_prev, 1.d-2)
    roughness = dmin1(dmax1(0.111*niu/velfrict_prev + &
    & 0.0144*velfrict_prev**2/g, 1.d-5),1.1d-1)
!    call SfcFlx_roughness (fetch, wind, velfrict_prev, 0.d0, &
!    & x1, x2,roughness, x4, x5)
!    roughness = 1.d-5
    aM = aMagw
    bM = bMagw
    if (varalb==0) then
      albedo = albedoofwater
    elseif (varalb==1) then
      albedo = WATER_ALBEDO( SINH0(year,month,day,hour,phi) )
      !albedo = (dirdif()*0.05/(sinh0()+0.15)+0.05)/ &
      !(1+dirdif())
    endif
    emissivity = emissivityofwater 
    relhums = 1.
end select  

END SUBROUTINE SURF_CHAR



FUNCTION DZETA(NL)
use ARRAYS, only : ddz
use DRIVING_PARAMS, only : M

!The function DZETA
!returns the dzeta coordinate [0..1]
!(non-dimensional coordinate in water layer)
!of the level nl [1..M+1].
!nl might be fractional number, however the function is opimized for
!integer or half-integer nl.

implicit none

real(8) :: DZETA

!Input variable
real(8), intent(in) :: nl

!local variables
real(8), allocatable, save :: dzeta_int(:)
real(8), allocatable, save :: dzeta_05int(:)

integer(4) :: nl_int
integer(4) :: i !Loop index

logical :: firstcall
data firstcall /.true./

if (firstcall) then
  allocate (dzeta_05int(1:M))
  dzeta_05int(1) = 0.5d0*ddz(1)
  do i = 2, M
    dzeta_05int(i) = dzeta_05int(i-1) + 0.5d0*(ddz(i-1) + ddz(i))
  enddo
  allocate (dzeta_int(1:M+1))
  dzeta_int(1) = 0.d0
  do i = 2, M+1
    dzeta_int(i) = dzeta_int(i-1) + ddz(i-1)
  enddo
endif

if (dmod(nl,1.d0) == 0.d0) then
  DZETA = dzeta_int(int(nl))
elseif (dmod(2.d0*nl,1.d0) == 0.d0) then
  DZETA = dzeta_05int(int(nl))
else
  nl_int = int(nl)
  if (nl_int>1) then
    DZETA = sum(ddz(1:nl_int-1))
  else
    DZETA = 0.
  endif
  if (nl<M+1) DZETA = DZETA + (nl-real(nl_int))*ddz(nl_int)
endif

if (firstcall) firstcall = .false.
END FUNCTION DZETA


FUNCTION DZETAI(NL)
use ARRAYS, only : ddzi
use DRIVING_PARAMS, only : Mice

!The function DZETAI
!returns the dzetai coordinate [0..1]
!(non-dimensional vertical coordinate in ice layers)
!of the level nl [1..M+1].
!nl might be fractional number, however the function is optimized
!for integer or half-integer nl.

implicit none

real(8) :: DZETAI

real(8) :: nl

real(8), allocatable, save :: dzetai_int(:)
real(8), allocatable, save :: dzetai_05int(:)

integer(4) :: nl_int
integer(4) :: i !Loop index

logical :: firstcall
data firstcall /.true./

if (firstcall) then
  allocate (dzetai_05int(1:Mice))
  dzetai_05int(1) = 0.5d0*ddzi(1)
  do i = 2, Mice
    dzetai_05int(i) = dzetai_05int(i-1) + 0.5d0*(ddzi(i-1) + ddzi(i))
  enddo
  allocate (dzetai_int(1:Mice+1))
  dzetai_int(1) = 0.d0
  do i = 2, Mice + 1
    dzetai_int(i) = dzetai_int(i-1) + ddzi(i-1)
  enddo
endif

if (dmod(nl,1.d0) == 0.d0) then
  DZETAI = dzetai_int(int(nl))
elseif (dmod(2.d0*nl,1.d0) == 0.d0) then
  DZETAI = dzetai_05int(int(nl))
else
  nl_int = int(nl)
  if (nl_int>1) then
    DZETAI = sum(ddzi(1:nl_int-1))
  else
    DZETAI = 0.
  endif
  if (nl<Mice+1) DZETAI = DZETAI + (nl-real(nl_int))*ddzi(nl_int)
endif

if (firstcall) firstcall = .false.
END FUNCTION DZETAI


REAL(8) FUNCTION VARMEAN(var,grid)
use arrays
use driving_params

!The function varmean computes
!the spatial average of the
!variable var.
!If grid = 1, than var is averaged as 
!if it has been placed in integer points of mesh  
!If grid = 2, than var is averaged as 
!if it has been placed in half-integer points of mesh  

implicit none

integer(4) grid,i
real(8), intent(in):: var(1:M+1)

varmean = 0.
if(grid==1) then
  varmean = varmean + var(1)  *0.5*ddz(1)
  varmean = varmean + var(M+1)*0.5*ddz(M)
  do i=2,M
    varmean = varmean + var(i)*0.5*(ddz(i-1)+ddz(i))
  enddo
elseif (grid==2) then
  do i=1,M
    varmean = varmean + var(i)*ddz(i)
  enddo
else
  print*, 'Illegal identifier GRID in the function varmean: STOP'
  STOP
endif

RETURN
END FUNCTION VARMEAN

