SUBROUTINE SNOWTEMP(ix,iy,nx,ny,year,month,day,hour,snowmass, &
& snowmass_init,a,b,c,d,Temp,phi,extwat,extice,fetch,dt)

! SNOWTEMP calculates temperature profile in the snow cover

use NUMERIC_PARAMS    
use ATMOS
use PHYS_CONSTANTS2
use DRIVING_PARAMS
use ARRAYS
implicit none


!-----------------------------MAIN VARIABLES-------------------------------
!   arrays: Tsn - snow temperature, C
!   lams - thermal conductivity of snow
!    
!   rofresh - density of fresh snow, kg/m**3     
!   snowmass - snow mass, kg/m**2

!   Boundary conditions:
!   t2 - air temperature (2 m), C  

real(8), parameter :: pi = 3.141592653589d0
integer(4), parameter :: vector_length = 350
 
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day

real(8)   , intent(in) :: hour
real(8)   , intent(in) :: phi
real(8)   , intent(in) :: extwat
real(8)   , intent(in) :: extice
real(8)   , intent(in) :: fetch


real(8) :: xx
real(8) :: dzsn(ms)
real(8) :: rofresh
real(8) :: t2
real(8) :: snowmass
real(8) :: pheat
real(8) :: pmelt
real(8) :: snowmass_init
real(8) :: Tf_old1
real(8) :: CCT(ML)
real(8) :: q(110)
real(8) :: lams(1:ms)
real(8) :: cs(1:ms)
real(8) :: Tsn(1:ms)
real(8) :: dt
real(8) :: dz

real(8) :: ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD, &
& ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW, &
& HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF, &
& ElatOld,HFold,PRSold,extinct
real(8) :: AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
real(8), dimension(1:vector_length) :: a, b, c, d, Temp

integer(4) :: itop
integer(4) :: i
integer(4) :: iyear
integer(4) :: imonth
integer(4) :: iday  

logical soil_ts_ext_iter

common /watericesnowarr/ lams, q
common /BL/ ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,&
     & ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,&
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,&
     & ElatOld,HFold,PRSold,extinct(ms)
common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),&
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
common /SOILDAT/ dz(ms),itop
common /ts_ext_iter/ Tf_old1,soil_ts_ext_iter 
common /snow_char/ Tsn,cs

SAVE

do i = itop,ms
  T(i) = Tsn(i)    
enddo

t2 = tempair
Erad = shortwave*(1-albedoofsnow)
Elatent = xlew

call addPrecipToSnow(t2, T(itop), snowmass, iyear, imonth, iday, &
& CCT, pmelt, pheat,dt)

call snow_calc(t2,Erad,T(itop),snowmass,iyear, &
& imonth, iday, CCt, pmelt, pheat,dt)
hs1 = hsnow

if (flag_snow_init == 1) then
! snowmass_init = snowmass  - (precip*dt*row-snmelt*dt*row-
! & Elatent/Liv*dt)
  flag_snow_init = 0
endif 

totalevaps = totalevaps + Elatent/(row0*Lwv)*dt
totalprecips = totalprecips + precip*dt
totalmelts = totalmelts + snmelt*dt

!if (hs1 == 0) hs1=0.00001

END SUBROUTINE SNOWTEMP


SUBROUTINE SNOW_COND_HEAT_COEF

use ARRAYS, only : &
& flag_snow_init, &
& hs1, &
& Ti1, &
& totalmelts, &
& totalprecips, &
& totalevaps
use ATMOS, only : &
& tempair
use PHYS_CONSTANTS2, only : &
& row0, &
& ci, &
& cw
use NUMERIC_PARAMS, only : &
& ms, &
& ML

implicit none

real(8) :: xx
real(8) :: dzsn(ms)
real(8) :: rofresh
real(8) :: CCT(ML)
real(8) :: q(110)
real(8) :: lams(1:ms)
real(8) :: cs(1:ms)
real(8) :: Tsn(1:ms)
real(8) :: dz

real(8) :: ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD, &
& ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW, &
& HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF, &
& ElatOld,HFold,PRSold,extinct
real(8) :: AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens

integer(4) :: itop
integer(4) :: i

common /watericesnowarr/ lams, q
common /BL/ ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,&
     & ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,&
     & HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,&
     & ElatOld,HFold,PRSold,extinct(ms)
common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),&
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
common /SOILDAT/ dz(ms),itop
common /snow_char/ Tsn,cs

if (flag_snow_init == 1) then
  itop = ms-2 
  hs1 = 0.02
  rofresh = 67.9 + 51.3*dexp(tempair/2.6)
  do i = itop, ms-1
    WL(i) = 0.0
    dens(i) = rofresh 
    dz(i) = 0.01
    T(i) = tempair
  end do
  wl(ms) = 0.
  dens(ms) = rofresh
  T(ms) = Ti1(1)
  dz(ms) = 0.5d0*dz(ms-1)
  do i=1,itop-1
    T(i)=0
    wl(i)=0
    dens(i)=0
    dz(i)=0
  enddo
  totalmelts=0.
  totalprecips=0.
  totalevaps=0.
end if


do i = itop, ms-2
  dzsn(i) = (dz(i)+dz(i+1))*0.5d0
enddo
do i = itop, ms-1
  !densav = (dens(i) + dens(i+1))/2
  !lams(i) =  0.419*(6.*(densav/row)**4+1.9*(densav/row)+0.05)
  lams(i) = 0.419*(6.*(dens(i)/row0)**4+1.9*(dens(i)/row0)+0.05)
  xx = wl(i)*row0/(dens(i)*dz(i)) 
  cs(i) = ci*(1-xx)+cw*xx
end do


END SUBROUTINE SNOW_COND_HEAT_COEF
