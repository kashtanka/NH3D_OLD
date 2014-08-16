SUBROUTINE S_DIFF(dt)
use ARRAYS
use DRIVING_PARAMS
use ATMOS, only: &
 & Sflux0
use PHYS_FUNC, only: &
 & w_sedim 
implicit none
      
! S_diff solves salinity (mineralization) diffusion equation  

! Parameters
integer(4), parameter :: vector_length = 350
integer(4), parameter :: water_salinity_indic = 6
integer(4), parameter :: soil_salinity_indic = 7

! Local variables
real(8) dt,Sflux1,Sflux_soil_bot
real(8), dimension(1:vector_length):: a,b,c,d,Sal
      
! Sflux1 --- salinity flux at the bottom boundary       
! Sflux0 --- salinity flux at the top boundary       

Sflux1 = 0.d0
Sflux_soil_bot = 0.d0

! It is allowed for atmospheric aerosol to immediately go to 
! water instead of first contaminating the snow cover.
! if (ice == 1) Sflux0 = 0.d0

! It is assumed, that the water freezing to ice, melting from ice and snow,
! the rain are all freshwater sinks/sources. It is not true in fact, as far as
! the acid rain may occur and contaminated snow thaw. For these cases
! the present scheme must be updated.

! Defines, if gravitational sedimentation of tracer is taken into account
! sedim = 1 it is taken into account
! sedim = 0 it is neglected
! sedim = 0

if (water==1) then
  call DIFF_COEF(a,b,c,d,2,M,2,M,water_salinity_indic,dt)
  c(1)   = -lamsal(1)/(h1*ddz(1))-dhw0/(2.d0*dt) - & ! - dhw0/(2*dt) is right sign!
  & ddz(1)*h1/(2.d0*dt)
  b(1)   = -lamsal(1)/(h1*ddz(1))+dhw0/(2.d0*dt)
  d(1)   = -ddz(1)*h1*Sal1(1)/(2.d0*dt) - Sflux0
  if (deepice==1) then
!   Case water,deepice and soil; upper ice and snow are allowed       
!   Salinity diffusion in water       
    c(M+1) = -lamsal(M)/(ddz(M)*h1)-ddz(M)*h1/(2.d0*dt) - &
    & (dhw-dhw0)/(2.d0*dt) ! - (dhw-dhw0)/(2*dt) is a right sign!
    a(M+1) = -lamsal(M)/(ddz(M)*h1)+(dhw-dhw0)/(2.d0*dt)
    d(M+1) = -Sal1(M+1)*ddz(M)*h1/(2.d0*dt)+Sflux1
    call PROGONKA (a,b,c,d,Sal,1,M+1)
    Sal(:)=dmax1(Sal(:),0.d0)
    Sal2(1:M+1)=Sal(1:M+1)
!   Salinity diffusion in soil
    call DIFF_COEF(a,b,c,d,2,ns-1,2,ns-1,soil_salinity_indic,dt)
    c(1) = -1.d0-dt*wsoil(1)/dzs(1)
    b(1) = dt*wsoil(1)/dzs(1)
    d(1) = -Sals1(1)
    c(ns) = -1.d0+wsoil(ns-1)*dt/(dzs(ns-1))
    a(ns) = -wsoil(ns-1)*dt/(dzs(ns-1))
    d(ns) = -Sals1(ns)+2.d0*dt*Sflux_soil_bot/dzs(ns-1)
    call PROGONKA (a,b,c,d,Sal,1,ns)
    Sal(:)=dmax1(Sal(:),0.d0)
    Sals2(1:ns)=Sal(1:ns)
  else
!   Case water, soil; upper ice and snow are allowed       
!------WATER-SOIL INTERFACE-------------------------
    a(M+1) = 0.5*(dhw-dhw0)/dt - lamsal(M)/(ddz(M)*h1)
    b(M+1) = 0.5*wsoil(1)
    c(M+1) = -( 0.5*(dzs(1)+ddz(M)*h1)/dt + 0.5*wsoil(1) + &
    & lamsal(M)/(ddz(M)*h1)+0.5*(dhw0-dhw)/dt )
    d(M+1) = -Sal1(M+1)*0.5*(dzs(1)+ddz(M)*h1)/dt
    call DIFF_COEF(a,b,c,d,2,ns-1,M+2,M+ns-1,soil_salinity_indic,dt)
    c(M+ns) = -1.d0+wsoil(ns-1)*dt/(dzs(ns-1))
    a(M+ns) = -wsoil(ns-1)*dt/(dzs(ns-1))
    d(M+ns) = -Sals1(ns)+2*dt*Sflux_soil_bot/dzs(ns-1)
    call PROGONKA (a,b,c,d,Sal,1,M+ns)
    Sal(:)=dmax1(Sal(:),0.d0)
    Sals2(1:ns)=Sal(M+1:M+ns)
    Sal2(1:M+1)=Sal(1:M+1)
  endif
  if (sedim == 1) call SAL_SEDIM(ddz,h1,dt,Sal2)
else
! Case of the soil and the ice above     
  call DIFF_COEF(a,b,c,d,2,ns-1,2,ns-1,soil_salinity_indic,dt)
  c(1) = -1.d0-dt*wsoil(1)/dzs(1)
  b(1) = dt*wsoil(1)/dzs(1)
  d(1) = -Sals1(1)
  c(ns) = -1.d0+wsoil(ns-1)*dt/(dzs(ns-1))
  a(ns) = -wsoil(ns-1)*dt/(dzs(ns-1))
  d(ns) = -Sals1(ns)+2.d0*dt*Sflux_soil_bot/dzs(ns-1)
  call PROGONKA (a,b,c,d,Sal,1,ns)
  Sal(:)=dmax1(Sal(:),0.d0)
  Sals2(1:ns)=Sal(1:ns)
endif 
          
END SUBROUTINE S_DIFF


SUBROUTINE SAL_SEDIM(ddz,h1,dt,Sal)

! The subroutine SAL_SEDIM updates the salinity profile
! due to gravitational sedimentation

use DRIVING_PARAMS, only: &
& M
use PHYS_FUNC, only: &
& W_SEDIM

implicit none

! Input variables
real(8), intent(in) :: ddz  (M) ! Spacing of dzeta-coordinate grid
      
real(8), intent(in) :: h1 ! Lake depth, m
real(8), intent(in) :: dt ! Timestep,   sec

! Input/output variables
real(8), intent(inout) :: Sal (M+1) ! Salinity at main levels, kg/kg

! Local variables
real(8) :: a(M+1) 
real(8) :: b(M+1) 
real(8) :: c(M+1) 
real(8) :: d(M+1) 

real(8) :: w_sediment (M) !Speed of gravitational sedimentation, m/s

real(8) x ! Help variable

integer i ! loop index

logical indstab
logical ind_bound

! External functions
logical, external:: CHECK_PROGONKA

! The speed of gravitational sedimentation, positive downwards
do i = 1, M
  w_sediment(i) = W_SEDIM()
enddo

! Boundary conditions at the top boundary (dzeta = 0)
x = 0.5d0*ddz(1)*h1/dt
c(1) = - 0.5d0*w_sediment(1) - x
b(1) =   0.5d0*w_sediment(1)
d(1) = - Sal(1)*x

! Boundary conditions at the bottom boundary (dzeta = 1)
x = 0.5d0*ddz(M)*h1/dt
a(M+1) = - 0.5d0*w_sediment(M)
c(M+1) =   0.5d0*w_sediment(M) - x
d(M+1) = - Sal(M+1)*x

! The coefficients of tridiagonal matrix
do i = 2, M
  x = 0.5d0*(ddz(i-1) + ddz(i))*h1/dt
  a(i) = - 0.5d0*w_sediment(i-1)
  c(i) =   0.5d0*w_sediment(i-1) - 0.5d0*w_sediment(i) - x
  b(i) =   0.5d0*w_sediment(i)
  d(i) = - Sal(i)*x
enddo
  
ind_bound = .true.
call IND_STAB_FACT_DB(a,b,c,1,M+1,indstab,ind_bound)
call PROGONKA(a,b,c,d,Sal,1,M+1)
if (.not.indstab) then
  print*,'Info: Unstable factorization method in SAL_SEDIM'
  print*,'The accuracy flag is', CHECK_PROGONKA(M+1,a,b,c,d,Sal)
  if (.not.CHECK_PROGONKA(M+1,a,b,c,d,Sal) ) STOP
endif

RETURN
END SUBROUTINE SAL_SEDIM
