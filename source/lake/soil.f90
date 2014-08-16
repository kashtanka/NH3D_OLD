SUBROUTINE COMSOILFORLAKE

!COMSOILFORLAKE specifies parameters of soil according to soil type

use NUMERIC_PARAMS
use DRIVING_PARAMS
use ARRAYS
implicit none

real(8) :: b(1:Num_Soil)
real(8) :: psi_max(1:Num_Soil)
real(8) :: Porosity(1:Num_Soil)
real(8) :: gamma_max(1:Num_Soil)
real(8) :: lambda_max(1:Num_Soil)
real(8) :: W_0(1:Num_Soil)
real(8) :: W_m(1:Num_Soil)
real(8) :: pow
  
integer(4) :: i
integer(4) :: j
logical :: firstcall

data firstcall /.true./

SAVE
if (firstcall) then
! allocate (AL(1:ns+2),DLT(1:ns+2),DVT(1:ns+2),DL(1:ns+2),
!& ALV(1:ns+2),DV(1:ns+2),Z(1:ns+2),T(1:ns+2),WL(1:ns+2),
!& WV(1:ns+2),WI(1:ns+2),dens(1:ns+2))
endif

!I=1 ! SAND
!I=2 ! LOAMY SAND
!I=3 ! SANDY LOAM
!I=4 ! LOAM
!I=5 ! SILT LOAM
!I=6 ! SANDY CLAY LOAM
!I=7 ! CLAY LOAM
!I=8 ! SILTY CLAY LOAM
!I=9 ! SANDY CLAY
!I=10! SILTY CLAY
!I=11! CLAY
!----------------------------------------------------------------
!  NN    bpsi_max  Por    gamma_max lambda_max W_0W_m
!----------------------------------------------------------------
!   1   4.05    3.50    0.395   0.01760   0.29800   0.02    0.01
!   2   4.38    1.78    0.410   0.01560   0.13900   0.05    0.02
!   3   4.90    7.18    0.435   0.00340   0.13700   0.08    0.03
!   4   5.39    14.6    0.451   0.00069   0.06190   0.20    0.08
!   5   5.30    56.6    0.485   0.00072   0.20400   0.18    0.07
!   6   7.12    8.63    0.420   0.00063   0.03070   0.13    0.06
!   7   8.52    36.1    0.476   0.00024   0.03720   0.27    0.12
!   8   7.75    14.6    0.477   0.00017   0.01240   0.24    0.11
!   9   10.4    6.16    0.426   0.00021   0.00641   0.23    0.10
!  10   10.4    17.4    0.492   0.00010   0.00755   0.30    0.15
!  11   11.4    18.6    0.482   0.00013   0.00926   0.40    0.20
!-----------------------------------------------------------------

zsoil(1) = 0.

if(UpperLayer > 0.) then
  pow = dlog(UpperLayer/depth)/dlog(1.d0/float(ns-1))
  do i=2,ns
    zsoil(i) = (1.*i/(ns-1))**pow*depth
  end do
else
  do i=2,ns
    zsoil(i) = (1.*i/(ns-1))*depth
  end do
end if

do i = 1, ns-1
  dzs(i)=zsoil(i+1)-zsoil(i)
end do

do i=1, ns
  if (i==1) then
    dzss(i) = 0.5d0*dzs(i)
  elseif (i==ns) then
    dzss(i) = 0.5d0*dzs(i-1)
  else
    dzss(i) = 0.5d0*(dzs(i-1)+dzs(i))
  endif 
end do 

!SoilType = 7

b(1) = 4.05
b(2) = 4.38
b(3) = 4.90
b(4) = 5.39
b(5) = 5.30
b(6) = 7.12
b(7) = 8.52
b(8) = 7.75
b(9) = 10.4
b(10)= 10.4
b(11)= 11.4

psi_max( 1) = 3.50/100. 
psi_max( 2) = 1.78/100. 
psi_max( 3) = 7.18/100. 
psi_max( 4) = 14.6/100. 
psi_max( 5) = 56.6/100. 
psi_max( 6) = 8.63/100. 
psi_max( 7) = 36.1/100. 
psi_max( 8) = 14.6/100. 
psi_max( 9) = 6.16/100. 
psi_max(10) = 17.4/100. 
psi_max(11) = 18.6/100. 

Porosity( 1) = 0.395
Porosity( 2) = 0.410
Porosity( 3) = 0.435
Porosity( 4) = 0.451
Porosity( 5) = 0.485
Porosity( 6) = 0.420
Porosity( 7) = 0.476
Porosity( 8) = 0.477
Porosity( 9) = 0.426
Porosity(10) = 0.492
Porosity(11) = 0.482

gamma_max( 1) = 0.01760/100. 
gamma_max( 2) = 0.01560/100. 
gamma_max( 3) = 0.00340/100. 
gamma_max( 4) = 0.00069/100. 
gamma_max( 5) = 0.00072/100. 
gamma_max( 6) = 0.00063/100. 
gamma_max( 7) = 0.00024/100. 
gamma_max( 8) = 0.00017/100. 
gamma_max( 9) = 0.00021/100. 
gamma_max(10) = 0.00010/100. 
gamma_max(11) = 0.00013/100. 

lambda_max( 1) = 0.29800/10000. 
lambda_max( 2) = 0.13900/10000.
lambda_max( 3) = 0.13700/10000.
lambda_max( 4) = 0.06190/10000.
lambda_max( 5) = 0.20400/10000.
lambda_max( 6) = 0.03070/10000.
lambda_max( 7) = 0.03720/10000.
lambda_max( 8) = 0.01240/10000.
lambda_max( 9) = 0.00641/10000.
lambda_max(10) = 0.00755/10000.
lambda_max(11) = 0.00926/10000.

W_0( 1) = 0.02
W_0( 2) = 0.05
W_0( 3) = 0.08
W_0( 4) = 0.20
W_0( 5) = 0.18
W_0( 6) = 0.13
W_0( 7) = 0.27
W_0( 8) = 0.24
W_0( 9) = 0.23
W_0(10) = 0.30
W_0(11) = 0.40

W_m( 1) = 0.01
W_m( 2) = 0.02
W_m( 3) = 0.03
W_m( 4) = 0.08
W_m( 5) = 0.07
W_m( 6) = 0.06
W_m( 7) = 0.12
W_m( 8) = 0.11
W_m( 9) = 0.10
W_m(10) = 0.15
W_m(11) = 0.20

if(SoilType>0) then
  do j = 1, ns
    BH(j)= b(SoilType)   ! PARAMETER B, DIMENSOINLESS
    PSIMAX(j) = - psi_max(SoilType) ! SAT. WATER POTENTIAL, M.    
    POR(j)    = Porosity(SoilType)  ! POROSITY, DIMENSIONLESS
    FLWMAX(j) = gamma_max(SoilType) ! SAT. HYDR. CONDUCTIVITY, M/S
    DLMAX(j)  = lambda_max(SoilType)! SAT. WATER DIFFUSIVITY, ?KG/(M*SEC)
    WLM0(j)   = W_0(SoilType) ! MAXIMAL UNFREEZING WATER AT 0C
    WLM7(j)   = W_m(SoilType) ! MAXIMAL UNFREEZING WATER AT T<<0C
  end do
else
!do j = 1, ns
! BH(j)= BH_soil! PARAMETER B, DIMENSOINLESS
! PSIMAX(j) = PSIMAX_soil ! SAT. WATER POTENTIAL, M.    
! POR(j)    = POR_soil    ! POROSITY, DIMENSIONLESS
! FLWMAX(j) = FLWMAX_soil ! SAT. HYDR. CONDUCTIVITY, M/S
! DLMAX(j)  = DLMAX_soil  ! SAT. WATER DIFFUSIVITY, ?KG/(M*SEC)
! WLM0(j)   = WLM0_soil   ! MAXIMAL UNFREEZING WATER AT 0C
! WLM7(j)   = WLM7_soil   ! MAXIMAL UNFREEZING WATER AT T<<0C
!end do
end if

if (firstcall) firstcall=.false.   
RETURN
END SUBROUTINE COMSOILFORLAKE


SUBROUTINE SOILFORLAKE( &
& ix,iy,nx,ny,year,month,day,hour,phi, &
& extwat,extice,fetch,dt,a,b,c,d,Temp)

!SOILFORLAKE calculates profiles of temperature, liquid and solid water 
!content in the soil under a lake

use PHYS_CONSTANTS2 
use DRIVING_PARAMS
use ARRAYS
use PHYS_FUNC, only : &
& MELTPNT, &
& UNFRWAT, &
& WL_MAX 
implicit none

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
    
real(8) :: a(1:vector_length)
real(8) :: b(1:vector_length)
real(8) :: c(1:vector_length)
real(8) :: d(1:vector_length)
real(8) :: Temp(1:vector_length)

real(8) :: Hflow
real(8) :: dt
real(8) :: dzmean
real(8) :: ma
real(8) :: mi
real(8) :: wflow
real(8) :: wfhigh
real(8) :: lammoist1
real(8) :: lammoist2
real(8) :: surplus
real(8) :: Potphenergy
real(8) :: watavailtofreeze
real(8) :: filtr_low
real(8) :: Tf_old
real(8) :: dhwfsoil1
real(8) :: dhwfsoil2
real(8) :: dhwfsoil3
real(8) :: dhwfsoil4
real(8) :: max1

real(8) :: alp(10)
real(8) :: psiit(10)

integer(4) :: i
integer(4) :: j
integer(4) :: k
integer(4) :: iter
integer(4) :: iter3

logical :: soil_ts_ext_iter
logical :: firstcall

common/ts_ext_iter/ Tf_old,soil_ts_ext_iter
data psiit/6,3,7,2,5,4,8,1,0,0/ !/3,2,4,1,0,0,0,0,0,0/  
data firstcall /.true./

SAVE

!open (133, file=dir_out//'\dhwfsoil.dat', status = 'unknown')   

!PARAMETERS OF ITERATIONAL PROCESS
ma=10.5 !5.5
mi=4.5
do i=1,8 
 alp(i)=2*(mi+ma-(ma-mi)*dcos((2*psiit(i)-1)/16*pi))**(-1)
enddo
   
!do i=1, ns
! wlmax_i = POR(i)*ROW0/(rosdry(i)*(1-por(i))) 
! BB1 = BH(i) + 2.
! csoil(i) = cr+WL1(i)*CW+WI1(i)*CI
! rosoil(i) = rosdry(i)*(1-por(i))*(1+wi1(i)+wl1(i)) 
! ARG = (WL1(i)+WI1(i))/wlmax_i
! ARG = dmax1(ARG,1.d-2)
! PSI = PSIMAX(i)*(ARG)**(-BH(i))
! PF = dlog(-PSI*100.)/dLOG(10.d0)
! IF(PF>5.1) THEN
!  lamsoil(i) = 4.1E-4*418.
! ELSE
!  lamsoil(i) = dexp(-PF-2.7)*418.
! END IF
! wlmax_i = POR(i)*ROW0/(rosdry(i)*(1-por(i))) - wi1(i)*row0/roi
! ARG = (WL1(i))/wlmax_i
! ARG = dmax1(ARG,1.d-2)
! lammoist(i) = dlmax(i)*ARG**BB1
!end do

Hflow = 0

!do i=1, ns-1
! wlmax_i=WL_MAX(por(i+1),rosdry(i+1),wi1(i+1))
! if (wl1(i+1)/wlmax_i>0.98.or.wlmax_i<0.01) then
!  filtr(i) = 0
! else
!  wlmax_i = (POR(i)+POR(i+1))*ROW0/((rosdry(i)+rosdry(i+1))*&
!&  (1-(por(i)+por(i+1))/2)) - (wi1(i)+wi1(i+1))/2*row0/roi
!  filtr(i) = (flwmax(i+1)+flwmax(i))/2*((wl1(i+1)+&
!&  wl1(i))/2/wlmax_i)**(bh(i+1)+bh(i)+3)
! endif
!enddo

!SPLITTING-UP METHOD FOR TEMPERATURE EQUATION:
!STEP 1: HEAT DIFFUSION

!do i = 1, ns-1
! lammoist1 = 0.5*(lammoist(i)+lammoist(i+1))
! wsoil(i) = ( - lammoist1*(wl1(i+1)-wl1(i))/dzs(i) + &
! & filtr(i) ) *0.5*(rosdry(i)+rosdry(i+1))* &
! & (1-0.5*(por(i)+por(i+1)) )/row0
!enddo


!SPLITTING-UP METHOD FOR MOISTURE EQUATION:
!STEP 1: MOISTURE DIFFUSION

Wflow = 0.
dhwfsoil=0.
dhwfsoil1=0.
dhwfsoil2=0.
dhwfsoil3=0.
dhwfsoil4=0.

if ((l1 /= 0.and. h1 == 0.).or.ls1/=0.or. &
  & (hs1/=0.and.l1==0.and.h1==0)) then
  Wfhigh = 0
!c(1)=1
!b(1)=1
!d(1)=0 
  c(1) = -1-dt*(lammoist(1)+lammoist(2))/(dzs(1)**2)
  b(1) =   -dt*(lammoist(1)+lammoist(2))/(dzs(1)**2)
  d(1) =   -wl1(1)
endif
   
if (h1 /= 0.and.ls1 == 0) then
  c(1)=1
  b(1)=0
  d(1)=WL_MAX(por(1),rosdry(1),wi1(1))
!dhwfsoil1 = dt*(wl1(2)-wl1(1))*rosdry(1)*(1-por(1))/row0
!& /dzs(1)*(lammoist(1)+lammoist(2))/2 !-
!&  (WL_MAX(por(1),rosdry(1),wi1(1))-wl1(1)) VS,23.09.07
!& *rosdry(1)*(1-por(1))*dzss(1)/row0  VS,23.09.07
endif    
  
lammoist1 = (lammoist(ns-1)+lammoist(ns))/2
!c(ns)=-lammoist1/dzs(ns-1)
!a(ns)=-lammoist1/dzs(ns-1)
!d(ns)=Wflow
c(ns) = -1-2*dt*lammoist1/dzs(ns-1)**2 
a(ns) = -2*dt*lammoist1/dzs(ns-1)**2 
d(ns) = -wl1(ns)

do i=2,ns-1
  lammoist2 = (lammoist(i)+lammoist(i+1))/2
  lammoist1 = (lammoist(i-1)+lammoist(i))/2
  dzmean = (dzs(i-1)+dzs(i))/2
  a(i)=-dt/dzmean*lammoist1/dzs(i-1)
  b(i)=-dt/dzmean*lammoist2/dzs(i)
  c(i)=-dt/dzmean*(lammoist1/dzs(i-1) + &
  & lammoist2/dzs(i))-1
  d(i)=-wl1(i)
 !wl2(i) = wl1(i)+dt*(lammoist2/dzs(i)*(wl1(i+1)-wl1(i))-
!& lammoist1/dzs(i-1)*(wl1(i)-wl1(i-1)))/dzmean
enddo   
    
call PROGONKA(a,b,c,d,wl2,1,ns)  
    
if (h1 /= 0.and.ls1 == 0) then 
  lammoist1 = (lammoist(1)+lammoist(2))/2.
  dhwfsoil1 = -(wl2(1)*(dzs(1)/(2*dt)+lammoist1/dzs(1)) - &
  & wl2(2)*lammoist1/dzs(1)-wl1(1)*dzs(1)/(2*dt)) * &
  & rosdry(1)*(1-por(1))/row0*dt
endif 
!STEP 2 OF SPLITTING-UP METHOD FOR MOISTURE EQUATION: 
!EVOLUTION OF MOISTURE DUE TO GRAVITATIONAL INFILTRATION  
 
do i=2,ns-1
  wl3(i) = wl2(i) + dt*(filtr(i-1)-filtr(i))/dzss(i)
enddo
 
if ((h1==0.and.l1/=0).or.ls1 /= 0) then
  wl3(1) = wl2(1) + dt*(-filtr(1)+0)/dzss(1)
endif

if (h1/=0.and.ls1 == 0) then
  wl3(1) = WL_MAX(por(1),rosdry(1),wi1(1))
  dhwfsoil2 = - filtr(1)*rosdry(1)*(1-por(1))/row0*dt 
endif

filtr_low = 0
wl3(ns) = wl2(ns) + dt*(-filtr_low+filtr(ns-1))/dzss(ns)

do i=1,ns
  if (wl3(i)>WL_MAX(por(i),rosdry(i),wi1(i))) then
    surplus=wl3(i)-WL_MAX(por(i),rosdry(i),wi1(i))
    cy1:do j=1,ns
      if (WL_MAX(por(j),rosdry(j),wi1(j))-wl3(j)>0) then
        if (WL_MAX(por(j),rosdry(j),wi1(j)) - wl3(j) > &
        & surplus*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl3(j)=wl3(j)+surplus*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          surplus=0
          exit cy1
        else
          surplus=surplus-(WL_MAX(por(j),rosdry(j),wi1(j))-wl3(j)) * &
          & (rosdry(j)*(1-por(j))*dzss(j)) / &
          & (rosdry(i)*(1-por(i))*dzss(i))
          wl3(j)=WL_MAX(por(j),rosdry(j),wi1(j)) 
        endif
      endif
    enddo cy1
    dhwfsoil3 = dhwfsoil3 + surplus*rosdry(i)*(1-por(i))*dzss(i)/row0
    wl3(i) = wl_max(por(i),rosdry(i),wi1(i))
  endif
  if(wl3(i)<0) then
    if (i/=1) then 
      do j=i-1,1,-1
        if (wl3(j)>-wl3(i)*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl3(j)=wl3(j)+wl3(i)*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          goto 1
        endif
      enddo
    endif
    if (i/=ns) then 
      do j=i+1,ns
        if (wl3(j)>-wl3(i)*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl3(j)=wl3(j)+wl3(i)*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          goto 1
        endif
      enddo 
    endif
 !dhwfsoil = dhwfsoil - abs(wl4(i))*rosdry(i)*(1-por(i))*dzss(i)/row
1   wl3(i)=0.
  endif
enddo

! STEP 3 OF SPLITTING-UP METHOD : PHASE PROCESSES

20 c2: do i=1, ns
  if (wi1(i)>0.and.Tsoil2(i)>Meltpnt(Sals2(i))+0.01) then
    Potphenergy = (Tsoil2(i)-(Meltpnt(Sals2(i))+0.01)) * &
    & csoil(i)*rosoil(i)*dzss(i)
    if (Potphenergy>=wi1(i)*rosdry(i)*dzss(i)*(1-por(i))*Lwi) then
      wl4(i) = wl3(i) + wi1(i)
      wi2(i) = 0.d0
      Tsoil3(i) = Meltpnt(Sals2(i))+0.01 + &
      & (Potphenergy-(wi1(i)*rosdry(i)*dzss(i)* &
      & (1-por(i))*Lwi))/(csoil(i)*rosoil(i)*dzss(i))
    else
      wl4(i) = wl3(i) + Potphenergy / &
      & (rosdry(i)*dzss(i)*(1-por(i))*Lwi)
      wi2(i) = wi1(i) - Potphenergy / &
      & (rosdry(i)*dzss(i)*(1-por(i))*Lwi)
      Tsoil3(i) = Meltpnt(Sals2(i)) + 0.01
    endif
  else
    if (wl3(i)>=0.and.Tsoil2(i)<Meltpnt(Sals2(i))-0.01) then
      Potphenergy = -(Tsoil2(i)-Meltpnt(Sals2(i))+0.01) * &
      & csoil(i)*rosoil(i)*dzss(i)
      watavailtofreeze = (wl3(i)-unfrwat(Tsoil2(i),i))*rosdry(i) * &
      & dzss(i)*(1-por(i))*Lwi
      if (watavailtofreeze<0) then 
        !cc = csoil(i)*rosoil(i)*dzss(i)
        !cc1 = rosdry(i)*dzss(i)*(1-por(i))*Lwi
        !bb = (-Potphenergy-0.1*cc-wl3(i)*cc1)/cc
        !aa = cc1/cc
        !call phase_iter(aa,bb,i,Tsoil3(i))
        Tsoil3(i) = Tsoil2(i) + &
        & (wl3(i)-unfrwat(Tsoil2(i),i))*rosdry(i)*dzss(i)*(1-por(i))*Lwi / &
        & (csoil(i)*rosoil(i)*dzss(i))
        wi2(i) = wi1(i) + (wl3(i)-unfrwat(Tsoil2(i),i))
        wl4(i) = unfrwat(Tsoil2(i),i)
        cycle c2  
      endif
      if (Potphenergy>=watavailtofreeze) then
        Tsoil3(i) = Meltpnt(Sals2(i)) -0.01 - &
        & (Potphenergy-watavailtofreeze) / &
        & (csoil(i)*rosoil(i)*dzss(i))
        wi2(i) = wi1(i) + (wl3(i)-unfrwat(Tsoil2(i),i))
        wl4(i) = unfrwat(Tsoil2(i),i)
      else
        wl4(i) = wl3(i)-Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
        wi2(i) = wi1(i)+Potphenergy/(rosdry(i)*dzss(i)*(1-por(i))*Lwi)
        Tsoil3(i) = Meltpnt(Sals2(i)) - 0.01
      endif
    else
      wl4(i) = wl3(i)
      wi2(i) = wi1(i)
      Tsoil3(i) = Tsoil2(i)
    endif
  endif
enddo c2

max1=0.
do i=1,ns
  if (Tsoil3(i)<Meltpnt(Sals2(i))-0.01.and. &
  &  dabs(wl4(i)-unfrwat(Tsoil3(i),i))>max1) then
    max1=dabs(wl4(i)-unfrwat(Tsoil3(i),i))
  endif
enddo    
  
if (max1>0.01.and.iter3<20) then
  !wl3=(wl4+wl3)/2.
  !wi1=(wi2+wi1)/2.
  !Tsoil2=(Tsoil3+Tsoil2)/2.
  wl3=wl3+alp(iter+1)*(wl4-wl3)
  wi1=wi1+alp(iter+1)*(wi2-wi1)
  Tsoil2=Tsoil2+alp(iter+1)*(Tsoil3-Tsoil2)
  iter = iter + 1
  iter3 = iter3 + 1
  if(iter==8) iter = 0
  goto 20
else
  iter = 0
  iter3 = 0 
endif   

!evol=0

do i=1,ns
  if (wi2(i)<0) then
    wl4(i)=wl4(i)+wi2(i)
    wi2(i)=0
  endif
enddo

do i=1,ns
  if (wl4(i)>POR(i)*ROW0/(rosdry(i)*(1-por(i)))-wi2(i)*row0/roi) then
    surplus = wl4(i)-(POR(i)*ROW0/(rosdry(i)*(1-por(i))) - wi2(i)*row0/roi)
    c3:do j=1,ns
      if (POR(j)*ROW0/(rosdry(j)*(1-por(j))) - wi2(j)*row0/roi - &
      & wl4(j)>0) then
        if (POR(j)*ROW0/(rosdry(j)*(1-por(j))) - wi2(j)*row0/roi - &
        & wl4(j)>surplus*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl4(j)=wl4(j)+surplus*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          surplus=0
          exit c3
        else
          surplus=surplus-(POR(j)*ROW0/(rosdry(j)*(1-por(j))) - &
          & wi2(j)*row0/roi-wl4(j))*(rosdry(j)*(1-por(j))*dzss(j)) / &
          & (rosdry(i)*(1-por(i))*dzss(i))
          wl4(j)=POR(j)*ROW0/(rosdry(j)*(1-por(j))) - wi2(j)*row0/roi 
        endif
      endif
    enddo c3
    dhwfsoil4 = dhwfsoil4 + surplus*rosdry(i)*(1-por(i))*dzss(i)/row0
    wl4(i) = POR(i)*ROW0/(rosdry(i)*(1-por(i))) - wi2(i)*row0/roi 
  endif
  if(wl4(i)<0) then
    if (i/=1) then 
      do j=i-1,1,-1
        if (wl4(j)>-wl4(i)*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl4(j)=wl4(j)+wl4(i)*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          goto 2
        endif
      enddo
    endif
    if (i/=ns) then 
      do j=i+1,ns
        if (wl4(j)>-wl4(i)*rosdry(i)*(1-por(i))*dzss(i) / &
        & (rosdry(j)*(1-por(j))*dzss(j))) then
          wl4(j)=wl4(j)+wl4(i)*rosdry(i)*(1-por(i))*dzss(i) / &
          & (rosdry(j)*(1-por(j))*dzss(j))
          goto 2
        endif
      enddo
    endif
    !dhwfsoil = dhwfsoil - abs(wl4(i))*rosdry(i)*(1-por(i))*dzss(i)/row
2   wl4(i)=0
  endif
enddo

  !!!!! II. FREEZING AND MELTING IN CASE OF TSOIL<MELTING POINT (O C)!!!!!

    !do i=1, ns
    ! if (Tsoil3(i)<-1.) then
    !   wi3(i) = wi2(i) + (wl4(i)-unfrwat(Tsoil3(i),i))
  !   wl5(i) = unfrwat(Tsoil3(i),i)
    !   Tsoil4(i) = Tsoil3(i) + & 
    !   (wl4(i)-unfrwat(Tsoil3(i),i))*rosdry(i)*dzss(i)*(1-por(i))*Lwi/&
    !   (csoil(i)*rosoil(i)*dzss(i))
    ! else
    !   wi3(i) = wi2(i)
    !   wl5(i) = wl4(i)
    !   Tsoil4(i) = Tsoil3(i)
  ! end if
    ! if (wl5(i)>POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi) then
    !  surplus = wl5(i)-(POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi)
    !  dhwfsoil = dhwfsoil + surplus*rosdry(i)*(1-por(i))*dzss(i)/row
    !  wl5(i) = POR(i)*ROW/(rosdry(i)*(1-por(i))) - wi3(i)*row/roi 
  ! endif
    ! if(wl5(i)<0) then 
    !  dhwfsoil = dhwfsoil - abs(wl5(i))*rosdry(i)*(1-por(i))*dzss(i)/row
    !  wl5(i)=0
    ! endif
  !enddo
   
   ! evol=0
!   evol1=0

!   do i=1,ns
!    evol=evol+(Tsoil3(i)-Tsoil1(i))*rosoil(i)*csoil(i)*dzss(i)
!    evol1=evol1+(wi2(i)-wi1(i))*rosdry(i)*(1-por(i))*dzss(i)*Lwi
!   enddo

dhwfsoil=dhwfsoil1+dhwfsoil2+dhwfsoil3+dhwfsoil4
  
Tsoil1 = Tsoil3
wl1 = wl4
wi1 = wi2

if (firstcall) firstcall=.false.
RETURN
END SUBROUTINE SOILFORLAKE


SUBROUTINE SOIL_COND_HEAT_COEF

use ARRAYS, only : &
& rosdry, &
& csoil, &
& rosoil, &
& lamsoil, &
& lammoist, &
& filtr, &
& wsoil, &
& por, &
& bh, &
& wl1, &
& wi1, &
& psimax, &
& flwmax, &
& dlmax, &
& dzs
use DRIVING_PARAMS, only : &
& ns
use PHYS_CONSTANTS2, only : &
& row0, &
& cw, &
& ci, &
& roi
use PHYS_FUNC, only : &
& WL_MAX


implicit none

real(8), parameter :: cr = 0.2*4180.d0

real(8) :: wlmax_i
real(8) :: bb1
real(8) :: arg
real(8) :: psi
real(8) :: pf
real(8) :: lammoist1

integer(4) :: i

rosdry(:) = 1200.
   
do i = 1, ns
  wlmax_i = POR(i)*ROW0/(rosdry(i)*(1-por(i))) 
  BB1 = BH(i) + 2.d0
  csoil(i) = cr+WL1(i)*CW+WI1(i)*CI
  rosoil(i) = rosdry(i)*(1-por(i))*(1+wi1(i)+wl1(i)) 
  ARG = (WL1(i)+WI1(i))/wlmax_i
  ARG = dmax1(ARG,1.d-2)
  PSI = PSIMAX(i)*(ARG)**(-BH(i))
  PF = dlog(-PSI*100.)/dLOG(10.d0)
  if(PF>5.1) then
    lamsoil(i) = 4.1E-4*418.
  else
    lamsoil(i) = dexp(-PF-2.7)*418.
  end if
  wlmax_i = POR(i)*ROW0/(rosdry(i)*(1-por(i))) - wi1(i)*row0/roi
  ARG = (WL1(i))/wlmax_i
  ARG = dmax1(ARG,1.d-2)
  lammoist(i) = dlmax(i)*ARG**BB1
end do

do i=1, ns-1
  wlmax_i=WL_MAX(por(i+1),rosdry(i+1),wi1(i+1))
  if (wl1(i+1)/wlmax_i>0.98.or.wlmax_i<0.01) then
    filtr(i) = 0
  else
    wlmax_i = (POR(i)+POR(i+1))*ROW0/((rosdry(i)+rosdry(i+1)) * &
    & (1-(por(i)+por(i+1))/2)) - (wi1(i)+wi1(i+1))/2*row0/roi
    filtr(i) = (flwmax(i+1)+flwmax(i))/2*((wl1(i+1) + &
    & wl1(i))/2/wlmax_i)**(bh(i+1)+bh(i)+3)
  endif
enddo

do i = 1, ns-1
  lammoist1 = 0.5*(lammoist(i)+lammoist(i+1))
  wsoil(i) = ( - lammoist1*(wl1(i+1)-wl1(i))/dzs(i) + &
  & filtr(i) ) *0.5*(rosdry(i)+rosdry(i+1))* &
  & (1-0.5*(por(i)+por(i+1)) )/row0 
enddo 

END SUBROUTINE SOIL_COND_HEAT_COEF