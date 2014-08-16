
      FUNCTION HEATBALANCE(Tsurf,surftyp,extwat,extice,fetch, dt)
      use ATMOS,  only : &
     & Radbal      , &
     & Radbal_surf , &
     & hflux       , &
     & Elatent     , &
     & eflux       , &
     & eflux0      , &
     & eflux0_kinem, &
     & hskin
      use DRIVING_PARAMS, only : &
     & skin
      use PHYS_CONSTANTS2, only : &
     & cw, &
     & row0

      implicit none

      real(8) :: HEATBALANCE

      integer(4), parameter :: water_indic = 4

      real(8), intent(in) :: Tsurf
      real(8), intent(in) :: dt
      real(8), intent(in) :: extwat
      real(8), intent(in) :: extice
      real(8), intent(in) :: fetch
      
      integer(4), intent(in) :: surftyp

      real(8) :: dHdt
      real(8) :: convect_flux

      real(8), external :: SHORTBAL_TOPHALFLAYER

      if     (skin == 1 .and. surftyp == water_indic) then
        call TSKIN_SOLVER(Tsurf,extwat,extice,fetch,dt,hskin)
        call HEATFLUXSURF(surftyp,dt,dHdt,convect_flux)

        HEATBALANCE = dHdt-(SHORTBAL_TOPHALFLAYER(extwat)+hskin) + &
     &   eflux + convect_flux
        eflux0_kinem = hskin /(cw*row0)
      else ! if (skin == 0) then
        call HEATFLUXSURF     (surftyp,dt,dHdt,convect_flux)
        call SENSLATMOM_FLUXES(Tsurf,fetch)
        call RADBALANCE       (Tsurf,extwat,extice,surftyp)
        
        HEATBALANCE = dHdt-(Radbal-hflux-Elatent) + &
     &   eflux + convect_flux
        eflux0      = Radbal_surf -hflux-Elatent
        eflux0_kinem = eflux0/(cw*row0)
      endif

      RETURN
      END FUNCTION HEATBALANCE


      SUBROUTINE HEATFLUXSURF(surftyp,dt,dHdt,convect_flux)
      use ATMOS, only: &
     & botflux, &
     & eflux
      use ARRAYS, only: &
     & lamsoil,Tsoil2,Tsoil1,dzs,csoil,rosoil, &
     & Ti2,Ti1,ddz,ddzi, &
     & lamw,Tw2,Tw1, &
     & h1,l1,hs1,ls1, &
     & dhw0,dhi0, &
     & PEMF, pt_down_f
      use PHYS_CONSTANTS2, only: &
     & lami,ci,roi,cw,row0
      use NUMERIC_PARAMS, only: &
     & ms,ML
      use DRIVING_PARAMS, only: &
     & M

      implicit none

!     Parameters
      integer(4), parameter :: soil_indic = 1
      integer(4), parameter :: ice_indic = 2
      integer(4), parameter :: snow_indic = 3
      integer(4), parameter :: water_indic = 4

!     Input variables
      integer(4), intent(in) :: surftyp
      real(8), intent(in) :: dt

!     Output variables
      real(8), intent(out) :: dHdt
      real(8), intent(out) :: convect_flux

      real(8), dimension(1:ms):: Tsn,cs
      common /snow_char/         Tsn,cs

      real(8), dimension(1:ms):: lams,q
      common /watericesnowarr/   lams,q

      integer(4) itop 
      real(8) AL,DLT,DVT,ALLL,DL, &
     & ALV,DV,Z,T,WL,WV,WI,dens,dz
      common /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML), &
     & ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),dens(ms)
      common /SOILDAT/ dz(ms),itop

      SELECT CASE (surftyp)   
       case (soil_indic)
        eflux = -(lamsoil(1)+lamsoil(2))*0.5d0*(Tsoil2(2)-Tsoil2(1))/dzs(1)
       case (ice_indic)
        eflux = -lami*(Ti2(2)-Ti2(1))/(ddzi(1)*l1) 
!        if (l1<0.1) eflux = -lami*(Ti2(2)-Ti2(1))/(ddz*0.1) 
       case (snow_indic)
        eflux = -(lams(itop)+lams(itop+1))*0.5d0*(Tsn(itop+1)-Tsn(itop))/dz(itop)
       case (water_indic)  
        eflux = -lamw(1)*(Tw2(2)-Tw2(1))/(ddz(1)*h1) 
!        if (h1<0.1) eflux = -lamw(1)*(Tw2(2)-Tw2(1))/(ddz*0.1) 
      END SELECT

      SELECT CASE (surftyp)   
       case (soil_indic)
        dHdt=csoil(1)*rosoil(1)*dzs(1)/2.*(Tsoil2(1)-Tsoil1(1))/dt
       case (ice_indic)
        dHdt=ci*roi*ddzi(1)*l1*0.5d0*(Ti2(1)-Ti1(1))/dt + &
     &   ci*roi*dhi0*(Ti2(2)-Ti2(1))/(2.d0*dt)
       case (snow_indic)
        dHdt=cs(itop)*dens(itop)*dz(itop)*0.5d0*(Tsn(itop)-T(itop))/dt
       case (water_indic)  
        dHdt=cw*row0*ddz(1)*h1*0.5d0*(Tw2(1)-Tw1(1))/dt + &
     &   cw*row0*dhw0*(Tw2(2)-Tw2(1))/(2.d0*dt)
      END SELECT

      SELECT CASE (surftyp)
        case(soil_indic)
          convect_flux = 0.d0
        case(ice_indic)
          convect_flux = 0.d0
        case(snow_indic)
          convect_flux = 0.d0
        case(water_indic)
          convect_flux = cw*row0*PEMF(1)*(pt_down_f(1)-0.5d0*(Tw2(2)+Tw2(1)) )
      END SELECT

      botflux = -lamw(M)*(Tw2(M+1)-Tw2(M))/(ddz(M)*h1)

      RETURN
      END SUBROUTINE HEATFLUXSURF


      SUBROUTINE SENSLATMOM_FLUXES(Tsurf,fetch)
      
!     Subroutine SENSLATMOM_FLUXES calculates
!     sensible, latent heat and momentum fluxes

      use ATMOS
      use DRIVING_PARAMS
      use PHYS_CONSTANTS2
      use ARRAYS
      use SFCFLX
      implicit none

      real(8), intent(in) :: Tsurf
      real(8), intent(in) :: fetch
      
      real(8) TET2,TET1,ro,esatsurf,humsurf,bx(7),bix(11), &
     & c_u,c_t,ksw,hwave,SHF,LHF
      real(8) wind10,urel,vrel,u,v,wr
      real(8) roughness,emissivity,albedo,aM,bM,relhums
      real(8) xx
      
      data ksw /2.d0/
     
      integer(4) itdrag

      common /wind/ wind10,urel,vrel,u,v
      common /surface/ roughness,emissivity,albedo,aM,bM,relhums
      
      SAVE
      
!    (urel,vrel) is wind, relative to lake currents

      if (relwind==2) then
        urel=uwind-u
        vrel=vwind-v
      elseif (relwind==1) then
        urel=uwind
        vrel=vwind
      endif

      wr       = dmax1(dsqrt(urel**2+vrel**2),0.1d0)
      wind10   = wr*log(10./0.01)/log(2./0.01) !0.01 m is water surface roughness
      TET2     = (tempair+273.15)*(100000./pressure)**0.286
      TET1     = (Tsurf  +273.15)*(100000./pressure)**0.286  
      ro       = pressure/(Rd*(tempair+273.15))  
      esatsurf = 610.7*10.**(aM*Tsurf/(bM+Tsurf))
      humsurf  = 0.622/pressure*esatsurf*relhums
      
      if (PBLpar==0) then
!     In this case the fluxes are set to constants, specified in input file
        hflux   = sensflux0
        Elatent = 0.d0
        tau     = momflux0
      endif

      if (PBLpar==1) then
      
!     Businger-Dayer, Beljaars parameterization for exchange coefficients

        itdrag= 10
        bx(1) = wr !wind
        bx(2) = TET2
        bx(3) = TET1
        bx(4) = humair 
        bx(5) = humsurf
        bx(6) = zref
        bx(7) = roughness 
        call dragvl (bx, bix, itdrag)
        c_u = bix(8)
        c_t = bix(9)
        hflux   = -cp*ro*c_u*c_t*wr*(TET2-TET1)
        if (l1 == 0) then
!     No ice and snow cover above lake: 
!     the latent heat of evaporation is used
          Elatent = -Lwv*ro*c_u*c_t*wr*(humair-humsurf) 
        else
!     There is ice and (probably) snow cover above lake: 
!     the latent heat of sublimation is used
          Elatent = -Liv*ro*c_u*c_t*wr*(humair-humsurf)
        endif
        tau     =  ro*c_u**2*wr**2
      endif

!     Parameterization of exchange coefficients acording to (Louis, 1979)

      if (PBLpar==2) then
        call RichPBL(hflux,Elatent,tau,Tsurf,humsurf,roughness)
      endif

      if (PBLpar==3) then
!     The parameterization of fluxes from lake model Flake (Mironov et al., 2006)
        call SfcFlx_momsenlat(zref    , zref    , fetch           , &
     &    wr       , tempair + 273.15 , humair  , Tsurf + 273.15  , &
     &    pressure , l1               , tau     , hflux           , &
     &    Elatent  , xx)
        tau = - tau
      endif 
      
!     SHALLOW WATER EFFECT ON HEAT, MOISTURE AND MOMENTUM FLUXES (PANIN ET. AL., 2006)
      
      if (l1==0.and.h1/=0.and.waveenh==1.and.(.not.PBLpar==0)) then
        hwave   = 0.07*wind10**2*(g*h1/dMAX1(wind10,1.d0)**2)**(3./5.)/g
        hflux   = hflux   + hflux   * ksw*hwave/h1
        Elatent = Elatent + Elatent * ksw*hwave/h1
        tau     = tau     + tau     * ksw*hwave/h1
      endif
      cdmw = tau/(ro*wr)
      xlew = Elatent
      hw   = hflux
          
      velfrict = dsqrt(tau/ro)    
      
! CALCULATED sensible and latent heat fluxes are used in heat balance

!     SHF=hflux
!     LHF=Elatent
      
! MEASURED sensible and latent heat fluxes are used in heat balance
!     SHF = Hm
!     LHF = LEm

      RETURN  
      END SUBROUTINE SENSLATMOM_FLUXES



      SUBROUTINE RADBALANCE(Tsurf,extwat,extice,surftyp)

      use ATMOS, only: &
     & longwave    , &
     & Radbal      , &
     & Radbal_surf , &
     & surfrad     , &
     & shortwave
      use DRIVING_PARAMS, only: &
     & skin        , &
     & ifrad
      use PHYS_CONSTANTS2, only: &
     & sigma       , &
     & sabs
      use ARRAYS, only: &
     & ddz,ddzi, &
     & h1,l1,hs1,ls1
      use PHYS_FUNC, only : &
     & EXTINCT_SNOW

      implicit none

!     Parameters
      integer(4), parameter :: soil_indic = 1
      integer(4), parameter :: ice_indic = 2
      integer(4), parameter :: snow_indic = 3
      integer(4), parameter :: water_indic = 4

      real(8),    intent(in) :: Tsurf
      real(8),    intent(in) :: extwat
      real(8),    intent(in) :: extice
      
      integer(4), intent(in) :: surftyp

      real(8)          roughness,emissivity,albedo,aM,bM,relhums
      common /surface/ roughness,emissivity,albedo,aM,bM,relhums

      if     (ifrad == 1) then
        surfrad = emissivity*sigma*(Tsurf+273.15)**4
      elseif (ifrad == 0) then
        surfrad = 0.d0
      endif

      SELECT CASE(surftyp)
      CASE(soil_indic)
        Radbal = shortwave*(1.d0-albedo)+longwave-surfrad
      CASE(ice_indic)
        Radbal = shortwave*(1.d0-albedo) * &
     &   (1.d0-dexp(-extice*0.5d0*ddzi(1)*l1))+longwave-surfrad
      CASE(snow_indic)
        Radbal = longwave - surfrad
!    &   + shortwave*(1-albedo)*
!    &   (1-EXTINCT_SNOW(dens(itop))**(0.5d0*dz(itop)) ) +
      CASE(water_indic)
        if     (skin==1) then
!     In this case Radbal is the net radiation balance of the
!     surface skin layer
         Radbal      = shortwave*(1.d0-albedo)*sabs + &
     &    longwave - surfrad
!     In this case Radbal_surf is the net radiation balance at the
!     top surface of skin layer
         Radbal_surf = shortwave*(1.d0-albedo) + &
     &    longwave - surfrad
        elseif (skin==0) then
!     In this case Radbal is the net radiation balance of the
!     top half layer of the water coloumn grid, assuming that
!     the sabs part of shortwave radiation is absorbed at the top of it
         Radbal = shortwave*(1.d0-albedo)*(1.d0-(1.d0-sabs) * &
     &    exp(-extwat*0.5d0*ddz(1)*h1))+longwave-surfrad
!     In this case Radbal_surf is the net radiation balance at the
!     outer boundary of top half layer
         Radbal_surf = shortwave*(1.d0-albedo)*sabs + &
     &    longwave - surfrad
        endif
      END SELECT

      RETURN
      END SUBROUTINE RADBALANCE

      
      FUNCTION SHORTBAL_TOPHALFLAYER(extwat)

!     The function SHORTBAL_TOPHALFLAYER computes
!     the shortwave radiative balance of the top half layer
!     of water coloumn, omitting the sabs fraction,
!     absorbed by the skin layer aloft

      use ATMOS,  only: &
     & shortwave
      use ARRAYS, only: &
     & ddz, &
     & h1
      use PHYS_CONSTANTS2, only: &
     & sabs

      implicit none
      
      real(8) :: SHORTBAL_TOPHALFLAYER
      
      real(8), intent(in) :: extwat

      real(8)          roughness,emissivity,albedo,aM,bM,relhums
      common /surface/ roughness,emissivity,albedo,aM,bM,relhums

      SHORTBAL_TOPHALFLAYER = shortwave*(1.d0-albedo)*(1.d0-sabs)* &
     & (1.d0-dexp(-extwat*0.5d0*ddz(1)*h1))

      END FUNCTION SHORTBAL_TOPHALFLAYER
