      MODULE phys_parameters
      REAL(8) z_boundary,AZ0,ALB,whc,hcond,SNCR,CCT_veg,TBEST,FFTMIN, &
      & FFQMIN,AA_veg,BB_veg,CC_veg,PSI1_veg,PSI2_veg,DTL_veg,ROOTSM_veg, &
      & CCF_veg,DZZ,HR,SoilDensity,BH_soil,PSIMAX_soil,POR_soil, &
      & FLWMAX_soil,DLMAX_soil,WLM0_soil,WLM7_soil
      INTEGER(4) VegetType
      SAVE

      contains

      SUBROUTINE DEFINE_phys_parameters()
      implicit none
!cccccccccccccccccccccc AIR cccccccccccccccccccccccc
      z_boundary = 1000        ! height of planetary boundary layer (m)
      AZ0 = 0.035              ! roughness lenght (m)
      ALB = 0.23               ! albedo (1)

!cccccccccccccccccccccc SNOW cccccccccccccccccccccccc
      whc = 0.04               ! water holding capacity of snow
      hcond = 0.01             ! hydraulic conductivity in snow (m/s)
      SNCR = 0.01              ! M OF WATER, FROM WHICH SNOW COVERS ALL THE CELL

!cccccccccccccccccccccc VEGETATION cccccccccccccccccccccccc
      CCT_veg = 0.0016     ! DEPEND. CANOPY RESIST. ON TEMPERATURE (1/K)
      TBEST = 298.15       ! LOWEST CANOPY RESIST. TEMPERATURE (K)
      FFTMIN = 0.01        ! MINIMAL RESIST. COEFF. FOR TEMPERATURE (1)
      FFQMIN = 0.01        ! MINIMAL RESIST. COEFF. FOR HUMIDITY (1)

      VegetType = 7        ! 1 -- 13
!cc  OR:
      AA_veg = 2582.0                    ! FOR LEAF RESISTANCE, J/M**3
      BB_veg = 1.09                      ! FOR LEAF RESISTANCE, W/M**2
      CC_veg = 110.0                     ! MINIMAL LEAF RESIST., S/M
      PSI1_veg = - 120.                  ! WATER POT. OF WILT BEGIN, M
      PSI2_veg = - 230.                  ! WATER POT. OF PERMANENT WILT, M
      DTL_veg = 3.0                      ! MAX. LAI, DIMENSIONLESS
      ROOTSM_veg = 0.5                   ! MAX. ROOTS LENGTH, M
      CCF_veg = 0.0238*1000./0.622       ! DEPEND. CANOPY RESIST. ON HUMIDITY (1/mb)

!cccccccccccccccccccccc SOIL cccccccccccccccccccccc
      DZZ = 5.               ! ACTUAL ST. DEV. OF TOPOGR. FOR RUNOFF CALC. (m)
      HR = 0.08              ! LAYER OF SOIL FROM WHICH EVAPORATION EXISTS (m)
      SoilDensity = 1200.    ! kg/m^3

!         SoilType = 0         ! 1 -- 11
!cc  OR:
      BH_soil     = 7.12             ! PARAMETER B, DIMENSIONLESS
      PSIMAX_soil = - 0.0863         ! SAT. WATER POTENTIAL, M
      POR_soil    = 0.401            ! POROSITY, DIMENSIONLESS
      FLWMAX_soil = 0.0000063        ! SAT. HYDR. CONDUCTIVITY, M/S
      FLWMAX_soil = 0.00002          ! SAT. HYDR. CONDUCTIVITY, M/S
      DLMAX_soil  = 0.000003070      ! SAT. WATER DIFFUSIVITY, M**2/S
      WLM0_soil   = 0.13             ! MAXIMAL UNFREEZING WATER AT 0C, kg/kg
      WLM7_soil   = 0.06             ! MAXIMAL UNFREEZING WATER AT T<<0C, kg/kg

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END SUBROUTINE DEFINE_phys_parameters
      END MODULE phys_parameters


      MODULE phys_constants
      REAL(8) grav,PL,PLI,GPA,SIGSB,APPA,R,CP,RHOW,RHOI,CR,CWAT,CICE,&
      & CSNOW,WMG,WLMMX,CEFF,CA,DZL,DZG,BMIN,DMIN,DMAX,D,WIINF,ZRM,&
      & ZRMM,TRM,TOMIN,FLXMIN,CK,cpa
      SAVE

      contains

      SUBROUTINE DEFINE_phys_constants()
      implicit none
!cccccccccccccccccccccc COMMON cccccccccccccccccccccccc
      grav = 9.81               ! acceleration due to gravity (m/s**2)
      PL = 600.*4180.           ! latent heat for evapor./condens. (cal/g)    
      PLI = 80.*4180.           ! latent heat for freezing/melting (cal/g)

!cccccccccccccccccccccc AIR cccccccccccccccccccccccc
      CPA = 1008.               ! air specific heat content (Dg /(kg deg))
      SIGSB = 5.67E-8           ! Stephan-Boltzman const.(Dg/(m**2 sec K**4))
      APPA = 0.286              ! kappa = c_v/c_p
      R = 287.                  ! universal gas constant
      CP = 1004.

!cccccccccccccccccccccc SOIL and SNOW cccccccccccccccccccccccc
      RHOW = 1000.              ! the liquid water density, kg/m**3
      RHOI = 916.8              ! the ice density, kg/m**3
      CR = 1008.                ! the soil specific heat content, Dg/(kg K)
      CWAT = 4200.              ! the water specific heat content, Dg/(kg K)
      CICE = 2100.              ! the ice specific heat content, Dg/(kg K)
      CSNOW = 2100.             ! the snow specific heat content, Dg/(kg K)
      WMG = 168.

!cccccccccccccccccccccc HYDROLOGY cccccccccccccccccccccc
      WLMMX = 5.E-4        ! M. OF WATER. MINIMAL SKIN LAYER CAPACITY
      CEFF = 0.5           ! 1 FRACTION OF INTERCEPTED BY SKIN LAYER PRECIP.
      CA = 1.              ! 1 FRACTION OF THE CELL, COVERED BY PRECIP.
      DZL = 100.           ! M. MINIMAL ST. DEV. OF TOPOGR. FOR RUNOFF CALC.
      DZG = 1000.          ! M. MAXIMAL ST. DEV. OF TOPOGR. FOR RUNOFF CALC.
      BMIN = 0.01          ! 1 MINIMAL B FOR RUNOFF CALC.
      DMIN = 2.8E-10        ! M/S. COEFF. FOR SLOW DRAINAGE
      DMAX = 2.8E-8        ! M/S. COEFF. FOR FAST DRAINAGE
      D = 1.5              ! 1 COEFF. FOR FAST DRAINAGE
      WIINF = 0.30         ! 1 SOIL ICE, FROM WHICH INFILTRATION BEGIN
      ZRM = 0.1           ! M. LAYER OF SOIL WITH LARGE AMOUNT OF ROOTS
      ZRMM = 0.2           ! M. LAYER OF SOIL FOR CALCULATION OF LAI AND VEG
      TRM = 5.0            ! C. CRITICAL SOIL TEMPERATURE FOR LAI AND VEG
      TOMIN = 0.1          ! C. CRITICAL TEMPERATURE FOR WATER ACCESS FOR CANOPY
      FLXMIN = 1.05        ! CAL/(M**2*S) ACCURACY OF BALANCE AT THE SURFACE
      CK = 0.9             ! 1 COEFF. FOR CANOPY RESIST.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      END SUBROUTINE DEFINE_phys_constants
      END MODULE phys_constants

 
      MODULE PHYS_CONSTANTS2
      
      real(8) :: g
      real(8) :: omega
      real(8) :: roa0
      real(8) :: row0
      real(8) :: roi
      real(8) :: cw
      real(8) :: ci
      real(8) :: Lwi
      real(8) :: Lwv
      real(8) :: Liv
      real(8) :: kappa
      real(8) :: roughnessofwater
      real(8) :: albedoofwater
      real(8) :: emissivityofwater 
      real(8) :: sigma
      real(8) :: cp
      real(8) :: aMag
      real(8) :: bMag
      real(8) :: Rd
      real(8) :: lami
      real(8) :: lamw0
      real(8) :: albedoofice
      real(8) :: albedoofsnow
      real(8) :: emissivityofice 
      real(8) :: sabs
      real(8) :: emissivityofsnow
      real(8) :: aMagw
      real(8) :: aMagi
      real(8) :: bMagw
      real(8) :: bMagi
      real(8) :: emissivityofsoil
      real(8) :: albedoofsoil
      real(8) :: alsal
      
!     Derived constants
      real(8) :: ci_m_roi
      real(8) :: cw_m_row0
      real(8) :: row0_m_Lwi
      real(8) :: roi_m_Lwi
      
      SAVE
      
      contains
      SUBROUTINE DEFINE_PHYS_CONSTANTS2 
      implicit none

      g     = 9.814    ! acceleration due to gravity,               m/s**2
      omega = 7.29d-5  ! angular velocity of the Earth's rotation,  1/s
      sigma = 5.67d-8  ! Stefan-Boltzman constant,                  W/(m**2*K**4)

      roa0  = 1.273    ! reference density of air,                 kg/m**3
      row0  = 1000.    ! reference density of water,                kg/m**3
      roi   = 917.     ! density of ice,                            kg/m**3

      cw    = 3990.    ! heat capacity of water,                    J/(kg*K)
      ci    = 2150.    ! heat capacity of ice,                      J/(kg*K)
      cp    = 1005.    ! heat capacity of air at constant pressure, J/(kg*K)

      lami  = 2.2      ! molecular conductivity of ice,             J/(m*s*K)
      lamw0 = 0.561    ! molecular conductivity of water,           J/(m*s*K)

      Lwi   = 267900.  ! latent heat of freezing,                   J/kg
      Lwv   = 2501000. ! latent heat of evaporation,                J/kg
      Liv   = 2834600. ! latent heat of sublimation,                J/kg

      kappa = 0.38     ! Karman constant,                           n/d 
      Rd    = 287.05   ! gas constant for dry air,                  J/(kg*K)        

      aMagw = 7.6326   ! coefficient for Magnus formula for water surface
      bMagw = 241.9    ! coefficient for Magnus formula for water surface
      aMagi = 9.5      ! coefficient for Magnus formula for ice   surface
      bMagi = 265.5    ! coefficient for Magnus formula for ice   surface

      albedoofwater = 0.07 ! albedo of water surface,               n/d ! 0.07
      albedoofice   = 0.50 ! albedo of ice   surface,               n/d ! 0.35
      albedoofsnow  = 0.80 ! albedo of snow  surface,               n/d ! 0.62
      albedoofsoil  = 0.3  ! albedo of soil  surface,               n/d

      emissivityofwater = 0.99 ! emissivity of water surface,       n/d ! 0.95
      emissivityofice   = 0.99 ! emissivity of ice   surface,       n/d ! 0.95
      emissivityofsnow  = 0.99 ! emissivity of snow  surface,       n/d
      emissivityofsoil  = 0.9  ! emissivity of soil  surface,       n/d

      sabs         = 0.01 !0.4 ! the part of solar radiation, absorbed at the water surface
      alsal        = 1.  ! the ratio of eddy diffusivity for salinity to those of heat
      
      ci_m_roi  = ci * roi
      cw_m_row0 = cw * row0
      row0_m_Lwi = row0 * Lwi
      roi_m_Lwi = roi * Lwi

      END SUBROUTINE DEFINE_PHYS_CONSTANTS2
      END MODULE PHYS_CONSTANTS2


      MODULE NUMERIC_PARAMS
       integer(4), PARAMETER:: KL = 1,ML = 131,MS = 100,NT = 52592, &
       & Num_Soil = 11,Num_Veget = 13 
       real(8), parameter:: dznorm = 0.05, dzmin = 0.01, &
       & UpperLayer = 0.008 
       
       real(8), parameter :: min_ice_thick = 0.01 
       real(8), parameter :: min_water_thick = 0.01
       real(8), parameter :: T_phase_threshold = 1.d-5

!         ML --- total number of layers
!         MS --- maximal number of layers in snow
!         dznorm --- standard layer depth in snow (m)
!         dzmin --- min layer depth in snow (m)
!         depth --- depth of the lowest layer in soil (m)
!         UpperLayer --- upper soil layer depth (m)

      END MODULE NUMERIC_PARAMS
      

      MODULE ATMOS

      real(8) :: tempair
      real(8) :: relhum
      real(8) :: humair
      real(8) :: pressure
      real(8) :: hw
      real(8) :: xlew
      real(8) :: cdmw
      real(8) :: uwind
      real(8) :: vwind
      real(8) :: wind
      real(8) :: longwave
      real(8) :: Radbal
	  real(8) :: Radbal_prev
      real(8) :: hflux
	  real(8) :: hflux_prev
      real(8) :: surfrad
      real(8) :: hbal
      real(8) :: zref
      real(8) :: botflux
      real(8) :: cdmw2
      real(8) :: shortwave
      real(8) :: precip
      real(8) :: Sflux0
      real(8) :: Elatent
	  real(8) :: elatent_prev
      real(8) :: Radbal_surf
      real(8) :: eflux
      real(8) :: eflux0
      real(8) :: eflux0_kinem
      real(8) :: hskin
      real(8) :: tau
      real(8) :: velfrict
      real(8) :: velfrict_prev
      real(8) :: cdm1
      real(8) :: turb_density_flux
     
      SAVE 

      END MODULE ATMOS


      MODULE DRIVING_PARAMS
      
! Switches: 1. PBL parameterization
!              PBLpar =1 (Businger-Dayer formulas (Monin-Obukhov theory) for exchange coefficients)
!              PBLpar =11(Businger-Dayer formulas with shallow water correction after Panin(19..))      
!              PBLpar =2 (formulation from NH3d)
!              PBLpar =21(formulation from NH3d with shallow water correction after Panin(19..))
!           2. Relative to water currents wind
!              relwind =1 (relative wind is off)
!              relwind =2 (relative wind is on)
!         3. Turbulent mixing parameterization
!              Turbpar =1 (analytyical profile of turbulent conductivity)
!              Turbpar =2 ("E-eps" parameterization of turbulent conductivity)        
!              Turbpar =3 : Nickuradze (NICK) formulation: Rodi (1993)    
!              Turbpar =4 : Parabolic (PARAB) formulation: Engelund (1976)
!              Turbpar =7 : RNG (re-normalization group) formulation: Simoes (1998)     
    
      real(8) :: dt_out
      real(8) :: sensflux0
      real(8) :: momflux0
      real(8) :: kwe
      real(8) :: d_surf
      real(8) :: d_bot
      real(8) :: depth

!     The group of tributaries characteristics

      integer(4) :: N_tribin
      integer(4) :: N_tribout
      integer(4) :: tribheat
      real(8), allocatable :: &
      & U_tribin     (:,:), &
      & U_tribout    (:,:), &
      & T_tribin     (:,:), &
      & width_tribin (:,:), &
      & width_tribout(:,:)  
      real(8), allocatable :: zTinitprof(:)
      real(8), allocatable :: Tinitprof(:)

!     In perspective, it is expected, that area_lake be an array,
!     as the horizontal section of a lake varies with depth

      integer, parameter :: nunit = 96
      
      integer(4) :: runmode
      integer(4) :: PBLpar
      integer(4) :: waveenh
      integer(4) :: relwind
      integer(4) :: Turbpar
      integer(4) :: stabfunc
      integer(4) :: varalb
      integer(4) :: skin
      integer(4) :: massflux
      integer(4) :: ifrad
      integer(4) :: sedim
      integer(4) :: M
      integer(4) :: Mice
      integer(4) :: ns
      integer(4) :: nstep_keps ! The number of timesteps of k-epsilon parmeterization per on model timestep
      integer(4) :: nscreen
      integer(4) :: SoilType
      integer(4) :: Tinitlength 

      integer(4) :: monthly_out
      integer(4) :: daily_out
      integer(4) :: hourly_out
      integer(4) :: time_series
      integer(4) :: everystep
      integer(4) :: turb_out
      integer(4) :: scale_output
      integer(4) :: assim
      integer(4) :: as_window
      integer(4) :: error_cov

      character(len=40) :: path
      character(len=60) :: dataname,setupfile
      
      real(8),    external :: getvarval
      integer(4), external :: igetvarval

      SAVE
      
      contains
      SUBROUTINE DEFINE_DRIVING_PARAMS( &
 & h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,&
     &         init_T_1	  )
	 
      
!     DEFINE_driving_params reads files with driving parameters

      implicit none
	  real(8),intent(out)::h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1 
	 integer(8),intent(out):: init_T_1
      character(len=200) :: line
      logical :: firstcall
      data firstcall, line/.true., 'begin file'/
      
      if (firstcall) then            
       open  (nunit,file='setup_lake.dat',status='old')
       read  (nunit,'(a)') line
       read  (nunit,'(a)') setupfile
       close (nunit)
       open  (nunit,file=setupfile,status='old')
       do while (line /= 'end')
        read (nunit,'(a)') line
        call readpar(line,&
 & h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,&
     &         init_T_1	  )
       enddo
       close (nunit)
	   close(nunit)
       firstcall = .false. 
      endif

      END SUBROUTINE DEFINE_driving_params


      SUBROUTINE READPAR(LINE,&
 & h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1,&
     &         init_T_1	  )
      implicit none
	   real(8),intent(out)::h10_1,l10_1,ls10_1,hs10_1,Ts0_1,&
     &          Tb0_1, &
     &          Tbb0_1,h_ML0_1,extwat_1,extice_1,kor_1,trib_inflow_1,&
     &         Sals0_1,Salb0_1,fetch_1,phi_1,lam_1,us0_1,vs0_1,Tm_1,&
     &         alphax_1,alphay_1,a_veg_1,c_veg_1,h_veg_1,area_lake_1 
	 integer(8),intent(out):: init_T_1
	  

      character(len=200) :: line
      integer(4) :: i
      integer(4) :: n2
      integer(4) :: n1
      integer(4) :: ncol
      real(8), allocatable :: work(:,:)
      
      SAVE

      i=1
      if (line(1:1)==' '.or.line(1:1)==char(1)) then
        do while (line(i:i)==' '.or.line(i:i)==char(1))
          i=i+1
        enddo
        n1=i
        do while (line(i:i)/=' '.and.line(i:i)/=char(1))
          i=i+1
        enddo
        n2=i-1
      else
        n1=1
        do while (line(i:i)/=' '.and.line(i:i)/=char(1))
          i=i+1
        enddo
        n2=i-1
      endif
      
      if (line(n1:n1)/='#') then

!     The group of controls for model physics

       SELECT CASE (line(n1:n2))
       CASE ('PBLpar')
        PBLpar   = igetvarval(n1,n2,line,'Param. of surf. fluxes   ')
       CASE ('runmode')
        runmode  = igetvarval(n1,n2,line,'Mode of model run        ')
       CASE ('relwind')
        relwind  = igetvarval(n1,n2,line,'Switch for relative wind ')
       CASE ('waveenh')
        waveenh  = igetvarval(n1,n2,line,'Switch for shal.wat.corr.')
       CASE ('kwe')
        kwe      = getvarval (n1,n2,line,'The wave break TKE-effect') 
       CASE ('Turbpar')
        Turbpar  = igetvarval(n1,n2,line,'Param. of turbulence     ')
       CASE ('stabfunc')
        stabfunc = igetvarval(n1,n2,line,'Stability function is    ') 
       CASE ('varalb')
        varalb   = igetvarval(n1,n2,line,'Switch for variab. albedo')
       CASE ('soiltype')
        soiltype = igetvarval(n1,n2,line,'Soil type                ')
       CASE ('soil_depth')
        depth    = getvarval (n1,n2,line,'Soil depth, m            ')
       CASE ('skin')
        skin     = igetvarval(n1,n2,line,'Switch for skin at surf. ')
       CASE ('sedim')
        sedim    = igetvarval(n1,n2,line,'Switch for sedimentation ') 
       CASE ('massflux')
        massflux = igetvarval(n1,n2,line,'Switch for massflux conv.')
       CASE ('ifrad')
        ifrad    = igetvarval(n1,n2,line,'Switch for radiation     ')
       CASE ('sensflux0')
        sensflux0= getvarval (n1,n2,line,'Const. sens. flux, W/m**2')
       CASE ('momflux0')
        momflux0 = getvarval (n1,n2,line,'Const. mom. flux,  N/m**2')

       CASE ('path')
        read (line((n2+1):100),*) path
        print*, 'path=', path

!     The group of numerical scheme parameters
       CASE ('nstep_keps')
        nstep_keps &
        &        = igetvarval(n1,n2,line,'Number of k-eps timesteps')
       CASE ('M')
        M        = igetvarval(n1,n2,line,'Number of layers in water')
       CASE ('Mice')
        Mice     = igetvarval(n1,n2,line,'Number of layers in ice  ')     
       CASE ('ns')
        ns       = igetvarval(n1,n2,line,'Number of layers in soil ')
       CASE ('d_surf')
        d_surf   = getvarval (n1,n2,line,'Surface zooming,     n/d ')
       CASE ('d_bot')
        d_bot    = getvarval (n1,n2,line,'Bottom zooming,      n/d ')

!     The group of initial conditions
       CASE ('h10')
        h10_1      = getvarval (n1,n2,line,'Init. depth of lake, m   ')
       CASE ('l10')
        l10_1      = getvarval (n1,n2,line,'Init. ice thickness, m   ')
       CASE ('ls10')
        ls10_1     = getvarval (n1,n2,line,'Init. deep ice thick., m ')
       CASE ('hs10')
        hs10_1     = getvarval (n1,n2,line,'Init. snow depth, m      ')
       CASE ('Ts0')
        Ts0_1      = getvarval (n1,n2,line,'Init. temp. at lake surf.')
       CASE ('Tb0')
        Tb0_1      = getvarval (n1,n2,line,'Init. temp. at lake bot. ')
       CASE ('Tbb0')
        Tbb0_1     = getvarval (n1,n2,line,'Init. temp. at soil bot. ')
       CASE ('h_ML0')
        h_ML0_1    = getvarval (n1,n2,line,'Init. mixed layer depth  ')
       CASE ('Tm')
        Tm_1       = getvarval (n1,n2,line,'Mean coloumn temperature ')
       CASE ('init_T')
        init_T_1   = igetvarval(n1,n2,line,'Switch for temp. init.   ')
       CASE ('Sals0')
        Sals0_1    = getvarval (n1,n2,line,'Init. surface salinity   ')
       CASE ('Salb0')
        Salb0_1    = getvarval (n1,n2,line,'Init. bottom salinity    ')
       CASE ('us0')
        us0_1      = getvarval (n1,n2,line,'Init. surface X-speed    ')
       CASE ('vs0')
        vs0_1      = getvarval (n1,n2,line,'Init. surface Y-speed    ')
   

       CASE ('T_profile')
        Tinitlength &
        &        =igetvarval (n1,n2,line,'Init. T-profile length   ') 
        allocate (Tinitprof  (1:Tinitlength) )
        allocate (zTinitprof (1:Tinitlength) )
        ncol = 2
        allocate (work(Tinitlength,ncol))
        call READPROFILE(nunit,Tinitlength,ncol,work) 
        zTinitprof(:) = work (:,1)
        Tinitprof (:) = work (:,2)
        deallocate(work)
          
!      The group of output controls

       CASE ('nscreen')
        nscreen  = igetvarval(n1,n2,line,'Screen output interval   ')        
       CASE ('monthly')
        monthly_out &
        &        = igetvarval(n1,n2,line,'Monthly output switch    ')
       CASE ('daily')
        daily_out= igetvarval(n1,n2,line,'Daily output switch      ')
       CASE ('hourly')
        hourly_out &
        &        = igetvarval(n1,n2,line,'Hourly output switch     ')
       CASE ('everystep')
        everystep &
        &        = igetvarval(n1,n2,line,'Everystep output switch  ')
       CASE ('dt_out')
        dt_out   = getvarval (n1,n2,line,'Output interval, hours   ')
       CASE ('time_series')
        time_series &
        &        = igetvarval(n1,n2,line,'Time series output switch')
       CASE ('turb_out')
        turb_out = igetvarval(n1,n2,line,'Switch turb. char. output')
       CASE ('scale_output')
        scale_output &
        &        = igetvarval(n1,n2,line,'Switch turb. out. scaling')

!     The group of physical properties
       CASE ('extwat')
        extwat_1   = getvarval (n1,n2,line,'The rad extinct in water ')
       CASE ('extice')
        extice_1   = getvarval (n1,n2,line,'The rad extinction in ice')
       CASE ('alphax')
        alphax_1   = getvarval (n1,n2,line,'The x-slope angle        ')
       CASE ('alphay')
        alphay_1   = getvarval (n1,n2,line,'The y-slope angle        ')
       CASE ('c_veg')
        c_veg_1    = getvarval (n1,n2,line,'c_veg                    ')
       CASE ('a_veg')
        a_veg_1    = getvarval (n1,n2,line,'a_veg                    ')
       CASE ('h_veg')
        h_veg_1    = getvarval (n1,n2,line,'The vegetation height    ')
       CASE ('kor')
        kor_1      = getvarval (n1,n2,line,'The Coriolis parameter   ')
       CASE ('phi')
        phi_1    = getvarval (n1,n2,line,'Latitude, deg            ')
       CASE ('lam')
        lam_1    = getvarval (n1,n2,line,'Longitude, deg           ')
       CASE ('fetch')
        fetch_1    = getvarval (n1,n2,line,'The fetch, m             ')

!     The group of data assimilation controls

       CASE ('assim')
        assim    = igetvarval(n1,n2,line,'The assimilation switch  ')
       CASE ('as_window')
        as_window= igetvarval(n1,n2,line,'Assimilation window      ')
       CASE ('error_cov')
        error_cov= igetvarval(n1,n2,line,'The switch for err_cov   ')

!     The group of parameters of inflows and outflows
       CASE ('trib_inflow')
        trib_inflow_1  = getvarval (n1,n2,line,'Tributary inflow, m/yr   ')
       CASE ('area_lake')
        area_lake_1= getvarval (n1,n2,line,'The area of lake, m**2   ')

       CASE ('tribheat')
        tribheat = igetvarval(n1,n2,line,'Switch for trib. heat    ')
       CASE ('inflowprof')
          N_tribin = 1
          allocate (width_tribin  (1:N_tribin ,1:M+1) )
          allocate (U_tribin      (1:N_tribin ,1:M+1) )
          allocate (T_tribin      (1:N_tribin ,1:M+1) )
          ncol = 4
          allocate (work(M+1,ncol))
          call READPROFILE(nunit,M+1,ncol,work) 
!         In current version there is only ONE inflow     
          width_tribin(1,:) = work (:,2)
          U_tribin    (1,:) = work (:,3)
          T_tribin    (1,:) = work (:,4)
          deallocate(work)
       CASE ('outflowprof')
          N_tribout = 1
          allocate (width_tribout  (1:N_tribout ,1:M+1) )
          allocate (U_tribout      (1:N_tribout ,1:M+1) )
          ncol = 3
          allocate (work(M+1,ncol))
          print*, 'The profiles in inflow are done'
          call READPROFILE(nunit,M+1,ncol,work) 
!         In current version there is only ONE outflow 
          width_tribout(1,:) = work (:,2)
          U_tribout    (1,:) = work (:,3)
          print*, 'The profiles in ouflow are done'
          deallocate(work)
        
       CASE DEFAULT
        if (.not.line(n1:n2)=='end') then
         print*,'Unknown keyword [',line(n1:n2),'] in setup file: STOP'
         STOP
        endif
       END SELECT

      endif

      END SUBROUTINE READPAR

	  


      SUBROUTINE READPROFILE(iunit,N_rows,N_coloumns,work)
      implicit none

!     Input variables:
      integer(4), intent(in) :: iunit
      integer(4), intent(in) :: N_rows
      integer(4), intent(in) :: N_coloumns

!     Output variables:
      real(8), intent(out) :: work(N_rows, N_coloumns)

!     Local variables
      integer(4) :: i
      integer(4) :: j

      do i=1, N_rows
        read(iunit,*) (work(i,j), j=1,N_coloumns)
      enddo

      END SUBROUTINE READPROFILE


      END MODULE DRIVING_PARAMS




      MODULE ARRAYS
      
!     MODULE ARRAYS contains allocatable arrays of the model  

      use DRIVING_PARAMS

      integer(4) snow,ice,water,deepice,nstep 
      real(8) time,tempairf,relhumf,humairf,pressuref, &
      & windf,shortwavef,longwavef,precipf,time_old,tempair_old, &
      & relhum_old,humair_old,pressure_old,wind_old, &
      & shortwave_old,longwave_old,precip_old,Erad, &
      & wdir_old,wdirf,dep_av
      real(8) zinv
      integer(4) time_int,time_int_old,man,par,flag_snow_init
      real(8) totalevaps,totalmelts,totalprecips
      real(8), allocatable :: Tskin(:), Tskin_2d(:,:)
      real(8) h1,l1,hs1,ls1
      real(8) dhw,dhw0,dhi,dhi0,dls,dls0,dhwfsoil
      real(8) Buoyancy0,H_mixed_layer, w_conv_scale,T_conv_scale
      real(8) :: SR_botsnow
      real(8) :: dt_keps
      real(8) :: deltaskin
      real(8), allocatable:: ddz    (:) , dz_full  (:), &
      &                      z_full (:) , z_half   (:)
      real(8), allocatable :: ddzi(:)
      real(8), allocatable :: Tw1(:),Tw2(:),Ti1(:),Ti2(:),Tis1(:), &
      & Tis2(:),RS(:),SR(:),lamw(:),lamsal(:),Sal1(:),Sal2(:), &
      & Sals1(:),Sals2(:),SRi(:),SRdi(:)
      real(8), allocatable:: Tsoil1(:),rosoil(:),csoil(:), &
      & lamsoil(:), dzs(:), dzss(:), wi1(:), wl1(:), Tsoil2(:), &
      & zsoil(:), wsoil(:) 
      real(8), allocatable:: S(:),Gen(:),F(:),TKE_turb_trans(:)
      real(8), allocatable:: KT(:),k2(:),u1(:),v1(:),w1(:),E1(:),eps1(:)
      real(8), allocatable:: WLM0(:),WLM7(:),bH(:),PSIMAX(:),POR(:),FLWMAX(:),DLMAX(:)
      real(8), allocatable:: KC(:),KLengT(:),RSR(:) 
      real(8), allocatable::l1_2d(:,:),h1_2d(:,:),ls1_2d(:,:),hs1_2d(:,:),time_2d(:,:) 
      real(kind(0.d0)), allocatable:: dep_2d(:,:),cdm2(:,:)
      real(8), allocatable::u_2d(:,:,:),v_2d(:,:,:),Tsoil1_2d(:,:,:), &
      & wi1_2d(:,:,:),wl1_2d(:,:,:),Tw1_2d(:,:,:),Ti1_2d(:,:,:), &
      & Tis1_2d(:,:,:),dz_2d(:,:,:),T_2d(:,:,:),wl_2d(:,:,:), &
      & dens_2d(:,:,:),E_2d(:,:,:),eps_2d(:,:,:),dhwfsoil_2d(:,:), &
      & Sal1_2d(:,:,:),Sals1_2d(:,:,:),Elatent_2d(:,:),dhw_2d(:,:), &
      & dhw0_2d(:,:),dhi_2d(:,:),dhi0_2d(:,:),dls0_2d(:,:), &
      & velfrict_2d(:,:),eflux0_kinem_2d(:,:),lamw_2d(:,:,:),&
	  & Radbal_2d(:,:),hflux_2d(:,:)
      real(8), allocatable :: snmelt_2d(:,:)
      integer(4), allocatable::fl_sn_2d(:,:),fl_sn_init_2d(:,:), &
      & itop_2d(:,:), init(:,:), num(:)
      integer(4), allocatable :: nstep_2d(:,:)
      real(8), allocatable:: Tsoil4(:), Tsoil3(:), &
      & wl2(:), wl3(:), wl4(:), wl5(:),  &
      & wi2(:), wi3(:), lammoist(:), rosdry(:), filtr(:)
      real(8), allocatable::E2(:),eps2(:),eps3(:), &
      & u2(:),v2(:),w2(:),k1(:),k3(:),k4(:),k5(:),u3(:), &
      & v3(:),C1aup(:),Re(:),row(:),Ri(:),uv(:),E_it1(:), &
      & E_it2(:),C_num(:),E_it3(:),Feps(:),E_it21(:), &
      & eps_it21(:),eps_it1(:),eps_it2(:), &
      & l(:),k2_mid(:),k3_mid(:),k4_mid(:),k5_mid(:), &
      & Um(:),Vm(:),E12(:),eps12(:),knum(:),k2t(:), &
      & Eeps(:),Ecent(:),Epscent(:),res_E(:), &
      & res_eps(:),dresd(:,:,:,:),dres2dE(:),dres2deps(:), &
      & veg_sw(:) 
      real(8), allocatable:: WU_(:),WV_(:),GAMT(:),GAMU(:),GAMV(:), TF(:),KLengu(:)
      real(8), allocatable:: PEMF    (:) , PDENS_DOWN (:), &
      &                       PT_DOWN (:) , PSAL_DOWN  (:), &
      &                       pt_down_f(:)
      real(8), allocatable:: k_turb_T_flux(:), T_massflux(:)

      data nstep /0/

      SAVE

      END MODULE ARRAYS

      

      MODULE TURB_CONST

!     PHYSICAL CONSTANTS
       real(8), parameter:: kar       = 0.4d0
       real(8), parameter:: niu       = 1.007d-6
       real(8), parameter:: roughness = 0.01d0 !water roughness

!     DIMENSIONLESS CONSTANTS IN E-EPS PARAMETERIZATION, according to Goudsmit et al., 2002
       real(8), parameter:: CE0       = 0.09d0 
       real(8), parameter:: CEt0      = 0.072d0
       real(8), parameter:: ceps1     = 1.44d0
       real(8), parameter:: ceps2     = 1.92d0
       real(8), parameter:: sigmaE    = 1.d0
       real(8), parameter:: sigmaeps  = 1.3d0
!      sigmaeps  = 1.111d0, according to Burchard, 2002

       real(8), parameter:: lam_T     = 1.d0 
       real(8), parameter:: lam_gen   = 1.d0 
       real(8), parameter:: lamTM     = 0.14d0

!	 real(8), parameter:: cmu_0     = 0.094d0
!      Other values for cmu_0, proposed in literature, are
!      0.121, 0.077 - see Burchard, 2002

       real(8), parameter:: lam_E0     = CE0/sigmaE   
       real(8), parameter:: lam_eps0   = CE0/sigmaeps 

!      The coefficients in formulas for stability functions (Canuto et al., 2001)
       real(8), parameter:: CEcoef01   = 0.127d0
       real(8), parameter:: CEcoef02   = 0.1526d-1
       real(8), parameter:: CEcoef03   = 0.16d-3

       real(8), parameter:: CEtcoef01  = 0.1190d0
       real(8), parameter:: CEtcoef02  = 0.4294d-2
       real(8), parameter:: CEtcoef03  = 0.66d-3

       real(8), parameter:: CEcoef11   = 1.d0
       real(8), parameter:: CEcoef12   = 0.2d0
       real(8), parameter:: CEcoef13   = 0.315d-1
       real(8), parameter:: CEcoef14   = 0.58d-2
       real(8), parameter:: CEcoef15   = 0.4d-2
       real(8), parameter:: CEcoef16   = 0.4d-4
 
!      The parameters of Mellor and Yamada model (1982)	 

       real(8), parameter:: A1  =  0.92d0
       real(8), parameter:: B1  =  1.66d+1
       real(8), parameter:: C1  =  0.8d-1
       real(8), parameter:: A2  =  0.74d0
       real(8), parameter:: B2  =  1.01d+1
       real(8), parameter:: CL_K_KL_MODEL = 2.d0**(1.5d0)/B1

!     Co is according to Satyanarayana et al., 1999
       real(8), parameter:: Co        = 1.9d0

!     CONSTANTS IN CONVECTION PARAMETERIZATION
       real(8), parameter:: KC0       = 0.d0 !1.d+5 !500.
       real(8), parameter:: dtdz0     = 2.d0

!     CONSTANTS OF SIMOES'S PARAMETERIZATION
       real(8), parameter:: Cs        = 1.5d0  !0.15
       real(8), parameter:: Cs1       = 100.d0    
      
!     lo=0.4*ddz*h1*3.                                                         VS,06.2007
       real(8), parameter:: L0        = 0.1d0 !0.1d0

       real(8), parameter:: CON0      = 0.05d0
       real(8), parameter:: CON1      = 0.11d0
       real(8), parameter:: CON2      = 1.56d0
       real(8), parameter:: CONUV     = 1.d0

      END MODULE TURB_CONST
      
      
      MODULE ASSIM_VAR
      use DRIVING_PARAMS
      integer n_modvar, n_obsvar, size_mod, size_obs, ntim
      
      real(8), allocatable:: y_obs(:,:), x_mod(:,:), b_mod(:,:), &
      & r_obs(:,:), h_obs(:,:), covar(:,:,:,:), c_diff(:,:,:), &
      & RpHBHT(:,:), dep_obs(:,:),y_obs_mean(:),x_mod_mean(:), &
      & delta(:)

      SAVE
     
      contains
      subroutine alloc_assim_var
      n_modvar = 1
      n_obsvar = 1
      size_mod = (M+1)*n_modvar
      size_obs = (M+1)*n_obsvar
      
      allocate (y_obs(1:ntim, 1:size_obs), &
      &          x_mod(1:ntim, 1:size_mod), &
      &          dep_obs(1:ntim, 1:size_obs))
      allocate (r_obs(1:size_obs, 1:size_obs), &
      &          RpHBHT(1:size_obs,1:size_obs), &
      &          b_mod(1:size_mod, 1:size_mod), &
      &          h_obs(1:size_obs, 1:size_mod))
      allocate (c_diff(1:M+1,1:n_obsvar,1:n_obsvar))
      allocate (covar(1:M+1, 1:M+1, 1:n_obsvar, 1:n_obsvar))
      allocate (y_obs_mean(1:size_obs), &
      &          x_mod_mean(1:size_mod), &
      &          delta     (1:size_obs))
      
      y_obs=0.
      x_mod=0.
      r_obs=0.
      b_mod=0. 
      h_obs=0.
      c_diff=0.
      covar=0.
      end subroutine alloc_assim_var
      
      END MODULE assim_var
	  
	  