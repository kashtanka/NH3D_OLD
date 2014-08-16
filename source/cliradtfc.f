 
      
      subroutine opticparset(m,np,albedo,emis,z11)
!-----Updates optic parameters, required for shortwave and longwave radiation 
!-----parametrizations;
!-----They are used as input to subroutines SORAD and IRRAD, called from THERMO      
      use optic_parameters
      
!-----The piece of the code, which calculates the parameters for shortwave code
!-----is adopted from:

C
C********************************************************************
C             re-written by Tarasova October 2005: 
C
C           New aerosol models with 100 levels CO2=360ppm
C
C-----The k-distributions of Fomin and Correa, (2005)
C-----are used 
C
C------------Aerosol properties of Continental, 
C              Maritime and Stratospheric aerosols
C                              are included
C      8 spectral intervals in 0.2-10 mcm:           
C UV+VIS: 
C  nband=3:  0.2-0.303, 0.303-0.323, 0.323-0.7  mcm
C NIR:
C  nband=5:  0.323-1.22, 0.70-1.22, 1.22-10.0, 1.22-2.27,  2.27-10.0
c
c
c-------  reading standard atm of 32 layers MLS.d 
c     
c*********************************************************************
c  Chou&Suarez   July, 1999 
c This is the source program for computing solar fluxes due to
c  absorption by water vapor, ozone, co2, o2, clouds, and aerosols
c  and due to scattering by clouds, aerosols, and gases.
c
c This is a vectorized code.  It computes fluxes simultaneously for
c  m soundings.
c
c The meaning, units and dimension of the input and output parameters
c  are given in the subroutine sorad.
c
c*********************************************************************

!-----The piece of the code, which calculates the parameters for longwave code
!-----is adopted from:

*************************  November 1, 2000  ***************************
c         ------- climate/user/ftp/pub/chou/radiation/ir -------
ctar
ctar  -----------  reading of TRA.d and MLS.d files for the comparison
ctar  -----------  with ICRCCM test cases
ctar
ctar
ctar

      implicit none

      integer m,np,i,j,k,ib,ict,icb
      integer NM           ! NEW
     
C
C---------NEW ARRAYS
C
ctar       REAL CO2(np+1),O2(np+1),co2a(m,np),o2a(m,np)       
ctar       REAL co2layer(np),o2layer(np)

       REAL CODW(np),CODI(np),CCF(np),CWCW(np)
       REAL CWCI(np),REFW(np),REFI(np)
       REAL CODW_m,CODI_m,CCF_m,CWCW_m
       REAL CWCI_m,REFW_m,REFI_m,TAUA_SUM
ctar       
       REAL TAU_S(8),SSA_S(8),ASYM_S(8)
       REAL TAU_M(8),SSA_M(8),ASYM_M(8)       
       REAL TAU_C(8),SSA_C(8),ASYM_C(8)       
       REAL TAUA_V(m,np)
       REAL SAUVB,SAUVD,SAIRB,SAIRD,CSZA 
       REAL FDS_C,FDS_A,FUS_A,FNS_A
       REAL ABS_C,ABS_A, FUT_C
       REAL, allocatable, dimension(:):: 
     &  Z,P,T,WV,OZ,OZm,olayer
      
c-----input parameters
      
      logical high,vege,trace,overcast_lw,cldwater_lw,aerosol

      character*6 hi,tr,cw,ov,ar

      CHARACTER C1*80, cMOD*13, aMOD*6, Catm*5      

c-----output parameters
     
      real(8) albedo(m), emis(m)
      real z11(m,np+1),z11_centr(m,np)
      logical firstcall,overcast,cldwater
      data firstcall /.true./
      
C-----SPECIFY AEROSOL PROPERTIES OF CONTINENTAL (C), MARITIME (M) AND
C-----STRATOSPHERIC (S) AEROSOLS 
 
C              CONTINENTAL (C) AEROSOLS,   8 intervals:
C
C     0.2-0.303 (0.252), 0.303-0.323 (0.313), 0.323-0.7 (0.512), 0.323-1.22, 
C     0.7-1.22,  1.22-10.0, 1.22-2.27, 2.27-10.0
      
      DATA TAU_C/2.09880,1.76676,1.08667,0.87224,0.56770,
     * 0.25307, 0.28098,0.12890/     
     
      DATA SSA_C/0.77936,0.89032,0.89711,0.86887,0.83504,
     * 0.75306,0.76509,0.69953/ 
	
      DATA ASYM_C/0.68684,0.65619,0.63822,0.63599,0.63209,
     * 0.68143,0.65978,0.77775/	    
C
C              MARITIME    (M)  AEROSOLS,   8 intervals
C
      DATA TAU_M/1.23760,1.16595,1.01222,0.97019,0.90512,
     * 0.76353,0.80001,0.60126/  
C
      DATA SSA_M/0.92032,0.97951,0.99000,0.98782,0.98620,
     * 0.96099,0.98319,0.86221/      
C
      DATA ASYM_M/0.75748,0.74395,0.74600,0.74975,0.75541,
     * 0.77743,0.77184,0.80229/

C              STRATOSPHERIC (S)  AEROSOLS,  8 intervals
C
      DATA TAU_S/1.55000,1.53946,1.09667,0.81293,0.45271,
     * 0.10020,0.11383,0.03955/    
C
      DATA SSA_S/5*0.9999,0.90023,0.99161,0.49372/
c       
      DATA ASYM_S/0.69380,0.71722,0.73244,0.68702,0.63582,
     * 0.39456,0.44335,0.17748/
      
      SAVE
    
      if1: IF (firstcall) THEN
             
       OPEN (unit=2001,file='optics/clrdtfc.in',status='old', 
     &  readonly)
       OPEN (unit=2009,file='optics/clrdtfc.out', status='unknown')    
       
       OPEN (unit=201, file='optics/mainlw.in', status='old', readonly)
       OPEN (unit=209, file='optics/mainlw.out', status='unknown')   
c
       READ (201,*) C1,       
!     * ts(1),tb(1),
     * C1,
     * high,    
     * C1,    
     * trace,         
     * C1,    
     * aerosol,
     * C1,    
     * vege,    
     * C1,    
     * overcast_lw,    
     * C1,        
     * cldwater_lw,
     * C1,    
     * co2,n2o,ch4,cfc11,cfc12,cfc22,    
     *  C1,    
!     * fs(1,1),fs(1,2),    
     *  C1,    
     * eg(1,1,1),ev(1,1,1),rvir(1,1,1),
     *  C1,
     *  Catm, 
     *  C1,
!     *  (taucl(1,k,2),k=1,np),   
     *  C1     
!     *  ,(fcld(1,k),k=1,np)                     

       CLOSE (201)
    
       READ (2001,*) C1,cMOD,
     * C1,      
     * C1, 
     * overcast,          
     * C1,
     * C1,     
     * cldwater,     
     * C1,
     * C1,
     * ict,icb,
     * C1,     
     * CODW_m,
     * C1,     
     * CODI_m,
     * C1,     
     * CWCW_m,     
     * C1,     
     * CWCI_m,
     * C1,     
     * REFW_m,     
     * C1,     
     * REFI_m,                        
     * C1,
     * CCF_m,     
     * C1,
     * SAUVB, SAUVD, SAIRB, SAIRD,
     * C1,
     * CSZA,
     * C1,
     * aMOD,
     * C1,
     * crel,
     * C1,   
     * co2,
     * C1,
     * s0  
                        
C  co2 = 300 ppv

       CLOSE (2001)
C
C------READING OF STANDARD ATMOSPHERIC MODEL:
C        MLS.d, TRA.d, SWA.d

       OPEN (unit=2002, file=cMOD, status='old', readonly)

       READ (2002,*) C1, NM
       READ (2002,*) C1

       allocate (Z(NM),P(NM),T(NM),WV(NM),OZ(NM),OZm(NM),
     &  olayer(NM-1))
     
       DO I=1,NM
        READ (2002,55) Z(I),P(I),T(I),WV(I),OZ(I)
       ENDDO

       CLOSE (2002)

       DO I=1,NM
        OZm(I)=OZ(I)*T(I)/(348.37*P(I))  ! kg/kg
       ENDDO        
        
       DO i=1,NM-1 
        olayer(i)=(OZm(NM+1-i)+OZm(NM-i))/2.
       ENDDO        

                     
!       WRITE (2009,*) 'AEROSOL MODEL=', aMOD
!       WRITE (2009,*) 'TAUA_SUM=',  TAUA_SUM
!       WRITE (2009,*) ' AEROSOL OPTICAL THICKNESS IN LAYERS'
!       DO I=1,np
!        WRITE (2009,*) TAUA_V(I)
!       ENDDO
!       WRITE (2009,*) 'ATMOSPHERIC MODEL=', cMOD 
                 
       ENDIF if1
       
      CODW = CODW_m
      CODI = CODI_m
      CCF  = CCF_m
      CWCW = CWCW_m
      CWCI = CWCI_m
      REFW = REFW_m
      REFI = REFI_m

!      Transforming units of z11 from m to km       
           
       do i = 1, m
        do j = 1, np
         z11_centr(i,j) = 0.5*(z11(i,j)+z11(i,j+1))
        enddo 
       enddo 
                     
       do i = 1, m
        do k = 1, np
         do j = 1, NM-1
          if (z(j+1)>=z11_centr(i,k).and.z(j)<=z11_centr(i,k)) then
           oa(i,k)=olayer(NM-j)
          endif 
         enddo 
        enddo
       enddo
          
       
c-----specify aerosol optical thickness (taual), single-scattering
c     albedo (ssaal), and asymmetry factor (asyal)

ctar------------ FORMATION OF AEROSOL VERTICAL PROFILE IN STRATOSPHERE

!-----------It is assumed, that the atmospheric model domain-----------------
!---------------------lies entirely in troposphere---------------------------

!      DO i=1, m
!       DO j=1, np-12 ! 12 is value, dependent on the height - to be corrected! 
!        IF (z11_centr(i,j) .GT. 12.0 .AND. z11_centr(i,j) .LE. 20.0)THEN
!         TAUA_V(i,j)=0.000218*(z11(i,j)-z11(i,j+1))
!        ENDIF
!        IF (z11_centr(i,j).GT.20.0) TAUA_V(i,j)=0.
!       ENDDO
!      ENDDO

!      IF (aMOD.EQ.'A-ZERO') TAUA_V = 0.
      
!      do ib = 1, 8
!       do k = 1, np-12
!         do i = 1, m
!          taual(i,k,ib) = TAU_S(ib)*TAUA_V(i,k)        
!          ssaal(i,k,ib) = SSA_S(ib)
!          asyal(i,k,ib) = ASYM_S(ib)
!         enddo
!        enddo
!      enddo

!-----------It is assumed, that the atmospheric model domain-----------------
!---------------------lies entirely in troposphere---------------------------

C---------IF CONT-I aerosol profile in troposphere

!--------In the case, when the upper part of the domain lies in stratosphere---
!--------the subsequent code should be rewritten in following manner------------

!stat      IF (aMOD.EQ.'CONT-I') THEN
!stat       DO i = 1, m
!stat        DO j = np-11, np ! 11 is value, dependent on the height - to be corrected! 
!stat         IF (z11_centr(i,j).GT.0.0.AND.z11_centr(i,j).LE.2.0) THEN
!stat          TAUA_V(i,j)=0.1*(z11(i,j)-z11(i,j+1))
!stat         ENDIF
!stat         IF (z11_centr(i,j).GT.2.0.AND.z11_centr(i,j).LE.12.0) THEN
!stat          TAUA_V(i,j)=0.0025*(z11(i,j)-z11(i,j+1))
!stat         ENDIF
!stat        ENDDO
!stat       ENDDO
       
       IF (aMOD.EQ.'CONT-I') THEN
       DO i = 1, m
        DO j = 1, np 
         IF (z11_centr(i,j).GT.0.0.AND.z11_centr(i,j).LE.2.0) THEN
          TAUA_V(i,j)=0.1*(z11(i,j)-z11(i,j+1))
         ENDIF
         IF (z11_centr(i,j).GT.2.0.AND.z11_centr(i,j).LE.12.0) THEN
          TAUA_V(i,j)=0.0025*(z11(i,j)-z11(i,j+1))
         ENDIF
        ENDDO
       ENDDO
       
       do ib = 1,8
        do k = 1, np
          do i = 1, m
           taual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)        
           ssaal(i,k,ib) = SSA_C(ib)
           asyal(i,k,ib) = ASYM_C(ib)
          enddo
         enddo
       enddo
      ENDIF

C---------IF CONT-A aerosol profile in troposphere
      IF (aMOD.EQ.'CONT-A') THEN
       DO i = 1, m
        DO j = 1, np ! 11 is value, dependent on the height - to be corrected! 
         IF (z11_centr(i,j).GT.0.0.AND.z11_centr(i,j).LE.2.0) THEN
          TAUA_V(i,j)=0.5*(z11(i,j)-z11(i,j+1))
         ENDIF
         IF (z11_centr(i,j).GT.2.0.AND.z11_centr(i,j).LE.12.0) THEN
          TAUA_V(i,j)=0.0025*(z11(i,j)-z11(i,j+1))
         ENDIF
        ENDDO
       ENDDO
       
       do ib=1,8
        do k = 1, np
          do i = 1, m
           taual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)        
           ssaal(i,k,ib) = SSA_C(ib)
           asyal(i,k,ib) = ASYM_C(ib)
          enddo
         enddo
       enddo
      ENDIF

C---------IF MAR-I aerosol profile in troposphere
      IF (aMOD.EQ.'MAR-I') THEN
       DO i=1, m
        DO j=1,np
         IF (z11_centr(i,j).GT.0.0.AND.z11_centr(i,j).LE.2.0) THEN
          TAUA_V(i,j)=0.025*(z11(i,j)-z11(i,j+1))
         ENDIF
         IF (z11_centr(i,j).GT.2.0.AND.z11_centr(i,j).LE.12.0) THEN
          TAUA_V(i,j)=0.0025*(z11(i,j)-z11(i,j+1))
         ENDIF
        ENDDO
       ENDDO

       do ib=1,8
        do k= 1, np-2
          do i= 1, m
           taual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)        
           ssaal(i,k,ib) = SSA_C(ib)
           asyal(i,k,ib) = ASYM_C(ib)
          enddo
         enddo
        enddo

        do ib=1,8
         do k= np-1, np
           do i= 1, m
            taual(i,k,ib) = TAU_M(ib)*TAUA_V(i,k)        
            ssaal(i,k,ib) = SSA_M(ib)
            asyal(i,k,ib) = ASYM_M(ib)
           enddo
          enddo
        enddo
      ENDIF

C--------IF MAR-II aerosol profile in troposphere
       IF (aMOD.EQ.'MAR-II') THEN
        DO i=1, m
         DO j=1,np
          IF (z11_centr(i,j).GT.0.0.AND.z11_centr(i,j).LE.2.0) THEN
           TAUA_V(i,j)=0.025*(z11(i,j)-z11(i,j+1))
          ENDIF
          IF (z11_centr(i,j).GT.2.0.AND.z11_centr(i,j).LE.6.0) THEN
           TAUA_V(i,j)=0.75*(z11(i,j)-z11(i,j+1))
          ENDIF
          IF (z11_centr(i,j).GT.6.0.AND.z11_centr(i,j).LE.12.0) THEN
           TAUA_V(i,j)=0.0025*(z11(i,j)-z11(i,j+1))
          ENDIF
         ENDDO
        ENDDO

        do ib=1,8
         do k= 1, np-2
           do i= 1, m
            taual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)        
            ssaal(i,k,ib) = SSA_C(ib)
            asyal(i,k,ib) = ASYM_C(ib)
           enddo
          enddo
        enddo

        do ib=1,8
         do k= np-1, np
           do i= 1, m
            taual(i,k,ib) = TAU_M(ib)*TAUA_V(i,k)        
            ssaal(i,k,ib) = SSA_M(ib)
            asyal(i,k,ib) = ASYM_M(ib)
           enddo
          enddo
        enddo
       ENDIF
C
       TAUA_SUM=0.
       DO i=1, m
        DO j=1,np
        TAUA_SUM=TAUA_SUM+TAUA_V(i,j)
        ENDDO
       ENDDO

c-----pl,ta,wa and oa are, respectively, the level pressure (mb), layer
c     temperature (k), layer water vapor amount (g/g), and layer ozone 
c     concentration g/g

c-----specify cloud optical thickness (taucld), amount (fcld), effective
c     particle size (reff).   cloud ice (index 1), liquid (index 2), and
c     rain (index 3) are allowed to co-exit in a layer.
c     cwc is the cloud ice/water concentration. If cldwater=.true.,
c     taucld is computed from cwc and reff.  If cldwater=.false.,
c     taucld is an input parameter, and cwc is irrelevant

       do k= 1, np
        do i= 1, m
         taucld(i,k,1)= CODI(k)  ! ICE PARTICLES
         taucld(i,k,2)= CODW(k)  ! WATER PARTICLES
         taucld(i,k,3)= 0.0      ! RAIN DROPS
         
!--------reff should be calculated in cloud micriphysics model!---------
         reff(i,k,1)  = REFI(k)  ! ICE PARTICLES
         reff(i,k,2)  = REFW(k)  ! WATER PARTICLES
         reff(i,k,3)  = 0.       ! ICE PARTICLES         
!            fcld(i,k)    = CCF(k)  ! CLOUD COVER 
!            cwc(i,k,1)   = CWCI(k)
!            cwc(i,k,2)   = CWCW(k)
!            cwc(i,k,3)   = 0.0
        enddo
       enddo
      

!-----specify surface albedos in solar ir and uv+par regions are assumed to be equal
!-----to integral albedo
       
         rsuvbm = albedo
         rsuvdf = albedo
         rsirbm = albedo
         rsirdf = albedo
         
!------Specification of longwave optic parameters-----------

 
!       if ( high        ) hi='true'
!       if ( .not. high  ) hi='false'
!       if ( trace       ) tr='true'
!       if ( .not. trace ) tr='false'
!       if ( overcast    ) ov='true'
!       if (.not.overcast) ov='false'
!       if ( cldwater    ) cw='true'
!       if (.not.cldwater) cw='false'
!       if ( aerosol     ) ar='true'
!       if (.not.aerosol)  ar='false'
 
c-----m x n is the number of atmospheres
 
c-----np is the number of atmospheric layers
 
c-----np+1 is the index for the surface level
 
c-----specify surface properties
 
        do ib=1,10
         do i=1,m
 
ctar           eg(i,1,ib)=1.0
ctar           ev(i,1,ib)=1.0
ctar           rv(i,1,ib)=0.0
	   
           eg(i,1,ib)=emis(i)
           ev(i,1,ib)=emis(i)
           rvir(i,1,ib)=rvir(1,1,1)
	   
ctar           eg(i,2,ib)=1.0
ctar           ev(i,2,ib)=1.0
ctar           rv(i,2,ib)=0.0
	   
           eg(i,2,ib)=emis(i)
           ev(i,2,ib)=emis(i)
           rvir(i,2,ib)=rvir(1,1,1)	   
 
         enddo
        enddo
 
c-----assign co2 amount. units are parts/part by volumn
 
ctar       co2=350.e-6
ctar       n2o=0.28e-6
ctar       ch4=1.75e-6
ctar       cfc11=0.3e-9
ctar       cfc12=0.5e-9
ctar       cfc22=0.2e-9
 
c-----specify cloud optical thickness (taucl), amount (fcld).  
c     cloud ice (index 1), water (index 2),and rain (index 3) are allowed to co-exit 
c     in a layer. cwc is the cloud ice/water concentration. If cldwater=.true.,
c     taucl will be computed from cwc.  If cldwater=.false.,
c     taucl is an input parameter, and cwc is irrelevent.
 
       do k = 1, np
        do i = 1, m
          taucl_lw(i,k,1) = 0.0
c          taucl(i,k,2) = 0.0
          taucl_lw(i,k,3) = 0.0
c          fcld(i,k)    = 0.0
!          cwc(i,k,1)   = 0.0
!          cwc(i,k,2)   = 0.0
!          cwc(i,k,3)   = 0.0
        enddo
       enddo
 
c----aerosols properties
 
      do j=1,na
       do ib=1,10
        do k=1,np
         do i=1,m
          taual_lw(i,k,ib,j)=0.00
          ssaal_lw(i,k,ib,j)=0.99
          asyal_lw(i,k,ib,j)=0.75
         enddo
        enddo
       enddo
      enddo
 
c----aerosol optical thickness
 
       do ib=1,10
        do k=np-5,np
         do i=1,m
          taual_lw(i,k,ib,1)=0.015/3.0
          taual_lw(i,k,ib,2)=0.015/3.0
          taual_lw(i,k,ib,3)=0.015/3.0
 
ctar          taual(i,k,ib,1)=0.00
ctar          taual(i,k,ib,2)=0.00
ctar          taual(i,k,ib,3)=0.00
 
          ssaal_lw(i,k,ib,1)=0.99
          ssaal_lw(i,k,ib,2)=0.99
          ssaal_lw(i,k,ib,3)=0.99
 
          asyal_lw(i,k,ib,1)=0.75
          asyal_lw(i,k,ib,2)=0.75
          asyal_lw(i,k,ib,3)=0.75
         enddo
        enddo
       enddo 
 
 
      if (firstcall) firstcall = .false.

c----- 
  100 format (/,'  ****  Sample mid-latitude summer atmosphere ****',/)
  101 format (/,'  s0    =', f7.1,'  w/m^2',/,'  cosz  =', f6.3,/,
     * '  ict   =',i3,/,'  icb   =',i3,/,
     * '  rsuvbm=',f5.2,/,'  rsuvdf=',f5.2,/,'  rsirbm=',f5.2,/,
     * '  rsirdf=',f5.2,/)
  102 format (/,' ******  INPUT DATA  ******',/)
  103 format (10x,'p',9x,'z',8x,'delz',5x,'T',8x,'q',
     * 8x,'o3',6x,'reff(1,2)',
     *        3x,'fcld',2x,'taual')
  104 format (9x,'(mb)',6x,'(km)',6x,'(km)',4x,'(K)',
     *5x,'(molec/cm2km)',5x, '(micron)',/)
  105 format (i4,f10.4,f10.4)
ctar
  305 format (i5,3E12.4)
  306 format (i4,3F10.2)    
  106 format (35x,f7.2,1p,2e10.3,1x,0p,1x,2f5.1,2x,f5.1,2x,f5.3)

  202 format (//,' ******  RESULTS  ******',/)
  203 format (10x,'p',8x,'flx',8x,'flc',5x,'all-sky heating')
  204 format (9x,'(mb)',4x,'(W/m^2)',4x,'(W/m^2)',7x,'(C/day)',/)
  205 format (i4,f10.4,f10.2,1x,f10.2)
  206 format (37x,f10.3)
  207 format (//,' fdiruv  = ',f7.2,/,' fdifuv  = ',f7.2,/,
     *           ' fdirpar = ',f7.2,/,' fdifpar = ',f7.2,/,
     *           ' fdirir  = ',f7.2,/,' fdifir  = ',f7.2)
  208 format (1x,'H,km',6x,'p(mb)',3x,'F all,d ',6x,
     * 'F all,u ',6x,'HR(K/d)')
ctar     * ,'F clr,d',2x,'F clr,u (W/m^2)',/)
  209 format (i5,f12.5,f12.6,1x,f12.6,1x,f10.2,1x,f10.2)
  210 format (2f10.4)   
ctar  55       FORMAT (F5.1,E12.4,F8.2,4E12.4)  
  55  FORMAT (1X,5E12.4) 
            
      end subroutine opticparset


c***********************       CLIRAD-SW      ************************
ctar 
      subroutine sorad (m,np,crel,pl,ta,wa,oa,co2,zdel,
     *                  overcast,cldwater,cwc,taucld,reff,fcld,ict,icb,
     *                  taual,ssaal,asyal,
     *                  cosz,rsuvbm,rsuvdf,rsirbm,rsirdf,
     *                  flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,
     *                  fdirir,fdifir,
     *                  flx_d,flx_u,flc_d,flc_u) 
c
c**********************************************************************
c 
c     The subroutine is re-written by Tarasova for using of
c              Fomin and Correa, 2005  parameterizations, August 2005
c   
c    wa,oa,co2- water vapor g/g, ozone g/g, co2=300 ppv 
c
c    swa(m,np),soa(m,np),sco2(m,np),so2(m,np) - H2O,O3,CO2,O2 concentration
c                          in molec/cm^2 km divided here by  e+21,e+17,e+19,e+22
c    wh(m,np),oh(m,np),co2h(m,np),o2h(m,np) - gaseous amount at the solar 
c    beam in molec/cm^2 divided by e+21,e+17,e+19,e+22
c    zdel((m,np) - layer height in km
c
c
c
c Following the NASA Technical Memorandum (NASA/TM-1999-104606, Vol. 15)
c  of Chou and Suarez (1999), this routine computes solar fluxes due to 
c  absorption by water vapor, ozone, co2, o2, clouds, and aerosols and 
c  due to scattering by clouds, aerosols, and gases.
c
c This code computes fluxes simultaneously for m soundings.
c
c Cloud ice, liquid, and rain particles are allowed to co-exist in a layer. 
c
c There is an option of providing either cloud ice/water mixing ratio 
c  (cwc) or optical thickness (taucld).  If the former is provided, set
c  cldwater=.true., and taucld is computed from cwc and reff as a
c  function of spectra band. Otherwise, set cldwater=.false., and
c  specify taucld, independent of spectral band.
c
c If no information is available for the effective particle size, reff, 
c  default values of 10 micron for liquid water and 75 micron for ice may be 
c  used. 
c  The size of raindrops, reff(3), is irrelevant in this code. It can be
c  set to any values.
c  For a clear layer, reff can be set to any values except zero.
c
c The maximum-random assumption is applied for treating cloud
c  overlapping. Clouds are grouped into high, middle, and low clouds
c  separated by the level indices ict and icb.  For detail, see
c  subroutine "cldscale".
c
c In a high spatial-resolution atmospheric model, fractional cloud cover
c  might be computed to be either 0 or 1.  In such a case, scaling of the
c  cloud optical thickness is not necessary, and the computation can be
c  made faster by setting overcast=.true.  Otherwise, set the option
c  overcast=.false.
c
c Aerosol optical thickness, single-scattering albedo, and asymmetry
c  factor can be specified as functions of height and spectral band.
c
c----- Input parameters:                           
c                                                   units      size
c  number of soundings (m)                          n/d         1
c  number of atmospheric layers (np)                n/d         1
c  level pressure (pl)                              mb        m*(np+1)
c  layer temperature (ta)                           k         m*np
c  layer specific humidity (wa)                     gm/gm     m*np
c  layer ozone concentration (oa)                   gm/gm     m*np
c  co2 mixing ratio by volume (co2)                 pppv        1
c  option for scaling cloud optical thickness       n/d         1
c        overcast="true" if scaling is NOT required
c        overcast="false" if scaling is required
c  option for cloud optical thickness               n/d         1
c        cldwater="true" if cwc is provided
c        cldwater="false" if taucld is provided
c  cloud water mixing ratio (cwc)                  gm/gm      m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d       m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  effective cloud-particle size (reff)          micrometer   m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud amount (fcld)                            fraction    m*np
c  level index separating high and middle           n/d         1
c        clouds (ict)
c  level index separating middle and low            n/d         1
c        clouds (icb)
c  aerosol optical thickness (taual)                n/d       m*np*11
c  aerosol single-scattering albedo (ssaal)         n/d       m*np*11
c  aerosol asymmetry factor (asyal)                 n/d       m*np*11
c        in the uv region :
c           index  1 for the 0.175-0.225 micron band
c           index  2 for the 0.225-0.245; 0.260-0.280 micron band
c           index  3 for the 0.245-0.260 micron band
c           index  4 for the 0.280-0.295 micron band
c           index  5 for the 0.295-0.310 micron band
c           index  6 for the 0.310-0.320 micron band
c           index  7 for the 0.325-0.400 micron band
c        in the par region :
c           index  8 for the 0.400-0.700 micron band
c        in the infrared region :
c           index  9 for the 0.700-1.220 micron band
c           index 10 for the 1.220-2.270 micron band
c           index 11 for the 2.270-10.00 micron band
c   cosine of solar zenith angle (cosz)              n/d      m
c   uv+visible sfc albedo for beam radiation
c        for wavelengths<0.7 micron (rsuvbm)       fraction   m
c   uv+visible sfc albedo for diffuse radiation
c        for wavelengths<0.7 micron (rsuvdf)       fraction   m
c   ir sfc albedo for beam radiation
c        for wavelengths>0.7 micron  (rsirbm)      fraction   m
c   ir sfc albedo for diffuse radiation (rsirdf)   fraction   m
c
c----- Output parameters
c
c   all-sky flux (downward minus upward) (flx)     fraction   m*(np+1)
c   clear-sky flux (downward minus upward) (flc)   fraction   m*(np+1)
c   all-sky direct downward uv (0.175-0.4 micron)
c                flux at the surface (fdiruv)      fraction   m
c   all-sky diffuse downward uv flux at
c                the surface (fdifuv)              fraction   m
c   all-sky direct downward par (0.4-0.7 micron)
c                flux at the surface (fdirpar)     fraction   m
c   all-sky diffuse downward par flux at
c                the surface (fdifpar)             fraction   m
c   all-sky direct downward ir (0.7-10 micron)
c                flux at the surface (fdirir)      fraction   m
c   all-sky diffuse downward ir flux at
c                the surface (fdifir)              fraction   m
c
c----- Notes:
c
c    (1) The unit of output fluxes (flx,flc,etc.) is fraction of the
c        insolation at the top of the atmosphere.  Therefore, fluxes
c        are the output fluxes multiplied by the extra-terrestrial solar
c        flux and the cosine of the solar zenith angle.
c    (2) pl( ,1) is the pressure at the top of the model, and
c        pl( ,np+1) is the surface pressure.
c    (3) the pressure levels ict and icb correspond approximately
c        to 400 and 700 mb.
c
c-----If coding errors are found, please notify Ming-Dah Chou at
c     chou@climate.gsfc.nasa.gov
c        
c*************************************************************************
 
      implicit none

c-----input parameters

      integer m,np,ict,icb
ctar      
      real pl(m,np+1),ta(m,np)
ctar      
      real so2(m,np),o2h(m,np),sco2(m,np),co2h(m,np)
      real swa(m,np), wh(m,np), soa(m,np), oh(m,np)
ctar      
      real wa(m,np),oa(m,np),co2            
ctar
      real zdel(m,np)      
      real cwc(m,np,3),taucld(m,np,3),reff(m,np,3),fcld(m,np)
      real taual(m,np,8),ssaal(m,np,8),asyal(m,np,8)
      real cosz(m),rsuvbm(m),rsuvdf(m),rsirbm(m),rsirdf(m)
      logical overcast,cldwater

c-----output parameters

      real flx(m,np+1),flc(m,np+1)
      REAL flx_d(m,np+1),flx_u(m,np+1)   ! NEW
      REAL flc_d(m,np+1),flc_u(m,np+1)   ! NEW      
C      
      real fdiruv (m),fdifuv (m)
      real fdirpar(m),fdifpar(m)
      real fdirir (m),fdifir (m)

c-----temporary array
 
      integer i,j,k,ntop,nctop(m)
      real cwp(m,np,3)
ctar      
      real dp(m,np)      
      real x
      real snt(m) 
      character crel*1     
c
      do i=1,m 
c-----snt is the secant of the solar zenith angle
         snt(i)=1.0/cosz(i)
      enddo
ctar
ctar
c----------wa, oa are in g/g, co2 is in ppv, o2 is 23.14%
c    Transformation of all gaseos concentrations into molec/cm^2 km
c             swa,    soa,    sco2,    so2   (molec/cm^2 km)
c  divided by e+21,   e+17,   e+19,    e+22
      
      do k=1,np
      do i=1,m
ctar
      swa(i,k)=wa(i,k)*3.346*348.37*0.5*
     * (pl(i,k)+pl(i,k+1))/ta(i,k)
ctar      
      soa(i,k)=oa(i,k)*1.255e+4*348.37*0.5*
     * (pl(i,k)+pl(i,k+1))/ta(i,k)
      so2(i,k)=0.007625e+3*(pl(i,k)+pl(i,k+1))/ta(i,k)      
      sco2(i,k)=0.0364*co2*(pl(i,k)+pl(i,k+1))/ta(i,k) 
ctar
      enddo 
      enddo       
      
c
!       write(2009,*) 'Concentration in molec/cm^2 km'
!       write(2009,115) 
!      do k=1,np 
!      write(2009,110) pl(1,k),swa(1,k),soa(1,k),so2(1,k),sco2(1,k)
!      enddo
!      write(2009,110) pl(1,np+1)
 110   format (5E12.4)
 115   format (3x,'P,mb',6x,'H2Oe+21',6x,'O3e+17',6x,
     * 'O2e+22', 6x, 'CO2e+19')
c 
ctar
c----Calculation of wh, oh, co2h, o2h - gaseous amount in molec/cm^2
c------------------at the solar beam    
c 
!       write(2009,*) 'snt=', snt(1)
!       write(2009,*) 'zdel(1,k)'
!       write(2009,*)  (zdel(1,k),k=1,np )     

      do i=1,m
      wh(i,1)=swa(i,1)*zdel(i,1)*snt(i)
      oh(i,1)=soa(i,1)*zdel(i,1)*snt(i)
      o2h(i,1)=so2(i,1)*zdel(i,1)*snt(i)       
      co2h(i,1)=sco2(i,1)*zdel(i,1)*snt(i)                  
      enddo       
c 
      do k=1,np-1
      do i=1,m
      wh(i,k+1)=wh(i,k)+swa(i,k+1)*zdel(i,k+1)*snt(i)
      oh(i,k+1)=oh(i,k)+soa(i,k+1)*zdel(i,k+1)*snt(i)
      o2h(i,k+1)= o2h(i,k)+so2(i,k+1)*zdel(i,k+1)*snt(i)      
      co2h(i,k+1)= co2h(i,k)+sco2(i,k+1)*zdel(i,k+1)*snt(i)                  
      enddo 
      enddo       
c
ctar          
c
!       write(2009,*) 'Amount at the solar beam in molec/cm^2 '
!       write(2009,115) 
!      do k=1,np 
!      write(2009,110) pl(1,k),wh(1,k),oh(1,k),o2h(1,k),co2h(1,k)
!      enddo
!      write(2009,110) pl(1,np+1)
c      
      do k=1,np
       do i=1,m
c
c-----compute layer thickness. indices for the surface level and
c     surface layer are np+1 and np, respectively.
 
          dp(i,k)=pl(i,k+1)-pl(i,k)
c 
c
c-----compute layer cloud water amount (gm/m**2)
c     the index is 1 for ice crystals, 2 for liquid drops, and
c     3 for rain drops

          x=1.02*10000.*dp(i,k)   
          cwp(i,k,1)=x*cwc(i,k,1)
          cwp(i,k,2)=x*cwc(i,k,2)
          cwp(i,k,3)=x*cwc(i,k,3)

       enddo
      enddo
c
c-----initialize fluxes for all-sky (flx), clear-sky (flc), and
c     flux reduction (df)
c
      do k=1,np+1
       do i=1,m
          flx(i,k)=0.
          flx_d(i,k)=0.    ! NEW
          flx_u(i,k)=0.    ! NEW
          flc(i,k)=0.                              
          flc_d(i,k)=0.    ! NEW
          flc_u(i,k)=0.    ! NEW                   
ctar         
       enddo
      enddo

c-----compute solar uv and par fluxes
C
ctar 
!      CONTINUE
      call soluv (m,np,crel,swa,wh,soa,oh,pl,ta,
     *            so2,o2h,dp,zdel,overcast,cldwater,
     *            cwp,taucld,reff,ict,icb,fcld,cosz,
     *            taual,ssaal,asyal,rsuvbm,rsuvdf,
     *            flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,
     *            flx_d,flx_u,flc_d,flc_u)  ! NEW   
     
!      continue
c-----compute and update solar ir fluxes
C
cc    
C
      call solir (m,np,crel,swa,wh,soa,oh,so2,o2h,sco2,co2h,
     * dp,zdel,pl,ta,
     * overcast,cldwater,cwp,taucld,reff,ict,icb,fcld,cosz,
     *            taual,ssaal,asyal,rsirbm,rsirdf,
     *            flx,flc,fdirir,fdifir,
     *            flx_d,flx_u,flc_d,flc_u)  ! NEW      
C
      return
      end  
c
c******************************************************************
ctar
      subroutine soluv (m,np,crel,swa,wh,soa,oh,pl,ta,
     *            so2,o2h,dp,zdel,overcast,cldwater,
     *            cwp,taucld,reff,ict,icb,fcld,cosz,
     *            taual,ssaal,asyal,rsuvbm,rsuvdf,
     *            flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,
     *            flx_d,flx_u,flc_d,flc_u)  ! NEW   

c******************************************************************
c     The subroutine is re-written by Tarasova for using of
c              Fomin and Correa, 2005  parameterizations, August 2005
c    wa(m,np),oa(m,np),co2,o2 - gaseous concentration in g/m^3, g/m^3, ppv, and %
c    swa(m,np),soa(m,np),so2(m,np) - H2O,O3,O2 concentration
c                          in molec/cm^2 km divided by  e+21,e+17,e+22
c    wh(m,np),oh(m,np),o2h(m,np) - gaseous amount at the solar beam in molec/cm^2
c    divided by e+21,e+17,e+22
c    zdel((m,np) - layer height in km
c 
c  nband=3:  0.2-0.303, 0.303-0.323, 0.323-0.7  mcm
c======================================================================
c  Chou and Suarez, 1999
c
c----- Input parameters:                            units      size
c
c  number of soundings (m)                          n/d         1
c  number of atmospheric layers (np)                n/d         1
c  layer scaled-water vapor content (wh)          gm/cm^2      m*np
c  layer ozone content (oh)                      (cm-atm)stp   m*np
c  layer pressure thickness (dp)                    mb         m*np
c  option for scaling cloud optical thickness       n/d         1
c        overcast="true" if scaling is NOT required
c        overcast="false" if scaling is required
c  input option for cloud optical thickness         n/d         1
c        cldwater="true" if taucld is provided
c        cldwater="false" if cwp is provided
c  cloud water amount (cwp)                        gm/m**2     m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d        m*np*3
c       index 1 for ice particles
c       index 2 for liquid drops
c       index 3 for rain drops
c  effective cloud-particle size (reff)          micrometer    m*np*3
c       index 1 for ice paticles
c       index 2 for liquid drops
c       index 3 for rain drops
c  level index separating high and                  n/d        m
c       middle clouds (ict)
c  level indiex separating middle and               n/d        m
c       low clouds (icb)
c  cloud amount (fcld)                            fraction     m*np
c  cosine of solar zenith angle (cosz)              n/d        m
c  aerosol optical thickness (taual)                n/d        m*np*11
c  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
c  aerosol asymmetry factor (asyal)                 n/d        m*np*11
c  uv+par surface albedo for beam                 fraction     m
c       radiation (rsuvbm)
c  uv+par surface albedo for diffuse              fraction     m
c       radiation (rsuvdf)
c
c---- temporary array
c
c  scaled cloud optical thickness                   n/d        m*np
c       for beam radiation (tauclb)
c  scaled cloud optical thickness                   n/d        m*np
c       for diffuse radiation  (tauclf)     
c
c----- output (updated) parameters:
c
c  all-sky net downward flux (flx)               fraction      m*(np+1)
c  clear-sky net downward flux (flc)             fraction      m*(np+1)
c  all-sky direct downward uv flux at
c       the surface (fdiruv)                     fraction      m
c  all-sky diffuse downward uv flux at
c       the surface (fdifuv)                     fraction      m
c  all-sky direct downward par flux at
c       the surface (fdirpar)                    fraction      m
c  all-sky diffuse downward par flux at
c       the surface (fdifpar)                    fraction      m
c
c***********************************************************************

      implicit none

c-----input parameters

      integer m,np,ict,icb
      real taucld(m,np,3),reff(m,np,3),fcld(m,np)
      real cwp(m,np,3),wh(m,np),oh(m,np),dp(m,np)
ctar
      real swa(m,np),soa(m,np),so2(m,np),o2h(m,np),
     *  sco2(m,np),co2h(m,np),zdel(m,np) 
      real pl(m,np+1),ta(m,np),pa(m,np)         
      real taual(m,np,11),ssaal(m,np,11),asyal(m,np,11)
      real rsuvbm(m),rsuvdf(m),cosz(m)
      logical overcast,cldwater

c-----output (updated) parameter

      real flx(m,np+1),flc(m,np+1)
      real fdiruv (m),fdifuv (m)
      real fdirpar(m),fdifpar(m)
      REAL flx_d(m,np+1),flx_u(m,np+1)   !NEW
      REAL flc_d(m,np+1),flc_u(m,np+1)   !NEW     

c-----static parameters

      integer nband
c   
ctar      
      parameter (nband=3)
      real acuv(m,np,nband)
ctar      
      real hk(nband),ry(nband)      
      real aig(3),awg(3),arg(3)
      real aib(2),awb(2),arb(2)

c-----temporary array

      integer i,k,ib
      integer ih1,ih2,im1,im2,is1,is2
      real dsm(m)
      real tauclb(m,np),tauclf(m,np),asycl(m,np)
ctar     
      real taurs(m,np),taugas(m,np)
      real tausto(m,np),ssatau(m,np),asysto(m,np)
      real tautob(m,np),ssatob(m,np),asytob(m,np)
      real tautof(m,np),ssatof(m,np),asytof(m,np)
      real taux,reff1,reff2,g1,g2,g3
      real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2),
     *     rs(m,np+1,2),ts(m,np+1,2)
      real fall(m,np+1),fclr(m,np+1),fsdir(m),fsdif(m)
      real asyclt(m),cc(m,3)
      real rrt(m,np),ttt(m,np),tdt(m,np),rst(m,np),tst(m,np)
      real dum1(m,np+1),dum2(m),dum(m,np)
      REAL fall_d(m,np+1),fall_u(m,np+1)   ! NEW
      REAL fclr_d(m,np+1),fclr_u(m,np+1)   ! NEW 
      character crel*1   
C      
c-----hk is the fractional extra-terrestrial solar flux in each
c     of the 3 bands.

      data hk/.01239,.009713,.44297/     
c
c-----ry is the extinction coefficient for Rayleigh scattering.
c     
c
      data ry /63.59,34.48,6.406/
c
c 
c-----coefficients for computing the extinction coefficients of ice, 
c     water, and rain particles, independent of spectral band. (Table 4)

      data aib/ 3.33e-4,2.52/
      data awb/-6.59e-3,1.65/
      data arb/ 3.07e-3,0.00/

c-----coefficients for computing the asymmetry factor of ice, water,
c     and rain particles, independent of spectral band. (Table 6)

      data aig/.74625,.0010541,-.00000264/
      data awg/.82562,.0052900,-.00014866/
      data arg/.883,0.0,0.0/

c-----initialize fdiruv, fdifuv, surface reflectances and transmittances.
c     the reflectance and transmittance of the clear and cloudy portions
c     of a layer are denoted by 1 and 2, respectively.
c     cc is the maximum cloud cover in each of the high, middle, and low
c     cloud groups.
c     1/dsm=1/cos(53) = 1.66
C
C
            
      do i=1,m                    
         dsm(i)=0.602
         fdiruv(i)=0.0
         fdifuv(i)=0.0
         rr(i,np+1,1)=rsuvbm(i)
         rr(i,np+1,2)=rsuvbm(i)
         rs(i,np+1,1)=rsuvdf(i)
         rs(i,np+1,2)=rsuvdf(i)
         td(i,np+1,1)=0.0
         td(i,np+1,2)=0.0
         tt(i,np+1,1)=0.0
         tt(i,np+1,2)=0.0
         ts(i,np+1,1)=0.0
         ts(i,np+1,2)=0.0
         cc(i,1)=0.0
         cc(i,2)=0.0
         cc(i,3)=0.0
      enddo

ctar
c
ctar-----Compute absorption coefficients of H2O,O3,O2
c
       call kuv (m,np,nband,soa,so2,swa,oh,o2h,
     * wh,acuv)
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.2-0.303, abs: O3,O2 '      
!      do k=1,np
!      write (2009,115) acuv(m,k,1)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.303-0.323, abs: O3 '      
!      do k=1,np
!      write (2009,115) acuv(m,k,2)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.323-0.7, abs: O3,H2O '      
!      do k=1,np
!      write (2009,115) acuv(m,k,3)       
!      enddo
c
c      write (9,*) 'gaseous absorption coefficients, km-1'
c      write (9,*) ' WL: 0.323-1.22, abs: O2 '      
c      do k=1,np
c      write (9,115) acuv(m,k,4)       
c      enddo
c
c      
 115   format (E12.4)      
c
c             
c-----Compute cloud optical thickness.  Eqs. (4.6) and (4.11)

      if (cldwater) then

       do k=1,np
        do i=1,m
          taucld(i,k,1)=cwp(i,k,1)*(aib(1)+aib(2)/reff(i,k,1))
          taucld(i,k,2)=cwp(i,k,2)*(awb(1)+awb(2)/reff(i,k,2))
          taucld(i,k,3)=cwp(i,k,3)* arb(1)
        enddo
       enddo

      endif

c-----options for scaling cloud optical thickness

      if (overcast) then

       do k=1,np
        do i=1,m
          tauclb(i,k)=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          tauclf(i,k)=tauclb(i,k)
        enddo
       enddo

       do k=1,3
        do i=1,m
           cc(i,k)=1.0
        enddo
       enddo

      else

c-----scale cloud optical thickness in each layer from taucld (with
c     cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
c     tauclb is the scaled optical thickness for beam radiation and
c     tauclf is for diffuse radiation (see section 7).
C
C
         call cldscale(m,np,cosz,fcld,taucld,ict,icb,
     *                 cc,tauclb,tauclf)

      endif
C       
C
c-----cloud asymmetry factor for a mixture of liquid and ice particles.
c     unit of reff is micrometers. Eqs. (4.8) and (6.4)

      do k=1,np
       do i=1,m

           asyclt(i)=1.0
           taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)

         if (taux.gt.0.02 .and. fcld(i,k).gt.0.01) then

           reff1=min(reff(i,k,1),130.)
           reff2=min(reff(i,k,2),20.0)

           g1=(aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucld(i,k,1)
           g2=(awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucld(i,k,2)
           g3= arg(1)*taucld(i,k,3)
           asyclt(i)=(g1+g2+g3)/taux

         endif
       enddo
         do i=1,m
           asycl(i,k)=asyclt(i)
         enddo
      enddo
c      
c-----pressure at middle levels in atm.    2026.5 = 1013.25*2.0   
c 
      do k=1,np
       do i=1,m   
      pa(i,k)=(pl(i,k)+pl(i,k+1))/2026.5
      enddo            
      enddo
c      
c-----integration over spectral bands

      do 100 ib=1,nband

       do k=1,np
        do i=1,m
ctar
c-----compute new rayleigh optical thicknesses

          taurs(i,k)=ry(ib)*pa(i,k)/ta(i,k)*zdel(i,k)

             if (crel.eq.'n')  taurs(i,k)=0.0
	
c-----compute gaseous optical thickness in each interval

          taugas(i,k)=acuv(i,k,ib)*zdel(i,k)

c-----compute clear-sky optical thickness, single scattering albedo,
c     and asymmetry factor (Eqs. 6.2-6.4)

ctar          tausto(i,k)=taurs+tauoz+tauwv+taual(i,k,ib)+1.0e-8
          tausto(i,k)=taurs(i,k)+taugas(i,k)+taual(i,k,ib)	  
          ssatau(i,k)=ssaal(i,k,ib)*taual(i,k,ib)+taurs(i,k)
          asysto(i,k)=asyal(i,k,ib)*ssaal(i,k,ib)*taual(i,k,ib)

c-----compute reflectance and transmittance of the clear portion of a layer

          tautob(i,k)=tausto(i,k)
          ssatob(i,k)=ssatau(i,k)/tautob(i,k)+1.0e-8
          ssatob(i,k)=min(ssatob(i,k),0.999999)
          asytob(i,k)=asysto(i,k)/(ssatob(i,k)*tautob(i,k))

        enddo
       enddo
       
!             WRITE (2009,*) '    taurs     taugas' 
C            
!             DO k=1,np
!             WRITE (2009,200) taurs(1,k), taugas(1,k)
!             ENDDO
C                             
 200    format (2E12.4)            
C
         call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)
C
C           
         call deledd (m,np,tautob,ssatob,asytob,dsm,rst,tst,dum)

       do k=1,np
        do i=1,m
           rr(i,k,1)=rrt(i,k)
           tt(i,k,1)=ttt(i,k)
           td(i,k,1)=tdt(i,k)
           rs(i,k,1)=rst(i,k)
           ts(i,k,1)=tst(i,k)
        enddo
       enddo

c-----compute reflectance and transmittance of the cloudy portion of a layer

       do k=1,np
        do i=1,m

c-----for direct incident radiation
c     The effective layer optical properties. Eqs. (6.2)-(6.4)

           tautob(i,k)=tausto(i,k)+tauclb(i,k)
           ssatob(i,k)=(ssatau(i,k)+tauclb(i,k))/tautob(i,k)+1.0e-8
           ssatob(i,k)=min(ssatob(i,k),0.999999)
           asytob(i,k)=(asysto(i,k)+asycl(i,k)*tauclb(i,k))
     *                /(ssatob(i,k)*tautob(i,k))

c-----for diffuse incident radiation

           tautof(i,k)=tausto(i,k)+tauclf(i,k)
           ssatof(i,k)=(ssatau(i,k)+tauclf(i,k))/tautof(i,k)+1.0e-8
           ssatof(i,k)=min(ssatof(i,k),0.999999)
           asytof(i,k)=(asysto(i,k)+asycl(i,k)*tauclf(i,k))
     *                /(ssatof(i,k)*tautof(i,k))

        enddo
       enddo

c-----for direct incident radiation
c     note that the cloud optical thickness is scaled differently for direct
c     and diffuse insolation, Eqs. (7.3) and (7.4).

         call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

         call deledd (m,np,tautof,ssatof,asytof,dsm,rst,tst,dum)

       do k=1,np
        do i=1,m

           rr(i,k,2)=rrt(i,k)
           tt(i,k,2)=ttt(i,k)
           td(i,k,2)=tdt(i,k)
           rs(i,k,2)=rst(i,k)
           ts(i,k,2)=tst(i,k)

        enddo
       enddo

c-----flux calculations

c     initialize clear-sky flux (fclr), all-sky flux (fall), 
c     and surface downward fluxes (fsdir and fsdif)

        do k=1,np+1
         do i=1,m
           fclr(i,k)=0.0
           fall(i,k)=0.0
           fclr_d(i,k)=0.0  ! NEW           
           fclr_u(i,k)=0.0  ! NEW           
           fall_d(i,k)=0.0  ! NEW
           fall_u(i,k)=0.0  ! NEW           
           dum1(i,k)=0.0           
         enddo
        enddo

        do i=1,m
           fsdir(i)=0.0
           fsdif(i)=0.0
	   dum2(i)=0.0
        enddo

        if (overcast) then

              ih1=1
              ih2=1
              im1=1
              im2=1
              is1=1
              is2=1

c-----for clear-sky fluxes only
C
C
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,dum1,dum2,dum2,
     *                fclr_d,fclr_u,dum1,dum1)   ! NEW
C
C             

           do k=1,np+1
            do i=1,m
              fall(i,k)=0.0         
           fall_d(i,k)=0.0  ! NEW
           fall_u(i,k)=0.0  ! NEW 
	   dum1(i,k)=0.0                  
C              
            enddo
           enddo

           do i=1,m
              fsdir(i)=0.0
              fsdif(i)=0.0
           enddo

              ih1=2
              ih2=2
              im1=2
              im2=2
              is1=2
              is2=2

c-----for cloudy-sky fluxes only

         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,dum1,fall,fsdir,fsdif,
     *                 dum1,dum1,fall_d,fall_u)  ! NEW
C
        else

              ih1=1
              ih2=2
              im1=1
              im2=2
              is1=1
              is2=2

c-----for clear- and all-sky fluxes
c     the all-sky flux, fall is the summation inside the brackets
c     of Eq. (7.11)

         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,
     *               fclr_d,fclr_u,fall_d,fall_u)   ! NEW

        endif

c-----flux integration, Eq. (6.1)

       do k=1,np+1
        do i=1,m
          flx(i,k)=flx(i,k)+fall(i,k)*hk(ib)
          flx_d(i,k)=flx_d(i,k)+fall_d(i,k)*hk(ib) !NEW         
          flx_u(i,k)=flx_u(i,k)+fall_u(i,k)*hk(ib) !NEW           
        enddo

        do i=1,m
          flc(i,k)=flc(i,k)+fclr(i,k)*hk(ib)
          flc_d(i,k)=flc_d(i,k)+fclr_d(i,k)*hk(ib) !NEW           
          flc_u(i,k)=flc_u(i,k)+fclr_u(i,k)*hk(ib) !NEW           
C          
        enddo
       enddo

c-----compute direct and diffuse downward surface fluxes in the UV
c     and par regions

       if(ib.le.2) then

        do i=1,m
          fdiruv(i)=fdiruv(i)+fsdir(i)*hk(ib)
          fdifuv(i)=fdifuv(i)+fsdif(i)*hk(ib)
        enddo

       else

        do i=1,m
          fdirpar(i)=fsdir(i)*hk(ib)
          fdifpar(i)=fsdif(i)*hk(ib)
        enddo

       endif
       
 100  continue
 
!        write (2009,*)  'pres        flx_d       flc_d'
!       do k=1,np+1
!       write (2009,210) pl(1,k),flx_d(1,k),flc_d(1,k)
!       enddo
 210   format (3E12.4)
c 
      return
      end
c
c
c***********************************************************************

      subroutine solir (m,np,crel,swa,wh,soa,oh,so2,o2h,sco2,co2h,
     * dp,zdel,pl,ta,
     * overcast,cldwater,cwp,taucld,reff,ict,icb,fcld,cosz,
     *            taual,ssaal,asyal,rsirbm,rsirdf,
     *            flx,flc,fdirir,fdifir,
     *            flx_d,flx_u,flc_d,flc_u)  ! NEW      
 
c************************************************************************
c 
c     The subroutine is re-written by Tarasova for using of
c              Fomin and Correa, 2005  parameterizations, August 2005
c    swa(m,np),soa(m,np),sco2(m,np),so2(m,np)-H2O,O3,CO2,O2 concentration
c                          in molec/cm^2 km divided by  e+21,e+17,e+19,e+22
c    wh(m,np),oh(m,np),co2h(m,np),o2h(m,np) - gaseous amount at the solar beam 
c    in molec/cm^2 divided by e+21,e+17,e+19,e+22
c    zdel((m,np) - layer height in km
c
c
c  compute solar flux in the infrared region. The spectrum is divided
c   into 5 bands:
c
c          band   wavenumber(/cm)  wavelength (micron)
c          1                       0.323-1.22
c          2     14280-8200        0.70-1.22
c          3     8200-1000         1.22-10.0
c          4     8200-4400         1.22-2.27
c          5     4400-1000         2.27-10.0
c        
c
c----- Input parameters:                            units      size
c
c  number of soundings (m)                          n/d         1
c  number of atmospheric layers (np)                n/d         1
c  layer scaled-water vapor content (wh)          gm/cm^2      m*np
c  option for scaling cloud optical thickness       n/d         1
c        overcast="true" if scaling is NOT required
c        overcast="false" if scaling is required
c  input option for cloud optical thickness         n/d         1
c        cldwater="true" if taucld is provided
c        cldwater="false" if cwp is provided
c  cloud water concentration (cwp)                gm/m**2      m*np*3
c        index 1 for ice particles
c        index 2 for liquid drops
c        index 3 for rain drops
c  cloud optical thickness (taucld)                 n/d        m*np*3
c        index 1 for ice paticles
c        index 2 for liquid drops
c        index 3 for rain drops
c  effective cloud-particle size (reff)           micrometer   m*np*3
c        index 1 for ice paticles
c        index 2 for liquid drops
c        index 3 for rain drops
c  level index separating high and                  n/d        m
c        middle clouds (ict)
c  level index separating middle and                n/d        m
c        low clouds (icb)
c  cloud amount (fcld)                            fraction     m*np
c  aerosol optical thickness (taual)                n/d        m*np*11
c  aerosol single-scattering albedo (ssaal)         n/d        m*np*11
c  aerosol asymmetry factor (asyal)                 n/d        m*np*11 
c  near ir surface albedo for beam                fraction     m
c        radiation (rsirbm)
c  near ir surface albedo for diffuse             fraction     m
c        radiation (rsirdf)
c
c---- temporary array
c
c  scaled cloud optical thickness                   n/d        m*np
c          for beam radiation (tauclb)
c  scaled cloud optical thickness                   n/d        m*np
c          for diffuse radiation  (tauclf)     
c
c----- output (updated) parameters:
c
c  all-sky flux (downward-upward) (flx)           fraction     m*(np+1)
c  clear-sky flux (downward-upward) (flc)         fraction     m*(np+1)
c  all-sky direct downward ir flux at
c          the surface (fdirir)                   fraction     m
c  all-sky diffuse downward ir flux at
c          the surface (fdifir)                   fraction     m
c
c**********************************************************************
      implicit none

c-----input parameters

      integer m,np,ict,icb
      integer ih1,ih2,im1,im2,is1,is2
      real cwp(m,np,3),taucld(m,np,3),reff(m,np,3)
      real fcld(m,np),cosz(m)
      real rsirbm(m),rsirdf(m)
      real taual(m,np,8),ssaal(m,np,8),asyal(m,np,8)
      real dp(m,np),wh(m,np)
      logical overcast,cldwater
ctar      
      real pl(m,np+1),ta(m,np),zdel(m,np),tsrf(m),pa(m,np)
      real so2(m,np),o2h(m,np),sco2(m,np),co2h(m,np)
      real swa(m,np), soa(m,np), oh(m,np)
      real acir(m,np,20)                        
ctar
c-----output (updated) parameters

      real flx(m,np+1),flc(m,np+1)
      real fdirir(m),fdifir(m)
      REAL flx_d(m,np+1),flx_u(m,np+1)   !NEW
      REAL flc_d(m,np+1),flc_u(m,np+1)   !NEW        
      
c-----static parameters

ctar     
      integer nband      
ctar  
      parameter (nband=5)
      integer nk(5),ichan      
ctar    
      real hk(nband,4),ry(nband)      
      real aib(nband,2),awb(nband,2),arb(nband,2)
      real aia(nband,3),awa(nband,3),ara(nband,3)
      real aig(nband,3),awg(nband,3),arg(nband,3) 
c-----temporary array
      integer ib,iv,ik,i,k
      real dsm(m)
      real tauclb(m,np),tauclf(m,np),cc(m,3)
      real ssacl(m,np),asycl(m,np)
      real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2),
     *     rs(m,np+1,2),ts(m,np+1,2)
      real fall(m,np+1),fclr(m,np+1),fsdir(m),fsdif(m)

ctar    
      real taurs(m,np),taugas(m,np)      
      real tausto(m,np),ssatau(m,np),asysto(m,np)
      real tautob(m,np),ssatob(m,np),asytob(m,np)
      real tautof(m,np),ssatof(m,np),asytof(m,np)
      real taux,reff1,reff2,w1,w2,w3,g1,g2,g3
      real ssaclt(m),asyclt(m)
      real rrt(m,np),ttt(m,np),tdt(m,np),rst(m,np),tst(m,np)
      real dum1(m,np+1),dum2(m),dum(m,np)
      REAL fall_d(m,np+1),fall_u(m,np+1)   ! NEW
      REAL fclr_d(m,np+1),fclr_u(m,np+1)   ! NEW  
      character crel*1  
      
c
c   number of k-distribution terms in each interval
c
      data nk/1,3,1,4,3/
      data hk/0.006887,0.20318,0.01032,0.09006,0.01381,
     *  0.,0.06211,0.,0.0187,0.009589,
     *  0.,0.05358,0.,0.01906,0.01336,
     *  0.,0.,0.,0.03427,0./
c
c
c-----ry is the extinction coefficient for Rayleigh scattering.
c 
      data ry /0.5382,0.5382,0.0,0.0598,0.0/      

c-----coefficients for computing the extinction coefficients of
c     ice, water, and rain particles 
ctar     
      data aib/
     1  .000333, .000333, .000333, .000333,.000333,
     2     2.52,   2.52,   2.52, 2.52, 2.52/     
ctar     
      data awb/
     1  -0.00806,-0.0101, -0.0198, -0.0166, -0.0339,
     2     1.68, 1.72,    1.91,    1.85,   2.16/     
ctar
      data arb/
     1   0.00307, 0.00307,0.00307, 0.00307, 0.00307,
     2   0.0    , 0.0    , 0.0    , 0.0  , 0.0  /     
     
c-----coefficients for computing the single-scattering co-albedo of
c     ice, water, and rain particles (Table 5)
ctar     
      data aia/
     1 -.00000260,-.00000260, .00215346, .00215346, .08938331,
     2  .00000746,  .00000746, .00073709, .00073709,  .00299387,
     3  .00000000, .00000000,-.00000134, -.00000134, -.00001038/     
ctar     
      data awa/
     1 .00000003, .00000007,.01209318, -.00019934, .01209318,
     2  .00000354,.00000845,.01784739, .00088757, .01784739,
     3 -.000000017, -.00000004,-.00036910, -.00000650, -.00036910/     
ctar
      data ara/
     1 .029, .029,      .342,  .342,  .466,
     2 .0000,  .0000,     .000,  .000,  .000,
     3 .0000, .0000,     .000,  .000,  .000/     

c-----coefficients for computing the asymmetry factor of 
c     ice, water, and rain particles (Table 6)
ctar  aig and arg are not recalculated (Tarasova)
      data aig/
     1 .74935228, .74935228, .84090400, .76098937, .84090400,
     2  .00119715, .00119715, .00126222 , .00141864, .00126222,
     3 -.00000367,-.00000367, -.00000385 , -.00000396, -.00000385/     
ctar     
      data awg/
     1  .79375035, .79375035, .83530748 ,.74513197, .83530748,
     2  .00832441, .00832441,.00257181 ,.01370071, .00257181,
     3  -.00023263, -.00023263,.00005519, -.00038203,.00005519/     
ctar     
      data arg/
     1 .891,  .891,      .948,  .948,    .971,
     2 .0000, .0000,     .000,  .000,     .000,
     3  .0000, .0000,     .000,   .000,    .000/     

c-----initialize surface fluxes, reflectances, and transmittances.
c     the reflectance and transmittance of the clear and cloudy portions
c     of a layer are denoted by 1 and 2, respectively.
c     cc is the maximum cloud cover in each of the high, middle, and low
c     cloud groups.
c     1/dsm=1/cos(53)=1.66

      do i=1,m
         dsm(i)=0.602
         fdirir(i)=0.0
         fdifir(i)=0.0
         rr(i,np+1,1)=rsirbm(i)
         rr(i,np+1,2)=rsirbm(i)
         rs(i,np+1,1)=rsirdf(i)
         rs(i,np+1,2)=rsirdf(i)
         td(i,np+1,1)=0.0
         td(i,np+1,2)=0.0
         tt(i,np+1,1)=0.0
         tt(i,np+1,2)=0.0
         ts(i,np+1,1)=0.0
         ts(i,np+1,2)=0.0
         cc(i,1)=0.0
         cc(i,2)=0.0
         cc(i,3)=0.0
      enddo

c-----Compute absorption coefficients of H2O,O3,O2,CO2 in 11 channels
c
         do i=1,m
          tsrf(i)=ta(i,np)
         enddo	  
c  
c-----pressure at middle levels (atm.)    2026.5 = 1013.25*2.0   
c 
      do k=1,np
       do i=1,m   
      pa(i,k)=(pl(i,k)+pl(i,k+1))/2026.5
      enddo            
      enddo      
c     
       call kir (m,np,so2,soa,sco2,swa,
     * o2h,oh,co2h,wh,ta,tsrf,pa,acir)       
c       
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.323-1.22, abs: O2 '      
!      do k=1,np
!      write (2009,115) acir(m,k,1)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.7-1.22, abs: H2O,O3,O2 '      
!      do k=1,np
!      write (2009,115) acir(m,k,5)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.7-1.22, abs: H2O,O3 '      
!      do k=1,np
!      write (2009,115) acir(m,k,6)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 0.7-1.22, abs: H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,7)       
!      enddo
c
c
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 1.22-10.0, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,9)       
!      enddo
c
c
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 1.22-2.27, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,13)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 1.22-2.27, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,14)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 1.22-2.27, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,15)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 1.22-2.27, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,16)       
!      enddo
c
c
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 2.27-10.0, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,17)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 2.27-10.0, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,18)       
!      enddo
c
!      write (2009,*) 'gaseous absorption coefficients, km-1'
!      write (2009,*) ' WL: 2.27-10.0, abs: CO2,H2O '      
!      do k=1,np
!      write (2009,115) acir(m,k,19)       
!      enddo
c     
 115   format (E12.4)      
c       
c-----integration over spectral bands

      do 100 ib=1,nband

ctar      
       iv=ib+3       

c-----Compute cloud optical thickness. Eqs. (4.6) and (4.11)

      if (cldwater) then

       do k=1,np
        do i=1,m
          taucld(i,k,1)=cwp(i,k,1)*(aib(ib,1)
     *                    +aib(ib,2)/reff(i,k,1))
          taucld(i,k,2)=cwp(i,k,2)*(awb(ib,1)
     *                    +awb(ib,2)/reff(i,k,2))
          taucld(i,k,3)=cwp(i,k,3)*arb(ib,1)
        enddo
       enddo

      endif

c-----options for scaling cloud optical thickness

      if (overcast) then

       do k=1,np
        do i=1,m
          tauclb(i,k)=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          tauclf(i,k)=tauclb(i,k)
        enddo
       enddo

       do k=1,3
        do i=1,m
          cc(i,k)=1.0
        enddo
       enddo

      else

c-----scale cloud optical thickness in each layer from taucld (with
c     cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
c     tauclb is the scaled optical thickness for beam radiation and
c     tauclf is for diffuse radiation.

       call cldscale(m,np,cosz,fcld,taucld,ict,icb,
     *              cc,tauclb,tauclf)

      endif

c-----compute cloud single scattering albedo and asymmetry factor
c     for a mixture of ice and liquid particles.
c     Eqs.(4.6)-(4.8), (6.2)-(6.4)

       do k=1,np

        do i=1,m

           ssaclt(i)=0.99999
           asyclt(i)=1.0

           taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)
          if (taux.gt.0.02 .and. fcld(i,k).gt.0.01) then

           reff1=min(reff(i,k,1),130.)
           reff2=min(reff(i,k,2),20.0)

           w1=(1.-(aia(ib,1)+(aia(ib,2)+
     *         aia(ib,3)*reff1)*reff1))*taucld(i,k,1)
           w2=(1.-(awa(ib,1)+(awa(ib,2)+
     *         awa(ib,3)*reff2)*reff2))*taucld(i,k,2)
           w3=(1.- ara(ib,1))*taucld(i,k,3)
           ssaclt(i)=(w1+w2+w3)/taux

           g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
           g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
           g3= arg(ib,1)*w3
           asyclt(i)=(g1+g2+g3)/(w1+w2+w3)

          endif

        enddo

         do i=1,m
           ssacl(i,k)=ssaclt(i)
         enddo
         do i=1,m
           asycl(i,k)=asyclt(i)
         enddo

       enddo
c
        
c-----integration over the nk k-distribution functions in each ib interval

       do 200 ik=1,nk(ib)

        do k=1,np
         do i=1,m
ctar
c-----compute new rayleigh optical thicknesses

         taurs(i,k)=ry(ib)*pa(i,k)/ta(i,k)*zdel(i,k)
	 
            if (crel.eq.'n')  taurs(i,k)=0.0	 

c-----compute gaseous optical thickness in each interval

          ichan=4*(ib-1)+ik
c	 
          taugas(i,k)=acir(i,k,ichan)*zdel(i,k)          
c
c-----compute clear-sky optical thickness, single scattering albedo,
c     and asymmetry factor. Eqs.(6.2)-(6.4)
 
ctar      tausto(i,k)=taurs(i,k)+taugas(i,k)+taual(i,k,iv)+1.0e-8
      tausto(i,k)=taurs(i,k)+taugas(i,k)+taual(i,k,iv)      
      ssatau(i,k)=ssaal(i,k,iv)*taual(i,k,iv)+taurs(i,k)
      asysto(i,k)=asyal(i,k,iv)*ssaal(i,k,iv)*taual(i,k,iv)
      
c       
c-----compute reflectance and transmittance of the clear portion of a layer

           tautob(i,k)=tausto(i,k)
           ssatob(i,k)=ssatau(i,k)/tautob(i,k)+1.0e-8
           ssatob(i,k)=min(ssatob(i,k),0.999999)
           asytob(i,k)=asysto(i,k)/(ssatob(i,k)*tautob(i,k))

         enddo
        enddo	
cc 
!      WRITE (2009,*) 'ichan=', ichan, 'iband=', ib, 'ik-dist=', ik
!             WRITE (2009,*)  'taurs           taugas'         
!             DO k=1,np
!             WRITE (2009,205) taurs(1,k), taugas(1,k)
!             ENDDO
cc                              
 205    format (2E12.4) 
       
c             	

c-----for direct incident radiation

         call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, Eqs. (6.5) and (6.6)

c    
         call deledd (m,np,tautob,ssatob,asytob,dsm,rst,tst,dum)
c        

        do k=1,np
         do i=1,m
            rr(i,k,1)=rrt(i,k)
            tt(i,k,1)=ttt(i,k)
            td(i,k,1)=tdt(i,k)
            rs(i,k,1)=rst(i,k)
            ts(i,k,1)=tst(i,k)
         enddo
        enddo

c-----compute reflectance and transmittance of the cloudy portion of a layer

        do k=1,np
         do i=1,m

c-----for direct incident radiation. Eqs.(6.2)-(6.4)

           tautob(i,k)=tausto(i,k)+tauclb(i,k)
           ssatob(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclb(i,k))
     *                /tautob(i,k)+1.0e-8
           ssatob(i,k)=min(ssatob(i,k),0.999999)
           asytob(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclb(i,k))
     *                /(ssatob(i,k)*tautob(i,k))

c-----for diffuse incident radiation

           tautof(i,k)=tausto(i,k)+tauclf(i,k)
           ssatof(i,k)=(ssatau(i,k)+ssacl(i,k)*tauclf(i,k))
     *                /tautof(i,k)+1.0e-8
           ssatof(i,k)=min(ssatof(i,k),0.999999)
           asytof(i,k)=(asysto(i,k)+asycl(i,k)*ssacl(i,k)*tauclf(i,k))
     *                /(ssatof(i,k)*tautof(i,k))

         enddo
        enddo
	
c        
c-----for direct incident radiation

          call deledd (m,np,tautob,ssatob,asytob,cosz,rrt,ttt,tdt)

c-----diffuse incident radiation is approximated by beam radiation with
c     an incident angle of 53 degrees, Eqs.(6.5) and (6.6)
ctar	 
          call deledd (m,np,tautof,ssatof,asytof,dsm,rst,tst,dum)
ctar	  	  
        do k=1,np
         do i=1,m

            rr(i,k,2)=rrt(i,k)
            tt(i,k,2)=ttt(i,k)
            td(i,k,2)=tdt(i,k)
            rs(i,k,2)=rst(i,k)
            ts(i,k,2)=tst(i,k)

         enddo
        enddo

c-----flux calculations

c-----initialize clear-sky flux (fclr), all-sky flux (fall), 
c     and surface downward fluxes (fsdir and fsdif)

        do k=1,np+1
         do i=1,m
           fclr(i,k)=0.0
           fall(i,k)=0.0
           fclr_d(i,k)=0.0  ! NEW           
           fclr_u(i,k)=0.0  ! NEW           
           fall_d(i,k)=0.0  ! NEW
           fall_u(i,k)=0.0  ! NEW
	   dum1(i,k)=0.0                  
         
         enddo
        enddo
ctar	
ctar     
        do i=1,m
           fsdir(i)=0.0
           fsdif(i)=0.0
	   dum2(i)=0.0
        enddo

        if (overcast) then

              ih1=1
              ih2=1
              im1=1
              im2=1
              is1=1
              is2=1

c-----for clear-sky fluxes only
	 
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,dum1,dum2,dum2,
     *                fclr_d,fclr_u,dum1,dum1)   ! NEW     

           do k=1,np+1
            do i=1,m
              fall(i,k)=0.0
           fall_d(i,k)=0.0  ! NEW
           fall_u(i,k)=0.0  ! NEW 
	   dum1(i,k)=0.0              
C              
            enddo
           enddo

           do i=1,m
              fsdir(i)=0.0
              fsdif(i)=0.0
           enddo
              ih1=2
              ih2=2
              im1=2
              im2=2
              is1=2
              is2=2
	      
c    
c-----for cloudy-sky fluxes only

         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,dum1,fall,fsdir,fsdif,
     *                 dum1,dum1,fall_d,fall_u)  ! NEW   

        else

              ih1=1
              ih2=2
              im1=1
              im2=2
              is1=1
              is2=2

c-----for clear- and all-sky fluxes
c     the all-sky flux, fall is the summation inside the brackets
c     of Eq. (7.11)
c	 
         call cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *                cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,
     *               fclr_d,fclr_u,fall_d,fall_u)   ! NEW     

        endif

c-----flux integration following Eq. (6.1)

       do k=1,np+1
        do i=1,m
          flx_d(i,k)=flx_d(i,k)+fall_d(i,k)*hk(ib,ik) !NEW         
          flx_u(i,k)=flx_u(i,k)+fall_u(i,k)*hk(ib,ik) !NEW           
          flx(i,k) = flx(i,k)+fall(i,k)*hk(ib,ik)
          
        enddo

        do i=1,m
          flc(i,k) = flc(i,k)+fclr(i,k)*hk(ib,ik)
          flc_d(i,k)=flc_d(i,k)+fclr_d(i,k)*hk(ib,ik) !NEW           
          flc_u(i,k)=flc_u(i,k)+fclr_u(i,k)*hk(ib,ik) !NEW        
        enddo
       enddo

c-----compute downward surface fluxes in the ir region

       do i=1,m
          fdirir(i) = fdirir(i)+fsdir(i)*hk(ib,ik)
          fdifir(i) = fdifir(i)+fsdif(i)*hk(ib,ik)
       enddo

  200 continue
  100 continue
c  
!        write (2009,*)  'pres        flx_d       flc_d'
!       do k=1,np+1
!       write (2009,210) pl(1,k),flx_d(1,k),flc_d(1,k)
!       enddo
 210   format (3E12.4)  
c 
      return
      end
c      
c*******************************************************************

      subroutine cldflx (m,np,ict,icb,ih1,ih2,im1,im2,is1,is2,
     *           cc,rr,tt,td,rs,ts,fclr,fall,fsdir,fsdif,
     *   fclr_d,fclr_u,fall_d,fall_u)      ! NEW

c*******************************************************************
c   Chou and Suarez 1999
c  compute upward and downward fluxes using a two-stream adding method
c  following equations (6.9)-(6.16).
c
c  clouds are grouped into high, middle, and low clouds which are assumed
c  randomly overlapped. It involves a maximum of 8 sets of calculations.
c  In each set of calculations, each atmospheric layer is homogeneous,
c  either totally filled with clouds or without clouds.

c  input parameters:
c
c   m:   number of soundings
c   np:  number of atmospheric layers
c   ict: the level separating high and middle clouds
c   icb: the level separating middle and low clouds
c   ih1,ih2,im1,im2,is1,is2: indices for three group of clouds
c   cc:  effective cloud covers for high, middle and low clouds
c   rr:  reflection of a layer illuminated by beam radiation
c   tt:  total diffuse transmission of a layer illuminated by beam radiation
c   td:  direct beam transmission
c   rs:  reflection of a layer illuminated by diffuse radiation
c   ts:  transmission of a layer illuminated by diffuse radiation
c
c  output parameters:
c
c     fclr:  clear-sky flux (downward minus upward)
c     fall:  all-sky flux (downward minus upward)
c     fsdir: surface direct downward flux
c     fsdif: surface diffuse downward flux
c
c*********************************************************************c

      implicit none

c-----input parameters

      integer m,np,ict,icb,ih1,ih2,im1,im2,is1,is2

      real rr(m,np+1,2),tt(m,np+1,2),td(m,np+1,2)
      real rs(m,np+1,2),ts(m,np+1,2)
      real cc(m,3)

c-----temporary array

      integer i,k,ih,im,is
      real rra(m,np+1,2,2),tta(m,np+1,2,2),tda(m,np+1,2,2)
      real rsa(m,np+1,2,2),rxa(m,np+1,2,2)
      real ch(m),cm(m),ct(m),flxdn(m,np+1)
      real fdndir(m),fdndif(m),fupdif
      real denm,xx,yy
      REAL fdn(m,np+1), fup(m,np+1)   !NEW      

c-----output parameters

      real fclr(m,np+1),fall(m,np+1)
      REAL fclr_u(m,np+1),fclr_d(m,np+1) !NEW
      REAL fall_u(m,np+1), fall_d (m,np+1)  !NEW    
      real fsdir(m),fsdif(m)

c-----compute transmittances and reflectances for a composite of
c     layers. layers are added one at a time, going down from the top.
c     tda is the composite transmittance illuminated by beam radiation
c     tta is the composite total transmittance illuminated by
c         beam radiation
c     rsa is the composite reflectance illuminated from below
c         by diffuse radiation
c     tta and rsa are computed from Eqs. (6.10) and (6.12)

c-----for high clouds
c     ih=1 for clear-sky condition, ih=2 for cloudy-sky condition

      do ih=ih1,ih2

       do i=1,m
          tda(i,1,ih,1)=td(i,1,ih)
          tta(i,1,ih,1)=tt(i,1,ih)
          rsa(i,1,ih,1)=rs(i,1,ih)
          tda(i,1,ih,2)=td(i,1,ih)
          tta(i,1,ih,2)=tt(i,1,ih)
          rsa(i,1,ih,2)=rs(i,1,ih)
       enddo

       do k=2,ict-1
        do i=1,m
          denm = ts(i,k,ih)/( 1.-rsa(i,k-1,ih,1)*rs(i,k,ih))
          tda(i,k,ih,1)= tda(i,k-1,ih,1)*td(i,k,ih)
          tta(i,k,ih,1)= tda(i,k-1,ih,1)*tt(i,k,ih)
     *          +(tda(i,k-1,ih,1)*rsa(i,k-1,ih,1)*rr(i,k,ih)
     *          +tta(i,k-1,ih,1)-tda(i,k-1,ih,1))*denm
          rsa(i,k,ih,1)= rs(i,k,ih)+ts(i,k,ih)
     *                  *rsa(i,k-1,ih,1)*denm
          tda(i,k,ih,2)= tda(i,k,ih,1)
          tta(i,k,ih,2)= tta(i,k,ih,1)
          rsa(i,k,ih,2)= rsa(i,k,ih,1)
        enddo

       enddo

c-----for middle clouds
c     im=1 for clear-sky condition, im=2 for cloudy-sky condition
      do im=im1,im2

       do k=ict,icb-1
        do i=1,m
          denm = ts(i,k,im)/( 1.-rsa(i,k-1,ih,im)*rs(i,k,im))
          tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,im)
          tta(i,k,ih,im)= tda(i,k-1,ih,im)*tt(i,k,im)
     *         +(tda(i,k-1,ih,im)*rsa(i,k-1,ih,im)*rr(i,k,im)
     *         +tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
          rsa(i,k,ih,im)= rs(i,k,im)+ts(i,k,im)
     *                  *rsa(i,k-1,ih,im)*denm
        enddo
       enddo

      enddo                 ! end im loop
      enddo                 ! end ih loop



c-----layers are added one at a time, going up from the surface.
c     rra is the composite reflectance illuminated by beam radiation
c     rxa is the composite reflectance illuminated from above
c         by diffuse radiation
c     rra and rxa are computed from Eqs. (6.9) and (6.11)
 
c-----for the low clouds
c     is=1 for clear-sky condition, is=2 for cloudy-sky condition

      do is=is1,is2

       do i=1,m
         rra(i,np+1,1,is)=rr(i,np+1,is)
         rxa(i,np+1,1,is)=rs(i,np+1,is)
         rra(i,np+1,2,is)=rr(i,np+1,is)
         rxa(i,np+1,2,is)=rs(i,np+1,is)
       enddo

       do k=np,icb,-1
        do i=1,m
          denm=ts(i,k,is)/( 1.-rs(i,k,is)*rxa(i,k+1,1,is) )
          rra(i,k,1,is)=rr(i,k,is)+(td(i,k,is)*rra(i,k+1,1,is)
     *        +(tt(i,k,is)-td(i,k,is))*rxa(i,k+1,1,is))*denm
          rxa(i,k,1,is)= rs(i,k,is)+ts(i,k,is)
     *        *rxa(i,k+1,1,is)*denm
          rra(i,k,2,is)=rra(i,k,1,is)
          rxa(i,k,2,is)=rxa(i,k,1,is)
        enddo
       enddo

c-----for middle clouds

      do im=im1,im2

       do k=icb-1,ict,-1
        do i=1,m
          denm=ts(i,k,im)/( 1.-rs(i,k,im)*rxa(i,k+1,im,is) )
          rra(i,k,im,is)= rr(i,k,im)+(td(i,k,im)*rra(i,k+1,im,is)
     *        +(tt(i,k,im)-td(i,k,im))*rxa(i,k+1,im,is))*denm
          rxa(i,k,im,is)= rs(i,k,im)+ts(i,k,im)
     *        *rxa(i,k+1,im,is)*denm
        enddo
       enddo

      enddo                 ! end im loop
      enddo                 ! end is loop

c-----integration over eight sky situations.
c     ih, im, is denotes high, middle and low cloud groups.

      do ih=ih1,ih2

c-----clear portion 

         if(ih.eq.1) then
           do i=1,m
             ch(i)=1.0-cc(i,1)
           enddo

          else

c-----cloudy portion

           do i=1,m
             ch(i)=cc(i,1)
           enddo

          endif

      do im=im1,im2

c-----clear portion

         if(im.eq.1) then

           do i=1,m
              cm(i)=ch(i)*(1.0-cc(i,2))
           enddo

         else

c-----cloudy portion

           do i=1,m
              cm(i)=ch(i)*cc(i,2) 
           enddo

         endif

      do is=is1,is2

c-----clear portion

         if(is.eq.1) then

           do i=1,m
             ct(i)=cm(i)*(1.0-cc(i,3)) 
           enddo

         else

c-----cloudy portion

           do i=1,m
             ct(i)=cm(i)*cc(i,3)
           enddo

         endif

c-----add one layer at a time, going down.

       do k=icb,np
        do i=1,m
          denm = ts(i,k,is)/( 1.-rsa(i,k-1,ih,im)*rs(i,k,is) )
          tda(i,k,ih,im)= tda(i,k-1,ih,im)*td(i,k,is)
          tta(i,k,ih,im)=  tda(i,k-1,ih,im)*tt(i,k,is)
     *         +(tda(i,k-1,ih,im)*rr(i,k,is)
     *         *rsa(i,k-1,ih,im)+tta(i,k-1,ih,im)-tda(i,k-1,ih,im))*denm
          rsa(i,k,ih,im)= rs(i,k,is)+ts(i,k,is)
     *         *rsa(i,k-1,ih,im)*denm
        enddo
       enddo

c-----add one layer at a time, going up.

       do k=ict-1,1,-1
        do i=1,m
          denm =ts(i,k,ih)/(1.-rs(i,k,ih)*rxa(i,k+1,im,is))
          rra(i,k,im,is)= rr(i,k,ih)+(td(i,k,ih)*rra(i,k+1,im,is)
     *        +(tt(i,k,ih)-td(i,k,ih))*rxa(i,k+1,im,is))*denm
          rxa(i,k,im,is)= rs(i,k,ih)+ts(i,k,ih)
     *        *rxa(i,k+1,im,is)*denm
        enddo
       enddo

c-----compute fluxes following Eq. (6.15) for fupdif and
c     Eq. (6.16) for (fdndir+fdndif)
 
c     fdndir is the direct  downward flux
c     fdndif is the diffuse downward flux
c     fupdif is the diffuse upward flux

      do k=2,np+1
       do i=1,m
         denm= 1./(1.-rsa(i,k-1,ih,im)*rxa(i,k,im,is))
         fdndir(i)= tda(i,k-1,ih,im)
         xx= tda(i,k-1,ih,im)*rra(i,k,im,is)
         yy= tta(i,k-1,ih,im)-tda(i,k-1,ih,im)
         fdndif(i)= (xx*rsa(i,k-1,ih,im)+yy)*denm
C
C--CALCULATION OF DOWNWARD (fdn(i,k)) and UPWARD (fupdif(i,k)) FLUXES
C  as well as net flux
C 
         fupdif= (xx+yy*rxa(i,k,im,is))*denm         
         flxdn(i,k)= fdndir(i)+fdndif(i)-fupdif  ! NET FLUX
C          flxdn(i,k)= fdndir(i)+fdndif(i)
          fup(i,k)= (xx+yy*rxa(i,k,im,is))*denm   !NEW         
          fdn(i,k)= fdndir(i)+fdndif(i)        !NEW 
       enddo
      enddo

       do i=1,m
          flxdn(i,1)=1.0-rra(i,1,im,is)
C           flxdn(i,1)=1.0
          fdn(i,1)=1.0 
          fup(i,1)=rra(i,1,im,is)                  
       enddo

c-----summation of fluxes over all sky situations;
c     the term in the brackets of Eq. (7.11)

       do k=1,np+1
        do i=1,m
           if(ih.eq.1 .and. im.eq.1 .and. is.eq.1) then
             fclr(i,k)=flxdn(i,k)
             fclr_d(i,k)=fdn(i,k)      ! NEW
             fclr_u(i,k)=fup(i,k)   ! NEW
           endif
             fall(i,k)=fall(i,k)+flxdn(i,k)*ct(i)
             fall_d(i,k)=fall_d(i,k)+fdn(i,k)*ct(i)      ! NEW             
             fall_u(i,k)=fall_u(i,k)+fup(i,k)*ct(i)   ! NEW            
        enddo
       enddo

        do i=1,m
            fsdir(i)=fsdir(i)+fdndir(i)*ct(i)
            fsdif(i)=fsdif(i)+fdndif(i)*ct(i)
        enddo

       enddo                 ! end is loop
       enddo                 ! end im loop
       enddo                 ! end ih loop

      return
      end
c
c********************************************************************

      subroutine cldscale (m,np,cosz,fcld,taucld,ict,icb,
     *                     cc,tauclb,tauclf)

c********************************************************************
c   Chou and Suarez 1999
c   This subroutine computes the high, middle, and low cloud
c    amounts and scales the cloud optical thickness (section 7)
c
c   To simplify calculations in a cloudy atmosphere, clouds are
c    grouped into high, middle and low clouds separated by the levels
c    ict and icb (level 1 is the top of the model atmosphere).
c
c   Within each of the three groups, clouds are assumed maximally
c    overlapped, and the cloud cover (cc) of a group is the maximum
c    cloud cover of all the layers in the group.  The optical thickness
c    (taucld) of a given layer is then scaled to new values (tauclb and
c    tauclf) so that the layer reflectance corresponding to the cloud
c    cover cc is the same as the original reflectance with optical
c    thickness taucld and cloud cover fcld.
c
c---input parameters
c
c    number of atmospheric soundings (m)
c    number of atmospheric layers (np)
c    cosine of the solar zenith angle (cosz)
c    fractional cloud cover (fcld)
c    cloud optical thickness (taucld)
c    index separating high and middle clouds (ict)
c    index separating middle and low clouds (icb)
c
c---output parameters
c
c    fractional cover of high, middle, and low cloud groups (cc)
c    scaled cloud optical thickness for direct  radiation (tauclb)
c    scaled cloud optical thickness for diffuse radiation (tauclf)
c
c********************************************************************

      implicit none

c-----input parameters

      integer m,np,ict,icb
      real cosz(m),fcld(m,np),taucld(m,np,3)

c-----output parameters

      real cc(m,3),tauclb(m,np),tauclf(m,np)

c-----temporary variables

      integer i,j,k,im,it,ia,kk
      real  fm,ft,fa,xai,taux

c-----pre-computed table
c     size of cosz-interval:         dm
c     size of taucld-interval:       dt
c     size of cloud amount-interval: da

      integer   nm,nt,na
      parameter (nm=11,nt=9,na=11) 
      real  dm,dt,da,t1,caib(nm,nt,na),caif(nt,na)
      parameter (dm=0.1,dt=0.30103,da=0.1,t1=-0.9031)

c-----include the pre-computed table of mcai for scaling the cloud optical
c     thickness under the assumption that clouds are maximally overlapped
c
c     caib is for scaling the cloud optical thickness for direct radiation
c     caif is for scaling the cloud optical thickness for diffuse radiation

      include "optics/mcai.data"

c-----clouds within each of the high, middle, and low clouds are assumed
c     to be maximally overlapped, and the cloud cover (cc) for a group
c     (high, middle, or low) is the maximum cloud cover of all the layers
c     within a group

      do i=1,m
         cc(i,1)=0.0
         cc(i,2)=0.0
         cc(i,3)=0.0
      enddo

      do k=1,ict-1
       do i=1,m
          cc(i,1)=max(cc(i,1),fcld(i,k))
       enddo
      enddo

      do k=ict,icb-1
       do i=1,m
          cc(i,2)=max(cc(i,2),fcld(i,k))
       enddo
      enddo

      do k=icb,np
       do i=1,m
          cc(i,3)=max(cc(i,3),fcld(i,k))
       enddo
      enddo

c-----scale the cloud optical thickness.
c     taucld(i,k,1) is the optical thickness for ice particles
c     taucld(i,k,2) is the optical thickness for liquid particles
c     taucld(i,k,3) is the optical thickness for rain drops
      
      do k=1,np

         if(k.lt.ict) then
            kk=1
         elseif(k.ge.ict .and. k.lt.icb) then
            kk=2
         else
            kk=3
         endif

       do i=1,m

         tauclb(i,k) = 0.0
         tauclf(i,k) = 0.0
         taux=taucld(i,k,1)+taucld(i,k,2)+taucld(i,k,3)

         if (taux.gt.0.02 .and. fcld(i,k).gt.0.01) then

c-----normalize cloud cover following Eq. (7.8)

           fa=fcld(i,k)/cc(i,kk)

c-----table look-up

           taux=min(taux,32.)

           fm=cosz(i)/dm
           ft=(log10(taux)-t1)/dt
           fa=fa/da
 
           im=int(fm+1.5)
           it=int(ft+1.5)
           ia=int(fa+1.5)
  
           im=max(im,2)
           it=max(it,2)
           ia=max(ia,2)
     
           im=min(im,nm-1)
           it=min(it,nt-1)
           ia=min(ia,na-1)

           fm=fm-float(im-1)
           ft=ft-float(it-1)
           fa=fa-float(ia-1)

c-----scale cloud optical thickness for beam radiation following Eq. (7.3)
c     the scaling factor, xai, is a function of the solar zenith
c     angle, optical thickness, and cloud cover.
 
           xai=    (-caib(im-1,it,ia)*(1.-fm)+
     *      caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)
         
           xai=xai+(-caib(im,it-1,ia)*(1.-ft)+
     *      caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

           xai=xai+(-caib(im,it,ia-1)*(1.-fa)+
     *     caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

           xai= xai-2.*caib(im,it,ia)
           xai=max(xai,0.0)
     
           tauclb(i,k) = taux*xai

c-----scale cloud optical thickness for diffuse radiation following Eq. (7.4)
c     the scaling factor, xai, is a function of the cloud optical
c     thickness and cover but not the solar zenith angle.

           xai=    (-caif(it-1,ia)*(1.-ft)+
     *      caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)

           xai=xai+(-caif(it,ia-1)*(1.-fa)+
     *      caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)

           xai= xai-caif(it,ia)
           xai=max(xai,0.0)
     
           tauclf(i,k) = taux*xai

         endif

       enddo
      enddo

      return
      end
c
c*********************************************************************

      subroutine deledd(m,np,tau,ssc,g0,cza,rr,tt,td)

c*********************************************************************
c     Chou and Suarez, 1999
c-----uses the delta-eddington approximation to compute the
c     bulk scattering properties of a single layer
c     coded following King and Harshvardhan (JAS, 1986)
c
c  inputs:
c       m:  number of soundings
c      np:  number of atmospheric layers
c     tau:  optical thickness
c     ssc:  single scattering albedo
c     g0:   asymmetry factor
c     cza:  cosine of the zenith angle
c
c  outputs:
c
c     rr:  reflection of the direct beam
c     tt:  total diffuse transmission of the direct beam
c     td:  direct transmission of the direct beam
c
c*********************************************************************
C
C
      implicit none

      real zero,one,two,three,four,fourth,seven,thresh
      parameter (one =1., three=3.)
      parameter (two =2., seven=7.)
      parameter (four=4., fourth=.25)
      parameter (zero=0., thresh=1.e-8)

c-----input parameters
      integer m,np
      real tau(m,np),ssc(m,np),g0(m,np),cza(m)

c-----output parameters
      real rr(m,np+1),tt(m,np+1),td(m,np+1)

c-----temporary parameters

      integer i,k
      real zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,
     *     all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4
 
c---------------------------------------------------------------------
C
C
      do k=1,np
       do i=1,m
C
          zth = cza(i)
 
c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  K & H eqs(27-29)

           ff  = g0(i,k)*g0(i,k)
           xx  = one-ff *ssc(i,k)
           taup= tau(i,k)*xx
           sscp= ssc(i,k)*(one-ff)/xx
           gp  = g0(i,k) /(one+g0(i,k))           
C            
c  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.

           xx  =  three*gp 
           gm1 =  (seven - sscp*(four+xx))*fourth
           gm2 = -(one   - sscp*(four-xx))*fourth

c  akk is k as defined in eq(25) of K & H
 
           akk = sqrt((gm1+gm2)*(gm1-gm2))
 
           xx  = akk * zth
           st7 = one - xx
           st8 = one + xx
           st3 = st7 * st8

           if (abs(st3) .lt. thresh) then
               zth = zth + 0.001
               xx  = akk * zth
               st7 = one - xx
               st8 = one + xx
               st3 = st7 * st8
           endif
C
C 
c  extinction of the direct beam transmission
 
           td(i,k)  = exp(-taup/zth)

c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H
 
           gm3  = (two - zth*three*gp)*fourth
           xx   = gm1 - gm2
           alf1 = gm1 - gm3 * xx
           alf2 = gm2 + gm3 * xx
 
c  all is last term in eq(21) of K & H
c  bll is last term in eq(22) of K & H
 
           xx  = akk * two
           all = (gm3 - alf2 * zth    )*xx*td(i,k)
           bll = (one - gm3 + alf1*zth)*xx
 
           xx  = akk * gm3
           cll = (alf2 + xx) * st7
           dll = (alf2 - xx) * st8
 
           xx  = akk * (one-gm3)
           fll = (alf1 + xx) * st8
           ell = (alf1 - xx) * st7
  
           st2 = exp(-akk*taup)
           st4 = st2 * st2
 
           st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)
 
c  rr is r-hat of eq(21) of K & H
c  tt is diffuse part of t-hat of eq(22) of K & H
 
           rr(i,k) =   ( cll-dll*st4         -all*st2)*st1
           tt(i,k) = - ((fll-ell*st4)*td(i,k)-bll*st2)*st1
 
           rr(i,k) = max(rr(i,k),zero)
           tt(i,k) = max(tt(i,k),zero)

           tt(i,k) = tt(i,k)+td(i,k)
C
C
        enddo
       enddo
C
C
      return
      end


      
      subroutine kuv (m,np,nband,uo3,uo2,uh2o,uo3s,uo2s,
     * uh2os,acuv)
!      
!  Parameterizations of Fomin and Correa, 2005:  
!  The subroutine is written by Tarasova Aug. 2005
!     uv and visible solar spectrum, 4 intervals
!                     input parameters:
!     m - number of model grid points
!     np - number of levels
!     nband - number of spectral intervals
!    uo3(m,np),uo2(m,np),uh2o(m,np) - concentration of O3,O2,H2O  in molec/(cm2 km)
!    divided by 1.0e17, 1.0e22, and 1.0e21, respectively
!    uo3s(m,np),uo2s(m,np),uh2os(m,np) - amount of O3,O2,H2O  in molec/(cm2)
!    at the solar beam divided by 1.0e17, 1.0e22, and 1.0e21, respectively
!                     output parameter: 
!    acuv(m,np,nband) - volume absorption coefficient (km-1)
!      
      implicit none             
      integer m, np, nband
      integer i,k,l
      real uo3(m,np), uo2(m,np), uh2o(m,np)
      real uo3s(m,np), uo2s(m,np), uh2os(m,np)       
      real acuv(m,np,nband)
      real u1o3(22),ac1o3l(22),u1o2(7),ac1o2l(7)
      real u2o3(8),ac2o3l(8),u4o2(8),ac4o2l(8)      
      real ac1o3,ac1o2,ac2o3,ac4o2      
!      
      data u1o3/1.0e-6,0.02,0.05,0.1,0.2,0.4,0.6,1.0,2.0,3.0,4.0, 
     * 5.0,6.0,7.0,9.0,15.0,20.0,30.0,100.0,200.0,250.,400.0/      
      data ac1o3l/0.365,0.365,0.362,0.355,0.345,0.32,0.3,0.26,
     * 0.185,0.15,0.12,0.105,0.095,0.09,0.08,0.063,0.057,0.047,0.032,
     * 0.026,0.023,0.005/
      data u1o2/1.0e-7,0.1,1.0,4.0,10.0,70.0,200.0/
      data ac1o2l/3.4e-3,3.4e-3,3.3e-3,3.05e-3,2.4e-3,0.6e-3,0.25e-3/
      data u2o3/1.0e-6,7.0e-4,0.4,2.0,7.0,20.0,50.0,100.0/
      data ac2o3l/7.35e-3,7.35e-3,8.35e-3,8.0e-3,7.5e-3,
     * 7.0e-3,6.0e-3,5.5e-3/
!ctar      data u4o2/0.01,0.05,0.1,0.3,1.0,3.0,10.0,100.0/
!ctar      data ac4o2l/1.9e-2,1.8e-2,1.5e-2,1.0e-2,0.7e-2,0.48e-2,
!ctar     * 0.25e-2,0.25e-2/            
!             
!             
      do i=1,m
      do k=1,np
! 
      ac1o3=ac1o3l(1)      
      if (uo3s(i,k).gt.u1o3(22)) ac1o3=ac1o3l(22)
!                 
      do l=1,21      
      if (uo3s(i,k).gt.u1o3(l).and.uo3s(i,k).le.u1o3(l+1)) then 
      ac1o3=ac1o3l(l)+(ac1o3l(l)-ac1o3l(l+1))*(uo3s(i,k)-
     * u1o3(l))/(u1o3(l)-u1o3(l+1)) 
      endif
      enddo
      
!      write (*,*) 'ac1o3='
!      write (*,30) ac1o3
!
      ac1o2=ac1o2l(1)
      if (uo2s(i,k).gt.u1o2(7)) ac1o2=ac1o2l(7) 
      do l=1,6      
      if (uo2s(i,k).gt.u1o2(l).and.uo2s(i,k).le.u1o2(l+1)) then 
      ac1o2=ac1o2l(l)+(ac1o2l(l)-ac1o2l(l+1))*(uo2s(i,k)-
     * u1o2(l))/(u1o2(l)-u1o2(l+1)) 
      endif      
      enddo
!
!      write (*,*) 'ac1o2='
!      write (*,30) ac1o2
      
      acuv(i,k,1)=ac1o3*uo3(i,k)+ac1o2*uo2(i,k)
!
      ac2o3=ac2o3l(1)
      if (uo3s(i,k).gt.u2o3(8)) ac2o3=ac2o3l(8)       
      do l=1,7      
      if (uo3s(i,k).gt.u2o3(l).and.uo3s(i,k).le.u2o3(l+1)) then 
      ac2o3=ac2o3l(l)+(ac2o3l(l)-ac2o3l(l+1))*(uo3s(i,k)-
     * u2o3(l))/(u2o3(l)-u2o3(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac2o3='
!      write (*,30) ac2o3      
      
!      
      acuv(i,k,2)=ac2o3*uo3(i,k) 
      acuv(i,k,3)=1.93e-4*uo3(i,k)+2.0e-5*uh2o(i,k)            
!
!tar      ac4o2=ac4o2l(1)
!tar      if (uo2s(i,k).gt.u4o2(8)) ac4o2=ac4o2l(8)
!tar      do l=1,7      
!tar      if (uo2s(i,k).gt.u4o2(l).and.uo2s(i,k).le.u4o2(l+1)) then 
!tar      ac4o2=ac4o2l(l)+(ac4o2l(l)-ac4o2l(l+1))*(uo2s(i,k)-
!tar     * u4o2(l))/(u4o2(l)-u4o2(l+1)) 
!tar      endif      
!tar      enddo
      
!      write (*,*) 'ac4o2='
!      write (*,30) ac4o2      
      
!      
!tar      acuv(i,k,4)=ac4o2*uo2(i,k)            
      
      enddo
      enddo
      
  30   format (1E12.4)
   
      return
      end
      
      
      subroutine kir (m,np,uo2,uo3,uco2,uh2o,
     * uo2s,uo3s,uco2s,uh2os,temp,ts,pres,acir)
!      
!  Parameterizations of Fomin and Correa, 2005:  
!   The subroutine is written by Tarasova Aug. 2005
!     infrared solar spectrum, 4 intervals, 11 channels
!               input parameters:
!     m - number of model grid points
!     np - number of levels
!    uco2(m,np),uh2o(m,np) - concentration of CO2,H2O  in molec/(cm2 km) 
!    divided by  1.0e+19 and 1.0e+21
!    uo3(m,np),uo2(m,np)  -  concentration of O3,O2  in molec/(cm2 km) 
!    divided by 1.0e+17 and 1.0e+22
!    uco2s(m,np),uh2os(m,np),uo3s(m,np,),uo2s(m,np), - amount of CO2,H2O,O3,O2 
!    at the solar beam in molec/cm2 divided by 1.0e+19, 1.0e+21, 1.0e+17 and 1.0e+22
!    temp(m,np),ts(m),pres(m,np) - temperature in K, surface temperature in K, pressure in atm
!               output parameter:
!    acir(m,np,20) - volume absorption coefficient (km-1) in 20 channels 5*4
!    
!        
        implicit none  
      integer m, np
      integer i,k,l
      real uco2(m,np),uh2o(m,np),uo2(m,np),uo3(m,np)
      real acir(m,np,20)     
      real uco2s(m,np),uh2os(m,np),uo2s(m,np),uo3s(m,np)      
      real temp(m,np), ts(m), pres(m,np), tmpr, tmpr1
!new
      real ac4o2,u4o2(8),ac4o2l(8)      
      real u7h2o(9),ac7h2ol(9),u8h2o(15),ac8h2ol(15)      
      real u8co2(10),ac8co2l(10),u10h2o(10),ac10h2ol(10)
      real u11h2o(8),ac11h2ol(8),u12h2o(12),ac12h2ol(12)
      real u14h2o(11),ac14h2ol(11),u15h2o(10),ac15h2ol(10)
      
      real ac7h2o,ac8h2o,ac8co2,ac10h2o,ac11h2o,ac12h2o,
     *  ac14h2o,ac15h2o           
!
      data u7h2o/0.001,0.01,0.03,0.1,1.0,
     * 3.0,20.0,100.0,200.0/
      data ac7h2ol/7.61e-2,6.85e-2,6.52e-2,5.98e-2,3.81e-2,
     * 2.72e-2,1.63e-2,1.09e-2,0.978e-2/
      
      data u8h2o/1.0e-6,5.0e-7,1.0e-5,4.0e-5,1.0e-4,5.0e-4,
     * 0.001,0.005,0.01,0.07,0.4,2.0,8.0,20.0,1.0e+2/
      data ac8h2ol/4.0,6.0,7.0,6.0,4.5,
     * 2.0,1.4,1.78,1.75,0.8,0.3,0.05,
     * 0.01,0.004,0.001/ 
         
      data u8co2/4.0e-5,1.0e-4,4.0e-4,0.002,0.01,0.1,
     * 1.0,10.0,100.0,1.0e+3/      
      data ac8co2l/0.8,0.5,0.25,0.1,0.06,
     * 0.02,0.008,0.0045,0.0025,0.0013/
     
      data u10h2o/0.01,0.1,1.0,4.0,10.0,40.0,
     * 70.0,100.0,200.0,1000.0/     
      data ac10h2ol/1.8e-3,2.34e-3,3.19e-3,4.03e-3,4.67e-3,
     * 6.05e-3,6.26e-3,6.37e-3,6.26e-3,5.84e-3/

      data u11h2o/0.01,0.1,1.0,10.0,20.0,
     * 40.0,100.0,1000.0/    
      data ac11h2ol/0.823e-2,1.65e-2,2.35e-2,3.29e-2,3.53e-2,
     * 3.47e-2,2.7e-2,1.76e-2/
     
      data u12h2o/1.0e-4,0.001,0.005,0.0105,0.02,
     * 0.2,0.7,2.0,5.0,20.0,100.0,1.0e+3/     
      data ac12h2ol/1.7,1.22,0.85,0.878,0.756,
     * 0.425,0.255,0.189,0.142,0.085,0.0283,9.45e-4/      

      data u14h2o/1.0e-3,0.01,0.1,1.0,10.0,20.0,
     * 40.0,100.0,150.0,200.0,1000.0/ 
      data ac14h2ol/0.209e-2,0.313e-2,0.584e-2,0.834e-2,1.3e-2,
     * 1.46e-2,1.49e-2,1.25e-2,1.1e-2,1.04e-2,0.834e-2/
      
      data u15h2o/1.0e-5,1.0e-4,0.001,0.01,0.1,1.0,
     * 3.0,10.0,40.0,100.0/     
      data ac15h2ol/7.49,4.16,1.17,1.66,0.749,
     * 0.25,0.125,0.0666,0.0499,0.025/ 
      data u4o2/0.01,0.05,0.1,0.3,1.0,3.0,10.0,100.0/
      data ac4o2l/1.9e-2,1.8e-2,1.5e-2,1.0e-2,0.7e-2,0.48e-2,
     * 0.25e-2,0.25e-2/            
!      
      do i=1,m
      do k=1,np
!
      ac4o2=ac4o2l(1)
      if (uo2s(i,k).gt.u4o2(8)) ac4o2=ac4o2l(8)
      do l=1,7      
      if (uo2s(i,k).gt.u4o2(l).and.uo2s(i,k).le.u4o2(l+1)) then 
      ac4o2=ac4o2l(l)+(ac4o2l(l)-ac4o2l(l+1))*(uo2s(i,k)-
     * u4o2(l))/(u4o2(l)-u4o2(l+1)) 
      endif      
      enddo
!tar
      acir(i,k,1)=ac4o2*uo2(i,k)
      acir(i,k,2)=0.      
      acir(i,k,3)=0.
      acir(i,k,4)=0.           
!tar  
        tmpr=3.2**((300.0-temp(i,k))/43.0)
	tmpr1=3.2**((300.0-temp(i,k))/25.0)
	
      acir(i,k,5)=1.35e-4*(1.0-0.0007*(300.0-temp(i,k)))*
     * (pres(i,k)**0.77)*uh2o(i,k)+1.19e-5*(pres(i,k)**0.59)*uo2(i,k)+
     * 1.7e-5*uo3(i,k) 
      
      acir(i,k,6)=2.1e-3*(pres(i,k)**0.55)*uh2o(i,k)+
     * 1.4e-5*uo3(i,k)	
!
      ac7h2o=ac7h2ol(1)                 
      if (uh2os(i,k).gt.u7h2o(9)) ac7h2o=ac7h2ol(9)
             
      do l=1,8      
      if (uh2os(i,k).gt.u7h2o(l).and.uh2os(i,k).le.u7h2o(l+1)) then 
      ac7h2o=ac7h2ol(l)+(ac7h2ol(l)-ac7h2ol(l+1))*(uh2os(i,k)-
     * u7h2o(l))/(u7h2o(l)-u7h2o(l+1)) 
      endif
      enddo
      
!      write (*,*) 'ac7h2o='
!      write (*,30) ac7h2o
      
      acir(i,k,7)=ac7h2o*uh2o(i,k)
      acir(i,k,8)=0.      
!
      ac8h2o=ac8h2ol(1)
      if (uh2os(i,k).gt.u8h2o(15)) ac8h2o=ac8h2ol(7) 
      do l=1,14      
      if (uh2os(i,k).gt.u8h2o(l).and.uh2os(i,k).le.u8h2o(l+1)) then 
      ac8h2o=ac8h2ol(l)+(ac8h2ol(l)-ac8h2ol(l+1))*(uh2os(i,k)-
     * u8h2o(l))/(u8h2o(l)-u8h2o(l+1)) 
      endif      
      enddo
!
!      write (*,*) 'ac8h2o='
!      write (*,30) ac8h2o
              
!      
      ac8co2=ac8co2l(1)
      if (uco2s(i,k).gt.u8co2(10)) ac8co2=ac8co2l(10) 
      do l=1,9      
      if (uco2s(i,k).gt.u8co2(l).and.uco2s(i,k).le.u8co2(l+1)) then 
      ac8co2=ac8co2l(l)+(ac8co2l(l)-ac8co2l(l+1))*(uco2s(i,k)-
     * u8co2(l))/(u8co2(l)-u8co2(l+1)) 
      endif      
      enddo
!
!      write (*,*) 'ac8co2='
!      write (*,30) ac8co2      
                  
       acir(i,k,9)=ac8h2o*uh2o(i,k)+ac8co2*uco2(i,k) 
       acir(i,k,10)=0. 
       acir(i,k,11)=0.
       acir(i,k,12)=0.                          
!
      acir(i,k,13)=3.77e-4*(1.0-0.00372*(300.0-temp(i,k)))*
     * (pres(i,k)**0.65)*uh2o(i,k)+6.5e-7*tmpr*
     * (uh2o(i,k)**2.0)+1.55e-5*(pres(i,k)**0.54)*uco2(i,k)
!
!
      ac10h2o=ac10h2ol(1)
      if (uh2os(i,k).gt.u10h2o(10)) ac10h2o=ac10h2ol(10)       
      do l=1,9      
      if (uh2os(i,k).gt.u10h2o(l).and.uh2os(i,k).le.u10h2o(l+1)) then 
      ac10h2o=ac10h2ol(l)+(ac10h2ol(l)-ac10h2ol(l+1))*(uh2os(i,k)-
     * u10h2o(l))/(u10h2o(l)-u10h2o(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac10h2o='
!      write (*,30) ac10h2o      
      
!      
      acir(i,k,14)=ac10h2o*((300.0/ts(i))**1.3)*uh2o(i,k)+6.5e-6*tmpr*
     * (uh2o(i,k)**2.0)+2.09e-5*(pres(i,k)**0.626)*uco2(i,k)
      
       ac11h2o=ac11h2ol(1)
      if (uh2os(i,k).gt.u11h2o(8)) ac11h2o=ac11h2ol(8)       
      do l=1,7      
      if (uh2os(i,k).gt.u11h2o(l).and.uh2os(i,k).le.u11h2o(l+1)) then 
      ac11h2o=ac11h2ol(l)+(ac11h2ol(l)-ac11h2ol(l+1))*(uh2os(i,k)-
     * u11h2o(l))/(u11h2o(l)-u11h2o(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac11h2o='
!      write (*,30) ac11h2o      
      
!      
      acir(i,k,15)= ac11h2o*((300.0/ts(i))**2.1)*uh2o(i,k)+2.5e-5*tmpr*
     * (uh2o(i,k)**2.0)+2.60e-5*(pres(i,k)**0.606)*uco2(i,k)    
 
 
        ac12h2o=ac12h2ol(1)
      if (uh2os(i,k).gt.u12h2o(12)) ac12h2o=ac12h2ol(12)       
      do l=1,11      
      if (uh2os(i,k).gt.u12h2o(l).and.uh2os(i,k).le.u12h2o(l+1)) then 
      ac12h2o=ac12h2ol(l)+(ac12h2ol(l)-ac12h2ol(l+1))*(uh2os(i,k)-
     * u12h2o(l))/(u12h2o(l)-u12h2o(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac12h2o='
!      write (*,30) ac12h2o      
      
!      
      acir(i,k,16)= ac12h2o*((300.0/ts(i))**1.2)*uh2o(i,k)+8.0e-5*tmpr*
     * (uh2o(i,k)**2.0)+6.85e-6*(pres(i,k)**0.518)*uco2(i,k)    
! 
      acir(i,k,17)= 1.3e-3*(1.0-0.00235*(300.0-temp(i,k)))*
     * (pres(i,k)**0.65)*uh2o(i,k)+6.0e-6*tmpr*(uh2o(i,k)**2.0)+
     * 1.99e-5*(pres(i,k)**0.782)*uco2(i,k)     
!

        ac14h2o=ac14h2ol(1)
      if (uh2os(i,k).gt.u14h2o(11)) ac14h2o=ac14h2ol(11)       
      do l=1,10      
      if (uh2os(i,k).gt.u14h2o(l).and.uh2os(i,k).le.u14h2o(l+1)) then 
      ac14h2o=ac14h2ol(l)+(ac14h2ol(l)-ac14h2ol(l+1))*(uh2os(i,k)-
     * u14h2o(l))/(u14h2o(l)-u14h2o(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac14h2o='
!      write (*,30) ac14h2o      
      
!      
      acir(i,k,18)= ac14h2o*uh2o(i,k)+2.6e-5*tmpr*
     * (uh2o(i,k)**2.0)+9.49e-6*(pres(i,k)**0.609)*uco2(i,k)                
!
        ac15h2o=ac15h2ol(1)
      if (uh2os(i,k).gt.u15h2o(10)) ac15h2o=ac15h2ol(10)       
      do l=1,9      
      if (uh2os(i,k).gt.u15h2o(l).and.uh2os(i,k).le.u15h2o(l+1)) then 
      ac15h2o=ac15h2ol(l)+(ac15h2ol(l)-ac15h2ol(l+1))*(uh2os(i,k)-
     * u15h2o(l))/(u15h2o(l)-u15h2o(l+1)) 
      endif      
      enddo
      
!      write (*,*) 'ac15h2o='
!      write (*,30) ac15h2o      
      
!      
      acir(i,k,19)= ac15h2o*uh2o(i,k)+2.0e-4*tmpr1*
     * (uh2o(i,k)**2.0)+3.64e-5*(pres(i,k)**0.799)*uco2(i,k)                
!        
       acir(i,k,20)=0.
!           
      enddo
      enddo
      
  30   format (1E12.4)
   
      return
      end
      
            

