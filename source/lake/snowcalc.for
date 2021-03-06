      SUBROUTINE snow_calc(t2,Erad_sc,ts,snow_sc,iyear,imonth,
     & iday,CCT,pmelt,pheat,dt)
      use phys_constants
      use phys_parameters
      use numeric_params
      use atmos
      use driving_params
      use arrays
      use phys_func, only:
     & EXTINCT_SNOW
      implicit none

CCCCCCCCCCCCCCC input parameters:
CC      e, cal/(cm2*s), --- surface heat balance: e = Radiation-SIGMA*(TS**4)-Latent-Sensible-InSoil
CC      Elatent, cal/(cm2*s), --- latent heat flux.
CC      ts, deg C, --- surface temperature.
CC      precip, cm --- precipitation during dt.
CC
CCCCCCCCC  some other parameters: CCCCCCCCCCC
CC      densny, g/cm3       --- density of fresh snow, depends on temperature at 2m
CC      dznorm, cm          --- initial thickness of snow layers
CC      dzmin, cm           --- minimal thickness of snow layers
CC      dt, sec             --- timestep
CC      ms                  --- max. number of layers in snowpack
CC      cwat, cal/(gr K)    --- the water specific heat content
CC      csnow, cal/(gr K)   --- the snow specific heat content
CC      rhow, gr/cm3        --- the liquid water density
CC      PLI,  cal/g         --- latent heat for freezing/melting
CC      PL,  cal/g          --- latent heat for evapor./condens.
CC      whc                 --- hydraulic constant
CC      hcond, cm/sec       --- hydraulic conductivity
CC
CCCCCCCCCCC  values to be calculate after every timestep:
CC      wl(i), cm       --- moisture content in snow layer No.i
CC      t(i), deg.C     --- mean temperature of snow layer No.i
CC      dz(i), cm       --- hickness of snow layer No.i
CC      dens(i), g/cm3  --- mean density of snow layer No.i
CC      hsnow, cm       --- snow depth
CC      swe, cm         --- water equivalent snow depth
CC      snmelt, cm/sec  --- mean snow melting speed during dt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      INTEGER(4) j,i,k,itop,iyear,imonth,iday
       
      integer(4) nmonth,nyear,nday,nhour
      
      common /time/ nyear,nmonth,nday,nhour
      COMMON /SOILDAT/ dz(ms),itop
      COMMON /BL/       ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     &                  ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     &                  HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     &                  ElatOld,HFold,PRSold,extinct(ms)
      COMMON /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     &                  ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),
     &                  dens(ms)
      COMMON /spin/ Cond,Phase

      common /watericesnowarr/ lams, q
      real(8) dt
      real(8) dz,q(110)
      real(8) lams(ms)
      real(8) ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD, 
     &ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,HS,ES,TGRNEW,QGRNEW,
     &WSNEW,SNNEW,RUNOFF,ElatOld,HFold,PRSold
      real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
      real(8) extinct,Erad_sc,ts,snow_sc,
     &evap,erain,countrad,qbase,dmass,densev,counter,qin,rf,
     &smelti,ein,dzi,ziz,porosityOfSnow,w0,fukt,dztop,smelt,eout,p1,q0,
     &dzmasi,swe,densny,Erad1,eflux1,Phase,Ei,sigma,Cond,Elatent1
      real(8) EInRad,T2,rho_dry,CCT(ML),a2,eta,P,dens_old,ti(ml)
      real(8) pheat,pmelt

      SAVE
       
      precip = precip*dt
      Elatent1 = Elatent

!     hsnow = 0.
!     do i = itop, ms-1
!       hsnow = hsnow + dz(i)
!     end do
!     eFlux = eFlux - Erad

      eflux1 = eflux !/2 
      Erad1 =  Erad_sc !/2

!     eflux1 = eflux1*(1.-exp(-0.02*snow))  ! not -0.02*snow*10.                  
      do i=itop,ms-1
!       wl(i) = wl(i)*dz(i)*dens(i)/rhow    
        ti(i) = t(i)
      end do
      
      snmelt = 0.
      swe = 0.0
      Phase = 0.
c     if(t2 .le. -10) then
c         densny = 0.05
c     elseif(t2 .le. -2) then
c         densny = 5./800.*(t2+18.)
c         densny = 0.1
c     else
c         densny = (t2+4.)/20.
c     end if
      densny = (67.9 + 51.3 * dexp(t2/2.6))
      if(ts .gt. 0.00001) then
          Elatent1 = 0
      elseif(snow_sc .lt. 0.000001 .and. Elatent .gt. 0.0) then    
          Elatent1 = 0
      end if

      dz(ms-1) = dz(ms-1) + dz(ms)
      evap = Elatent1/(PLI+PL) 
!     totalevap=totalevap+evap/rhow*dt
      hsnow=0

!     extinct(itop)=0
      do i = itop,ms-1
!       extinct(i) = (0.9 - 0.4/0.45*(dens(i) - 0.05))*100.
!       extinct(i) = dexp(-(0.13*dens(i)+3.4))
        extinct(i) = EXTINCT_SNOW(dens(i))
!       extinct(i) = 0.
      end do

      if(itop.eq.ms) then

      else
CC        --- first top snow layer ---
          i=itop
          dztop=dz(i)
          Erain = pheat                           
          if(extinct(itop).eq.0.0) then
              Eflux1 = Eflux1 + Erad1 + Erain                            
          else    
              if(itop.eq.ms-1) then
                 Eflux1 = Eflux1 + Erain + Erad1            
              else
                 Eflux1 = Eflux1 + Erain + Erad1*(1 - extinct(i)**dz(i))      
              end if
              countrad = Erad1*(1 - extinct(i)**dz(i))
          end if               
CC rain   heat from rain is added as pheat
CC rain   the rain water itself is assumed to freeze giving heat to top layer
CC rain   then this heat will melt snow
          if(Eflux1.gt.0.0)then
              if(itop.eq.ms-1) then
                  if(Ts.ge.-0.0001)then
                    if(Eflux1/PLI*dt .le. DENS(i)*DZ(i)-WL(i)*rhow) then
                          smelt = Eflux1/PLI/rhow                          
                          Eout=0.
                      else
                          smelt = (DZ(i)*DENS(i)/rhow - WL(i))/dt
                          Eout = Eflux1 - smelt*rhow*PLI    
                      end if
                  else
                      smelt = 0.0
                      Eout = 0.
                  end if
              else
                  smelt = 0.0
                  Eout = Eflux1
              end if
          else
              smelt = 0.
              Eout = 0.
          end if

      
Cc
Cc    PERCOLATION
Cc
          smelt = smelt + wl(i)/dt
          wl(i) = 0.
          qbase = smelt

CC
CC    TOP  LAYER THICKNESS AND DENSITY
CC
          dmass=dens(i)*dztop-(evap+smelt*rhow)*dt     
          if(evap.gt.0.0)then
              densev=dens(i)
          else
              densev=rhoi
          endif
          dztop=dztop-(smelt*rhow/dens(i)+evap/densev)*dt   
        if(dztop.eq.0.) dztop = 0.000001  
          dens(i)= dmax1(dmass/dztop,0.d0)
          if(dens(i).gt.rhoi) then
              dztop = dztop*dens(i)/rhoi
              dens(i)=rhoi
              dmass=dens(i)*dztop
          endif
          dz(i)=dztop
          if(itop .eq.ms-1) snmelt=qbase
      end if 


      counter = dz(itop)
      if(itop.lt.ms-1) then
CC
CC INTERNAL LAYERS
CC
          do i=itop+1,ms-1
              qin = qbase                                     
              rf = 0.0                                        
              smelti = 0.0                                    
              Ein=Eout                                        
              T(i) = T(i)*dz(i)*dens(i)/(qbase*rhow*dt+dz(i)*dens(i))  
              dzi=dz(i)                                       
              wl(i) = wl(i) + qbase*dt
              if(wl(i).lt.0.) wl(i) = 0.
              w0=wl(i)
              if(extinct(i).eq.0.0) then
                  Ei = Ein                
              else
                  if(i.eq.ms-1) then
                      Ei=Ein + (Erad1 - countrad)
                  else
                      Ei=Ein + Erad1*(extinct(i)**(counter) -
     &                    extinct(i)**(counter+dz(i)) )
                      EInRad = Erad1*(extinct(i)**(counter) -
     &                    extinct(i)**(counter+dz(i)) )
                  end if                            
                  countrad = countrad + Erad1*(extinct(i)**(counter)-
     &              extinct(i)**(counter+dz(i)) )
                  counter = counter+dz(i)
              end if

              dz(i)=dz(i)+qbase*dt  
              ziz=dzi/dz(i)
              dens(i)=ziz*dens(i)+(1.0-ziz)*rhow

            
              if(Ei.gt.0.0) then
                  if(T(i).ge.0.0001) then  ! 0.0001
                      Ei = Ei + csnow*dens(i)*T(i)*dz(i)/dt
                      T(i) = 0.0
                      if(Ei/PLI*dt .le. (DZ(i)*DENS(i)-WL(i)*rhow)) then
                          smelti = Ei/PLI/rhow 
                          Eout = 0.
                      else
                          smelti = (DZ(i)*DENS(i)/rhow - WL(i))/dt
                          Eout = Ei - smelti*rhow*PLI 
                          dzi = 0.
                      end if
                      wl(i) = wl(i) + smelti*dt
                  else
Cc        T<0
                      T(i) = T(i) + Ei*dt/(csnow*dens(i)*dz(i))
                      if(T(i).lt.-0.0001) then   !0.0001
                      if(wl(i) .gt. -T(i)*csnow*dz(i)*dens(i)/PLI/rhow) 
     &                            then
                            rf = -T(i)*csnow*dz(i)*dens(i)/PLI/dt/rhow 
                              T(i) = 0.0
                              wl(i) = wl(i) - rf*dt
                          else
                              rf = wl(i)/dt
                     T(i) = T(i) + rf*rhow*dt*PLI/(csnow*dz(i)*dens(i))    
                              wl(i) = 0.
                          end if
                          Eout = 0.
                      end if
                      if(T(i).ge.0.0001)then  !0.0001
                          if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &                        dz(i)*dens(i)) then
                         smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow 
                              Eout = 0.
                          else
                              smelti = dz(i)*dens(i)/dt/rhow  
                              Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                              - smelti*PLI*rhow 
                              dzi = 0.
                          end if
                          T(i) = 0.0
                          wl(i) = wl(i) + smelti*dt
                      endif
                  endif
              else

Cc    E = 0 (E<0 is impossible case)
                  if(T(i).ge.0.0001) then     !0.0001
                      if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &                        dz(i)*dens(i)) then
                          smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow  
                          T(i) = 0.0
                          Eout = 0.
                          wl(i) = wl(i) + smelti*dt
                      else
                          smelti = dz(i)*dens(i)/dt/rhow  
                          Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                            - smelti*PLI*rhow  
                          T(i) = 0.0
                          dzi = 0.
                          wl(i) = wl(i) + smelti*dt
                      end if
                  else
Cc    T<0
                  if(wl(i) .gt. -csnow*dens(i)*T(i)*dz(i)/PLI/rhow) then  
                          rf = -csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow   
                          T(i) = 0.0
                          wl(i) = wl(i) - rf*dt
                      else
                          rf = wl(i)/dt
                          wl(i) = 0.
                      T(i) = T(i) + rf*rhow*dt*PLI/(csnow*dz(i)*dens(i))
                      end if
                      Eout = 0.
                      if(T(i).ge.0.0001) then   !0.0001
                          if(csnow*dens(i)*T(i)*dz(i)/PLI .le.
     &                         dz(i)*dens(i)) then
                          smelti = csnow*dens(i)*T(i)*dz(i)/dt/PLI/rhow 
                              Eout = 0.0
                          else
                              smelti = dz(i)*dens(i)/dt/rhow
                              Eout = csnow*dens(i)*T(i)*dz(i)/dt
     &                            - smelti*PLI*rhow 
                              dzi = 0.
                          end if
                          T(i) = 0.0
                          wl(i) = wl(i) + smelti*dt
                      endif
                  endif
              endif
      
Cc
Cc
Cc    internal q
Cc
              rho_dry = (dens(i)-w0/(dz(i))*RHOW)/
     *            (1-w0/(dz(i)))
              if(rho_dry.lt.0.) WRITE(*,*) 'rho_dry = ', rho_dry

!     PorosityOfSnow = 1 - rho_dry/rhoi - wl(i)/dz(i)*(1-rhow/rhoi) wl(i)/rhow/dz(i)
             PorosityOfSnow = 1 - rho_dry/rhoi - wl(i)/dz(i)
              PorosityOfSnow = dmax1(PorosityOfSnow,whc + 0.1)
              p1 = PorosityOfSnow - whc   ! whc = water holding capacity
              fukt = wl(i)/dz(i)
              if(fukt.gt.whc)then
                  fukt = (fukt-whc)/p1
                  q0 = hcond*fukt**3.0
      
                  qbase = dmin1(q0,wl(i)/dt)
                  wl(i) = wl(i) - qbase*dt
              else
                  qbase = 0.0
                  q0 = 0.0
              end if
              if(dzi .le. 0.000001) then   
                  qbase = qbase + wl(i)/dt
                  dz(i)=dzi
              else
                  dzmasi=dz(i)*dens(i)-qbase*dt*rhow   
                  dzi=dz(i)+rhow*(-qbase/RHOW-rf/RHOW+rf/rhoi-  
     &                    smelti/rho_dry+smelti/RHOW)*dt     
                  dens(i)=dmax1(dzmasi/dzi,0.d0)
                  dz(i)=dzi
                  if(dens(i).gt.rhoi)then
                      dz(i)=dz(i)*dens(i)/rhoi
                      dens(i)=rhoi
                  end if
              end if
              Phase = Phase                                 !diagnostic variable
     &                + (t(i)-ti(i))*dens(i)*dz(i)/snow_sc/dt
      
              if(i.eq.ms-1) snmelt=qbase
          end do
      end if
      if (itop==ms-1.and.dz(itop)<0.000001) then
        dz(itop)=0
        dens(itop)=0
        goto 123
      endif
      if (itop .lt. ms .and. dz(itop) .le. 0.000001) then 
        dz(itop)=0.
        dens(itop)=0.
        itop=itop+1
      end if
          
222   do i=ms-1,itop,-1
          if(dz(i).le.0.000001) then   
              do k=i-1,itop,-1
                  if(dz(k).gt.0.000001) goto 111   
              end do
      
111           do j=i,itop,-1
                  dz(j)=dz(k+j-i)
                  dens(j)=dens(k+j-i)
                  t(j)=t(k+j-i)
                  wl(j)=wl(k+j-i)
              end do
              itop=itop-k+i
              goto 222
          end if
      end do
123   continue
       
333   do i=itop,ms-1,1
          if(dz(i).lt.dzmin.and.itop.lt.ms-1) then
          if(i.lt.ms-1) then
              dens(i+1)=dens(i+1)*dz(i+1)/(dz(i)+dz(i+1))+
     &                dens(i)*dz(i)/(dz(i)+dz(i+1))
              t(i+1)=t(i+1)*dz(i+1)/(dz(i)+dz(i+1))+
     &                t(i)*dz(i)/(dz(i)+dz(i+1))
              wl(i+1)=wl(i+1)+wl(i)
              dz(i+1)=dz(i+1)+dz(i)
              do j=i,itop+1,-1
                  dens(j)=dens(j-1)
                  t(j)=t(j-1)
                  wl(j)=wl(j-1)
                  dz(j)=dz(j-1)
              end do
              itop=itop + 1
              goto 333
          end if
          end if
      end do
      
      swe = 0.
      hsnow = 0.
      do j=ms-1,itop,-1
          swe = swe + dz(j)*dens(j)
          hsnow = hsnow + dz(j)
      end do
       
cccccccccccccc snow densification due to gravity and metamorphism cccccccccccccc
      a2 = 0.00000066
      sigma = 75.                             ! metamorphic processes, Pa
      do j = itop+1,ms-1
          if(dens(j).lt.3*densny) then
              dens_old = dens(j)
              eta =                           ! compactive viscosity of snow 
     &            a2*dexp(19.3*dens(j)/rhoi)*
     &            dexp(67300/8.314/(t(j)+273.15))
              P = 0.0                         ! gravity, Pa
              do i = j,itop,-1
                  !P = P + dens(i)/1000.*grav*(dz(i))*10000.
                  P = P + dens(i)*grav*dz(i)
              end do
              dens(j) = dens(j) + 0.001*(dt*dens(j)*1000.*(sigma+P)/eta)
                          dens(j) = dmin1(dens(j),rhoi)
              dz(j) = dz(j) * dens_old/dens(j)
          end if
      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      swe = 0.
      hsnow = 0.
      do j=ms-1,itop,-1
          swe = swe + dz(j)*dens(j)
          hsnow = hsnow + dz(j)
      end do
!     do m = itop,ms-1
!       por(m) = 1 - dens(m)/rhoi
!     end do
      
      snow_sc = swe

        !! dz(ms) !!
      if(itop.ne.ms-1) then
          dz(ms) = dz(ms-1)/3.
          dz(ms-1) = dz(ms-1) - dz(ms)
      else
          dz(ms) = dz(ms-1)/2.
          dz(ms-1) = dz(ms-1) - dz(ms)
      end if
      if(itop.eq.ms-1) then
          z(ms-1) = z(ms) - dz(ms-1) - dz(ms)
      else
          z(ms-1) = z(ms) - dz(ms-1)/2. - dz(ms)
      end if
      DO j = MS-2,itop+1,-1
          Z(j) = z(j+1)-dz(j)/2.-dz(j+1)/2.
      END DO
      if(itop.lt.ms-1) z(itop) = z(itop+1) - dz(itop+1)/2. - dz(itop)

!     do i=itop,ms-1
!       wl(i) = wl(i)/dz(i)/dens(i)*rhow
!     end do


      
      do i = 1, itop-1
          dz(i) = 0.
      end do

      
      precip = precip/dt
      !eflux = eflux+Erad
      
      return
      END


      SUBROUTINE addPrecipToSnow(T2,ts,snow,iyear,imonth,
     &    iday,CCT,pmelt,pheat,dt)
      use phys_constants    
      use phys_parameters
      use numeric_params
      use driving_params
      use atmos
      use arrays, only:
     & h1,l1,hs1,ls1
      
      implicit none
CCCCCCCCCCCCCCC input parameters:
CC      e, cal/(cm2*s), --- surface heat balance: e = Radiation-SIGMA*(TS**4)-Latent-Sensible-InSoil
CC      Elatent, cal/(cm2*s), --- latent heat flux.
CC      ts, deg C, --- surface temperature.
CC      precip, cm --- precipitation during dt.
CC
CCCCCCCCC  some other parameters: CCCCCCCCCCC
CC      densny, g/cm3       --- density of fresh snow, depends on surface temperature
CC      dznorm, cm          --- initial thickness of snow layers
CC      dzmin, cm           --- minimal thickness of snow layers
CC      dt, sec             --- timestep
CC      ms                  --- max. number of layers in snowpack
CC      cwat, cal/(gr K)    --- the water specific heat content
CC      csnow, cal/(gr K)   --- the snow specific heat content
CC      rhow, gr/cm3        --- the liquid water density
CC      PLI,  cal/g          --- latent heat for freezing/melting
CC      PL,  cal/g          --- latent heat for evapor./condens.
CC      whc                 --- hydraulic constant
CC      hcond, cm/sec       --- hydraulic conductivity
CC
CCCCCCCCCCC  values to be calculate after every timestep:
CC      wl(i), cm       --- moisture content in snow layer No.i
CC      t(i), deg.C     --- mean temperature of snow layer No.i
CC      dz(i), cm       --- thickness of snow layer No.i
CC      dens(i), g/cm3  --- mean density of snow layer No.i
CC      hsnow, cm       --- snow depth
CC      swe, cm         --- water equivalent snow depth
CC      snmelt, cm/sec  --- mean snow melting speed during dt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      COMMON /SOILDAT/ dz(ms),itop
      COMMON /BL/       ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     &                  ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,
     &                  HS,ES,TGRNEW,QGRNEW,WSNEW,SNNEW,RUNOFF,
     &                  ElatOld,HFold,PRSold,extinct(ms)
      COMMON /SOILSOL/ AL(ML),DLT(ML),DVT(ML),ALLL(ML),DL(ML),
     &                  ALV(ML),DV(ML),Z(ML),T(ML),WL(ML),WV(ML),WI(ML),
     &                  dens(ms)
      COMMON /spin/ Cond,Phase
      real(8) dz
      real(8) ST,PGR,TGROLD,QGROLD,RADIAT,WSOLD,SNOLD,
     &ZS,thsoil,whsoil,BOLD,RF1,RF2,SNMELT,HSNOW,HS,ES,TGRNEW,QGRNEW,
     &WSNEW,SNNEW,RUNOFF,ElatOld,HFold,PRSold
      real(8) AL,DLT,DVT,ALLL,DL,ALV,DV,Z,T,WL,WV,WI,dens
      real(8) extinct,ts,snow,dz0,pmelt,pheat,
     & z0z,dmass, swe,densny,Cond,Phase
      INTEGER(4) j,i,k,mn,itop,iyear,imonth,iday
      real(8) T2,CCT(ML)
      real(8) dt 

      SAVE
      
      precip = precip*dt
          
            
!     if (snow > 0.) then
!       write (*, *) snow
!       STOP
!     end if
!     T2 = T2 - 273.15

c     if(T2 .le. -10) then
c         densny = 0.05
c     elseif(ts .le. -2) then
c         densny = 5./800.*(t2+18.)
cc        densny = 0.1
c     else
c         densny = (t2+4.)/20.
c     end if

      densny = (67.9 + 51.3 * dexp(t2/2.6))
      if(itop.lt.ms)  dz(ms-1) = dz(ms-1) + dz(ms) 
      
      snow = 0.
      do mn = itop,ms-1
          snow = snow + dz(mn)*dens(mn)
      end do
      
      
      j=itop
555   continue
      if(dz(itop).lt.dznorm+dzmin) goto 444
      dz(j-1) = dz(j) - dznorm
      dz(j) = dznorm
      z(j-1)= z(itop) 
      z(j)= z(j-1) + dz(j-1)*2.
      wl(j-1)=0.0
      dens(j-1)=dens(j)
      t(j-1)=t(j)
      j=j-1
      itop=j
      goto 555
444   continue
      
      if(T2.gt.0.) then
          pmelt = precip/dt                       
!         pmelt = 0.                       
          pheat = T2*precip*cwat*rhow/dt      
      else       
          pmelt=0.0
          pheat=0.0
          if(itop.eq.ms .and. precip.gt.0.0) then
              itop=itop-1
              z(itop) = z(itop+1)
              t(itop) = t(itop+1)
              dens(itop) = densny
              wl(itop) = 0
          end if
          if(itop.ne.ms) then
              dz0 = dz(itop)
              z(itop) = z(itop) - precip*rhow/densny   
              dz(itop) = dz(itop) + precip*rhow/densny   
              if(dz(itop) .lt. dznorm+dzmin) then
                  z0z=dz0/dz(itop)
                  dens(itop)=z0z*dens(itop)+(1.0-z0z)*densny
              end if
              j=itop
55            continue
              if(dz(itop) .lt. dznorm+dzmin) goto 54
              dz(j-1) = dz(j) - dznorm
              dz(j) = dznorm
              z(j-1)= z(itop) 
              z(j)= z(itop) - dz(itop)*0.5
              z0z=dz0/(dz(itop)+dz(itop-1))
              dens(itop)=z0z*dens(itop)+(1.0-z0z)*densny
              wl(j-1)=0.0
              dz0 = 0.0
              dens(j-1)=densny
              t(j-1)=t(j)
              j=j-1
              itop=j
              goto 55
54            continue
          end if
      end if
          
      if(pheat.gt.0.0) then
          if(itop.eq.ms-1) then
              dmass = dens(itop)*dz(itop) + precip*rhow    
              dz(itop) = dz(itop) + precip*rhow/rhoi   
              if(dz(itop).eq.0.) dz(itop) = 0.000001   
              dens(itop)= dmax1(dmass/dz(itop),0.d0)
              if(dens(itop).gt.rhoi) then
                  dz(itop) = dz(itop)*dens(itop)/rhoi
                  dens(itop)=rhoi
              end if
          else
              if(itop.lt.ms-1) then
                  mn = itop+1
                  dmass = dens(mn)*dz(mn) + precip*rhow  
                  dz(mn) = dz(mn) + precip*rhow/rhoi 
                  if(dz(mn).eq.0.) dz(mn) = 0.000001    
                  dens(mn)= dmax1(dmass/dz(mn),0.d0)
                  if(dens(mn).gt.rhoi) then
                      dz(mn) = dz(mn)*dens(mn)/rhoi
                      dens(mn)=rhoi
                  end if
                  wl(mn) = wl(mn) + precip
          T(mn) = T(mn)*dz(mn)*dens(mn)/(precip*rhow+dz(mn)*dens(mn)) 
              end if
          end if
      end if

        !!  dz(ms) !!!
      if(itop.ne.ms-1) then
          dz(ms) = dz(ms-1)/3.
          dz(ms-1) = dz(ms-1) - dz(ms)
      else
          dz(ms) = dz(ms-1)/2.
          dz(ms-1) = dz(ms-1) - dz(ms)
      end if
      if(itop.eq.ms-1) then
          z(ms-1) = z(ms) - dz(ms-1) - dz(ms)
      else 
          z(ms-1) = z(ms) - dz(ms-1)/2. - dz(ms)
      end if
      DO Mn = MS-2,itop+1,-1
          Z(Mn) = z(mn+1)-dz(mn)/2.-dz(mn+1)/2.
      END DO
      if(itop.lt.ms-1) z(itop) = z(itop+1) - dz(itop+1)/2. - dz(itop)

      !if (dz(itop).lt.dznorm+dzmin) goto 777 
      
      swe = 0.0
      do mn=ms-1,itop,-1
          swe=swe+dz(mn)*dens(mn)
      end do
      swe = swe + dz(ms)*dens(ms-1)
      snow = swe
      
      hsnow = 0.
      do i = itop, ms
        hsnow = hsnow + dz(i)
      end do
             
      precip = precip/dt
      
      return
      END
