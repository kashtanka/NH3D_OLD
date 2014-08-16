      SUBROUTINE RichPBL(h,le,tau,Tsurf1,humsurf,z0)
      
!     RichPBL calculates heat fluxes according to Louis parameterization (Louis, 1079)    
      
      use atmos, only:
     & tempair,
     & humair,
     & pressure,
     & zref
      implicit none
      
      real(8),parameter:: gg = 9.80665
      real(8),parameter:: ra = 287.05
      real(8),parameter:: rv = 461.51
      real(8),parameter:: cp = 1005.46
      real(8),parameter:: cpv = 1869.46
      real(8),parameter:: xl = 2500800.
      real(8),parameter:: stefan = 5.6697d-08
      real(8),parameter:: xpi = 3.141592654
      real(8),parameter:: too = 273.16
      real(8),parameter:: xpoo = 100000.
      real(8),parameter:: vkarmn = .4
      real(8),parameter:: xli = 2834500.
      real(8),parameter:: tcdil = 86400.
      real(8),parameter:: rascp = ra / cp
      real(8),parameter:: xlsg = xl / gg
      real(8),parameter:: xlscp = xl / cp
      real(8),parameter:: cpsl = cp / xl
      real(8),parameter:: cpsg = cp / gg
      real(8),parameter:: gscp = gg / cp
      real(8),parameter:: unsg = 1. / gg
      real(8),parameter:: gsra = gg / ra
      real(8),parameter:: rasl = ra / xl
      real(8),parameter:: rascp2 = ra / (2.*cp)
      real(8),parameter:: rasrv = ra / rv
      real(8),parameter:: rvsra = rv / ra
      real(8),parameter:: etv = rvsra - 1.
      real(8),parameter:: ecph = cpv / cp - 1.
      real(8),parameter:: xlsrv = xl / rv
      real(8),parameter:: unscp = 1. / cp
      real(8),parameter:: cpvmcp = cpv - cp
      real(8),parameter:: etvq = 1. - rasrv
      real(8),parameter:: xlf = xli - xl
      real(8),parameter:: xliscp = xli / cp
      real(8),parameter:: xlisg = xli / gg

      real(8) z0hz0,Tsurf1,humsurf,z,xmu,xfh,ts,h,le,leg,
     & urel,vrel,ta,qa,ps,ua,z0h,z0,hu,veg,ztvi,rhoa,zdepl,zua,ztvis,Ri,
     & valnon,zeps,zch,zsta,iyra,zdi,rra,zdim,cdm,zds,chstar2,qs,ph2,
     & cdh,cmstar2,pm2,zdsm,xhu,tau,wind10,u,v
      common /wind/ wind10,urel,vrel,u,v

      SAVE
      
      ta=tempair+273.15
      qa=humair
      ps=pressure
      ua=dsqrt(urel**2+vrel**2)
      ts=Tsurf1+273.15
      qs=humsurf

      z0hz0=0.1
      z0h=z0*z0hz0

c
c**********************************************************************
c           delta calculatlon
c
!      delta=0.
!      if (veg.gt.0.) delta=(wr(li)/wrmax)**(2./3.)
c
c**********************************************************************
c                         hv calculation
c                         ______________
c
      ! hv = 1. - dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs *
      !&      (1. - delta) / (rra + rs )

        hu=1. !air above water is saturated
        veg=0. !no vegetation on the lake
c
c?    rra ‚ a resistencia aerodinamica = 1/ChVa  -volta
c        a ser calculada … frente (correctamente

c        Implicit resolution of the surface temperature equation
c        linearization of Ts

        ztvi = ta * (1. + etv * qa)
        rhoa = ps / (ra * ztvi)
        !qs = (hu * (1. - veg) + veg * hv) * qsat (ts(li),ps) +
     &  !     (1. - hv) * veg * qa
c
c**********************************************************************
c                CMWF, 1981/11/25-27, pp. 59-79)
c**********************************************************************
c
        zdepl = 0.
        !write(*,*) 'nhli9904zref',zref
        zua=zref + zdepl
        ztvis = ts * (1. + etv * qs)
c______________________________________________________________________
c       calculo do numero de richardson e da resistencia aerodinamica
c
      !write(*,*) 'zref',zref
      ua=dmax1(1.d0,ua)
      ri=2.*gg*zua*zua*(ztvi-ztvis+gscp*zref)/(ztvi+ztvis)/(ua*ua)/zref
c
c       calculo do comprimento de Monin -Obukhov,de ustar,tstar
c
        valnon=-999.
c        call calclmon(z,zO,zOh,ua,ztvi,ztvls,xlmon,ustar,tstar,valnon)
        zeps= 1.e-6
        !write(* ,*) 'zua,z0',zua,z0
        zch=(vkarmn/dlog(dabs(zua/z0)))**2
        !write(*,*) 'zch',zch
        zsta=ri*ua*ua
c
c  iyra escolhe entre a nova e a velha versao da rotina yra
c         de NP89. Quando iyra (ficheiro de entrada(?!)) e maior
c         do que 0 utiliza-se a nova versao que destingue entre
c         os comprimentos de rugosidade
c
        iyra=1
        if(iyra.eq.0) then
c
          if(ri.lt.0.) then
           zdi=1./(ua+75.*zch*dsqrt(-zref*zsta/z0h))
           rra=zch*(ua-15.*zsta*zdi)
           zdim=1./(1.+75.*zch*dsqrt(-zref*ri/z0))
           cdm=zch*(1.-10.*ri*zdim)*ua
          else
           zds=dsqrt(ua * ua + 5. * zsta + zeps)
           rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)
           zdim=(1.+5.*ri)**2.
           cdm=zch*ua/zdim
          endif
c
        else
c
c nova rotina yramng
          ! write(*,*) 'estez0h',z0h
          !write(*,*) 'z,zref',z,zref
           z=zref
           xmu=dlog(z0/z0h)
           xfh=dlog(z/z0)/dlog(z/z0h)
          ! write(*,*) '--------------------'
          ! write(*,*) 'z0,z0h',z0,z0h
           if(ri.lt.0.) then
                zdi=1./(ua+chstar2(xmu)*zch*15.*(z/z0h)**
     &          ph2(xmu)*xfh
     s              *dsqrt(-zsta))
                rra=zch*(ua-15.*zsta*zdi)*xfh
c
c alteracoes para o calculo do coeficiente de dragg_______________
c        cdm=Cd*uvsuf
c
                cdh=rra
                zdim=1./(ua+cmstar2(xmu)*zch*10.*(z/z0)**
     &          pm2(xmu)
     s              *dsqrt(-zsta))
                cdm=zch*(ua-10.*zsta*zdim)
c_______________________________________________________________________
           else
                zds = dsqrt(ua * ua + 5. * zsta + zeps)
                rra = zch*ua/(1.+15.*zsta*zds/ua/ua/ua)*xfh
c
                cdh=rra
                zdsm= dsqrt(ua * ua + 5. * zsta + zeps)
                cdm = zch*ua/(1.+10.*zsta/zdsm/ua)
      
           endif
c
        endif
c
        rra=1./rra
c----------------------------------------------------------------------
c
!      zrsra = rhoa / rra
!      if (iclay.ge.0) then
!         ct = 1. / ((1. - veg) / cso + veg / cv)
!         hv=1.-dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa)) * rs*
!     &    (1. - delta) / (rra + rs)
c
!         za=1. / dt + ct * (4. * emis * stefan * (ts(li)**3) +
!     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu)
!     &      +zrsra*cp)+2.*xpi/tcdil
!         zb=1./dt+ ct*(3. * emis * stefan * (ts(li)** 3) +
!     &      zrsra*xl*dqsat(ts(li),ps)*(veg*hv+(1.-veg)*hu*xhu))
!         zc=2.*xpi*t2(li)/tcdil+ct*(zrsra*cp*ta+rg*
!     &      (1.-alb)+emis*rat-zrsra*xl*(veg*hv*(qsat(ts(li),
!     &      ps)-qa)+(1.-veg)*xhu*(hu*qsat(ts(li),ps)-qa)))
c
!         ts(lf) = (ts(li) * zb + zc) / za
c
c       resolution of t2 equation
c
!         t2(lf) = (t2(li) + dt * ts(lf) / tcdil) / (1. + dt / tcdil)
!      else
!          ts(lf)=tsfunc
!          t2(lf)=t2(li)
!      endif
c
c**********************************************************************
c
      xhu=1.
      !hqs2=hu*qsat(ts(lf),ps)
      !rn = rg * (1. - alb) + emis * (rat - stefan * (ts(lf)**4))
      h = rhoa * cp * (ts - ta) / rra
      leg = rhoa * xl*(1.-veg) *xhu*(qs - qa)/ rra 
      !lev = rhoa * xl * veg * hv * (qsat(ts(lf), ps) - qa) / rra
      !zzhv = dmax1(0.d0,dsign(1.d0,qsat(ts(li), ps) - qa))
      !letr = zzhv * (1. - delta) * rhoa * xl * veg * (qsat(ts(lf), ps)
     &!       - qa) / (rra + rs)
      le = leg !+ lev
      tau = rhoa*cdm*ua
      !hw=h;xlew=le;cdmw=cdm

      return
      END SUBROUTINE RichPBL
          
      function chstar2(xmu)

      implicit real*8(a-h,o-z)
      chstar2=3.2165+4.3431*xmu+.536*xmu*xmu-.0781*xmu*xmu*xmu
      return
      end

      function cmstar2(xmu)

      implicit real*8(a-h,o-z)
      cmstar2=6.8741+2.6933*xmu+.3601*xmu*xmu-.0154*xmu*xmu*xmu
      return
      end

      function ph2(xmu)

      implicit real*8(a-h,o-z)
      ph2    =0.5802-0.1571*xmu+.0327*xmu*xmu-.0026*xmu*xmu*xmu
      return
      end
c
c******************************************************************
      function pm2(xmu)
c****************************************************************
c
      implicit real*8(a-h,o-z)
      pm2    =0.5233-0.0815*xmu+.0135*xmu*xmu-.001*xmu*xmu*xmu
      return
      end
