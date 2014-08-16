      SUBROUTINE PBLDAT

!     PBLDAT assigns values to parameters used in surface layer parameterization
!     (see subr. dragvl and subr., called by dragvl)

      IMPLICIT real (A-H,O-Z)
C
C*=====================================================================
C*       INITIALIZATION OF PBL COMMON BLOCKS (FBL,DRAG)               =
C*    FBL(1)     :  FOR LINEAR EXTRAPOLATION OF WIND ONTO SURFACE     =
C*    FBL(2)     :  ---------------------------------------------     =
C*    FBL(3)     :  PRESCRIBED MINIMAL VALUE OF SURFACE WIND MODULE   =
C*    FBL(4)     :  HEIGHT OF CONSTANT FLUX LAYER ( M )               =
C*    FBL(5)     :  CRITICAL VALUE OF RELATIVE HUMIDITY IN            =
C*               :  CALCULATION OF EQUIVALENT POTENTIAL TEMPERATURE   =
C*    FBL(6)     :  SEA WATER TEMPERATURE UNDER ICE ( DEG. K )        =
C*    FBL(7)     :  PRECRIBED MINIMAL VALUE OF SEA ICE SURFACE        =
C*                  TEMPERATURE ( DEG. K )                            =
C*    FBL(8)     :  RESERVED                                          =
C*    FBL(9)     :  PBL WIND TURNING ANGLE OVER OCEAN ( DEG )         =
C*    FBL(10)    :  --------------------------- ICE ----------------- =
C*    FBL(11)    :  --------------------------- SNOW COVERED SOIL --- =
C*    FBL(12)    :  --------------------------- BARED SOIL ---------- =
C*    FBL(13)    :  --------------------------- IN TROPICS ---------- =
C*    FBL(14)    :  BOUNDARIES OF TROPICS ( RADIANS )                 =
C*    FBL(15)    :  RESERVED
C*    FBL(16)    :  (SNOW COVERED SOIL HEAT CONDUCTIVITY)/DEPTH ----- =
C*    FBL(17)    :  (SEA ICE HEAT CONDUCTIVITY)/(SEA ICE DEPTH) ----- =
C*    FBL(18)    :  (SOIL HEAT CONDUCTIVITY)/DEPTH, CAL/(DEG.*M**2*S) =
C*    FBL(19)    :  BULK SOIL HEAT CAPACITY, CAL/(DEG.*M**2)          =
C*    FBL(20)    :  (SOIL MOISTURE CONDUCTIVITY)/DEPTH,SEC**-1        =
C*=====================================================================
      COMMON /FBL/ FBL(20)
      COMMON /BL1/ AZ0P,VEG,WL,DZZ
      COMMON /HYDR/ SNCR,WLMMX,WSMAX,WSL,WSG,CEFF,
     &              CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF,
     &              ZRM,ZRMM,TRM,FLXMIN,TOMIN,HSNold
      COMMON /VEGSW/ CCQ,CK,SWW,TL
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     &                zeta_st,a_m,a_t,b_m,b_t
      real FBL,AZ0P,VEG,WL,DZZ,SNCR,WLMMX,WSMAX,WSL,WSG,CEFF,
     &CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF,ZRM,ZRMM,TRM,
     &FLXMIN,TOMIN,CCQ,CK,SWW,TL,
     &vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,alpha_m,alpha_t,beta_m,
     &beta_t,gamma_m,gamma_t,zeta_st,a_m,a_t,b_m,b_t

      FBL(1) = 1.
      FBL(2) = 1.
      FBL(3) = 1.
      FBL(4) = 70.
      FBL(5) = 1.
      FBL(6) = 271.5
      FBL(7) = 243.2
      FBL(8) = 0.
      FBL(9) = 20.
      FBL(10) = 10.
      FBL(11) = 30.
      FBL(12) = 30.
      FBL(13) = 0.
      FBL(14) = 20.*3.1415926536/180.
C     FBL(8)...FBL(14) WILL BE NOT USED IN THE MULTILEVEL PBL
      FBL(15) = 0.
      FBL(16) = 0.04/2.
      FBL(17) = 0.5/3.
      FBL(18) = 0.6/2.
      FBL(19) = 4.0E+04
      FBL(20) = 1.0E-6/5.
c------------------------------------------------------------------
c*    subroutine dragvl constants
c------------------------------------------------------------------
      vkc = 0.4
      anu = 0.000015
      z0min = 1.5e-5
      chc = 0.0132
      zeta_st = -2.
      p_m = - 0.25
      p_t = - 0.5
      q_m = - 1./3.
      q_t = - 1./3.
      alpha_m = 1.
      alpha_t = 1.
      beta_m = 4.7
      beta_t = 4.7
      gamma_m = 16.
      gamma_t = 16.
      a_m = (1. + gamma_m*(p_m/q_m - 1.)*zeta_st)
     &      *(1. - gamma_m*zeta_st)**(p_m/q_m - 1.)
      b_m = a_m*gamma_m*(p_m/q_m)/(1. + gamma_m*(p_m/q_m - 1.)*zeta_st)
      a_t = (1. + gamma_t*(p_t/q_t - 1.)*zeta_st)
     &      *(1. - gamma_t*zeta_st)**(p_t/q_t - 1.)
      b_t = a_t*gamma_t*(p_t/q_t)/(1. + gamma_t*(p_t/q_t - 1.)*zeta_st)

      RETURN
      END


      SUBROUTINE dragvl(bx2,bix2,itdrag)
      
!     dragvl calculates exchange coefficients in aerodynamic formulas for surface
!     sensible heat, latent heat and momentum fluxes at the surface
!     following Businger-Dayer interpolation formulas, Beljaars parameterization and others      
      
      real bx(7),bix(11)

      real(8), intent(in) :: bx2(7)
      real(8), intent(out):: bix2(11)
          
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t,grav
      real u_a,tpot_a,tpot_s,q_a,q_s,z,z0m,z0t,tvir_a,
     &tvir_s,dtvir,alam,c_u,c_t,ust,hfl,re,xx,zeta,zeta0m,
     &zeta0t,smo,psi_m,psi_t,phi_m,phi_t
      INTEGER iter,itdrag
     
         
      grav = 9.81
      
      bx=sngl(bx2)  

      u_a = amax1(bx(1),0.1) ! speed 
      tpot_a = bx(2) ! potential temperature - air
      tpot_s = bx(3) ! potential temperature - surface
      q_a = bx(4)
      q_s = bx(5)
      z = bx(6)
      if(bx(7) .gt. 0.) then
          z0m = bx(7)
      else
          z0m = z0min
      end if
      tvir_a = tpot_a * (1. + 0.61 * q_a) ! virtual temp - air
      tvir_s = tpot_s * (1. + 0.61 * q_s) ! virtual temp - surf
      dtvir = tvir_a - tvir_s  ! for smo
      alam = grav/(0.5*(tpot_a + tpot_s))  ! for smo
      c_u = vkc/log(z/z0m) ! for smo
      c_t = c_u
      do iter = 1,itdrag
          ust = c_u * u_a  ! for smo
          hfl =  - c_u * c_t * u_a * dtvir ! for smo Hflux
          if(bx(7) .lt. 0.) z0m = amax1(chc*ust**2/grav, z0min)
          re = ust*z0m/anu
          if(re .le. 0.111) then
              xx = - 2.43
          else
              if(0.111 .lt. re .and. re .le. 16.3) then
                  xx = 0.83*log(re) - 0.6
              else
                  xx = 0.49 * re**0.45
              end if
          end if
          z0t = amax1(z0m*exp(-xx), z0min)
c         z0t = z0m
          if(abs(dtvir) .lt. 1.e-3) then
              zeta = 1.
              zeta0m = 1.
              zeta0t = 1.
          else
              smo = - ust**3/(alam*vkc*hfl) ! vkc = 0.4
              smo = sign(1.,smo)*amax1(abs(smo),1.e-3) !
              zeta = z/smo
              zeta0m = z0m/smo
              zeta0t = z0t/smo
          end if
          c_u = vkc/amax1( (log(z/z0m) - psi_m(zeta,zeta0m)), 0.1)
          c_t = vkc/amax1( (log(z/z0t) - psi_t(zeta,zeta0t)), 0.1)
      end do
CC         c_t = c_u
      bix(1) = zeta
      bix(2) = alam*z*dtvir/u_a**2
      bix(3) = re
      bix(4) = log(z0m/z0t)
      bix(5) = z0m
      bix(6) = z0t
      bix(7) = 0.
      bix(8) = c_u
      bix(9) = c_t
      bix(10) = vkc*ust*z/phi_m(zeta)
      bix(11) = phi_m(zeta) / phi_t(zeta)

      bix2=dble(bix)

      return
      END     
                        

      real function phi_m(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
      a=1.0

c*   universal function for momentum: businger et al. (1971) formulation
c*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = alpha_m + beta_m*zeta     !  businger et al.
          y = 1. + zeta*(a+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
          if(zeta. ge. zeta_st) then
              y = (1. - gamma_m*zeta)**(p_m)
          else
              y = (a_m - b_m*zeta)**(q_m)
          end if
      end if
      phi_m = y

      return
      end

      real function phi_t(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
      a=1.0

c*   universal function for heat: businger et al. (1971) formulation
c*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = alpha_t + beta_t*zeta
          y = 1. + zeta*(
     &    a*sqrt(1.+2.*zeta*a/3.)+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
          if(zeta .ge. zeta_st) then
              y = (1. - gamma_t*zeta)**(p_t)
          else
              y = (a_t - b_t*zeta)**(q_t)
          end if
      end if
      phi_t = y

      return
      end

      real function psi_m1(x)

      real x,pi
      pi = 3.1415926
c*   integrated universal function for momentum based on businger et al.
c*   (1971) formulation for zeta_st .le. zeta .le. 0.
c*    p_m = - 1./4. (paulson, 1970)
      psi_m1 = log(0.125*(1. + x**2)*(1. + x)**2) -2.*atan(x) + pi/2.
      return
      end

      real function psi_m2(a,x)
      real x,a
c*   integrated universal function for momentum based on large et al.
c*   (1994) formulation for zeta .le. zeta_st
c*    only for q_m = -1./3.
      psi_m2 = (0.5*log(x**2 + x + 1.)
     &         + sqrt(3.)*atan((2.*x + 1.)/sqrt(3.)))
     &         / sign(1.,a)*abs(a)**(1./3.)
      return
      end

      real function psi_m(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,zeta0,x,x_st,psi_m2,psi_m1,x0,x0_st,a
      
      a=1.0
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = (1.-alpha_m)*log(zeta/zeta0) - beta_m*(zeta-zeta0)
          y = - (a*zeta+
     &      0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+(0.667*5./0.35))  ! b&h
      else
          if(zeta .ge. zeta_st) then
              x = (1. - gamma_m*zeta)**(- p_m)
              y = psi_m1(x)
          else
              x0 = 1. - (b_m/a_m)*zeta
              x = sign(1.,x0)*abs(x0)**(- q_m)
              x0_st = 1. - (b_m/a_m)*zeta_st
              x_st = sign(1.,x0_st)*abs(x0_st)**(- q_m)
              y = psi_m1(zeta_st) + psi_m2(a_m,x) - psi_m2(a_m,x_st)
          end if
      end if
      psi_m = y
      return
      end

      real function psi_t1(x)
      real x
c*   integrated universal function for heat based on businger et al.
c*   (1971) formulation for zeta_st .le. zeta .le. 0.
c*    p_t = - 1./2. (paulson, 1970)
      psi_t1 = log(0.25*(1. + x)**2) !in (paulsen, 1970) log(0.25*(1. + x**2)**2)!
      return
      end

      real function psi_t(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,
     *                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t,
     *                zeta_st,a_m,a_t,b_m,b_t
      real zeta,zeta0,y,x,x0,psi_t1,x0_st,x_st,psi_m2,a
      a=1.0
c*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
c         y = (1.-alpha_t)*log(zeta/zeta0) - beta_t*(zeta-zeta0)
          y = - (a*(1.+(2./3.)*zeta*a)**(1.5)
     &   + 0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+((0.667*5./0.35)-1.))
      else
          if(zeta .ge. zeta_st) then
              x = (1. - gamma_t*zeta)**(- p_t)
              y = psi_t1(x)
          else
              x0 = 1. - (b_t/a_t)*zeta
              x = sign(1.,x0)*abs(x0)**(- q_t)
              x0_st = 1. - (b_t/a_t)*zeta_st
              x_st = sign(1.,x0_st)*abs(x0_st)**(- q_t)
              y = psi_t1(x_st) + psi_m2(a_t,x) - psi_m2(a_t,x_st)
          end if
      end if
      psi_t = y
      return
      end
