      SUBROUTINE TSKIN_SOLVER(Tsurf,extwat,extice,fetch,dt,hskin)

!     The subroutine TSKIN_SOLVER performs iterations to
!     find the skin temperature from heat balance equation
      use DRIVING_PARAMS, only: &
      & path
      use ARRAYS, only: &
      & Tskin , &
      & deltaskin
      use PHYS_CONSTANTS2, only : &
      & sabs      
      implicit none

      real(8), intent(in) :: Tsurf
      real(8), intent(in) :: extwat
      real(8), intent(in) :: extice
      real(8), intent(in) :: fetch
      real(8), intent(in) :: dt
      
      real(8), intent(out) :: hskin

      real(8) :: Tskin1,Tskin2,Tskin3,Balskin1,Balskin2,Balskin3
      real(8) :: interval2
      
      integer iter, maxiter
      data          maxiter /500/

      real(8), external :: HEATBALANCESKIN
      real(8), external :: SKIN_THICKNESS
      real(8), external :: FRAC_SOLAR_ABSORB_SKIN

      logical firstcall
      data    firstcall /.true./

      if (firstcall) then
        open (998,file=path(1:len_trim(path))//'results/debug/iter_Tskin.dat')
      endif

      deltaskin = SKIN_THICKNESS()
	  sabs = FRAC_SOLAR_ABSORB_SKIN(deltaskin)

      interval2 = dmax1((dabs(Tsurf-Tskin(1)) + 0.1d0),10.d0)

      Tskin1 = Tsurf - interval2 !10.
      Tskin2 = Tsurf + interval2 !10.

      Balskin1 = HEATBALANCESKIN(Tskin1,Tskin(1),Tsurf, &
      & extwat,extice,fetch,deltaskin,dt,hskin)
      Balskin2 = HEATBALANCESKIN(Tskin2,Tskin(1),Tsurf, &
      & extwat,extice,fetch,deltaskin,dt,hskin)
      if (Balskin1*Balskin2>0) then
        STOP 'Chorde method is not applicable for skin temperature equation: STOP'
      endif

!     print*, Balskin1, Balskin2

      iter=0
      Balskin3=1.d0
      
!     CHORDE METHOD TO FIND SURFACE TEMPERATURE

      cyskin: do while (dabs(Balskin3)>0.1d0)
        iter=iter+1
        Balskin1 = HEATBALANCESKIN(Tskin1,Tskin(1),Tsurf, &
        & extwat,extice,fetch,deltaskin,dt,hskin)
        Balskin2 = HEATBALANCESKIN(Tskin2,Tskin(1),Tsurf, &
        & extwat,extice,fetch,deltaskin,dt,hskin)
        
        Tskin3=(Tskin1*Balskin2-Tskin2*Balskin1)/(Balskin2-Balskin1)
        
        Balskin3 = HEATBALANCESKIN(Tskin3,Tskin(1),Tsurf, &
        & extwat,extice,fetch,deltaskin,dt,hskin)
        if (iter>maxiter) then
          if (dabs(Tskin1-Tskin2)<0.01d0) then
            write (998,*) Balskin3
            exit cyskin
          else
            STOP 'Skin temperature iterations do not converge'
          endif        
        endif
        if     (Balskin1*Balskin3<0.) then
          Tskin2=Tskin3
        elseif (Balskin2*Balskin3<0.) then
          Tskin1=Tskin3
        endif
      enddo cyskin

!     print*, 'Skin disbalance', Balskin3, Tskin3, hskin
!     read*

      Tskin(2) = Tskin3

      if (firstcall) firstcall=.false.
      RETURN
      END SUBROUTINE TSKIN_SOLVER



 !     FUNCTION HEATBALANCESKIN(Tskin_iter,Tskin_prev,Tsurf, &
 !     & extwat,extice,fetch,deltaskin,dt,hskin)
      FUNCTION HEATBALANCESKIN(Tskin_iter,Tskin_prev,Tsurf, &
      & extwat,extice,fetch,deltaskin,dt,hskin)
      
!     The function HEATBALANCESKIN computes the heat balance
!     of the skin at the top of water coloumn

      use PHYS_CONSTANTS2, only: &
      & cw, &
      & row0, &
      & lamw0, &
	  & sabs
      use ATMOS, only: &
      & hflux, &
      & Elatent, &
      & Radbal
!	  use ARRAYS, only:&
!	  & deltaskin   

      implicit none
      
      real(8) :: HEATBALANCESKIN

      real(8), intent(in) :: Tskin_iter
      real(8), intent(in) :: Tskin_prev
      real(8), intent(in) :: Tsurf
      real(8), intent(in) :: extwat
      real(8), intent(in) :: extice
      real(8), intent(in) :: fetch
      real(8), intent(in) :: deltaskin      
      real(8), intent(in) :: dt
      
      real(8), intent(out):: hskin

	  real(8), external :: FRAC_SOLAR_ABSORB_SKIN
      real(8), external :: SKIN_THICKNESS
      real(8), external :: HEATCONTENTSKIN

      integer(4), save :: surftyp
      data                surftyp /4/
  

      call SENSLATMOM_FLUXES(Tskin_iter,fetch)

!     deltaskin = SKIN_THICKNESS()
!	  sabs = FRAC_SOLAR_ABSORB_SKIN(deltaskin)
!      sabs=0.03
      call RADBALANCE       (Tskin_iter,extwat,extice,surftyp)
!      deltaskin = SKIN_THICKNESS()
      hskin = - lamw0*(Tsurf-Tskin_iter)/deltaskin

      HEATBALANCESKIN = HEATCONTENTSKIN(Tskin_iter,Tskin_prev,deltaskin,dt) - &
      & Radbal + hflux + Elatent + hskin
      
      END FUNCTION HEATBALANCESKIN
      
      
      FUNCTION HEATCONTENTSKIN(Tskin_iter,Tskin_prev,deltaskin,dt)
      
!     Function HEATCONTENTSKIN calculates the rate of change of
!     skin heat content

      use PHYS_CONSTANTS2, only: &
      & cw_m_row0
      
      implicit none
      real(8) :: HEATCONTENTSKIN
      
      real(8), intent(in) :: Tskin_iter
      real(8), intent(in) :: Tskin_prev
      real(8), intent(in) :: deltaskin
      real(8), intent(in) :: dt
      
      HEATCONTENTSKIN = deltaskin*cw_m_row0*(Tskin_iter-Tskin_prev)/dt
      
      END FUNCTION HEATCONTENTSKIN
      
      
      FUNCTION SKIN_THICKNESS()
      
!     Function SKIN_THICKNESS calculates the thickness of cool skin
!     at the water surface
      
	  use TURB_CONST, only : &
      & niu
	  use PHYS_CONSTANTS2, only : &
	  & roa0, &
	  & row0, &
	  & g,    &
	  & cw,   &
	  & lamw0
	  use ATMOS, only : &
	  & velfrict,       &
	  & velfrict_prev , &
	  & hflux_prev,     &
      & elatent_prev,   &
      & Radbal_prev
	  

      implicit none

	  real(8), parameter :: const_Saunders = 6.d0 ! Dimensionless
      real(8), parameter :: SKIN_THICKNESS_min = 0.0001d0
      real(8), parameter :: SKIN_THICKNESS_max = 0.005d0
	  real(8), parameter :: sw_alfa = 0.00021d0     ! thermal expansion coeff, C -1

      real(8) :: SKIN_THICKNESS
	  real(8) :: const_Fairall
      real(8) :: big_c
	  real(8) :: potok
	  big_c = 16.*cw*row0*g*sw_alfa*niu**3/(((roa0/row0)**2)*(lamw0**2))
	  potok=Radbal_prev-hflux_prev-elatent_prev
	  
	  const_Fairall = const_Saunders/(1+(big_c*(abs(potok))/velfrict_prev**4)**(3./4.))**(1./3.)

      SKIN_THICKNESS = 0.001d0
	  if (velfrict.ne.0) then
        SKIN_THICKNESS = const_Saunders*niu/(sqrt(roa0/row0)*velfrict_prev)
	  else
        SKIN_THICKNESS = 0.005d0
	  endif
      
      continue

      SKIN_THICKNESS = dmax1(SKIN_THICKNESS, SKIN_THICKNESS_min)
	  SKIN_THICKNESS = dmin1(SKIN_THICKNESS, SKIN_THICKNESS_max)
      
      END FUNCTION SKIN_THICKNESS



      FUNCTION FRAC_SOLAR_ABSORB_SKIN(skin_thickness)
      
!     Function FRAC_SOLAR_ABSORB_SKIN calculates the fraction of solar radiation
!     absorbed by cool skin at the water surface according to Fairall et al., 1996
      
      implicit none
      
      real(8) :: FRAC_SOLAR_ABSORB_SKIN
      
      real(8), parameter :: c1 = 6.5d-2
      real(8), parameter :: c2 = 1.1d+1
      real(8), parameter :: c3 = 6.6d-5
      real(8), parameter :: c4 = 8.0d-4
      
!     Input variables      
      real(8), intent(in) :: skin_thickness
      
      FRAC_SOLAR_ABSORB_SKIN = c1 + c2*skin_thickness - c3/skin_thickness * &
      & (1 - dexp(-skin_thickness/c4))
      
      END FUNCTION FRAC_SOLAR_ABSORB_SKIN