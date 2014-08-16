SUBROUTINE K_EPSILON(ix, iy, nx, ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, dt, &
& b0, tau_air, tau_i, tau_gr)
      
! KT_eq calculates eddy diffusivity (ED) in water coloumn following parameterizations:

! 1) empirical profile of ED;
! 2) Semi-empirical formulations for ED;
! 2) k-epsilon parameterization for ED, including kolmogorov relation.
      
use ATMOS
use DRIVING_PARAMS
use ARRAYS
use PHYS_CONSTANTS2, only: &
& row0, roa0, &
& lamw0, &
& cw, &
& g

use TURB_CONST
use WATER_DENSITY, only: &
& water_dens_t_unesco, &
& water_dens_ts

implicit none

integer(4), parameter :: vector_length = 350

! Input variables

! Reals
real(8), intent(in) :: dt
real(8), intent(in) :: kor
real(8), intent(in) :: h_veg
real(8), intent(in) :: a_veg
real(8), intent(in) :: c_veg
real(8), intent(in) :: hour
real(8), intent(in) :: b0
real(8), intent(in) :: tau_air
real(8), intent(in) :: tau_i
real(8), intent(in) :: tau_gr

! Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day
        
! Local variables     
! Reals
real(8) :: C_eps1(M)
real(8) :: C_eps2(M)
real(8) :: C_eps3(M)
real(8) :: CE(M)
real(8) :: CEt(M)
real(8) :: lam_E(M+1)
real(8) :: lam_eps(M+1)

real(8) :: a(vector_length)
real(8) :: b(vector_length)
real(8) :: c(vector_length)
real(8) :: d(vector_length)

!real(8) :: am(vector_length,2,2)
!real(8) :: bm(vector_length,2,2)
!real(8) :: cm(vector_length,2,2)
!real(8) :: ym(vector_length,2)
!real(8) :: dm(vector_length,2)

real(8) :: work(1:M+1,1:2) ! Work array

real(8) :: AG(5)

real(8) :: wind10
real(8) :: urel
real(8) :: vrel
real(8) :: u
real(8) :: v
real(8) :: wr

real(8) :: pi
real(8) :: ufr
real(8) :: dist_surf
real(8) :: ACC
real(8) :: ACC2
real(8) :: ACCk
real(8) :: dE_it
real(8) :: deps_it
real(8) :: E_mid
real(8) :: eps_mid
real(8) :: xx
real(8) :: yy
real(8) :: lm
real(8) :: ext_lamw
real(8) :: month_sec
real(8) :: day_sec
real(8) :: hour_sec
real(8) :: al_it
real(8) :: FS
real(8) :: FB
real(8) :: al
real(8) :: alft
real(8) :: GAMUN
real(8) :: urels
real(8) :: vrels
real(8) :: minwind
real(8) :: ceps3
real(8) :: lambda_M
real(8) :: lambda_N
real(8) :: E_min
real(8) :: eps_min
real(8) :: vdamp
real(8) :: dE_it_max
real(8) :: deps_it_max
     
real(8) :: GRADU
real(8) :: GRADV
real(8) :: GRADT
real(8) :: GAMVN
real(8) :: TFR
real(8) :: COGRTN
real(8) :: COGRUN
real(8) :: COGRVN
real(8) :: DTDZH

! Integers
integer(4), parameter :: maxiter = 35 !35
integer(4) :: iter
integer(4) :: iter_t
integer(4) :: nl
integer(4) :: i
integer(4) :: j
integer(4) :: k
integer(4) :: numkeps
integer(4) :: keps_coef
! integer(4) :: stabfunc
integer(4) :: time_deriv

! Logicals
logical, allocatable :: init_keps(:,:)
logical :: indstab
logical :: ind_bound
logical :: iterat
logical :: firstcall
logical :: cycle_keps
logical :: next_tstep_keps
logical :: smooth
logical :: perdam
logical :: semi_implicit_scheme

! Characters
character :: tp_name*10
character :: numt*1

common /wind/  wind10,urel,vrel,u,v

data firstcall /.true./

! Externals
real(8), external :: KRON
real(8), external :: DZETA
real(8), external :: CE_CANUTO 
real(8), external :: CEt_CANUTO
real(8), external :: SMOMENT_GALPERIN 
real(8), external :: SHEAT_GALPERIN
real(8), external :: INTERPOLATE_TO_INTEGER_POINT

SAVE

if (firstcall) then
  if (turb_out==1) then
    write (numt,'(i1)') Turbpar
!    open(2114,file=path(1:len_trim(path))//'results/'// &
!    & 'err_progon.dat',  status = 'unknown')
!    open(2115,file=path(1:len_trim(path))// &
!    & 'results/'//'E-eps.dat', &
!    & status = 'unknown')
!    open(2116,file=path(1:len_trim(path))// &
!    & 'results/'//'E-eps1.dat', &
!    & status = 'unknown')
!    open (2117,file=path(1:len_trim(path))//'results/' &
!    & //'turb'//numt//'.dat',      status = 'unknown')
!    open (2118,file=path(1:len_trim(path))//'results/' &
!    & //'temp'//numt//'.dat',      status = 'unknown')
!    open (2119,file=path(1:len_trim(path))//'results/' &
!    & //'uv'//numt//'.dat',        status = 'unknown')
!    open (2120,file=path(1:len_trim(path))//'results/' &
!    & //'KT'//numt//'.dat',        status = 'unknown') 
!    open (2121,file=path(1:len_trim(path))//'results/' &
!    & //'dE1'//numt//'.dat',       status = 'unknown') 
!    open (2122,file=path(1:len_trim(path))//'results/' &
!    & //'dE2'//numt//'.dat',       status = 'unknown') 
    open (2123,file=path(1:len_trim(path))//'results/&
    &debug/err_keps_iter.dat',      status = 'unknown')  
  endif
  
  allocate (init_keps(1:nx, 1:ny) )
  init_keps(:,:) = .false.
       
  month_sec   = 30.*24.*60.*60.
  day_sec     = 24*60.*60.
  hour_sec    = 60.*60.

  pi          = 4.*dtan(1.d0)

  AL    = g/row0
  
! Parameters of numerical scheme
  semi_implicit_scheme = .false. ! No iterations in k-epsilon, semi-implicit scheme
  al_it   = 0.2 !0.8
  ACC     = 1.d-20 !1.d-20
  ACC2    = 1.d-20 !1.d-20
  ACCk    = 1.d-20 !1.d-20
  knum    = 0.
  minwind = 1.0d0

  eps_min = 1.d-30 ! This value of eps_min needs to be checked on consistency
                   ! with values observed in nature
  E_min   = dsqrt(eps_min*lamw0*(cw*row0*CEt0)**(-1) ) ! This expression ensures that
                                                       ! eddy diffusivity in the decaying
                                                       ! turbulence (E -> Emin, eps -> eps_min)
                                                       ! is of the order of molecular diffusivity

  dE_it_max = 1.d-9
  deps_it_max = 1.d-12

  smooth = .false.
  perdam = .false.
  vdamp  = 1.d-3 !1.d-3
 
  AG(1)   = 0.d0
  AG(2)   = 0.d0
  AG(3)   = 0.d0
! AG(4)=1. !0.
! AG(5)=1. !0.

  numkeps   = 1
  keps_coef = 1

! keps_coef = 1 - standard empirical coefficients in epsilon-equation
! keps_coef = 2 - the coefficient by (Aupoix et al., 1989) is used

! stabfunc  = 3
  
! stabfunc = 1  - stability functions CE and CEt are constants
! stabfunc = 2  - stability functions CE and CEt are dependent 
! on shear and buoyancy (Canuto et al., 2001)
! stabfunc = 3  - stability functions CE and CEt are dependent
! only on buoyancy (Galperin et al., 1988)

       
! FRICTION: VEGETATION

  if (h_veg>0) then
    print*, 'Vegetation friction parameterization &
    & is not operational currently: STOP'
    STOP
  endif

  veg_sw = 0.
  do i = 1, M+1
    if (h1-DZETA(dble(i))*h1 < h_veg) veg_sw(i) = 1.
  enddo

endif

!u=u1(1)
!v=v1(1)
!
!if (relwind==2) then
!  urel=uwind-u
!  vrel=vwind-v
!elseif (relwind==1) then
!  urel=uwind
!  vrel=vwind
!else
!  print*, 'relwind=',relwind,' - incorrect meaning: STOP' 
!  stop
!endif
!
!if (dabs(urel)<minwind) then
!  urels=sign(minwind,urel)
!else
!  urels=urel
!endif
!
!if (dabs(vrel)<minwind) then
!  vrels=sign(minwind,vrel)
!else   
!  vrels=vrel
!endif
!
!wr   = dsqrt(urels**2+vrels**2)
!tau_air = roa0*cdmw2*wr
!ufr  = dsqrt(tau_air/row0)    

!tau_wav = tau_air !*kw
!tau_sbl = tau_air !*(1.-kw)

! Calculation of Shezy coefficient
! rp = 0.1 !height of roughness elements
!if (l1/=0) then
!  Cz_gr=7.7*dsqrt(g)*(h1/2/rp)**(1./6.) !50. !ground Shezy coefficient
!else
!  Cz_gr=7.7*dsqrt(g)*(h1/rp)**(1./6.)
!endif
!rp = 0.01 
!Cz_i = 7.7*dsqrt(g)*(h1/2/rp)**(1./6.)
!
!! Friction: ground
!tau_gr = row0*g*Cz_gr**(-2)*(u1(M)**2+v1(M)**2)
!
!! Friction: ice
!tau_i = row0*g*Cz_i**(-2)*(u1(1)**2+v1(1)**2)
     
      
ifkeps: if (Turbpar==2) then

! Implementation of the Krank-Nikolson numerical scheme
! for k-epsilon parameterization using simple iterations to
! solve the set of nonlinear finite-difference equations

  Um = 0.5d0 * (u2 + u1)
  Vm = 0.5d0 * (v2 + v1)
  
  next_tstep_keps = .false.

  keps_mode: do while (.not.next_tstep_keps)
  
  if (.not.init_keps(ix,iy)) then
! Initialization mode of k-epsilon parameterization: time derivatives are neglected  
    time_deriv = 0
  else
! Evolution mode of k-epsilon parameterization: full equations are solved
    time_deriv = 1
  endif
  
  E_it1   = E1
  eps_it1 = eps1
  iter_t = 0
  k2_mid = 0.d0
  k2     = 0.d0

  DTDZH=(Tw1(2)-Tw1(1))/(h1*ddz(1))
  ALFT=1.d0 

  cycle_keps = .true.
  
  if (semi_implicit_scheme) then
    iterat = .false.
  else
    iterat = .true.
  endif  

! do i = 1, M
!   if (stabfunc == 1) then
!     lam_eps(i) = lam_eps0
!     lam_E  (i) = lam_E0
!     CE(i) = CE0
!     CEt    (i) = CEt0
!   else
!     lambda_N = E1(i)*E1(i)/((eps1(i)+ACC2)*(eps1(i)+ACC2))*
!     & ((row(i+1)-row(i))/(h1*ddz(i)))*AL
!     lambda_M = E1(i)*E1(i)/((eps1(i)+ACC2)*(eps1(i)+ACC2))*
!     & ( (Um(i+1)-Um(i))*(Um(i+1)-Um(i)) +
!     &   (Vm(i+1)-Vm(i))*(Vm(i+1)-Vm(i)) )/
!     &   (h1*h1*ddz(i)*ddz(i))
!     CE(i)  = CE_CANUTO (lambda_M, lambda_N)
!     CEt(i) = CEt_CANUTO(lambda_M, lambda_N)
!     lam_eps(i) = CE(i)/sigmaeps
!     lam_E  (i) = CE(i)/sigmaE
!   endif
! enddo
          
! 10

  iterat_keps: do while (cycle_keps)

! The cycle implements iterations to
! solve the set of nonlinear finite-difference equations
! of k-epsilon parameterization

   E12   = 0.5d0*(E1   + E_it1)
   eps12 = 0.5d0*(eps1 + eps_it1)
 
!    E12   =  E_it1
!    eps12 =  eps_it1

    do i = 1, M
      if (stabfunc == 1) then
        lam_eps(i) = lam_eps0
        lam_E  (i) = lam_E0
        CE     (i) = CE0
        CEt    (i) = CEt0
      else
        lambda_N = E12(i)*E12(i)/((eps12(i)+ACC2)*(eps12(i)+ACC2))* & 
        & ((row(i+1)-row(i))/(h1*ddz(i)))*AL
        lambda_M = E12(i)*E12(i)/((eps12(i)+ACC2)*(eps12(i)+ACC2))* &
        & ( (Um(i+1)-Um(i))*(Um(i+1)-Um(i)) + &
        &   (Vm(i+1)-Vm(i))*(Vm(i+1)-Vm(i)) ) / &
        &   (h1*h1*ddz(i)*ddz(i))
        if (stabfunc == 2) then
          CE(i)  = CE_CANUTO (lambda_M, lambda_N)
          CEt(i) = CEt_CANUTO(lambda_M, lambda_N)
        elseif (stabfunc == 3) then
          CE(i)  = dsqrt(2.d0)*CL_K_KL_MODEL*SMOMENT_GALPERIN(lambda_N)
          CEt(i) = dsqrt(2.d0)*CL_K_KL_MODEL*SHEAT_GALPERIN(lambda_N)
        endif
        lam_eps(i) = CE(i)/sigmaeps
        lam_E  (i) = CE(i)/sigmaE
      endif
    enddo
              
    do i = 2, M
!      E_mid=0.5d0*(E12(i-1)+E12(i))
!      eps_mid=0.5d0*(eps12(i-1)+eps12(i))
      E_mid = INTERPOLATE_TO_INTEGER_POINT(E12(i-1),E12(i),ddz(i-1),ddz(i))
      eps_mid = INTERPOLATE_TO_INTEGER_POINT(eps12(i-1),eps12(i),ddz(i-1),ddz(i))
      if (iter_t==0) then
        k2_mid(i) = dmax1(E_mid*E_mid/(eps_mid+ACCk),1.d-8)
      else
        k2_mid(i)=k2_mid(i)*al_it + &
        & dmax1(E_mid*E_mid/(eps_mid+ACCk),1.d-8)*(1-al_it)
      endif
      k5_mid(i) = niu + k2_mid(i)* &
      & INTERPOLATE_TO_INTEGER_POINT(lam_eps(i-1),lam_eps(i),ddz(i-1),ddz(i))
      k3_mid(i) = niu + k2_mid(i)* & 
      & INTERPOLATE_TO_INTEGER_POINT(lam_E(i-1),lam_E(i),ddz(i-1),ddz(i))
    enddo
    
!    print*, iter_t,k2_mid
!    read*

    do i = 1, M
      if (iter_t==0) then
        k2(i) = dmax1(E12(i)*E12(i)/(eps12(i)+ACCk),1.d-8)
      else
        k2(i)=k2(i)*al_it + &
        & dmax1(E12(i)*E12(i)/(eps12(i)+ACCk),1.d-8)*(1-al_it)
      endif
      k5(i) = niu + k2(i)*lam_eps(i)
      l(i)=E12(i)**(3./2.)/(eps12(i)+ACC)
      TF(i)=dSQRT(E12(i))/(l(i)+ACC)
    enddo

    do i=2,M-1
      WU_(i)=(E_it1(i+1)-E_it1(i-1))*(U2(i)-U2(i-1))/ &
      & (h1**2*(0.5*ddz(i-1)+0.5*ddz(i+1)+ddz(i))*ddz(i-1) )
      WV_(i)=(E_it1(i+1)-E_it1(i-1))*(V2(i)-V2(i-1))/ &
      & (h1**2*(0.5*ddz(i-1)+0.5*ddz(i+1)+ddz(i))*ddz(i-1) )
    enddo
       
    do i=2,M
      GRADT = (Tw1(i+1)-Tw1(i-1))/(h1*(ddz(i-1)+ddz(i) ) )
      GRADU = (U2(i)-U2(i-1))/(h1*ddz(i-1) )
      GRADV = (V2(i)-V2(i-1))/(h1*ddz(i-1) )
      GAMUN = - CONUV*(WU_(i+1)-WU_(i-1))/(h1*(ddz(i-1)+ddz(i) ))
      GAMVN = - CONUV*(WV_(i+1)-WV_(i-1))/(h1*(ddz(i-1)+ddz(i) ))
      COGRTN = DTDZH * 0.1
      TFR=E_it1(i)/(l(i)*l(i)+ACC) + acc
      COGRUN = GAMUN/TFR
      COGRVN = GAMVN/TFR 
      if(gradT < 0.) then
        GAMT(i) = - AG(3)*cogrtn
        GAMU(i) = - AG(1)*cogrun
        GAMV(i) = - AG(2)*cogrvn
      else
        GAMT(i) = - AG(3)*(AL*GRADT**2+CON0*TFR*COGRTN) / &
        &            (AL*GRADT+CON0*TFR)
        GAMU(i) = - AG(1)*(AL*(CON2*GRADT+(CON2-1.)*GAMT(i))*GRADU + &
        &            CON1*TFR*COGRUN)/(AL*GRADT+CON1*TFR)
        GAMV(i) = - AG(2)*(AL*(CON2*GRADT+(CON2-1.)*GAMT(i))*GRADV + &
        &            CON1*TFR*COGRVN)/(AL*GRADT+CON1*TFR)
      end if
    enddo
    GAMT(1)=0.
    GAMU(1)=0.
    GAMV(1)=0.

    do i=1,M
      Gen(i) = ((Um(i+1)-Um(i))*(Um(i+1)-Um(i)+h1*ddz(i)*GAMU(i)) + &
      &        (Vm(i+1)-Vm(i))*(Vm(i+1)-Vm(i) +h1*ddz(i)*GAMV(i))) / &
      &        (h1**2*ddz(i)**2)*CE(i)
      S(i) = -((row(i+1)-row(i))/(h1*ddz(i))+GAMT(i))*ALFT*AL*CEt(i)
      F(i) = Gen(i)+S(i)
    enddo
      
! Solution of the equation for TKE (turbulent kinetic energy)

! Boundary condition at the top
    if (l1==0) then 
      FS=kwe*(tau_air/row0)**(3./2.)
    elseif (l1>0.and.l1<=L0) then
      FS=((b0*tau_i+(1-b0)*tau_air)/row0)**(3./2.)
    elseif (l1>L0) then
      FS=(tau_i/row0)**(3./2.)
    endif
    xx=0.5*(k3_mid(2)*dt/(h1**2*ddz(1)*0.5*(ddz(1)+ddz(2)) ) + &
    &  time_deriv*DZETA(1.5d0)*dhw/(h1*0.5*(ddz(1)+ddz(2)) ) - &
    &  time_deriv*dhw0/(h1*0.5*(ddz(1)+ddz(2)) ) )
    yy=dt/(h1*ddz(1))
    b(1)=-xx
    c(1)=-(xx+time_deriv*1.d0)-(dABS(S(1))-S(1))/(TF(1)+ACC)*dt/2-TF(1)*dt     
    d(1)=-time_deriv*E1(1)-xx*(E1(2)-E1(1))-yy*FS - &
    & k2(1)*(S(1)+dabs(S(1)))*dt/2-k2(1)*Gen(1)*dt !+eps_it1(i)*dt

! Boundary condition at the bottom       
    FB=-(tau_gr/row0)**(3./2.)
    xx=0.5*(k3_mid(M)*dt/(h1**2*ddz(M)*0.5*(ddz(M)+ddz(M-1)) ) - &
    & time_deriv*DZETA(dble(M+0.5))*dhw/(h1*0.5*(ddz(M)+ddz(M-1)) ) + &
    & time_deriv*dhw0/(h1*0.5*(ddz(M)+ddz(M-1)) ) )
    yy=dt/(h1*ddz(M))
    a(M)=-xx
    c(M)=-(xx+time_deriv*1.d0)-(dABS(S(M))-S(M))/(TF(M)+ACC)*dt/2 - &
    & TF(M)*dt
    d(M)=-time_deriv*E1(M)-xx*(E1(M-1)-E1(M))+yy*FB - &
    & k2(M)*(S(M)+dabs(S(M)))*dt/2-k2(M)*Gen(M)*dt !+eps_it1(i)*dt

! The coefficients of linear algebraic equations in the internal points of the mesh
    do i=2,M-1
      a(i)=time_deriv*dzeta(dble(i+0.5))*dhw/(2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & time_deriv*dhw0/(2.*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & k3_mid(i)*dt/(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )

      b(i)=-time_deriv*dzeta(dble(i+0.5))*dhw/(2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) + &
      & time_deriv*dhw0/(2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & k3_mid(i+1)*dt/(h1**2*ddz(i)*(ddz(i+1)+ddz(i)) )

      c(i)=a(i)+b(i)-time_deriv*1.d0-(dABS(S(i))-S(i))/(TF(i)+ACC)*dt/2-TF(i)*dt

      d(i)=-time_deriv*E1(i)-k2(i)*(S(i)+dabs(S(i)))*dt/2-k2(i)*Gen(i)*dt !+eps_it1(i)*dt
      d(i)=d(i)-dt*(h1**2*ddz(i)*(ddz(i)+ddz(i+1)) )**(-1)* &
      & (k3_mid(i+1)*(E1(i+1)-E1(i)))   
      d(i)=d(i)+dt*(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )**(-1)* &
      & (k3_mid(i)*(E1(i)-E1(i-1)))
      d(i)=d(i)-time_deriv*dzeta(dble(i+0.5))*dhw* &
      & (2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1)**(-1)* &
      & (E1(i+1)-E1(i-1)) 
      d(i)=d(i)+time_deriv*dhw0* &
      & (2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1)**(-1)* &
      & (E1(i+1)-E1(i-1)) 
    enddo
    ind_bound=.false.
    call ind_stab_fact_db (a,b,c,1,M,indstab,ind_bound)
    if (indstab==.false.) then
      do i=2,M-1
        a(i)=-k3_mid(i)*dt/(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )
        b(i)=-k3_mid(i+1)*dt/(h1**2*ddz(i)*(ddz(i+1)+ddz(i)) )
      enddo
    endif 
    call PROGONKA (a,b,c,d,E_it2,1,M)
    E_it2(M+1)=0.
                
    call CHECK_MIN_VALUE(E_it2, M, E_min) 

    dE_it=maxval(dabs(E_it1-E_it2))

! C1 is given following Satyanarayana et al., 1999; Aupoix et al., 1989

    if (keps_coef==1) then
      do i = 1, M
        if (S(i)<0.) then
          ceps3=-0.4d0
!         ceps3 < -0.4 should be used for stable stratification (Burchard, 2002)
        else
          ceps3= 1.14d0 ! Baum and Capony (1992)
!         Other values for ceps3:
!         0.<ceps3<0.29 for stable and ceps3 = 1.44 for unstable conditions (Rodi, 1987);
!         ceps3 >= 1 (Kochergin and Sklyar, 1992)
        endif
        C_eps1(i) = ceps1
        C_eps2(i) = ceps2
        C_eps3(i) = ceps3
      enddo
    elseif (keps_coef==2) then
      do i = 1, M 
        Re(i)=ACC+(2*E_it1(i)/3.)**2/(eps_it1(i)*niu+ACC)
        C1aup(i)=Co/(1.d0+0.69d0*(2.d0-Co)/dsqrt(Re(i)))
        C_eps1(i) = C1aup(i)
        C_eps2(i) = C1aup(i)
        C_eps3(i) = C1aup(i)
      enddo  
    endif

! The solution of the equation for dissipation rate

! dist_surf is the distance from the surface, 
! where the b.c. for epsilon are formulated

    dist_surf = 0.d0 !0.5d0*ddz(1)*h1

! Boundary condition at the top
    if (l1==0) then
      FS = CE0**(0.75d0)*k5(1)*E_it1(1)**(3./2.)/kar* &
      & (roughness + dist_surf)**(-2)
    else
      FS = CE0**(0.75d0)*k5(1)*E_it1(1)**(3./2.)/kar* &
      & (0.01d0 + dist_surf)**(-2) !0.01 m - roughness of ice 
    endif 
    xx=0.5*(k5_mid(2)*dt/(h1**2*ddz(1)*0.5*(ddz(1)+ddz(2)) ) + &
    & time_deriv*DZETA(1.5d0)*dhw/(h1*0.5*(ddz(1)+ddz(2)) ) - &
    & time_deriv*dhw0/(h1*0.5*(ddz(1)+ddz(2)) ) )
    yy=dt/(h1*ddz(1))
    b(1)=-xx
    c(1)=-(xx+time_deriv*1.d0)-C_eps2(1)*TF(1)*dt + &
    & C_eps3(1)*0.5*(-dABS(S(1))+S(1))/(TF(1)+ACC)*dt
    d(1)=-time_deriv*eps1(1)-xx*(eps1(2)-eps1(1))-yy*FS-0.5*C_eps3(1) * &
    & (dABS(S(1))+S(1))*TF(1)*k2(1)*dt - &
    & C_eps1(1)*TF(1)*k2(1)*Gen(1)*dt !+eps_it1(i)*dt
 
 ! Boundary condition at the bottom
    dist_surf = 0.d0 ! ddz(M)/2*h1
    FB = - CE0**(0.75d0)*k5(M)* &
    & E_it1(M)**(3./2.)/kar*(0.01d0 + dist_surf)**(-2) 
    xx=0.5*(k5_mid(M)*dt/(h1**2*ddz(M)*0.5*(ddz(M)+ddz(M-1)) ) - &
    & time_deriv*DZETA(dble(M+0.5))*dhw/(h1*0.5*(ddz(M)+ddz(M-1)) ) + &
    & time_deriv*dhw0/(h1*0.5*(ddz(M)+ddz(M-1)) ))
    yy=dt/(h1*ddz(M))
    a(M)=-xx
    c(M)=-(xx+time_deriv*1.d0)-C_eps2(M)*TF(M)*dt + &
    & C_eps3(M)*0.5*(-dABS(S(M))+S(M))/(TF(M)+ACC)*dt
    d(M)=-time_deriv*eps1(M)-xx*(eps1(M-1)-eps1(M))+yy*FB-0.5*C_eps3(M)* &
    & (dABS(S(M))+S(M))*TF(M)*k2(M)*dt - &
    & C_eps1(M)*TF(M)*k2(M)*Gen(M)*dt

    do i=2,M-1
      a(i)=time_deriv*dzeta(dble(i+0.5))*dhw/(2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & time_deriv*dhw0/(2.*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & k5_mid(i)*dt/(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )

      b(i)=-time_deriv*dzeta(dble(i+0.5))*dhw/(2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) + &
      & time_deriv*dhw0/(2.*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1) - &
      & k5_mid(i+1)*dt/(h1**2*ddz(i)*(ddz(i+1)+ddz(i)) )

      c(i)=a(i)+b(i)-time_deriv*1.d0-C_eps2(i)*TF(i)*dt + &
      & C_eps3(i)*0.5*(-dABS(S(i))+S(i))/(TF(i)+ACC)*dt

      d(i)=-time_deriv*eps1(i)-0.5*C_eps3(i)*(dABS(S(i))+S(i))*TF(i)*k2(i)*dt- &
      & C_eps1(i)*TF(i)*k2(i)*Gen(i)*dt
      d(i)=d(i)-dt*(h1**2*ddz(i)*(ddz(i)+ddz(i+1)) )**(-1)* &
      & (k5_mid(i+1)*(eps1(i+1)-eps1(i)))
      d(i)=d(i)+dt*(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )**(-1)* &
      & (k5_mid(i)*(eps1(i)-eps1(i-1)))
      d(i)=d(i)-time_deriv*dzeta(dble(i+0.5))*dhw* &
      & (2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1)**(-1)* &
      & (eps1(i+1)-eps1(i-1))
      d(i)=d(i)+time_deriv*dhw0* &
      & (2*(ddz(i)+0.5*ddz(i-1)+0.5*ddz(i+1))*h1)**(-1)* &
      & (eps1(i+1)-eps1(i-1))
    enddo
    ind_bound=.false.
    call ind_stab_fact_db (a,b,c,1,M,indstab,ind_bound)
    if (indstab==.false.) then
      do i=2,M-1
        a(i)=-k5_mid(i)*dt/(h1**2*ddz(i)*(ddz(i)+ddz(i-1)) )
        b(i)=-k5_mid(i+1)*dt/(h1**2*ddz(i)*(ddz(i+1)+ddz(i)) )
      enddo
    endif  

    call PROGONKA (a,b,c,d,eps_it2,1,M)
    eps_it2(M+1)=0.
 
    call CHECK_MIN_VALUE(eps_it2, M, eps_min)

    deps_it=maxval(dabs(eps_it1-eps_it2))
    
    continue   

    if (iterat) then
      if ((dE_it>dE_it_max .or. deps_it>deps_it_max).and.iter_t < maxiter) then
        E_it1   = E_it1   * al_it + E_it2   * (1-al_it)
        eps_it1 = eps_it1 * al_it + eps_it2 * (1-al_it)
        iter=iter+1
        iter_t=iter_t + 1
        if (iter==8) iter=0
      else
        if (dE_it<1.E-6.and.deps_it<1.E-7) then
          E2=E_it2
          eps2=eps_it2
          iter=0
          iter_t=0
          nl=2
          cycle_keps = .false.
        else
          eps_it1=eps1
          E_it1=E1
          iterat=.false.
          iter_t = 0
          nl=1
        endif
      endif
    else
      E2=E_it2
      eps2=eps_it2
!      E2   = E_it1   * al_it + E_it2   * (1-al_it)
!      eps2 = eps_it1 * al_it + eps_it2 * (1-al_it)
      iter=0
      iter_t=0
      cycle_keps = .false.
    endif
!    if (dE_it<1.E-10.or.deps_it<1.E-13) nl=3
              
  enddo iterat_keps
  
  if (.not.init_keps(ix,iy)) then
! K-epsilon solver is implemented in initialization mode, i.e.
! to obtain initial profiles of E1 and eps1. 
! Time derivatives in k-epslilon equations are neglected in this case.
    init_keps(ix,iy) = .true.
    E1 = E2
    eps1 = eps2
  else
! The case when k-epsilon solver is implemented in evolution mode, i.e.
! E2 and epsilon2 are obtained at the next time step.
    next_tstep_keps = .true.
  endif
  
  enddo keps_mode

! Smoothing of the k and epsilon profiles
  if (smooth) then
    call VSMOP_LAKE(E2,E2,lbound(E2,1),ubound(E2,1),1,M,vdamp,perdam)
    call CHECK_MIN_VALUE(E2, M, E_min)
    call VSMOP_LAKE(eps2,eps2,lbound(eps2,1),ubound(eps2,1),1,M,vdamp,perdam)
    call CHECK_MIN_VALUE(eps2, M, eps_min)
  endif

  E12=0.5d0*(E1+E2)
  eps12=0.5d0*(eps1+eps2)

  do i=1,M
    KT(i)=lam_T*dmax1(CEt(i)*E2(i)**2/(eps2(i)+ACCk),1.d-8)*row0*cw
  enddo

! Diagnostic calculation
  do i=2, M-1
    TKE_turb_trans(i) = 1.d0/(h1*h1*ddz(i))* &
    & (k3_mid(i+1)*(E2(i+1) - E2(i)  ) / &
    & (0.5d0*(ddz(i)   + ddz(i+1)))  - &
    & k3_mid(i)  *(E2(i)   - E2(i-1)) / &
    & (0.5d0*(ddz(i-1) + ddz(i)  )) )
  enddo
      
! In the top and bottom layers TKE_turb_trans is currently set to zero,
! while actually it is not zero.

  TKE_turb_trans(1) = 0.d0
  TKE_turb_trans(M) = 0.d0

! Diagnostic calculations
  S  (:) = S  (:)*k2(:)
  Gen(:) = Gen(:)*k2(:)

!  do i = 2, M-1
!    work(i,1) = dhw/(dt*h1)*DZETA(dble(i)+0.5d0)*(E2(i+1)-E2(i-1))/ &
!    & (0.5d0*ddz(i-1)+ddz(i)+0.5d0*ddz(i+1))/Buoyancy0
!    work(i,2) = dhw0/(dt*h1)*(E2(i+1)-E2(i-1))/ &
!    & (0.5d0*ddz(i-1)+ddz(i)+0.5d0*ddz(i+1))/Buoyancy0
!  enddo
  
! Output
  if (turb_out==1) then
!  do i=1,M+1
!    write (2115,7) time/hour_sec+12.-3*24, dzeta(dble(i+0.5))*h1, &
!    & E2(i)*1.E9,eps2(i)*1.E9, KT(i),tau_air*1.E3,Tw1(i)
!  enddo
!  write (2116,8) time/hour_sec+12.-3*24, E2(1)*1.E9, &
!  & eps2(1)*1.E9, KT(1)
    if (dE_it > 1.d-7 .or. deps_it > 1.d-9) then
      write (2123, *) nstep, nl, dE_it, deps_it
    endif
  endif
7 format (f10.4,f8.2,2f12.1,3f8.2)
8 format (f10.4,2f12.1,f8.2)      
                  
       
  E1=E2
  eps1=eps2

endif ifkeps


SELECT CASE (Turbpar)
  CASE (2)
  CASE (1)
    do i=1,M+1
      KT(i)=(k2(i)-niu)*cw*row0 !*0.01
    enddo
  CASE DEFAULT
    KT=lamTM*k2*cw*row0 !*0.01
    if (turb_out==1) then
      if (firstcall) write (2117,'(a)') tp_name
      do i=1,M+1
        write (2117,'(F8.2, F5.2, E10.3, F8.1, F5.2, F6.2)') &
        & time/hour_sec+12.-3*24, dzeta(dble(i+0.5))*h1, Ri(i), &
        & KT(i),dzeta(dble(i))*h1,uv(i)
      enddo
    endif
ENDSELECT 
     
if (firstcall) firstcall=.false.

END SUBROUTINE K_EPSILON


!SUBROUTINE KEPS_INIT


!END SUBROUTINE KEPS_INIT


SUBROUTINE CHECK_MIN_VALUE(E,N,threshold)
implicit none

! Input variables
! Integers
integer(4), intent(in) :: N
! Reals
real(8), intent(in) :: threshold

! Input/output variables
real(8), intent(inout) :: E(1:N)

integer(4) :: i ! Loop index

do i = 1, N 
  E(i) = dmax1(E(i),threshold)
enddo
END SUBROUTINE CHECK_MIN_VALUE


REAL(8) FUNCTION CE_CANUTO(lambda_M, lambda_N)

! The function CE_CANUTO calculates stability function for momentum
! according to Canuto et al., 2001

use TURB_CONST, only: &
& CEcoef01, &
& CEcoef02, &
& CEcoef03, &

& CEcoef11, &
& CEcoef12, &
& CEcoef13, &
& CEcoef14, &
& CEcoef15, &
& CEcoef16
implicit none

! Upper limit for stability function for momentum (Galperin et al., 1988)
real(8), parameter :: upper_limit = 0.46d0 
real(8), parameter :: small_number = 1.d-20 
real(8), parameter :: lambda_N_low_limit = - 5.9d0

real(8), intent(in) :: lambda_M 
real(8), intent(in) :: lambda_N

real(8) :: lambda_N1

! Putting lower limit for lambda_N, that is defined from the condition
! d(CE_CANUTO)/d(lambda_N)<0

lambda_N1 = dmax1(lambda_N,lambda_N_low_limit)

CE_CANUTO = &
& (CEcoef01 + CEcoef02*lambda_N1 - CEcoef03*lambda_M)/ &
& (CEcoef11 + CEcoef12*lambda_N1 + CEcoef13*lambda_M + &
&  CEcoef14*lambda_N1**2 + CEcoef15*lambda_N1*lambda_M - & 
&  CEcoef16*lambda_M**2)

CE_CANUTO = dmax1(small_number,CE_CANUTO)
CE_CANUTO = dmin1(upper_limit,CE_CANUTO)

END FUNCTION CE_CANUTO



REAL(8) FUNCTION CEt_CANUTO(lambda_M, lambda_N)

! The function CEt_CANUTO calculates stability function for scalars
! according to Canuto et al., 2001

use TURB_CONST, only: &
& CEtcoef01, &
& CEtcoef02, &
& CEtcoef03, &

& CEcoef11, &
& CEcoef12, &
& CEcoef13, &
& CEcoef14, &
& CEcoef15, &
& CEcoef16

implicit none
! Upper limit for stability function for scalars (Galperin et al., 1988)
real(8), parameter:: upper_limit = 0.61d0
real(8), parameter:: small_number = 1.d-20 
real(8), parameter:: lambda_N_low_limit = - 5.9d0   

real(8), intent(in):: lambda_M
real(8), intent(in):: lambda_N
 
real(8) lambda_N1
 
! Putting lower limit for lambda_N, that is defined from the condition
! d(CEt_CANUTO)/d(lambda_N)<0 
lambda_N1 = dmax1(lambda_N,lambda_N_low_limit)

CEt_CANUTO = &
& (CEtcoef01 + CEtcoef02*lambda_N1 + CEtcoef03*lambda_M) / &
& (CEcoef11 + CEcoef12*lambda_N1 + CEcoef13*lambda_M + &
& CEcoef14*lambda_N1**2 + CEcoef15*lambda_N1*lambda_M - &
& CEcoef16*lambda_M**2)

CEt_CANUTO = dmax1(small_number,CEt_CANUTO)
CEt_CANUTO = dmin1(upper_limit,CEt_CANUTO)
 
END FUNCTION CEt_CANUTO
      
      
REAL(8) FUNCTION SMOMENT_GALPERIN(lambda_N)
 
! The function SMOMENT_GALPERIN calculates the equilibrium stability function
! for momentum according to Galperin et al., 1988
 
use TURB_CONST, only: &
& A1, &
& B1, &
& C1, &
& A2, &
& B2, &
& CL_K_KL_MODEL 

implicit none

! Upper limit for stability function (Galperin et al., 1988) 
! real(8), parameter:: upper_limit = 0.46d0 ! This value is in terms of CE
real(8), parameter:: small_number = 1.d-2
! Below this limit for lambda_N stability function has negative values    
!real(8), parameter:: lambda_N_lower_limit = - 1.8d0 
real(8), parameter:: GH_MAX = (A2*(12.d0*A1 + B1 + 3.d0*B2))**(-1) ! Galperin et al., 1988
real(8), parameter:: GH_MIN = - 0.28d0 ! Galperin et al., 1988
 
real(8), intent(in):: lambda_N

real(8) :: GH
 
GH = - lambda_N ! According to Mellor and Yamada, 1982 GH has the opposite sign to lambda_N
! The scaling of lambda_N according to Mellor and Yamada, 1982 
GH = CL_K_KL_MODEL**2/2.d0*GH
! Implementing limitations Galperin et al., 1988
GH = dmin1(GH,GH_MAX)
GH = dmax1(GH,GH_MIN)
 
SMOMENT_GALPERIN = 1.d0 - 3.d0 * C1 - 6.d0 * A1 / B1 - &
& 3.d0 * A2 * GH * ( (B2 - 3.d0*A2) * (1.d0 - 6.d0*A1 / B1) - &
& 3.d0 * C1 * (B2 + 6.d0 * A1) )
SMOMENT_GALPERIN = SMOMENT_GALPERIN / &
& ( (1.d0 - 3.d0 * A2 * GH * (6.d0 * A1 + B2) ) * &
&   (1.d0 - 9.d0 * A1 * A2 * GH) )
SMOMENT_GALPERIN = SMOMENT_GALPERIN * A1
 
SMOMENT_GALPERIN = dmax1(small_number,SMOMENT_GALPERIN)
!SMOMENT_GALPERIN = dmin1(upper_limit, SMOMENT_GALPERIN)

END FUNCTION SMOMENT_GALPERIN
      
      
REAL(8) FUNCTION SHEAT_GALPERIN(lambda_N)
 
! The function SHEAT_GALPERIN calculates the equilibrium stability function
! for heat according to Galperin et al., 1988
 
use TURB_CONST, only: &
& A1, &
& B1, &
& C1, &
& A2, &
& B2, &
& CL_K_KL_MODEL

implicit none

! Upper limit for stability function (Galperin et al., 1988)
! real(8), parameter:: upper_limit = 0.61d0 ! This value is in terms of CEt
real(8), parameter:: small_number = 1.d-2 
! Below this limit for lambda_N stability function has negative values
! real(8), parameter:: lambda_N_lower_limit = - 1.8d0 
real(8), parameter:: GH_MAX = (A2*(12.d0*A1 + B1 + 3.d0*B2))**(-1) ! Galperin et al., 1988
real(8), parameter:: GH_MIN = - 0.28d0 ! Galperin et al., 1988

real(8), intent(in):: lambda_N
 
real(8) :: GH
 
! lambda_N1 = - dmax1(lambda_N,lambda_N_lower_limit)
! The scaling of lambda_N according to Mellor and Yamada, 1982
GH = - lambda_N ! According to Mellor and Yamada, 1982 GH has the opposite sign to lambda_N
! The scaling of lambda_N according to Mellor and Yamada, 1982 
GH = CL_K_KL_MODEL**2/2.d0*GH
! Implementing limitations Galperin et al., 1988
GH = dmin1(GH,GH_MAX)
GH = dmax1(GH,GH_MIN)
 
SHEAT_GALPERIN = 1.d0 - 6.d0 * A1 / B1
SHEAT_GALPERIN = SHEAT_GALPERIN / &
& (1.d0 - 3.d0 * A2 * GH * (6.d0 * A1 + B2) )
SHEAT_GALPERIN = SHEAT_GALPERIN * A2
 
SHEAT_GALPERIN = dmax1(small_number,SHEAT_GALPERIN)
! SHEAT_GALPERIN = dmin1(upper_limit, SHEAT_GALPERIN)
 
END FUNCTION SHEAT_GALPERIN


SUBROUTINE VSMOP_LAKE(f,fs,nmin,nmax,k0,k1,gama,perdam)

! only the perturbation is smoothed (if perdam is .true.)
!
! 4th-order vertical smoothing. 2nd order for grid-points adjacent
! to the boundary. no smoothing for boundary points.
!
! the second order coeficient is calculated from the fourth order one,
! imposing similar damping for 2 grid-length waves
!
! fourth order damping: response function
! r(lx,ly=)1-16*gama*sin(pi*dx/lx)**4

implicit none

! Input variables
! Integers
integer(4), intent(in) :: nmin
integer(4), intent(in) :: nmax
integer(4), intent(in) :: k0
integer(4), intent(in) :: k1

! Reals
real(8), intent(in) :: fs(nmin:nmax)
real(8), intent(in) :: gama

! Logicals
logical, intent(in) :: perdam

! Input/output variables
! Reals
real(8), intent(inout) :: f(nmin:nmax)

! Local variables
! Reals
real(8) :: work(nmin:nmax)
real(8) :: um6g
real(8) :: gasc
real(8) :: gasc05
real(8) :: umgasc

! Integers
integer(4) :: i ! Loop index

um6g = 1. - 6.*gama
gasc = 8.*gama
gasc05 = 0.5*gasc
umgasc = 1. - gasc

if (perdam) then
  do i = k0, k1
    f(i) = f(i) - fs(i)
  enddo
endif

do i = k0+2, k1-2
  work(i) = um6g*f(i) - gama*(f(i+2) + f(i-2) - 4.*(f(i+1) + f(i-1)) )
enddo

do i = k0+1, k1-1, k1-k0-2
  work(i) = umgasc*f(i) + gasc05*(f(i-1) + f(i+1))
enddo

do i = k0+1, k1-1
  f(i) = work(i)
enddo

if (perdam) then
  do i = k0, k1
    f(i) = f(i) + fs(i)
  enddo
endif

RETURN
END SUBROUTINE VSMOP_LAKE


FUNCTION INTERPOLATE_TO_INTEGER_POINT(fim05,fip05,ddzim1,ddzi)

! Interpolates linearly the grid function 
! from two half-integer points to integer point between them

implicit none

real(8) :: INTERPOLATE_TO_INTEGER_POINT

! Input variables
real(8), intent(in) :: fim05  ! The function at (i-0.5)-th level
real(8), intent(in) :: fip05  ! The function at (i+0.5)-th level
real(8), intent(in) :: ddzim1 ! The thickness of (i-1)-th layer
real(8), intent(in) :: ddzi   ! The thickness of i-th layer

INTERPOLATE_TO_INTEGER_POINT = fim05 + (fip05 - fim05)*ddzim1/(ddzim1 + ddzi)

END FUNCTION INTERPOLATE_TO_INTEGER_POINT
