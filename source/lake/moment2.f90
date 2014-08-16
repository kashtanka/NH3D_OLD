SUBROUTINE MOMENT_SOLVER(ix, iy, nx, ny, year, month, day, hour, kor, a_veg, c_veg, h_veg, &
 & alphax, alphay, dt, b0, tau_air, tau_i, tau_gr)

! Subroutine MOMENT_SOLVER solves the momentum equations
! for two horizontal components od speed

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
real(8), intent(in) :: alphax
real(8), intent(in) :: alphay
real(8), intent(in) :: hour

! Integers
integer(4), intent(in) :: ix
integer(4), intent(in) :: iy
integer(4), intent(in) :: nx
integer(4), intent(in) :: ny
integer(4), intent(in) :: year
integer(4), intent(in) :: month
integer(4), intent(in) :: day

! Output variables

real(8), intent(out) :: b0
real(8), intent(out) :: tau_air
real(8), intent(out) :: tau_gr
real(8), intent(out) :: tau_i
        
! Local variables     
! Reals

real(8) :: CE(M)

real(8) :: a(vector_length)
real(8) :: b(vector_length)
real(8) :: c(vector_length)
real(8) :: d(vector_length)

real(8) :: am(vector_length,2,2)
real(8) :: bm(vector_length,2,2)
real(8) :: cm(vector_length,2,2)
real(8) :: ym(vector_length,2)
real(8) :: dm(vector_length,2)

real(8) :: work(1:M+1,1:2) ! Work array

real(8) :: wind10
real(8) :: urel
real(8) :: vrel
real(8) :: u
real(8) :: v
real(8) :: wr

real(8) :: pi
real(8) :: ufr
real(8) :: ACC2
real(8) :: ACCk
real(8) :: dudz2dvdz2
real(8) :: taux
real(8) :: tauy
real(8) :: kor2
real(8) :: Cz_gr
real(8) :: coef1
real(8) :: coef2
real(8) :: rp
real(8) :: Cz_i
real(8) :: tau_sbl
real(8) :: tau_wav
real(8) :: xx
real(8) :: yy
real(8) :: lm
real(8) :: ext_lamw
real(8) :: month_sec
real(8) :: day_sec
real(8) :: hour_sec
real(8) :: AL
real(8) :: zup
real(8) :: urels
real(8) :: vrels
real(8) :: minwind
real(8) :: lambda_M
real(8) :: lambda_N

! Integers
integer(4) :: i
integer(4) :: j
integer(4) :: k

! Logicals
logical :: indstab
logical :: ind_bound
logical :: firstcall

! Characters
character :: tp_name*10
character :: numt*1

common /wind/  wind10,urel,vrel,u,v

data firstcall /.true./

! Externals
real(8), external :: KRON
real(8), external :: DZETA
real(8), external :: CE_CANUTO 
real(8), external :: SMOMENT_GALPERIN 

SAVE

firstcall_if : if (firstcall) then
       
  month_sec   = 30.*24.*60.*60.
  day_sec     = 24*60.*60.
  hour_sec    = 60.*60.

  pi          = 4.*dtan(1.d0)

  AL    = g/row0

! Parameters of numerical scheme
  ACC2    = 1.d-20 !1.d-20
  ACCk    = 1.d-20 !1.d-20
  knum    = 0.
  minwind = 1.d-1 !1.0d0
       
! Friction by vegetation

  if (h_veg>0) then
    print*, 'Vegetation friction parameterization &
    & is not operational currently: STOP'
    STOP
  endif

  veg_sw = 0.
  do i = 1, M+1
    if (h1-DZETA(dble(i))*h1 < h_veg) veg_sw(i) = 1.
  enddo

endif firstcall_if

u=u1(1)
v=v1(1)

if (relwind==2) then
  urel=uwind-u
  vrel=vwind-v
elseif (relwind==1) then
  urel=uwind
  vrel=vwind
else
  print*, 'relwind=',relwind,' - incorrect meaning: STOP' 
  stop
endif

if (dabs(urel)<minwind) then
  urels=sign(minwind,urel)
else
  urels=urel
endif

if (dabs(vrel)<minwind) then
  vrels=sign(minwind,vrel)
else   
  vrels=vrel
endif

kor2  = 0.5d0*kor*dt 

wr   = dsqrt(urels**2+vrels**2)
tau_air = roa0*cdmw2*wr
ufr  = dsqrt(tau_air/row0)    

!if (month == 5 .and. day == 17) then
!  write(*,*) tau_air
!  read*
!endif

! PARAMETERIZATION OF WIND STRESS SPLIT-UP INTO TWO PARTS:
! tau=tau_sbl+tau_wave
! if (wind10>5.) then
!   C10=1.26*10**(-3.)
! else
!   C10=0.0044*wind10**(-1.15)
! endif
! a.PARAMETERIZATION USING SIGNIFICANT WAVE HEIGHT (NB1, pp. 8)
! co1=4.*kws/sqrt(C10*1000.*g) !1000 m is "average fetch" of lake Syrdakh
! kw=min(max(co1*wind10,kws),0.75)
! b.PARAMETERIZATION USING WAVE AGE (agw) (NB1, pp. 9)
! agw=0.14*(g*1000.)**(1./3.)*C10**(1./6.)*wind10**(-2./3.) 
! kw=min(max(kws*1.15/agw,kws),0.75)
! KW=0.75

! tau_wav = tau_air !*kw
tau_sbl = tau_air !*(1.-kw)

! Calculation of Shezy coefficient
rp = 0.1 !height of roughness elements
if (l1/=0) then
  Cz_gr=7.7*dsqrt(g)*(h1/2/rp)**(1./6.) !50. !ground Shezy coefficient
else
  Cz_gr=7.7*dsqrt(g)*(h1/rp)**(1./6.)
endif
rp = 0.01 
Cz_i = 7.7*dsqrt(g)*(h1/2/rp)**(1./6.)

! Friction: ground
 tau_gr = row0*g*Cz_gr**(-2)*(u1(M)**2+v1(M)**2)
! tau_grx=row0*g*Cz_gr**(-2)*u1(M)*dsqrt((u1(M)**2+v1(M)**2)) !*0.005
! tau_gry=row0*g*Cz_gr**(-2)*v1(M)*dsqrt((u1(M)**2+v1(M)**2)) !*0.005

! Friction: ice
 tau_i = row0*g*Cz_i**(-2)*(u1(1)**2+v1(1)**2)
! tau_ix=row0*g*Cz_i**(-2)*u1(1)*dsqrt(u1(1)**2+v1(1)**2)
! tau_iy=row0*g*Cz_i**(-2)*v1(1)*dsqrt(u1(1)**2+v1(1)**2)
coef1 = -g*Cz_i**(-2)*dsqrt(u1(1)**2+v1(1)**2)
coef2 = g*Cz_gr**(-2)*dsqrt(u1(M)**2+v1(M)**2) !*0.005 !*row0
      
! Eddy viscosity parameterization
do i=1,M+1
  uv(i)=dsqrt(u1(i)**2+v1(i)**2)
enddo
do i=1,M
! Ri(i)=dmin1(dmax1(g/row0/(((uv(i+1)-uv(i))/
! & (ddz*h1))**2)*(row(i+1)-row(i))/(ddz*h1),-10.d0),10.d0)
  dudz2dvdz2 = dmax1( ( (u1(i+1)-u1(i))/(ddz(i)*h1) )**2 + &
  & ( (v1(i+1)-v1(i))/(ddz(i)*h1) )**2, 1.d-1)
  Ri(i) = g/row0*(row(i+1)-row(i))/(ddz(i)*h1)/dudz2dvdz2
enddo

SELECT CASE (Turbpar)
       
! 1. "Empirical" parametrization: Stepanenko, Lykossov (2005)
  CASE (1)
    tp_name='SL_Empir' 
    k2(1) =(lamw0*10.+(wind*2./zref/20.)*(lamw0*1000.-lamw0*10.))/ &
    & (cw*row0)
    ext_lamw = dlog(k2(1)*cw*row0/(lamw0*10.))/h1
    do i=2,M+1
      k2(i) = k2(1)*dexp(-ext_lamw*dzeta(dble(i+0.5))*h1)+niu
    enddo
    k2(1)=k2(1)+niu

! 2. "E-epsilon" parameterization: k=CE*E**2/eps with 
!    prognostic equations for E and eps
  CASE (2)
    do i=1,M+1
      if (E1(i)<=0) E1(i)=10.**(-16)
      if (eps1(i)<=0) eps1(i)=10.**(-18)
    enddo
    do i = 1, M
      if (stabfunc == 1) then
        CE(i) = CE0
      else
        lambda_N = E1(i)*E1(i)/((eps1(i)+ACC2)*(eps1(i)+ACC2))* &
        & ((row(i+1)-row(i))/(h1*ddz(i)))*AL
        lambda_M = E1(i)*E1(i)/((eps1(i)+ACC2)*(eps1(i)+ACC2))* &
        & ( (u1(i+1)-u1(i))*(u1(i+1)-u1(i)) + &
        &   (v1(i+1)-v1(i))*(v1(i+1)-v1(i)) ) / &
        &    (h1*h1*ddz(i)*ddz(i))
        if (stabfunc == 2) then
          CE(i)  = CE_CANUTO (lambda_M, lambda_N)
        elseif (stabfunc == 3)  then
          CE(i)  = dsqrt(2.d0)*CL_K_KL_MODEL*SMOMENT_GALPERIN (lambda_N)
        endif
      endif
    enddo
    do i = 1, M
      k2(i)=dmax1(CE(i)*E1(i)**2/(eps1(i)+ACCk),1.d-8)+niu+knum(i)
    end do

! 3. Nickuradze (NICK) formulation: Rodi (1993)
  CASE (3)
    tp_name='NICK'
    do i=1,M
      zup=h1-dzeta(dble(i+0.5))*h1
      lm=h1*(0.14-0.08*(1-zup/h1)**2-0.06*(1-zup/h1)**4)
      k2(i)=lm**2*dabs((uv(i+1)-uv(i))/(ddz(i)*h1))*dexp(-Cs*Ri(i)) &
      & + niu
    enddo
      
! 4. Parabolic (PARAB) formulation: Engelund (1976)
  CASE (4)
    tp_name='PARAB'
    do i=1,M
      zup=h1-dzeta(dble(i+0.5))*h1
      k2(i)=kar*ufr*zup*(1-zup/h1)*dexp(-Cs*Ri(i)) &
      & +niu
    enddo 

! 5. W2 (used in Version 2 of CE-QUAL-W2 model): Cole and Buchak (1995)   
  CASE (5)
    print*, 'Turbpar = 5 is not operational setting: STOP'
    STOP

! 6. W2N (W2 with mixing length of Nickuradze): Cole and Buchak (1995) and Rodi (1993)
  CASE (6)
    print*, 'Turbpar = 6 is not operational setting: STOP'
    STOP
   
! 7. RNG (re-normalization group) formulation: Simoes (1998)
  CASE (7)
    tp_name='RNG'
    do i=1,M
      zup=h1-dzeta(dble(i+0.5))*h1
      k2(i)=niu*(1+dmax1(3*kar*(zup/niu*ufr)**3* &
      & (1-zup/h1)**3-Cs1,0.d0))**(1./3.)*dexp(-Cs*Ri(i))+niu  
    enddo
END SELECT
        
      
! Equations for horizontal velocities (u and v)

xx=0.5*(1.-dhw0*h1*ddz(1)/(k2(1)*dt))
yy=(ddz(1)*h1)**2/(k2(1)*dt)

! The top boundary conditions for u and v
! 1-st case: water surface is free from ice:
! momentum flux = momentum flux from atmosphere

if (l1==0.d0) then
  taux=urels/wr*tau_sbl/row0
  tauy=vrels/wr*tau_sbl/row0
  do i=1,2
    do j=1,2
      cm(1,i,j)=KRON(i,j)*(xx+yy+ & 
      & 0.5*yy*dt*a_veg*c_veg*dsqrt(u1(1)**2+v1(1)**2)*veg_sw(1))
      bm(1,i,j)=KRON(i,j)*xx
    enddo
  enddo

  cm(1,1,2)=-kor*(h1*ddz(1))**2/(2*k2(1))
  cm(1,2,1)=kor*(h1*ddz(1))**2/(2*k2(1))

  dm(1,1)=(taux*ddz(1)*h1/k2(1)+yy*u1(1)) 
  dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx + &
  & kor*(h1*ddz(1))**2*v1(1)/(2*k2(1)) + &
  & yy*dt*g*dtan(pi*alphax/180.)
  dm(1,2)=(tauy*ddz(1)*h1/k2(1)+yy*v1(1))
  dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx - &
  & kor*(h1*ddz(1))**2*u1(1)/(2*k2(1)) + &
  & yy*dt*g*dtan(pi*alphay/180.)
endif

! 2-nd case: thin ice:
! momentum flux = 
! weighted momentum flux from atmosphere plus weighted friction 
! at the water-ice interface

if (l1>0.and.l1<=L0) then 
  b0=l1/L0
  taux=urels/wr*tau_sbl/row0 !wdir+pi
  tauy=vrels/wr*tau_sbl/row0 !wdir+pi
  do i=1,2
    do j=1,2
      cm(1,i,j)=KRON(i,j)*(xx+yy-b0*coef1*ddz(1)*h1/k2(1) + &
      & 0.5*yy*dt*a_veg*c_veg*dsqrt(u1(1)**2+v1(1)**2)*veg_sw(1))
      bm(1,i,j)=KRON(i,j)*xx
    enddo
  enddo

  cm(1,1,2)=-kor*(h1*ddz(1))**2/(2*k2(1))
  cm(1,2,1)=kor*(h1*ddz(1))**2/(2*k2(1))

  dm(1,1)=yy*u1(1) !-tau_ix*ddz*h1/k2(1)
  dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx+2*taux*ddz(1)*h1/k2(1)*(1-b0) + &
  & kor*(h1*ddz(1))**2*v1(1)/(2*k2(1))+yy*dt*g*dtan(pi*alphax/180.)
  dm(1,2)=yy*v1(1) !-tau_iy*ddz*h1/k2(1)
  dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx+2*tauy*ddz(1)*h1/k2(1)*(1-b0) - &
  & kor*(h1*ddz(1))**2*u1(1)/(2*k2(1))+yy*dt*g*dtan(pi*alphay/180.)
endif

! 3-rd case: thick ice:
! momentum flux = water-ice friction 

if (l1>L0) then
  do i=1,2
    do j=1,2
      cm(1,i,j)=kron(i,j)*(xx+yy-coef1*ddz(1)*h1/k2(1) + &
      & 0.5*yy*dt*a_veg*c_veg*dsqrt(u1(1)**2+v1(1)**2)*veg_sw(1))
      bm(1,i,j)=kron(i,j)*xx
   enddo
  enddo

  cm(1,1,2)=-kor*(h1*ddz(1))**2/(2*k2(1))
  cm(1,2,1)=kor*(h1*ddz(1))**2/(2*k2(1))

  dm(1,1)=yy*u1(1) !-tau_ix*ddz*h1/k2(1)
  dm(1,1)=dm(1,1)+(u1(2)-u1(1))*xx + &
  & kor*(h1*ddz(1))**2*v1(1)/(2*k2(1)) + &
  & yy*dt*g*tan(pi*alphax/180.)
  dm(1,2)=yy*v1(1) !-tau_iy*ddz*h1/k2(1)
  dm(1,2)=dm(1,2)+(v1(2)-v1(1))*xx - &
  & kor*(h1*ddz(1))**2*u1(1)/(2*k2(1)) + &
  & yy*dt*g*tan(pi*alphay/180.)
endif

! Boundary conditions for u and v at the lake bottom

xx=0.5*(1-ddz(M)*h1*(dhw-dhw0)/(k2(M)*dt))
yy=(ddz(M)*h1)**2/(k2(M)*dt)
do i=1,2
  do j=1,2
    cm(M+1,i,j)=KRON(i,j)*(xx+yy+coef2*ddz(M)*h1/k2(M) + &
    & 0.5*yy*dt*a_veg*c_veg*dsqrt(u1(M+1)**2+v1(M+1)**2)*veg_sw(M+1))
    am(M+1,i,j)=KRON(i,j)*xx
  enddo
enddo

cm(M+1,1,2)=-kor*(h1*ddz(M))**2/(2*k2(M))
cm(M+1,2,1)=kor*(h1*ddz(M))**2/(2*k2(M))

dm(M+1,1)=yy*u1(M+1) !-tau_grx*ddz*h1/k2(M)
dm(M+1,1)=dm(M+1,1)-(u1(M+1)-u1(M))*xx + &
& kor*(h1*ddz(M))**2*v1(M+1)/(2*k2(M)) + &
& yy*dt*g*tan(pi*alphax/180.)
dm(M+1,2)=yy*v1(M+1) !-tau_gry*ddz*h1/k2(M)
dm(M+1,2)=dm(M+1,2)-(v1(M+1)-v1(M))*xx - &
& kor*(h1*ddz(M))**2*u1(M+1)/(2*k2(M)) + &
& yy*dt*g*tan(pi*alphay/180.)

! The coefficients of equation for u and v in the internal points
! of the mesh

do k=2,M
  do i=1,2
    do j=1,2
      am(k,i,j)=KRON(i,j)*(dzeta(dble(k))*dhw/(2.*(ddz(k-1)+ddz(k))*h1) - &
      & dhw0/(2.*(ddz(k-1)+ddz(k))*h1) - &
      & k2(k-1)*dt/(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) ))
      bm(k,i,j)=-KRON(i,j)*(dzeta(dble(k))*dhw/(2*(ddz(k-1)+ddz(k))*h1) + &
      & dhw0/(2*(ddz(k-1)+ddz(k))*h1) + &
      & k2(k)*dt/(h1**2*(ddz(k)+ddz(k-1))*ddz(k) ))
    enddo
  enddo

  cm(k,1,1)=-(k2(k)/ddz(k)+k2(k-1)/ddz(k-1))*dt/ &
  & (h1**2*(ddz(k)+ddz(k-1)) ) - 1. - &
  & 0.5*dt*c_veg*a_veg*dsqrt(u1(k)**2+v1(k)**2)*veg_sw(k)
  cm(k,1,2)=kor2 !Krank-Nickolson
  cm(k,2,1)=-cm(k,1,2)
  cm(k,2,2)=cm(k,1,1)

  dm(k,1)=-u1(k)-kor2*v1(k)-dt*g*tan(pi*alphax/180.) !Krank-Nickolson
  dm(k,1)=dm(k,1)-dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k) )**(-1)* &
  & (k2(k)*(u1(k+1)-u1(k)))
  dm(k,1)=dm(k,1)+dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) )**(-1)* &
  & (k2(k-1)*(u1(k)-u1(k-1)))
  dm(k,1)=dm(k,1)-dzeta(dble(k))*dhw*(2*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(u1(k+1)-u1(k-1))
  dm(k,1)=dm(k,1)+dhw0*(2.*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(u1(k+1)-u1(k-1))

  dm(k,2)=-v1(k)+kor2*u1(k)-dt*g*tan(pi*alphay/180.) !Krank-Nickolson
  dm(k,2)=dm(k,2)-dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k) )**(-1)* &
  & (k2(k)*(v1(k+1)-v1(k)))
  dm(k,2)=dm(k,2)+dt*(h1**2*(ddz(k)+ddz(k-1))*ddz(k-1) )**(-1)* &
  & (k2(k-1)*(v1(k)-v1(k-1)))
  dm(k,2)=dm(k,2)-dzeta(dble(k))*dhw*(2*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(v1(k+1)-v1(k-1))
  dm(k,2)=dm(k,2)+dhw0*(2.*(ddz(k)+ddz(k-1) )*h1)** &
  & (-1)*(v1(k+1)-v1(k-1))
enddo

! ind_bound=.false.
! call ind_stab_fact_db (a,b,c,1,M+1,indstab,ind_bound)
! if (indstab==.false.) then
! do i=2,M
!   a(i)=-k2(i)*dt/(h1**2*ddz**2)
!   b(i)=-k2(i+1)*dt/(h1**2*ddz**2)
! enddo
! endif  

call MATRIXPROGONKA(am,bm,cm,dm,ym,M+1)
do k=1,M+1
  u2(k)=ym(k,1)
  v2(k)=ym(k,2)
enddo

!print*, '(us, vs)', u2(1), v2(1)
     
if (firstcall) firstcall=.false.

END SUBROUTINE MOMENT_SOLVER