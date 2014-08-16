       program  gamma
	   implicit none

	   real*8 w_2,w_2s,ust_s,dzits,hbl,wstar,z,rho,cp,h,karm,deltaz,tst_s,difk,dift,gammah,g,thv,dift_tm,difk_tm
	   integer i

	   ust_s=0.3
	   hbl=1000.
	   dzits=-1.
	   h=400.
	   cp=1005.
	   rho=1.4
	   z=200.
       g=9.81
	   thv=246.
	   deltaz=14.
	   karm=0.4
      
	   open(10,file='gamma.txt')
	   open(20,file='fim.txt')

	   do i=14,1000,10

	   z=i


	   tst_s=h/rho/cp/ust_s

	   wstar=(g/thv*h/rho/cp*hbl)**(1./3.)

	   w_2s=1.44*ust_s**2.*(1.-1.5*dzits)**(2./3.)

	w_2=(1.6*ust_s**2.*(1-z/hbl))**(3./2.)+ &
     & 1.2*wstar**3.*z/hbl*(1-0.9*z/hbl)**(3./2.)
	w_2=w_2**(2./3.)
	
	gammah=2.*wstar*h/cp/rho/w_2/hbl
	

	dift=karm*ust_s*deltaz/((1-16.*dzits)**(-0.5)-karm*deltaz/tst_s*gammah)* &
	& ((hbl-z)/(hbl-deltaz))**2.*(ust_s*karm*z+wstar*hbl* &
     & (z/hbl)**(4./3.))/(ust_s*karm*deltaz+ &
     & wstar*hbl*(deltaz/hbl)**(4./3.))
!      dift=wstar(ix,iy)*hbl(ix,iy)*(z/hbl(ix,iy))**(4./3.)
 !    :	*(1.-z/hbl(ix,iy))**2.*(1.+(-0.2)*z/hbl(ix,iy))

	difk=dift*((1-16.*dzits)**(-0.25)+2.*wstar* &
     & ust_s*karm*deltaz/hbl  &
     & /((1-16.*dzits)**(-0.25))/w_2s)


!	 difk_tm=ust_s*(1-7.*dzits)**(1./3.)*karm*z*(1-z/hbl)**2.
    difk_tm=ust_s*(1-16.*dzits)**(0.25)*karm*z*(1-z/hbl)**2.
	   dift_tm=difk_tm/(1./(1.-16.*dzits)**0.25+karm*deltaz/hbl*6.5)

	 write(10,'(2f9.2,f9.4)')dift_tm,difk_tm,z/hbl


	 enddo

	 do i=1,500
	 dzits=-i/100.
	 write(20,'(f6.2,2f9.3)')dzits,(1-7.*dzits)**(-1./3.),(1-16.*dzits)**(-0.25)
	 enddo

	   end