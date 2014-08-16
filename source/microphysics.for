      subroutine iceini
c--------------------------------------------------------------------------------------------------------------c      
c
c      initiation of ice crystals   (water vapor => cloud ice)                            DC,10.2009
c      after Dudhia,1989
c---------------------------------------------------------------------------------------------------------------c
	use alloc    
      
      implicit none
	integer ix,iy,is
	real(8) xMo,p,t,xNice
	real(8),external::qsati,qsat
	
	!xMo=1.e-12
	xMo=1.e-12
      
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
            p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	      if (t.le.273.15) then
	      xNice=min(10.**5,0.01*exp(0.6*(273.15-t)))                 !ice crystals concentration after Fletcher
	      else
	      xNice=0.
	      endif
	      if(t.le.273.15) then
	      vini(ix,iy,is)=max((xMo*xNice-qci(ix,iy,is,2))/dtl,0.)
	      else
	      vini(ix,iy,is)=0.
	      endif
		  if (qv(ix,iy,is,3)-qs(ix,iy,is).gt.0) then
            vini(ix,iy,is)=min(vini(ix,iy,is),
     :	  (qv(ix,iy,is,3)-qs(ix,iy,is))/dtl)
	      else
	      vini(ix,iy,is)=0.
	      endif

	    enddo
	  enddo
	enddo
	write(0,*) 'iceini done', sum(vini)
      return
	end


	subroutine icemelt
c-----------------------------------------------------------------------------------------------------------------c
c     cloud ice melting rate if T>0 Celcium  (cloud ice => cloud water)
c
c-----------------------------------------------------------------------------------------------------------------c
      use alloc
	implicit none
	integer ix,iy,is
	real*8 p,t    
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if (t.gt.273.15) then
	    imlt(ix,iy,is)=max(0.,qci(ix,iy,is,2))/dtl
	    else
	    imlt(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'icemelt done', sum(imlt)
	return
	end

	subroutine homofreeze
c-------------------------------------------------------------------------------------------------------------------c
c      homogeneous freezing of cloud water at T < Too = -40. Celcium
c      cloud water => cloud ice
c-------------------------------------------------------------------------------------------------------------------c
      use alloc
	implicit none
	integer ix,iy,is
	real*8 p,t,Too
	Too=233.15    
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if(t.lt.Too) then
	    hmfrz(ix,iy,is)=max(0.,qc(ix,iy,is,2))/dtl
	    else
          hmfrz(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'homofreeze done',sum(hmfrz)
	return
	end

	subroutine accretcwsn
c------------------------------------------------------------------------------------------------------------------c
c       calculation of accretion rate of cloud water by snow  (after Lin et al,1983)
c       (cloud water => snow)
c       if t>0 then it contributes to rain water growth 
c------------------------------------------------------------------------------------------------------------------c
      use alloc
	implicit none
	integer ix,iy,is
	real*8 p,t,gamma
	real(8),external::lamsnow,gamf
	gamma=gamf(3.+dsnow)
	do is=0,ns
        do iy=0,ny1
        do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if (qsn(ix,iy,is,2).gt.0) then
	    sacrw(ix,iy,is)=pi*nso*csnow*qc(ix,iy,is,2)*gamma
     :	    *sqrt(rho0/rho(ix,iy,is))/
     :      (4.*lamsnow(rho(ix,iy,is),qsn(ix,iy,is,2))**(3.+dsnow))
	    else
		sacrw(ix,iy,is)=0. 
	    endif
	    if(t.ge.273.15.and.qc(ix,iy,is,2).gt.0.) then
          sacrwr(ix,iy,is)=min(sacrw(ix,iy,is),qc(ix,iy,is,2)/dtl)
	    else
          sacrwr(ix,iy,is)=0.
	    endif
	    if (t.lt.273.15.and.qc(ix,iy,is,2).gt.0.) then
	    sacrw(ix,iy,is)=min(sacrw(ix,iy,is),qc(ix,iy,is,2)/dtl)
	    else
	    sacrw(ix,iy,is)=0.
	    endif
	enddo
	enddo
	enddo
	write(0,*) 'accretcwsn done',sum(sacrw),sum(sacrwr)
	return
	end

	subroutine bergeron
c--------------------------------------------------------------------------------------------------------------------c
c        bergeron process (deposition and riming)  after Lin et al, 1983, Hsie et al, 1980 
c        transfer of cloud water to form snow     (sbercw)
c        and transfer from cloud ice to snow      (sberci)
c        and transfer from cloud water to cloud ice (ibercw)
c---------------------------------------------------------------------------------------------------------------------c

      use alloc

      implicit none

	real*8 mi50,ri50,ui50,ni50,qi50,dt1,Ein,p,t,mn,
     : xNice,mi40,ri60,mi60,ni60,qi60,ui60,mn2
      real(8), external :: a1koenig,a2koenig
	integer ix,iy,is

	mi50=4.8e-10
	mi40=2.46e-10
	ri50=50.e-6
	ri60=60.e-6
	ui50=1.
	ui60=1.
        Ein=1.
	mn=1.05e-18
	mi60=4.52e-10
	mi40=1.34e-10
      
	
	do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if (t.lt.273.15) then
	    dt1=1/(a1koenig(t)*(1.-a2koenig(t)))*
     :	(mi60**(1.-a2koenig(t))-mi40**(1.-a2koenig(t)))
	    else
	    dt1=0.
	    qi60=0.
	    endif
         if(dt1.gt.200) dt1=10.     !  dt1 should be chosen so that sberci*dt = 0.05-0.1*qci, which gives dt1=10...20
	    qi60=qci(ix,iy,is,2)*dtl/dt1
	    ni60=qi60/mi60
	    sbercw(ix,iy,is)=ni60*(a1koenig(t)*mi60**a2koenig(t)+
     :    pi*Ein*rho(ix,iy,is)*qc(ix,iy,is,2)*ri60**2.*ui60)
          sbercw(ix,iy,is)=min(sbercw(ix,iy,is),qc(ix,iy,is,2)/dtl)
		if (qc(ix,iy,is,2).gt.0.and.qci(ix,iy,is,2).gt.0) then		      
	    sberci(ix,iy,is)=qci(ix,iy,is,2)/dt1
	    sberci(ix,iy,is)=min(qci(ix,iy,is,2)/dtl,sberci(ix,iy,is))
	    else
	    sberci(ix,iy,is)=0.
	    endif
	    if (t.lt.273.15) then 
		  xNice=min(10.**5,0.01*exp(0.6*(273.15-t)))
	!      xNice=1000.*exp(12.96*Sice-0.639)
	      else
	      xNice=0.
	      endif
          if(qci(ix,iy,is,2).eq.0) then
		mn2=mn
		else  
	    mn2=qci(ix,iy,is,2)*rho(ix,iy,is)/xNice
	    endif
	    if(t.lt.273.15.and.qc(ix,iy,is,2).gt.0.
     :	and.qci(ix,iy,is,2).gt.0) then
	    ibercw(ix,iy,is)=xNice/(1000.*rho(ix,iy,is))
     :      *a1koenig(t)*mn**a2koenig(t)
	    ibercw(ix,iy,is)=min(ibercw(ix,iy,is),qc(ix,iy,is,2)/dtl)
	    else
	    ibercw(ix,iy,is)=0.
	    endif
	    if (t.ge.273.15) then
	    sbercw(ix,iy,is)=0.
	    sberci(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'begeron done',sum(sberci),sum(sbercw),sum(ibercw)
	return
	end

	subroutine saggregation
c------------------------------------------------------------------------------------------------------------------c
c         ice crystal aggregation rate to form snow 
c         after Lin,1983 and Kessler,1969 
c         cloud ice => snow
c------------------------------------------------------------------------------------------------------------------c
       use alloc

	implicit none
	real(8) alfa,qci0,t,p,t0
	integer ix,iy,is
        t0=273.15
	!qci0=1.e-3          !Lin,1983
	qci0=1.2e-4        !BEM
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
          p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    alfa=1.e-3*exp(0.025*(t-t0))
	    sagg(ix,iy,is)=alfa*(qci(ix,iy,is,2)-qci0)
	    if(t.lt.t0.and.qci(ix,iy,is,2).gt.0.) then
	    sagg(ix,iy,is)=max(0.,sagg(ix,iy,is))
	    else
	    sagg(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'saggregation done',sum(sagg)
	return
	end

	subroutine accretcisn
c------------------------------------------------------------------------------------------------------------------c
c         accretion rate of cloud ice by snow
c         after Lin,1983
c         cloud ice => snow
c------------------------------------------------------------------------------------------------------------------c
      use alloc

	implicit none
	real(8) p,t,t0,esi,gamma
	integer ix,iy,is
	real(8), external :: gamf,lamsnow
	gamma=gamf(dsnow+3.)
	t0=273.15
	do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    esi=exp(0.025*(t-t0))
	    if (qsn(ix,iy,is,2).gt.0) then
	    saci(ix,iy,is)=pi*esi*nso*csnow*qci(ix,iy,is,2)*gamma
     :	    *sqrt(rho0/rho(ix,iy,is))*0.25
     :      /lamsnow(rho(ix,iy,is),qsn(ix,iy,is,2))**(dsnow+3.)
	    else
	    saci(ix,iy,is)=0.
	    endif     
	    if(t.lt.t0)
     :    then
	    saci(ix,iy,is)=min(saci(ix,iy,is),qci(ix,iy,is,2)/dtl)
	    else
	    saci(ix,iy,is)=0.
	    endif

          enddo
	  enddo
	enddo
      write(0,*) 'accretcisn done',sum(saci)
	return
	end

	subroutine accret3comp
c------------------------------------------------------------------------------------------------------------------c
c       Three-component process: cloud ice accreted by rain water to form snow
c       raci - sink term for cloud ice
c       iacr - sink term for rain
c       raci+iacr - source term for snow
c       after Lin, 1983
c
c------------------------------------------------------------------------------------------------------------------c
      use alloc

	implicit none

	real(8) p,t,eri,Mi,gamma1,gamma2
	integer ix,iy,is
	real(8), external :: gamf,lamsnow,lamrain
	gamma1=gamf(3.+xb)
	gamma2=gamf(6.+xb)
	eri=1.
	mi=4.19e-13
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if(qr(ix,iy,is,2).gt.0) then
	    raci(ix,iy,is)=pi*eri*xnor*xa*qci(ix,iy,is,2)*gamma1
     :	    *sqrt(rho0/rho(ix,iy,is))*0.25
     :      /lamrain(rho(ix,iy,is),qr(ix,iy,is,2))**(3.+xb)
	    else
	    raci(ix,iy,is)=0.
	    endif
	    if(t.lt.273.15) then
	    raci(ix,iy,is)=min(raci(ix,iy,is),qci(ix,iy,is,2)/dtl)
	    else
	    raci(ix,iy,is)=0.
	    endif
	    if(qr(ix,iy,is,2).gt.0) then
	    iacr(ix,iy,is)=pi**2*eri*xnor*xa*qci(ix,iy,is,2)*rhol*
     :	    gamma2*sqrt(rho0/rho(ix,iy,is))/
     :      (24.*Mi*lamrain(rho(ix,iy,is),qr(ix,iy,is,2))**(6.+xb))
	    else
	    iacr(ix,iy,is)=0.
	    endif
	    if(t.lt.273.15) then
	    iacr(ix,iy,is)=min(iacr(ix,iy,is),qr(ix,iy,is,2)/dtl)
	    else
	    iacr(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'accret3comp done',sum(raci),sum(iacr)
	return
	end

	subroutine accretsr
c------------------------------------------------------------------------------------------------------------------c
c      accretion rate of snow at the expense of rain water 
c      after Lin et al, 1983
c------------------------------------------------------------------------------------------------------------------c
      use alloc

	implicit none
	integer is,iy,ix
	real(8),external :: lamsnow,lamrain,Usn,Urn,gamf
	real(8) ls,lr,p,t,gamma1,gamma2
	gamma1=gamf(4.+dsnow)
	gamma2=gamf(4.+brain)
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if(qr(ix,iy,is,2).gt.0.and.qsn(ix,iy,is,2).gt.0) then
	    lr=lamrain(rho(ix,iy,is),qr(ix,iy,is,2))
	    ls=lamsnow(rho(ix,iy,is),qsn(ix,iy,is,2))
	    sacrr(ix,iy,is)=pi**2*nso*xnor
     :      *abs(Usn(rho(ix,iy,is),qsn(ix,iy,is,2),
     :	    gamma1)-Urn(rho(ix,iy,is),qr(ix,iy,is,2),gamma2))
     :      *(rhol/rho(ix,iy,is))*
     :      (5./(lr**6.*ls)+2./(lr**5.*ls**2.)+0.5/(lr**4.*ls**3.))
	    else
	    sacrr(ix,iy,is)=0.
	    endif
	    if(t.lt.273.15) then
	    sacrr(ix,iy,is)=min(sacrr(ix,iy,is),qr(ix,iy,is,2)/dtl)
	    else
	    sacrr(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
	write(0,*) 'accretsr done',sum(sacrr)
	return
	end
	

	subroutine snowmelt
c------------------------------------------------------------------------------------------------------------------c
c          melting rate of snow to form rain 
c          after Lin et al, 1983
c------------------------------------------------------------------------------------------------------------------c
      use alloc

	implicit none
	real(8) p,t,lfus,xmyu,xKa,diffus,cw,Sc,gamma
	integer ix,iy,is
	real(8), external :: gamf,lamsnow,lvap,qsati
	gamma=gamf((dsnow+5.)/2.)
	lfus=3.336e5
	cw=4.187e3
	do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    diffus=8.794/(10.**5)*(t**1.81)/p
	    xmyu=1.496/(10.**6)*(t**1.5)/(t+120.)
	    xKa=1.414*10.**3*xmyu
	    Sc=xmyu/(rho(ix,iy,is)*diffus)
	    if(qsn(ix,iy,is,2).gt.0) then
	    smlt(ix,iy,is)=-((-2.*pi/rho(ix,iy,is)/lfus)
     :      *(xKa*(t-273.15)-hlat*diffus
     :	    *rho(ix,iy,is)*(qsati(t,p)-qv(ix,iy,is,2)))*nso*(0.78/
     :      lamsnow(rho(ix,iy,is),qsn(ix,iy,is,2))**2
     :      +0.31*Sc**(1./3.)
     :      *gamma*sqrt(csnow)*(rho0/rho(ix,iy,is))**0.25/
     :      sqrt(xmyu/rho(ix,iy,is))/
     :      lamsnow(rho(ix,iy,is),qsn(ix,iy,is,2))**((dsnow+5.)/2.))
     :      -cw*(t-273.15)/lfus*(sacrw(ix,iy,is)+sacrr(ix,iy,is)))
	    else
	    smlt(ix,iy,is)=0.
	    endif
	    if (t.gt.273.15) then
	    smlt(ix,iy,is)=min(smlt(ix,iy,is),qsn(ix,iy,is,2)/dtl)
	    else
	    smlt(ix,iy,is)=0.
	    endif
	    if(smlt(ix,iy,is).lt.0) smlt(ix,iy,is)=0.
	    enddo
	  enddo
	enddo
      write(0,*) 'snowmelt done',sum(smlt)
      return
	end

      subroutine snowprec
c------------------------------------------------------------------------------c
c           preciptating snow
c           after Lin et al, 1983
c------------------------------------------------------------------------------c
	use alloc
	implicit none
	real(8) p,t,gamma
	real(8), external:: Usn,gamf
	integer ix,iy,is
	gamma=gamf(4.+dsnow)
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	    p=sigma0(is)*pp(ix,iy,2)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
	    if(qsn(ix,iy,is,2).gt.0) then
	    vsnow(ix,iy,is)=Usn(rho(ix,iy,is),qsn(ix,iy,is,2),gamma)
	    else
	    vsnow(ix,iy,is)=0.
	    endif
	    enddo
	  enddo
	enddo
      
	is=ns
      do iy=0,ny1
        do ix=0,nx1
          p=sigma0(is)*pp(ix,iy,2)+ptop
          t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
          precsnow(ix,iy)=(0.5*(vsnow(ix,iy,is)
     :      +vsnow(ix,iy,is-1))-w(ix,iy,is,2))
     :      *rho(ix,iy,is)*qsn(ix,iy,is,2)
        enddo
      enddo
      snowacc=snowacc+precsnow
      write(0,*) 'snowprec done'
	return
	end
!---------------------------------------------------------------------!
!      SATURATION ADJUSTMENT after W.-K. TAO, 1989
!        deposition/sublimation   (qv <---> qci)
!    and condensation/evaporation (qv <---> qc)
!        if not enough of qci and qc to remove subsaturation,
!        sublimation of snow takes place (qsn ---> qv)
!---------------------------------------------------------------------!

      subroutine adjustment
      use alloc
	implicit none
	real(8) p,t,CND,DEP,Too,To,dq,dqi,dqc
	real(8) r1,r2,A1,A2,A3,th,dqs
      integer ix,iy,is
	real,external:: qsat,qsati
	Too=233.15
	To=273.15
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	      
	      p=sigma0(is)*pp(ix,iy,3)+ptop
              t=(pts(ix,iy,is)+pt(ix,iy,is,3))*(p/p00)**akapa
	      A1=(237.3*17.27*(p/p00)**akapa)/(t-35.5)**2
	      A2=(265.5*21.88*(p/p00)**akapa)/(t-7.5)**2

	      if (qc(ix,iy,is,3).ne.0.or.qci(ix,iy,is,3).ne.0) then
	      r2=(A1*qc(ix,iy,is,3)*qsat(t,p)+
     :	  A2*qci(ix,iy,is,3)*qsati(t,p))/
     :      (qc(ix,iy,is,3)+qci(ix,iy,is,3))
	      if (t.gt.Too.and.t.le.To) then
            CND=(t-Too)/(To-Too)
	      DEP=(To-t)/(To-Too)
	      elseif (t.gt.To) then
	      CND=1.
	      DEP=0.
	      elseif (t.le.Too) then
	      CND=0.
	      DEP=1.
	      endif
	      endif
	      if(qc(ix,iy,is,3).eq.0.and.qci(ix,iy,is,3).eq.0) then
	      if (t.gt.To) then 
		  r2=A1*qsat(t,p)
            CND=1.
	      DEP=0.
	      endif
	      if (t.le.To.and.t.gt.Too) then 
		  r2=0.5*(A2*qsati(t,p)+A1*qsat(t,p))
            CND=(t-Too)/(To-Too)
	      DEP=(To-t)/(To-Too)
	      endif
            if (t.le.Too) then
            r2=A2*qsati(t,p)
            CND=0.
	      DEP=1.
	      endif
	      endif 


	     
	      A3=(hlat*CND+hsub*DEP)/(cp*(p/p00)**akapa)
	      r1=qv(ix,iy,is,3)-qs(ix,iy,is)
	      
		  

	      if(r1.gt.0) then
	      dq=-r1/(1+r2*A3)
            qv(ix,iy,is,3)=qv(ix,iy,is,3)+dq
		  endif

	      if (r1.lt.0) then
	      if(qc(ix,iy,is,3).gt.0.and.qci(ix,iy,is,3).gt.0) then
	      dq=min(-r1/(1+r2*A3),
     :	  qc(ix,iy,is,3)+qci(ix,iy,is,3))
	      qv(ix,iy,is,3)=qv(ix,iy,is,3)+dq
	      elseif (qc(ix,iy,is,3).eq.0.and.qci(ix,iy,is,3).gt.0) then
	      dq=min(-r1/(1+r2*A3),qci(ix,iy,is,3))
            qv(ix,iy,is,3)=qv(ix,iy,is,3)+dq
	      CND=0.
	      DEP=1.
	      elseif (qc(ix,iy,is,3).gt.0.and.qci(ix,iy,is,3).eq.0) then
	      dq=min(-r1/(1+r2*A3),qc(ix,iy,is,3))
            qv(ix,iy,is,3)=qv(ix,iy,is,3)+dq
	      DEP=0.
	      CND=1.
	      elseif(qc(ix,iy,is,3).eq.0.and.qci(ix,iy,is,3).eq.0) then
	      dq=0.
	      endif
	      endif

	      if (r1.lt.0) then
	      if (dq.lt.-r1/(1+r2*A3).and.qsn(ix,iy,is,3).gt.0.) then
	      dqs=-r1/(1+r2*A3)-dq
		  dqs=min(dqs,qsn(ix,iy,is,3))
            qv(ix,iy,is,3)=qv(ix,iy,is,3)+dqs
            qsn(ix,iy,is,3)=qsn(ix,iy,is,3)-dqs
	      else
	      dqs=0.
	      endif
	      else 
	      dqs=0.
	      endif
	
	      dqc=dq*CND
	      dqi=dq*DEP

	     if(dqi.gt.0.) dqi=min(dqi,qci(ix,iy,is,3))
	     if(dqc.gt.0.) dqc=min(dqc,qc(ix,iy,is,3))
	
	      qci(ix,iy,is,3)=qci(ix,iy,is,3)-dqi
	      qc(ix,iy,is,3)=qc(ix,iy,is,3)-dqc	
	      if(dq.ne.0) then
	      pt(ix,iy,is,3)=pt(ix,iy,is,3)-hsub/cp*(p00/p)**akapa*
     :	  (dqi+dqs)-hlatcp*(p00/p)**akapa*dqc
		  endif  
          enddo
        enddo
      enddo
	return
	end

	subroutine microphys_update
	use alloc
	implicit none
	real(8) condens,transform,p,dqc,dqr,dqci,dqsn
	integer is,iy,ix
	do is=0,ns
        do iy=0,ny1
          do ix=0,nx1

      if (ifqi.ne.0.) then
      dqc=dtl*(auto(ix,iy,is)
     :    +col(ix,iy,is)-imlt(ix,iy,is)+hmfrz(ix,iy,is)                                       
     :    +sacrw(ix,iy,is)+sacrwr(ix,iy,is)+sbercw(ix,iy,is)                                   
     :    +ibercw(ix,iy,is))
      dqr=dtl*(-auto(ix,iy,is)-col(ix,iy,is)
     :      -smlt(ix,iy,is)-sacrwr(ix,iy,is)+sacrr(ix,iy,is)+             
     :      iacr(ix,iy,is)                                                  
     :      +evap(ix,iy,is))
	dqci=dtl*(-vini(ix,iy,is)-hmfrz(ix,iy,is)
     :    -ibercw(ix,iy,is)+imlt(ix,iy,is)+sberci(ix,iy,is)
     :    +sagg(ix,iy,is)+saci(ix,iy,is)+raci(ix,iy,is))
	dqsn=dtl*(
     :      -sbercw(ix,iy,is)-sberci(ix,iy,is)-                  
     :      sacrw(ix,iy,is)-sagg(ix,iy,is)-saci(ix,iy,is)-                     
     :      raci(ix,iy,is)-iacr(ix,iy,is)-sacrr(ix,iy,is)+                     
     :      smlt(ix,iy,is))
       if(dqc.gt.0.and.dqc.gt.qc(ix,iy,is,3)) dqc=qc(ix,iy,is,3)
	if(dqr.gt.0.and.dqr.gt.qr(ix,iy,is,3)) dqr=qr(ix,iy,is,3)
      if(dqci.gt.0.and.dqci.gt.qci(ix,iy,is,3)) dqci=qci(ix,iy,is,3)
	if(dqsn.gt.0.and.dqsn.gt.qsn(ix,iy,is,3)) dqsn=qsn(ix,iy,is,3) 
	
	elseif(ifqi.eq.0.and.ifqr.ne.0) then
      dqc=dtl*(-cond(ix,iy,is)+auto(ix,iy,is)+col(ix,iy,is))
      dqr=dtl*(-auto(ix,iy,is)-col(ix,iy,is)
     :      +evap(ix,iy,is))
       if(dqc.gt.0.and.dqc.gt.qc(ix,iy,is,3)) dqc=qc(ix,iy,is,3)
	if(dqr.gt.0.and.dqr.gt.qr(ix,iy,is,3)) dqr=qr(ix,iy,is,3)
 
      elseif(ifqi.eq.0.and.ifqr.eq.0) then
      dqc=dtl*(-cond(ix,iy,is))

	if(dqc.gt.0.and.dqc.gt.qc(ix,iy,is,3)) dqc=qc(ix,iy,is,3)
	endif

	


	p=sigma0(is)*pp(ix,iy,3)+ptop
	! --------potential temperature-------------!
      if(ifqc.ne.0.and.ifqi.eq.0.and.ifqr.ne.0) then
             condens=hlatcp*(p00/p)**akapa*(cond(ix,iy,is)
     :         -evap(ix,iy,is))
      elseif (ifqc.ne.0.and.ifqi.eq.0.and.ifqr.eq.0) then
       condens=hlatcp*(p00/p)**akapa*(cond(ix,iy,is))
           endif
	if(ifqc.ne.0.and.ifqi.ne.0) then
  !           condens=hlatcp*(p00/p)**akapa*(
  !   :         -evap(ix,iy,is))
           endif
       if(ifqc.eq.0) then
	condens=0.
	endif


	     if(ifqi.ne.0.) then
             transform=
     :         hsub/cp*(p00/p)**akapa*(vini(ix,iy,is))
     :         +hfus/cp*(p00/p)**akapa*(                        
     :         hmfrz(ix,iy,is)+sacrw(ix,iy,is)+iacr(ix,iy,is)+                  
     :         sbercw(ix,iy,is)+ibercw(ix,iy,is)+sacrr(ix,iy,is)-               
     :         smlt(ix,iy,is)-imlt(ix,iy,is))                                   
           else
             transform=0.
           endif
	     pt(ix,iy,is,3)=pt(ix,iy,is,3)
     :       -dtl*(-condens-transform)
	 
!------------water vapor--------------------------!
      if(ifqi.ne.0.) then 
          qv(ix,iy,is,3)=qv(ix,iy,is,3)
     :      -dtl*(-evap(ix,iy,is)
     :      +vini(ix,iy,is))      
	    else
	      qv(ix,iy,is,3)=qv(ix,iy,is,3)
     :      -dtl*(cond(ix,iy,is)-evap(ix,iy,is))
	    endif
!----------could water----------------------------!
         qc(ix,iy,is,3)=qc(ix,iy,is,3)-dqc 
!---------cloud ice-------------------------------!
        if (ifqi.ne.0.) then
          qci(ix,iy,is,3)=qci(ix,iy,is,3)             
     :    -dqci
	  endif
!----------rain------------------------------------!
       if (ifqr.ne.0.) then
	 qr(ix,iy,is,3)=qr(ix,iy,is,3)-dqr
	 endif
!---------snow--------------------------------------!
       if (ifqi.ne.0.) then
        qsn(ix,iy,is,3)=qsn(ix,iy,is,3)-dqsn
	 endif
	    enddo
        enddo
      enddo
      
	return
	end

	subroutine saturation
	use alloc
	implicit none
	real(8) p,t
	integer is,iy,ix
	real, external:: qsat,qsati
      do is=0,ns
        do iy=0,ny1
          do ix=0,nx1
	      p=sigma0(is)*pp(ix,iy,3)+ptop
            t=(pts(ix,iy,is)+pt(ix,iy,is,3))*(p/p00)**akapa
	if(qc(ix,iy,is,3).ne.0.or.qci(ix,iy,is,3).ne.0) then
	      qs(ix,iy,is)=(qc(ix,iy,is,3)*qsat(t,p)+
     :		  qci(ix,iy,is,3)*qsati(t,p))
     :      /(qc(ix,iy,is,3)+qci(ix,iy,is,3))
	 else
	if(t.gt.273.15)	qs(ix,iy,is)=qsat(t,p)
	if(t.le.273.15.and.t.gt.233.15)	
     :qs(ix,iy,is)=0.5*(qsati(t,p)+qsat(t,p))
	if(t.le.233.15)qs(ix,iy,is)=qsati(t,p)
	 endif
	if(t.gt.273.15)	qs(ix,iy,is)=qsat(t,p)
	if(t.le.233.15)qs(ix,iy,is)=qsati(t,p)

	    enddo
        enddo
      enddo
	return
	end

	
c==================================================================================================================c
c         FUNCTIONS
c
c==================================================================================================================c
      function lamsnow(rho,qsnow)
	implicit none
	real(8) rho,qsnow,pi,rosn,lamsnow,xnso
c------------------------------------------------------------------------------------------------------------------c
c        slope parameter for snow distribution, after Lin,1983
c------------------------------------------------------------------------------------------------------------------c
	
	pi=3.1421
	rosn=1000.
	xnso=3.*10**6
	lamsnow=(pi*rosn*xnso/(rho*qsnow))**0.25
	end

	function lamrain(rho,qrain)
c------------------------------------------------------------------------------------------------------------------c
c        slope parameter for rain distribution, after Lin,1983
c------------------------------------------------------------------------------------------------------------------c
	implicit none
	real(8) pi,xnro,rhol,lamrain,rho,qrain
	pi=3.1421
	xnro=8.*10**6
	rhol=1.e3
	lamrain=(pi*rhol*xnro/(rho*qrain))**0.25
	end

	function Usn(rho,qsnow,gamma)
c------------------------------------------------------------------------------------------------------------------c
c     mass-weight mean terminate velocity of precipitating snow, after Lin, 1983
c------------------------------------------------------------------------------------------------------------------c
      implicit none
	real(8) csnow,dsnow,rho0,Usn,qsnow,rho,gamma
	real(8),external:: lamsnow
	csnow=11.72
	dsnow=0.25
      rho0=1.23
	Usn=csnow*gamma*sqrt(rho0/rho)
     :/(6.*lamsnow(rho,qsnow)**dsnow)
	end

	function Urn(rho,qrain,gamma)
c------------------------------------------------------------------------------------------------------------------c
c     mass-weight mean terminate velocity of precipitating rain, after Lin, 1983
c------------------------------------------------------------------------------------------------------------------c
      implicit none
	real(8) arain,brain,rho0,Urn,qrain,rho,gamma
	real(8),external:: gamf,lamrain
	arain=842.
	brain=0.8
      rho0=1.23
	Urn=arain*gamma*sqrt(rho0/rho)
     :/(6.*lamrain(rho,qrain)**brain)
	end

	function gamf(arg)
c-----------------------------------------------------------------------------------------------------------------c
c	! Gamma function calculation                                        DC,10.2009
c------------------------------------------------------------------------------------------------------------------c
	implicit double precision(a-h,o-z)
	real*8 ag(8)
	 data ag/-0.577191652,0.988205891,-0.897056937,0.918206857,
     :    -0.756704078,0.482199394,-0.193527818,0.035868343/

	 do while(arg.gt.2)
	 arg=arg-1.
	 xn=xn+1
	 enddo
	 do i=1,8
	 slagaemoe=ag(i)*(arg-1)**i
	 summa=summa+slagaemoe
	 enddo
	 gamf12=1.+summa
	 xmnojitel=arg+xn
	 do i=1,xn
       xmnojitel=xmnojitel*(arg+xn-i)
	 enddo
	 gamf=xmnojitel*gamf12
	 summa=0.
	 slagaemoe=0.
	 xmnojitel=0.
	 gamf12=0.
	 xn=0.
	end
      
	function a1koenig(t)
	implicit double precision (a-h,o-z)
	t0=273.15
	if(t.gt.t0) a1koenig=0.0
	if(t.le.t0.and.t.gt.t0-1.5) a1koenig=0.7939e-7
	if(t.le.t0-1.5.and.t.gt.t0-2.5) a1koenig=0.7841e-6
	if(t.le.t0-2.5.and.t.gt.t0-3.5) a1koenig=0.3369e-5
	if(t.le.t0-3.5.and.t.gt.t0-4.5) a1koenig=0.4336e-5
	if(t.le.t0-4.5.and.t.gt.t0-5.5) a1koenig=0.5285e-5
	if(t.le.t0-5.5.and.t.gt.t0-6.5) a1koenig=0.3728e-5
	if(t.le.t0-6.5.and.t.gt.t0-7.5) a1koenig=0.1852e-5
	if(t.le.t0-7.5.and.t.gt.t0-8.5) a1koenig=0.2991e-6
	if(t.le.t0-8.5.and.t.gt.t0-9.5) a1koenig=0.4248e-6
	if(t.le.t0-9.5.and.t.gt.t0-10.5) a1koenig=0.7434e-6
	if(t.le.t0-10.5.and.t.gt.t0-11.5) a1koenig=0.1812e-5
	if(t.le.t0-11.5.and.t.gt.t0-12.5) a1koenig=0.4394e-5
	if(t.le.t0-12.5.and.t.gt.t0-13.5) a1koenig=0.9145e-5
	if(t.le.t0-13.5.and.t.gt.t0-14.5) a1koenig=0.1725e-6
	if(t.le.t0-14.5.and.t.gt.t0-15.5) a1koenig=0.3348e-4
	if(t.le.t0-15.5.and.t.gt.t0-16.5) a1koenig=0.1725e-4
	if(t.le.t0-16.5.and.t.gt.t0-17.5) a1koenig=0.9175e-5
	if(t.le.t0-17.5.and.t.gt.t0-18.5) a1koenig=0.4412e-5
	if(t.le.t0-18.5.and.t.gt.t0-19.5) a1koenig=0.2252e-5
	if(t.le.t0-19.5.and.t.gt.t0-20.5) a1koenig=0.9115e-6
	if(t.le.t0-20.5.and.t.gt.t0-21.5) a1koenig=0.4876e-6
	if(t.le.t0-21.5.and.t.gt.t0-22.5) a1koenig=0.3473e-6
	if(t.le.t0-22.5.and.t.gt.t0-23.5) a1koenig=0.4758e-6
	if(t.le.t0-23.5.and.t.gt.t0-24.5) a1koenig=0.6306e-6
	if(t.le.t0-24.5.and.t.gt.t0-25.5) a1koenig=0.8573e-6
	if(t.le.t0-25.5.and.t.gt.t0-26.5) a1koenig=0.7868e-6
	if(t.le.t0-26.5.and.t.gt.t0-27.5) a1koenig=0.7192e-6
	if(t.le.t0-27.5.and.t.gt.t0-28.5) a1koenig=0.6513e-6
	if(t.le.t0-28.5.and.t.gt.t0-29.5) a1koenig=0.5956e-6
	if(t.le.t0-29.5.and.t.gt.t0-30.5) a1koenig=0.5333e-6
	if(t.le.t0-30.5) a1koenig=0.4834e-6	
	end

	function a2koenig(t)
	implicit double precision (a-h,o-z)
	t0=273.15
	if(t.gt.t0) a2koenig=0.0
	if(t.le.t0.and.t.gt.t0-1.5) a2koenig=0.4006
	if(t.le.t0-1.5.and.t.gt.t0-2.5) a2koenig=0.4831
	if(t.le.t0-2.5.and.t.gt.t0-3.5) a2koenig=0.5320
	if(t.le.t0-3.5.and.t.gt.t0-4.5) a2koenig=0.5307
	if(t.le.t0-4.5.and.t.gt.t0-5.5) a2koenig=0.5319
	if(t.le.t0-5.5.and.t.gt.t0-6.5) a2koenig=0.5249
	if(t.le.t0-6.5.and.t.gt.t0-7.5) a2koenig=0.4888
	if(t.le.t0-7.5.and.t.gt.t0-8.5) a2koenig=0.3894
	if(t.le.t0-8.5.and.t.gt.t0-9.5) a2koenig=0.4047
	if(t.le.t0-9.5.and.t.gt.t0-10.5) a2koenig=0.4318
	if(t.le.t0-10.5.and.t.gt.t0-11.5) a2koenig=0.4771
	if(t.le.t0-11.5.and.t.gt.t0-12.5) a2koenig=0.5183
	if(t.le.t0-12.5.and.t.gt.t0-13.5) a2koenig=0.5463
	if(t.le.t0-13.5.and.t.gt.t0-14.5) a2koenig=0.5651
	if(t.le.t0-14.5.and.t.gt.t0-15.5) a2koenig=0.5813
	if(t.le.t0-15.5.and.t.gt.t0-16.5) a2koenig=0.5655
	if(t.le.t0-16.5.and.t.gt.t0-17.5) a2koenig=0.5478
	if(t.le.t0-17.5.and.t.gt.t0-18.5) a2koenig=0.5203
	if(t.le.t0-18.5.and.t.gt.t0-19.5) a2koenig=0.4906
	if(t.le.t0-19.5.and.t.gt.t0-20.5) a2koenig=0.4447
	if(t.le.t0-20.5.and.t.gt.t0-21.5) a2koenig=0.4126
	if(t.le.t0-21.5.and.t.gt.t0-22.5) a2koenig=0.3960
	if(t.le.t0-22.5.and.t.gt.t0-23.5) a2koenig=0.4149
	if(t.le.t0-23.5.and.t.gt.t0-24.5) a2koenig=0.4320
	if(t.le.t0-24.5.and.t.gt.t0-25.5) a2koenig=0.4506
	if(t.le.t0-25.5.and.t.gt.t0-26.5) a2koenig=0.4483
	if(t.le.t0-26.5.and.t.gt.t0-27.5) a2koenig=0.4460
	if(t.le.t0-27.5.and.t.gt.t0-28.5) a2koenig=0.4433
	if(t.le.t0-28.5.and.t.gt.t0-29.5) a2koenig=0.4413
	if(t.le.t0-29.5.and.t.gt.t0-30.5) a2koenig=0.4382
	if(t.le.t0-30.5) a2koenig=0.4361
	end
