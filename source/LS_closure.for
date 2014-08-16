	subroutine diffu_LS(fngrd
     :  ,def11,h1,q1
     :  ,def12,h2,q2
     :  ,def13,h3,q3
     :  ,def22
     :  ,def23
     :  ,def33
     :  ,dif000
     :  ,dif010
     :  ,dif001,dif100)

      use alloc
	implicit none
	character*80 fngrd
      real(kind(0.d0)),dimension(0:nx1,0:ny1,0:ns1)::
     :  def11,def12,def13,def22,def23,def33
     :  ,h1,h2,h3,dif000,dif010,dif001,dif100,dift000,
     :   q1,q2,q3
	real(kind(0.d0)),dimension(0:nx1,0:ny1):: lmax,
     :  wth_h
	real*8 z,mixl,def,rich,difk,xnsq,qsat1,qsat2,
     :  qsatur,tp1,tp2,t,p1,p2,p,drag,ztop,tk,dift,
     :  w_2,w_2s,dift2,difk2,tstar_f,qstar_f,wm3
	integer is,iy,ix,ixm,iym
	real,parameter:: karm=0.4
	real,parameter:: b_hm=3.
	real,parameter:: B2=5.
        real,parameter:: A=4.5
	real,parameter:: b=0.
	real,external:: qsat
        do iy=1,ny
        do ix=1,nx 
	if (h(ix,iy).le.20.) then
           lmax(ix,iy)=40.
	else
	   lmax(ix,iy)=max(50.,0.15*hbl(ix,iy))
        endif
        wm3=wstar(ix,iy)+B2*ust_s(ix,iy)
        wth_h(ix,iy)=-A*wm3/hbl(ix,iy)
	enddo
	enddo
	do ix=1,nx1
	w(ix,ny1,ns+1,2)=w(ix,ny1-1,ns+1,2)
	enddo

	do iy=1,ny1
	w(nx1,iy,ns+1,2)=w(nx1-1,iy,ns+1,2)
	enddo
      
	
      do is=1,ns-1
      do iy=1,ny1
      do ix=1,nx1
        def33(ix,iy,is)=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
	enddo
      enddo
      enddo
c
c boundary conditions d/ds=0
c
      do is=0,ns,ns
      do iy=1,ny1
      do ix=1,nx1
        def33(ix,iy,is)=-(s(ix,iy,is+1)+s(ix,iy,is))
     :    *(w(ix,iy,is+1,2)-w(ix,iy,is,2))/ds12(is)
      enddo
      enddo
      enddo

	do is=1,ns
      do iy=1,ny1
      do ix=1,nx
        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix+1,iy,is,2)-w(ix,iy,is,2))/dx-s1ds4a(is)*ppdx10(ix,iy,2)
     :    *(w(ix+1,iy,is+1,2)+w(ix,iy,is+1,2)-w(ix+1,iy,is-1,2)
     :    -w(ix,iy,is-1,2))
	 enddo
      enddo
      enddo
c
      do is=1,ns
      do iy=1,ny
      do ix=1,nx1
        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
     :    +(w(ix,iy+1,is,2)-w(ix,iy,is,2))/dy-s1ds4a(is)*ppdy01(ix,iy,2)
     :    *(w(ix,iy+1,is+1,2)+w(ix,iy,is+1,2)-w(ix,iy+1,is-1,2)
     :    -w(ix,iy,is-1,2))
      enddo
      enddo
      enddo
      
c-
c calculate deformation/richardson number dependent diffusion coefficient
c-
        if(prt.and.driout) allocate(ri(0:nx1,0:ny1,0:ns1))
	  if(prt.and.driout) allocate(bosa(0:nx1,0:ny1,0:ns1))

	
	
	  do is=1,ns-1
        do iy=2,ny
        do ix=2,nx
c
          def=dsqrt((0.25*(s(ix,iy,is+1)+s(ix,iy,is))*((u(ix,iy,is-1,2)
     :	-u(ix,iy,is+1,2))/ds02(is)+(u(ix-1,iy,is-1,2)
     :    -u(ix-1,iy,is+1,2))
     :    /ds02(is)))**2+(0.25*(s(ix,iy,is+1)+s(ix,iy,is))*
     :    ((v(ix,iy,is-1,2)-v(ix,iy,is+1,2))/ds02(is)+
     :    (v(ix,iy-1,is-1,2)-v(ix,iy-1,is+1,2))/ds02(is)))**2)


	    if(prt.and.driout) bosa(ix,iy,is)=def
	    
!	    if(qif.eq.0) then
		xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :    +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :    /(pt(ix,iy,is,2)+pts(ix,iy,is))/ds04a(is)
 !         else
 !         xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*((pt(ix,iy,is+1,2)
 !    :    +pts(ix,iy,is+1))*(1.+0.61*qv(ix,iy,is+1,2))
 !    :    -(pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
 !    :    *(1.+0.61*qv(ix,iy,is-1,2)))
 !    :    /(pt(ix,iy,is,2)+pts(ix,iy,is))
 !    :    *(1.+0.61*qv(ix,iy,is+1,2))/ds04a(is)
 !         endif

	    rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich

	    z=(phis(ix,iy,is)+phi(ix,iy,is))/g
	    mixl=karm*z/(1+karm*z/lmax(ix,iy))
          
	    if (h(ix,iy).le.20) then !(Fv(ix,iy).le.0) then
	      if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift=difk
	        else
	          difk=mixl**2*def*b
	          dift=difk
	        endif  
	      else
	        difk=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift=difk*(min(3.,(1.-16.*rich)**0.25))
	      endif
	    else

	    if(z.gt.hbl(ix,iy)) then
	
	      if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift=difk
	        else
	          difk=mixl**2*def*b
	          dift=difk
	        endif  
	      else
	        difk=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift=difk*(min(3.,(1.-16.*rich)**0.25))
	      endif
	    
	    else      
	
           if(rich.ge.0) then
	        if(rich.le.0.2) then
                difk2=mixl**2*def*(max(b,(1.-5.*rich))**2)
	          dift2=difk2
	        else
	          difk2=mixl**2*def*b
	          dift2=difk2
	        endif  
	      else
	        difk2=mixl**2*def*(min(9.,sqrt(1.-16.*rich)))
	        dift2=difk2*(min(3.,(1.-16.*rich)**0.25))
	      endif



      w_2s=(1.6*ust_s(ix,iy)**2.*(1-deltaz(ix,iy)/hbl(ix,iy)))**(3./2.)+
     :1.2*wstar(ix,iy)**3.*deltaz(ix,iy)/hbl(ix,iy)*
     :(1-0.9*deltaz(ix,iy)/hbl(ix,iy))**(3./2.)
	w_2s=w_2s**(2./3.)
!      w_2s=1.44*ust_s(ix,iy)**2.*(1.-1.5*dzits(ix,iy))**(2./3.)

	w_2=(1.6*ust_s(ix,iy)**2.*(1-z/hbl(ix,iy)))**(3./2.)+
     :1.2*wstar(ix,iy)**3.*z/hbl(ix,iy)*(1-0.9*z/hbl(ix,iy))**(3./2.)
	w_2=w_2**(2./3.)

	tstar_f=-ust_s(ix,iy)*tst_s(ix,iy)/wstar(ix,iy)
	
	gammah(ix,iy)=b_hm*wstar(ix,iy)**2*tstar_f/w_2s/hbl(ix,iy)
	

	dift=karm*ust_s(ix,iy)*deltaz(ix,iy)/((1-16.*dzits(ix,iy))**(-0.5)
     :-karm*deltaz(ix,iy)/tst_s(ix,iy)*gammah(ix,iy))*
     :((hbl(ix,iy)-z)/(hbl(ix,iy)-deltaz(ix,iy)))**2.*
     :(ust_s(ix,iy)*karm*z+wstar(ix,iy)*hbl(ix,iy)*
     :(z/hbl(ix,iy))**(4./3.))/(ust_s(ix,iy)*karm*deltaz(ix,iy)+
     :wstar(ix,iy)*hbl(ix,iy)*(deltaz(ix,iy)/hbl(ix,iy))**(4./3.))


	difk=dift*((1-16.*dzits(ix,iy))**(-0.25)+b_hm*wstar(ix,iy)*
     :ust_s(ix,iy)*karm*deltaz(ix,iy)/hbl(ix,iy)
     :/((1-16.*dzits(ix,iy))**(-0.25))/w_2s)
       

!	 dift=max(dift,dift2)
!	 difk=max(difk,difk2)

	    endif
	    endif
            
	    dif000(ix,iy,is)=difk+difl
	    dift000(ix,iy,is)=dift
	enddo    
	enddo
	enddo
!	do iy=2,ny
 !       do ix=2,nx
 !       dif000(ix,iy,ns-1)=0.5*(cdm(ix,iy)*deltaz(ix,iy)
 !    :   +dif000(ix,iy,ns-1))
 !       dift000(ix,iy,ns-1)=
 !    :   0.5*(tst_s(ix,iy)*ust_s(ix,iy)*deltaz(ix,iy)/
 !    :         (pts(ix,iy,ns-1)+pt(ix,iy,ns-1,3)-tsurf(ix,iy)+0.0000001)
 !    :   +dift000(ix,iy,ns-1))  
 !     dif000(ix,iy,ns-1)=cdm(ix,iy)*deltaz(ix,iy)
 !     dift000(ix,iy,ns-1)=
 !    :   tst_s(ix,iy)*ust_s(ix,iy)*deltaz(ix,iy)/
 !    :         (pts(ix,iy,ns-1)+pt(ix,iy,ns-1,3)-tsurf(ix,iy)+0.0000001)
    !      write(0,*) ix,iy,dif000(6,18,ns-1),dift000(6,18,ns-1) 
 !       enddo
 !     enddo
      
	call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
	call extra3(nx1,ny1,ns1,dift000,1,nx1,1,ny1,0,ns)
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(dift000,'kt',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	  else
            tk=0.
            call wri3ar(ri,2,nx,2,ny,2,ns-1,'ri      ',tk,iomod,nx,ny)
            call wri3ar(dif000,2,nx,2,ny,2,ns-1,'dif000  ',tk,iomod
     :        ,nx,ny)
          endif
          deallocate(ri)
	    deallocate(bosa)
        endif
       
c
c top boundary enhanced diffusion:
c
      if(olddra) then
        ixm=nx1/2
        iym=ny1/2
        ztop=phis(ixm,iym,1)/g
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          z=(phi(ix,iy,is)+phis(ix,iy,is))/g
          drag=(dragm+(dragt-dragm)*cos(pi05*abs((z-ztop)/(zm-ztop)))
     :      **2)*dim(z,zm)/(z-zm+zero0)
          dif000(ix,iy,is)=dif000(ix,iy,is)+drag
	    dift000(ix,iy,is)=dift000(ix,iy,is)+drag
        enddo
        enddo
        enddo
      elseif(.not.raylei) then
        do is=0,idrmax
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=dif000(ix,iy,is)+endrag(is)
	    dift000(ix,iy,is)=dift000(ix,iy,is)+endrag(is)
        enddo
        enddo
        enddo
      endif

	do is=0,ns
      do iy=1,ny
      do ix=1,nx
        dif010(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix,iy+1,is))
        dif100(ix,iy,is)=0.5*(dif000(ix,iy,is)+dif000(ix+1,iy,is))
      enddo
      enddo
      enddo
c
      call extrah(nx1,ny1,dif010,0,nx1,0,ny1,0,ns)
      call extrah(nx1,ny1,dif100,0,nx1,0,ny1,0,ns)
c
      do is=0,ns
         do iy=1,ny1
         do ix=1,nx1
           def33(ix,iy,is)=def33(ix,iy,is)*dif000(ix,iy,is)
         enddo
         enddo
         enddo
c
c
         do is=1,ns
         do iy=1,ny1
         do ix=1,nx
           def13(ix,iy,is)=def13(ix,iy,is)*0.5*(dif100(ix,iy,is)
     :       +dif100(ix,iy,is-1))
         enddo
         enddo
         enddo
c
         do is=1,ns
         do iy=1,ny
         do ix=1,nx1
           def23(ix,iy,is)=def23(ix,iy,is)*0.5*(dif010(ix,iy,is)
     :       +dif010(ix,iy,is-1))
         enddo
         enddo
         enddo

	do iy=1,ny1
	do ix=1,nx
	def13(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix+1,iy))*u(ix,iy,ns-1,2)
	enddo
	enddo
	do iy=1,ny
         do ix=1,nx1
      def23(ix,iy,ns)=0.5*(cdm(ix,iy)+cdm(ix,iy+1))*v(ix,iy,ns-1,2)
	enddo
	enddo
c--
c difunu
c--
      do is=1,ns-1
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,is)=(
     :    -(s(ix+1,iy,is+1)+s(ix+1,iy,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def13(ix,iy,is+1)-def13(ix,iy,is))/ds14(is))*pp10(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=1,nx
        difunu(ix,iy,0)=difunu(ix,iy,1)
        difunu(ix,iy,ns)=difunu(ix,iy,ns-1)
	!difunu(ix,iy,ns)=0.
      enddo
      enddo

c--
c difunv
c--
      do is=1,ns-1
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,is)=(
     :    -(s(ix,iy+1,is+1)+s(ix,iy+1,is)+s(ix,iy,is+1)+s(ix,iy,is))
     :    *(def23(ix,iy,is+1)-def23(ix,iy,is))/ds14(is))*pp01(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=1,ny
      do ix=2,nx
        difunv(ix,iy,0)=difunv(ix,iy,1)
        difunv(ix,iy,ns)=difunv(ix,iy,ns-1)
	!difunv(ix,iy,ns)=0.
	enddo
      enddo


c--
c  difunw
c--
      do is=2,ns-1
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,is)=(
     :    -s(ix,iy,is)*(def33(ix,iy,is)-def33(ix,iy,is-1))/ds0(is))
     :    *pp(ix,iy,2)
      enddo
      enddo
      enddo
c
      do iy=2,ny
      do ix=2,nx
        difunw(ix,iy,ns)=difunw(ix,iy,ns-1)
	!difunw(ix,iy,ns)=0.
	  difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
      
!----------------------------------------------------------------!
!            DIFUNT                                              !
!----------------------------------------------------------------!
	do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
	  h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*tdif
!	  z=(phi(ix,iy,is)+phis(ix,iy,is))/g
	   z=0.5*(phi(ix,iy,is)+phis(ix,iy,is)+
     :   phi(ix,iy,is+1)+phis(ix,iy,is+1))/g
	  if(z.le.hbl(ix,iy).and.h(ix,iy).gt.20) then
	w_2=(1.6*ust_s(ix,iy)**2.*(1.-z/hbl(ix,iy)))**(3./2.)+
     :1.2*wstar(ix,iy)**3.*z/hbl(ix,iy)*(1.-0.9*z/hbl(ix,iy))**(3./2.)
	w_2=w_2**(2./3.)
	tstar_f=-ust_s(ix,iy)*tst_s(ix,iy)/wstar(ix,iy)
	gammah(ix,iy)=b_hm*wstar(ix,iy)**2*tstar_f/w_2/hbl(ix,iy)
	  h3c(ix,iy,is)=-gammah(ix,iy)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))
!	  h3e(ix,iy,is)=(10.*B2*ust_s(ix,iy))/hbl(ix,iy)*A*
!     :    (z/hbl(ix,iy))**3.
         h3e(ix,iy,is)=-0.2*ust_s(ix,iy)*tst_s(ix,iy)*
     :  (z/hbl(ix,iy))**3.
!         h3e(ix,iy,is)=-wth_h(ix,iy)*(z/hbl(ix,iy))**3.
	else
        h3c(ix,iy,is)=0.
	h3e(ix,iy,is)=0.
	endif
	enddo
      enddo
      enddo
	
c
      do iy=1,ny1
      do ix=1,nx1
        h3(ix,iy,1)=0.
        h3c(ix,iy,1)=0.
	h3e(ix,iy,1)=0.
        !h3(ix,iy,ns)=0.
        h3(ix,iy,ns)=ust_s(ix,iy)*tst_s(ix,iy)
	h3c(ix,iy,ns)=0.
	h3e(ix,iy,ns)=0.
      enddo
      enddo

      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunt(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(h3(ix,iy,is+1)+h3c(ix,iy,is+1)
     :    +h3e(ix,iy,is+1)-h3(ix,iy,is)-h3c(ix,iy,is)-h3e(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
c
	if(prt.and.driout) then
          if(grdout) then
            call wgrids(h3,'hl',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(h3c,'hc',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(h3e,'he',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(h3+h3c+h3e,'ht',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
	
!--------------------------------------------------------------------------------!
!                            DIFUNQ                                              !
!--------------------------------------------------------------------------------!
	if(qif.gt.0.) then

	do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qv(ix,iy,is,2)-qv(ix,iy,is-1,2)
     :    +qsdif*(qvs(ix,iy,is)-qvs(ix,iy,is-1)))/ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*qdif
         z=0.5*(phi(ix,iy,is)+phis(ix,iy,is)+
     :   phi(ix,iy,is-1)+phis(ix,iy,is-1))/g
	  if(z.le.hbl(ix,iy).and.h(ix,iy).gt.20) then
	    w_2=(1.6*ust_s(ix,iy)**2.*(1.-z/hbl(ix,iy)))**(3./2.)+
     :      1.2*wstar(ix,iy)**3.*z/hbl(ix,iy)*
     :      (1.-0.9*z/hbl(ix,iy))**(3./2.)
	    w_2=w_2**(2./3.)
	    qstar_f=-ust_s(ix,iy)*qst_s(ix,iy)/wstar(ix,iy)
	    gammaq(ix,iy)=b_hm*wstar(ix,iy)**2*qstar_f/w_2/hbl(ix,iy)
	    q3c(ix,iy,is)=-gammaq(ix,iy)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))
	  else
          q3c(ix,iy,is)=0.
	  endif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.
        q3c(ix,iy,1)=0.
        q3c(ix,iy,ns)=0.
        !q3(ix,iy,ns)=0.0
        q3(ix,iy,ns)=qst_s(ix,iy)*ust_s(ix,iy)
      enddo
      enddo
      
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunq(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)+q3c(ix,iy,is+1)
     :    -q3(ix,iy,is)-q3c(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
     
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(q3,'qf',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(q3c,'qk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(q3+q3c,'qg',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
      endif
      

	if(iodify.eq.0) then
        call extray(nx1,ny1,difunu,1,nx,2,ny,0,ns)
        call extray(nx1,ny1,difunv,2,nx,2,ny-1,0,ns)
        call extray(nx1,ny1,difunw,1,nx1,2,ny,0,ns)
      endif
c
c these extrapolations get rid of unused boundary values for diffusion:
c
      call extrah(nx1,ny1,difunu,1,nx,1,ny1,0,ns)
      call extrah(nx1,ny1,difunv,1,nx1,1,ny,0,ns)
      call extrah(nx1,ny1,difunw,1,nx1,1,ny1,1,ns) 
      

	return
	end
