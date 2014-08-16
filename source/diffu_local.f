	subroutine diffu_local(fngrd
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
     :  ,h1,h2,h3,q1,q2,q3,dif000,dif010,dif001,dif100,dift000
	real(kind(0.d0)),dimension(0:nx1,0:ny1):: lmax
	real*8 z,mixl,def,rich,difk,xnsq,qsat1,qsat2,
     :  qsatur,tp1,tp2,t,p1,p2,p,drag,ztop,tk,dift,
     :  difk_h
	integer is,iy,ix,ixm,iym
	real,parameter:: karm=0.4
	real,external:: qsat
!-----------------------------------------------------------------------------------------------------!
!      write(0,*) 'BROCKEN SOCIAL SCENE 2010'
      do ix=1,nx1
	do iy=1,ny1
      if (ust_s(ix,iy)*tst_s(ix,iy).gt.0) then
!	  lmax(ix,iy)=max(50.,0.007*ust_s(ix,iy)/fcor)
!        lmax(ix,iy)=max(50.,0.15*hbl(ix,iy))
        lmax(ix,iy)=40.
	else
	lmax(ix,iy)= 200. !max(50.,0.15*hbl(ix,iy))
	endif
!      lmax(ix,iy)=150.
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

!	do is=1,ns
!     do iy=1,ny1
!      do ix=1,nx
!        def13(ix,iy,is)=-(s(ix+1,iy,is)+s(ix,iy,is))
!     :    *(u(ix,iy,is,2)-u(ix,iy,is-1,2))/ds02(is)
!      enddo
!      enddo
!      enddo
c
!      do is=1,ns
!      do iy=1,ny
!      do ix=1,nx1
!        def23(ix,iy,is)=-(s(ix,iy,is)+s(ix,iy+1,is))
!     :    *(v(ix,iy,is,2)-v(ix,iy,is-1,2))/ds02(is)
!      enddo
!      enddo
!      enddo

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

	if(iodif.eq.1) then
c-
c calculate linear diffusion coefficients.
c-
        do is=0,ns
        do iy=1,ny1
        do ix=1,nx1
          dif000(ix,iy,is)=difcof
        enddo
        enddo
        enddo
      else
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

	 !
	 ! if(is.eq.35) write(0,*) is,def
	    if(prt.and.driout) bosa(ix,iy,is)=def
c
    !       if(qif.ne.0. .and. ifqc.ne.0) then
              
           if(qif.ne.0.and.ifqc.ge.0) then
             p=sigma0(is)*pp(ix,iy,2)+ptop
              t=(pts(ix,iy,is)+pt(ix,iy,is,2))*(p/p00)**akapa
             if(qv(ix,iy,is,2).ge.0.9*qsat(t,p)) then
            
 !           p2=sigma0(is+1)*pp(ix,iy,2)+ptop
 !           p1=sigma0(is-1)*pp(ix,iy,2)+ptop
            
 !           tp2=(pts(ix,iy,is+1)+pt(ix,iy,is+1,2))*(p2/p00)**akapa
 !           tp1=(pts(ix,iy,is-1)+pt(ix,iy,is-1,2))*(p1/p00)**akapa
 !           qsatur=qsat(t,p)
 !           qsat2=qsat(tp2,p2)
 !           qsat1=qsat(tp1,p1)
 !           xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*
 !    :      ((1+hlat*qsatur/(rv*t))/(1+0.622*hlat*hlatcp*qsatur/
 !    :      (rv*t*t))
 !    :      *((pt(ix,iy,is+1,2)
 !    :      +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
 !    :      /(pt(ix,iy,is,2)+pts(ix,iy,is))+hlatcp/t*(qsat2-qsat1))-
 !    :      (qv(ix,iy,is+1,2)+qc(ix,iy,is+1,2)-qv(ix,iy,is-1,2)-
 !    :      qc(ix,iy,is-1,2)))
 !    :      /ds04a(is)
               
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)+hlat/cp*qv(ix,iy,is+1,2)
     :        -pt(ix,iy,is-1,2)-pts(ix,iy,is-1)
     :        -hlat/cp*qv(ix,iy,is-1,2))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is)+
     :          +hlat/cp*qv(ix,iy,is,2))
     :        /ds04a(is)
               
          else
            xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
           else
             xnsq=-g*(s(ix,iy,is+1)+s(ix,iy,is))*(pt(ix,iy,is+1,2)
     :        +pts(ix,iy,is+1)-pt(ix,iy,is-1,2)-pts(ix,iy,is-1))
     :        /(pt(ix,iy,is,2)+pts(ix,iy,is))
     :        /ds04a(is)
          endif
          
          rich=xnsq/(def*def+zero0)
          if(prt.and.driout) ri(ix,iy,is)=rich
c

          
		  z=(phis(ix,iy,is)+phi(ix,iy,is))/g
          mixl=karm*z/(1.+karm*z/lmax(ix,iy))
	    if(rich.ge.0) then
          difk=mixl**2.*def*(max(0.,(1.-5.*rich))**2.)
	    dift=difk
	    endif
	    if(rich.lt.0) then
	    difk=mixl**2.*def*(min(9.,sqrt(1.-16.*rich)))
	    dift=difk*(min(3.,(1.-16.*rich)**0.25))
	 
	     !difk=delka2(is)*def*dsqrt(dim(1.0d0,rich*rikey))
	     !dift=difk
	
          endif
  !        z=(phi(ix,iy,is)+phis(ix,iy,is))/g
  !         if(z.gt.hbl(ix,iy)) then
  !         difk=0
  !         dift=0
  !         endif
 	   dif000(ix,iy,is)=difk+difl
	    dift000(ix,iy,is)=dift
        enddo
        enddo
        enddo
	call extra3(nx1,ny1,ns1,dif000,1,nx1,1,ny1,0,ns)
	call extra3(nx1,ny1,ns1,dift000,1,nx1,1,ny1,0,ns)
c
        if(prt.and.driout) then
          if(grdout) then
            call wgrids(ri,'ri',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
            call wgrids(dif000,'kk',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	      call wgrids(bosa,'de',2,nx,2,ny,2,ns-1,iomod
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
        difunw(ix,iy,1)=difunw(ix,iy,2)
      enddo
      enddo
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(def13,'ux',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(def23,'vx',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(def33,'wx',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif


     
 !     do is=1,ns-1
 !     do iy=1,ny1
 !     do ix=1,nx-1
 !       umfl(ix,iy,is)=0.25*(def13(ix,iy,is)+def13(ix+1,iy,is)+
 !    :  def13(ix,iy,is+1)+def13(ix+1,iy,is+1))
 !     enddo
 !     enddo
 !     enddo
c
 !     do is=1,ns-1
 !     do iy=1,ny-1
 !     do ix=1,nx1
 !       vmfl(ix,iy,is)=0.25*(def23(ix,iy,is)+def23(ix,iy+1,is)+
 !    :  def23(ix,iy,is+1)+def23(ix,iy+1,is+1))
 !     enddo
 !     enddo
 !     enddo
!--------------------------------------------------------------------------------!
!                            DIFUNT                                              !
!--------------------------------------------------------------------------------!
	if (qif.ne.0.and.ifqc.ge.0) then
	  do is=2,ns-1
        do iy=1,ny
        do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +hlat/cp*qv(ix,iy,is,2)
     :    -hlat/cp*qv(ix,iy,is-1,2)   
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*tdif	
        enddo
        enddo
        enddo
	else
	  do is=2,ns-1
        do iy=1,ny
        do ix=1,nx
        h3(ix,iy,is)=-s(ix,iy,is)*(pt(ix,iy,is,2)-pt(ix,iy,is-1,2)
     :    +tsdif*(pts(ix,iy,is)-pts(ix,iy,is-1)))/ds0(is)
        h3(ix,iy,is)=h3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))*tdif	   
        enddo
        enddo
        enddo
      endif
	
	!write(0,*)pts(18,18,47)+pt(18,18,47,2),u(18,18,47,2),v(18,18,47,2)
c
      do iy=1,ny1
      do ix=1,nx1
        h3(ix,iy,1)=0.0
        !h3(ix,iy,ns)=0.0
	  if(qif.ge.0.and.ifqc.ge.0) then
	    h3(ix,iy,ns)=(ust_s(ix,iy)*tst_s(ix,iy)+
     :    hlat/cp*qst_s(ix,iy)*ust_s(ix,iy))
!	    hfl(ix,iy,ns)=h(ix,iy)+le(ix,iy)
	  else
	    h3(ix,iy,ns)=ust_s(ix,iy)*tst_s(ix,iy)
!	    hfl(ix,iy,ns)=h(ix,iy)
	  endif
	  
      enddo
      enddo

	do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunt(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(h3(ix,iy,is+1)-h3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
      
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(h3,'he',2,nx,2,ny,2,ns-1,iomod
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
     :    )/ds0(is)   !
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :    +dift000(ix,iy,is-1))
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        !q3(ix,iy,ns)=0.0
        q3(ix,iy,ns)=qst_s(ix,iy)*ust_s(ix,iy)
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunq(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo
      endif
      
      if(prt.and.driout) then
          if(grdout) then
            call wgrids(q3,'qf',2,nx,2,ny,2,ns-1,iomod
     :        ,0,0,0,fngrd)
	    endif
	endif
!--------------------------------------------------------------------------------!
!                            DIFUNQC                                              !
!--------------------------------------------------------------------------------!
      if(ifqc.gt.0.) then
      do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qc(ix,iy,is,2)-qc(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :     +dift000(ix,iy,is-1))
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqc(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo

      endif
!--------------------------------------------------------------------------------!
!                            DIFUNQCI                                              !
!--------------------------------------------------------------------------------!		
      if(ifqi.gt.0.) then
	do is=2,ns-1                                                              
      do iy=1,ny                                                                
      do ix=1,nx                                                                
        q3(ix,iy,is)=-s(ix,iy,is)*(qci(ix,iy,is,2)-qci(ix,iy,is-1,2))           
     :     /ds0(is)                                                            
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)                         
     :     +dift000(ix,iy,is-1))*qdif                                          
      enddo                                                                    
      enddo                                                                  
      enddo                                                                
c                                          
      do iy=1,ny1                                                        
      do ix=1,nx1                                                           
        q3(ix,iy,1)=0.0                                                     
        q3(ix,iy,ns)=0.0                                                   
      enddo                                                                  
      enddo                                                                   
c
      do is=1,ns-1                                                         
      do iy=2,ny                                                               
      do ix=2,nx                                                              
        difunqci(ix,iy,is)=                                                 
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))           
     :    /ds12(is)                                                            
      enddo                                                                   
      enddo                                                                   
      enddo                                                                   

      endif         
!--------------------------------------------------------------------------------!
!                            DIFUNQSN                                              !
!--------------------------------------------------------------------------------!

      if(ifqi.gt.0.) then
	do is=2,ns-1                                                           
      do iy=1,ny                                                           
      do ix=1,nx                                                              
        q3(ix,iy,is)=-s(ix,iy,is)*(qsn(ix,iy,is,2)-qsn(ix,iy,is-1,2))         
     :     /ds0(is)                                                            
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)                      
     :     +dift000(ix,iy,is-1))*qdif                                           
      enddo                                                                  
      enddo                                                                 
      enddo                                                                   
c                                          
      do iy=1,ny1                                                             
      do ix=1,nx1                                                             
        q3(ix,iy,1)=0.0                                                       
        q3(ix,iy,ns)=0.0                                                    
      enddo                                                                   
      enddo                                                                    
c
      do is=1,ns-1                                                              
      do iy=2,ny                                                               
      do ix=2,nx                                                              
        difunqsn(ix,iy,is)=                                                    
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))            
     :    /ds12(is)                                                             
      enddo                                                                     
      enddo                                                                  
      enddo                                                                    

      endif 
!--------------------------------------------------------------------------------!
!                            DIFUNQR                                              !
!--------------------------------------------------------------------------------!
	
	if(ifqr.gt.0.) then
	do is=2,ns-1
      do iy=1,ny
      do ix=1,nx
        q3(ix,iy,is)=-s(ix,iy,is)*(qr(ix,iy,is,2)-qr(ix,iy,is-1,2))
     :     /ds0(is)
        q3(ix,iy,is)=q3(ix,iy,is)*0.5*(dift000(ix,iy,is)
     :     +dift000(ix,iy,is-1))*qdif
      enddo
      enddo
      enddo
c
      do iy=1,ny1
      do ix=1,nx1
        q3(ix,iy,1)=0.0
        q3(ix,iy,ns)=0.0
      enddo
      enddo
c
      do is=1,ns-1
      do iy=2,ny
      do ix=2,nx
        difunqr(ix,iy,is)=
     :    -(s(ix,iy,is+1)+s(ix,iy,is))*(q3(ix,iy,is+1)-q3(ix,iy,is))
     :    /ds12(is)
      enddo
      enddo
      enddo	
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
