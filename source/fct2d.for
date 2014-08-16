      subroutine fct2d(f,fm1,fm2,u,um1,v,vm1,wsig,wsigm1,nx,ny,ns
     :   ,dtxs,dtys,dtxy,dxys,dxdt,dydt,dt,ds0
     :   ,flx,fly,fls,sx,sy,ss,wk7,wk8,prt,grdout,iomod,fngrd)
     
      !  FCT advection scheme in horizontal after Zalesak'79
      !  f - the advected field      
      implicit none
      real*8 zero0
      parameter (zero0=1.d-30)
      integer ix,iy,is,nx,ny,ns
      logical prt,grdout
      integer iomod
      character*80 fngrd
      real*8 f(0:nx+1,0:ny+1,0:ns+1),fm1(0:nx+1,0:ny+1,0:ns+1)
     :   ,fm2(0:nx+1,0:ny+1,0:ns+1)
      real*8 u(0:nx+1,0:ny+1,0:ns+1),um1(0:nx+1,0:ny+1,0:ns+1)
      real*8 v(0:nx+1,0:ny+1,0:ns+1),vm1(0:nx+1,0:ny+1,0:ns+1)
      real*8 wsig(0:nx+1,0:ny+1,0:ns+1),wsigm1(0:nx+1,0:ny+1,0:ns+1)
      real*8 flx(0:nx+1,0:ny+1,0:ns+1),fly(0:nx+1,0:ny+1,0:ns+1)
     :   ,fls(0:nx+1,0:ny+1,0:ns+1)
      real*8 sx(0:nx+1,0:ny+1,0:ns+1),sy(0:nx+1,0:ny+1,0:ns+1)
     :   ,ss(0:nx+1,0:ny+1,0:ns+1)
      real*8 wk7(0:nx+1,0:ny+1,0:ns+1),wk8(0:nx+1,0:ny+1,0:ns+1)
      real*8 dtxs(0:ns+1),dtys(0:ns+1),dxys(0:ns+1),ds0(0:ns+1)
      real*8 dtxy,dxdt,dydt,dt
      !local variables
      real*8 pplus(0:nx+1,0:ny+1,0:ns+1),rplus(0:nx+1,0:ny+1,0:ns+1)
     :   ,qplus
      real*8 pminus(0:nx+1,0:ny+1,0:ns+1)
     :   ,rminus(0:nx+1,0:ny+1,0:ns+1)
     :   ,qminus
      real*8 sm(0:nx+1,0:ny+1,0:ns+1)
      real*8 cx(0:nx+1,0:ny+1,0:ns+1)
      real*8 cy(0:nx+1,0:ny+1,0:ns+1)
      real*8 a,b,c
      
c
c      equivalence(wk1,flx),(wk2,fly),(wk3,fls)
c      equivalence(wk4,sx),(wk5,sy),(wk6,ss)
c
      
      !step 1: calculate fluxes at the cells interfaces with higher order scheme, 
      !        not monotonic
      
      do 10 is=1,ns-1
      do 10 iy=2,ny
      do 10 ix=1,nx
      flx(ix,iy,is)=0.5*dtys(is)*um1(ix,iy,is)*(fm1(ix+1,iy,is)
     :   +fm1(ix,iy,is))/dt
10    continue
c
      do 20 is=1,ns-1
      do 20 iy=1,ny
      do 20 ix=2,nx
      fly(ix,iy,is)=0.5*dtxs(is)*vm1(ix,iy,is)*(fm1(ix,iy+1,is)
     :   +fm1(ix,iy,is))/dt
20    continue

      ! step 2: calculate fluxes (sx and sy) with lower order scheme
      ! but monotonic, and make advection step on the field with this 
      ! low order scheme:
      ! step 3: calculate antidiffusive fluxes
      do 60 is=1,ns-1
      do 60 iy=2,ny
      do 60 ix=1,nx
      sx(ix,iy,is)=dtys(is)*(0.5*
     :   ((um1(ix,iy,is)+abs(um1(ix,iy,is)))*fm1(ix,iy,is)
     :   +(um1(ix,iy,is)-abs(um1(ix,iy,is)))*fm1(ix+1,iy,is)))/dt
      flx(ix,iy,is)=flx(ix,iy,is)-sx(ix,iy,is)
60    continue
      
      do 70 is=1,ns-1
      do 70 iy=1,ny
      do 70 ix=2,nx
      sy(ix,iy,is)=dtxs(is)*(0.5*
     :   ((vm1(ix,iy,is)+abs(vm1(ix,iy,is)))*fm1(ix,iy,is)
     :   +(vm1(ix,iy,is)-abs(vm1(ix,iy,is)))*fm1(ix,iy+1,is)))/dt
      fly(ix,iy,is)=fly(ix,iy,is)-sy(ix,iy,is)
70    continue
      
      do 90 is=1,ns-1
      do 90 iy=2,ny
      do 90 ix=2,nx
      f(ix,iy,is)=fm2(ix,iy,is)-2.*dt*(sx(ix,iy,is)-sx(ix-1,iy,is)
     :   +sy(ix,iy,is)-sy(ix,iy-1,is))
     :   /dxys(is)
90    continue

      call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
      call extra3(nx+1,ny+1,ns+1,f,0,nx+1,0,ny+1,0,ns+1)
      
      !step 4: setting some antidiffusive fluxes to zero
      
         do 110 is=1,ns-1
         do 110 iy=2,ny
         do 110 ix=1,nx
         sm(ix,iy,is)=(f(ix+1,iy,is)-f(ix,iy,is))
110       continue
c
c neuman b.c. in the horizontal
c
         do 111 ix=0,nx+1,nx+1
         do 111 is=1,ns-1
         do 111 iy=2,ny
         sm(ix,iy,is)=0.
111       continue
c
         do 120 is=1,ns-1
         do 120 iy=2,ny
         do 120 ix=1,nx
         a=flx(ix,iy,is)*sm(ix,iy,is)
         b=flx(ix,iy,is)*sm(ix+1,iy,is)
         c=flx(ix,iy,is)*sm(ix-1,iy,is)
c         if(a.lt.0. .and. (b.lt.0. .or. c.lt.0.)) flx(ix,iy,is)=0.
         flx(ix,iy,is)=flx(ix,iy,is)+flx(ix,iy,is)*(sign(0.5d0,a)-0.5)
     :      *min(1.d0,-sign(0.5d0,b)-sign(0.5d0,c)+1.d0)
120       continue

           do 210 is=1,ns-1
         do 210 iy=1,ny
         do 210 ix=2,nx
         sm(ix,iy,is)=(f(ix,iy+1,is)-f(ix,iy,is))
210      continue
c
c neuman b.c. in the horizontal
c
         do 211 iy=0,ny+1,ny+1
         do 211 is=1,ns-1
         do 211 ix=2,nx
         sm(ix,iy,is)=0.
211      continue
c
         do 220 is=1,ns-1
         do 220 iy=1,ny
         do 220 ix=2,nx
         a=fly(ix,iy,is)*sm(ix,iy,is)
         b=fly(ix,iy,is)*sm(ix,iy+1,is)
         c=fly(ix,iy,is)*sm(ix,iy-1,is)
c         if(a.lt.0. .and. (b.lt.0. .or. c.lt.0.)) fly(ix,iy,is)=0.
         fly(ix,iy,is)=fly(ix,iy,is)+fly(ix,iy,is)*(sign(0.5d0,a)-0.5)
     :      *min(1.d0,-sign(0.5d0,b)-sign(0.5d0,c)+1.d0)
220      continue

      ! step 5: compute the sum of all antidiffusive flux into grid cell
      call extra3(nx+1,ny+1,ns+1,flx,0,nx+1,1,ny+1,0,ns)
      call extra3(nx+1,ny+1,ns+1,fly,1,nx+1,0,ny+1,0,ns)
      do 410 is=1,ns
      do 410 iy=1,ny+1
      do 410 ix=1,nx+1
      pplus(ix,iy,is)=max(0.d0,flx(ix-1,iy,is))-min(0.d0,flx(ix,iy,is))
     :   +max(0.d0,fly(ix,iy-1,is))-min(0.d0,fly(ix,iy,is))
410   continue
      ! step 6: compute the limitation on net antidif flux into grid cell:
      do 420 is=1,ns
      do 420 iy=1,ny
      do 420 ix=1,nx
      qplus=(max(f(ix,iy,is),f(ix-1,iy,is),f(ix+1,iy,is)
     :   ,f(ix,iy-1,is),f(ix,iy+1,is)
     :   ,fm1(ix,iy,is),fm1(ix-1,iy,is),fm1(ix+1,iy,is)
     :   ,fm1(ix,iy-1,is),fm1(ix,iy+1,is))
     :   -f(ix,iy,is))*dxys(is)/(2.*dt)
         if(pplus(ix,iy,is).gt.0) then
           rplus(ix,iy,is)=min(1.d0,qplus/(pplus(ix,iy,is)+zero0))
         elseif (pplus(ix,iy,is).eq.0) then
           rplus(ix,iy,is)=0.
         endif
420   continue
      
      call extra3(nx+1,ny+1,ns+1,rplus,1,nx+1,1,ny+1,0,ns)
      
      do 440 is=1,ns-1
      do 440 iy=2,ny
      do 440 ix=1,nx
      cx(ix,iy,is)=max(-sign(rplus(ix,iy,is),flx(ix,iy,is))
     :   ,sign(rplus(ix+1,iy,is),flx(ix,iy,is)))
440   continue
c
      do 450 is=1,ns-1
      do 450 iy=1,ny
      do 450 ix=2,nx
      cy(ix,iy,is)=max(-sign(rplus(ix,iy,is),fly(ix,iy,is))
     :   ,sign(rplus(ix,iy+1,is),fly(ix,iy,is)))
450   continue

      do 530 is=1,ns-1
      do 530 iy=2,ny
      do 530 ix=1,nx
      sm(ix,iy,is)=cx(ix,iy,is)*flx(ix,iy,is)
530   continue
c
c calculate flux directing away  from grid point (ix,iy,is) p(ix,iy,is)- .
c
      do 540 is=1,ns-1
      do 540 iy=2,ny
      do 540 ix=2,nx
      pminus(ix,iy,is)=max(0.d0,sm(ix,iy,is))-min(0.d0,sm(ix-1,iy,is))
540   continue
c
c apply flux correction on fly.
c
      do 550 is=1,ns-1
      do 550 iy=1,ny
      do 550 ix=2,nx
      sm(ix,iy,is)=cy(ix,iy,is)*fly(ix,iy,is)
550   continue

      do 560 is=1,ns-1
      do 560 iy=2,ny
      do 560 ix=2,nx
      pminus(ix,iy,is)=pminus(ix,iy,is)+max(0.d0,sm(ix,iy,is))
     :   -min(0.d0,sm(ix,iy-1,is))
560   continue

      do 600 is=1,ns-1
      do 600 iy=2,ny
      do 600 ix=2,nx
      qminus=(f(ix,iy,is)-min(f(ix,iy,is)
     :   ,f(ix-1,iy,is),f(ix+1,iy,is)
     :   ,f(ix,iy-1,is),f(ix,iy+1,is)
     :   ,fm1(ix,iy,is),fm1(ix-1,iy,is),fm1(ix+1,iy,is)
     :   ,fm1(ix,iy-1,is),fm1(ix,iy+1,is)))*dxys(is)/(2.*dt)
        if(pminus(ix,iy,is).gt.0) then
           rminus(ix,iy,is)=min(1.d0,qminus
     :    /(pminus(ix,iy,is)+zero0))
        elseif (pminus(ix,iy,is).eq.0) then
           rminus(ix,iy,is)=0.
        endif
600   continue
c
c determine boundary values of r- .
c
      call extra3(nx+1,ny+1,ns+1,rminus,1,nx+1,1,ny+1,0,ns)
     
      do 620 is=1,ns-1
      do 620 iy=2,ny
      do 620 ix=1,nx
      cx(ix,iy,is)=cx(ix,iy,is)*
     :   max(sign(rminus(ix,iy,is),flx(ix,iy,is)),
     :   -sign(rminus(ix+1,iy,is),flx(ix,iy,is)))
620   continue
c
      do 630 is=1,ns-1
      do 630 iy=1,ny
      do 630 ix=2,nx
      cy(ix,iy,is)=cy(ix,iy,is)*
     :   max(sign(rminus(ix,iy,is),fly(ix,iy,is)),
     :   -sign(rminus(ix,iy+1,is),fly(ix,iy,is)))
630   continue

c final flux corrections
c
      do 700 is=1,ns-1
      do 700 iy=2,ny
      do 700 ix=1,nx
      flx(ix,iy,is)=flx(ix,iy,is)*cx(ix,iy,is)
700   continue
c
      do 710 is=1,ns-1
      do 710 iy=1,ny
      do 710 ix=2,nx
      fly(ix,iy,is)=fly(ix,iy,is)*cy(ix,iy,is)
710   continue
    
     
!      do 430 is=1,ns
!      do 430 iy=1,ny+1
!      do 430 ix=1,nx+1
!      pminus(ix,iy,is)=max(0.d0,flx(ix,iy,is))-min(0.d0,flx(ix-1,iy,is))
!     :   +max(0.d0,fly(ix,iy,is))-min(0.d0,fly(ix,iy-1,is))
!430   continue

!      do 440 is=1,ns
!      do 440 iy=1,ny+1
!      do 440 ix=1,nx+1
!      qminus=(f(ix,iy,is)-min(f(ix,iy,is),f(ix-1,iy,is),f(ix+1,iy,is)
!     :   ,f(ix,iy-1,is),f(ix,iy+1,is)
!     :   ,fm1(ix,iy,is),fm1(ix-1,iy,is),fm1(ix+1,iy,is)
!     :   ,fm1(ix,iy-1,is),fm1(ix,iy+1,is))
!     :   )*dxys(is)/(2.*dt)
!         if(pminus(ix,iy,is).gt.0) then
!             rminus(ix,iy,is)=min(1.d0,qminus/(pminus(ix,iy,is)+zero0))
!         elseif (pminus(ix,iy,is).eq.0) then
!             rminus(ix,iy,is)=0.
!         endif
!440   continue

     
      
      
      ! antidiffusive step - limited antidiffusive fluxes are used
      do 51 is=1,ns-1
      do 51 iy=2,ny
      do 51 ix=2,nx
      f(ix,iy,is)=f(ix,iy,is)-2.*dt*(flx(ix,iy,is)-flx(ix-1,iy,is)
     :   +fly(ix,iy,is)-fly(ix,iy-1,is))
     :   /dxys(is)
51    continue
c
      call extra3(nx+1,ny+1,ns+1,f,1,nx+1,1,ny+1,0,ns)
      
      
!      if(prt) then
!          if(grdout) then
!            call wgrids(flx,'a1',2,nx,2,ny,2,ns-1,iomod
!     :        ,0,0,0,fngrd)
!	      call wgrids(fly,'a2',2,nx,2,ny,2,ns-1,iomod
!     :        ,0,0,0,fngrd)
!	    endif

 !     endif
c
      return

      
      
      end
