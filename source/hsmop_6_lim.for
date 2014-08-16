      subroutine hsmop_6_lim(f,fs,work,nxm,nym,nsm,i0,i1,j0,j1,k0,k1,
     : gama,perdam)
      implicit real*8(a-h,o-z)
      dimension f(0:nxm,0:nym,0:nsm),fs(0:nxm,0:nym,0:nsm)
     :   ,work(0:nxm,0:nym)
      logical perdam
      dimension gama(0:*)
      real(kind(0.d0)),allocatable,dimension(:)::
     :  gama05,um6g,gasco4,umgasc,gama6_05
      real(8) fx1,fx2,fy1,fy2,grx1,grx2,gry1,gry2
c     dimension gama05(0:maxns),um6g(0:maxns),gasco4(0:maxns)
c    :   ,umgasc(0:maxns)
c
      if(.not.allocated(gama05)) then
        allocate(gama05(0:nsm))
        allocate(um6g(0:nsm))
        allocate(gasco4(0:nsm))
        allocate(umgasc(0:nsm))
        allocate(gama6_05(0:nsm))
      endif

      do is=k0,k1
        gama05(is)=0.5*gama(is)
        um6g(is)=1.-6.*gama(is)
        gasco4(is)=2.*gama(is)
        umgasc(is)=1.-8.*gama(is)
        gama6_05(is)=0.25*gama05(is)
      enddo
      
      if(perdam) then
         do 11 is=k0,k1
         do 11 iy=j0,j1
         do 11 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)-fs(ix,iy,is)
11       continue
      endif
      
      do is=k0,k1
      
      do iy=j0+1,j1-1
        do ix=i0+3,i1-3
          fx1=f(ix+3,iy,is)-5.*f(ix+2,iy,is)+
     :     10.*f(ix+1,iy,is)-10.*f(ix,iy,is)+
     :     5.*f(ix-1,iy,is)-f(ix-2,iy,is)
          fx2=f(ix+2,iy,is)-5.*f(ix+1,iy,is)+
     :     10.*f(ix,iy,is)-10.*f(ix-1,iy,is)+
     :     5.*f(ix-2,iy,is)-f(ix-3,iy,is)
          grx1=f(ix+1,iy,is)-f(ix,iy,is)
          grx2=f(ix,iy,is)-f(ix-1,iy,is)
          if(fx1*grx1.le.0) fx1=0.
          if(fx2*grx2.le.0) fx2=0.
          work(ix,iy)=f(ix,iy,is)+gama6_05(is)*(fx1-fx2)
        enddo
        do ix=i0+2,i1-2,i1-i0-4
          fx1=f(ix+2,iy,is)-3.*f(ix+1,iy,is)
     :    +3.*f(ix,iy,is)-f(ix-1,iy,is)
          fx2=f(ix+1,iy,is)-3.*f(ix,iy,is)
     :    +3.*f(ix-1,iy,is)-f(ix-2,iy,is)
          grx1=f(ix+1,iy,is)-f(ix,iy,is)
          grx2=f(ix,iy,is)-f(ix-1,iy,is)
          if(-fx1*grx1.le.0) fx1=0.
          if(-fx2*grx2.le.0) fx2=0.
          work(ix,iy)=f(ix,iy,is)-gama05(is)*(fx1-fx2)
        enddo
        do ix=i0+1,i1-1,i1-i0-2
          work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix-1,iy,is)
     :    +f(ix+1,iy,is)+f(ix,iy,is)+f(ix,iy,is))
        enddo    
      enddo
      
      do iy=j0+1,j1-1
        do ix=i0+1,i1-1
          f(ix,iy,is)=work(ix,iy)
        enddo
      enddo
      
      do ix=i0+1,i1-1
        do iy=j0+3,j1-3
          fy1=f(ix,iy+3,is)-5.*f(ix,iy+2,is)+
     :     10.*f(ix,iy+1,is)-10.*f(ix,iy,is)+
     :     5.*f(ix,iy-1,is)-f(ix,iy-2,is)
          fy2=f(ix,iy+2,is)-5.*f(ix,iy+1,is)+
     :     10.*f(ix,iy,is)-10.*f(ix,iy-1,is)+
     :     5.*f(ix,iy-2,is)-f(ix,iy-3,is)
          gry1=f(ix,iy+1,is)-f(ix,iy,is)
          gry2=f(ix,iy,is)-f(ix,iy-1,is)
          if(fy1*gry1.le.0) fy1=0.
          if(fy2*gry2.le.0) fy2=0.
          work(ix,iy)=f(ix,iy,is)+gama6_05(is)*(fy1-fy2)
        enddo
        do iy=j0+2,j1-2,j1-j0-4
          fy1=f(ix,iy+2,is)-3.*f(ix,iy+1,is)
     :    +3.*f(ix,iy,is)-f(ix,iy-1,is)
          fy2=f(ix,iy+1,is)-3.*f(ix,iy,is)
     :    +3.*f(ix,iy-1,is)-f(ix,iy-2,is)
          gry1=f(ix,iy+1,is)-f(ix,iy,is)
          gry2=f(ix,iy,is)-f(ix,iy-1,is)
          if(-fy1*gry1.le.0) fy1=0.
          if(-fy2*gry2.le.0) fy2=0.
          work(ix,iy)=f(ix,iy,is)-gama05(is)*(fy1-fy2)
        enddo
        do iy=j0+1,j1-1,j1-j0-2
          work(ix,iy)=umgasc(is)*f(ix,iy,is)+gasco4(is)*(f(ix,iy,is)
     :    +f(ix,iy,is)+f(ix,iy-1,is)+f(ix,iy+1,is))
        enddo
      enddo
      
      do iy=j0+1,j1-1
        do ix=i0+1,i1-1
          f(ix,iy,is)=work(ix,iy)
        enddo
      enddo
      
      enddo
      
      if(perdam) then
         do 12 is=k0,k1
         do 12 iy=j0,j1
         do 12 ix=i0,i1
         f(ix,iy,is)=f(ix,iy,is)+fs(ix,iy,is)
12       continue
      endif
      
      return
      end