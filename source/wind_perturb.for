      subroutine iniperturb
      use alloc, only :
     : hbl,nx1,ny1,ns1,upert,vpert,
     : dx, f_chord, xm_pert,phi,phis,g
      
      implicit none
      integer ix,iy,is
      
      real c1,c2,p,z
      real x(0:nx1),y(0:ny1),amp(0:nx1,0:ny1),
     : wave(0:nx1),wave1(0:nx1),wave2(0:nx1)
      
c     defining parabola parameters

      p=0.5*f_chord
 
c     setting the coordinate axes     
 
      x(0)=-(nx1+1)*dx*0.001/2
      y(0)=-(ny1+1)*dx*0.001/2
      do ix=1,nx1
        x(ix)=x(ix-1)+dx*0.001
        write(0,*) x(ix)
      enddo
      do iy=1,ny1
        y(iy)=y(iy-1)+dx*0.001
      enddo
      
      write(0,*) 'scale=',xm_pert


c     calculating the amplitude modulation
      
      do iy=0,ny1
      do ix=0,nx1
         amp(ix,iy)=(-(x(ix))**2./(4.*p)-(y(iy))**2./(4.*p)+p)/p
     :              *xm_pert
         if (amp(ix,iy).le.0) amp(ix,iy)=0.
         if (amp(ix,iy).gt.0) write(0,*) ix,iy,amp(ix,iy)
      enddo
      enddo
      
c     parabolic wave in east-west direction
      
      do ix=0,nx1
        wave1(ix)=((x(ix)-p)**2./(2.*p)-p/2)/(0.5*p)
        if(wave1(ix).ge.0) wave1(ix)=0.
        wave2(ix)=(-(x(ix)+p)**2./(2.*p)+p/2)/(0.5*p)
        if(wave2(ix).le.0) wave2(ix)=0.
        wave(ix)=wave1(ix)+wave2(ix)
      enddo
      
      do is=0,ns1
      do iy=0,ny1
      do ix=0,nx1
       z=(phi(ix,iy,is)+phis(ix,iy,is))/g
       upert(ix,iy,is)=0.
       vpert(ix,iy,is)=wave(ix)*amp(ix,iy)*exp(-0.001*z)
      enddo
      enddo
      enddo
      
      end
      
      subroutine pert_nudging
      use alloc, only:
     : us,vs,u,v,upert,vpert,tau_nudg,nx1,ny1,ns1
     
      implicit none
      integer ix,iy,is
      
      do is=0,ns1
      do iy=0,ny1
      do ix=0,nx1
      if(upert(ix,iy,is).ne.0)
     : u(ix,iy,is,3)=u(ix,iy,is,3)+tau_nudg*upert(ix,iy,is)
      if(vpert(ix,iy,is).ne.0)
     : v(ix,iy,is,3)=v(ix,iy,is,3)+tau_nudg*vpert(ix,iy,is)
      enddo
      enddo
      enddo
      
      end