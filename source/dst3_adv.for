      subroutine dst3adv (pu,pv,pold,pnew,u,v,phi,phiold,
     : fx,fy,phinew)
      use alloc, only :
     : nx,ny,ns,nx1,ny1,ns1,dtl,dx,dy,dt
      implicit none
      real*8,intent(in):: phi(0:nx1,0:ny1,0:ns1)
      real*8,intent(in):: phiold(0:nx1,0:ny1,0:ns1)
      real*8,intent(in):: pu(0:nx1,0:ny1)
      real*8,intent(in):: pv(0:nx1,0:ny1)
      real*8,intent(in):: pold(0:nx1,0:ny1)
      real*8,intent(in):: pnew(0:nx1,0:ny1)
      real*8,intent(in):: u(0:nx1,0:ny1,0:ns1)
      real*8,intent(in):: v(0:nx1,0:ny1,0:ns1)
      real*8,intent(out):: fx(0:nx1,0:ny1,0:ns1)
      real*8,intent(out):: fy(0:nx1,0:ny1,0:ns1)
      real*8,intent(inout):: phinew(0:nx1,0:ny1,0:ns1)
!      real*8,intent(out):: tendx(0:nx1,0:ny1,0:ns1)
!      real*8,intent(out):: tendy(0:nx1,0:ny1,0:ns1)
      
      real*8 c,d0,d1
      integer ix,iy,is
      
 !     do is = 1,ns-1
 !       do iy = 1,ny
 !         do ix = 1,nx-1
 !         c=u(ix,iy,is)*dtl/dx
   
 !         d0=1./6.*(2.-abs(c))*(1.-abs(c))
 !         d1=1./6.*(1.-abs(c))*(1.+abs(c))
 !         if(u(ix,iy,is).ge.0) then
 !           fx(ix,iy,is)=pu(ix,iy)*u(ix,iy,is)*(phi(ix+1,iy,is)+d0
 !    :       *(phi(ix+1,iy,is)-phi(ix,iy,is))+d1*(phi(ix,iy,is)
 !    :       -phi(ix-1,iy,is)))
 !         else
 !           fx(ix,iy,is)=pu(ix,iy)*u(ix,iy,is)*(phi(ix,iy,is)+d0
 !    :      *(phi(ix+1,iy,is)-phi(ix,iy,is))+d1*(phi(ix+2,iy,is)
 !    :      -phi(ix+1,iy,is)))
 !         endif          
  !        enddo
  !      enddo
  !    enddo
      
  !    do is = 1,ns
  !      do iy = 1,ny
  !           fx(nx,iy,is)=0.5*u(nx,iy,is)*pu(nx,iy)
  !   :       *(phi(nx+1,iy,is)+phi(nx,iy,is))
  !      enddo
  !    enddo
      
  !    call extra3(nx+1,ny+1,ns+1,fx,0,nx+1,0,ny+1,0,ns)
      
      
  !    do is = 1,ns
  !      do iy = 1,ny
  !        do ix = 1,nx
  !           phinew(ix,iy,is)=(phiold(ix,iy,is)*pold(ix,iy)
  !   :        -dt*(fx(ix,iy,is)
  !   :        -fx(ix-1,iy,is))/dx)/pnew(ix,iy)
  !        enddo
  !      enddo
  !    enddo
      
      
      
      do is = 1,ns
        do iy = 1,ny-1
          do ix = 1,nx
          c=v(ix,iy,is)*dtl/dy
          d0=1./6.*(2.-abs(c))*(1.-abs(c))
          d1=1./6.*(1.-abs(c))*(1.+abs(c))
          if(v(ix,iy,is).ge.0) then
           fy(ix,iy,is)=pv(ix,iy)*v(ix,iy,is)*(phi(ix,iy,is)+d0
     :     *(phi(ix,iy+1,is)-phi(ix,iy,is))+d1*(phi(ix,iy,is)
     :     -phi(ix,iy-1,is)))
          else
            fy(ix,iy,is)=pv(ix,iy)*v(ix,iy,is)*(phi(ix,iy+1,is)-d0
     :     *(phi(ix,iy+1,is)-phi(ix,iy,is))-d1*(phi(ix,iy+2,is)
     :      -phi (ix,iy+1,is)))
          endif          
          enddo
        enddo
      enddo
      
      do is = 1,ns
        do ix = 1,nx
             fy(ix,ny,is)=0.5*v(ix,ny,is)*pv(ix,ny)
     :       *(phi(ix,ny+1,is)+phi(ix,iy,is))
        enddo
      enddo
      
      call extra3(nx+1,ny+1,ns+1,fy,0,nx+1,0,ny+1,0,ns)
      
      do is = 1,ns
        do iy = 1,ny
          do ix = 1,nx
             phinew(ix,iy,is)=(phiold(ix,iy,is)*pold(ix,iy)
     :        -dtl*(fy(ix,iy,is)
     :        -fy(ix,iy-1,is))/dy)/pnew(ix,iy)
            
          enddo
        enddo
      enddo
      
      return
      end