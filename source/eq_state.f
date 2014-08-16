      subroutine eq_state
      use alloc, only: rho,sigma0,pp,ptop,
     :                 r,pts,pt,p00,akapa,qv,qc,qci,
     :                 ns1,ny1,nx1,qif,ifqc,ifqi,
     :                 ns,nx,ny
      implicit none
      real p
      integer is,iy,ix
      
      if (qif.eq.0) then
      do is=0,ns1
         do iy=0,ny1
            do ix=0,nx1
               p=sigma0(is)*pp(ix,iy,2)+ptop
               rho(ix,iy,is)=p/(r*(pts(ix,iy,is)
     :              +pt(ix,iy,is,2))*(p/p00)**akapa)
            enddo
         enddo
      enddo
      elseif(qif.ne.0.and.ifqc.eq.0) then
      do is=0,ns1
         do iy=0,ny1
            do ix=0,nx1
               p=sigma0(is)*pp(ix,iy,2)+ptop
               rho(ix,iy,is)=p/(r*(pts(ix,iy,is)
     :              +pt(ix,iy,is,2))*(p/p00)**akapa*
     :              (1.+0.61*qv(ix,iy,is,2)))
            enddo
         enddo
      enddo
      elseif (ifqc.ne.0.and.ifqi.eq.0) then
      do is=0,ns1
         do iy=0,ny1
            do ix=0,nx1
               p=sigma0(is)*pp(ix,iy,2)+ptop
               rho(ix,iy,is)=p/(r*(pts(ix,iy,is)
     :              +pt(ix,iy,is,2))*(p/p00)**akapa*
     :              (1.+0.61*qv(ix,iy,is,2)-qc(ix,iy,is,2)))
            enddo
         enddo
      enddo
      elseif (ifqi.ne.0) then
      do is=0,ns1
         do iy=0,ny1
            do ix=0,nx1
               p=sigma0(is)*pp(ix,iy,2)+ptop
               rho(ix,iy,is)=p/(r*(pts(ix,iy,is)
     :              +pt(ix,iy,is,2))*(p/p00)**akapa*
     :              (1.+0.61*qv(ix,iy,is,2)-qc(ix,iy,is,2)
     :              -qci(ix,iy,is,2)))
            enddo
         enddo
      enddo
      endif
      end
