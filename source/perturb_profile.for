      subroutine perturb_profile

      use alloc

      implicit real*8 (a-h,o-z)
      
      do iy=0,ny1
        do ix=0,nx1
          do is=0,ns-1
          z=phis(ix,iy,is)/g-
     :      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1))/g
            do iz=1,ndat
              if(zthdat(iz,1).gt.z.and.zthdat(iz-1,1).le.z) then
              
                pt(ix,iy,is,3)=
     :           (zthdat(iz,1)-z)/(zthdat(iz,1)-zthdat(iz-1,1))
     :           *per_dat(iz-1,1)
     :           +(z-zthdat(iz-1,1))/(zthdat(iz,1)-zthdat(iz-1,1))
     :           *per_dat(iz,1)
                !write(0,*)zthdat(iz,1),z,pt(ix,iy,is,3)
              endif
            enddo
          enddo
        enddo
      enddo
      
        do is=0,ns1
        ini_pts(is)=pt(nx/2,ny/2,is,3)
        enddo
      
      do is=0,ns-1
      write(0,*)
     : (phis(nx/2,ny/2,is)-
     : 0.5*(phis(nx/2,ny/2,ns)+phis(nx/2,ny/2,ns-1)))/g,
     :pts(nx/2,ny/2,is)+pt(nx/2,ny/2,is,3) 
      enddo
      !stop
      end