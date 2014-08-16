   subroutine iniperturb
   use alloc
   use refstate
   implicit none
   integer ip,i,is,iy,ix
   real*8 zsia
   integer izia0,izia
   do ip=1,nprof
     call inipro(perdat(0,ip),zperdat(0,ip),ndat,fpro(0),zpro,npro)
   
      do i=0,npro
        dperdz(i,ip)=(fpro(i+1)-fpro(i))/dzpro
        perpro0(i,ip)=fpro(i)-dperdz(i,ip)*zpro(i)
        if (i.lt.100) write(0,*)i,fpro(i)
      enddo
      
      do is=0,ns
      do iy=0,ny1
      do ix=0,nx1
         zsia=phis(ix,iy,is)/g
            
         zsia=phis(ix,iy,is)/g- &
     &      0.5*(phis(ix,iy,ns)+phis(ix,iy,ns-1))/g
            
         izia0=zsia/dzpro
         izia=max(0,min(izia0,npro))

         thpert(ix,iy,is)=perpro0(izia,ip)+dperdz(izia,ip)*zsia
         pt(ix,iy,is,3)=thpert(ix,iy,is)-pts(ix,iy,is)
         if(iy.eq.43.and.ix.eq.24) write(0,*)is,pts(ix,iy,is),pt(ix,iy,is,3)
         
       enddo
       enddo
       inipert(is)=pt(nx,ny,is,3)
       enddo
   enddo
   return
   end