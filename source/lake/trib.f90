      SUBROUTINE TRIBTEMP(dt,area_lake,Tw)

!     Subroutine TRIBUTARIES computes the change of the
!     horizontally average lake temperature due to
!     tributaries (inflows and outflows) heat advection

!     INPUT variables:
!     dt --- timestep, sec

!     INTPUT/OUTPUT variables:
!     Tw --- the temperature profile in lake, C

      use driving_params, only : &
     & N_tribin, & 
     & N_tribout, &
     & U_tribin, &
     & U_tribout, &
     & T_tribin, &
     & width_tribin, &
     & width_tribout, &
     & M

      implicit none

      real(8), intent(in)   :: dt
      real(8), intent(in)   :: area_lake
      real(8), intent(inout):: Tw(1:M+1)

      real(8) :: heatinflow
      real(8) :: heatoutflow
      real(8) :: invdt

      integer(4) :: i
      integer(4) :: k

      invdt = 1./dt
     
      do i = 1, M+1

        heatinflow = 0.
        do k = 1, N_tribin
          heatinflow = heatinflow + U_tribin(k,i)*T_tribin(k,i)*width_tribin(k,i)
        enddo
        heatinflow = heatinflow/area_lake

        heatoutflow = 0.
        do k = 1, N_tribout
          heatoutflow = heatoutflow + U_tribout(k,i)*width_tribout(k,i)
        enddo
        heatoutflow = heatoutflow/area_lake

        Tw(i) = (invdt*Tw(i) + heatinflow)/(invdt + heatoutflow)

      enddo

!     print*, 'I am in TRIBUTARIES!' ! Debug string

      RETURN
      END SUBROUTINE TRIBTEMP
