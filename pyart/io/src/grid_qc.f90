!  Module: grid_qc.f90


subroutine despeckle(data, window_size, noise_ratio, proc, &
                     fill_value, nx, ny, nz, data_qc)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   integer(kind=4), intent(in)                    :: window_size, proc
   real(kind=8), intent(in)                       :: noise_ratio, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: data
   real(kind=8), intent(out), dimension(nz,ny,nx) :: data_qc


!  Local variables ===========================================================

   real(kind=8)                 :: n_total, n_noise

   logical, dimension(nz,ny,nx) :: m_data

   integer(kind=4)              :: window_x, window_y, window_z, &
                                   i, j, k, imin, imax, jmin, jmax, &
                                   kmin, kmax

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: window_size, proc
   !f2py real(kind=8), intent(in)              :: noise_ratio, fill_value
   !f2py real(kind=8), intent(in)              :: data
   !f2py real(kind=8), intent(out)             :: data_qc

!  ===========================================================================


!  The window size specified by the user should be the total length of the
!  window for each dimension. Therefore we compute half of this length in
!  order to know how far to each side of a grid point the window extends
   window_x = int(0.5d0 * window_size, kind=4)
   window_y = int(0.5d0 * window_size, kind=4)
   window_z = int(0.5d0 * window_size, kind=4)

!  Filter the gridded data using a 3-D box filter. Here we look at all
!  the surrounding grid points within a 3-D window and count how many of
!  these points are noise or missing data. If the number of surrounding grid
!  points that are noise is greater than the user-defined threshold, we flag
!  that point as noise
   data_qc = data
   m_data = data == fill_value

   !$omp parallel num_threads(proc)

   !$omp do
   do i = 1, nx
   
!     Get stencil of x grid points within window
      imin = max(1, i - window_x)
      imax = min(nx, i + window_x)
      
      do j = 1, ny
      
!        Get stencil of y grid points within window
         jmin = max(1, j - window_y)
         jmax = min(ny, j + window_y)
         
         do k = 1, nz
            
!           Get stencil of z grid points within window
            kmin = max(1, k - window_z)
            kmax = min(nz, k + window_z)
            
!           Count the total number of grid points within the 3-D window and
!           the number of grid points that are noise or missing data. Note
!           how the total number of grid points within the 3-D window will
!           be the same value for most grid points, except for grid points
!           that are near the boundaries of the grid, where the window will
!           naturally be smaller
            n_total = (imax-imin+1.d0) * (jmax-jmin+1.d0) * (kmax-kmin+1.d0)
            n_noise = dble(count(m_data(kmin:kmax,jmin:jmax,imin:imax)))
            
            if (100.d0 * n_noise / n_total > noise_ratio) then
               data_qc(k,j,i) = fill_value
            endif
            
         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine despeckle
