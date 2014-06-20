!  Module: echo_steiner.f90


subroutine convective_radius(ze_bkg, area_relation, conv_rad)

   implicit none

   real(kind=8), intent(in)      :: ze_bkg
   character(len=16), intent(in) :: area_relation
   real(kind=8), intent(out)     :: conv_rad


!  F2PY directives ===========================================================

   !f2py character(len=16), intent(in) :: area_relation
   !f2py real(kind=8), intent(in)      :: ze_bkg
   !f2py real(kind=8), intent(out)     :: conv_rad

!  ===========================================================================


!  Given a mean background reflectivity value, we determine via a step
!  function what the corresponding convective radius would be
!
!  Higher background reflectivitives are expected to have larger convective
!  influence on surrounding areas, so a larger convective radius would be
!  prescribed

   if (area_relation == 'small') then

      if (ze_bkg < 30.d0) then
         conv_rad = 1000.d0
      elseif (ze_bkg >= 30.d0 .and. ze_bkg < 35.d0) then
         conv_rad = 2000.d0
      elseif (ze_bkg >= 35.d0 .and. ze_bkg < 40.d0) then
         conv_rad = 3000.d0
      elseif (ze_bkg >= 40.d0 .and. ze_bkg < 45.d0) then
         conv_rad = 4000.d0
      else
         conv_rad = 5000.d0
      endif
      

   elseif (area_relation == 'medium') then

      if (ze_bkg < 25.d0) then
         conv_rad = 1000.d0
      elseif (ze_bkg >= 25.d0 .and. ze_bkg < 30.d0) then
         conv_rad = 2000.d0
      elseif (ze_bkg >= 30.d0 .and. ze_bkg < 35.d0) then
         conv_rad = 3000.d0
      elseif (ze_bkg >= 35.d0 .and. ze_bkg < 40.d0) then
         conv_rad = 4000.d0
      else
         conv_rad = 5000.d0
      endif
      
   elseif (area_relation == 'large') then

      if (ze_bkg < 20.d0) then
         conv_rad = 1000.d0
      elseif (ze_bkg >= 20.d0 .and. ze_bkg < 25.d0) then
         conv_rad = 2000.d0
      elseif (ze_bkg >= 25.d0 .and. ze_bkg < 30.d0) then
         conv_rad = 3000.d0
      elseif (ze_bkg >= 30.d0 .and. ze_bkg < 35.d0) then
         conv_rad = 4000.d0
      else
         conv_rad = 5000.d0
      endif
      
   elseif (area_relation == 'sgp') then

      if (ze_bkg < 40.d0) then
         conv_rad = 0.d0
      elseif (ze_bkg >= 40.d0 .and. ze_bkg < 45.d0) then
         conv_rad = 1000.d0
      elseif (ze_bkg >= 45.d0 .and. ze_bkg < 50.d0) then
         conv_rad = 2000.d0
      elseif (ze_bkg >= 50.d0 .and. ze_bkg < 55.d0) then
         conv_rad = 6000.d0
      else
         conv_rad = 8000.d0
      endif
      
   else

      stop

   endif

   return

 end subroutine convective_radius


 subroutine peakedness(ze_bkg, peak_relation, peak)

   implicit none

   real(kind=8), intent(in)      :: ze_bkg
   character(len=16), intent(in) :: peak_relation
   real(kind=8), intent(out)     :: peak


!  F2PY directives ===========================================================

   !f2py character(len=16), intent(in) :: peak_relation
   !f2py real(kind=8), intent(in)     :: ze_bkg
   !f2py real(kind=8), intent(out)    :: peak

!  ===========================================================================


!  Given a background reflectivity value, we determine what the necessary
!  peakedness (or difference) has to be between a grid point's reflectivity
!  and the background reflectivity in order for that grid point to be labeled
!  convective

   if (peak_relation == 'default') then

      if (ze_bkg < 0.d0) then
         peak = 10.d0
      elseif (ze_bkg >= 0.d0 .and. ze_bkg < 42.43d0) then
         peak = 10.d0 - ze_bkg**2 / 180.d0
      else
         peak = 0.d0
      endif
      
   elseif (peak_relation == 'sgp') then

      if (ze_bkg < 0.d0) then
         peak = 14.d0
      elseif (ze_bkg >= 0.d0 .and. ze_bkg < 42.43d0) then
         peak = 14.d0 - ze_bkg**2 / 180.d0
      else
         peak = 4.d0
      endif
      
   else

      stop

   endif

   return

end subroutine peakedness


subroutine classify(ze, x, y, z, dx, dy, intense, bkg_rad, work_level, &
                    area_relation, peak_relation, use_intense, fill_value, &
                    nx, ny, nz, sclass)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   logical, intent(in)                            :: use_intense
   character(len=16), intent(in)                  :: area_relation, &
                                                     peak_relation
   real(kind=8), intent(in)                       :: dx, dy, intense, &
                                                     bkg_rad, work_level, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nx)        :: x
   real(kind=8), intent(in), dimension(ny)        :: y
   real(kind=8), intent(in), dimension(nz)        :: z
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: ze
   integer(kind=4), intent(out), dimension(ny,nx) :: sclass


!  Define local variables  ===================================================

   real(kind=8)                   :: conv_rad, peak, ze_bkg, &
                                     rad, n, sum_ze

   real(kind=8), dimension(ny,nx) :: ze_s

   logical, dimension(ny,nx)      :: m_ze

   integer(kind=4)                :: i, j, k, l, m, &
                                     imin, imax, jmin, jmax, &
                                     lmin, lmax, mmin, mmax

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py logical, intent(in)                   :: use_intense
   !f2py character(len=16), intent(in)         :: area_relation, peak_relation
   !f2py real(kind=8), intent(in)              :: x, y, z, dx, dy
   !f2py real(kind=8), intent(in)              :: intense, bkg_rad
   !f2py real(kind=8), intent(in)              :: ze, work_level, fill_value
   !f2py integer(kind=4), intent(out)          :: sclass

!  ===========================================================================


!  Get index of working level and grab constant height cross section of
!  reflectivity

   k = minloc(abs(z - work_level), dim=1)

   ze_s = ze(k,:,:)

   m_ze = ze_s /= fill_value

!  We perform the Steiner et al. (1995) algorithm for echo classification
!  using only the reflectivity field in order to classify each grid point
!  as either convective, stratiform or undefined. Grid points are classified
!  as follows,
!
!  0 = Undefined
!  1 = Stratiform
!  2 = Convective

   sclass = 0

   do i = 1, nx
   
!     Get stencil of x grid points within the background radius

      imin = max(1, int(i - bkg_rad / dx))
      imax = min(nx, int(i + bkg_rad / dx))
      
      do j = 1, ny
      
         if (m_ze(j,i)) then

!        First make sure that the current grid point has not already been
!        classified. This can happen when grid points within the convective
!        radius of a previous grid point have also been classified

         if (sclass(j,i) == 0) then
         
!           Get stencil of y grid points within the background radius

            jmin = max(1, int(j - bkg_rad / dy))
            jmax = min(ny, int(j + bkg_rad / dy))
            
            n = 0.d0
            sum_ze = 0.d0

!           Calculate the mean background reflectivity for the current grid
!           point, which will be used to determine the convective radius and
!           the required peakedness
            
            do l = imin, imax
               do m = jmin, jmax
               
                  if (m_ze(m,l)) then

                  rad = sqrt((x(l) - x(i))**2 + (y(m) - y(j))**2) ! (m)

!                 The mean background reflectivity will first be computed in
!                 linear units, i.e. mm^6/m^3, then converted to decibel units

                  if (rad <= bkg_rad) then
                     n = n + 1.d0
                     sum_ze = sum_ze + 10.d0**(ze_s(m,l) / 10.d0)
                  endif

                  endif

               enddo
            enddo
            
            ze_bkg = 10.d0 * log10(sum_ze / n) ! (dBZ)

!           Now get the corresponding convective radius knowing the mean
!           background reflectivity

            call convective_radius(ze_bkg, area_relation, conv_rad)
            
!           Now we want to investigate the points surrounding the current
!           grid point that are within the convective radius, and whether
!           they too are convective, stratiform or undefined
!
!           Get stencil of x and y grid points within the convective radius

            lmin = max(1, int(i - conv_rad / dx))
            lmax = min(nx, int(i + conv_rad / dx))
            mmin = max(1, int(j - conv_rad / dy))
            mmax = min(ny, int(j + conv_rad / dy))
            
!           First we check whether the current grid point meets the intensity
!           criteria, which simply means the current grid point has a large
!           enough reflectivity that it must necessarily be convective, as a
!           stratiform grid point could not possibly have this large a
!           reflectivity

            if (use_intense .and. ze_s(j,i) >= intense) then
            
               sclass(j,i) = 2
               
               do l = lmin, lmax
                  do m = mmin, mmax
                  
                     if (m_ze(m,l)) then

                     rad = sqrt((x(l) - x(i))**2 + (y(m) - y(j))**2) ! (m)
                     
                     if (rad <= conv_rad) then
                        sclass(m,l) = 2
                     endif

                     endif

                  enddo
               enddo

            else
            
!              Now get the corresponding peakedness value required for a
!              grid point within the convective radius to be classified as
!              convective

               call peakedness(ze_bkg, peak_relation, peak)
               
               if (ze_s(j,i) - ze_bkg >= peak) then

                  sclass(j,i) = 2
                  
                  do l = imin, imax
                     do m = jmin, jmax
                     
                        if (m_ze(m,l)) then

                        rad = sqrt((x(l) - x(i))**2 + (y(m) - y(j))**2) ! (m)
                        
                        if (rad <= conv_rad) then
                           sclass(m,l) = 2
                        endif

                        endif

                     enddo
                  enddo
                  
               else

!                 If by now the current grid point has not been classified as
!                 convective by either the intensity criteria or the
!                 peakedness criteria, then it must be stratiform

                  sclass(j,i) = 1
                  
               endif

            endif

         endif

         endif

      enddo
   enddo

   return

end subroutine classify
