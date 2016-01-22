!  ----------------------------------------------------------------------------
!  pyart.retrieve.kdp_brute
!  ========================
!
!  A Fortran module for various routines involving specific differential phase
!  retrievals when more brute force approaches are required, e.g., second-order
!  derivatives for each radar range gate.
!
!  ----------------------------------------------------------------------------

! Necessary and/or potential future improvements to this module:
!
! * High order finite difference schemes still need to be added. Note that the
!   value of using a higher order scheme has yet to be determined and may be
!   fruitless in the end.


subroutine lowpass_maesaka_cost(d2kdr2, Clpf, fill_value, proc, nr, ng, Jlpf)
!  ----------------------------------------------------------------------------
!  Compute the value of the low-pass filter cost functional similar to equation
!  (15) in Maesaka et al. (2012).
!
!  Parameters
!  ----------
!  d2kdr2 : array, dim(r,g), float64
!     Second-order derivative of the control variable k with respect to range.
!     The control variable k is proportional to the square root of specific
!     differential phase.
!  Clpf : float64
!     The low-pass filter (radial smoothness) constraint weight.
!  fill_value : float64
!     Value indicating missing or bad data in second-order derivative input.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
!  Returns
!  -------
!  Jlpf : float64
!     Value of the low-pass filter (radial smoothness) cost functional.
!
!  ----------------------------------------------------------------------------

   implicit none

!  Define inputs and outputs
   integer(kind=4), intent(in)                :: nr, ng
   real(kind=8), intent(in), dimension(nr,ng) :: d2kdr2
   real(kind=8), intent(in)                   :: Clpf, fill_value
   integer(kind=4), intent(in)                :: proc
   real(kind=8), intent(out)                  :: Jlpf

!  Define local variables
   integer(kind=4) :: r, g

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nr, ng
!  f2py real(kind=8), intent(in)              :: d2kdr2, Clpf, fill_value
!  f2py integer(kind=4), intent(in)           :: proc
!  f2py real(kind=4), intent(out)             :: Jlpf

!  The low-pass filter cost computed in this subroutine is given by,
!
!  Jlpf = 0.5 * Clpf * sum[ (d2k/dr2)**2 ] ,
!
!  where the sum is over all range gates for all rays.

   Jlpf = 0.d0

   do g = 1, ng
   do r = 1, nr

      Jlpf = Jlpf + 0.5d0 * Clpf * (d2kdr2(r,g))**2

   enddo
   enddo

   return

end subroutine lowpass_maesaka_cost


subroutine lowpass_maesaka_jac(d2kdr2, dr, Clpf, finite_order, fill_value, &
                               proc, nr, ng, dJlpfdk)
!  ----------------------------------------------------------------------------
!  Compute the Jacobian of the low-pass filter cost functional similar to
!  equation (18) in Maesaka et al. (2012). This subroutine does not currently
!  support radars with variable range resolution.
!
!  Parameters
!  ----------
!  d2kdr2 : array, dim(r,g), float64
!     Second-order derivative of the control variable k with respect to range.
!     The control variable k is proportional to the square root of specific
!     differential phase.
!  dr : float64
!     The range resolution in meters.
!  Clpf : float64
!     The low-pass filter (radial smoothness) constraint weight.
!  finite_order :  'low' or 'high', character
!     The finite difference accuracy used to compute the second-order range
!     derivative of the control variable k.
!  fill_value : float64
!     Value indicating missing or bad data in the derivative input.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
!  Returns
!  -------
!  dJlpfdk : array, dim(nr,ng), float64
!     The Jacobian of the low-pass filter cost functional with respect to the
!     control variable k.
!
!  ----------------------------------------------------------------------------

   implicit none

!  Define inputs and outputs
   integer(kind=4), intent(in)                 :: nr, ng
   real(kind=8), intent(in), dimension(nr,ng)  :: d2kdr2
   real(kind=8), intent(in)                    :: dr, Clpf, fill_value
   character(len=16), intent(in)               :: finite_order
   integer(kind=4), intent(in)                 :: proc
   real(kind=8), intent(out), dimension(nr,ng) :: dJlpfdk

!  Define local variables
   integer(kind=4) :: r, g
   real(kind=8)    :: dr2

!  Define F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nr, ng
!  f2py real(kind=8), intent(in)              :: d2kdr2, dr, Clpf, fill_value
!  f2py character(len=16), intent(in)         :: finite_order
!  f2py integer(kind=4), intent(in)           :: proc
!  f2py real(kind=4), intent(out)             :: dJlpfdk

!  The low-pass filter cost is defined as,
!
!  Jlpf = 0.5 * Clpf * sum[ (d2k/dr2)**2 ] ,
!
!  where the sum is over all range gates for all rays.

!  $omp parallel num_threads(proc)

   dr2 = dr**2

!  The Jacobian of Jlpf when a low finite order has been used to compute the
!  second-order range derivative of the control variable k
   if (finite_order == 'low') then

!     $omp do
      do g = 1, ng
      do r = 1, nr

      if (g > 3 .and. g < ng - 2) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g-1) - 2.d0 * d2kdr2(r,g) + &
                                d2kdr2(r,g+1)) / dr2

      elseif (g == 3) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g-2) + d2kdr2(r,g-1) - 2.d0 * &
                                d2kdr2(r,g) + d2kdr2(r,g+1)) / dr2

      elseif (g == 2) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g+1) - 2.d0 * d2kdr2(r,g) - 2.d0 * &
                                d2kdr2(r,g-1)) / dr2

      elseif (g == 1) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g) + d2kdr2(r,g+1)) / dr2

      elseif (g == ng - 2) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g+2) + d2kdr2(r,g+1) - 2.d0 * &
                                d2kdr2(r,g) + d2kdr2(r,g-1)) / dr2

      elseif (g == ng - 1) then
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g-1) - 2.d0 * d2kdr2(r,g) - 2.d0 * &
                                d2kdr2(r,g+1)) / dr2

      else
         dJlpfdk(r,g) = Clpf * (d2kdr2(r,g) + d2kdr2(r,g-1)) / dr2

      endif

      enddo
      enddo
!     $omp end do

!  The Jacobian of Jlpf when a high finite order has been used to compute the
!  second-order range derivative of the control variable k
   elseif (finite_order == 'high') then

!     $omp do
      do g = 1, ng
      do r = 1, nr

      enddo
      enddo
!     $omp end do

   else

      stop

   endif

!  $omp end parallel

   return

end subroutine lowpass_maesaka_jac


subroutine lowpass_maesaka_term(k, dr, finite_order, fill_value, proc, &
                                nr, ng, d2kdr2)
!  ----------------------------------------------------------------------------
!  Compute the low-pass filter term found in Maesaka et al. (2012). This term
!  represents the second-order derivative of the control variable k with
!  respect to range. This subroutine does not currently support radars with
!  variable range resolution.
!
!  Parameters
!  ----------
!  k : array, dim(nr,ng), float64
!     Control variable k defined in Maesaka et al. (2012). This variable is
!     proportional to the square root of specific differential phase.
!  dr : float64
!     The range resolution in meters.
!  finite_order : 'low' or 'high', character
!     The finite difference accuracy to use when computing the second-order
!     range derivative of the control variable k.
!  fill_value : float64
!     Value indicating missing or bad data.
!  proc : int32
!     Number of parallel threads (CPUs) to use.
!
!  Returns
!  -------
!  d2kdr2 : array, dim(nr,ng), float64
!     Second-order derivative of k with respect to range.
!
!  ----------------------------------------------------------------------------

   implicit none

!  Define local variables
   integer(kind=4), intent(in)                 :: nr, ng
   real(kind=4), intent(in), dimension(nr,ng)  :: k
   real(kind=4), intent(in)                    :: dr, fill_value
   character(len=16), intent(in)               :: finite_order
   integer(kind=4), intent(in)                 :: proc
   real(kind=8), intent(out), dimension(nr,ng) :: d2kdr2

!  Define local variables
   integer(kind=4) :: r, g
   real(kind=8)    :: dr2

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nr, ng
!  f2py real(kind=8), intent(in)              :: k, dr, fill_value
!  f2py character(len=16), intent(in)         :: accuracy
!  f2py integer(kind=4), intent(in)           :: proc
!  f2py real(kind=4), intent(out)             :: d2kdr2

!  $omp parallel num_threads(proc)

   dr2 = dr**2

!  Use a low order finite difference scheme to compute the second-order range
!  derivative
   if (finite_order == 'low') then

!     $omp do
      do g = 1, ng
      do r = 1, nr

!     For interior range gates, i.e.,
!
!     g = [2, ng-1]
!
!     use a centered difference scheme where p = 2. When at ray boundaires,
!     i.e.,
!
!     g = 1 or ng
!
!     use either a forward or backward difference scheme where p = 1 depending
!     on which boundary
!
!     Computing --> d2k/dr2
      if (g > 1 .and. g < ng) then
         d2kdr2(r,g) = (k(r,g+1) - 2.d0 * k(r,g) + k(r,g-1)) / dr2

      elseif (g == 1) then
         d2kdr2(r,g) = (k(r,g) - 2.d0 * k(r,g+1) + k(r,g+2)) / dr2

      else
         d2kdr2(r,g) = (k(r,g) - 2.d0 * k(r,g-1) + k(r,g-2)) / dr2

      endif

      enddo
      enddo
!     $omp end do

!  Use a high order finite difference scheme to compute the second-order range
!  derivative
   elseif (finite_order == 'high') then

      do g = 1, ng
      do r = 1, nr

      enddo
      enddo

   else

      stop

   endif

!  $omp end parallel

   return

end subroutine lowpass_maesaka_term


subroutine fill_psidp(psidp, fill_value, proc, nr, ng)
!  ----------------------------------------------------------------------------
!  Fill in missing range gates with accumulated differential phase data, e.g.,
!  if a range gate is missing, fill it in with the differential phase of the
!  previous range gate. If the entire ray is missing, all range gates are
!  filled in with 0 deg.
!
!  Parameters
!  ----------
!  psidp : array, dim(r,g), float64
!     Differential phase data in degrees.
!  fill_value : float64
!     Value indicating missing or bad data in the differential phase array.
!  proc : int32
!     The number of parallel threads (CPUs) to use.
!
!  ----------------------------------------------------------------------------

   implicit none

!  Define inputs and outputs
   integer(kind=4), intent(in)                   :: nr, ng
   real(kind=8), intent(in)                      :: fill_value
   integer(kind=4), intent(in)                   :: proc
   real(kind=8), intent(inout), dimension(nr,ng) :: psidp

!  Define local variables
   real(kind=8), parameter   :: atol=1.d-5
   integer(kind=4)           :: r, g, g0
   logical, dimension(nr,ng) :: mask

!  F2PY directives
!  f2py integer(kind=4), optional, intent(in) :: nr, ng
!  f2py real(kind=8), intent(in)              :: fill_value
!  f2py real(kind=8), intent(inout)           :: psidp

!  $omp parallel num_threads(proc)

!  Determine all radar gates which have missing or bad values
   mask = abs(psidp - fill_value) <= atol

!  $omp do
   do r = 1, nr
   do g = 1, ng

      if (g > 1 .and. mask(r,g)) then
         psidp(r,g) = psidp(r,g-1)

!     Special case where no differential phase measurements are available in
!     the entire ray, in which case set the differential phase to zero degrees
!     at all range gates
      elseif (all(mask(r,:))) then
         psidp(r,:) = 0.d0
         exit

!     Special case where the first radar gate is missing a differential phase
!     measurement but there are differential phase measurements at further
!     range gates along the ray
      elseif (g == 1 .and. mask(r,g)) then
         do g0 = 2, ng
            if (.not. mask(r,g0)) then
               psidp(r,g) = psidp(r,g0)
            endif
         enddo

      endif

   enddo
   enddo
!  $omp end do

!  $omp end parallel

   return

end subroutine fill_psidp
