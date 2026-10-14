! Triggered component of the ETAS conditional intensity, evaluated at the
! observed events.
!
! Adapted from subroutine `etasfull8newserial` of the R package etasFLP
! (version 2.3.0), Copyright (C) Marcello Chiodi, Giada Adelfio,
! licensed under GPL (>= 2). The numerical computation is unchanged; the
! only modifications are explicit variable declarations (implicit none) and
! the subroutine name.
!
! For i = 2, ..., n:
!   l(i) = k * sum_{j < i, t(j) < t(i)} (t(i) - t(j) + c)^(-p) * exp(predictor(j))
!          * [ (|s_i - s_j|^2 / exp(g * m(j)) + d)^(-q)   if tflag == 0 ]
! l(1) is left untouched (the caller initialises it to zero).
!
! Input magnitudes m must be expressed as m(j) - m0.
! The argument mu is not used; it is kept for compatibility with etasFLP.

subroutine pp_etas_trigger(tflag, n, mu, k, c, p, g, d, q, x, y, t, m, &
                           predictor, l)
  implicit none
  integer, intent(in) :: tflag, n
  double precision, intent(in) :: mu, k, c, p, g, d, q
  double precision, intent(in) :: x(n), y(n), t(n), m(n), predictor(n)
  double precision, intent(inout) :: l(n)

  integer :: i, j
  double precision :: dx, dy, ds, dt, xi, yi, ti, etas, inc, dum

  dum = mu

  do i = 2, n
    inc = 0.0d0
    ti = t(i)
    xi = x(i)
    yi = y(i)
    do j = 1, i - 1
      dt = ti - t(j)
      etas = 0.0d0
      if (dt > 0.0d0) then
        dx = xi - x(j)
        dy = yi - y(j)
        ds = dx * dx + dy * dy
        if (tflag > 0) then
          etas = ((dt + c)**(-p)) * exp(predictor(j))
        else
          etas = ((dt + c)**(-p)) * exp(predictor(j)) * &
                 (ds / exp(g * m(j)) + d)**(-q)
        end if
      end if
      inc = inc + etas
    end do
    l(i) = inc * k
  end do

end subroutine pp_etas_trigger
