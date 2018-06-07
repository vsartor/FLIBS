!! math.f90
!! Copyright (C) 2018 Victhor Sartorio
!! * This Source Code Form is part of the 'dynbayes' project.
!! * This Source Code Form is part of the 'flibs' project.
!! * This Source Code Form is subject to the terms of the Mozilla Public
!!   License, v. 2.0. If a copy of the MPL was not distributed with this
!!   file, You can obtain one at http://mozilla.org/MPL/2.0/.
!! * This Source Code Form may use code originally licensed under possibly
!!   different open software licenses. For copying, read through subroutine
!!   and function description comments to check for those instances.

module math
    implicit none

contains

!! Function: log1p
!!
!! Approximation of log(1+x) for small x.
!!
!! Input:
!!   * x: A small real value.
!! Output:
!!   * Approximation of log(1+x)
function log1p(x)
    implicit none

    real*8, intent(in) :: x
    real*8             :: log1p

    real*8             :: z

    z     = 1d0 + x
    log1p = log(z) - ((z - 1d0)-x)/z
end function log1p

!! Function: expm1
!!
!! Approximation of exp(x) - 1 for small x.
!!
!! Input:
!!   * x: A small real value
!! Output:
!!   * Approximation of exp(x) - 1
function expm1(x)
    implicit none

    real*8, intent(in) :: x
    real*8             :: expm1

    real*8             :: a

    a = abs(x)

    if (a .lt. 1d-300) then
        expm1 = x
        goto 600
    end if

    if (a .gt. 0.697d0) then
        expm1 = exp(x) - 1d0
        goto 600
    end if

    if (a .gt. 1d-8) then
        expm1 = exp(x) - 1d0
    else
        expm1 = (0.5d0 * x + 1d0) * x
    end if

    expm1 = expm1 - (1d0 + expm1) * (log1p(expm1) - x)

600 expm1 = expm1
end function expm1

end module math
