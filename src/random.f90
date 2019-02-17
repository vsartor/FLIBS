!! random.f90
!! Copyright (C) 2018 Victhor Sartorio
!! * This Source Code Form is part of the FLIBS project.
!! * This Source Code Form is subject to the terms of the Mozilla Public
!!   License, v. 2.0. If a copy of the MPL was not distributed with this
!!   file, You can obtain one at http://mozilla.org/MPL/2.0/.
!! * This Source Code Form may use code originally licensed under possibly
!!   different open software licenses. For copying, read through subroutine
!!   and function description comments to check for those instances.

module random
    use math
    implicit none

    ! Unit ID to be used for opening the `urandom` device
    integer, parameter, private :: rdunit = 155

    ! Internal state of the Xoroshiro128+, 2 integers of 64 bits
    integer*8, save, private :: state0 = 0
    integer*8, save, private :: state1 = 0

contains


!! Subroutine: init_seed
!!
!! Initializes the state of the random number generator with random values, read
!! from the system's random device. Relies on '/dev/urandom' being a readable
!! file. If, for example, on Windows, do not call this function and use `set_seed`
!! instead.
subroutine init_seed()
    implicit none

    open  (rdunit, file='/dev/urandom', access='stream', form='unformatted')
    read  (rdunit) state0, state1
    close (rdunit)
end subroutine init_seed


!! Subroutine: set_seed
!!
!! Initializes the state based on a specific input for reproducible results.
!! The transformation of the seed into state is based on transformations using
!! prime numbers and bit rotations, followed by jumping three steps into the
!! process.
!!
!! Parameters:
!!   * seed [in]: A 64bit integer from which a state will be determined.
subroutine set_seed(seed)
    implicit none

    integer*8, intent(in) :: seed
    integer*8 :: x

    x = 2860486313_8 * seed + 5463458053_8
    x = ishftc(x, 55)
    state0 = x
    x = 3267000013_8 * x + 9576890767_8
    x = ishftc(x, 36)
    state1 = x

    x = rand64()
    x = rand64()
    x = rand64()
end subroutine set_seed


!! Function: rand64
!!
!! Returns a pseudo-random sequence of 64 bits as an integer.
!!
!! This function uses the xoroshiro128+ algorithm, described in
!! http://xoroshiro.di.unimi.it/, authored by David Blackman and Sebastiano
!! Vigna in 2016. The implementation is a translation from the C code available
!! at http://xoroshiro.di.unimi.it/xoroshiro128plus.c.
!!
!! Input:
!!   * None
!! Output:
!!   * A 64-bit pseudo-random integer.
function rand64()
    implicit none

    integer*8 :: rand64, s0, s1

    ! The random value to be returned
    rand64 = state0 + state1

    ! Advancing the state one step
    s0 = state0
    s1 = ieor(state1, s0)
    state0 = ieor(ieor(ishftc(s0, 55), s1), ishft(s1, 14))
    state1 = ishftc(S1, 36)
end function rand64


!! Function: randf
!!
!! Returns a random value between 0 and 1.
!!
!! Input:
!!   * None
!! Output:
!!   * A random value between 0 and 1.
function randf()
    implicit none

    real*8 :: randf

    randf = (real(rand64(), 8) / 18446744073709551615d0) + 0.5d0
end function randf


!! Function: rnorm1
!!
!! Returns single normally distributed random value.
!!
!! Input:
!!   * mu:    The mean of the distribution.
!!   * sigma: The scale of the distribution.
!! Output:
!!   * A random value.
function rnorm1(mu, sigma)
    implicit none

    real*8,    intent(in) :: mu, sigma
    real*8                :: rnorm1

    real*8                :: u, v, s

    s = 1.1d0
    do while (s .ge. 1)
        u = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
        v = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
        s = u*u + v*v
    end do
    rnorm1 = u * sqrt(-2 * log(s) / s) * sigma + mu
end function rnorm1


!! Function: rnorm
!!
!! Returns a vector of normally distributed random values.
!!
!! Input:
!!   * n:     The number of values to be generated.
!!   * mu:    The mean of the distribution.
!!   * sigma: The scale of the distribution.
!! Output:
!!   * A vector with random values.
function rnorm(n, mu, sigma)
    implicit none

    integer*8, intent(in) :: n
    real*8,    intent(in) :: mu, sigma
    real*8                :: rnorm(n)

    real*8                :: u, v, s, nfac
    integer*8             :: i

    if (mod(n, 2) .eq. 1) then
        s = 1.1d0
        do while (s .ge. 1d0)
            u = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
            v = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
            s = u*u + v*v
        end do
        rnorm(n) = u * sqrt(-2 * log(s) / s) * sigma + mu
    end if

    do i = 2, n, 2
        s = 1.1d0
        do while (s .ge. 1d0)
            u = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
            v = (real(rand64(), 8) / 18446744073709551615d0) * 2d0
            s = u*u + v*v
        end do
        nfac       = sqrt(-2 * log(s) / s)
        rnorm(i-1) = u * nfac * sigma + mu
        rnorm(i)   = v * nfac * sigma + mu
    end do
end function rnorm


!! Subourtine: rmvnorm
!!
!! Returns a matrix with row vectors obtained from a
!! multivariate normal distribution.
!!
!! Parameters:
!!   * n [in]:     The number of samples desired.
!!   * m [--]:     The dimension of mu and order of sigma.
!!   * mu [in]:    The mean vector of the distribution.
!!   * sigma [in]: The covariance matrix of the distribution.
!!   * y [out]:    The matrix with the random values.
subroutine rmvnorm(n, m, mu, sigma, y)
    implicit none

    integer,   intent(in)  :: m
    integer*8, intent(in)  :: n
    real*8,    intent(in)  :: mu(m), sigma(m,m)
    real*8,    intent(out) :: y(n,m)

    integer                :: errcode
    integer*8              :: i
    real*8                 :: z(m)

    ! Cholesky decomposition of the covariance matrix into lower
    ! triangular form, stored overriding the sigma matrix.
    call dpotrf2('L', m, sigma, m, errcode)
    if (errcode .ne. 0) then
        error stop 'Error code from dpotrf2 on rmvnorm.'
    end if

    do i = 1, n
        z = rnorm(int(m, 8), 0d0, 1d0)
        ! Multiply lower triangular matrix to vector. The result
        ! overwrites the vector z.
        call dtrmv('L', 'N', 'N', m, sigma, m, z, 1)
        y(i,:) = z + mu
    end do
end subroutine rmvnorm


!! Function: rexp1
!!
!! Generates a single value from the standard exponential
!! distribution.
!!
!! Input:
!!   * None.
!! Output:
!!   * A random positive value.
function rexp1()
    implicit none

    real*8 :: rexp1

    rexp1 = -log(randf())
end function rexp1


!! Function: rgamma1
!!
!! Returns a value from a gamma distribution.
!! Based on implementation from Rmath project (GPL v2).
!!
!! Input:
!!   * a:    The shape parameter
!!   * rate: The rate parameter
!! Output:
!!   * A positive random value.
function rgamma1(a, rate) result(ret_val)
    implicit none

    real*8, parameter  :: sqrt32 =  5.656854d0
    real*8, parameter  :: exp_m1 =  0.36787944117144232159d0
    real*8, parameter  :: q1     =  0.04166669d0
    real*8, parameter  :: q2     =  0.02083148d0
    real*8, parameter  :: q3     =  0.00801191d0
    real*8, parameter  :: q4     =  0.00144121d0
    real*8, parameter  :: q5     = -7.3880d-5
    real*8, parameter  :: q6     =  2.4511d-4
    real*8, parameter  :: q7     =  2.4240d-4
    real*8, parameter  :: a1     =  0.3333333d0
    real*8, parameter  :: a2     = -0.2500030d0
    real*8, parameter  :: a3     =  0.2000062d0
    real*8, parameter  :: a4     = -0.1662921d0
    real*8, parameter  :: a5     =  0.1423657d0
    real*8, parameter  :: a6     = -0.1367177d0
    real*8, parameter  :: a7     =  0.1233795d0

    real*8, intent(in) :: a, rate
    real*8             :: ret_val

    real*8             :: scale, s, s2, d, q0, b, si, c
    real*8             :: e, p, q, r, t, u, v, w, x

    if (a .le. 0 .or. rate .le. 0) then
        ret_val = 0d0
        goto 600
    end if

    scale = 1d0 / rate

    ! -- GS algorithm for parameters a < 1 --
    if (a .lt. 1) then
        e = 1d0 + exp_m1 * a
        do
            p = e * randf()
            if (p .ge. 1d0) then
                x = -log((e - p) / a)
                if (rexp1() .ge. (1d0 - a) * log(x)) then
                    exit
                end if
            else
                x = exp(log(p) / a)
                if (rexp1() .ge. x) then
                    exit
                end if
            end if
        end do
        ret_val = scale * x
        goto 600
    end if

    ! -- GD algorithm for parameters a >= 1 --

    s2 = a - 0.5d0
    s  = sqrt(s2)
    d  = sqrt32 - s * 12d0

    t       = rnorm1(0d0, 1d0)
    x       = s + 0.5d0 * t
    ret_val = x * x
    if (t .ge. 0) then
        ret_val = scale * ret_val
        goto 600
    end if

    u = randf()
    if (d * u .le. t * t * t) then
        ret_val = scale * ret_val
        goto 600
    end if

    r  = 1d0 / a
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
    if (a .le. 3.686d0) then
        b  = 0.463d0 + s + 0.178d0 * s2
        si = 1.235d0
        c  = 0.195d0 / s - 0.079d0 + 0.16d0 * s
    else if (a .le. 13.022) then
        b  = 1.654d0 + 0.0076d0 * s2
        si = 1.68d0 / s + 0.275d0
        c  = 0.062d0 / s + 0.024d0
    else
        b  = 1.77d0
        si = 0.75d0
        c  = 0.1515d0 / s
    end if

    if (x .gt. 0) then
        v = t / (s + s)
        if (abs(v) .le. 0.25d0) then
            q = q0+0.5d0*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
        else
            q = q0-s*t+0.25d0*t*t+(s2+s2)*log(1d0+v)
        end if

        if (log(1d0-u) .le. q) then
            ret_val = scale * ret_val
            goto 600
        end if
    end if

    do
        e = rexp1()
        u = randf()
        u = u + u - 1d0
        if (u .lt. 0) then
            t = b - si * e
        else
            t = b + si * e
        end if

        if (t .ge. -0.71874483771719d0) then
            v = t / (s + s)
            if (abs(v) .le. 0.25d0) then
                q = q0+0.5d0*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
            else
                q = q0-s*t+0.25d0*t*t+(s2+s2)*log(1d0 + v)
            end if
            if (q .gt. 0) then
                w = expm1(q)
                if (c * abs(u) .le. w * exp(e - 0.5d0 * t * t)) then
                    exit
                end if
            end if
        end if
    end do

    x = s + 0.5d0 * t
    ret_val = scale * x * x
    goto 600

600 continue
end function rgamma1


!! Function: rgamma
!!
!! Returns a vector with values from a gamma distribution.
!! Based on implementation from Rmath project (GPL v2).
!!
!! Input:
!!   * n:    The sample size required.
!!   * a:    The shape parameter
!!   * rate: The rate parameter
!! Output:
!!   * The vector with positive random values.
function rgamma(n, a, rate) result(ret_val)
    implicit none

    real*8, parameter     :: sqrt32 =  5.656854d0
    real*8, parameter     :: exp_m1 =  0.36787944117144232159d0
    real*8, parameter     :: q1     =  0.04166669d0
    real*8, parameter     :: q2     =  0.02083148d0
    real*8, parameter     :: q3     =  0.00801191d0
    real*8, parameter     :: q4     =  0.00144121d0
    real*8, parameter     :: q5     = -7.3880d-5
    real*8, parameter     :: q6     =  2.4511d-4
    real*8, parameter     :: q7     =  2.4240d-4
    real*8, parameter     :: a1     =  0.3333333d0
    real*8, parameter     :: a2     = -0.2500030d0
    real*8, parameter     :: a3     =  0.2000062d0
    real*8, parameter     :: a4     = -0.1662921d0
    real*8, parameter     :: a5     =  0.1423657d0
    real*8, parameter     :: a6     = -0.1367177d0
    real*8, parameter     :: a7     =  0.1233795d0

    integer*8, intent(in) :: n
    real*8,    intent(in) :: a, rate
    real*8                :: ret_val(n)

    integer*8             :: i
    real*8                :: scale, s, s2, d, q0, b, si, c
    real*8                :: e, p, q, r, t, u, v, w, x

    if (a .le. 0 .or. rate .le. 0) then
        ret_val(1:n) = 0d0
        goto 600
    end if

    scale = 1d0 / rate

    ! -- GS algorithm for parameters a < 1 --
    if (a .lt. 1) then
        do i = 1, n
            e = 1d0 + exp_m1 * a
            do
                p = e * randf()
                if (p .ge. 1d0) then
                    x = -log((e - p) / a)
                    if (rexp1() .ge. (1d0 - a) * log(x)) then
                        exit
                    end if
                else
                    x = exp(log(p) / a)
                    if (rexp1() .ge. x) then
                        exit
                    end if
                end if
            end do
            ret_val(i) = scale * x
        end do
        goto 600
    end if

    ! -- GD algorithm for parameters a >= 1 --

    ! Pre-compute constant values (Block 1 and 2)

    ! Block 1
    s2 = a - 0.5d0
    s  = sqrt(s2)
    d  = sqrt32 - s * 12d0

    ! Block 2
    r  = 1d0 / a
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
    if (a .le. 3.686d0) then
        b  = 0.463d0 + s + 0.178d0 * s2
        si = 1.235d0
        c  = 0.195d0 / s - 0.079d0 + 0.16d0 * s
    else if (a .le. 13.022) then
        b  = 1.654d0 + 0.0076d0 * s2
        si = 1.68d0 / s + 0.275d0
        c  = 0.062d0 / s + 0.024d0
    else
        b  = 1.77d0
        si = 0.75d0
        c  = 0.1515d0 / s
    end if

    ! Actual generation loop

    outer: do i = 1, n
        t          = rnorm1(0d0, 1d0)
        x          = s + 0.5d0 * t
        ret_val(i) = x * x
        if (t .ge. 0) then
            ret_val(i) = scale * ret_val(i)
            cycle outer
        end if

        u = randf()
        if (d * u .le. t * t * t) then
            ret_val(i) = scale * ret_val(i)
            cycle outer
        end if

        if (x .gt. 0) then
            v = t / (s + s)
            if (abs(v) .le. 0.25d0) then
                q = q0+0.5d0*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
            else
                q = q0-s*t+0.25d0*t*t+(s2+s2)*log(1d0+v)
            end if

            if (log(1d0-u) .le. q) then
                ret_val(i) = scale * ret_val(i)
                cycle outer
            end if
        end if

        inner: do
            e = rexp1()
            u = randf()
            u = u + u - 1d0
            if (u .lt. 0) then
                t = b - si * e
            else
                t = b + si * e
            end if

            if (t .ge. -0.71874483771719d0) then
                v = t / (s + s)
                if (abs(v) .le. 0.25d0) then
                    q = q0+0.5d0*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
                else
                    q = q0-s*t+0.25d0*t*t+(s2+s2)*log(1d0 + v)
                end if
                if (q .gt. 0) then
                    w = expm1(q)
                    if (c * abs(u) .le. w * exp(e - 0.5d0 * t * t)) then
                        exit inner
                    end if
                end if
            end if
        end do inner

        x = s + 0.5d0 * t
        ret_val(i) = scale * x * x
    end do outer

    ! Quick exit statement
600 continue
end function rgamma

end module random
