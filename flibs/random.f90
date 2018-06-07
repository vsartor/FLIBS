!! random.f90
!! Copyright (C) 2018 Victhor Sartorio
!! * This Source Code Form is part of the 'dynbayes' project.
!! * This Source Code Form is part of the 'flibs' project.
!! * This Source Code Form is subject to the terms of the Mozilla Public
!!   License, v. 2.0. If a copy of the MPL was not distributed with this
!!   file, You can obtain one at http://mozilla.org/MPL/2.0/.
!! * This Source Code Form may use code originally licensed under possibly
!!   different open software licenses. For copying, read through subroutine
!!   and function description comments to check for those instances.

module random
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
!! file. If, for example, on Windows, do not call this function and look towards
!! `set_seed` instead.
subroutine init_seed()
    implicit none

    open  (rdunit, file = '/dev/urandom', access = 'stream', form = 'unformatted')
    read  (rdunit) state0, state1
    close (rdunit)
end subroutine init_seed

!! Subroutine: set_seed
!!
!! Initializes the state based on a specific input for reproducible results.
!! The transformation of the seed into state is based on transformations using
!! prime numbers and bit rotations followed by jumping three steps.
!!
!! Parameters:
!!   * seed [in]: A 64bit integer from which a state will be determined.
subroutine set_seed(seed)
    implicit none

    integer*8, intent(in) :: seed
    integer*8 :: x

    x = seed
    x = 2860486313_8 * x + 5463458053_8
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

end module random
