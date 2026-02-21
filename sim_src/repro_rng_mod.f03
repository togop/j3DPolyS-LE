! Portable RNG (xorshift64*) for GPU/CPU reproducibility.
! Same seed produces the same sequence on host and device (OpenACC).
module repro_rng_mod
    implicit none
    private
    public :: repro_rng_set_base_seed, repro_rng_clear_base_seed, repro_rng_init_from_trajectory, repro_rng_next

    integer, parameter :: ik8 = selected_int_kind(18)

    ! Base seed for trajectory-based init (set by main when --seed:N is used)
    integer(ik8) :: repro_rng_base_seed = 0_ik8
    logical, public :: repro_rng_use_portable = .false.

contains

    subroutine repro_rng_set_base_seed(seed_value)
        integer, intent(in) :: seed_value
        repro_rng_base_seed = int(seed_value, ik8)
        repro_rng_use_portable = .true.
    end subroutine repro_rng_set_base_seed

    subroutine repro_rng_clear_base_seed()
        repro_rng_base_seed = 0_ik8
        repro_rng_use_portable = .false.
    end subroutine repro_rng_clear_base_seed

    ! Initialize state from (base_seed, trajectory_id). Deterministic. Host-only.
    subroutine repro_rng_init_from_trajectory(trajectory_id, state)
        integer, intent(in) :: trajectory_id
        integer(ik8), intent(out) :: state
        integer(ik8) :: s0, s1
        s0 = repro_rng_base_seed
        s1 = int(trajectory_id, ik8)
        ! Combine into one 64-bit state (must not be 0)
        state = ior(ishft(iand(s0, 4294967295_ik8), 32), iand(s1, 4294967295_ik8))
        if (state == 0_ik8) state = 1_ik8
    end subroutine repro_rng_init_from_trajectory

    ! Advance state and return next value in (0, 1]. xorshift64 (no multiply to avoid overflow)
    function repro_rng_next(state) result(r)
        !$acc routine seq
        integer(ik8), intent(inout) :: state
        real(8) :: r
        integer(ik8) :: x
        x = state
        x = ieor(x, ishft(x, -12))
        x = ieor(x, ishft(x, 25))
        x = ieor(x, ishft(x, -27))
        state = x
        ! Map to (0,1] using lower 32 bits (deterministic, no overflow)
        r = real(iand(x, 4294967295_ik8), 8) / 4294967296.0_8
        if (r <= 0.0_8) r = 1.0e-16_8
    end function repro_rng_next

end module repro_rng_mod
