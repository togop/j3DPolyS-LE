! fortran version of the ran_uniform.c
! Note: This function is NOT used in GPU code paths
! The GPU code uses repro_rng_mod.f03 for reproducible RNG

function randomnumber() result(r)
    ! Note: Removed !$acc routine seq because random_number() intrinsic
    ! is not supported in NVHPC OpenACC device code
    real*8 :: r    ! same: (KIND=8)
    call random_number(r)
end function randomnumber