! fortran version of the ran_uniform.c

function randomnumber() result(r)
    real*8 :: r    ! same: (KIND=8)
    call random_number(r)
end function randomnumber