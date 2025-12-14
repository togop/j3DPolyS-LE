! GPU-accelerated random number generator using OpenACC
! This module provides GPU-safe random number generation for parallel execution

module randomnumber_gpu_mod
    implicit none
    private
    public :: randomnumber_gpu, init_random_gpu, finalize_random_gpu
    
    ! GPU state for random number generation
    !$acc declare create(random_state)
    integer, allocatable :: random_state(:)
    integer :: random_state_size = 0
    logical :: gpu_initialized = .false.

contains

    subroutine init_random_gpu(seed, nthreads)
        ! Initialize GPU random number generator
        ! seed: random seed
        ! nthreads: number of parallel threads (for state array size)
        implicit none
        integer, intent(in) :: seed
        integer, intent(in) :: nthreads
        integer :: i
        
        random_state_size = nthreads
        allocate(random_state(nthreads))
        
        ! Initialize each thread's random state with different seeds
        do i = 1, nthreads
            random_state(i) = seed + i
        end do
        
        !$acc enter data create(random_state)
        !$acc update device(random_state)
        
        gpu_initialized = .true.
    end subroutine init_random_gpu

    subroutine finalize_random_gpu()
        ! Clean up GPU random number generator
        implicit none
        
        if (gpu_initialized) then
            !$acc exit data delete(random_state)
            deallocate(random_state)
            gpu_initialized = .false.
        end if
    end subroutine finalize_random_gpu

    function randomnumber_gpu(thread_id) result(r)
        ! GPU-safe random number generator
        ! thread_id: thread index (0-based or 1-based depending on usage)
        implicit none
        integer, intent(in) :: thread_id
        real*8 :: r
        integer :: state_idx, state_val
        
        ! Map thread_id to state array index (1-based)
        state_idx = mod(thread_id, random_state_size) + 1
        if (state_idx < 1) state_idx = 1
        if (state_idx > random_state_size) state_idx = random_state_size
        
        ! Simple LCG (Linear Congruential Generator) for GPU
        ! This is a simple implementation - for production, consider cuRAND
        !$acc routine seq
        state_val = random_state(state_idx)
        state_val = mod(16807 * state_val, 2147483647)
        if (state_val == 0) state_val = 1  ! Avoid zero state
        random_state(state_idx) = state_val
        r = real(state_val, kind=8) / 2147483647.0d0
    end function randomnumber_gpu

end module randomnumber_gpu_mod

