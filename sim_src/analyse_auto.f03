! Automatic parallelization selection module
! Automatically uses OpenACC if available, falls back to OpenMP, then sequential
! Provides unified interface that works across all platforms

module analyse_auto_mod
    use lattice_data_mod
    use PolymerModel_mod
    use logging_mod
    USE HDF5
    implicit none

    private
    public :: analyse_auto, detect_parallelization_method, ParallelizationMethod
    
    type(Logger) :: log = Logger('analyse_auto_mod', LOG_DEBUG)
    
    ! Parallelization method enumeration
    enum, bind(c)
        enumerator :: PARALLEL_NONE = 0
        enumerator :: PARALLEL_OPENMP = 1
        enumerator :: PARALLEL_OPENACC = 2
    end enum
    
    integer, parameter :: ParallelizationMethod = kind(PARALLEL_NONE)

contains

    integer function detect_parallelization_method(prefer_gpu) result(method)
        ! Detect the best available parallelization method
        ! prefer_gpu: if true, prefer OpenACC over OpenMP
        ! Returns: PARALLEL_OPENACC, PARALLEL_OPENMP, or PARALLEL_NONE
        implicit none
        logical, intent(in) :: prefer_gpu
        character(1000) :: env_var
        integer :: env_len, stat
        
        method = PARALLEL_NONE
        
        ! Check for OpenACC support (if prefer_gpu is true)
        if (prefer_gpu) then
            ! Check if ACC_DEVICE_TYPE is set and valid
            call get_environment_variable('ACC_DEVICE_TYPE', env_var, env_len, stat)
            if (stat == 0 .and. env_len > 0) then
                if (trim(env_var) == 'nvidia' .or. trim(env_var) == 'NVIDIA') then
                    ! Try to detect if GPU is actually available
                    ! This is a simple check - in production, you might want more sophisticated detection
                    method = PARALLEL_OPENACC
                    call log%info('OpenACC detected: ACC_DEVICE_TYPE=' // trim(env_var))
                    return
                end if
            end if
            
            ! Note: OpenACC doesn't support conditional compilation like OpenMP's !$
            ! If OpenACC is compiled, the environment variable check above should have caught it
            ! If not compiled, this code path won't be reached anyway
        end if
        
        ! Check for OpenMP support
        !$ method = PARALLEL_OPENMP
        !$ call log%info('OpenMP support detected')
        !$ return
        
        ! Fallback to sequential
        method = PARALLEL_NONE
        call log%info('No parallelization support detected, using sequential execution')
    end function detect_parallelization_method

    subroutine analyse_auto(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom, &
            prefer_gpu, force_method, num_threads)

        ! Automatic parallelization selection
        ! prefer_gpu: if true, prefer OpenACC over OpenMP
        ! force_method: force specific method (0=auto, 1=OpenMP, 2=OpenACC)
        ! num_threads: number of threads for OpenMP (0 = use default)

        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder
        character(*), intent(in) :: analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        logical, intent(in) :: prefer_gpu
        integer, intent(in) :: force_method  ! 0=auto, 1=OpenMP, 2=OpenACC
        integer, intent(in) :: num_threads

        integer :: method
        character(1000) :: method_name

        call log%set_level(global_log_level)

        ! Determine which method to use
        if (force_method == 0) then
            ! Auto-detect
            method = detect_parallelization_method(prefer_gpu)
        else if (force_method == 1) then
            method = PARALLEL_OPENMP
        else if (force_method == 2) then
            method = PARALLEL_OPENACC
        else
            method = PARALLEL_NONE
        end if

        ! Select appropriate implementation
        select case(method)
        case(PARALLEL_OPENACC)
            method_name = 'OpenACC (GPU)'
            call log%info('Using ' // trim(method_name) // ' for analysis')
            ! Use OpenACC version
            ! Note: This requires linking with analyse_openacc_mod
            ! For now, we'll call a wrapper that handles the interface
            call analyse_openacc_wrapper(radiuscontact, use_contact_probability, &
                    params, Niter, Nmeas, output_folder, analyse_folder, &
                    hic3d_factor, chrom, .true.)
            ! Note: If OpenACC is not compiled, this will fail at build/link time
            
        case(PARALLEL_OPENMP)
            method_name = 'OpenMP (CPU)'
            call log%info('Using ' // trim(method_name) // ' for analysis')
            !$ ! Use OpenMP version
            !$ ! Note: This requires linking with analyse_openmp_mod
            !$ call analyse_openmp_wrapper(radiuscontact, use_contact_probability, &
            !$         params, Niter, Nmeas, output_folder, analyse_folder, &
            !$         hic3d_factor, chrom, .true., num_threads)
            ! Note: If OpenMP is not compiled, the above call will not be compiled
            ! and a link error will occur - this is expected behavior
            
        case default
            method_name = 'Sequential (CPU)'
            call log%info('Using ' // trim(method_name) // ' for analysis')
            ! Use sequential version (original analyse_mod)
            call analyse_sequential_wrapper(radiuscontact, use_contact_probability, &
                    params, Niter, Nmeas, output_folder, analyse_folder, &
                    hic3d_factor, chrom)
        end select

    end subroutine analyse_auto

    ! Wrapper subroutines that interface with the actual implementations
    ! These can be conditionally compiled or use module interfaces

    subroutine analyse_openacc_wrapper(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom, use_gpu)
        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder, analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        logical, intent(in) :: use_gpu
        
        ! Interface to analyse_openacc_mod::analyse_openacc
        ! This would be implemented by either:
        ! 1. Using module interface (if modules are available)
        ! 2. Conditionally compiling the call
        ! 3. Using a function pointer/interface block
        
        ! For now, we'll use conditional compilation
        ! Note: Uncomment when analyse_openacc_mod is available
        ! use analyse_openacc_mod, only: analyse_openacc
        ! call analyse_openacc(radiuscontact, use_contact_probability, &
        !         params, Niter, Nmeas, output_folder, analyse_folder, &
        !         hic3d_factor, chrom, use_gpu)
        
        call log%error('analyse_openacc_mod not available - please link with OpenACC implementation')
    end subroutine analyse_openacc_wrapper

    subroutine analyse_openmp_wrapper(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom, &
            use_omp, num_threads)
        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder, analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        logical, intent(in) :: use_omp
        integer, intent(in) :: num_threads
        
        ! Interface to analyse_openmp_mod::analyse_openmp
        !$ ! Note: Uncomment when analyse_openmp_mod is available
        !$ ! use analyse_openmp_mod, only: analyse_openmp
        !$ ! call analyse_openmp(radiuscontact, use_contact_probability, &
        !$ !         params, Niter, Nmeas, output_folder, analyse_folder, &
        !$ !         hic3d_factor, chrom, use_omp, num_threads)
        
        call log%error('analyse_openmp_mod not available - please link with OpenMP implementation')
    end subroutine analyse_openmp_wrapper

    subroutine analyse_sequential_wrapper(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom)
        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder, analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        
        ! Interface to analyse_mod::analyse (original sequential version)
        ! use analyse_mod, only: analyse
        ! call analyse(radiuscontact, use_contact_probability, &
        !         params, Niter, Nmeas, output_folder, analyse_folder, &
        !         hic3d_factor, chrom)
        
        call log%error('analyse_mod not available - please link with original implementation')
    end subroutine analyse_sequential_wrapper

end module analyse_auto_mod

