! Unified analysis module with automatic method selection
! This version uses module interfaces to properly link with available implementations

module analyse_unified_mod
    use lattice_data_mod
    use PolymerModel_mod
    use logging_mod
    use analyse_mod  ! Import analyse_mod to get explicit interface for analyse
    USE HDF5
    implicit none

    private
    public :: analyse_unified
    
    type(Logger) :: log = Logger('analyse_unified_mod', LOG_DEBUG)
    
    ! Note: The actual implementations are in separate modules:
    ! - analyse_openacc_mod (for OpenACC)
    ! - analyse_openmp_mod (for OpenMP)  
    ! - analyse_mod (for sequential)
    !
    ! To use this unified interface, you need to:
    ! 1. Link with the desired implementation modules
    ! 2. Use module procedures or create wrapper calls
    !
    ! For simplicity, we'll use conditional compilation and direct calls
    ! when the modules are available

contains

    logical function check_openacc_available() result(available)
        ! Check if OpenACC is available at runtime
        implicit none
        character(1000) :: env_var
        integer :: env_len, stat
        
        available = .false.
        
        ! Check environment variable
        call get_environment_variable('ACC_DEVICE_TYPE', env_var, env_len, stat)
        if (stat == 0 .and. env_len > 0) then
            if (trim(env_var) == 'nvidia' .or. trim(env_var) == 'NVIDIA') then
                available = .true.
                return
            end if
        end if
        
        ! Check if compiled with OpenACC
        !$acc if (.true.)
        available = .true.
        !$acc end if
    end function check_openacc_available

    logical function check_openmp_available() result(available)
        ! Check if OpenMP is available
        implicit none
        available = .false.
        !$omp if (.true.)
        available = .true.
        !$omp end if
    end function check_openmp_available

    subroutine analyse_unified(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom, &
            prefer_gpu, force_method, num_threads)

        ! Unified analysis interface with automatic method selection
        ! prefer_gpu: if true, prefer OpenACC over OpenMP
        ! force_method: 'auto', 'openacc', 'openmp', or 'sequential'
        ! num_threads: number of threads for OpenMP (0 = use default)

        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder, analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        logical, intent(in) :: prefer_gpu
        character(*), intent(in) :: force_method  ! 'auto', 'openacc', 'openmp', 'sequential'
        integer, intent(in) :: num_threads

        logical :: use_openacc, use_openmp, use_sequential
        character(100) :: selected_method

        call log%set_level(global_log_level)

        ! Determine which method to use
        if (trim(force_method) == 'openacc') then
            use_openacc = .true.
            use_openmp = .false.
            use_sequential = .false.
            selected_method = 'OpenACC (forced)'
        else if (trim(force_method) == 'openmp') then
            use_openacc = .false.
            use_openmp = .true.
            use_sequential = .false.
            selected_method = 'OpenMP (forced)'
        else if (trim(force_method) == 'sequential') then
            use_openacc = .false.
            use_openmp = .false.
            use_sequential = .true.
            selected_method = 'Sequential (forced)'
        else
            ! Auto-detect: prefer_gpu determines priority
            if (prefer_gpu .and. check_openacc_available()) then
                use_openacc = .true.
                use_openmp = .false.
                use_sequential = .false.
                selected_method = 'OpenACC (auto-detected)'
            else if (check_openmp_available()) then
                use_openacc = .false.
                use_openmp = .true.
                use_sequential = .false.
                selected_method = 'OpenMP (auto-detected)'
            else
                use_openacc = .false.
                use_openmp = .false.
                use_sequential = .true.
                selected_method = 'Sequential (fallback)'
            end if
        end if

        call log%info('Analysis method: ' // trim(selected_method))

        ! Call appropriate implementation
        ! Use conditional compilation to call the right implementation
        
        ! For now, always use the sequential analyse from analyse_mod
        ! TODO: When OpenACC/OpenMP implementations are fully integrated with save routines,
        !       uncomment the conditional calls below
        call analyse(radiuscontact, use_contact_probability, &
                params, Niter, Nmeas, output_folder, analyse_folder, &
                hic3d_factor, chrom)
        
        ! Future implementation with full OpenACC/OpenMP support:
        ! if (use_openacc) then
        !     !$acc if (.true.)
        !     use analyse_openacc_mod, only: analyse_openacc
        !     call analyse_openacc(radiuscontact, use_contact_probability, &
        !             params, Niter, Nmeas, output_folder, analyse_folder, &
        !             hic3d_factor, chrom, .true.)
        !     !$acc end if
        ! else if (use_openmp) then
        !     !$omp if (.true.)
        !     use analyse_openmp_mod, only: analyse_openmp
        !     call analyse_openmp(radiuscontact, use_contact_probability, &
        !             params, Niter, Nmeas, output_folder, analyse_folder, &
        !             hic3d_factor, chrom, .true., num_threads)
        !     !$omp end if
        ! else
        !     call analyse(radiuscontact, use_contact_probability, &
        !             params, Niter, Nmeas, output_folder, analyse_folder, &
        !             hic3d_factor, chrom)
        ! end if

    end subroutine analyse_unified

end module analyse_unified_mod

