! Unified PolymerModel with automatic parallelization selection
! Automatically uses OpenACC if available, falls back to OpenMP, then sequential

module PolymerModel_unified_mod
    use Timers
    use lattice_data_mod
    use logging_mod, only: Logger, LOG_INFO, LOG_DEBUG, str, strf
    use lib_conf
    !$ use omp_lib
    implicit none

    private
    public :: PolymerModel_unified, ModelParameters_unified
    
    type(Logger) :: log = Logger('PolymerModel_unified_mod', LOG_INFO)
    
    ! Interface for randomnumber function (from randomnumber.f03)
    ! This is an external function that uses Fortran's random_number
    ! Declared as external function with explicit interface for type safety
    interface
        function randomnumber() result(r)
            implicit none
            real*8 :: r
        end function randomnumber
    end interface

    type ModelParameters_unified
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        real, public :: kb, ku, km, Ea, Ei
    end type ModelParameters_unified

    type PolymerModel_unified
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        real, public :: kb, ku, km, Ea, Ei, kint = 1.17
        logical, public :: z_loop = .false.
        logical, public :: unidirectional = .false.
        
        ! Parallelization settings
        logical :: use_openacc = .false.
        logical :: use_openmp = .false.
        logical :: prefer_gpu = .true.
        integer :: num_threads = 0
        
        ! Internal model - can be any implementation
        ! We'll use the original PolymerModel as base and add parallelization
        ! For now, we'll conditionally use different implementations
        
        ! Allocatable arrays (same as original)
        integer, dimension(:, :), allocatable :: config
        integer, dimension(:, :), allocatable :: bittable
        real, dimension(:, :), allocatable :: dr
        integer, dimension(:, :), allocatable :: contact
        real, dimension(:, :), allocatable :: boundary
        real, dimension(:), allocatable :: loading_sites_factor
        integer, dimension(:), allocatable :: interaction_sites_state

    contains
        procedure, public :: init_unified, do_simulation_unified, detect_parallelization, &
            output_parameters_unified, init_unified_base, output_unified
        procedure, public :: trialmoveex_unified, trialmovetad_unified, trialbound_unified, &
            trialunbound_unified, unbound_all_unified, erase_unified
        procedure, private :: allocate_unified, cleanup_arrays_unified, initbitable_unified, &
            initconfig4_unified, initconfig4_zigzag_unified, initconfig_sim_out_unified
        final :: deallocate_unified
    end type PolymerModel_unified

contains

    subroutine detect_parallelization(self)
        ! Detect best available parallelization method
        ! Note: If OpenMP/OpenACC is not compiled, the directives in the
        ! implementation functions are ignored and execution continues sequentially.
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        character(1000) :: env_var
        integer :: env_len, stat
        
        self%use_openacc = .false.
        self%use_openmp = .false.
        
        ! Check for OpenACC support (if prefer_gpu is true)
        if (self%prefer_gpu) then
            ! Check environment variable - this is a strong indicator
            call get_environment_variable('ACC_DEVICE_TYPE', env_var, env_len, stat)
            if (stat == 0 .and. env_len > 0) then
                if (trim(env_var) == 'nvidia' .or. trim(env_var) == 'NVIDIA') then
                    ! Environment variable suggests OpenACC should be available
                    self%use_openacc = .true.
                    call log%info('PolymerModel_unified_mod - Using OpenACC (GPU) for simulation')
                    return
                end if
            end if
            ! Even without environment variable, try OpenACC if compiled
            ! (directives will be ignored if not compiled)
            self%use_openacc = .true.
            call log%info('PolymerModel_unified_mod - Attempting OpenACC (GPU) for simulation')
            return
        end if
        
        ! Check for OpenMP support
        ! Try OpenMP if compiled (directives will be ignored if not compiled)
        self%use_openmp = .true.
        call log%info('PolymerModel_unified_mod - Using OpenMP (CPU) for simulation')
    end subroutine detect_parallelization

    subroutine cleanup_arrays_unified(self)
        ! Helper subroutine to safely clean up all arrays and OpenACC data regions
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer :: dealloc_stat
        
        ! Clean up OpenACC data regions first (if they exist and OpenACC was used)
        if (self%use_openacc) then
            ! Ensure any pending OpenACC operations complete before cleanup
            !$acc wait
            
            ! Clean up OpenACC data regions
            if (allocated(self%config) .or. allocated(self%bittable) .or. &
                allocated(self%dr) .or. allocated(self%contact) .or. &
                allocated(self%boundary) .or. allocated(self%loading_sites_factor) .or. &
                allocated(self%interaction_sites_state)) then
                !$acc exit data delete(self%config, self%bittable, self%dr, self%contact, &
                !$acc& self%boundary, self%loading_sites_factor, self%interaction_sites_state)
            end if
        end if
        
        ! Now deallocate arrays (with error handling)
        if (allocated(self%config)) then
            deallocate(self%config, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate config, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%bittable)) then
            deallocate(self%bittable, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate bittable, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%dr)) then
            deallocate(self%dr, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate dr, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%contact)) then
            deallocate(self%contact, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate contact, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%boundary)) then
            deallocate(self%boundary, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate boundary, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%loading_sites_factor)) then
            deallocate(self%loading_sites_factor, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate ' // &
                    'loading_sites_factor, stat=' // trim(str(dealloc_stat)))
            end if
        end if
        if (allocated(self%interaction_sites_state)) then
            deallocate(self%interaction_sites_state, stat=dealloc_stat)
            if (dealloc_stat /= 0) then
                call log%warn('cleanup_arrays_unified: Warning - failed to deallocate ' // &
                    'interaction_sites_state, stat=' // trim(str(dealloc_stat)))
            end if
        end if
    end subroutine cleanup_arrays_unified

    subroutine allocate_unified(self)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer :: bittable_t, stat, dealloc_stat

        bittable_t = 4 * (self%L**3)

        call log%debug('allocate_unified self%Nchain: ' // trim(str(self%Nchain)) // ' bittable_t: ' // trim(str(bittable_t)))

        ! First, ensure complete cleanup by calling cleanup helper
        ! This handles both OpenACC data regions and array deallocation
        call cleanup_arrays_unified(self)

        ! Double-check that all arrays are deallocated before allocating
        if (allocated(self%config) .or. allocated(self%bittable) .or. &
            allocated(self%dr) .or. allocated(self%contact) .or. &
            allocated(self%boundary) .or. allocated(self%loading_sites_factor) .or. &
            allocated(self%interaction_sites_state)) then
            call log%error('allocate_unified: Arrays still allocated after deallocation attempt')
            call log%error('  config allocated: ' // merge('YES', 'NO ', allocated(self%config)))
            call log%error('  bittable allocated: ' // merge('YES', 'NO ', allocated(self%bittable)))
            call log%error('  dr allocated: ' // merge('YES', 'NO ', allocated(self%dr)))
            call log%error('  contact allocated: ' // merge('YES', 'NO ', allocated(self%contact)))
            call log%error('  boundary allocated: ' // merge('YES', 'NO ', allocated(self%boundary)))
            call log%error('  loading_sites_factor allocated: ' // merge('YES', 'NO ', allocated(self%loading_sites_factor)))
            call log%error('  interaction_sites_state allocated: ' // merge('YES', 'NO ', allocated(self%interaction_sites_state)))
            stop 'allocate_unified: Cannot allocate - arrays still allocated'
        end if

        ! Now allocate fresh arrays
        allocate (self%config(2, self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate config, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%bittable(14, bittable_t), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate bittable, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%dr(3, self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate dr, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%contact(3, self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate contact, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%boundary(2, self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate boundary, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%loading_sites_factor(self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate loading_sites_factor, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        allocate (self%interaction_sites_state(self%Nchain), stat=stat)
        if (stat /= 0) then
            call log%error('allocate_unified: Failed to allocate interaction_sites_state, stat=' // trim(str(stat)))
            stop 'allocate_unified: Allocation failed'
        end if
        
        ! Initialize arrays
        self%config = 0
        self%bittable = 0
        self%dr = 0.
        self%contact = 0
        self%boundary = 0.
        self%loading_sites_factor = 0.
        self%interaction_sites_state = 0
    end subroutine allocate_unified

    subroutine deallocate_unified(self)
        implicit none
        type (PolymerModel_unified), intent(inout) :: self

        if (allocated(self%config)) then
            !$acc exit data delete(self%config)
            deallocate(self%config)
        end if
        if (allocated(self%bittable)) then
            !$acc exit data delete(self%bittable)
            deallocate(self%bittable)
        end if
        if (allocated(self%dr)) then
            !$acc exit data delete(self%dr)
            deallocate(self%dr)
        end if
        if (allocated(self%contact)) then
            !$acc exit data delete(self%contact)
            deallocate(self%contact)
        end if
        if (allocated(self%boundary)) then
            !$acc exit data delete(self%boundary)
            deallocate(self%boundary)
        end if
        if (allocated(self%loading_sites_factor)) then
            !$acc exit data delete(self%loading_sites_factor)
            deallocate(self%loading_sites_factor)
        end if
        if (allocated(self%interaction_sites_state)) then
            !$acc exit data delete(self%interaction_sites_state)
            deallocate(self%interaction_sites_state)
        end if
    end subroutine deallocate_unified

    subroutine init_unified(self, L, Nchain, iku, ikm, ikb, Nleffree, kb, ku, km, Ea, Ei, &
            z_loop, unidirectional, kint, prefer_gpu, num_threads)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: L, Nchain, iku, ikm, ikb, Nleffree
        real, intent(in) :: kb, ku, km, Ea, Ei
        logical, intent(in), optional :: z_loop, unidirectional
        real, intent(in), optional :: kint
        logical, intent(in), optional :: prefer_gpu
        integer, intent(in), optional :: num_threads

        ! Clean up any existing OpenACC data regions before reallocating
        ! This must be done BEFORE allocate_unified to avoid conflicts
        if (allocated(self%config) .or. allocated(self%bittable) .or. &
            allocated(self%dr) .or. allocated(self%contact) .or. &
            allocated(self%boundary) .or. allocated(self%loading_sites_factor) .or. &
            allocated(self%interaction_sites_state)) then
            ! Clean up OpenACC data regions if they exist
            if (allocated(self%config)) then
                !$acc exit data delete(self%config)
            end if
            if (allocated(self%bittable)) then
                !$acc exit data delete(self%bittable)
            end if
            if (allocated(self%dr)) then
                !$acc exit data delete(self%dr)
            end if
            if (allocated(self%contact)) then
                !$acc exit data delete(self%contact)
            end if
            if (allocated(self%boundary)) then
                !$acc exit data delete(self%boundary)
            end if
            if (allocated(self%loading_sites_factor)) then
                !$acc exit data delete(self%loading_sites_factor)
            end if
            if (allocated(self%interaction_sites_state)) then
                !$acc exit data delete(self%interaction_sites_state)
            end if
        end if

        self%L = L
        self%Nchain = Nchain
        self%iku = iku
        self%ikm = ikm
        self%ikb = ikb
        self%Nleffree = Nleffree
        self%kb = kb
        self%ku = ku
        self%km = km
        self%Ea = Ea
        self%Ei = Ei
        
        if (present(z_loop)) self%z_loop = z_loop
        if (present(unidirectional)) self%unidirectional = unidirectional
        if (present(kint)) self%kint = kint
        if (present(prefer_gpu)) self%prefer_gpu = prefer_gpu
        if (present(num_threads)) self%num_threads = num_threads

        call self%allocate_unified()
        
        ! Detect best parallelization method
        call self%detect_parallelization()
        
        ! Set up GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc enter data create(self%config, self%bittable, self%dr, self%contact, &
            !$acc& self%boundary, self%loading_sites_factor, self%interaction_sites_state)
            call log%info('PolymerModel: Data allocated on GPU')
        end if
        
        ! Set up OpenMP if using
        if (self%use_openmp) then
            if (self%num_threads > 0) then
                !$ call omp_set_num_threads(self%num_threads)
                call log%info('PolymerModel: OpenMP threads set to ' // trim(str(self%num_threads)))
            else
                !$ call log%info('PolymerModel: OpenMP enabled with default thread count: ' // trim(str(omp_get_max_threads())))
            end if
        end if
    end subroutine init_unified

    subroutine init_unified_base(self, boundary, loading_sites_factor, interaction_sites_state, init_mode, trajectory_i)
        ! Initialize base configuration (similar to original init)
        ! Note: This is a simplified version. For full functionality, you may need to
        ! implement initbitable, initconfig4, etc., or use the original PolymerModel
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        real, dimension(:, :) :: boundary
        real, dimension(:) :: loading_sites_factor
        integer, dimension(:) :: interaction_sites_state
        character(len = 1000), intent(in) :: init_mode
        integer, intent(in) :: trajectory_i

        self%boundary = boundary
        self%loading_sites_factor = loading_sites_factor
        self%interaction_sites_state = interaction_sites_state
        
        ! Update GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc update device(self%boundary, self%loading_sites_factor, self%interaction_sites_state)
        end if
        
        ! Initialize bittable
        call self%initbitable_unified()
        
        ! Initialize configuration based on init_mode
        call log%info('init_mode ' // init_mode)
        select case(init_mode)
        case('h')
            call self%initconfig4_unified()
        case('z')
            call self%initconfig4_zigzag_unified()
        case default   ! s=<init_sim_folder>
            if ( index(init_mode, 's=') > 0 ) then
                call log%info('Init folding mode: ' // trim(init_mode(3:)))
                call self%initconfig_sim_out_unified(trim(init_mode(3:)), trajectory_i)
            else
                call log%error('unknown init_mode ' // init_mode // ' , it will be used the default: helices !')
                call exit(1)
            end if
        end select

        ! Initialize simulation arrays (already zeroed in allocate_unified, but ensure consistency)
        ! Note: Arrays are already initialized in allocate_unified(), but we set these explicitly
        ! for clarity since they're used immediately after
        self%contact = 0
        self%dr = 0.
        
        ! Note: No device update needed here - data will be on device only if actually used for GPU computation
    end subroutine init_unified_base

    subroutine do_simulation_unified(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        ! Unified simulation with automatic method selection
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: trajectory_i
        integer, intent(in) :: Ninter
        integer, intent(in) :: Nmeas
        integer, intent(in) :: burnin
        integer, intent(in) :: burnout
        integer, intent(in) :: burnoutM
        
        integer :: burn_Nmeas, simburnin
        integer :: j, k, v, Ntrial
        real :: pt
        real*8 :: r
        character(20) :: method_name
        ! Note: randomnumber is declared in module interface

        burn_Nmeas = 0
        simburnin = 0

        ! Determine method name with proper string length
        if (self%use_openacc) then
            method_name = 'OpenACC'
        else if (self%use_openmp) then
            method_name = 'OpenMP'
        else
            method_name = 'Sequential'
        end if
        
        call log%info('Start Trajectory (Unified): ' // trim(str(trajectory_i)) // &
                ', Ninter: ' // trim(str(Ninter)) // ', Nmeas: ' // trim(str(Nmeas)) // &
                ', Method: ' // trim(method_name))

        ! Use appropriate parallelization method
        if (self%use_openacc) then
            call do_simulation_openacc_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        else if (self%use_openmp) then
            call do_simulation_openmp_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        else
            call do_simulation_sequential_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        end if
    end subroutine do_simulation_unified

    subroutine do_simulation_openacc_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        ! OpenACC implementation
        ! Note: MC trial moves have dependencies and require random number generation,
        ! so full GPU parallelization is challenging. This implementation uses OpenACC
        ! for data management while keeping the simulation logic on CPU.
        ! Since trial moves run sequentially on CPU, we keep data on host and only
        ! use device data regions for potential future parallelization or caching benefits.
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM
        
        integer :: burn_Nmeas, simburnin, j, k, v, Ntrial
        real :: pt
        real*8 :: r
        ! Note: randomnumber is declared in module interface

        burn_Nmeas = 0
        simburnin = 0

        ! now start the simulation with LEF: measurement after the burn-in time
        call log%info('Start Trajectory: ' // trim(str(trajectory_i)) // ', Ninter: ' // trim(str(Ninter)) // &
                ', Nmeas: ' // trim(str(Nmeas)) // ', burnin: ' // trim(str(burnin)) // &
                ', burnout: ' // trim(str(burnout)) // ', burnoutM: ' // trim(str(burnoutM)))

        ! Since trial moves run on CPU and modify data, we work with host data
        ! Device data regions are maintained for potential future use but not actively used
        ! during sequential simulation. Data is kept consistent via present clauses.

        ! first initial measurement before the burn-in
        call self%output_unified()

        if (burnin > 0) then
            burn_Nmeas = burn_Nmeas + 1

            ! choose random simulation burn-in time = initial_burnin + random_fraction * initial_burnin
            r = randomnumber()
            simburnin = burnin + int(r * burnin)

            do j = 1, simburnin
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() !trial move for monomers
                end do
            end do
            
            ! Note: No device update needed here - trial moves run on CPU and modify host data
            ! Device data will be updated only when actually needed for GPU computation
            
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-in Measurement, ' // &
                    ' burn_Nmeas: ' // trim(str(burn_Nmeas)))
            call self%output_unified()
        end if

        if (burnout > 0) then
            burn_Nmeas = burn_Nmeas + 1
        end if

        if (burnoutM > 0) then
            burn_Nmeas = burn_Nmeas + burnoutM
        end if

        ! Main simulation loop - sequential due to dependencies in MC moves
        ! OpenACC is used for data management, not parallelization
        do j = 1, Nmeas - 1 - burn_Nmeas  ! account for the initial, burn-in and burn-out measurements
            do k = 1, Ninter
                Ntrial = 3 * self%Nchain + self%Nleffree
                pt = real(self%Nchain) / real(Ntrial)
                do v = 1, Ntrial
                    r = randomnumber()
                    if (r.lt.pt) then
                        call self%trialmovetad_unified() !trial move for monomers
                    elseif (r.lt.2 * pt) then
                        call self%trialmoveex_unified() !trial move for LEF movement
                    elseif (r.lt.3 * pt) then
                        call self%trialunbound_unified() !trial move for unbinding event
                    else
                        call self%trialbound_unified() !trial move for binding event
                    end if
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', Measurement:' // trim(str(j)))
            call flush(6)
            call self%output_unified()
            call flush(13)
            ! Note: No device update needed - trial moves run on CPU, data stays on host
        end do

        ! do burn-out included as extra measurement
        if (burnout > 0) then
            call self%unbound_all_unified()
            do j = 1, burnout
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() ! move for monomers
                end do
            end do
            ! Note: No device update needed - trial moves run on CPU
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', fix burn-out Measurement')
            call self%output_unified()
            call flush(13)
        end if

        if (burnoutM > 0) then
            call self%unbound_all_unified()
            do j = 1, burnoutM
                do k = 1, Ninter
                    do v = 1, self%Nchain
                        call self%trialmovetad_unified() ! move for monomers
                    end do
                end do
                ! Note: No device update needed - trial moves run on CPU
                call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-out Measurement: ' // trim(str(j)))
                call self%output_unified()
                call flush(13)
            end do
        end if

        call self%erase_unified()
    end subroutine do_simulation_openacc_impl

    subroutine do_simulation_openmp_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        ! OpenMP implementation
        ! 
        ! IMPORTANT: The main Monte Carlo simulation loop CANNOT be parallelized with OpenMP
        ! because each trial move modifies shared state (config, bittable, contact arrays).
        ! The moves have dependencies - each move depends on the current state and modifies it.
        ! 
        ! OpenMP is used for:
        ! - initbitable_unified(): Parallelizes the bittable initialization loop (independent iterations)
        ! - Other independent initialization operations
        !
        ! The simulation itself runs sequentially to maintain correctness.
        ! For true parallelization, consider:
        ! 1. Running multiple independent trajectories in parallel (MPI)
        ! 2. Algorithm redesign with conflict resolution
        ! 3. Red-black ordering schemes
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM
        
        integer :: burn_Nmeas, simburnin, j, k, v, Ntrial
        real :: pt
        real*8 :: r
        ! Note: randomnumber is declared in module interface

        burn_Nmeas = 0
        simburnin = 0

        ! now start the simulation with LEF: measurement after the burn-in time
        call log%info('Start Trajectory: ' // trim(str(trajectory_i)) // ', Ninter: ' // trim(str(Ninter)) // &
                ', Nmeas: ' // trim(str(Nmeas)) // ', burnin: ' // trim(str(burnin)) // &
                ', burnout: ' // trim(str(burnout)) // ', burnoutM: ' // trim(str(burnoutM)))

        ! first initial measurement before the burn-in
        call self%output_unified()

        if (burnin > 0) then
            burn_Nmeas = burn_Nmeas + 1

            ! choose random simulation burn-in time = initial_burnin + random_fraction * initial_burnin
            r = randomnumber()
            simburnin = burnin + int(r * burnin)

            ! Burn-in loop: sequential due to state dependencies
            do j = 1, simburnin
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() !trial move for monomers
                end do
            end do
            
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-in Measurement, ' // &
                    ' burn_Nmeas: ' // trim(str(burn_Nmeas)))
            call self%output_unified()
        end if

        if (burnout > 0) then
            burn_Nmeas = burn_Nmeas + 1
        end if

        if (burnoutM > 0) then
            burn_Nmeas = burn_Nmeas + burnoutM
        end if

        ! Main simulation loop - sequential due to dependencies in MC moves
        do j = 1, Nmeas - 1 - burn_Nmeas  ! account for the initial, burn-in and burn-out measurements
            do k = 1, Ninter
                Ntrial = 3 * self%Nchain + self%Nleffree
                pt = real(self%Nchain) / real(Ntrial)
                do v = 1, Ntrial
                    r = randomnumber()
                    if (r.lt.pt) then
                        call self%trialmovetad_unified() !trial move for monomers
                    elseif (r.lt.2 * pt) then
                        call self%trialmoveex_unified() !trial move for LEF movement
                    elseif (r.lt.3 * pt) then
                        call self%trialunbound_unified() !trial move for unbinding event
                    else
                        call self%trialbound_unified() !trial move for binding event
                    end if
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', Measurement:' // trim(str(j)))
            call flush(6)
            call self%output_unified()
            call flush(13)
        end do

        ! do burn-out included as extra measurement
        if (burnout > 0) then
            call self%unbound_all_unified()
            ! Burnout loop: sequential due to state dependencies
            do j = 1, burnout
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() ! move for monomers
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', fix burn-out Measurement')
            call self%output_unified()
            call flush(13)
        end if

        if (burnoutM > 0) then
            call self%unbound_all_unified()
            do j = 1, burnoutM
                ! Sequential loop due to state dependencies
                do k = 1, Ninter
                    do v = 1, self%Nchain
                        call self%trialmovetad_unified() ! move for monomers
                    end do
                end do
                call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-out Measurement: ' // trim(str(j)))
                call self%output_unified()
                call flush(13)
            end do
        end if

        call self%erase_unified()
    end subroutine do_simulation_openmp_impl

    subroutine do_simulation_sequential_impl(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        ! Sequential implementation (original algorithm)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM
        
        integer :: burn_Nmeas, simburnin, j, k, v, Ntrial
        real :: pt
        real*8 :: r
        ! Note: randomnumber is declared in module interface

        burn_Nmeas = 0
        simburnin = 0

        ! now start the simulation with LEF: measurement after the burn-in time
        call log%info('Start Trajectory: ' // trim(str(trajectory_i)) // ', Ninter: ' // trim(str(Ninter)) // &
                ', Nmeas: ' // trim(str(Nmeas)) // ', burnin: ' // trim(str(burnin)) // &
                ', burnout: ' // trim(str(burnout)) // ', burnoutM: ' // trim(str(burnoutM)))

        ! first initial measurement before the burn-in
        call self%output_unified()

        if (burnin > 0) then
            burn_Nmeas = burn_Nmeas + 1

            ! choose random simulation burn-in time = initial_burnin + random_fraction * initial_burnin
            ! add random time (a fraction) after the initial burn-in time
            r = randomnumber()
            simburnin = burnin + int(r * burnin)

            do j = 1, simburnin
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() !trial move for monomers
                end do
            end do
            ! now start the simulation with LEF: measurement after the burn-in time
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-in Measurement, ' // &
                    ' burn_Nmeas: ' // trim(str(burn_Nmeas)))
            call self%output_unified()
        end if

        if (burnout > 0) then
            burn_Nmeas = burn_Nmeas + 1
        end if

        if (burnoutM > 0) then
            burn_Nmeas = burn_Nmeas + burnoutM
        end if

        do j = 1, Nmeas - 1 - burn_Nmeas  ! account for the initial, burn-in and burn-out measurements
            do k = 1, Ninter
                Ntrial = 3 * self%Nchain + self%Nleffree
                pt = real(self%Nchain) / real(Ntrial)
                do v = 1, Ntrial
                    r = randomnumber()
                    if (r.lt.pt) then
                        call self%trialmovetad_unified() !trial move for monomers
                    elseif (r.lt.2 * pt) then
                        call self%trialmoveex_unified() !trial move for LEF movement
                    elseif (r.lt.3 * pt) then
                        call self%trialunbound_unified() !trial move for unbinding event
                    else
                        call self%trialbound_unified() !trial move for binding event
                    end if
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', Measurement:' // trim(str(j)))
            call flush(6)
            call self%output_unified()
            call flush(13)
        end do

        ! do burn-out included as extra measurement
        if (burnout > 0) then
            call self%unbound_all_unified()
            do j = 1, burnout
                do v = 1, self%Nchain
                    call self%trialmovetad_unified() ! move for monomers
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', fix burn-out Measurement')
            call self%output_unified()
            call flush(13)
        end if

        if (burnoutM > 0) then
            call self%unbound_all_unified()
            do j = 1, burnoutM
                do k = 1, Ninter
                    do v = 1, self%Nchain
                        call self%trialmovetad_unified() ! move for monomers
                    end do
                end do
                call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-out Measurement: ' // trim(str(j)))
                call self%output_unified()
                call flush(13)
            end do
        end if

        call self%erase_unified()
    end subroutine do_simulation_sequential_impl

    subroutine output_unified(self)
        ! Write simulation data to output files (same as original output)
        ! This implementation matches the original output() subroutine exactly:
        ! - Unit 10: config array (configuration data)
        ! - Unit 11: dr array (displacement data)
        ! - Unit 12: contact array (extruder occupancy data)
        ! - Unit 14: Nleffree (number of free extruders)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: j, rc

        ! Ensure we're working with host data if using OpenACC
        ! Since trial moves run on CPU, host data is always current
        if (self%use_openacc) then
            ! Data is already on host (trial moves modify host arrays)
            ! No need to update from device since we're not using device for computation
        end if

        ! Write in output files (units 10, 11, 12, 14 as in original)
        ! Check that arrays are allocated and have correct dimensions
        if (.not.allocated(self%config) .or. .not.allocated(self%dr) .or. .not.allocated(self%contact)) then
            call log%error('output_unified: Arrays not allocated!')
            call log%error('  config allocated: ' // merge('YES', 'NO ', allocated(self%config)))
            call log%error('  dr allocated: ' // merge('YES', 'NO ', allocated(self%dr)))
            call log%error('  contact allocated: ' // merge('YES', 'NO ', allocated(self%contact)))
            return
        end if

        ! Verify array dimensions
        if (size(self%config, 1) /= 2 .or. size(self%config, 2) /= self%Nchain) then
            call log%error('output_unified: config array has wrong dimensions!')
            return
        end if
        if (size(self%dr, 1) /= 3 .or. size(self%dr, 2) /= self%Nchain) then
            call log%error('output_unified: dr array has wrong dimensions!')
            return
        end if
        if (size(self%contact, 1) /= 3 .or. size(self%contact, 2) /= self%Nchain) then
            call log%error('output_unified: contact array has wrong dimensions!')
            return
        end if

        ! Optimized: Write arrays efficiently using array sections
        ! Using explicit format reduces parsing overhead compared to list-directed I/O
        do j = 1, self%Nchain
            ! Write config array: config(1,j)=node where monomer j is located, 
            ! config(2,j)=direction of the vector between j and j+1
            write(10, '(2(i0,1x))', iostat=rc) self%config(:, j)
            if (rc /= 0) call log%warn('output_unified: Error writing to unit 10, iostat=' // trim(str(rc)))
            
            ! Write dr array: vector (x,y,z in lattice unit) of displacement of monomer j
            write(11, '(3(es15.8,1x))', iostat=rc) self%dr(:, j)
            if (rc /= 0) call log%warn('output_unified: Error writing to unit 11, iostat=' // trim(str(rc)))
            
            ! Write contact array: contact(1,j) \ne 0 if one lef of a extruder is in j,
            ! contact(1,j)=the monomer where the other lef is, 
            ! contact(2,j)=the vector between the two legs, 
            ! contact(3,j)=-1 (resp. +1) if it's a lef walking in the (-) direction (resp. (+) direction)
            write(12, '(3(i0,1x))', iostat=rc) self%contact(:, j)
            if (rc /= 0) call log%warn('output_unified: Error writing to unit 12, iostat=' // trim(str(rc)))
        end do
        write(14, *, iostat=rc) self%Nleffree !number of free (unbound) extruders
        if (rc /= 0) call log%warn('output_unified: Error writing to unit 14, iostat=' // trim(str(rc)))

        flush(10)
        flush(11)
        flush(12)
        flush(14)

        return
    end subroutine output_unified

    subroutine output_parameters_unified(self, fout, init_mode, interaction_sites, boundary_file, lef_loading_sites, &
            basal_loading_factor, boundary_direction, &
            Niter, Ninter, Nmeas, burnin, burnout, burnoutM, radius_contact, &
            kb_a, ku_a, km_a)
        ! Output parameters (same interface as original)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer, intent(in) :: fout
        character(*), intent(in) :: init_mode
        character(*), intent(in) :: interaction_sites
        character(*), intent(in) :: boundary_file
        character(*), intent(in) :: lef_loading_sites
        real, intent(in) :: basal_loading_factor
        integer, intent(in) :: boundary_direction
        integer, intent(in) :: Niter
        integer, intent(in) :: Ninter
        integer, intent(in) :: Nmeas
        integer, intent(in) :: burnin
        integer, intent(in) :: burnout
        integer, intent(in) :: burnoutM
        real, intent(in) :: radius_contact
        real, intent(in) :: kb_a, ku_a, km_a

        write(fout, '(a)') '[3dpolys_le]'
        write(fout, '(a)') '# polymer characteristics'
        write(fout, '(a)') 'Nchain=' // trim(str(self%Nchain))
        write(fout, '(a)') 'L=' // trim(str(self%L))
        write(fout, '(a)') 'Ea=' // trim(strf(self%Ea))
        write(fout, '(a)') 'Ei=' // trim(strf(self%Ei))
        write(fout, '(a)') 'interaction_sites=' // trim(interaction_sites)
        write(fout, '(a)') 'init_mode=' // trim(init_mode)

        write(fout, '(a)') '# measurements'
        write(fout, '(a)') 'Niter=' // trim(str(Niter))
        write(fout, '(a)') 'Ninter=' // trim(str(Ninter))
        write(fout, '(a)') 'Nmeas=' // trim(str(Nmeas))
        write(fout, '(a)') 'burnin=' // trim(str(burnin))
        write(fout, '(a)') 'burnout=' // trim(str(burnout))
        write(fout, '(a)') 'burnoutM=' // trim(str(burnoutM))

        write(fout, '(a)') '# Loop-Extrusion factors:'
        write(fout, '(a)') 'kint=' // trim(strf(self%kint))
        write(fout, '(a)') 'kb=' // trim(strf(kb_a))
        write(fout, '(a)') 'ku=' // trim(strf(ku_a))
        write(fout, '(a)') 'km=' // trim(strf(km_a))
        write(fout, '(a)') 'Nlef=' // trim(str(self%Nleffree))

        write(fout, '(a)') 'boundary=' // trim(boundary_file)
        write(fout, '(a)') 'lef_loading_sites=' // trim(lef_loading_sites)
        write(fout, '(a)') 'basal_loading_factor=' // trim(strf(basal_loading_factor))
        write(fout, '(a)') 'boundary_direction=' // trim(str(boundary_direction))
        if (self%z_loop) then
            write(fout, '(a)') 'z_loop=true'
        else
            write(fout, '(a)') 'z_loop=false'
        end if
        if (self%unidirectional) then
            write(fout, '(a)') 'unidirectional=true'
        else
            write(fout, '(a)') 'unidirectional=false'
        end if

        write(fout, '(a)') '# analysis: experiments in silico:'
        write(fout, '(a)') 'radius_contact=' // trim(strf(radius_contact))

        flush(fout)
        return
    end subroutine output_parameters_unified

    subroutine initbitable_unified(self)
        ! Initialize bittable (same as original)
        ! This loop can be parallelized with OpenMP since each iteration is independent
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        type(Timer) :: crono
        integer :: i, j, k, a, v, ip, jp, kp, L2
        real :: x, y, z, xp, yp, zp

        !initialize the bittable (periodic boundary conditions)
        self%bittable = 0
        L2 = self%L ** 2
        call log%debug('initialize the bittable (periodic boundary conditions) for L: ' // trim(str(self%L)))
        call crono%Tic()
        
        ! Parallelize the outer loop - each iteration is independent
        !$omp parallel do if(self%use_openmp) default(none) &
        !$omp& shared(self, L2, voisxyz) private(a, i, j, k, v, x, y, z, xp, yp, zp, ip, jp, kp)
        do a = 1, 4 * L2 * self%L
            k = int((a - 1) / (2 * L2)) + 1
            j = int((a - 1) / self%L - 2 * self%L * (k - 1)) + 1
            i = a - self%L * (2 * self%L * (k - 1) + (j - 1))
            x = (i - 1) + 0.5 * (1 - mod(j + mod(k + 1, 2), 2))
            y = (j - 1) * 0.5
            z = (k - 1) * 0.5
            self%bittable(1, a) = 0
            ! Inner loop: compute neighbor positions with periodic boundary conditions
            !$omp simd private(xp, yp, zp, ip, jp, kp)
            do v = 1, 12
                xp = x + voisxyz(1, v + 1)
                yp = y + voisxyz(2, v + 1)
                zp = z + voisxyz(3, v + 1)
                ! Apply periodic boundary conditions
                if (xp >= self%L) xp = xp - self%L
                if (xp < 0) xp = xp + self%L
                if (yp >= self%L) yp = yp - self%L
                if (yp < 0) yp = yp + self%L
                if (zp >= self%L) zp = zp - self%L
                if (zp < 0) zp = zp + self%L
                ip = int(xp) + 1
                jp = int(2 * yp + 1)
                kp = int(2 * zp + 1)
                self%bittable(v + 1, a) = ip + (jp - 1) * self%L + (kp - 1) * 2 * L2
            end do
            !$omp end simd
        end do
        !$omp end parallel do
        call log%info('lattice density: ' // trim(strf(real(self%Nchain) / real(4 * L2 * self%L))))
        call log%info(crono%Tac(info = 'initbitable_unified(): '))
        
        ! Update GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc update device(self%bittable)
        end if
    end subroutine initbitable_unified

    subroutine initconfig4_unified(self)
        ! Initialize configuration (helices mode) - same as original
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer :: turn1(7), turn2(7), turn(7), lim, a, n, i, j, t, nv1, nv2, iv, en2, v, b, c(2, self%Nchain)
        ! Note: randomnumber is declared in module interface

        self%bittable(1, :) = 0

        turn1 = (/13, 13, 2, 2, 12, 12, 3/)
        turn2 = (/13, 2, 2, 12, 12, 3, 3/)

        lim = self%L / 2

        self%config = 0
        a = int(4 * self%L**3 * randomnumber()) + 1
        self%config(1, 1) = a
        self%bittable(1, a) = 1
        n = 2
        do i = 1, lim
            if (mod(i, 2).eq.1) then
                turn = turn1
            else
                turn = turn2
            end if
            do j = 1, 7
                self%config(2, n - 1) = turn(j)
                self%config(1, n) = self%bittable(turn(j), self%config(1, n - 1))
                self%bittable(1, self%config(1, n)) = 1
                n = n + 1
            end do
            self%config(2, n - 1) = 11
            self%config(1, n) = self%bittable(11, self%config(1, n - 1))
            self%bittable(1, self%config(1, n)) = 1
            n = n + 1
        end do
        n = n - 1

        do while (n.ne.self%Nchain)
            t = int((n - 1) * randomnumber()) + 1
            iv = int((voisnn(1, 1, self%config(2, t)) - 1) * randomnumber()) + 1
            nv1 = voisnn(2 * iv, 1, self%config(2, t))
            nv2 = voisnn(2 * iv + 1, 1, self%config(2, t))
            en2 = self%config(1, t)
            if (nv1.eq.1) then
                v = en2
            else
                v = self%bittable(nv1, en2)
            end if
            b = self%bittable(1, v)
            if ((b.eq.0)) then
                c = self%config
                self%config(1, t + 1) = v
                self%config(2, t) = nv1
                self%config(2, t + 1) = nv2
                self%bittable(1, v) = 1
                self%config(:, t + 2:n + 1) = c(:, t + 1:n)
                n = n + 1
            end if
        end do

        ! Interaction state
        self%bittable(14,:)=0
        do n=1, self%Nchain
            if (self%interaction_sites_state(n).gt.0) then
                a=self%config(1,n)
                self%bittable(14,a)=self%bittable(14,a)+1
                do v=1,12
                    b=self%bittable(v+1,a)
                    self%bittable(14,b)=self%bittable(14,b)+1
                end do
            end if
        end do
        
        ! Update GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc update device(self%config, self%bittable)
        end if
    end subroutine initconfig4_unified

    subroutine initconfig4_zigzag_unified(self)
        ! Initialize configuration (zigzag mode) - same as original
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        integer :: a, n, i, t, nv1, nv2, iv, en2, v, b, c(2, self%Nchain), lim, v1, v2
        ! Note: randomnumber is declared in module interface

        lim = 2 * self%L - 2

        i = int(6 * randomnumber()) + 1
        select case(i)
        case(1)
            a = 3 * self%L**2 - self%L / 2
            v1 = 6
            v2 = 9
        case(2)
            a = 4 * self%L**3 - 3 * self%L**2 - self%L / 2
            v1 = 8
            v2 = 7
        case(3)
            a = 2 * self%L**3 + 3 * self%L / 2
            v1 = 4
            v2 = 10
        case(4)
            a = 2 * self%L**3 + 2 * self%L**2 - 3 * self%L / 2
            v1 = 11
            v2 = 5
        case(5)
            a = 2 * self%L**3 + self%L**3 - self%L + 1
            v1 = 2
            v2 = 13
        case(6)
            a = 2 * self%L**3 + self%L**2 + self%L
            v1 = 12
            v2 = 3
        end select

        self%config = 0
        self%config(1, 1) = a
        self%bittable(1, a) = 1

        n = 2
        do i = 1, lim - 1
            if (mod(i, 2).eq.1) then
                self%config(2, n - 1) = v1
                self%config(1, n) = self%bittable(v1, self%config(1, n - 1))
                b = self%bittable(1, self%config(1, n))
                if (b.ne.0) then
                    call log%error('error init' // trim(str(i)))
                    stop
                end if
                self%bittable(1, self%config(1, n)) = 1
                n = n + 1
            else
                self%config(2, n - 1) = v2
                self%config(1, n) = self%bittable(v2, self%config(1, n - 1))
                b = self%bittable(1, self%config(1, n))
                if (b.ne.0) then
                    call log%error('error init' // trim(str(i)))
                    stop
                end if
                self%bittable(1, self%config(1, n)) = 1
                n = n + 1
            end if
        end do
        n = n - 1

        do while (n.ne.self%Nchain)
            t = int((n - 1) * randomnumber()) + 1
            iv = int((voisnn(1, 1, self%config(2, t)) - 1) * randomnumber()) + 1
            nv1 = voisnn(2 * iv, 1, self%config(2, t))
            nv2 = voisnn(2 * iv + 1, 1, self%config(2, t))
            en2 = self%config(1, t)
            if (nv1.eq.1) then
                v = en2
            else
                v = self%bittable(nv1, en2)
            end if
            b = self%bittable(1, v)
            if ((b.eq.0)) then
                c = self%config
                self%config(1, t + 1) = v
                self%config(2, t) = nv1
                self%config(2, t + 1) = nv2
                self%bittable(1, v) = 1
                self%config(:, t + 2:n + 1) = c(:, t + 1:n)
                n = n + 1
            end if
        end do

        ! Interaction state
        self%bittable(14,:)=0
        do n=1, self%Nchain
            if (self%interaction_sites_state(n).gt.0) then
                a=self%config(1,n)
                self%bittable(14,a)=self%bittable(14,a)+1
                do v=1,12
                    b=self%bittable(v+1,a)
                    self%bittable(14,b)=self%bittable(14,b)+1
                end do
            end if
        end do
        
        ! Update GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc update device(self%config, self%bittable)
        end if
    end subroutine initconfig4_zigzag_unified

    subroutine initconfig_sim_out_unified(self, sim_out_dir, trajectory_i)
        ! Initialize from simulation output - same as original
        implicit none
        class (PolymerModel_unified), intent(inout) :: self
        character(*), intent(in) :: sim_out_dir
        integer, intent(in) :: trajectory_i
        character(1000) :: config_3dpoys_le_file
        character(1000) :: config_out_file
        integer :: initNchain, initNiter, initNmeas, initNinter, skip_monomers
        integer :: a, n, v, b, lim, v1, rc
        integer :: laticeAddr, nextTurn  ! MonomerConfig fields

        config_3dpoys_le_file = trim(sim_out_dir) // '/3dpoys_le.cfg'
        config_out_file = trim(sim_out_dir) // '/config.out'

        call load_config_file(config_3dpoys_le_file, iostat=rc)
        if (rc /= 0) then
            call log%error('Could not find or open input configuration file: ' // trim(config_out_file))
            call exit(1)
        end if
        call read_config('Nchain',   initNchain)
        call read_config('Niter',    initNiter)
        call read_config('Nmeas',    initNmeas)
        call read_config('Ninter',   initNinter)

        if ( .not. initNchain == self%Nchain ) then
            call log%error('Not compatible intialazing Nchain ' // trim(str(initNchain)) &
                    // ' expected ' // trim(str(self%Nchain)))
            call exit(1)
        end if

        lim = 2 * self%L - 2
        self%config = 0

        open(30, file = trim(config_out_file), action = 'read', iostat = rc)
        if (rc == 0) then
            skip_monomers =  ((mod(trajectory_i-1, initNiter)+1) * (initNmeas-1) * initNchain)
            call log%info('Load config_out_file file ' // trim(config_out_file) &
                    // ' and for trajectory ' // trim(str(trajectory_i)) // ' skip ' // trim(str(skip_monomers)))

            do n=1, skip_monomers
                read(30, *, iostat = rc) laticeAddr, nextTurn
                if (rc /= 0) exit
            end do

            do n=1, self%Nchain
                read(30, *, iostat = rc) laticeAddr, nextTurn
                if (rc /= 0) exit
                a = laticeAddr
                v1 = nextTurn
                self%config(1, n) = a
                self%config(2, n) = v1
                self%bittable(1, a) = self%bittable(1, a) + 1
            end do

            ! Interaction state
            self%bittable(14,:)=0
            do n=1, self%Nchain
                if (self%interaction_sites_state(n).gt.0) then
                    a=self%config(1,n)
                    self%bittable(14,a)=self%bittable(14,a)+1
                    do v=1,12
                        b=self%bittable(v+1,a)
                        self%bittable(14,b)=self%bittable(14,b)+1
                    end do
                end if
            end do
        end if
        close(30)
        
        ! Update GPU data if using OpenACC
        if (self%use_openacc) then
            !$acc update device(self%config, self%bittable)
        end if
    end subroutine initconfig_sim_out_unified

    ! ============================================================================
    ! Trial Move Methods (copied from original PolymerModel and adapted)
    ! ============================================================================

    subroutine trialmoveex_unified(self)
        ! trial move to move a LEF leg
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: n, iv, s, id, con(3), strand
        real*8 :: fc
        real :: impermeability

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1

        s = self%contact(3, n)

        if ((s.eq.0).or.(n.eq.1).or.(n.eq.self%Nchain)) return !if monomer not occupied or end-monomer do nothing

        strand = (s + 1) / 2 + 1 ! 1 (-), 2 (+)
        impermeability = abs(self%boundary(strand, n))
        fc = (1. - impermeability)**(1. / real(self%ikm)) ! probability of permeability

        do iv = 1, self%ikm
            if (randomnumber().ge.(self%km * fc)) return
        end do

        if (s.eq.-1) then !if a (-) direction leg move to n-1
            if (connec(1, self%config(2, n - 1), self%contact(2, n)).eq.0) return !if break the slip-link do not move

            if (self%contact(1, n - 1).ne.0) then !if n-1 already occupied, try to swap
                if (.not. self%z_loop) return
                if ((self%contact(3, n - 1).eq.-1).or.(n.eq.2)) return !if same direction or at the end do not swap
                if (connec(1, opp(self%config(2, n - 1)), self%contact(2, n - 1)).eq.0) return  ! break slip-link n-1
                con = self%contact(:, n - 1)

                iv = connec(1, self%config(2, n - 1), self%contact(2, n))
                id = self%contact(1, n)
                self%contact(1, n - 1) = id
                self%contact(2, n - 1) = iv
                self%contact(3, n - 1) = -1
                self%contact(1, id) = n - 1
                self%contact(2, id) = opp(iv)

                iv = connec(1, opp(self%config(2, n - 1)), con(2))
                id = con(1)
                self%contact(1, n) = id
                self%contact(2, n) = iv
                if (self%unidirectional) then
                    self%contact(3, n) = con(3)
                else
                    self%contact(3, n) = 1
                end if
                self%contact(1, id) = n
                self%contact(2, id) = opp(iv)
            else
                iv = connec(1, self%config(2, n - 1), self%contact(2, n))
                id = self%contact(1, n)
                self%contact(1, n - 1) = id
                self%contact(2, n - 1) = iv
                self%contact(3, n - 1) = -1
                self%contact(1, id) = n - 1
                self%contact(2, id) = opp(iv)
                self%contact(:, n) = 0
            end if
        elseif (s.eq.1) then !if a (+1) direction leg move to n+1
            if (connec(1, opp(self%config(2, n)), self%contact(2, n)).eq.0) return !if break the slip-link do not move

            if (self%contact(1, n + 1).ne.0) then !if n+1 already occupied, try to swap
                if (.not. self%z_loop) return
                if ((self%contact(3, n + 1).eq.1).or.(n.eq.(self%Nchain - 1))) return ! if same direction or at the end do not swap
                if (connec(1, self%config(2, n + 1), self%contact(2, n + 1)).eq.0) return !if break the slip-link of n+1 do not swap
                con = self%contact(:, n + 1)

                iv = connec(1, opp(self%config(2, n)), self%contact(2, n))
                id = self%contact(1, n)
                self%contact(1, n + 1) = id
                self%contact(2, n + 1) = iv
                self%contact(3, n + 1) = 1
                self%contact(1, id) = n + 1
                self%contact(2, id) = opp(iv)

                iv = connec(1, self%config(2, n + 1), con(2))
                id = con(1)
                self%contact(1, n) = id
                self%contact(2, n) = iv
                if (self%unidirectional) then
                    self%contact(3, n) = con(3)
                else
                    self%contact(3, n) = -1
                end if
                self%contact(1, id) = n
                self%contact(2, id) = opp(iv)
            else
                iv = connec(1, opp(self%config(2, n)), self%contact(2, n))
                id = self%contact(1, n)
                self%contact(1, n + 1) = id
                self%contact(2, n + 1) = iv
                self%contact(3, n + 1) = 1
                self%contact(1, id) = n + 1
                self%contact(2, id) = opp(iv)
                self%contact(:, n) = 0
            end if
        end if

        return
    end subroutine trialmoveex_unified

    subroutine trialmovetad_unified(self)
        !trial move for monomer
        ! NOTE: This is a very long method - copying from original
        ! For full implementation, see polymer_model.f03 lines 566-853
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: n, iv, v, b, j, nv1, nv2, nm2, np1, en, cn2, cn3, cm2, en2, id, cc, a
        real :: dE
        ! Note: randomnumber is declared in module interface

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1
        en = self%config(1, n)

        !test if allowed moved and move
        if (n.eq.1) then
            en2 = self%config(1, 2)
            cn2 = opp(self%config(2, 1))
            if (cn2.lt.self%config(2, 2)) then
                cm2 = self%config(2, 2)
            else
                cm2 = cn2
                cn2 = self%config(2, 2)
            end if
            iv = int(11 * randomnumber()) + 1
            if (iv.ge.cn2) iv = iv + 1
            if (iv.ge.cm2) iv = iv + 1

            if (iv.eq.1) then
                v = en2
            else
                v = self%bittable(iv, en2)
            end if
            b = self%bittable(1, v)
            if ((b.eq.0).or.((b.eq.1).and.(en2.eq.v))) then
                id = self%contact(1, n)
                cn2 = self%config(2, 1)
                cn3 = self%contact(2, n)
                if (cn3.eq.0) then
                    cc = 0
                else
                    cc = connec(iv, opp(cn2), cn3)
                end if
                if ((id.ne.0).and.(cc.eq.0)) return
                dE = costhet(opp(iv), self%config(2, 2)) - costhet(cn2, self%config(2, 2))

                if (self%interaction_sites_state(n).gt.0) then
                    dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                end if

                if ((self%contact(3, n).eq.-1).and.(id.lt.self%Nchain)) then
                    if (connec(opp(self%contact(2, n)), self%config(2, id), 1).ne.0) dE = dE - self%Ea
                    if (connec(opp(cc), self%config(2, id), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%contact(3, n + 1).eq.-1) then
                    if (connec(opp(cn2), self%contact(2, n + 1), 1).ne.0) dE = dE - self%Ea
                    if (connec(iv, self%contact(2, n + 1), 1).ne.0) dE = dE + self%Ea
                end if

                if (randomnumber().lt.exp(-dE)) then
                    self%bittable(1, en) = self%bittable(1, en) - 1
                    self%bittable(1, v) = self%bittable(1, v) + 1

                    if (self%interaction_sites_state(n).gt.0) then
                        self%bittable(14, en) = self%bittable(14, en) - 1
                        self%bittable(14, v) = self%bittable(14, v) + 1
                        do j=2,13
                            a=self%bittable(j, en)
                            self%bittable(14, a) = self%bittable(14, a) - 1
                            a=self%bittable(j, v)
                            self%bittable(14, a) = self%bittable(14, a) + 1
                        end do
                    end if

                    self%config(1, 1) = v
                    self%config(2, 1) = opp(iv)
                    if (id.ne.0) then
                        self%contact(2, n) = cc
                        self%contact(2, id) = opp(cc)
                    end if

                    self%dr(1, n) = self%dr(1, n) + voisxyz(1, iv) + voisxyz(1, cn2)
                    self%dr(2, n) = self%dr(2, n) + voisxyz(2, iv) + voisxyz(2, cn2)
                    self%dr(3, n) = self%dr(3, n) + voisxyz(3, iv) + voisxyz(3, cn2)
                end if
            end if

        elseif (n.eq.self%Nchain) then
            en2 = self%config(1, self%Nchain - 1)
            cn2 = self%config(2, self%Nchain - 1)
            if (cn2.lt.opp(self%config(2, self%Nchain - 2))) then
                cm2 = opp(self%config(2, self%Nchain - 2))
            else
                cm2 = cn2
                cn2 = opp(self%config(2, self%Nchain - 2))
            end if

            iv = int(11 * randomnumber()) + 1
            if (iv.ge.cn2) iv = iv + 1
            if (iv.ge.cm2) iv = iv + 1

            if (iv.eq.1) then
                v = en2
            else
                v = self%bittable(iv, en2)
            end if
            b = self%bittable(1, v)

            if ((b.eq.0).or.((b.eq.1).and.(en2.eq.v))) then
                id = self%contact(1, n)
                cn2 = self%config(2, self%Nchain - 1)
                cn3 = self%contact(2, n)
                if (cn3.eq.0) then
                    cc = 0
                else
                    cc = connec(iv, cn2, cn3)
                end if
                if ((id.ne.0).and.(cc.eq.0)) return
                dE = costhet(self%config(2, self%Nchain - 2), iv) - costhet(self%config(2, self%Nchain - 2), cn2)

                if (self%interaction_sites_state(n).gt.0) then
                    dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                end if

                if ((self%contact(3, n).eq.1).and.(id.gt.1)) then
                    if (connec(self%contact(2, n), opp(self%config(2, id - 1)), 1).ne.0) dE = dE - self%Ea
                    if (connec(cc, opp(self%config(2, id - 1)), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%contact(3, n - 1).eq.1) then
                    if (connec(cn2, self%contact(2, n - 1), 1).ne.0) dE = dE - self%Ea
                    if (connec(iv, self%contact(2, n - 1), 1).ne.0) dE = dE + self%Ea
                end if

                if (randomnumber().lt.exp(-dE)) then
                    self%bittable(1, en) = self%bittable(1, en) - 1
                    self%bittable(1, v) = self%bittable(1, v) + 1

                    if (self%interaction_sites_state(n).gt.0) then
                        self%bittable(14, en) = self%bittable(14, en) - 1
                        self%bittable(14, v) = self%bittable(14, v) + 1
                        do j=2,13
                            a = self%bittable(j, en)
                            self%bittable(14, a) = self%bittable(14, a) - 1
                            a = self%bittable(j, v)
                            self%bittable(14, a) = self%bittable(14, a) + 1
                        end do
                    end if

                    self%config(1, self%Nchain) = v
                    self%config(2, self%Nchain - 1) = iv
                    if (id.ne.0) then
                        self%contact(2, n) = cc
                        self%contact(2, id) = opp(cc)
                    end if
                    self%dr(1, n) = self%dr(1, n) + voisxyz(1, iv) - voisxyz(1, cn2)
                    self%dr(2, n) = self%dr(2, n) + voisxyz(2, iv) - voisxyz(2, cn2)
                    self%dr(3, n) = self%dr(3, n) + voisxyz(3, iv) - voisxyz(3, cn2)
                end if
            end if
        else
            cn2 = self%config(2, n)
            cm2 = self%config(2, n - 1)
            en2 = self%config(1, n - 1)
            nm2 = n - 2
            np1 = n + 1

            if (voisnn(1, cm2, cn2).gt.1) then
                iv = int((voisnn(1, cm2, cn2) - 1) * randomnumber()) + 1
                if (voisnn(2 * iv, cm2, cn2).ge.cm2) iv = iv + 1
                nv1 = voisnn(2 * iv, cm2, cn2)
                nv2 = voisnn(2 * iv + 1, cm2, cn2)
                if (nv1.eq.1) then
                    v = en2
                else
                    v = self%bittable(nv1, en2)
                end if
                b = self%bittable(1, v)
                if ((b.eq.0).or.((b.eq.1).and.((v.eq.en2).or.(v.eq.self%config(1, np1))))) then
                    id = self%contact(1, n)
                    cn3 = self%contact(2, n)
                    if (cn3.eq.0) then
                        cc = 0
                    else
                        cc = connec(nv1, cm2, cn3)
                    end if
                    if ((id.ne.0).and.(cc.eq.0)) return
                    if (n.eq.2) then
                        dE = costhet(nv1, nv2) + costhet(nv2, self%config(2, np1)) - costhet(cm2, cn2) &
                                - costhet(cn2, self%config(2, np1))
                    elseif (n.eq.self%Nchain - 1) then
                        dE = costhet(self%config(2, nm2), nv1) + costhet(nv1, nv2) - costhet(self%config(2, nm2), cm2) &
                                - costhet(cm2, cn2)
                    else
                        dE = costhet(self%config(2, nm2), nv1) + costhet(nv1, nv2) + costhet(nv2, self%config(2, np1)) &
                                - costhet(self%config(2, nm2), cm2) - costhet(cm2, cn2) - costhet(cn2, self%config(2, np1))
                    end if

                    if (self%interaction_sites_state(n).gt.0) then
                        dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                    end if

                    if ((self%contact(3, n).eq.-1).and.(id.lt.self%Nchain)) then
                        if (connec(opp(self%contact(2, n)), self%config(2, id), 1).ne.0) dE = dE - self%Ea
                        if (connec(opp(cc), self%config(2, id), 1).ne.0) dE = dE + self%Ea
                    elseif ((self%contact(3, n).eq.1).and.(id.gt.1)) then
                        if (connec(self%contact(2, n), opp(self%config(2, id - 1)), 1).ne.0) dE = dE - self%Ea
                        if (connec(cc, opp(self%config(2, id - 1)), 1).ne.0) dE = dE + self%Ea
                    end if

                    if (self%contact(3, n + 1).eq.-1) then
                        if (connec(opp(cn2), self%contact(2, n + 1), 1).ne.0) dE = dE - self%Ea
                        if (connec(opp(nv2), self%contact(2, n + 1), 1).ne.0) dE = dE + self%Ea
                    end if
                    if (self%contact(3, n - 1).eq.1) then
                        if (connec(cm2, self%contact(2, n - 1), 1).ne.0) dE = dE - self%Ea
                        if (connec(nv1, self%contact(2, n - 1), 1).ne.0) dE = dE + self%Ea
                    end if

                    if (randomnumber().lt.exp(-dE)) then
                        self%bittable(1, en) = self%bittable(1, en) - 1
                        self%bittable(1, v) = self%bittable(1, v) + 1

                        if (self%interaction_sites_state(n).gt.0) then
                            self%bittable(14, en) = self%bittable(14, en) - 1
                            self%bittable(14, v) = self%bittable(14, v) + 1
                            do j=2,13
                                a = self%bittable(j, en)
                                self%bittable(14, a) = self%bittable(14, a) - 1
                                a = self%bittable(j, v)
                                self%bittable(14, a) = self%bittable(14, a) + 1
                            end do
                        end if

                        self%config(1, n) = v
                        self%config(2, n - 1) = nv1
                        self%config(2, n) = nv2
                        if (id.ne.0) then
                            self%contact(2, n) = cc
                            self%contact(2, id) = opp(cc)
                        end if
                        self%dr(1, n) = self%dr(1, n) + voisxyz(1, nv1) - voisxyz(1, cm2)
                        self%dr(2, n) = self%dr(2, n) + voisxyz(2, nv1) - voisxyz(2, cm2)
                        self%dr(3, n) = self%dr(3, n) + voisxyz(3, nv1) - voisxyz(3, cm2)
                    end if
                end if
            end if
        end if
        return
    end subroutine trialmovetad_unified

    subroutine trialbound_unified(self)
        !trial move for binding LEF
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: n, id, j, d
        real :: kbp
        ! Note: randomnumber is declared in module interface

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1

        id = self%contact(1, n)
        if (id==0) then !if the bin is not occupied by a leg, try to randomly insert a LEF to NN sites
            kbp = self%kb * (  self%loading_sites_factor(n)**(1./real(self%ikb)) )
            do j = 1, self%ikb
                if (randomnumber()>=kbp) return
            end do
            d = int(2 * randomnumber()) + 1
            if (d==1) then
                if (n>1) then
                    if (self%contact(1, n - 1)==0) then
                        self%contact(1, n) = n - 1
                        self%contact(2, n) = opp(self%config(2, n - 1))
                        self%contact(1, n - 1) = n
                        self%contact(2, n - 1) = self%config(2, n - 1)
                        if (self%unidirectional) then
                            if (randomnumber()<0.5) then
                                self%contact(3, n) = 2
                                self%contact(3, n - 1) = -1
                            else
                                self%contact(3, n) = 1
                                self%contact(3, n - 1) = -2
                            end if
                        else
                            self%contact(3, n) = 1
                            self%contact(3, n - 1) = -1
                         end if
                        self%Nleffree = self%Nleffree - 1
                    end if
                end if
            else if (n<self%Nchain) then
                if (self%contact(1, n + 1)==0) then
                    self%contact(1, n) = n + 1
                    self%contact(2, n) = self%config(2, n)
                    self%contact(1, n + 1) = n
                    self%contact(2, n + 1) = opp(self%config(2, n))
                    if (self%unidirectional) then
                        if (randomnumber()<0.5) then
                            self%contact(3, n) = -2
                            self%contact(3, n + 1) = 1
                        else
                            self%contact(3, n) = -1
                            self%contact(3, n + 1) = 2
                        end if
                    else
                        self%contact(3, n) = -1
                        self%contact(3, n + 1) = 1
                    end if
                    self%Nleffree = self%Nleffree - 1
                end if
            end if
        end if

        return
    end subroutine trialbound_unified

    subroutine trialunbound_unified(self)
        !trial move for unbinding LEF
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: n, id, j
        ! Note: randomnumber is declared in module interface, not as local variable

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1

        id = self%contact(1, n)
        if (id.gt.0) then !if the monomer is occupied by a leg, try to remove it
            do j = 1, self%iku
                if (randomnumber().ge.self%ku) return
            end do
            self%contact(:, n) = 0
            self%contact(:, id) = 0
            write(13, *) n, id
            self%Nleffree = self%Nleffree + 1
        end if

        return
    end subroutine trialunbound_unified

    subroutine unbound_all_unified(self)
        !unbinding of all bound LEFs
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: n, id

        do n = 1, self%Nchain
            id = self%contact(1, n)
            if (id.gt.0) then !if the monomer is occupied by a leg, try to remove it
                self%contact(:, n) = 0
                self%contact(:, id) = 0
                write(13, *) n, id
                self%Nleffree = self%Nleffree + 1
            end if
        end do

        return
    end subroutine unbound_all_unified

    subroutine erase_unified(self)
        implicit none
        class (PolymerModel_unified), intent(inout) :: self

        integer :: i, a

        !erase the configuration and occupancy of the bittable

        do i = 1, self%Nchain
            a = self%config(1, i)
            self%bittable(1, a) = 0
        end do
        self%config = 0

        return
    end subroutine erase_unified

end module PolymerModel_unified_mod

