character(len = 1) function path_separator()
    character(len = 10000) :: path
    call get_environment_variable('PATH', path)
    path_separator = path(1:1)
    ! print*, 'path=', trim(path), ' path_separator=', path_separator
end function path_separator

character(len = 1000) function find_path_program()
    character(len = 10000) :: path
    integer :: i, j
    character(*), parameter :: program_bin = '3dpolys_le/bin' !not optimal to look for: if instalation folder changes
    !character(*), parameter :: program_bin = 'anaconda3/bin' ! used just for local test

    call get_environment_variable('PATH', path)
    j = index(path, program_bin, .true.)
    if (j>0) then
        ! print*, 'j=', j, ' path(1:j)=', path(1:j)
        i = index(path(1:j), ':', .true.)

        find_path_program = path(i + 1:j + len(program_bin) - 1)
    end if
end function find_path_program

subroutine print_version()
    character(*), parameter :: VERSION = '2022.5'
    character(1000) :: program_location = './', find_path_program
    character(1000) :: program_folder
    character(25) :: var_name, program_name = '3dpolys_le', program__version = '2022.5'
    character(2) :: eq_sign = '='
    character(1) :: path_separator, path_sep
    logical :: file_exists
    integer :: i, rc
    character(1000) :: init__py_file = '__version__.py'

    path_sep = path_separator()

    CALL get_command_argument(0, program_location)
    if (len_trim(program_location) > 0) then
        ! print*, ' program_location ', trim(program_location)

        i = index(program_location, path_sep, .true.)
        if (i==0) then
            ! probably in PATH environment variable: try to find it by folder name
            program_folder = find_path_program()
        else
            ! used the whome program path or simlink
            program_folder = program_location(:i-1)
        end if
        if (len(trim(program_folder))>0) then
            ! print*, ' program_folder ', trim(program_folder)

            init__py_file = adjustl(trim(program_folder)) // path_sep // '..' // path_sep // '__init__.py'
            ! print*, trim(init__py_file)

            inquire(file = trim(init__py_file), exist = file_exists)
            if (file_exists) then
                open(10, file = trim(init__py_file), action = 'read', iostat = rc)
                if (rc == 0) then
                    read(10, *, iostat = rc) var_name, eq_sign, program_name
                    read(10, *, iostat = rc) var_name, eq_sign, program__version
                end if
            end if
        end if

        print*, program_name, ' version ', VERSION

    end if
end subroutine print_version

subroutine print_help()
    call print_version()
    print*, 'Usage: 3dpolys_le [options] [<3dpolys_le.cfg file>]'
    print*, 'Options:'
    print*, '   -h|--help   Display this information.'
    print*, '   --log:<log level> Sets the log output level: OFF, FATAL, ERROR, WARN, INFO, DEBUG, TRACE. Default: INFO'
    print*, '   --hic3d:<hic3d_factor>    Produce hic3d-map with the given factor (resolution) on the last measurement &
            & otherwise do not. Deafault: not'
    print*, '   -o|--output_folder:<output folder> Path to a simulation output folder. Default: the folder of the &
            & 3dpolys_le.cfg file.'
    print*, '   -a|--analyse:<analyse folder> Perform analyse step on an already done simulation&
            & and store results in a given folder.'
    print*, '   -c|--chrom:Chromosome name.'
    print*, '   -b|--boundary:<boundary sites file> Boundary sites file in a csv format with the following columns:&
            & name,midpoint,impermeability. Default: no bondaries'
    print*, '   -lls|--lef_loading_sites:<loop extrusion loading sites file> LEFs loadding sites file in a csv format&
            & with the following columns: name,position,length,probability. Default: if not given, the whole polymer'
    print*, '   -r|--radius_contact:<radius> contact radius in lattice unit (1=70nm) for extracting Hi-C matrixes, Default: 1.42'
    print*, '   -cp|--contact_probability Together with the contact radius use contact probability with the formula&
            & (1 - <actual_radius_contact>^2 / <radius_contact param>^2),&
            & Default: false - uniform distributed.'
    print*, '   -l|-nlef:<Nlef_val>: Nlef value, overwriting the one from the inpiut.dat'
    print*, '   -m|--km:<km_val>: km value, overwriting the one from the inpiut.dat'
    print*, '   -bd|--boundary_direction:<boundary_direction>: impermeability direction applied to all boundaries:&
            & -1:opposite direction, 0:both, 1:same direction. Default: 0'
    print*, '   -im|--init_mode:<init_mode>: Initial folding mode: h for helices-like, z for zigzag-like polymer state.'
    print*, '   -z|--z_loop : Allow z_loop for LEFs move, where LEFs can traverse one another. Default: false'
    print*, '   -u|--unidirectional : Unidirectional mode for LEFs move otherwise bidirectional. Default: false=bidirectional'
    print*, '<3dpolys_le.cfg file>: path to the inpit.dat file. Default: ./3dpolys_le.cfg'
    print*, '<output folder>: path to output folder. Default: the folder of the <3dpolys_le.cfg file>'
    print*, 'inpiut.cfg format:'
    print*, '<name>=<value>'
    print*, 'parameters in a configuration file (3dpolys_le.cfg):'
    print*, 'Nchain     Polymer chain length in monomers of 2kb.'
    print*, 'L          Polymer compartment box size L (choose L so that Nchain/(4*L^3) ~ 0.5).'
    print*, 'Niter      Number of iterations, aka number of independent trajectories as polymer replicas.'
    print*, 'Nmeas      Number of measures >=3 (initial, burin-in, n*simulation steps, burn-out), aka number of snapshots.'
    print*, 'Ninter     Interval between measures, aka number of Monte Carlo steps (MCS) between two snapshots. &
            & 1 min = 12000 interaction steps'
    print*, 'kint       Bending energy of the polymer (do not change).'
    print*, 'kb         LEFs binding rate (do not change).'
    print*, 'ku         LEFs half-unbinding rate (do not change).'
    print*, 'km         LEFs rate of movement, aka velocity (2.7e-3 = 50kb/min=833bp/s).'
    print*, 'Ea         Energy of extrusion (passive: Ea=0, facilitated Ea<0).'
    print*, 'Nlef       Maximal number of bound LEFs (average number of bound extruders=Nlef kb/(kb+2ku).'
    print*, 'chrom      Chromosome name to be used for output files. Default: chrS'
    print*, 'burnin     Number of iterations steps before introducing LEFs = burnin + random(0,1)*burnin.'
    print*, 'burnout    Number of iterations steps after the last simulation measurement with removed LEFs, &
            & aka relaxation time in steps.'
    print*, 'burnoutM   Number of measurements at the end with removed LEFs, aka relaxation time in measurements.'
end subroutine print_help

subroutine check_iostat(io_unit, iostat)
    integer, intent(in) :: io_unit
    integer, intent(in) :: iostat
    if (iostat /= 0) then
        close(io_unit)
        call print_help()
        call exit(1)
    end if
end subroutine check_iostat

program mainprogram
    use Timers
    use lattice_data_mod
    use PolymerModel_mod
    use analyse_mod
    use mpi
    use logging_mod
    use lib_conf

    implicit none

    character(1000) :: input_options
    character(1000) :: input_dat_file
    character(1000) :: save_input_cfg_file
    character(1000) :: input_folder
    character(1000) :: output_folder = ''
    character(1000) :: boundary_file = ''
    character(1000) :: lef_loading_sites = ''
    character(20) :: opt_s = ''
    character(20) :: chrom = ''
    integer :: ai = 1  ! input argument position the 3dpolys_le.cfg in the CLI, the ai+1 is the output folder
    character(1) :: path_separator, path_sep
    character(1) :: init_mode = ''

    integer :: L, Nchain, Niter, Nmeas, Ninter, iku, ikm, ikb, Nlef = 0, burnin = 0, burnout = 0, burnoutM = 0
    !integer :: simburnin = 0, burn_Nmeas = 0
    real :: kint, kb, ku, km = 0., Ea, kb_factor

    integer :: i, time
    !real :: pt
    !real*8 :: randomnumber, r
    ! real :: seed
    type(Timer) :: crono
    real, dimension(:, :), allocatable :: boundary  ! (strand -/+, permeability)
    real, dimension(:), allocatable :: loading_sites_factor
    integer :: loading_sites_count = 0 ! deafault: there is no lef_loading_site file so the whole polymer is loading site
    real :: basal_loading_factor = -1. ! not defined
    integer :: boundary_direction = -9 ! not defned direction
    integer, parameter :: resolution_factor = 2000
    type(PolymerModel) :: model
    type(ModelParameters) :: params
    real :: radius_contact = 0 !in lattice unit (recall: 1 lattice unit=70nm)
    logical :: use_contact_probability = .false.
    logical :: do_analyse = .false.
    character(1000) :: analyse_folder = ''

    integer :: rank, Nproc, ierr
    integer :: rank_Niter, mod_Niter
    integer :: trajectory_i
    !real*8 :: start, finish
    integer*4 :: status = 0
    integer :: rc
    character(len = 20) :: col1, col2, col3, col4
    type(Logger) :: log = Logger(source = '3dpolys_le', level = LOG_INFO)
    integer :: hic3d_factor = 0
    logical :: z_loop = .false.
    logical :: unidirectional = .false.
    logical :: file_exists

    type :: BoundarySite
        character(len = 20) :: name
        integer :: midpoint
        real :: impermeability
    end type BoundarySite
    type(BoundarySite) :: boundary_site

    type :: LoadingSite
        character(len = 20) :: name
        integer :: position
        integer :: length
        real :: factor
    end type LoadingSite
    type(LoadingSite) :: loading_site

    real :: impermeability_prev
    integer :: strand

    path_sep = path_separator()

    call MPI_Init(ierr)
    if (ierr /= 0) error stop 'mpi init error'

    !>  Get the number of processes.
    call MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    if (rank == 0) then
        call print_version()

        ! print *, compiler_version()
        ! wtime = MPI_Wtime()
        call log%info('number of processes: ' // trim(str(Nproc)))
        !print *, 'MPI LIBRARY_VERSION: , ', len
    end if

    !load parameters
    global_log_level = LOG_INFO  ! default
    if (command_argument_count() > 0) then
        CALL get_command_argument(ai, input_options)
        i = index(input_options, '-')
        do while (i == 1)
            ! call log%debug('argument i: ' // trim(str(i)))
            if (index(input_options, '--log:') > 0) then
                i = index(input_options, ':')
                global_log_level = str2loglevel(trim(input_options(i + 1:)))
                !print*, 'log_level', global_log_level
                call log%set_level(global_log_level)
            elseif (index(input_options, '--hic3d:') > 0) then
                i = index(input_options, ':')
                opt_s = trim(input_options(i + 1:))
                READ(opt_s, *) hic3d_factor
                if (rank == 0) then
                    call log%info('Produce hic3d-map on the last measurement with factor: ' // trim(str(hic3d_factor)))
                end if
            elseif ((index(input_options, '--analyse:') > 0).or.(index(input_options, '-a:') > 0)) then
                do_analyse = .true.
                i = index(input_options, ':')
                analyse_folder = trim(input_options(i + 1:))
                ! add last / if missing
                if (analyse_folder(len(analyse_folder):) /= path_sep) then
                    analyse_folder = trim(analyse_folder) // path_sep
                end if
                if (rank == 0) then
                    call log%info('Do only analyse on already done simulation and store results in ' // trim(analyse_folder))
                end if
            elseif ((index(input_options, '--boundary:') > 0).or.(index(input_options, '-b:') > 0)) then
                i = index(input_options, ':')
                boundary_file = trim(input_options(i + 1:))
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(boundary_file))
                end if
            elseif ((index(input_options, '--chrom:') > 0).or.(index(input_options, '-c:') > 0)) then
                i = index(input_options, ':')
                chrom = trim(input_options(i + 1:))
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(chrom))
                end if
            elseif ((index(input_options, '--lef_loading_sites:') > 0).or.(index(input_options, '-lls:') > 0)) then
                i = index(input_options, ':')
                lef_loading_sites = trim(input_options(i + 1:))
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(lef_loading_sites))
                end if
            elseif ((index(input_options, '--radius_contact:') > 0).or.(index(input_options, '-r:') > 0)) then
                i = index(input_options, ':')
                opt_s = trim(input_options(i + 1:))
                READ(opt_s, *) radius_contact
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(strf(radius_contact)))
                end if
            elseif ((index(input_options, '--contact_probability') > 0).or.(index(input_options, '-cp') > 0)) then
                use_contact_probability = .true.
                if (rank == 0) then
                    call log%info('Together with the monomer contact radius use contact probability: (1 - r^2 / max_r^2)')
                end if
            elseif ((index(input_options, '--output_folder:') > 0).or.(index(input_options, '-o:') > 0)) then
                i = index(input_options, ':')
                output_folder = trim(input_options(i + 1:))
                ! add last / if missing
                if (output_folder(len(output_folder):) /= path_sep) then
                    output_folder = trim(output_folder) // path_sep
                end if
                if (rank == 0) then
                    call log%info('Store all results in output folder: ' // trim(output_folder))
                end if
            elseif ((index(input_options, '--nlef:') > 0).or.(index(input_options, '-l:') > 0)) then
                i = index(input_options, ':')
                opt_s = trim(input_options(i + 1:))
                READ(opt_s, *) Nlef
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(str(Nlef)))
                end if
            elseif ((index(input_options, '--km:') > 0).or.(index(input_options, '-m:') > 0)) then
                i = index(input_options, ':')
                opt_s = trim(input_options(i + 1:))
                READ(opt_s, *) km
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(strf(km)))
                end if
            elseif ((index(input_options, '--boundary_direction:') > 0).or.(index(input_options, '-bd:') > 0)) then
                i = index(input_options, ':')
                opt_s = trim(input_options(i + 1:))
                READ(opt_s, *) boundary_direction
                if (rank == 0) then
                    call log%info('input ' // trim(input_options(:i)) // trim(str(boundary_direction)))
                end if
            elseif ((index(input_options, '--z_loop') > 0).or.(index(input_options, '-z') > 0)) then
                z_loop = .true.
                if (rank == 0) then
                    call log%info('Allow z_loop for LEFs move')
                end if
            elseif ((index(input_options, '--unidirectional') > 0).or.(index(input_options, '-u') > 0)) then
                unidirectional = .true.
                if (rank == 0) then
                    call log%info('Unidirectional mode for LEFs move')
                end if
            elseif ((index(input_options, '--init_mode') > 0).or.(index(input_options, '-im') > 0)) then
                i = index(input_options, ':')
                init_mode = trim(input_options(i + 1:i + 2))
                if (rank == 0) then
                    call log%info('Init folding mode: ' // init_mode)
                end if
            elseif ((index(input_options, '--help') > 0).or.(index(input_options, '-h') > 0)) then
                call print_help()
                call exit(0)
            else
                if (rank == 0) then
                    call log%info('unrecognized option: ' // trim(input_options))
                end if
                call print_help()
                call exit(0)
            end if
            ai = ai + 1
            if ((command_argument_count() >= ai)) then
                CALL get_command_argument(ai, input_options)
            else
                exit
            end if
            i = index(input_options, '-')   ! next argument
        end do
        CALL get_command_argument(ai, input_dat_file)
    else
        input_dat_file = './3dpolys_le.cfg'
    end if

    call log%debug('MPI running rank: ' // trim(str(rank)))
    !start = MPI_Wtime()

    !initialize random generator
    call SYSTEM_CLOCK(time)
    call srand(time)
    call random_seed()  ! use PUT=seed_int
    ! call random_number(seed)
    ! seed = 7 ! used only for synchronize compare with the original code
    !call initializerandomnumbergenerator(dble(seed))

    i = index(input_dat_file, path_sep, .true.)
    if (i>=0) then
        input_folder = input_dat_file(:i)
    else
        input_folder = './'
    end if

    if (rank == 0) then
        call log%info('Load initial parameters from ' // trim(input_dat_file) // &
                ' file and input_folder: ' // trim(input_folder))
    end if

    if (len(output_folder) == 0) then
        output_folder = input_folder
    end if

    if (rank == 0) then
        call log%info('Output folder: ' // trim(output_folder))
    end if

    call load_config_file(input_dat_file, iostat=ierr)
    if (ierr /= 0) then
        if (rank == 0) then
            call log%error('Could not find or open input configuration file: ' // trim(input_dat_file))
            call print_help()
        end if
        call exit(1)
    end if
    call read_config('Nchain',    Nchain)
    call read_config('L',    L)
    call read_config('Niter',    Niter)
    call read_config('Nmeas',    Nmeas)
    call read_config('Ninter',   Ninter)
    ! call read_config('kint',     kint) ! not used
    call read_config('kb',       kb)
    call read_config('ku',       ku)
    if (km == 0.) then
        call read_config('km',   km)
    end if
    call read_config('Ea',       Ea)
    if (Nlef == 0.) then
        call read_config('Nlef', Nlef)
    end if
    call read_config('burnin',   burnin)
    call read_config('burnout',  burnout)
    call read_config('burnoutM', burnoutM)
    if (init_mode == '') then
        call read_config('init_mode', init_mode)
    end if
    if (boundary_file == '') then
        call read_config('boundary', boundary_file, default='')
    end if
    if (chrom == '') then
        call read_config('chrom', chrom, default='chrS')
    end if
    if (lef_loading_sites == '') then
        call read_config('lef_loading_sites', lef_loading_sites, default='')
    end if
    if (basal_loading_factor == -1) then
        call read_config('basal_loading_factor', basal_loading_factor, default=1.)
    end if
    if (hic3d_factor == 0) then
        call read_config('hic3d_factor', hic3d_factor, default=0)
    end if
    if (boundary_direction < -1) then
        call read_config('boundary_direction', boundary_direction)
    end if
    if (.not.z_loop) then
        call read_config('z_loop', z_loop)
    end if
    if (.not.unidirectional) then
        call read_config('unidirectional', unidirectional)
    end if
    if (radius_contact == 0.) then
        call read_config('radius_contact', radius_contact, default=1.42)
    end if

    ! loaded input parameters:
    if (rank == 0) then
        call log%info('Nchain=' // trim(str(Nchain)  ))
        call log%info('L=' // trim(str(L)  ))
        call log%info('Niter=' // trim(str(Niter)  ))
        call log%info('Nmeas=' // trim(str(Nmeas)  ))
        call log%info('Ninter=' // trim(str(Ninter) ))
        call log%info('kint=' // trim(strf(kint)   ))
        call log%info('kb=' // trim(strf(kb)     ))
        call log%info('ku=' // trim(strf(ku)     ))
        call log%info('km=' // trim(strf(km)     ))
        call log%info('Ea=' // trim(strf(Ea)     ))
        call log%info('Nlef=' // trim(str(Nlef)   ))
        call log%info('burnin=' // trim(str(burnin) ))
        call log%info('burnout=' // trim(str(burnout)))
        call log%info('burnoutM=' // trim(str(burnoutM)))

        call log%info('init_mode=' // trim(init_mode))
        call log%info('boundary=' // trim(boundary_file))
        call log%info('lef_loading_sites=' // trim(lef_loading_sites))
        call log%info('basal_loading_factor=' // trim(strf(basal_loading_factor)))
        call log%info('boundary_direction=' // trim(str(boundary_direction)))
        if (z_loop) then
            call log%info('z_loop=true')
        else
            call log%info('z_loop=false')
        end if
        if (unidirectional) then
            call log%info('unidirectional=true')
        else
            call log%info('unidirectional=false')
        end if
        call log%info('radius_contact=' // trim(strf(radius_contact)))
        call log%info('chrom=' // trim(chrom))
    end if

    !load the local state of the 2kbp-bins
    ! used to init a polymermodel.boundary
    allocate (boundary(2, Nchain))
    boundary = 0.

    if (trim(boundary_file) /= '') then
        open(10, file = trim(boundary_file), action = 'read', iostat = rc)
        if (rc == 0) then
            if (rank == 0) then
                call log%info('Load boundaries file ' // trim(boundary_file) // ' ...')
            end if
            col1 = ''
            col2 = ''
            col3 = ''
            read(10, *) col1, col2, col3
            if ((trim(col1) == 'name').and.(trim(col2)=='midpoint').and.(trim(col3)=='impermeability')) then
                do
                    read(10, *, iostat = rc) boundary_site
                    if (rc /= 0) exit
                    !print*, boundary_site
                    i = 1 + boundary_site%midpoint / resolution_factor ! indexing starts from 1
                    if ((boundary_direction * boundary_site%impermeability) < 0) then
                        strand = 1
                    else if ((boundary_direction * boundary_site%impermeability) > 0) then
                        strand = 2
                    else
                        strand = 0  ! both
                    end if
    !                call log%info('Nchain: ' // trim(str(Nchain)) // ', boundary_site%impermeability: ' &
    !                        // trim(strf(boundary_site%impermeability)) // ', strand: ' // trim(str(strand)) &
    !                        // ', i: ' // trim(str(i)) &
    !                        // ', boundary_site%midpoint: ' // trim(str(boundary_site%midpoint)))
                    if ((i <= Nchain) .and. (boundary_site%impermeability /= 0)) then
                        impermeability_prev = 0
                        ! print*, 'boundary(' // trim(str(strand)) // ', ' // trim(str(i)) // ')=' // &
                        !        trim(strf(boundary(max(1, strand), i)))
                        if (boundary(max(1, strand), i) /= 0) then
                            impermeability_prev = boundary(max(1, strand), i)
                            if (rank == 0) then
                                call log%info('update boundary-site[' // trim(str(strand)) // ',' // trim(str(i)) // ']=' // &
                                        trim(boundary_site%name) // ', prev.impermeability: ' // &
                                        trim(strf(impermeability_prev)) // ' with ' // trim(strf(boundary_site%impermeability)))
                            end if
                        end if

                        ! print*, 'old impermeability_prev: ' // trim(strf(impermeability_prev))
                        impermeability_prev = 1 - (1 - impermeability_prev) * (1 - abs(boundary_site%impermeability))
                        if (strand /= 0) then
                            boundary(strand, i) = impermeability_prev
                        else
                            boundary(1, i) = impermeability_prev
                            boundary(2, i) = impermeability_prev
                        end if
                        if (rank == 0) then
                            call log%info('set boundary-site[' // trim(str(strand)) // ',' // trim(str(i)) // ']=' // &
                                    trim(boundary_site%name) // ', impermeability: ' // &
                                    trim(strf(boundary(max(1, strand), i))))
                        end if
                    else
                        if (rank == 0) then
                            call log%warn('skip boundary-site[' // trim(str(strand)) // ',' // trim(str(i)) // ']=' // &
                                    trim(boundary_site%name) // ') out of boundary size or impermeability = 0')
                        end if
                    end if
                end do
            else
                if (rank == 0) then
                    call log%error('wrong format for boundary sites: ' // col1 // col2 // col3 // ' but expected: ' // &
                            'name midpoint impermeability. STOP PROCEEDING. Please provide a correct boundary.csv file.')

                end if
                call MPI_FINALIZE(ierr)
                if (ierr /= 0) error stop 'mpi finalize error'

                call exit(1)
            end if
        else
            if (rank == 0) then
                call log%error('Cound not find ' // trim(boundary_file) // ' file!')
            end if
            call MPI_FINALIZE(ierr)
            if (ierr /= 0) error stop 'mpi finalize error'

            call exit(1)
        end if
        close(10)
    end if

    allocate (loading_sites_factor(Nchain))
    loading_sites_factor = basal_loading_factor

    if (trim(lef_loading_sites) /= '') then
        open(10, file = trim(lef_loading_sites), action = 'read', iostat = rc)
        if (rc == 0) then
            if (rank == 0) then
                call log%info('Load lef_loading_sites file ' // trim(lef_loading_sites) // ' ...')
            end if
            col1 = ''
            col2 = ''
            col3 = ''
            col4 = ''
            read(10, *) col1, col2, col3, col4
            ! print "(A, A, A)", 'col1[', trim(col1), ']'
            ! TODO add (trim(col1) == "name").and.
            if ((trim(col2)=="position").and.(trim(col3)=="length")) then
                !&
                !.and.(trim(col4)=='factor')) then

                loading_sites_count = 0
                kb_factor = 0
                do
                    read(10, *, iostat = rc) loading_site
                    if (rc /= 0) exit
                    do i = (1 + loading_site%position/resolution_factor), &
                            (1 + loading_site%position/resolution_factor + (loading_site%length - 1)/resolution_factor)
                        loading_sites_count = loading_sites_count + 1
                        loading_sites_factor(i) = loading_site%factor
                        kb_factor = kb_factor + loading_site%factor
                        if (rank == 0) then
                            call log%info('Added  LEFs loading site: loading_sites_factor(' // trim(str(i)) &
                                    // ')= ' // trim(strf(loading_site%factor)))
                        end if
                    end do
                end do
                ! correct kb
                kb_factor = real(Nchain) / (kb_factor + real(Nchain-loading_sites_count)*basal_loading_factor)
                ! should be same == real(Nchain) / SUM(loading_sites_factor)
                kb = kb * kb_factor
                if (rank == 0) then
                    call log%info('Added ' // trim(str(loading_sites_count)) // ' LEFs binding sites.')
                    call log%info('Loading kb_factor= ' // trim(strf(kb_factor)))
                    call log%info('New kb= ' // trim(strf(kb)))
                end if
            else
                if (rank == 0) then
                    call log%error('wrong format for lef binding sites: ' // col1 // col2 // col3 // col4 // ' but expected: ' // &
                            'name position length probability. STOP PROCEEDING. Please provide a correct loading_sites.csv file.')

                end if
                call MPI_FINALIZE(ierr)
                if (ierr /= 0) error stop 'mpi finalize error'

                call exit(1)
            end if
        else
            if (rank == 0) then
                call log%error('Cound not find ' // trim(lef_loading_sites) // ' file!')
            end if
            call MPI_FINALIZE(ierr)
            if (ierr /= 0) error stop 'mpi finalize error'

            call exit(1)
        end if
        close(10)
    end if

    !the first and last monomer should be impenetrable
    boundary(1, 1) = 1.
    boundary(2, 1) = 1.
    boundary(1, Nchain) = 1.
    boundary(2, Nchain) = 1.
    do i = 2, Nchain - 1
        boundary(1, i) = min(1., max(0., boundary(1, i))) !(should be between 0 and 1)
        boundary(2, i) = min(1., max(0., boundary(2, i))) !(should be between 0 and 1)
    end do

    !tricks to effectively reduce the rate
    ikm = 1
    if (km.gt.0.) then
        do while (km**(1. / real(ikm)).lt.0.001)
            ikm = ikm + 1
        end do
    end if
    ikb = 1
    if (kb.gt.0.) then
        do while (kb**(1. / real(ikb)).lt.0.001)
            ikb = ikb + 1
        end do
    end if
    iku = 1
    if (ku.gt.0.) then
        do while (ku**(1. / real(iku)).lt.0.001)
            iku = iku + 1
        end do
    end if
    km = km**(1. / real(ikm))
    ku = ku**(1. / real(iku))
    kb = kb**(1. / real(ikb))

    !load and predefine lattice properties
    call lattice_init(kint)

    if (rank >= Niter) then
        rank_Niter = 0
    elseif (Niter <= Nproc) then
        rank_Niter = 1
    else
        rank_Niter = Niter / Nproc ! run only rank_Niter Trajectory for a rank
        ! take care of the last rank: probably not the same rank_Niter
        mod_Niter = mod(Niter, Nproc)
        if ((mod_Niter /= 0).and.(rank < mod_Niter)) then
            rank_Niter = rank_Niter + 1
        end if
    end if

    call log%info('Running simulations for rank:' // trim(str(rank)) // ' #trajectories:' // trim(str(rank_Niter)) // ' ...')

    params = ModelParameters(L = L, Nchain = Nchain, iku = iku, ikm = ikm, ikb = ikb, Nleffree = Nlef, &
            kb = kb, ku = ku, km = km, Ea = Ea)

    if ((rank_Niter > 0).and.(.not.do_analyse)) then
        status = SYSTEM('mkdir -p ' // trim(output_folder))

        !open output files
        ! TODO consistent open/write/flush/close files
        open(10, file = trim(output_folder) // 'config' // '_' // trim(str(rank)) // '.out', &
                action = 'write', status = 'new', iostat = rc)
        open(11, file = trim(output_folder) // 'dr' // '_' // trim(str(rank)) // '.out', &
                action = 'write', status = 'new', iostat = rc)
        open(12, file = trim(output_folder) // 'contact' // '_' // trim(str(rank)) // '.out', &
                action = 'write', status = 'new', iostat = rc)
        open(13, file = trim(output_folder) // 'process' // '_' // trim(str(rank)) // '.out', &
                action = 'write', status = 'new', iostat = rc)
        open(14, file = trim(output_folder) // 'Nlef' // '_' // trim(str(rank)) // '.out', &
                action = 'write', status = 'new', iostat = rc)

        call log%debug('L:' // trim(str(L)) // ', Nchain:' // trim(str(Nchain)) // ' Nchain:' // trim(str(Niter)) &
                // ' Nmeas:' // trim(str(Nmeas)))

        do i = 1, rank_Niter
            trajectory_i = rank * rank_Niter + i
            !call crono%Tic()

            !generate initial configuration
            model = PolymerModel(L = L, Nchain = Nchain, iku = iku, ikm = ikm, ikb = ikb, Nleffree = Nlef, &
                    kb = kb, ku = ku, km = km, Ea = Ea, z_loop = z_loop, unidirectional = unidirectional)

            if ((rank == 0).and.(i == 1)) then                ! do it only once
                save_input_cfg_file = trim(trim(output_folder) // '3dpoys_le.cfg')
                call log%info('Save parameters in file: ' // save_input_cfg_file)

                open(20, file = save_input_cfg_file, action = 'write', status = 'new', iostat = rc)
                call model%output_parameters(20, init_mode, boundary_file, lef_loading_sites, &
                        basal_loading_factor, boundary_direction, &
                        Niter, Ninter, Nmeas, burnin, burnout, burnoutM, radius_contact)
                close(20)
            end if

            call model%init(boundary, loading_sites_factor, init_mode)

            call model%do_simulation(trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)

            !call log%info(crono%Tac(info = ' tajectory ' // trim(str(trajectory_i)) // ' for rank ' // trim(str(rank))))
        end do

        close(10)
        close(11)
        close(12)
        close(13)
        close(14)

        analyse_folder = output_folder
    endif

    call log%debug('WAIT MPI_Barrier for rank ' // trim(str(rank)) // '...')
    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (ierr /= 0) then
        call log%error('FAILED MPI_Barrier for rank ' // trim(str(rank)))
        ! error stop 'mpi barrier error'
    end if
    call log%debug('PASSED MPI_Barrier for rank ' // trim(str(rank)))

    if (rank == 0) then
        call crono%Tic()
        ! TODO take care of old files or remove them before generating the new one
        inquire(file = trim(trim(output_folder) // 'config.out'), exist = file_exists)
        if (.NOT.file_exists) then
            call log%info('Merge files to collect all results together from ' // trim(output_folder) // ' ...')
            status = SYSTEM('cat ' // trim(output_folder) // 'config_*.out > ' // trim(output_folder) // 'config.out')
            if (status .eq. 0) CALL SYSTEM('rm ' // trim(output_folder) // 'config_*.out')   ! maybe keep config files but in one file are again more then one!
            status = SYSTEM('cat ' // trim(output_folder) // 'dr_*.out > ' // trim(output_folder) // 'dr.out')
            if (status == 0) CALL SYSTEM('rm ' // trim(output_folder) // 'dr_*.out')
            status = SYSTEM('cat ' // trim(output_folder) // 'contact_*.out > ' // trim(output_folder) // 'contact.out')
            if (status == 0) CALL SYSTEM('rm ' // trim(output_folder) // 'contact_*.out')
            status = SYSTEM('cat ' // trim(output_folder) // 'process_*.out > ' // trim(output_folder) // 'process.out')
            if (status == 0) CALL SYSTEM('rm ' // trim(output_folder) // 'process_*.out')
            status = SYSTEM('cat ' // trim(output_folder) // 'Nlef_*.out > ' // trim(output_folder) // 'Nlef.out')
            if (status == 0) CALL SYSTEM('rm ' // trim(output_folder) // 'Nlef_*.out')
            call log%info(crono%Tac('Merged files'))
        end if
    end if

    if ((rank == 0).and.do_analyse) then
        call log%info('MPI root rank final analyse ...')

        status = SYSTEM('mkdir -p ' // trim(analyse_folder))

        call crono%Tic()
        call analyse(radiuscontact = radius_contact, use_contact_probability = use_contact_probability, &
                params = params, Niter = Niter, Nmeas = Nmeas, &
                output_folder = output_folder, analyse_folder = analyse_folder, &
                hic3d_factor = hic3d_factor, chrom = chrom)
        call log%info(crono%Tac('Finished analyse'))
    end if

    !finish = MPI_Wtime()
    !write(6, '(i3,f7.4)') rank, finish - start
    call log%debug('Finished MPI rank ' // trim(str(rank)))

    call MPI_FINALIZE(ierr)
    if (ierr /= 0) error stop 'mpi finalize error'

end program
