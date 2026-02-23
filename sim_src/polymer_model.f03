! Created by  on 19.09.19.

module PolymerModel_mod
    use Timers
    use lattice_data_mod
    use logging_mod
    use lib_conf
    use repro_rng_mod, only: repro_rng_use_portable, repro_rng_init_from_trajectory, repro_rng_next

    implicit none
    private
    interface
        function randomnumber()
            real(8) :: randomnumber
        end function randomnumber
    end interface
    public :: PolymerModel, ModelParameters, MonomerConfig
    type(Logger) :: log = Logger('PolymerModel_mod', LOG_INFO)

    type ModelParameters
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        real, public :: kb, ku, km, Ea, Ei
        ! integer, public :: Niter, Nmeas, Ninter, kint  ! not used
    end type ModelParameters

    type MonomerConfig
        integer :: laticeAddr
        integer :: nextTurn  ! TODO rename: link
    end type MonomerConfig

    type PolymerModel
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        real, public :: kb, ku, km, Ea, Ei, kint = 1.17
        logical, public :: z_loop = .false.
        logical, public :: unidirectional = .false.
        ! allocatable
        integer, dimension(:, :), allocatable :: config
        integer, dimension(:, :), allocatable :: bittable
        real, dimension(:, :), allocatable :: dr
        integer, dimension(:, :), allocatable :: contact
        real, dimension(:, :), allocatable :: boundary
        real, dimension(:), allocatable :: loading_sites_factor
        integer, dimension(:), allocatable :: interaction_sites_state
        ! Per-replica RNG: portable (seed set) or compiler (no seed)
        integer, dimension(:), allocatable, private :: rng_state
        integer(8), dimension(:), allocatable, private :: rng_state_rep

    contains
        procedure, public :: init, do_simulation, advance_one_measurement_interval, trialmoveex, trialmovetad, &
                trialbound, trialunbound, unbound_all, erase, output, output_parameters
        procedure, public :: run_burnin_phase, run_burnout_only_phase, run_burnoutM_one_phase
        procedure, public :: export_state_flat, import_state_flat
        procedure, private :: get_random
        procedure :: allocate, initbitable, initconfig4, initconfig4_zigzag, initconfig_sim_out
        final :: deallocate
    end type PolymerModel

contains

    subroutine init(self, boundary, loading_sites_factor, interaction_sites_state, init_mode, trajectory_i)
        implicit none
        class (PolymerModel), intent(inout) :: self
        real, dimension(:, :) :: boundary
        real, dimension(:) :: loading_sites_factor
        integer, dimension(:) :: interaction_sites_state
        character(len = 1000), intent(in) :: init_mode
        integer, intent(in) :: trajectory_i
        character(len = 1000) :: init_sim_folder

        call log%set_level(global_log_level)

        call self%allocate()

        self%boundary = boundary
        self%loading_sites_factor = loading_sites_factor
        self%interaction_sites_state = interaction_sites_state
        !print*, 'initialize boundary ', shape(boundary), shape(self%boundary)

        !display interaction states
        !do i = 1, self%Nchain
        !    if (self%interaction_sites_state(i) > 0) then
        !        call log%info('interaction_site['// trim(str(i)) // ']='// trim(str(self%interaction_sites_state(i))))
        !    end if
        !end do

        call self%initbitable()
        call log%info('init_mode ' // init_mode)
        select case(init_mode)
        case('h')
            call self%initconfig4()
        case('z')
            call self%initconfig4_zigzag()
        case default   ! s=<init_sim_folder>
            if ( index(init_mode, 's=') > 0 ) then
                init_sim_folder = trim(init_mode(3:))
                call log%info('Init folding mode: ' // init_sim_folder)
                call self%initconfig_sim_out(init_sim_folder, trajectory_i)
            else
                call log%error('unknown init_mode ' // init_mode // ' , it will be used the default: helices !')
                ! call self%initconfig4()
                call exit(1)
            end if
        end select

        !make simulations
        self%contact = 0
        self%dr = 0.

        ! Per-replica RNG: portable (GPU/CPU match when seed set) or compiler state
        if (repro_rng_use_portable) then
            if (allocated(self%rng_state)) deallocate(self%rng_state)
            if (allocated(self%rng_state_rep)) deallocate(self%rng_state_rep)
            allocate(self%rng_state_rep(1))
            call repro_rng_init_from_trajectory(trajectory_i, self%rng_state_rep(1))
        else
            block
                integer :: sz
                call random_seed(size=sz)
                if (allocated(self%rng_state_rep)) deallocate(self%rng_state_rep)
                if (allocated(self%rng_state)) deallocate(self%rng_state)
                allocate(self%rng_state(sz))
                call random_seed(get=self%rng_state)
            end block
        end if

        return
    end subroutine

    subroutine allocate(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer :: bittable_t

        ! include 'global.f03'

        ! Deallocate if already allocated (e.g. re-init for next trajectory in single process)
        if (allocated(self%config)) deallocate(self%config)
        if (allocated(self%bittable)) deallocate(self%bittable)
        if (allocated(self%dr)) deallocate(self%dr)
        if (allocated(self%contact)) deallocate(self%contact)
        if (allocated(self%boundary)) deallocate(self%boundary)
        if (allocated(self%loading_sites_factor)) deallocate(self%loading_sites_factor)
        if (allocated(self%interaction_sites_state)) deallocate(self%interaction_sites_state)
        if (allocated(self%rng_state)) deallocate(self%rng_state)
        if (allocated(self%rng_state_rep)) deallocate(self%rng_state_rep)

        bittable_t = 4 * (self%L**3)

        call log%debug('allocate self%Nchain: ' // trim(str(self%Nchain)) // ' bittable_t: ' // trim(str(bittable_t)))

        allocate (self%config(2, self%Nchain))
        allocate (self%bittable(14, bittable_t))
        allocate (self%dr(3, self%Nchain))
        allocate (self%contact(3, self%Nchain))
        allocate (self%boundary(2, self%Nchain))
        allocate (self%loading_sites_factor( self%Nchain))
        allocate (self%interaction_sites_state( self%Nchain))

        return
    end subroutine allocate

    subroutine deallocate(self)
        implicit none
        type (PolymerModel), intent(inout) :: self

        if (allocated(self%config)) deallocate(self%config)
        if (allocated(self%bittable)) deallocate(self%bittable)
        if (allocated(self%dr)) deallocate(self%dr)
        if (allocated(self%contact)) deallocate(self%contact)
        if (allocated(self%boundary)) deallocate(self%boundary)
        if (allocated(self%loading_sites_factor)) deallocate(self%loading_sites_factor)
        if (allocated(self%interaction_sites_state)) deallocate(self%interaction_sites_state)
        if (allocated(self%rng_state)) deallocate(self%rng_state)
        if (allocated(self%rng_state_rep)) deallocate(self%rng_state_rep)
    end subroutine deallocate

    subroutine initbitable(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        type(Timer) :: crono
        integer :: i, j, k, a, v, ip, jp, kp, L2
        real :: x, y, z, xp, yp, zp

        ! include 'global.f03'

        ! initialize bittable (periodic BC): bittable(1,a)=monomers on node a, bittable(2:13,a)=adresses of 12 NN
        self%bittable = 0
        L2 = self%L ** 2
        call log%debug('initialize the bittable (periodic boundary conditions) for L: ' // trim(str(self%L)))
        call crono%Tic()
        do a = 1, 4 * L2 * self%L
            k = int((a - 1) / (2 * L2)) + 1
            j = int((a - 1) / self%L - 2 * self%L * (k - 1)) + 1
            i = a - self%L * (2 * self%L * (k - 1) + (j - 1))
            x = (i - 1) + 0.5 * (1 - mod(j + mod(k + 1, 2), 2))
            y = (j - 1) * 0.5
            z = (k - 1) * 0.5
            self%bittable(1, a) = 0
            do v = 1, 12
                xp = x + voisxyz(1, v + 1)
                yp = y + voisxyz(2, v + 1)
                zp = z + voisxyz(3, v + 1)
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
        end do
        call log%info('lattice density: ' // trim(strf(real(self%Nchain) / real(4 * L2 * self%L))))
        call log%info(crono%Tac(info = 'initbitable(): '))

        return
    end subroutine initbitable

    subroutine initconfig4(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        ! include 'global.f03'
        integer :: turn1(7), turn2(7), turn(7), lim, a, n, i, j, t, nv1, nv2, iv, en2, v, b, c(2, self%Nchain)

        !generate initial configuration without knots (do not change)

        self%bittable(1, :) = 0

        turn1 = (/13, 13, 2, 2, 12, 12, 3/)
        turn2 = (/13, 2, 2, 12, 12, 3, 3/)

        lim = self%L / 2

        self%config = 0
        a = int(4 * self%L**3 * self%get_random()) + 1
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
            t = int((n - 1) * self%get_random()) + 1
            iv = int((voisnn(1, 1, self%config(2, t)) - 1) * self%get_random()) + 1
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

        !!! NEW: Interaction state
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
        !!!!!!!!!!!!!!!!!!!!!!!!!

        return

    end subroutine initconfig4

    subroutine initconfig4_zigzag(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer :: a, n, i, t, nv1, nv2, iv, en2, v, b, c(2, self%Nchain), lim, v1, v2

        lim = 2 * self%L - 2

        i = int(6 * self%get_random()) + 1
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
            t = int((n - 1) * self%get_random()) + 1
            iv = int((voisnn(1, 1, self%config(2, t)) - 1) * self%get_random()) + 1
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

        !!! NEW: Interaction state
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
        !!!!!!!!!!!!!!!!!!!!!!!!!

        return

    end subroutine initconfig4_zigzag

    subroutine initconfig_sim_out(self, sim_out_dir, trajectory_i)
        implicit none
        class (PolymerModel), intent(inout) :: self
        character(*), intent(in) :: sim_out_dir
        integer, intent(in) :: trajectory_i
        character(1000) :: config_3dpoys_le_file
        character(1000) :: config_out_file
        integer :: initNchain, initNiter, initNmeas, initNinter, skip_monomers
        type(MonomerConfig) :: monomer_config
        integer :: rc
        integer :: a, n, v, b, lim, v1

        config_3dpoys_le_file = trim(sim_out_dir) // '/3dpoys_le.cfg'
        config_out_file = trim(sim_out_dir) // '/config.out'

        call load_config_file(config_3dpoys_le_file, iostat=rc)
        if (rc /= 0) then
            call log%error('Could not find or open input configuration file: ' // trim(config_out_file))
            call print_help()
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
            ! jump to the last measure, shift by trajectory index
            skip_monomers =  ((mod(trajectory_i-1, initNiter)+1) * (initNmeas-1) * initNchain)
            call log%info('Load config_out_file file ' // trim(config_out_file) &
                    // ' and for trajectory ' // trim(str(trajectory_i)) // ' skip ' // trim(str(skip_monomers)))

            do n=1, skip_monomers
                read(30, *, iostat = rc) monomer_config
                if (rc /= 0) exit
            end do

            do n=1, self%Nchain
                read(30, *, iostat = rc) monomer_config
                if (rc /= 0) exit
                a = monomer_config%laticeAddr
                v1 = monomer_config%nextTurn
                self%config(1, n) = a
                self%config(2, n) = v1
                self%bittable(1, a) = self%bittable(1, a) + 1
            end do

            !!! NEW: Interaction state
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
            !!!!!!!!!!!!!!!!!!!!!!!!!
        end if
        close(30)

        return

    end subroutine initconfig_sim_out

    subroutine trialmoveex(self)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        ! trial move to move a LEF leg
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, iv, s, id, con(3), strand
        real*8 :: fc
        real :: impermeability

        !choose randomly a monomer
        n = int(self%Nchain * self%get_random()) + 1

        s = self%contact(3, n)

        if ((s.eq.0).or.(n.eq.1).or.(n.eq.self%Nchain)) return !if monomer not occupied or end-monomer do nothing

        strand = (s + 1) / 2 + 1 ! 1 (-), 2 (+)
        impermeability = abs(self%boundary(strand, n))
        ! if impermeability is very low only the directionality is taken into account: the sign o
        fc = (1. - impermeability)**(1. / real(self%ikm)) ! probability of permeability
        ! check the singn of the boundary : -1 to decide to move or not

        !if (impermeability /= 0) then
        !    print*, 'fc[n:', n,'] = ', fc, ', impermeability=', impermeability
        !end if
        do iv = 1, self%ikm
            if (self%get_random().ge.(self%km * fc)) return
        end do

        if (s.eq.-1) then !if a (-) direction leg move to n-1
            if (connec(1, self%config(2, n - 1), self%contact(2, n)).eq.0) return !if break the slip-link do not move

            if (self%contact(1, n - 1).ne.0) then !if n-1 already occupied, try to swap
                if (.not. self%z_loop) return
                !if n-1 already occupied & z-loop, try to swap
                if ((self%contact(3, n - 1).eq.-1).or.(n.eq.2)) return ! same dir or end: no swap
                ! if break slip-link of n-1 do not swap
                if (connec(1, opp(self%config(2, n - 1)), self%contact(2, n - 1)).eq.0) &
                    return
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
        elseif (s.eq.1) then ! TODO Ask Daniel if this is only unidirectional!!! !if a (+1) direction leg move to n+1
            if (connec(1, opp(self%config(2, n)), self%contact(2, n)).eq.0) return !if break the slip-link do not move

            if (self%contact(1, n + 1).ne.0) then !if n+1 already occupied, try to swap
                if (.not. self%z_loop) return
                !if n+1 already occupied & z-loop, try to swap
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
        !elseif (.not. self%unidirectional) then
        !    call log%error('missied case s: ' // trim(str(s)) // ' in trialmoveex for monomer: ' // trim(str(n)))
        end if

        return

    end subroutine trialmoveex

    ! *************

    subroutine trialmovetad(self)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        !trial move for monomer
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, iv, v, b, j, nv1, nv2, nm2, np1, en, cn2, cn3, cm2, en2, id, cc, a
        real :: dE

        !choose randomly a monomer
        n = int(self%Nchain * self%get_random()) + 1
        en = self%config(1, n)

        !display interaction states
        !do iv = 1, self%Nchain
        !    if (self%interaction_sites_state(iv) > 0) then
        !        call log%info('interaction_site['// trim(str(iv)) // ']='// trim(str(self%interaction_sites_state(iv))))
        !    end if
        !end do
        !if (self%interaction_sites_state(n) > 0) then
        !    call log%info('interaction_site['// trim(str(n)) // ']='// trim(str(self%interaction_sites_state(n))))
        !end if

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
            iv = int(11 * self%get_random()) + 1
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

                !!! NEW: Interaction state
                if (self%interaction_sites_state(n).gt.0) then
                    dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                    !call log%info('NEW: Interaction state[' // trim(str(n)) // ']: Ei: ' // trim(strf(real(self%Ei))) &
                    !        // ', dE: ' // trim(strf(dE)))
                !else
                !    call log%info('NO: Interaction state[' // trim(str(n)) // ']')
                end if
                !!!!!!!!!!!!

                if ((self%contact(3, n).eq.-1).and.(id.lt.self%Nchain)) then
                    if (connec(opp(self%contact(2, n)), self%config(2, id), 1).ne.0) dE = dE - self%Ea
                    if (connec(opp(cc), self%config(2, id), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%contact(3, n + 1).eq.-1) then
                    if (connec(opp(cn2), self%contact(2, n + 1), 1).ne.0) dE = dE - self%Ea
                    if (connec(iv, self%contact(2, n + 1), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%get_random().lt.exp(-dE)) then
                    self%bittable(1, en) = self%bittable(1, en) - 1
                    self%bittable(1, v) = self%bittable(1, v) + 1

                    !!! NEW: Interaction state
                    self%bittable(14, en) = self%bittable(14, en) - 1
                    self%bittable(14, v) = self%bittable(14, v) + 1
                    do j=2,13
                        a=self%bittable(j, en)
                        self%bittable(14, a) = self%bittable(14, a) - 1
                        a=self%bittable(j, v)
                        self%bittable(14, a) = self%bittable(14, a) + 1
                    end do
                    !!!!!!!!!!!!

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

            iv = int(11 * self%get_random()) + 1
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
                ! TODO ckeck with Daniel if this is a propper bugfix
                if (cn3.eq.0) then
                    cc = 0
                else
                    cc = connec(iv, cn2, cn3)
                end if
                if ((id.ne.0).and.(cc.eq.0)) return
                dE = costhet(self%config(2, self%Nchain - 2), iv) - costhet(self%config(2, self%Nchain - 2), cn2)

                !!! NEW: Interaction state
                if (self%interaction_sites_state(n).gt.0) then
                    dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                    !call log%info('NEW2: Interaction state[' // trim(str(n)) // ']: Ei: ' // trim(strf(real(self%Ei))) &
                    !        // ', dE: ' // trim(strf(dE)))
                !else
                !    call log%info('NO: Interaction state[' // trim(str(n)) // ']')
                end if
                !!!!!!!!!!!!

                if ((self%contact(3, n).eq.1).and.(id.gt.1)) then
                    if (connec(self%contact(2, n), opp(self%config(2, id - 1)), 1).ne.0) dE = dE - self%Ea
                    if (connec(cc, opp(self%config(2, id - 1)), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%contact(3, n - 1).eq.1) then
                    if (connec(cn2, self%contact(2, n - 1), 1).ne.0) dE = dE - self%Ea
                    if (connec(iv, self%contact(2, n - 1), 1).ne.0) dE = dE + self%Ea
                end if

                if (self%get_random().lt.exp(-dE)) then
                    self%bittable(1, en) = self%bittable(1, en) - 1
                    self%bittable(1, v) = self%bittable(1, v) + 1

                    !!! NEW: Interaction state
                    self%bittable(14, en) = self%bittable(14, en) - 1
                    self%bittable(14, v) = self%bittable(14, v) + 1
                    do j=2,13
                        a = self%bittable(j, en)
                        self%bittable(14, a) = self%bittable(14, a) - 1
                        a = self%bittable(j, v)
                        self%bittable(14, a) = self%bittable(14, a) + 1
                    end do
                    !!!!!!!!!!!!

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

            !print*, 'size voisnn, cm2, cn2, voisnn(1, 1, 1) ', shape(voisnn), cm2, cn2, voisnn(1, 1, 1)

            if (voisnn(1, cm2, cn2).gt.1) then
                iv = int((voisnn(1, cm2, cn2) - 1) * self%get_random()) + 1
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
                    !print*, 'id, n, nv1, cm2, contact(2, n): ', id, n, nv1, cm2, cn3
                    !print*, 'shape(connec): ', shape(connec)
                    ! why self%contact(2, n)=0 and is a problem: only if you add -fbounds-check
                    if (cn3.eq.0) then
                        cc = 0
                    else
                        cc = connec(nv1, cm2, cn3)
                    end if
                    !print *, 'cc: ', cc
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

                    !!! NEW: Interaction state
                    if (self%interaction_sites_state(n).gt.0) then
                        dE = dE + self%Ei*(self%bittable(14, v) - self%bittable(14, en))
                        !call log%info('NEW3: Interaction state[' // trim(str(n)) // ']Ei:' // trim(strf(real(self%Ei))) &
                        !        // ', dE: ' // trim(strf(dE)))
                    end if
                    !!!!!!!!!!!!

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

                    if (self%get_random().lt.exp(-dE)) then
                        self%bittable(1, en) = self%bittable(1, en) - 1
                        self%bittable(1, v) = self%bittable(1, v) + 1

                        !!! NEW: Interaction state
                        self%bittable(14, en) = self%bittable(14, en) - 1
                        self%bittable(14, v) = self%bittable(14, v) + 1
                        do j=2,13
                            a = self%bittable(j, en)
                            self%bittable(14, a) = self%bittable(14, a) - 1
                            a = self%bittable(j, v)
                            self%bittable(14, a) = self%bittable(14, a) + 1
                        end do
                        !!!!!!!!!!!!

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

    end subroutine

    ! ****************

    subroutine trialbound(self)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        !trial move for binding LEF
        implicit none
        class (PolymerModel), intent(inout) :: self

        integer :: n, id, j, d
        real :: kbp

        !choose randomly a monomer
        n = int(self%Nchain * self%get_random()) + 1

        id = self%contact(1, n)
        if (id==0) then !if the bin is not occupied by a leg, try to randomly insert a LEF to NN sites
            kbp = self%kb * (  self%loading_sites_factor(n)**(1./real(self%ikb)) )
            do j = 1, self%ikb
                if (self%get_random()>=kbp) return
            end do
            d = int(2 * self%get_random()) + 1
            if (d==1) then
                if (n>1) then
                    if (self%contact(1, n - 1)==0) then
                        self%contact(1, n) = n - 1
                        self%contact(2, n) = opp(self%config(2, n - 1))
                        self%contact(1, n - 1) = n
                        self%contact(2, n - 1) = self%config(2, n - 1)
                        if (self%unidirectional) then
                            if (self%get_random()<0.5) then
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
                        if (self%get_random()<0.5) then
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

    end subroutine trialbound

    ! ******

    subroutine trialunbound(self)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        !trial move for unbinding LEF
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, id, j

        !choose randomly a monomer
        n = int(self%Nchain * self%get_random()) + 1

        id = self%contact(1, n)
        if (id.gt.0) then !if the monomer is occupied by a leg, try to remove it
            do j = 1, self%iku
                if (self%get_random().ge.self%ku) return
            end do
            self%contact(:, n) = 0
            self%contact(:, id) = 0
            write(13, *) n, id
            self%Nleffree = self%Nleffree + 1
        end if

        return

    end subroutine

    subroutine unbound_all(self)
        !unbinding of all bound LEFs
        implicit none
        class (PolymerModel), intent(inout) :: self

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

    end subroutine

    subroutine erase(self)
        implicit none
        class (PolymerModel), intent(inout) :: self

        integer :: i, a

        !erase the configuration and occupancy of the bittable

        do i = 1, self%Nchain
            a = self%config(1, i)
            self%bittable(1, a) = 0
        end do
        self%config = 0

        return
    end subroutine

    subroutine output(self)
        implicit none
        class (PolymerModel), intent(inout) :: self

        !write in output files

        integer :: j

        do j = 1, self%Nchain
            ! TODO maybe possible to write thw whole chain in one round
            write(10, *) self%config(:, j) ! config(1,j)=node of monomer j, config(2,j)=dir j->j+1
            write(11, *) self%dr(:, j) ! displacement (x,y,z) of monomer j from t=0 to now
            write(12, *) self%contact(:, j) ! contact(1,j)/=0: LEG at j; (1,j)=other LEG; (2,j)=vec; (3,j)=+/-1 dir
        end do
        write(14, *) self%Nleffree !number of free (unbound) extruders

        flush(10)
        flush(11)
        flush(12)
        flush(14)

        return
    end subroutine

    subroutine output_parameters(self, fout, init_mode, interaction_sites, boundary_file, lef_loading_sites, &
            basal_loading_factor, boundary_direction, &
            Niter, Ninter, Nmeas, burnin, burnout, burnoutM, radius_contact, &
            kb_a, ku_a, km_a)
        implicit none
        class (PolymerModel), intent(inout) :: self
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
        real, intent(in) :: kb_a, ku_a, km_a ! original argument's values for logging

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
        ! write(fout, '(a)') '# _k? = k?**(1. / real(ik?)) > 0.001' ! if needed to log for debug purpose
        ! write(fout, '(a)') '_kb=' // trim(strf(self%kb))
        ! write(fout, '(a)') '_ku=' // trim(strf(self%ku))
        ! write(fout, '(a)') '_km=' // trim(strf(self%km))
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
    end subroutine

    subroutine export_state_flat(self, config, bittable, contact, dr, boundary, &
            interaction_sites_state, loading_sites_factor, rng_state, Nleffree)
        implicit none
        class (PolymerModel), intent(in) :: self
        integer, intent(out), dimension(2, *) :: config
        integer, intent(out), dimension(14, *) :: bittable
        integer, intent(out), dimension(3, *) :: contact
        real, intent(out), dimension(3, *) :: dr
        real, intent(out), dimension(2, *) :: boundary
        integer, intent(out), dimension(*) :: interaction_sites_state
        real, intent(out), dimension(*) :: loading_sites_factor
        integer(8), intent(out), dimension(*) :: rng_state
        integer, intent(out) :: Nleffree
        integer :: nch, nbt
        nch = self%Nchain
        nbt = size(self%bittable, 2)
        config(1:2, 1:nch) = self%config
        bittable(1:14, 1:nbt) = self%bittable
        contact(1:3, 1:nch) = self%contact
        dr(1:3, 1:nch) = self%dr
        boundary(1:2, 1:nch) = self%boundary
        interaction_sites_state(1:nch) = self%interaction_sites_state
        loading_sites_factor(1:nch) = self%loading_sites_factor
        Nleffree = self%Nleffree
        if (allocated(self%rng_state_rep)) rng_state(1) = self%rng_state_rep(1)
    end subroutine export_state_flat

    subroutine import_state_flat(self, config, bittable, contact, dr, boundary, &
            interaction_sites_state, loading_sites_factor, rng_state, Nleffree)
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer, intent(in), dimension(2, *) :: config
        integer, intent(in), dimension(14, *) :: bittable
        integer, intent(in), dimension(3, *) :: contact
        real, intent(in), dimension(3, *) :: dr
        real, intent(in), dimension(2, *) :: boundary
        integer, intent(in), dimension(*) :: interaction_sites_state
        real, intent(in), dimension(*) :: loading_sites_factor
        integer(8), intent(in), dimension(*) :: rng_state
        integer, intent(in) :: Nleffree
        integer :: nch, nbt
        nch = self%Nchain
        nbt = size(self%bittable, 2)
        self%config = config(1:2, 1:nch)
        self%bittable = bittable(1:14, 1:nbt)
        self%contact = contact(1:3, 1:nch)
        self%dr = dr(1:3, 1:nch)
        self%boundary = boundary(1:2, 1:nch)
        self%interaction_sites_state = interaction_sites_state(1:nch)
        self%loading_sites_factor = loading_sites_factor(1:nch)
        self%Nleffree = Nleffree
        if (allocated(self%rng_state_rep)) self%rng_state_rep(1) = rng_state(1)
    end subroutine import_state_flat

    real(8) function get_random(self)
        ! Note: Removed !$acc routine seq because random_seed() intrinsic
        ! is not supported in NVHPC OpenACC device code
        class (PolymerModel), intent(inout) :: self
        if (allocated(self%rng_state_rep)) then
            get_random = repro_rng_next(self%rng_state_rep(1))
        else
            if (allocated(self%rng_state)) call random_seed(put=self%rng_state)
            get_random = randomnumber()
            if (allocated(self%rng_state)) call random_seed(get=self%rng_state)
        end if
        !get_random = randomnumber()
    end function get_random

    subroutine advance_one_measurement_interval(self, Ninter)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        ! Advance state by one measurement interval (Ninter * Ntrial steps). No I/O.
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer, intent(in) :: Ninter
        integer :: k, v, Ntrial
        real :: pt
        real*8 :: r

        Ntrial = 3 * self%Nchain + self%Nleffree
        pt = real(self%Nchain) / real(Ntrial)
        do k = 1, Ninter
            do v = 1, Ntrial
                r = self%get_random()
                if (r.lt.pt) then
                    call self%trialmovetad()
                elseif (r.lt.2 * pt) then
                    call self%trialmoveex()
                elseif (r.lt.3 * pt) then
                    call self%trialunbound()
                else
                    call self%trialbound()
                end if
            end do
        end do
        return
    end subroutine advance_one_measurement_interval

    subroutine run_burnin_phase(self, burnin)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        ! Run burn-in: simburnin = burnin + int(r*burnin), then simburnin*Nchain trialmovetad. No I/O.
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer, intent(in) :: burnin
        integer :: j, v, simburnin

        simburnin = burnin + int(self%get_random() * burnin)
        do j = 1, simburnin
            do v = 1, self%Nchain
                call self%trialmovetad()
            end do
        end do
        return
    end subroutine run_burnin_phase

    subroutine run_burnout_only_phase(self, burnout)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        ! Run burnout block: unbound_all + burnout*Nchain trialmovetad. No I/O.
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer, intent(in) :: burnout
        integer :: j, v

        if (burnout <= 0) return
        call self%unbound_all()
        do j = 1, burnout
            do v = 1, self%Nchain
                call self%trialmovetad()
            end do
        end do
        return
    end subroutine run_burnout_only_phase

    subroutine run_burnoutM_one_phase(self, Ninter)
        ! Note: Removed !$acc routine seq - not used in GPU flat array path
        ! Run one burnoutM step: unbound_all + Ninter*Nchain trialmovetad. No I/O.
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer, intent(in) :: Ninter
        integer :: k, v

        call self%unbound_all()
        do k = 1, Ninter
            do v = 1, self%Nchain
                call self%trialmovetad()
            end do
        end do
        return
    end subroutine run_burnoutM_one_phase

    subroutine do_simulation(self, trajectory_i, Ninter, Nmeas, burnin, burnout, burnoutM)
        implicit none
        class (PolymerModel), intent(inout) :: self
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

        burn_Nmeas = 0
        simburnin = 0

        ! now start the simulation with LEF: measurement after the burn-in time
        call log%info('Start Trajectory: ' // trim(str(trajectory_i)) // ', Ninter: ' // trim(str(Ninter)) // &
                ', Nmeas: ' // trim(str(Nmeas)) // ', burnin: ' // trim(str(burnin)) // &
                ', burnout: ' // trim(str(burnout)) // ', burnoutM: ' // trim(str(burnoutM)))

        ! first initial measurement before the burn-in
        call self%output()

        if (burnin > 0) then
            burn_Nmeas = burn_Nmeas + 1

            ! choose rondome simulation burn-in time = initial_burnin + random_fraction * initial_burnin
            ! add random time (a fraction) after the initial burn-in time
            r = self%get_random()
            simburnin = burnin + int(r * burnin)

            do j = 1, simburnin
                do v = 1, self%Nchain
                    call self%trialmovetad() !trial move for monomers
                end do
            end do
            ! now start the simulation with LEF: measurement after the burn-in time
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-in Measurement, ' // &
                    ' burn_Nmeas: ' // trim(str(burn_Nmeas)))
            call self%output()
        end if

        if (burnout > 0) then
            burn_Nmeas = burn_Nmeas + 1
        end if

        if (burnoutM > 0) then
            burn_Nmeas = burn_Nmeas + burnoutM
        end if

        do j = 1, Nmeas - 1 - burn_Nmeas  ! account for the initial, burn-in and burn-out measurements
            call self%advance_one_measurement_interval(Ninter)
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', Measurement:' // trim(str(j)))
            call flush(6)
            call self%output()
            call flush(13)
        end do

        ! do burn-out included as extra measurement
        if (burnout > 0) then
            call self%unbound_all()
            do j = 1, burnout
                do v = 1, self%Nchain
                    call self%trialmovetad() ! move for monomers
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', fix burn-out Measurement')
            call self%output()
            call flush(13)
        end if

        if (burnoutM > 0) then
            call self%unbound_all()
            do j = 1, burnoutM
                do k = 1, Ninter
                    do v = 1, self%Nchain
                        call self%trialmovetad() ! move for monomers
                    end do
                end do
                call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', burn-out Measurement: ' // trim(str(j)))
                call self%output()
                call flush(13)
            end do
        end if

        call self%erase()
    end subroutine do_simulation

end module PolymerModel_mod