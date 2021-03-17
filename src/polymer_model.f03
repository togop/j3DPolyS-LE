! Created by  on 19.09.19.

module PolymerModel_mod
    use Timers
    use lattice_data_mod
    use logging_mod

    implicit none
    private
    public :: PolymerModel, ModelParameters
    type(Logger) :: log = Logger('PolymerModel_mod', LOG_INFO)

    type ModelParameters
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        real, public :: kb, ku, km, Ea
        ! integer, public :: Niter, Nmeas, Ninter, kint  ! not used
    end type ModelParameters

    type PolymerModel
        private
        integer, public :: L, Nchain, iku, ikm, ikb, Nleffree
        ! GLOBAL: integer ::voisnn(13,13,13),opp(13),connec(13,13,13)
        ! GLOBAL: real,dimension(3,13)	::voisxyz
        ! GLOBAL: real,dimension(13,13)	::costhet
        real, public :: kb, ku, km, Ea
        logical, public :: z_loop = .false.
        logical, public :: unidirectional = .false.
        ! allocatable
        integer, dimension(:, :), allocatable :: config
        integer, dimension(:, :), allocatable :: bittable
        real, dimension(:, :), allocatable :: dr
        integer, dimension(:, :), allocatable :: contact
        real, dimension(:, :), allocatable :: boundary

    contains
        procedure, public :: init, do_simulation, trialmoveex, trialmovetad, trialbound, trialunbound, unbound_all, &
                erase, output
        procedure :: allocate, initbitable, initconfig4, initconfig4_zigzag
        final :: deallocate
    end type PolymerModel

contains

    subroutine init(self, boundary, init_mode)
        implicit none
        class (PolymerModel), intent(inout) :: self
        real, dimension(:, :) :: boundary
        character(len = 1), intent(in) :: init_mode

        call log%set_level(global_log_level)

        call self%allocate()

        self%boundary = boundary
        !print*, 'initialize boundary ', shape(boundary), shape(self%boundary)

        call self%initbitable()
        call log%info('init_mode ' // init_mode)
        select case(init_mode)
        case('h')
            call self%initconfig4()
        case('z')
            call self%initconfig4_zigzag()
        case default   ! h: helices
            call log%warn('unknown init_mode ' // init_mode // ' , it will be used the default: helices !')
            call self%initconfig4()
        end select

        !make simulations
        self%contact = 0
        self%dr = 0.

        return
    end subroutine

    subroutine allocate(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer :: bittable_t

        ! include 'global.f03'

        bittable_t = 4 * (self%L**3)

        call log%debug('allocate self%Nchain: ' // trim(str(self%Nchain)) // ' bittable_t: ' // trim(str(bittable_t)))

        ! call deallocate(self) ! make sure it's free
        allocate (self%config(2, self%Nchain))
        allocate (self%bittable(13, bittable_t))
        allocate (self%dr(3, self%Nchain))
        allocate (self%contact(3, self%Nchain))
        allocate (self%boundary(2, self%Nchain))

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
    end subroutine deallocate

    subroutine initbitable(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        type(Timer) :: crono
        integer :: i, j, k, a, v, ip, jp, kp, L2
        real :: x, y, z, xp, yp, zp

        ! include 'global.f03'

        !initialize the bittable (periodic boundary conditions) bittable(1,a)=number of monomers on lattice node a, bittable(2:13,a)=adresses of the 12 NN
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
                if (xp.ge.self%L) xp = xp - self%L
                if (xp.lt.0) xp = xp + self%L
                if (yp.ge.self%L) yp = yp - self%L
                if (yp.lt.0) yp = yp + self%L
                if (zp.ge.self%L) zp = zp - self%L
                if (zp.lt.0) zp = zp + self%L
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
        real*8 :: randomnumber

        !generate initial configuration without knots (do not change)

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

        return

    end subroutine initconfig4

    subroutine initconfig4_zigzag(self)
        implicit none
        class (PolymerModel), intent(inout) :: self
        integer :: a, n, i, t, nv1, nv2, iv, en2, v, b, c(2, self%Nchain), lim, v1, v2
        real*8 :: randomnumber

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

        return

    end subroutine initconfig4_zigzag

    subroutine trialmoveex(self)
        ! trial move to move a LEF leg
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, iv, s, id, con(3), strand
        real*8 :: randomnumber, fc
        real :: impermeability

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1

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
            if (randomnumber().ge.(self%km * fc)) return
        end do

        if (s.eq.-1) then !if a (-) direction leg move to n-1
            if (connec(1, self%config(2, n - 1), self%contact(2, n)).eq.0) return !if break the slip-link do not move

            if (self%contact(1, n - 1).ne.0) then !if n-1 already occupied, try to swap
                if (.not. self%z_loop) return
                !if n-1 already occupied & z-loop, try to swap
                if ((self%contact(3, n - 1).eq.-1).or.(n.eq.2)) return !if same direction or at the end do not swap
                if (connec(1, opp(self%config(2, n - 1)), self%contact(2, n - 1)).eq.0) return !if break the slip-link of n-1 do not swap
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
        !trial move for monomer
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, iv, v, b, nv1, nv2, nm2, np1, en, cn2, cn3, cm2, en2, id, cc
        real :: dE
        real*8 :: randomnumber

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
                ! TODO ckeck with Daniel if this is a propper bugfix
                if (cn3.eq.0) then
                    cc = 0
                else
                    cc = connec(iv, opp(cn2), cn3)
                end if
                if ((id.ne.0).and.(cc.eq.0)) return
                dE = costhet(opp(iv), self%config(2, 2)) - costhet(cn2, self%config(2, 2))

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
                ! TODO ckeck with Daniel if this is a propper bugfix
                if (cn3.eq.0) then
                    cc = 0
                else
                    cc = connec(iv, cn2, cn3)
                end if
                if ((id.ne.0).and.(cc.eq.0)) return
                dE = costhet(self%config(2, self%Nchain - 2), iv) - costhet(self%config(2, self%Nchain - 2), cn2)

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
                    !print*, 'id, n, nv1, cm2, contact(2, n): ', id, n, nv1, cm2, cn3
                    !print*, 'shape(connec): ', shape(connec)
                    ! why self%contact(2, n)=0 and is a problem: only if you add -fbounds-check
                    ! TODO ckeck with Daniel if this is a propper bugfix
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
        !trial move for binding LEF
        implicit none
        class (PolymerModel), intent(inout) :: self

        integer :: n, id, j, d
        real*8 :: randomnumber

        !choose randomly a monomer
        n = int(self%Nchain * randomnumber()) + 1

        id = self%contact(1, n)
        if (id==0) then !if the bin is not occupied by a leg, try to randomly insert a LEF to NN sites
            do j = 1, self%ikb
                if (randomnumber()>=self%kb) return
            end do
            d = int(2 * randomnumber()) + 1
            if (d==1) then
                ! bugfix IFs and ands  ! TODO maybe remind Daniel
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

    end subroutine trialbound

    ! ******

    subroutine trialunbound(self)
        !trial move for unbinding LEF
        implicit none
        class (PolymerModel), intent(inout) :: self

        ! include 'global.f03'

        integer :: n, id, j
        real*8 :: randomnumber

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
            write(10, *) self%config(:, j) !the configuration: config(1,j)=node where monomer j is located, config(2,j)=direction of the vector between j and j+1
            write(11, *) self%dr(:, j) !the displacement: vector (x,y,z in lattice unit) of displacement of monomer j between time 0 and current time
            write(12, *) self%contact(:, j) !the extruder occupancy and information: contact(1,j) \ne 0 if one lef of a extruder is in j, contact(1,j)=the monomer where the other lef is, contact(2,j)=the vector between the two legs, contact(3,j)=-1 (resp. +1) if it's a lef walking in the (-) direction (resp. (+) direction)
        end do
        write(14, *) self%Nleffree !number of free (unbound) extruders

        flush(10)
        flush(11)
        flush(12)
        flush(14)

        return
    end subroutine

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
        real*8 :: randomnumber, r

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

            ! choose rondome simulation burn-in timw = initial_burnin + random_fraction * initial_burnin
            ! add random time (a fraction) after the initial burn-in time
            r = randomnumber()
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
            do k = 1, Ninter
                Ntrial = 3 * self%Nchain + self%Nleffree
                pt = real(self%Nchain) / real(Ntrial)
                do v = 1, Ntrial
                    r = randomnumber()
                    if (r.lt.pt) then
                        call self%trialmovetad() !trial move for monomers
                    elseif (r.lt.2 * pt) then
                        call self%trialmoveex() !trial move for LEF movement
                    elseif (r.lt.3 * pt) then
                        call self%trialunbound() !trial move for unbinding event
                    else
                        call self%trialbound() !trial move for binding event
                    end if
                end do
            end do
            call log%info('Trajectory: ' // trim(str(trajectory_i)) // ', Measurement:' // trim(str(j)))
            call flush(6)
            call self%output()
            flush(13)
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
            flush(13)
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
                flush(13)
            end do
        end if

        call self%erase()
    end subroutine do_simulation

end module PolymerModel_mod