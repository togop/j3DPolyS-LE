! Flattened GPU replica path: state in plain arrays so OpenACC can actually offload.
! Use when use_gpu is true for reliable GPU execution (avoids array-of-derived-type deep copy).
! See sim_src/GPU_DESIGN.md section 10.

module gpu_flat_mod
    use PolymerModel_mod
    use lattice_data_mod
    use logging_mod
    use repro_rng_mod, only: repro_rng_next
    implicit none
    private
    public :: do_simulation_replicas_flat

    ! Flat arrays (allocated in do_simulation_replicas_flat, then copyin to device)
    integer, allocatable :: config_f(:, :, :), bittable_f(:, :, :), contact_f(:, :, :), &
            interaction_sites_state_f(:, :)
    integer(8), allocatable :: rng_state_f(:, :)   ! (1, nrep) portable RNG state
    real, allocatable :: dr_f(:, :, :), boundary_f(:, :, :), loading_sites_factor_f(:, :)
    integer, allocatable :: Nleffree_f(:)
    ! RNG buffering for reduced function call overhead
    real*8, allocatable :: rng_buffer_f(:, :)  ! (BUFFER_SIZE, nrep) pre-generated random numbers
    integer, allocatable :: rng_buffer_idx_f(:) ! (nrep) current index in buffer
    integer, parameter :: RNG_BUFFER_SIZE = 5000
    !$acc declare create(config_f, bittable_f, contact_f, dr_f, boundary_f, &
    !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f, &
    !$acc   rng_buffer_f, rng_buffer_idx_f)
    ! Scalars (same for all replicas; used in device routines)
    integer :: nrep_f, Nchain_f, bittable_t_f, L_f, iku_f, ikm_f, ikb_f
    integer :: burnin_f, burnout_f, Ninter_f
    real :: Ea_f, Ei_f, ku_f, km_f, kb_f
    logical :: unidirectional_f, z_loop_f
    !$acc declare copyin(nrep_f, Nchain_f, bittable_t_f, L_f, iku_f, ikm_f, ikb_f, &
    !$acc   burnin_f, burnout_f, Ninter_f, Ea_f, Ei_f, ku_f, km_f, kb_f, unidirectional_f, z_loop_f)

contains

    subroutine do_simulation_replicas_flat(replicas, nrep, Ninter, Nmeas, burnin, burnout, burnoutM, batch_size)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: nrep, Ninter, Nmeas, burnin, burnout, burnoutM
        integer, intent(in), optional :: batch_size
        integer :: burn_Nmeas, j, r, output_interval, batch, actual_batch_size
        type(Logger) :: logger

        ! Initialize logger
        logger%source = 'gpu_flat_mod'
        logger%level = LOG_INFO

        ! Default batch size: 10 measurement intervals before transfer (reduces transfers by 10x)
        output_interval = 10
        if (present(batch_size)) then
            if (batch_size > 0) output_interval = batch_size
        end if
        call logger%info('GPU batch size (measurements per transfer): ' // trim(str(output_interval)))
        burn_Nmeas = 0
        if (burnin > 0) burn_Nmeas = burn_Nmeas + 1
        if (burnout > 0) burn_Nmeas = burn_Nmeas + 1
        burn_Nmeas = burn_Nmeas + burnoutM

        nrep_f = nrep
        Nchain_f = replicas(1)%Nchain
        L_f = replicas(1)%L
        bittable_t_f = 4 * (L_f**3)
        iku_f = replicas(1)%iku
        ikm_f = replicas(1)%ikm
        ikb_f = replicas(1)%ikb
        Ea_f = replicas(1)%Ea
        Ei_f = replicas(1)%Ei
        ku_f = replicas(1)%ku
        km_f = replicas(1)%km
        kb_f = replicas(1)%kb
        unidirectional_f = replicas(1)%unidirectional
        z_loop_f = replicas(1)%z_loop
        !$acc update device(nrep_f, Nchain_f, bittable_t_f, L_f, iku_f, ikm_f, ikb_f, &
        !$acc   Ea_f, Ei_f, ku_f, km_f, kb_f, unidirectional_f, z_loop_f)

        ! Allocate flat arrays (rng_state_f: 1 integer(8) per replica for portable RNG)
        allocate(config_f(2, Nchain_f, nrep), bittable_f(14, bittable_t_f, nrep), &
                contact_f(3, Nchain_f, nrep), dr_f(3, Nchain_f, nrep), boundary_f(2, Nchain_f, nrep), &
                interaction_sites_state_f(Nchain_f, nrep), loading_sites_factor_f(Nchain_f, nrep), &
                rng_state_f(1, nrep), Nleffree_f(nrep), &
                rng_buffer_f(RNG_BUFFER_SIZE, nrep), rng_buffer_idx_f(nrep))

        ! Initial output (from replicas on host)
        do r = 1, nrep
            call replicas(r)%output()
        end do

        ! Pack host replicas -> flat
        call pack_replicas_to_flat(replicas, nrep)

        ! Initialize RNG buffers (will be filled on first use on device)
        rng_buffer_idx_f = RNG_BUFFER_SIZE + 1  ! Force initial refill

        !$acc enter data copyin(config_f, bittable_f, contact_f, dr_f, boundary_f, &
        !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f, &
        !$acc   rng_buffer_f, rng_buffer_idx_f)

        ! Burn-in
        if (burnin > 0) then
            burnin_f = burnin
            !$acc update device(burnin_f)
            !$acc parallel loop gang vector_length(128) &
            !$acc   present(config_f, bittable_f, contact_f, dr_f, boundary_f, &
            !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f) &
            !$acc   copyin(nrep_f, Nchain_f, bittable_t_f,  Ea_f, Ei_f, L_f, iku_f, ikm_f, ikb_f, &
            !$acc   ku_f, km_f, kb_f, unidirectional_f, z_loop_f, burnin_f) &
            !$acc   async(1)
            do r = 1, nrep_f
                call run_burnin_phase_flat(r)
            end do
            !$acc end parallel loop
            ! Optimized transfer: skip bittable_f as it's not needed for output
            !$acc update host(config_f, contact_f, dr_f, Nleffree_f, rng_state_f) async(1)
            !$acc wait(1)
            call unpack_flat_to_replicas(replicas, nrep)
            do r = 1, nrep
                call replicas(r)%output()
            end do
        end if

        ! Main measurement loop with batching, async operations, and optimized transfers
        Ninter_f = Ninter
        !$acc update device(Ninter_f) async(1)
        j = 1
        do while (j <= Nmeas - 1 - burn_Nmeas)
            ! Determine actual batch size for this iteration
            actual_batch_size = min(output_interval, Nmeas - burn_Nmeas - j)

            ! Run multiple measurement intervals on GPU without transfer (batch processing)
            do batch = 1, actual_batch_size
                !$acc parallel loop gang vector_length(128) &
                !$acc   present(config_f, bittable_f, contact_f, dr_f, boundary_f, &
                !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f) &
                !$acc   copyin(nrep_f, Nchain_f, bittable_t_f,  Ea_f, Ei_f, L_f, iku_f, ikm_f, ikb_f, &
                !$acc   ku_f, km_f, kb_f, unidirectional_f, z_loop_f, Ninter_f) &
                !$acc   async(1)
                do r = 1, nrep_f
                    call advance_one_measurement_interval_flat(r)
                end do
                !$acc end parallel loop
            end do

            ! Optimized transfer: only data needed for output (skip bittable_f - 70% size reduction)
            !$acc update host(config_f, contact_f, dr_f, Nleffree_f, rng_state_f) async(1)
            !$acc wait(1)

            call unpack_flat_to_replicas(replicas, nrep)
            call logger%info('Replica batch (flat GPU), Measurement: ' // trim(str(j + actual_batch_size - 1)) &
                // ' (batched ' // trim(str(actual_batch_size)) // ' intervals)')
            call flush(6)
            do r = 1, nrep
                call replicas(r)%output()
            end do
            call flush(13)

            j = j + actual_batch_size
        end do

        ! Burnout block
        if (burnout > 0) then
            burnout_f = burnout
            !$acc update device(burnout_f)
            !$acc parallel loop gang vector_length(128) &
            !$acc   present(config_f, bittable_f, contact_f, dr_f, boundary_f, &
            !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f) &
            !$acc   copyin(nrep_f, Nchain_f, bittable_t_f,  Ea_f, Ei_f, L_f, iku_f, ikm_f, ikb_f, &
            !$acc   ku_f, km_f, kb_f, unidirectional_f, z_loop_f, burnout_f) &
            !$acc   async(1)
            do r = 1, nrep_f
                call run_burnout_only_phase_flat(r)
            end do
            !$acc end parallel loop
            ! Optimized transfer: skip bittable_f as it's not needed for output
            !$acc update host(config_f, contact_f, dr_f, Nleffree_f, rng_state_f) async(1)
            !$acc wait(1)
            call unpack_flat_to_replicas(replicas, nrep)
            do r = 1, nrep
                call replicas(r)%output()
            end do
            call flush(13)
        end if

        ! BurnoutM steps
        if (burnoutM > 0) then
            Ninter_f = Ninter
            !$acc update device(Ninter_f) async(1)
            do j = 1, burnoutM
                !$acc parallel loop gang vector_length(128) &
                !$acc   present(config_f, bittable_f, contact_f, dr_f, boundary_f, &
                !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f) &
                !$acc   copyin(nrep_f, Nchain_f, bittable_t_f,  Ea_f, Ei_f, L_f, iku_f, ikm_f, ikb_f, &
                !$acc   ku_f, km_f, kb_f, unidirectional_f, z_loop_f, Ninter_f) &
                !$acc   async(1)
                do r = 1, nrep_f
                    call run_burnoutM_one_phase_flat(r)
                end do
                !$acc end parallel loop
                ! Optimized transfer: skip bittable_f as it's not needed for output
                !$acc update host(config_f, contact_f, dr_f, Nleffree_f, rng_state_f) async(1)
                !$acc wait(1)
                call unpack_flat_to_replicas(replicas, nrep)
                call logger%info('Replica batch (flat GPU), burn-out Measurement: ' // trim(str(j)))
                do r = 1, nrep
                    call replicas(r)%output()
                end do
                call flush(13)
            end do
        end if

        !$acc exit data copyout(config_f, bittable_f, contact_f, dr_f, boundary_f, &
        !$acc   interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f) &
        !$acc   delete(rng_buffer_f, rng_buffer_idx_f)

        deallocate(config_f, bittable_f, contact_f, dr_f, boundary_f, &
                interaction_sites_state_f, loading_sites_factor_f, rng_state_f, Nleffree_f, &
                rng_buffer_f, rng_buffer_idx_f)

        do r = 1, nrep
            call replicas(r)%erase()
        end do
    end subroutine do_simulation_replicas_flat

    subroutine pack_replicas_to_flat(replicas, nrep)
        type(PolymerModel), intent(in) :: replicas(:)
        integer, intent(in) :: nrep
        integer :: r
        do r = 1, nrep
            call replicas(r)%export_state_flat(config_f(:, :, r), bittable_f(:, :, r), contact_f(:, :, r), &
                    dr_f(:, :, r), boundary_f(:, :, r), interaction_sites_state_f(:, r), loading_sites_factor_f(:, r), &
                    rng_state_f(:, r), Nleffree_f(r))
        end do
    end subroutine pack_replicas_to_flat

    subroutine unpack_flat_to_replicas(replicas, nrep)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: nrep
        integer :: r
        do r = 1, nrep
            call replicas(r)%import_state_flat(config_f(:, :, r), bittable_f(:, :, r), contact_f(:, :, r), &
                    dr_f(:, :, r), boundary_f(:, :, r), interaction_sites_state_f(:, r), loading_sites_factor_f(:, r), &
                    rng_state_f(:, r), Nleffree_f(r))
        end do
    end subroutine unpack_flat_to_replicas

    ! RNG buffering helper: refill buffer for a replica
    subroutine refill_rng_buffer_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: i
        do i = 1, RNG_BUFFER_SIZE
            rng_buffer_f(i, irep) = repro_rng_next(rng_state_f(1, irep))
        end do
        rng_buffer_idx_f(irep) = 1
    end subroutine refill_rng_buffer_flat

    ! RNG buffering helper: get next random number (with auto-refill)
    function get_random_flat(irep) result(rv)
        !$acc routine seq
        integer, intent(in) :: irep
        real*8 :: rv
        if (rng_buffer_idx_f(irep) > RNG_BUFFER_SIZE) then
            call refill_rng_buffer_flat(irep)
        end if
        rv = rng_buffer_f(rng_buffer_idx_f(irep), irep)
        rng_buffer_idx_f(irep) = rng_buffer_idx_f(irep) + 1
    end function get_random_flat

    subroutine advance_one_measurement_interval_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: k, v, Ntrial
        real :: pt
        real*8 :: rv

        Ntrial = 3 * Nchain_f + Nleffree_f(irep)
        pt = real(Nchain_f) / real(Ntrial)
        do k = 1, Ninter_f
            do v = 1, Ntrial
                rv = repro_rng_next(rng_state_f(1, irep))
                if (rv < pt) then
                    call trialmovetad_flat(irep)
                else if (rv < 2 * pt) then
                    call trialmoveex_flat(irep)
                else if (rv < 3 * pt) then
                    call trialunbound_flat(irep)
                else
                    call trialbound_flat(irep)
                end if
            end do
        end do
    end subroutine advance_one_measurement_interval_flat

    subroutine run_burnin_phase_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: j, v
        real*8 :: rv
        integer :: simburnin
        rv = repro_rng_next(rng_state_f(1, irep))
        simburnin = burnin_f + int(rv * burnin_f)
        do j = 1, simburnin
            do v = 1, Nchain_f
                call trialmovetad_flat(irep)
            end do
        end do
    end subroutine run_burnin_phase_flat

    subroutine run_burnout_only_phase_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: j, v
        if (burnout_f <= 0) return
        call unbound_all_flat(irep)
        do j = 1, burnout_f
            do v = 1, Nchain_f
                call trialmovetad_flat(irep)
            end do
        end do
    end subroutine run_burnout_only_phase_flat

    subroutine run_burnoutM_one_phase_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: k, v
        call unbound_all_flat(irep)
        do k = 1, Ninter_f
            do v = 1, Nchain_f
                call trialmovetad_flat(irep)
            end do
        end do
    end subroutine run_burnoutM_one_phase_flat

    subroutine unbound_all_flat(irep)
        !$acc routine seq
        ! On device we skip write(13) (no I/O). Unbind log not written for flat GPU path.
        integer, intent(in) :: irep
        integer :: n, id
        do n = 1, Nchain_f
            id = contact_f(1, n, irep)
            if (id > 0) then
                contact_f(:, n, irep) = 0
                contact_f(:, id, irep) = 0
                Nleffree_f(irep) = Nleffree_f(irep) + 1
            end if
        end do
    end subroutine unbound_all_flat

    subroutine trialmovetad_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: n, iv, v, b, j, nv1, nv2, nm2, np1, en, cn2, cn3, cm2, en2, id, cc, a
        real :: dE
        real*8 :: rv
        rv = repro_rng_next(rng_state_f(1, irep))
        n = int(Nchain_f * rv) + 1
        en = config_f(1, n, irep)
        if (n == 1) then
            en2 = config_f(1, 2, irep)
            cn2 = opp(config_f(2, 1, irep))
            if (cn2 < config_f(2, 2, irep)) then
                cm2 = config_f(2, 2, irep)
            else
                cm2 = cn2
                cn2 = config_f(2, 2, irep)
            end if
            iv = int(11 * repro_rng_next(rng_state_f(1, irep))) + 1
            if (iv >= cn2) iv = iv + 1
            if (iv >= cm2) iv = iv + 1
            if (iv == 1) then
                v = en2
            else
                v = bittable_f(iv, en2, irep)
            end if
            b = bittable_f(1, v, irep)
            if ((b == 0) .or. ((b == 1) .and. (en2 == v))) then
                id = contact_f(1, n, irep)
                cn2 = config_f(2, 1, irep)
                cn3 = contact_f(2, n, irep)
                if (cn3 == 0) then
                    cc = 0
                else
                    cc = connec(iv, opp(cn2), cn3)
                end if
                if ((id /= 0) .and. (cc == 0)) return
                dE = costhet(opp(iv), config_f(2, 2, irep)) - costhet(cn2, config_f(2, 2, irep))
                if (interaction_sites_state_f(n, irep) > 0) &
                    dE = dE + Ei_f * (bittable_f(14, v, irep) - bittable_f(14, en, irep))
                if ((contact_f(3, n, irep) == -1) .and. (id < Nchain_f)) then
                    if (connec(opp(contact_f(2, n, irep)), config_f(2, id, irep), 1) /= 0) dE = dE - Ea_f
                    if (connec(opp(cc), config_f(2, id, irep), 1) /= 0) dE = dE + Ea_f
                end if
                if (contact_f(3, n + 1, irep) == -1) then
                    if (connec(opp(cn2), contact_f(2, n + 1, irep), 1) /= 0) dE = dE - Ea_f
                    if (connec(iv, contact_f(2, n + 1, irep), 1) /= 0) dE = dE + Ea_f
                end if
                if (repro_rng_next(rng_state_f(1, irep)) < exp(-dE)) then
                    bittable_f(1, en, irep) = bittable_f(1, en, irep) - 1
                    bittable_f(1, v, irep) = bittable_f(1, v, irep) + 1
                    bittable_f(14, en, irep) = bittable_f(14, en, irep) - 1
                    bittable_f(14, v, irep) = bittable_f(14, v, irep) + 1
                    do j = 2, 13
                        a = bittable_f(j, en, irep)
                        bittable_f(14, a, irep) = bittable_f(14, a, irep) - 1
                        a = bittable_f(j, v, irep)
                        bittable_f(14, a, irep) = bittable_f(14, a, irep) + 1
                    end do
                    config_f(1, 1, irep) = v
                    config_f(2, 1, irep) = opp(iv)
                    if (id /= 0) then
                        contact_f(2, n, irep) = cc
                        contact_f(2, id, irep) = opp(cc)
                    end if
                    dr_f(1, n, irep) = dr_f(1, n, irep) + voisxyz(1, iv) + voisxyz(1, cn2)
                    dr_f(2, n, irep) = dr_f(2, n, irep) + voisxyz(2, iv) + voisxyz(2, cn2)
                    dr_f(3, n, irep) = dr_f(3, n, irep) + voisxyz(3, iv) + voisxyz(3, cn2)
                end if
            end if
        else if (n == Nchain_f) then
            en2 = config_f(1, Nchain_f - 1, irep)
            cn2 = config_f(2, Nchain_f - 1, irep)
            if (cn2 < opp(config_f(2, Nchain_f - 2, irep))) then
                cm2 = opp(config_f(2, Nchain_f - 2, irep))
            else
                cm2 = cn2
                cn2 = opp(config_f(2, Nchain_f - 2, irep))
            end if
            iv = int(11 * repro_rng_next(rng_state_f(1, irep))) + 1
            if (iv >= cn2) iv = iv + 1
            if (iv >= cm2) iv = iv + 1
            if (iv == 1) then
                v = en2
            else
                v = bittable_f(iv, en2, irep)
            end if
            b = bittable_f(1, v, irep)
            if ((b == 0) .or. ((b == 1) .and. (en2 == v))) then
                id = contact_f(1, n, irep)
                cn2 = config_f(2, Nchain_f - 1, irep)
                cn3 = contact_f(2, n, irep)
                if (cn3 == 0) then
                    cc = 0
                else
                    cc = connec(iv, cn2, cn3)
                end if
                if ((id /= 0) .and. (cc == 0)) return
                dE = costhet(config_f(2, Nchain_f - 2, irep), iv) - costhet(config_f(2, Nchain_f - 2, irep), cn2)
                if (interaction_sites_state_f(n, irep) > 0) &
                    dE = dE + Ei_f * (bittable_f(14, v, irep) - bittable_f(14, en, irep))
                if ((contact_f(3, n, irep) == 1) .and. (id > 1)) then
                    if (connec(contact_f(2, n, irep), opp(config_f(2, id - 1, irep)), 1) /= 0) dE = dE - Ea_f
                    if (connec(cc, opp(config_f(2, id - 1, irep)), 1) /= 0) dE = dE + Ea_f
                end if
                if (contact_f(3, n - 1, irep) == 1) then
                    if (connec(cn2, contact_f(2, n - 1, irep), 1) /= 0) dE = dE - Ea_f
                    if (connec(iv, contact_f(2, n - 1, irep), 1) /= 0) dE = dE + Ea_f
                end if
                if (repro_rng_next(rng_state_f(1, irep)) < exp(-dE)) then
                    bittable_f(1, en, irep) = bittable_f(1, en, irep) - 1
                    bittable_f(1, v, irep) = bittable_f(1, v, irep) + 1
                    bittable_f(14, en, irep) = bittable_f(14, en, irep) - 1
                    bittable_f(14, v, irep) = bittable_f(14, v, irep) + 1
                    do j = 2, 13
                        a = bittable_f(j, en, irep)
                        bittable_f(14, a, irep) = bittable_f(14, a, irep) - 1
                        a = bittable_f(j, v, irep)
                        bittable_f(14, a, irep) = bittable_f(14, a, irep) + 1
                    end do
                    config_f(1, Nchain_f, irep) = v
                    config_f(2, Nchain_f - 1, irep) = iv
                    if (id /= 0) then
                        contact_f(2, n, irep) = cc
                        contact_f(2, id, irep) = opp(cc)
                    end if
                    dr_f(1, n, irep) = dr_f(1, n, irep) + voisxyz(1, iv) - voisxyz(1, cn2)
                    dr_f(2, n, irep) = dr_f(2, n, irep) + voisxyz(2, iv) - voisxyz(2, cn2)
                    dr_f(3, n, irep) = dr_f(3, n, irep) + voisxyz(3, iv) - voisxyz(3, cn2)
                end if
            end if
        else
            cn2 = config_f(2, n, irep)
            cm2 = config_f(2, n - 1, irep)
            en2 = config_f(1, n - 1, irep)
            nm2 = n - 2
            np1 = n + 1
            if (voisnn(1, cm2, cn2) > 1) then
                iv = int((voisnn(1, cm2, cn2) - 1) * repro_rng_next(rng_state_f(1, irep))) + 1
                if (voisnn(2 * iv, cm2, cn2) >= cm2) iv = iv + 1
                nv1 = voisnn(2 * iv, cm2, cn2)
                nv2 = voisnn(2 * iv + 1, cm2, cn2)
                if (nv1 == 1) then
                    v = en2
                else
                    v = bittable_f(nv1, en2, irep)
                end if
                b = bittable_f(1, v, irep)
                if ((b == 0) .or. ((b == 1) .and. ((v == en2) .or. (v == config_f(1, np1, irep))))) then
                    id = contact_f(1, n, irep)
                    cn3 = contact_f(2, n, irep)
                    if (cn3 == 0) then
                        cc = 0
                    else
                        cc = connec(nv1, cm2, cn3)
                    end if
                    if ((id /= 0) .and. (cc == 0)) return
                    if (n == 2) then
                        dE = costhet(nv1, nv2) + costhet(nv2, config_f(2, np1, irep)) - costhet(cm2, cn2) &
                                - costhet(cn2, config_f(2, np1, irep))
                    else if (n == Nchain_f - 1) then
                        dE = costhet(config_f(2, nm2, irep), nv1) + costhet(nv1, nv2) - costhet(config_f(2, nm2, irep), cm2) &
                                - costhet(cm2, cn2)
                    else
                        dE = costhet(config_f(2, nm2, irep), nv1) + costhet(nv1, nv2) + costhet(nv2, config_f(2, np1, irep)) &
                                - costhet(config_f(2, nm2, irep), cm2) - costhet(cm2, cn2) - costhet(cn2, config_f(2, np1, irep))
                    end if
                    if (interaction_sites_state_f(n, irep) > 0) &
                        dE = dE + Ei_f * (bittable_f(14, v, irep) - bittable_f(14, en, irep))
                    if ((contact_f(3, n, irep) == -1) .and. (id < Nchain_f)) then
                        if (connec(opp(contact_f(2, n, irep)), config_f(2, id, irep), 1) /= 0) dE = dE - Ea_f
                        if (connec(opp(cc), config_f(2, id, irep), 1) /= 0) dE = dE + Ea_f
                    else if ((contact_f(3, n, irep) == 1) .and. (id > 1)) then
                        if (connec(contact_f(2, n, irep), opp(config_f(2, id - 1, irep)), 1) /= 0) dE = dE - Ea_f
                        if (connec(cc, opp(config_f(2, id - 1, irep)), 1) /= 0) dE = dE + Ea_f
                    end if
                    if (contact_f(3, n + 1, irep) == -1) then
                        if (connec(opp(cn2), contact_f(2, n + 1, irep), 1) /= 0) dE = dE - Ea_f
                        if (connec(opp(nv2), contact_f(2, n + 1, irep), 1) /= 0) dE = dE + Ea_f
                    end if
                    if (contact_f(3, n - 1, irep) == 1) then
                        if (connec(cm2, contact_f(2, n - 1, irep), 1) /= 0) dE = dE - Ea_f
                        if (connec(nv1, contact_f(2, n - 1, irep), 1) /= 0) dE = dE + Ea_f
                    end if
                    if (repro_rng_next(rng_state_f(1, irep)) < exp(-dE)) then
                        bittable_f(1, en, irep) = bittable_f(1, en, irep) - 1
                        bittable_f(1, v, irep) = bittable_f(1, v, irep) + 1
                        bittable_f(14, en, irep) = bittable_f(14, en, irep) - 1
                        bittable_f(14, v, irep) = bittable_f(14, v, irep) + 1
                        do j = 2, 13
                            a = bittable_f(j, en, irep)
                            bittable_f(14, a, irep) = bittable_f(14, a, irep) - 1
                            a = bittable_f(j, v, irep)
                            bittable_f(14, a, irep) = bittable_f(14, a, irep) + 1
                        end do
                        config_f(1, n, irep) = v
                        config_f(2, n - 1, irep) = nv1
                        config_f(2, n, irep) = nv2
                        if (id /= 0) then
                            contact_f(2, n, irep) = cc
                            contact_f(2, id, irep) = opp(cc)
                        end if
                        dr_f(1, n, irep) = dr_f(1, n, irep) + voisxyz(1, nv1) - voisxyz(1, cm2)
                        dr_f(2, n, irep) = dr_f(2, n, irep) + voisxyz(2, nv1) - voisxyz(2, cm2)
                        dr_f(3, n, irep) = dr_f(3, n, irep) + voisxyz(3, nv1) - voisxyz(3, cm2)
                    end if
                end if
            end if
        end if
    end subroutine trialmovetad_flat

    subroutine trialmoveex_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: n, iv, s, id, con(3), strand
        real*8 :: rv, fc
        real :: impermeability
        rv = repro_rng_next(rng_state_f(1, irep))
        n = int(Nchain_f * rv) + 1
        s = contact_f(3, n, irep)
        if ((s == 0) .or. (n == 1) .or. (n == Nchain_f)) return
        strand = (s + 1) / 2 + 1
        impermeability = abs(boundary_f(strand, n, irep))
        fc = (1. - impermeability)**(1. / real(ikm_f))
        do iv = 1, ikm_f
            if (repro_rng_next(rng_state_f(1, irep)) >= (km_f * fc)) return
        end do
        if (s == -1) then
            if (connec(1, config_f(2, n - 1, irep), contact_f(2, n, irep)) == 0) return
            if (contact_f(1, n - 1, irep) /= 0) then
                if (.not. z_loop_f) return
                if ((contact_f(3, n - 1, irep) == -1) .or. (n == 2)) return
                if (connec(1, opp(config_f(2, n - 1, irep)), contact_f(2, n - 1, irep)) == 0) return
                con = contact_f(:, n - 1, irep)
                iv = connec(1, config_f(2, n - 1, irep), contact_f(2, n, irep))
                id = contact_f(1, n, irep)
                contact_f(1, n - 1, irep) = id
                contact_f(2, n - 1, irep) = iv
                contact_f(3, n - 1, irep) = -1
                contact_f(1, id, irep) = n - 1
                contact_f(2, id, irep) = opp(iv)
                iv = connec(1, opp(config_f(2, n - 1, irep)), con(2))
                id = con(1)
                contact_f(1, n, irep) = id
                contact_f(2, n, irep) = iv
                if (unidirectional_f) then
                    contact_f(3, n, irep) = con(3)
                else
                    contact_f(3, n, irep) = 1
                end if
                contact_f(1, id, irep) = n
                contact_f(2, id, irep) = opp(iv)
            else
                iv = connec(1, config_f(2, n - 1, irep), contact_f(2, n, irep))
                id = contact_f(1, n, irep)
                contact_f(1, n - 1, irep) = id
                contact_f(2, n - 1, irep) = iv
                contact_f(3, n - 1, irep) = -1
                contact_f(1, id, irep) = n - 1
                contact_f(2, id, irep) = opp(iv)
                contact_f(:, n, irep) = 0
            end if
        else if (s == 1) then
            if (connec(1, opp(config_f(2, n, irep)), contact_f(2, n, irep)) == 0) return
            if (contact_f(1, n + 1, irep) /= 0) then
                if (.not. z_loop_f) return
                if ((contact_f(3, n + 1, irep) == 1) .or. (n == Nchain_f - 1)) return
                if (connec(1, config_f(2, n + 1, irep), contact_f(2, n + 1, irep)) == 0) return
                con = contact_f(:, n + 1, irep)
                iv = connec(1, opp(config_f(2, n, irep)), contact_f(2, n, irep))
                id = contact_f(1, n, irep)
                contact_f(1, n + 1, irep) = id
                contact_f(2, n + 1, irep) = iv
                contact_f(3, n + 1, irep) = 1
                contact_f(1, id, irep) = n + 1
                contact_f(2, id, irep) = opp(iv)
                iv = connec(1, config_f(2, n + 1, irep), con(2))
                id = con(1)
                contact_f(1, n, irep) = id
                contact_f(2, n, irep) = iv
                if (unidirectional_f) then
                    contact_f(3, n, irep) = con(3)
                else
                    contact_f(3, n, irep) = -1
                end if
                contact_f(1, id, irep) = n
                contact_f(2, id, irep) = opp(iv)
            else
                iv = connec(1, opp(config_f(2, n, irep)), contact_f(2, n, irep))
                id = contact_f(1, n, irep)
                contact_f(1, n + 1, irep) = id
                contact_f(2, n + 1, irep) = iv
                contact_f(3, n + 1, irep) = 1
                contact_f(1, id, irep) = n + 1
                contact_f(2, id, irep) = opp(iv)
                contact_f(:, n, irep) = 0
            end if
        end if
    end subroutine trialmoveex_flat

    subroutine trialbound_flat(irep)
        !$acc routine seq
        integer, intent(in) :: irep
        integer :: n, id, j, d
        real :: kbp
        real*8 :: rv
        rv = repro_rng_next(rng_state_f(1, irep))
        n = int(Nchain_f * rv) + 1
        id = contact_f(1, n, irep)
        if (id == 0) then
            kbp = kb_f * (loading_sites_factor_f(n, irep)**(1. / real(ikb_f)))
            do j = 1, ikb_f
                if (repro_rng_next(rng_state_f(1, irep)) >= kbp) return
            end do
            d = int(2 * repro_rng_next(rng_state_f(1, irep))) + 1
            if (d == 1) then
                if (n > 1) then
                    if (contact_f(1, n - 1, irep) == 0) then
                        contact_f(1, n, irep) = n - 1
                        contact_f(2, n, irep) = opp(config_f(2, n - 1, irep))
                        contact_f(1, n - 1, irep) = n
                        contact_f(2, n - 1, irep) = config_f(2, n - 1, irep)
                        if (unidirectional_f) then
                            if (repro_rng_next(rng_state_f(1, irep)) < 0.5) then
                                contact_f(3, n, irep) = 2
                                contact_f(3, n - 1, irep) = -1
                            else
                                contact_f(3, n, irep) = 1
                                contact_f(3, n - 1, irep) = -2
                            end if
                        else
                            contact_f(3, n, irep) = 1
                            contact_f(3, n - 1, irep) = -1
                        end if
                        Nleffree_f(irep) = Nleffree_f(irep) - 1
                    end if
                end if
            else if (n < Nchain_f) then
                if (contact_f(1, n + 1, irep) == 0) then
                    contact_f(1, n, irep) = n + 1
                    contact_f(2, n, irep) = config_f(2, n, irep)
                    contact_f(1, n + 1, irep) = n
                    contact_f(2, n + 1, irep) = opp(config_f(2, n, irep))
                    if (unidirectional_f) then
                        if (repro_rng_next(rng_state_f(1, irep)) < 0.5) then
                            contact_f(3, n, irep) = -2
                            contact_f(3, n + 1, irep) = 1
                        else
                            contact_f(3, n, irep) = -1
                            contact_f(3, n + 1, irep) = 2
                        end if
                    else
                        contact_f(3, n, irep) = -1
                        contact_f(3, n + 1, irep) = 1
                    end if
                    Nleffree_f(irep) = Nleffree_f(irep) - 1
                end if
            end if
        end if
    end subroutine trialbound_flat

    subroutine trialunbound_flat(irep)
        !$acc routine seq
        ! On device we skip write(13); unbind log not written for flat GPU path.
        integer, intent(in) :: irep
        integer :: n, id, j
        n = int(Nchain_f * repro_rng_next(rng_state_f(1, irep))) + 1
        id = contact_f(1, n, irep)
        if (id > 0) then
            do j = 1, iku_f
                if (repro_rng_next(rng_state_f(1, irep)) >= ku_f) return
            end do
            contact_f(:, n, irep) = 0
            contact_f(:, id, irep) = 0
            Nleffree_f(irep) = Nleffree_f(irep) + 1
        end if
    end subroutine trialunbound_flat

end module gpu_flat_mod
