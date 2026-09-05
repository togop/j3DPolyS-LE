! Optional CUDA Monte-Carlo backend (ISO_C_BINDING).
! When compiled without -DUSE_CUDA the public API is a no-op stub.

module mc_cuda_mod
    use iso_c_binding
    use PolymerModel_mod
    use logging_mod
#ifdef USE_CUDA
    use lattice_data_mod
#endif
    implicit none
    private
    public :: mc_cuda_is_compiled, mc_cuda_has_device, run_cuda_rank

    type(Logger) :: log = Logger('mc_cuda_mod', LOG_INFO)

contains

#ifndef USE_CUDA

    integer function mc_cuda_is_compiled()
        mc_cuda_is_compiled = 0
    end function mc_cuda_is_compiled

    integer function mc_cuda_has_device()
        mc_cuda_has_device = 0
    end function mc_cuda_has_device

    subroutine run_cuda_rank(models, ntraj, rank, Ninter, Nmeas, burnin, burnout, burnoutM, seed, success)
        type(PolymerModel), intent(inout) :: models(:)
        integer, intent(in) :: ntraj, rank, Ninter, Nmeas, burnin, burnout, burnoutM, seed
        logical, intent(out) :: success
        success = .false.
        if (ntraj < 0 .or. rank < 0 .or. Ninter < 0 .or. Nmeas < 0) return
        if (burnin < 0 .or. burnout < 0 .or. burnoutM < 0 .or. seed == 0) return
        if (size(models) < 0) return
    end subroutine run_cuda_rank

#else

    type, bind(C) :: McCudaParamsC
        integer(c_int) :: T, N, naddr, nlef0
        integer(c_int) :: iku, ikm, ikb
        integer(c_int) :: z_loop, unidirectional, seed
        real(c_float) :: kb, ku, km, Ea, Ei
    end type McCudaParamsC

    interface
        function mc_cuda_compiled() bind(C, name='mc_cuda_compiled')
            import
            integer(c_int) :: mc_cuda_compiled
        end function mc_cuda_compiled

        function mc_cuda_available() bind(C, name='mc_cuda_available')
            import
            integer(c_int) :: mc_cuda_available
        end function mc_cuda_available

        function mc_cuda_create(out, p, opp, voisxyz, costhet, voisnn, connec, &
                bittable, boundary, loadfac, istate) bind(C, name='mc_cuda_create')
            import
            integer(c_int) :: mc_cuda_create
            type(c_ptr) :: out
            type(McCudaParamsC), intent(in) :: p
            integer(c_int), intent(in) :: opp(*)
            real(c_float), intent(in) :: voisxyz(*)
            real(c_float), intent(in) :: costhet(*)
            integer(c_int), intent(in) :: voisnn(*)
            integer(c_int), intent(in) :: connec(*)
            integer(c_int), intent(in) :: bittable(*)
            real(c_float), intent(in) :: boundary(*)
            real(c_float), intent(in) :: loadfac(*)
            integer(c_int), intent(in) :: istate(*)
        end function mc_cuda_create

        function mc_cuda_upload_traj(s, t, config, contact, bittable, dr, nfree, traj_index) &
                bind(C, name='mc_cuda_upload_traj')
            import
            integer(c_int) :: mc_cuda_upload_traj
            type(c_ptr), value :: s
            integer(c_int), value :: t
            integer(c_int), intent(in) :: config(*)
            integer(c_int), intent(in) :: contact(*)
            integer(c_int), intent(in) :: bittable(*)
            real(c_float), intent(in) :: dr(*)
            integer(c_int), value :: nfree
            integer(c_int), value :: traj_index
        end function mc_cuda_upload_traj

        function mc_cuda_set_nsteps_t(s, nsteps_t) bind(C, name='mc_cuda_set_nsteps_t')
            import
            integer(c_int) :: mc_cuda_set_nsteps_t
            type(c_ptr), value :: s
            type(c_ptr), value :: nsteps_t
        end function mc_cuda_set_nsteps_t

        function mc_cuda_run(s, mode, nsteps, step0, coherent) bind(C, name='mc_cuda_run')
            import
            integer(c_int) :: mc_cuda_run
            type(c_ptr), value :: s
            integer(c_int), value :: mode, nsteps, step0, coherent
        end function mc_cuda_run

        function mc_cuda_download_traj(s, t, config, contact, bittable, dr, nfree) &
                bind(C, name='mc_cuda_download_traj')
            import
            integer(c_int) :: mc_cuda_download_traj
            type(c_ptr), value :: s
            integer(c_int), value :: t
            integer(c_int), intent(out) :: config(*)
            integer(c_int), intent(out) :: contact(*)
            type(c_ptr), value :: bittable
            real(c_float), intent(out) :: dr(*)
            integer(c_int), intent(inout) :: nfree
        end function mc_cuda_download_traj

        function mc_cuda_download_unbind(s, t, max_events, n_out, id_out, count) &
                bind(C, name='mc_cuda_download_unbind')
            import
            integer(c_int) :: mc_cuda_download_unbind
            type(c_ptr), value :: s
            integer(c_int), value :: t, max_events
            integer(c_int), intent(out) :: n_out(*)
            integer(c_int), intent(out) :: id_out(*)
            integer(c_int), intent(out) :: count
        end function mc_cuda_download_unbind

        subroutine mc_cuda_destroy(s) bind(C, name='mc_cuda_destroy')
            import
            type(c_ptr), value :: s
        end subroutine mc_cuda_destroy
    end interface

    integer function mc_cuda_is_compiled()
        mc_cuda_is_compiled = int(mc_cuda_compiled())
    end function mc_cuda_is_compiled

    integer function mc_cuda_has_device()
        mc_cuda_has_device = int(mc_cuda_available())
    end function mc_cuda_has_device

    subroutine run_cuda_rank(models, ntraj, rank, Ninter, Nmeas, burnin, burnout, burnoutM, seed, success)
        type(PolymerModel), intent(inout) :: models(:)
        integer, intent(in) :: ntraj, rank, Ninter, Nmeas, burnin, burnout, burnoutM, seed
        logical, intent(out) :: success

        type(c_ptr) :: st
        type(McCudaParamsC) :: p
        integer :: rc, t, j, k, meas, n_mixed, burn_Nmeas, step0
        integer :: N, naddr, coherent
        integer, allocatable, target :: burn_n(:)
        integer, allocatable :: snap_cfg(:, :, :, :), snap_con(:, :, :, :), snap_nlef(:, :)
        real, allocatable :: snap_dr(:, :, :, :)
        integer, parameter :: evcap = 262144
        integer, parameter :: evbuf = 65536
        integer, allocatable :: all_ev_n(:, :), all_ev_id(:, :), all_ev_cnt(:)
        integer :: bufn(evbuf), bufid(evbuf)
        integer :: traj_index
        real*8, external :: randomnumber
        real*8 :: r

        success = .false.
        call log%set_level(global_log_level)

        if (ntraj <= 0) return
        N = models(1)%Nchain
        naddr = size(models(1)%bittable, 2)
        coherent = 1

        p%T = ntraj
        p%N = N
        p%naddr = naddr
        p%nlef0 = models(1)%Nleffree
        p%iku = models(1)%iku
        p%ikm = models(1)%ikm
        p%ikb = models(1)%ikb
        p%z_loop = 0
        if (models(1)%z_loop) p%z_loop = 1
        p%unidirectional = 0
        if (models(1)%unidirectional) p%unidirectional = 1
        p%seed = seed
        p%kb = models(1)%kb
        p%ku = models(1)%ku
        p%km = models(1)%km
        p%Ea = models(1)%Ea
        p%Ei = models(1)%Ei

        st = c_null_ptr
        rc = mc_cuda_create(st, p, opp, voisxyz, costhet, voisnn, connec, &
                models(1)%bittable, models(1)%boundary, models(1)%loading_sites_factor, &
                models(1)%interaction_sites_state)
        if (rc /= 0 .or. (.not. c_associated(st))) then
            call log%warn('mc_cuda_create failed')
            return
        end if

        allocate(snap_cfg(2, N, Nmeas, ntraj))
        allocate(snap_con(3, N, Nmeas, ntraj))
        allocate(snap_dr(3, N, Nmeas, ntraj))
        allocate(snap_nlef(Nmeas, ntraj))
        allocate(all_ev_n(evcap, ntraj))
        allocate(all_ev_id(evcap, ntraj))
        allocate(all_ev_cnt(ntraj))
        all_ev_cnt = 0

        do t = 1, ntraj
            traj_index = rank * ntraj + t
            rc = mc_cuda_upload_traj(st, t - 1, models(t)%config, models(t)%contact, &
                    models(t)%bittable, models(t)%dr, models(t)%Nleffree, traj_index)
            if (rc /= 0) then
                call log%warn('mc_cuda_upload_traj failed')
                call mc_cuda_destroy(st)
                return
            end if
            snap_cfg(:, :, 1, t) = models(t)%config
            snap_con(:, :, 1, t) = models(t)%contact
            snap_dr(:, :, 1, t) = models(t)%dr
            snap_nlef(1, t) = models(t)%Nleffree
        end do

        call log%info('CUDA batch start: T=' // trim(str(ntraj)) // &
                ', N=' // trim(str(N)) // ', coherent_moves=true')

        meas = 1
        burn_Nmeas = 0
        if (burnin > 0) burn_Nmeas = burn_Nmeas + 1
        if (burnout > 0) burn_Nmeas = burn_Nmeas + 1
        if (burnoutM > 0) burn_Nmeas = burn_Nmeas + burnoutM

        if (burnin > 0) then
            allocate(burn_n(ntraj))
            do t = 1, ntraj
                r = randomnumber()
                burn_n(t) = burnin + int(r * burnin)
            end do
            rc = mc_cuda_set_nsteps_t(st, c_loc(burn_n(1)))
            if (rc /= 0) then
                call mc_cuda_destroy(st)
                return
            end if
            rc = mc_cuda_run(st, 1, maxval(burn_n), 0, coherent)
            rc = rc + mc_cuda_set_nsteps_t(st, c_null_ptr)
            deallocate(burn_n)
            if (rc /= 0) then
                call log%warn('CUDA burn-in kernel failed')
                call mc_cuda_destroy(st)
                return
            end if
            meas = meas + 1
            call download_all(st, models, ntraj, N, meas, snap_cfg, snap_con, snap_dr, &
                    snap_nlef, all_ev_n, all_ev_id, all_ev_cnt, evcap, bufn, bufid, evbuf)
            call log%info('CUDA burn-in measurement stored for T=' // trim(str(ntraj)))
        end if

        n_mixed = Nmeas - 1 - burn_Nmeas
        step0 = 0
        do j = 1, n_mixed
            rc = mc_cuda_run(st, 0, Ninter, step0, coherent)
            if (rc /= 0) then
                call log%warn('CUDA mixed kernel failed')
                call mc_cuda_destroy(st)
                return
            end if
            step0 = step0 + Ninter
            meas = meas + 1
            call download_all(st, models, ntraj, N, meas, snap_cfg, snap_con, snap_dr, &
                    snap_nlef, all_ev_n, all_ev_id, all_ev_cnt, evcap, bufn, bufid, evbuf)
            call log%info('CUDA measurement ' // trim(str(j)) // ' stored for T=' // trim(str(ntraj)))
        end do

        if (burnout > 0) then
            rc = mc_cuda_run(st, 2, 0, 0, coherent)
            rc = rc + mc_cuda_run(st, 1, burnout, 0, coherent)
            if (rc /= 0) then
                call log%warn('CUDA burnout kernel failed')
                call mc_cuda_destroy(st)
                return
            end if
            meas = meas + 1
            call download_all(st, models, ntraj, N, meas, snap_cfg, snap_con, snap_dr, &
                    snap_nlef, all_ev_n, all_ev_id, all_ev_cnt, evcap, bufn, bufid, evbuf)
            call log%info('CUDA fix burn-out measurement stored')
        end if

        if (burnoutM > 0) then
            rc = mc_cuda_run(st, 2, 0, 0, coherent)
            if (rc /= 0) then
                call mc_cuda_destroy(st)
                return
            end if
            do j = 1, burnoutM
                rc = mc_cuda_run(st, 1, Ninter, 0, coherent)
                if (rc /= 0) then
                    call mc_cuda_destroy(st)
                    return
                end if
                meas = meas + 1
                call download_all(st, models, ntraj, N, meas, snap_cfg, snap_con, snap_dr, &
                        snap_nlef, all_ev_n, all_ev_id, all_ev_cnt, evcap, bufn, bufid, evbuf)
                call log%info('CUDA burn-out measurement ' // trim(str(j)) // ' stored')
            end do
        end if

        if (meas /= Nmeas) then
            call log%warn('CUDA snapshot count ' // trim(str(meas)) // ' != Nmeas ' // trim(str(Nmeas)))
        end if

        do t = 1, ntraj
            do j = 1, min(meas, Nmeas)
                models(t)%config = snap_cfg(:, :, j, t)
                models(t)%contact = snap_con(:, :, j, t)
                models(t)%dr = snap_dr(:, :, j, t)
                models(t)%Nleffree = snap_nlef(j, t)
                call models(t)%output()
            end do
            do k = 1, all_ev_cnt(t)
                write(13, *) all_ev_n(k, t), all_ev_id(k, t)
            end do
            call flush(13)
            call models(t)%erase()
        end do

        call mc_cuda_destroy(st)
        success = .true.
        call log%info('CUDA batch finished, wrote snapshots in trajectory-major order')
    end subroutine run_cuda_rank

    subroutine download_all(st, models, ntraj, N, meas, snap_cfg, snap_con, snap_dr, &
            snap_nlef, all_ev_n, all_ev_id, all_ev_cnt, evcap, bufn, bufid, evbuf)
        type(c_ptr), intent(in) :: st
        type(PolymerModel), intent(in) :: models(:)
        integer, intent(in) :: ntraj, N, meas, evcap, evbuf
        integer, intent(inout) :: snap_cfg(:, :, :, :), snap_con(:, :, :, :), snap_nlef(:, :)
        real, intent(inout) :: snap_dr(:, :, :, :)
        integer, intent(inout) :: all_ev_n(:, :), all_ev_id(:, :), all_ev_cnt(:)
        integer, intent(inout) :: bufn(:), bufid(:)
        integer :: t, rc, nfree, evc, i
        integer, allocatable :: tmp_cfg(:, :), tmp_con(:, :)
        real, allocatable :: tmp_dr(:, :)

        allocate(tmp_cfg(2, N), tmp_con(3, N), tmp_dr(3, N))
        do t = 1, ntraj
            nfree = models(t)%Nleffree
            rc = mc_cuda_download_traj(st, t - 1, tmp_cfg, tmp_con, c_null_ptr, tmp_dr, nfree)
            if (rc == 0) then
                snap_cfg(:, :, meas, t) = tmp_cfg
                snap_con(:, :, meas, t) = tmp_con
                snap_dr(:, :, meas, t) = tmp_dr
                snap_nlef(meas, t) = nfree
            end if
            evc = 0
            rc = mc_cuda_download_unbind(st, t - 1, evbuf, bufn, bufid, evc)
            if (rc == 0) then
                do i = 1, evc
                    if (all_ev_cnt(t) < evcap) then
                        all_ev_cnt(t) = all_ev_cnt(t) + 1
                        all_ev_n(all_ev_cnt(t), t) = bufn(i)
                        all_ev_id(all_ev_cnt(t), t) = bufid(i)
                    end if
                end do
            end if
        end do
        deallocate(tmp_cfg, tmp_con, tmp_dr)
    end subroutine download_all

#endif

end module mc_cuda_mod
