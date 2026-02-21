! GPU replica-parallel execution for do_simulation.
!
! Run multiple trajectories (replicas) in parallel. Each replica has its own
! state and RNG stream so Monte Carlo stochasticity is preserved.
! See sim_src/GPU_DESIGN.md for full design.
!
! Build with OpenACC (e.g. nvfortran -acc) to enable GPU. Without -acc the
! !$acc directives are treated as comments and the loop runs on host.
! The main program can disable GPU at runtime via --no-gpu (use_gpu=.false.).

module gpu_replica_mod
    use PolymerModel_mod
    use logging_mod
    implicit none
    private
    public :: run_measurement_interval_replicas, run_burnin_replicas, run_burnout_only_replicas, &
            run_burnoutM_one_replicas, do_simulation_replicas

contains

    ! Run one measurement interval (Ninter * Ntrial steps) for all replicas.
    subroutine run_measurement_interval_replicas(replicas, Ninter, use_gpu)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: Ninter
        logical, intent(in) :: use_gpu
        integer :: r, nrep

        nrep = size(replicas)
        if (nrep <= 0) return

        if (use_gpu) then
            !$acc parallel loop present(replicas) copyin(Ninter)
            do r = 1, nrep
                call replicas(r)%advance_one_measurement_interval(Ninter)
            end do
            !$acc end parallel loop
        else
            !$omp parallel do schedule(static) default(none) shared(replicas, Ninter, nrep)
            do r = 1, nrep
                call replicas(r)%advance_one_measurement_interval(Ninter)
            end do
            !$omp end parallel do
        end if
    end subroutine run_measurement_interval_replicas

    subroutine run_burnin_replicas(replicas, burnin, use_gpu)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: burnin
        logical, intent(in) :: use_gpu
        integer :: r, nrep

        nrep = size(replicas)
        if (nrep <= 0 .or. burnin <= 0) return

        if (use_gpu) then
            !$acc parallel loop present(replicas) copyin(burnin)
            do r = 1, nrep
                call replicas(r)%run_burnin_phase(burnin)
            end do
            !$acc end parallel loop
        else
            !$omp parallel do schedule(static) default(none) shared(replicas, burnin, nrep)
            do r = 1, nrep
                call replicas(r)%run_burnin_phase(burnin)
            end do
            !$omp end parallel do
        end if
    end subroutine run_burnin_replicas

    subroutine run_burnout_only_replicas(replicas, burnout, use_gpu)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: burnout
        logical, intent(in) :: use_gpu
        integer :: r, nrep

        nrep = size(replicas)
        if (nrep <= 0 .or. burnout <= 0) return

        if (use_gpu) then
            !$acc parallel loop present(replicas) copyin(burnout)
            do r = 1, nrep
                call replicas(r)%run_burnout_only_phase(burnout)
            end do
            !$acc end parallel loop
        else
            !$omp parallel do schedule(static) default(none) shared(replicas, burnout, nrep)
            do r = 1, nrep
                call replicas(r)%run_burnout_only_phase(burnout)
            end do
            !$omp end parallel do
        end if
    end subroutine run_burnout_only_replicas

    subroutine run_burnoutM_one_replicas(replicas, Ninter, use_gpu)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: Ninter
        logical, intent(in) :: use_gpu
        integer :: r, nrep

        nrep = size(replicas)
        if (nrep <= 0) return

        if (use_gpu) then
            !$acc parallel loop present(replicas) copyin(Ninter)
            do r = 1, nrep
                call replicas(r)%run_burnoutM_one_phase(Ninter)
            end do
            !$acc end parallel loop
        else
            !$omp parallel do schedule(static) default(none) shared(replicas, Ninter, nrep)
            do r = 1, nrep
                call replicas(r)%run_burnoutM_one_phase(Ninter)
            end do
            !$omp end parallel do
        end if
    end subroutine run_burnoutM_one_replicas

    ! Full replica-parallel do_simulation: initial output, burn-in, main loop, burnout, then erase.
    ! Output is written on host after each phase/measurement (measurement-major order in files).
    subroutine do_simulation_replicas(replicas, nrep, Ninter, Nmeas, burnin, burnout, burnoutM, use_gpu)
        type(PolymerModel), intent(inout) :: replicas(:)
        integer, intent(in) :: nrep, Ninter, Nmeas, burnin, burnout, burnoutM
        logical, intent(in) :: use_gpu
        integer :: burn_Nmeas, j, r
        type(Logger) :: log

        log = Logger('gpu_replica_mod', LOG_INFO)
        burn_Nmeas = 0
        if (burnin > 0) burn_Nmeas = burn_Nmeas + 1
        if (burnout > 0) burn_Nmeas = burn_Nmeas + 1
        burn_Nmeas = burn_Nmeas + burnoutM

        ! Initial output (replicas on host)
        do r = 1, nrep
            call replicas(r)%output()
        end do

        if (use_gpu) then
            ! Keep replica data on device; copy to host only when we need to output (reduces transfer vs copy every kernel)
            !$acc enter data copyin(replicas)
        end if

        ! Burn-in
        if (burnin > 0) then
            call run_burnin_replicas(replicas, burnin, use_gpu)
            if (use_gpu) then
                !$acc update host(replicas)
            end if
            do r = 1, nrep
                call replicas(r)%output()
            end do
        end if

        ! Main measurement loop
        do j = 1, Nmeas - 1 - burn_Nmeas
            call run_measurement_interval_replicas(replicas, Ninter, use_gpu)
            if (use_gpu) then
                !$acc update host(replicas)
            end if
            call log%info('Replica batch, Measurement: ' // trim(str(j)))
            call flush(6)
            do r = 1, nrep
                call replicas(r)%output()
            end do
            call flush(13)
        end do

        ! Burnout block
        if (burnout > 0) then
            call run_burnout_only_replicas(replicas, burnout, use_gpu)
            if (use_gpu) then
                !$acc update host(replicas)
            end if
            do r = 1, nrep
                call replicas(r)%output()
            end do
            call flush(13)
        end if

        ! BurnoutM steps
        if (burnoutM > 0) then
            do j = 1, burnoutM
                call run_burnoutM_one_replicas(replicas, Ninter, use_gpu)
                if (use_gpu) then
                    !$acc update host(replicas)
                end if
                call log%info('Replica batch, burn-out Measurement: ' // trim(str(j)))
                do r = 1, nrep
                    call replicas(r)%output()
                end do
                call flush(13)
            end do
        end if

        if (use_gpu) then
            !$acc exit data copyout(replicas)
        end if

        do r = 1, nrep
            call replicas(r)%erase()
        end do
    end subroutine do_simulation_replicas

end module gpu_replica_mod
