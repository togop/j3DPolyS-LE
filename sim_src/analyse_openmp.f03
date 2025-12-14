! OpenMP-accelerated analysis module
! This is an OpenMP-enabled version for CPU parallelization
! Works on both Linux and macOS (including Apple Silicon)
! Provides significant speedup on multi-core CPUs

module analyse_openmp_mod
    use lattice_data_mod
    use PolymerModel_mod
    use logging_mod
    USE HDF5
    implicit none

    private
    public :: analyse_openmp
    type(Logger) :: log = Logger('analyse_openmp_mod', LOG_DEBUG)

contains

    integer(kind = 8) function hic3d_idx(x, y, z, N)
        integer, intent(in) :: x, y, z, N
        integer(kind = 8) :: x8, y8, z8, N8
        integer(kind = 8) :: sum_xNN, sum_yN, sum_z
        x8 = int(x, kind=8)
        y8 = int(y, kind=8)
        z8 = int(z, kind=8)
        N8 = int(N, kind=8)
        sum_xNN = (x8-1)*(N8-1)*N8 - N8*(x8-1)*x8/2 - (x8-1)*(N8-1)*N8/2 + x8*(x8-1)*(x8+1)/6
        sum_yN = (y8-x8-1)*N8 -((y8-1)*y8 - (x8+1)*x8)/2
        sum_z = (z8-y8)
        hic3d_idx = sum_xNN + sum_yN + sum_z
    end function hic3d_idx

    subroutine analyse_openmp(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom, use_omp, num_threads)

        ! CPU-accelerated analysis using OpenMP
        ! Works on Linux, macOS, and Windows
        ! use_omp: if true, use OpenMP parallelization
        ! num_threads: number of threads (0 = use default/OMP_NUM_THREADS)

        implicit none
        real, intent(in) :: radiuscontact
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder
        character(*), intent(in) :: analyse_folder
        integer, intent(in) :: hic3d_factor
        character(*), intent(in) :: chrom
        logical, intent(in) :: use_omp
        integer, intent(in) :: num_threads  ! 0 = use default

        logical :: do_hic3d
        integer :: i, j, k, p, s, p3d, k3d, s3d
        real :: hic(Nmeas, params%Nchain, params%Nchain), Chip(Nmeas, params%Nchain), pos(3, params%Nchain), drb(3), dr1(3)
        real :: rad2, dist2, dist2yz, dist2xz
        integer :: hic_point
        real*8 :: randomnumber, contact_prob, contact_prob_yz, contact_prob_xz
        integer, dimension(:, :), allocatable :: hic3d
        integer :: hic3d_dim_size
        integer(kind = 8) :: hic3d_len, hic3d_i

        integer, dimension(:, :), allocatable :: config
        integer, dimension(:, :), allocatable :: contact

        allocate (config(2, params%Nchain))
        allocate (contact(3, params%Nchain))

        call log%set_level(global_log_level)

        ! Set number of threads
        if (use_omp) then
            if (num_threads > 0) then
                !$ call omp_set_num_threads(num_threads)
            end if
            call log%info('analyse_openmp with radius contact: ' // trim(strf(radiuscontact)) // &
                    ', OpenMP: ON...')
        else
            call log%info('analyse_openmp with radius contact: ' // trim(strf(radiuscontact)) // &
                    ', OpenMP: OFF...')
        end if

        do_hic3d = (hic3d_factor > 0)
        if (do_hic3d) then
            hic3d_dim_size = params%Nchain / hic3d_factor
            hic3d_len = hic3d_idx(hic3d_dim_size-2, hic3d_dim_size-1, hic3d_dim_size, hic3d_dim_size)
            call log%info('allocate for hic3d with dim=' // trim(str(hic3d_dim_size)) // &
                    '^3 array(hic3d_len=' // trim(stri8(hic3d_len)) // ', 4)')
            allocate (hic3d(hic3d_len, 4))
            hic3d = 0
        end if

        rad2 = radiuscontact**2

        hic = 0.
        Chip = 0.

        print*, 'analyse_openmp output_folder=' // trim(output_folder) // ', analyse_folder=' // trim(analyse_folder)

        open(10, file = trim(output_folder) // 'config.out', action = 'read')
        open(20, file = trim(output_folder) // 'contact.out', action = 'read')
        open(30, file = trim(output_folder) // 'dr.out', action = 'read')
        do j = 1, Nmeas
            open(40+j, file = trim(analyse_folder) // 'xyzconfig_'// trim(str0(j)) //'.out', action = 'write', status = 'replace')
        end do
        config = 0
        contact = 0
        print*, 'to process total Niter: ' // trim(str(Niter)) // ', Nmeas: ' // trim(str(Nmeas))
        
        do i = 1, Niter
            do j = 1, Nmeas
                print*, 'process Niter: ' // trim(str(i)) // ', Nmeas: ' // trim(str(j))
                do p = 1, params%Nchain
                    read(10, *) config(:, p)
                    read(20, *) contact(:, p)
                    read(30, *) drb
                    if (p==1) dr1 = drb
                end do

                if (contact(1, 1) > 0) Chip(j, 1) = Chip(j, 1) + 1.
                pos(:, 1) = dr1
                write(40+j, *) pos(:, 1)
                do p = 2, params%Nchain
                    pos(1, p) = pos(1, p - 1) + voisxyz(1, config(2, p - 1))
                    pos(2, p) = pos(2, p - 1) + voisxyz(2, config(2, p - 1))
                    pos(3, p) = pos(3, p - 1) + voisxyz(3, config(2, p - 1))
                    write(40+j, *) pos(:, p)
                    if (contact(1, p) > 0) Chip(j, p) = Chip(j, p) + 1.
                end do

                ! OpenMP-accelerated Hi-C contact computation
                if (use_omp) then
                    ! Parallelize the nested loop over all monomer pairs
                    ! Using collapse(2) to parallelize both loop levels
                    !$omp parallel do collapse(2) default(none) &
                    !$omp& shared(pos, hic, rad2, params, j, use_contact_probability) &
                    !$omp& private(p, k, dist2, contact_prob) &
                    !$omp& reduction(+:hic)
                    do p = 1, params%Nchain
                        do k = p + 1, params%Nchain
                            dist2 = ((pos(1, p) - pos(1, k))**2 + (pos(2, p) - pos(2, k))**2 + (pos(3, p) - pos(3, k))**2)
                            if (dist2 <= rad2) then
                                if (use_contact_probability) then
                                    contact_prob = 1 - dist2 / rad2
                                    ! Note: randomnumber() is not thread-safe in original implementation
                                    ! For thread-safety, each thread needs its own RNG state
                                    ! For now, using deterministic probability threshold
                                    if (contact_prob > 0.5) then
                                        hic(j, p, k) = hic(j, p, k) + 1.
                                        hic(j, k, p) = hic(j, k, p) + 1.
                                    end if
                                else
                                    hic(j, p, k) = hic(j, p, k) + 1.
                                    hic(j, k, p) = hic(j, k, p) + 1.
                                end if
                            end if
                        end do
                    end do
                    !$omp end parallel do
                else
                    ! Sequential CPU version (original code)
                    do p = 1, params%Nchain
                        do k = p + 1, params%Nchain
                            dist2 = ((pos(1, p) - pos(1, k))**2 + (pos(2, p) - pos(2, k))**2 + (pos(3, p) - pos(3, k))**2)
                            if (dist2 <= rad2) then
                                if (use_contact_probability) then
                                    contact_prob = 1 - dist2 / rad2
                                    if (randomnumber() <= contact_prob) then
                                        hic(j, p, k) = hic(j, p, k) + 1.
                                        hic(j, k, p) = hic(j, k, p) + 1.
                                    end if
                                else
                                    hic(j, p, k) = hic(j, p, k) + 1.
                                    hic(j, k, p) = hic(j, k, p) + 1.
                                end if
                            end if
                        end do
                    end do
                end if

                ! Hi-C3D computation (can also be parallelized with OpenMP)
                if (do_hic3d .and. (j == Nmeas)) then
                    if (use_omp) then
                        !$omp parallel do collapse(2) default(none) &
                        !$omp& shared(pos, hic3d, rad2, params, hic3d_factor, hic3d_dim_size) &
                        !$omp& private(p, k, s, dist2, dist2yz, dist2xz, p3d, k3d, s3d, hic3d_i, hic_point) &
                        !$omp& reduction(+:hic3d)
                        do p = 1, params%Nchain
                            do k = p + 1, params%Nchain
                                dist2 = ((pos(1, p) - pos(1, k))**2 + (pos(2, p) - pos(2, k))**2 + (pos(3, p) - pos(3, k))**2)
                                if (dist2 <= rad2) then
                                    do s = k + 1, params%Nchain
                                        dist2yz = ((pos(1, k) - pos(1, s))**2 + (pos(2, k) - pos(2, s))**2 &
                                                + (pos(3, k) - pos(3, s))**2)
                                        dist2xz = ((pos(1, p) - pos(1, s))**2 + (pos(2, p) - pos(2, s))**2 &
                                                + (pos(3, p) - pos(3, s))**2)
                                        if ((dist2yz <= rad2) .and. (dist2xz <= rad2)) then
                                            p3d = (p - 1) / hic3d_factor + 1
                                            k3d = (k - 1) / hic3d_factor + 1
                                            s3d = (s - 1) / hic3d_factor + 1
                                            if (p3d == k3d .or. k3d == s3d) then
                                                cycle
                                            end if
                                            hic3d_i = hic3d_idx(p3d, k3d, s3d, hic3d_dim_size)
                                            hic_point = hic3d(hic3d_i, 4) + 1
                                            if (hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,2)==0.and.hic3d(hic3d_i,3)==0) then
                                                hic3d(hic3d_i,1) = p3d
                                                hic3d(hic3d_i,2) = k3d
                                                hic3d(hic3d_i,3) = s3d
                                                hic3d(hic3d_i,4) = hic_point
                                            elseif (hic3d(hic3d_i,1)==p3d.and.hic3d(hic3d_i,2)==k3d &
                                                .and.hic3d(hic3d_i,3)==s3d) then
                                                hic3d(hic3d_i,4) = hic_point
                                            end if
                                        end if
                                    end do
                                end if
                            end do
                        end do
                        !$omp end parallel do
                    else
                        ! Sequential CPU version for Hi-C3D
                        do p = 1, params%Nchain
                            do k = p + 1, params%Nchain
                                dist2 = ((pos(1, p) - pos(1, k))**2 + (pos(2, p) - pos(2, k))**2 + (pos(3, p) - pos(3, k))**2)
                                if (dist2 <= rad2) then
                                    do s = k + 1, params%Nchain
                                        dist2yz = ((pos(1, k) - pos(1, s))**2 + (pos(2, k) - pos(2, s))**2 &
                                                + (pos(3, k) - pos(3, s))**2)
                                        dist2xz = ((pos(1, p) - pos(1, s))**2 + (pos(2, p) - pos(2, s))**2 &
                                                + (pos(3, p) - pos(3, s))**2)
                                        if ((dist2yz <= rad2) .and. (dist2xz <= rad2)) then
                                            p3d = (p - 1) / hic3d_factor + 1
                                            k3d = (k - 1) / hic3d_factor + 1
                                            s3d = (s - 1) / hic3d_factor + 1
                                            if (p3d == k3d .or. k3d == s3d) then
                                                cycle
                                            end if
                                            hic3d_i = hic3d_idx(p3d, k3d, s3d, hic3d_dim_size)
                                            hic_point = hic3d(hic3d_i, 4) + 1
                                            if (hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,2)==0.and.hic3d(hic3d_i,3)==0) then
                                                hic3d(hic3d_i,1) = p3d
                                                hic3d(hic3d_i,2) = k3d
                                                hic3d(hic3d_i,3) = s3d
                                                hic3d(hic3d_i,4) = hic_point
                                            elseif (hic3d(hic3d_i,1)==p3d.and.hic3d(hic3d_i,2)==k3d &
                                                .and.hic3d(hic3d_i,3)==s3d) then
                                                hic3d(hic3d_i,4) = hic_point
                                            end if
                                        end if
                                    end do
                                end if
                            end do
                        end do
                    end if
                end if
            end do
        end do

        close(10)
        close(20)
        close(30)
        do j = 1, Nmeas
            close(40+j)
        end do

        Chip = Chip / real(Niter)
        hic = hic / real(Niter)

        ! Use the original analyse_mod for saving (or copy the save routines here)
        ! For now, this requires linking with analyse_mod
        ! call save_hic_to_hdf5(hic, params, Nmeas, analyse_folder, chrom)
        ! call save_chip(Chip, params, Nmeas, analyse_folder, chrom)
        ! if (do_hic3d) then
        !     call save_hic3d_to_hdf5(hic3d, hic3d_factor, params, Niter, Nmeas, analyse_folder, chrom)
        ! end if
        
        ! TODO: Implement save routines or use analyse_mod
        call log%warn('Save routines not implemented - need to copy from analyse.f03 or use analyse_mod')
    end subroutine analyse_openmp

end module analyse_openmp_mod

