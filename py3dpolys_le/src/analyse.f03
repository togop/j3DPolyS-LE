module analyse_mod
    use lattice_data_mod
    use PolymerModel_mod
    use logging_mod
    !    USE ISO_C_BINDING
    USE HDF5
    implicit none

    private
    public :: str, analyseradius
    type(Logger) :: log = Logger('analyse_mod', LOG_DEBUG)

contains


    integer function hic3d_idx(x, y, z, N)
        integer, intent(in) :: x, y, z, N
        hic3d_idx = (x-1)*(N-1)*N - N*(x-1)*x/2 - (x-1)*(N-1)*N/2 + x*(x-1)*(x+1)/6 &  ! sum_xNN
                + (y-x-1)*N -((y-1)*y - (x+1)*x)/2 &  ! sum_yN
                + (z-y)  ! sum_z
        !call log%debug('hic3d_idx: x=' // trim(str(x)) // ', y=' // trim(str(y)) // ', z='// trim(str(z)) // &
        !       ', N=' // trim(str(N)) // ')=' // trim(str(hic3d_idx)) )
    end function hic3d_idx

    subroutine analyseradius(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor)

        !analyse data to estimate the Hi-C map for a given radius of contact and the Chip-profile of LEF legs
        ! also compute the xyz coordinates

        implicit none
        real, intent(in) :: radiuscontact   ! default = 1. !in lattice unit (recall: 1 lattice unit=70nm)
        logical, intent(in) :: use_contact_probability
        class (ModelParameters), intent(in) :: params
        integer, intent(in) :: Niter, Nmeas
        character(*), intent(in) :: output_folder
        character(*), intent(in) :: analyse_folder
        integer, intent(in) :: hic3d_factor

        logical :: do_hic3d
        integer :: i, j, k, p, s, p3d, k3d, s3d
        real :: hic(Nmeas, params%Nchain, params%Nchain), Chip(Nmeas, params%Nchain), pos(3, params%Nchain), drb(3), dr1(3)
        real :: rad2, dist2, dist2yz, dist2xz
        integer :: hic_point
        real*8 :: randomnumber, contact_prob, contact_prob_yz, contact_prob_xz
        integer, dimension(:, :), allocatable :: hic3d  ! TODO mybe use integer and if want normalized do it when saving
        integer :: hic3d_dim_size, hic3d_len, hic3d_i

        integer, dimension(:, :), allocatable :: config
        integer, dimension(:, :), allocatable :: contact

        allocate (config(2, params%Nchain))
        allocate (contact(3, params%Nchain))

        call log%set_level(global_log_level)

        call log%info('analyse with radius contact: ' // trim(strf(radiuscontact)) // '...')

        do_hic3d = (hic3d_factor > 0)
        if (do_hic3d) then
            hic3d_dim_size = params%Nchain / hic3d_factor
            hic3d_len = hic3d_idx(hic3d_dim_size-2, hic3d_dim_size-1, hic3d_dim_size, hic3d_dim_size)  ! hic3d_dim_size * (hic3d_dim_size + 1) * (hic3d_dim_size + 2) / 6
            call log%info('allocate for hic3d with dim=' // trim(str(hic3d_dim_size)) // &
                    '^3 array(hic3d_len=' // trim(str(hic3d_len)) // ', 4)')
            allocate (hic3d(hic3d_len, 4))
            hic3d = 0
            ! allocate (hic3d(hic3d_dim_size, hic3d_dim_size, hic3d_dim_size))
        end if

        rad2 = radiuscontact**2

        hic = 0.
        Chip = 0.

        print*, 'analyseradius output_folder=' // trim(output_folder) // ', analyse_folder=' // trim(analyse_folder)

        open(10, file = trim(output_folder) // 'config.out', action = 'read')
        open(20, file = trim(output_folder) // 'contact.out', action = 'read')
        open(30, file = trim(output_folder) // 'dr.out', action = 'read')
        open(40, file = trim(analyse_folder) // 'xyzconfig.out', action = 'write', status = 'replace')
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
                write(40, *) pos(:, 1)
                do p = 2, params%Nchain
                    pos(1, p) = pos(1, p - 1) + voisxyz(1, config(2, p - 1))
                    pos(2, p) = pos(2, p - 1) + voisxyz(2, config(2, p - 1))
                    pos(3, p) = pos(3, p - 1) + voisxyz(3, config(2, p - 1))
                    write(40, *) pos(:, p)
                    if (contact(1, p) > 0) Chip(j, p) = Chip(j, p) + 1.
                end do
                do p = 1, params%Nchain
                    do k = p + 1, params%Nchain
                        dist2 = ((pos(1, p) - pos(1, k))**2 + (pos(2, p) - pos(2, k))**2 + (pos(3, p) - pos(3, k))**2)
                        if (dist2 <= rad2) then
                            if (use_contact_probability) then
                                contact_prob = 1 - dist2 / rad2   ! alt: (rad2 - dist2) / (rad2 + start1_contact_prob2)
                                if (randomnumber() <= contact_prob) then
                                    hic(j, p, k) = hic(j, p, k) + 1.
                                    hic(j, k, p) = hic(j, k, p) + 1.
                                end if
                                if (do_hic3d .and. (j == Nmeas)) then
                                    do s = k + 1, params%Nchain
                                        dist2yz = ((pos(1, k) - pos(1, s))**2 + (pos(2, k) - pos(2, s))**2 &
                                                + (pos(3, k) - pos(3, s))**2)
                                        dist2xz = ((pos(1, p) - pos(1, s))**2 + (pos(2, p) - pos(2, s))**2 &
                                                + (pos(3, p) - pos(3, s))**2)
                                        if ((dist2yz <= rad2) .and. (dist2xz <= rad2)) then
                                            contact_prob_yz = 1 - dist2yz / rad2   ! alt: (rad2 - dist2) / (rad2 + start1_contact_prob2)
                                            contact_prob_xz = 1 - dist2xz / rad2
                                            if ((randomnumber() <= contact_prob_yz) .and. (randomnumber() <= contact_prob_xz)) then
                                                p3d = (p - 1) / hic3d_factor + 1
                                                k3d = (k - 1) / hic3d_factor + 1
                                                s3d = (s - 1) / hic3d_factor + 1
                                                if (p3d == k3d .or. k3d == s3d) then
                                                    cycle  ! skip the diagonal: leave it 0
                                                end if
                                                ! hic3d cool like format
                                                hic3d_i = hic3d_idx(p3d, k3d, s3d, hic3d_dim_size)
                                                hic_point = hic3d(hic3d_i, 4) + 1
                                                if (hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,1)==0) then
                                                    hic3d(hic3d_i,1) = p3d
                                                    hic3d(hic3d_i,2) = k3d
                                                    hic3d(hic3d_i,3) = s3d
                                                    hic3d(hic3d_i,4) = hic_point
                                                    !call log%debug('init hic3d(hic3d_i:' // trim(str(hic3d_i)) // &
                                                    !        ')= hic(p):' // trim(str(hic3d(hic3d_i,1))) // &
                                                    !        ', hic(k):' // trim(str(hic3d(hic3d_i,2))) // &
                                                    !        ', hic(s):' // trim(str(hic3d(hic3d_i,3))) // &
                                                    !         ', hic(v):' // trim(str(hic_point)))
                                                elseif (hic3d(hic3d_i,1)==p3d.and.hic3d(hic3d_i,2)==k3d &
                                                    .and.hic3d(hic3d_i,3)==s3d) then
                                                    hic3d(hic3d_i,4) = hic_point
                                                else
                                                    call log%error('not mathching p3d:' // trim(str(p3d)) // &
                                                            ',k3d:' // trim(str(k3d)) // ', s3d:' // trim(str(s3d)))
                                                    call log%error('with hic3d(hic3d_i:' // trim(str(hic3d_i)) // &
                                                            ')= hic(p3d):' // trim(str(hic3d(hic3d_i,1))) // &
                                                            ', hic(k3d):' // trim(str(hic3d(hic3d_i,2))) // &
                                                            ', hic(s3d):' // trim(str(hic3d(hic3d_i,3))))
                                                    call exit(1)
                                                end if
                                            end if
                                        end if
                                    end do
                                end if
                            else
                                hic(j, p, k) = hic(j, p, k) + 1.
                                hic(j, k, p) = hic(j, k, p) + 1.
                                if (do_hic3d .and. (j == Nmeas)) then
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
                                                cycle  ! skip the diagonal: leave it 0
                                            end if
                                            ! hic3d cool like format
                                            hic3d_i = hic3d_idx(p3d, k3d, s3d, hic3d_dim_size)
                                            ! call log%info('init hic3d(hic3d_i:' // trim(str(hic3d_i)))
                                            hic_point = hic3d(hic3d_i, 4) + 1
                                            if (hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,1)==0.and.hic3d(hic3d_i,1)==0) then
                                                hic3d(hic3d_i,1) = p3d
                                                hic3d(hic3d_i,2) = k3d
                                                hic3d(hic3d_i,3) = s3d
                                                hic3d(hic3d_i,4) = hic_point
                                                !call log%debug('init hic3d(hic3d_i:' // trim(str(hic3d_i)) // &
                                                !        ')= hic(p):' // trim(str(hic3d(hic3d_i,1))) // &
                                                !        ', hic(k):' // trim(str(hic3d(hic3d_i,2))) // &
                                                !        ', hic(s):' // trim(str(hic3d(hic3d_i,3))) // &
                                                !        ', hic(v):' // trim(str(hic_point)))
                                            elseif (hic3d(hic3d_i,1)==p3d.and.hic3d(hic3d_i,2)==k3d &
                                                .and.hic3d(hic3d_i,3)==s3d) then
                                                hic3d(hic3d_i,4) = hic_point
                                                !call log%debug('update hic3d(hic3d_i:' // trim(str(hic3d_i)) // &
                                                !        ')= hic(p):' // trim(str(hic3d(hic3d_i,1))) // &
                                                !        ', hic(k):' // trim(str(hic3d(hic3d_i,2))) // &
                                                !        ', hic(s):' // trim(str(hic3d(hic3d_i,3))) // &
                                                !        ', hic(v):' // trim(str(hic_point)))
                                            else
                                                call log%error('not mathching p3d:' // trim(str(p3d)) // &
                                                        ',k3d:' // trim(str(k)) // ', s3d:' // trim(str(s3d)))
                                                call log%error('with hic3d(hic3d_i:' // trim(str(hic3d_i)) // &
                                                        ')= hic(p3d):' // trim(str(hic3d(hic3d_i,1))) // &
                                                        ', hic(k3d):' // trim(str(hic3d(hic3d_i,2))) // &
                                                        ', hic(s3d):' // trim(str(hic3d(hic3d_i,3))))
                                                call exit(1)
                                            end if
                                        end if
                                    end do
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end do

        close(10)
        close(20)
        close(30)
        close(40)

        Chip = Chip / real(Niter)
        hic = hic / real(Niter)

        call save_hic_to_hdf5(hic, params, Nmeas, analyse_folder)
        call save_chip(Chip, params, Nmeas, analyse_folder)

        if (do_hic3d) then
            ! to Normlize myaybe do this
            ! hic3d(:,4) = hic3d(:,4) / real(Niter * hic3d_factor**2)
            call save_hic3d_to_hdf5(hic3d, hic3d_factor, params, Niter, Nmeas, analyse_folder)
        end if
    end subroutine analyseradius

    subroutine save_hic_to_hdf5(hic, params, Nmeas, analyse_folder)
        implicit none
        class (ModelParameters), intent(in) :: params
        real, intent(in), TARGET :: hic(Nmeas, params%Nchain, params%Nchain)    ! or: hic(*)
        integer, intent(in) :: Nmeas
        character(*), intent(in) :: analyse_folder

        integer :: rc, i

        !character(1000) :: filename

        ! HDF5 support
        character(1000) :: filename_h5
        INTEGER, PARAMETER :: real_kind_7 = SELECTED_REAL_KIND(6, 37) !should map to REAL*4 on most modern processors
        INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(15, 307) !should map to REAL*8 on most modern processors
        CHARACTER(LEN = 7), PARAMETER :: dataset = "hic_map"  ! Dataset name
        INTEGER(HSIZE_T), DIMENSION(2) :: data_dims   ! Dataset dimensions =  (/params%Nchain, params%Nchain/)
        INTEGER(HID_T) :: file_id       ! File identifier
        INTEGER(HID_T) :: dset_id       ! Dataset identifier
        INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
        INTEGER :: rank = 2                            ! Dataset rank
        INTEGER :: error ! Error flag

        ! TODO split output in different method

        ! HDF5 support
        ! Initialize FORTRAN interface.
        CALL h5open_f(rc)

        data_dims(1) = params%Nchain
        data_dims(2) = params%Nchain

        !allocate(dims(2, params%Nchain))

        do i = 1, Nmeas
            ! output in plain text tsv file
            !            filename = trim(output_folder)//'hic' // '_' // trim(str(i)) // '.out'
            !            open(10, file = filename, action='write', status='replace', iostat = rc)
            !            do p = 1, params%Nchain
            !                write(10, *) hic(i, p, :)
            !            end do
            !            close(10)

            ! output in HDF5 file
            filename_h5 = trim(analyse_folder) // 'hic' // '_' // trim(str0(i)) // '.hdf5'

            ! Create a new file
            CALL h5fcreate_f(filename_h5, H5F_ACC_TRUNC_F, file_id, rc)

            ! Create the dataspace.
            CALL h5screate_simple_f(rank, data_dims, dspace_id, error)

            ! Create the dataset with default properties.
            CALL h5dcreate_f(file_id, dataset, h5kind_to_type(real_kind_7, H5_REAL_KIND), dspace_id, dset_id, error)

            !hic_j = hic(i, :, :)
            !f_ptr = C_LOC(hic(i, :, :))
            CALL h5dwrite_f(dset_id, h5kind_to_type(real_kind_7, H5_REAL_KIND), hic(i, :, :), data_dims, error) !, xfer_prp = plist_id)

            ! Close the dataset.
            CALL h5dclose_f(dset_id, error)

            ! Terminate access to the data space.
            CALL h5sclose_f(dspace_id, error)

            ! Close the file.
            CALL h5fclose_f(file_id, error)

            call log%info('Saved ' // trim(filename_h5))
        end do

        ! Close FORTRAN interface.
        CALL h5close_f(error)
    end subroutine save_hic_to_hdf5

    subroutine save_hic3d_to_hdf5(hic3d, factor, params, Niter, Nmeas, analyse_folder)
        ! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
        implicit none
        class (ModelParameters), intent(in) :: params
        ! real, intent(in), TARGET :: hic3d(:, :, :)
        integer, intent(in), TARGET :: hic3d(:, :)    ! cool like format
        integer, intent(in) :: factor, Niter, Nmeas
        character(*), intent(in) :: analyse_folder

        integer :: hic3d_shape(2)
        integer, parameter :: sim_resolution = 2000
        integer :: resolution

        ! HDF5 support
        integer :: rc
        character(1000) :: filename_h5
        INTEGER, PARAMETER :: real_kind_7 = SELECTED_REAL_KIND(6, 37) !should map to REAL*4 on most modern processors
        INTEGER, PARAMETER :: real_kind_15 = SELECTED_REAL_KIND(15, 307) !should map to REAL*8 on most modern processors
        INTEGER, PARAMETER :: int_kind_1 = SELECTED_INT_KIND(2)  !should map to INTEGER*1 on most modern processors
        INTEGER, PARAMETER :: int_kind_4 = SELECTED_INT_KIND(4)  !should map to INTEGER*2 on most modern processors
        INTEGER, PARAMETER :: int_kind_8 = SELECTED_INT_KIND(9)  !should map to INTEGER*4 on most modern processors
        INTEGER, PARAMETER :: int_kind_16 = SELECTED_INT_KIND(18) !should map to INTEGER*8 on most modern processors
        CHARACTER(LEN = 6), PARAMETER :: group_pixels = "pixels"
        CHARACTER(LEN = 14), PARAMETER :: dataset_bin1_id = "bin1_id"
        CHARACTER(LEN = 14), PARAMETER :: dataset_bin2_id = "bin2_id"
        CHARACTER(LEN = 14), PARAMETER :: dataset_bin3_id = "bin3_id"
        CHARACTER(LEN = 12), PARAMETER :: dataset_count = "count"
        CHARACTER(LEN = 4), PARAMETER :: group_bins = "bins"
        CHARACTER(LEN = 5), PARAMETER :: dataset_chorm = "chrom"
        CHARACTER(LEN = 5), PARAMETER :: dataset_start = "start"
        CHARACTER(LEN = 3), PARAMETER :: dataset_end = "end"
        CHARACTER(LEN = 6), PARAMETER :: group_chroms = "chroms"
        CHARACTER(LEN = 4), PARAMETER :: dataset_name = "name"
        CHARACTER(LEN = 6), PARAMETER :: dataset_length = "length"
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: grp_pixels_id
        INTEGER(HID_T) :: did_bin1_id
        INTEGER(HID_T) :: did_bin2_id
        INTEGER(HID_T) :: did_bin3_id
        INTEGER(HID_T) :: did_count
        INTEGER(HID_T) :: grp_bins_id
        INTEGER(HID_T) :: did_chrom
        INTEGER(HID_T) :: did_start
        INTEGER(HID_T) :: did_end
        INTEGER(HID_T) :: grp_chroms_id
        INTEGER(HID_T) :: did_name
        INTEGER(HID_T) :: did_length
        INTEGER(HID_T) :: dspace_bin1_id
        INTEGER(HID_T) :: dspace_bin2_id
        INTEGER(HID_T) :: dspace_bin3_id
        INTEGER(HID_T) :: dspace_count
        INTEGER(HID_T) :: dspace_chrom
        INTEGER(HID_T) :: dspace_start
        INTEGER(HID_T) :: dspace_end
        INTEGER(HID_T) :: dspace_name
        INTEGER(HID_T) :: dspace_length
        INTEGER :: rank_1 = 1               ! Dataset rank
        INTEGER :: error ! Error flag
        INTEGER :: bins_size, hic3d_len, I
        INTEGER(HSIZE_T), DIMENSION(1) :: pixels_dims, bins_dim, chroms_dim = 1
        logical, dimension(:), allocatable :: mask
        integer, dimension(:), allocatable :: bin1_id, bin2_id, bin3_id
        integer, dimension(:), allocatable :: count  ! TODO normalize it real: / Niter
        integer, dimension(:), allocatable :: chrom, start, end
        CHARACTER(LEN = 4), dimension(1) :: name = "chrS"
        integer, dimension(1) :: length

        ! HDF5 support
        ! Initialize FORTRAN interface.
        CALL h5open_f(rc)

        hic3d_shape = shape(hic3d)
        hic3d_len = hic3d_shape(1)

        call log%info('Calling save_hic3d_to_hdf5: analyse_folder=' // trim(analyse_folder) // &
                ', Niter=' // trim(str0(Niter)) // ', Nmeas=' // trim(str0(Nmeas)) // ', factor=' // trim(str(factor)))
        call log%info('in save_hic3d_to_hdf5: hic3d_len=' // trim(str(hic3d_len)))
        call log%info('in save_hic3d_to_hdf5: params%Nchain=' // trim(str(params%Nchain)))

        allocate (mask(hic3d_len))

        ! mask = .false.
        mask = hic3d(:,4) /= 0  ! take only non zero count pixels
        ! base bind_ids to start from 0: -1
        bin1_id = pack(hic3d(:,1), mask) - 1
        pixels_dims = size(bin1_id)
        call log%info('in  save_hic3d_to_hdf5: pixels_dims=' // trim(str(size(bin1_id))))
        bin2_id = pack(hic3d(:,2), mask) - 1
        bin3_id = pack(hic3d(:,3), mask) - 1
        count = pack(hic3d(:,4), mask)

        deallocate (mask)

        resolution = sim_resolution * factor

        bins_size = params%Nchain / factor
        bins_dim = params%Nchain / factor

        allocate(chrom(bins_size))
        allocate(start(bins_size))
        allocate(end(bins_size))

        chrom = 0  ! we have only one chromosome
        start = [0, (I, I = resolution, (bins_size-1)*resolution, resolution)]
        end = [resolution, (I, I = 2*resolution, bins_size*resolution, resolution)]
        length = bins_size*resolution

        ! output in HDF5 file
        filename_h5 = trim(analyse_folder) // 'hic3d_' // trim(str0(Nmeas)) // '_' // &
                trim(str(resolution/1000)) // 'k_cool.hdf5'

        ! Create a new file
        CALL h5fcreate_f(filename_h5, H5F_ACC_TRUNC_F, file_id, rc)

        CALL h5gcreate_f(file_id, group_pixels, grp_pixels_id, rc)
        CALL h5gcreate_f(file_id, group_bins, grp_bins_id, rc)
        CALL h5gcreate_f(file_id, group_chroms, grp_chroms_id, rc)

        ! Create the dataspace.
        !
        CALL h5screate_simple_f(rank_1, pixels_dims, dspace_bin1_id, error)
        CALL h5screate_simple_f(rank_1, pixels_dims, dspace_bin2_id, error)
        CALL h5screate_simple_f(rank_1, pixels_dims, dspace_bin3_id, error)
        CALL h5screate_simple_f(rank_1, pixels_dims, dspace_count, error)

        CALL h5screate_simple_f(rank_1, pixels_dims, dspace_count, error)

        CALL h5screate_simple_f(rank_1, bins_dim, dspace_chrom, error)
        CALL h5screate_simple_f(rank_1, bins_dim, dspace_start, error)
        CALL h5screate_simple_f(rank_1, bins_dim, dspace_end, error)

        CALL h5screate_simple_f(rank_1, chroms_dim, dspace_name, error)
        CALL h5screate_simple_f(rank_1, chroms_dim, dspace_length, error)

        ! Create the pixels group datasets
        CALL h5dcreate_f(grp_pixels_id, dataset_bin1_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_bin1_id, did_bin1_id, error)
        CALL h5dcreate_f(grp_pixels_id, dataset_bin2_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_bin2_id, did_bin2_id, error)
        CALL h5dcreate_f(grp_pixels_id, dataset_bin3_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_bin3_id, did_bin3_id, error)
        CALL h5dcreate_f(grp_pixels_id, dataset_count, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_count, did_count, error)

        CALL h5dcreate_f(grp_bins_id, dataset_chorm, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_chrom, did_chrom, error)
        CALL h5dcreate_f(grp_bins_id, dataset_start, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_start, did_start, error)
        CALL h5dcreate_f(grp_bins_id, dataset_end, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
                dspace_end, did_end, error)

!        CALL h5dcreate_f(grp_chroms_id, dataset_name, h5kind_to_type(SELECTED_CHAR_KIND, H5_CSET_UTF8), &
!                dspace_name, did_name, error)
!        CALL h5dcreate_f(grp_chroms_id, dataset_length, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), &
!                dspace_length, did_length, error)

        ! Write data
        CALL h5dwrite_f(did_bin1_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), bin1_id, pixels_dims, error)
        CALL h5dwrite_f(did_bin2_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), bin2_id, pixels_dims, error)
        CALL h5dwrite_f(did_bin3_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), bin3_id, pixels_dims, error)
        CALL h5dwrite_f(did_count, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), count, pixels_dims, error)

        CALL h5dwrite_f(did_chrom, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), chrom, bins_dim, error)
        CALL h5dwrite_f(did_start, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), start, bins_dim, error)
        CALL h5dwrite_f(did_end, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), end, bins_dim, error)

!        CALL h5dwrite_f(did_name, h5kind_to_type(SELECTED_CHAR_KIND, H5_CSET_UTF8), name, chroms_dim, error)
!        CALL h5dwrite_f(did_length, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), length, chroms_dim, error)

        ! Close the dataset.
        CALL h5dclose_f(did_bin1_id, error)
        CALL h5dclose_f(did_bin2_id, error)
        CALL h5dclose_f(did_bin3_id, error)
        CALL h5dclose_f(did_count, error)

        CALL h5dclose_f(did_chrom, error)
        CALL h5dclose_f(did_start, error)
        CALL h5dclose_f(did_end, error)

        ! Terminate access to the data space.
        CALL h5sclose_f(dspace_bin1_id, error)
        CALL h5sclose_f(dspace_bin2_id, error)
        CALL h5sclose_f(dspace_bin3_id, error)
        CALL h5sclose_f(dspace_count, error)

        CALL h5sclose_f(dspace_chrom, error)
        CALL h5sclose_f(dspace_start, error)
        CALL h5sclose_f(dspace_end, error)

        CALL h5gclose_f(grp_pixels_id, error)
        CALL h5gclose_f(grp_bins_id, error)
        CALL h5gclose_f(grp_chroms_id, error)

        ! Close the file.
        CALL h5fclose_f(file_id, error)

        call log%info('in save_hic3d_to_hdf5: Saved ' // trim(filename_h5))

        ! Close FORTRAN interface.
        CALL h5close_f(error)
    end subroutine save_hic3d_to_hdf5

    subroutine save_chip(chip, params, Nmeas, analyse_folder)
        implicit none
        class (ModelParameters), intent(in) :: params
        real :: chip(Nmeas, params%Nchain)
        integer, intent(in) :: Nmeas
        character(*), intent(in) :: analyse_folder

        integer :: rc, i

        open(10, file = trim(analyse_folder) // 'Chip.out', action = 'write', status = 'replace', iostat = rc)
        do i = 1, Nmeas
            write(10, *) chip(i, :)
        end do
        close(10)

        call log%info('Saved ' // trim(analyse_folder) // 'Chip.out')
    end subroutine save_chip

end module analyse_mod