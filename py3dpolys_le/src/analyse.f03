module analyse_mod
    use lattice_data_mod
    use PolymerModel_mod
    use logging_mod
    !    USE ISO_C_BINDING
    USE HDF5
    implicit none

    private
    public :: str, analyse
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

    subroutine analyse(radiuscontact, use_contact_probability, &
            params, Niter, Nmeas, output_folder, analyse_folder, hic3d_factor, chrom)

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
        character(*), intent(in) :: chrom

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

        print*, 'analyse output_folder=' // trim(output_folder) // ', analyse_folder=' // trim(analyse_folder)

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
        do j = 1, Nmeas
            close(40+j)
        end do

        Chip = Chip / real(Niter)
        hic = hic / real(Niter)

        call save_hic_to_hdf5(hic, params, Nmeas, analyse_folder, chrom)
        call save_chip(Chip, params, Nmeas, analyse_folder, chrom)

        if (do_hic3d) then
            ! to Normlize myaybe do this
            ! hic3d(:,4) = hic3d(:,4) / real(Niter * hic3d_factor**2)
            call save_hic3d_to_hdf5(hic3d, hic3d_factor, params, Niter, Nmeas, analyse_folder, chrom)
        end if
    end subroutine analyse

    subroutine save_hic_to_hdf5(hic, params, Nmeas, analyse_folder, chrom)
        implicit none
        class (ModelParameters), intent(in) :: params
        real, intent(in), TARGET :: hic(Nmeas, params%Nchain, params%Nchain)    ! or: hic(*)
        integer, intent(in) :: Nmeas
        character(*), intent(in) :: analyse_folder
        character(*), intent(in) :: chrom
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

            CALL h5_add_attr_str(file_id, "chromosome", chrom)
            CALL h5_add_attr_str(file_id, "format", "HDF5:hic_matrix")
            CALL h5_add_attr_str(file_id, "format-url", "https://gitlab.com/togop/3DPolyS-LE")
            CALL h5_add_attr_str(file_id, "format-version", "1")
            CALL h5_add_attr_str(file_id, "generated", "3DPolyS-LEv2022.7")

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

    subroutine h5_add_attr_int(dest_id, name, value)
        INTEGER(HID_T) :: dest_id
        character(*), intent(in) :: name
        integer, intent(in) :: value
        INTEGER, PARAMETER :: int_kind_8 = SELECTED_INT_KIND(9)
        INTEGER :: rank = 1
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsf = 1
        INTEGER(HID_T) :: dspace_attr
        INTEGER(HID_T) :: attr_id
        INTEGER :: status

        CALL h5screate_simple_f(rank, dimsf, dspace_attr, status)
        CALL h5screate_f(H5S_SCALAR_F, dspace_attr, status)
        CALL h5acreate_f(dest_id, name, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), dspace_attr, attr_id, status)
        CALL h5awrite_f(attr_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), value, dimsf, status)

        CALL h5sclose_f(dspace_attr, status)
        !CALL h5sclose_f(attr_id, status)

    end subroutine h5_add_attr_int

    subroutine h5_add_attr_str(dest_id, name, value)
        INTEGER(HID_T) :: dest_id
        character(*), intent(in) :: name
        character(*), intent(in) :: value
        INTEGER, PARAMETER :: int_kind_8 = SELECTED_INT_KIND(9)
        !INTEGER :: rank = 1
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsf = 1
        INTEGER(HID_T) :: dspace_attr
        INTEGER(HID_T) :: attr_id, attr_type
        INTEGER :: status

        INTEGER          , PARAMETER :: dim0      = 1
        INTEGER(SIZE_T) :: str_len
        INTEGER(HID_T)  :: memtype

        str_len = len(value)
        ! sdim = str_len
        CALL H5Tcopy_f(H5T_C_S1, attr_type, status)
        CALL H5Tset_size_f(attr_type, str_len+1, status)

        CALL H5Tcopy_f( H5T_FORTRAN_S1, memtype, status)
        CALL H5Tset_size_f(memtype, str_len, status)

        CALL H5Screate_f(H5S_NULL_F, dspace_attr, status)

        CALL h5screate_simple_f(1, dimsf, dspace_attr, status)
        CALL H5Acreate_f(dest_id, name, attr_type, dspace_attr, attr_id, status)

        CALL h5awrite_f(attr_id, memtype, value, dimsf, status)

        CALL h5aclose_f(attr_id, status)
        CALL h5sclose_f(dspace_attr, status)
        CALL H5tclose_f(attr_type, status)
        CALL H5tclose_f(memtype, status)

    end subroutine h5_add_attr_str

    subroutine h5_add_dataset_int(dest_id, name, dataset)
        INTEGER(HID_T) :: dest_id
        character(*), intent(in) :: name
        integer, dimension(:), intent(in) :: dataset
        INTEGER, PARAMETER :: int_kind_8 = SELECTED_INT_KIND(9)
        INTEGER(HSIZE_T), DIMENSION(1) :: dims
        INTEGER(HID_T) :: dspace_id
        INTEGER(HID_T) :: dataset_id
        INTEGER :: status

        dims = size(dataset)

        CALL h5screate_simple_f(1, dims, dspace_id, status)
        CALL h5dcreate_f(dest_id, name, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), dspace_id, dataset_id, status)
        CALL h5dwrite_f(dataset_id, h5kind_to_type(int_kind_8, H5_INTEGER_KIND), dataset, dims, status)

        CALL h5dclose_f(dataset_id, status)
        CALL h5sclose_f(dspace_id, status)
    end subroutine h5_add_dataset_int

    subroutine h5_add_dataset_real(dest_id, name, dataset)
        INTEGER(HID_T) :: dest_id
        character(*), intent(in) :: name
        real, dimension(:), intent(in) :: dataset
        INTEGER, PARAMETER :: real_kind_7 = SELECTED_REAL_KIND(6, 37) !should map to REAL*4 on most modern processors
        INTEGER(HSIZE_T), DIMENSION(1) :: dims
        INTEGER(HID_T) :: dspace_id
        INTEGER(HID_T) :: dataset_id
        INTEGER :: status

        dims = size(dataset)

        CALL h5screate_simple_f(1, dims, dspace_id, status)
        CALL h5dcreate_f(dest_id, name, h5kind_to_type(real_kind_7,H5_REAL_KIND), dspace_id, dataset_id, status)
        CALL h5dwrite_f(dataset_id, h5kind_to_type(real_kind_7,H5_REAL_KIND), dataset, dims, status)

        CALL h5dclose_f(dataset_id, status)
        CALL h5sclose_f(dspace_id, status)
    end subroutine h5_add_dataset_real

    subroutine save_hic3d_to_hdf5(hic3d, factor, params, Niter, Nmeas, analyse_folder, chrom)
        ! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
        implicit none
        class (ModelParameters), intent(in) :: params
        integer, intent(in), TARGET :: hic3d(:, :)    ! cool like format
        integer, intent(in) :: factor, Niter, Nmeas
        character(*), intent(in) :: analyse_folder
        character(*), intent(in) :: chrom
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
        CHARACTER(LEN = 4), PARAMETER :: group_bins = "bins"
        CHARACTER(LEN = 6), PARAMETER :: group_chroms = "chroms"
        INTEGER(HID_T) :: file_id
        INTEGER(HID_T) :: grp_pixels_id
        INTEGER(HID_T) :: grp_bins_id
        INTEGER(HID_T) :: grp_chroms_id
        INTEGER :: error ! Error flag
        INTEGER :: bins_size, hic3d_len, I
        logical, dimension(:), allocatable :: mask
        integer, dimension(:), allocatable :: bin1_id, bin2_id, bin3_id
        real, dimension(:), allocatable :: count
        integer, dimension(:), allocatable :: chrom_id, start_pos, end_pos
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
        ! base bind_ids to start_pos from 0: -1
        bin1_id = pack(hic3d(:,1), mask) - 1
        call log%info('in  save_hic3d_to_hdf5: pixels_dims=' // trim(str(size(bin1_id))))
        bin2_id = pack(hic3d(:,2), mask)
        bin2_id = bin2_id - 1
        bin3_id = pack(hic3d(:,3), mask)
        bin3_id = bin3_id - 1
        count = real(pack(hic3d(:,4), mask)) / real(Niter)

        deallocate (mask)

        resolution = sim_resolution * factor

        bins_size = params%Nchain / factor

        allocate(chrom_id(bins_size))
        allocate(start_pos(bins_size))
        allocate(end_pos(bins_size))

        chrom_id = 0  ! we have only one chromosome
        start_pos = [0, (I, I = resolution, (bins_size-1)*resolution, resolution)]
        end_pos = [resolution, (I, I = 2*resolution, bins_size*resolution, resolution)]
        length(1) = bins_size*resolution

        ! output in HDF5 file
        filename_h5 = trim(analyse_folder) // 'hic3d_' // trim(str0(Nmeas)) // '_' // &
                trim(str(resolution/1000)) // 'k.cool3d'

        ! Create a new file
        CALL h5fcreate_f(filename_h5, H5F_ACC_TRUNC_F, file_id, rc)

        CALL h5gcreate_f(file_id, group_pixels, grp_pixels_id, rc)
        CALL h5gcreate_f(file_id, group_bins, grp_bins_id, rc)
        CALL h5gcreate_f(file_id, group_chroms, grp_chroms_id, rc)

        CALL h5_add_dataset_int(grp_pixels_id, "bin1_id", bin1_id)
        CALL h5_add_dataset_int(grp_pixels_id, "bin2_id", bin2_id)
        CALL h5_add_dataset_int(grp_pixels_id, "bin3_id", bin3_id)
        CALL h5_add_dataset_real(grp_pixels_id, "count", count)

        CALL h5_add_dataset_int(grp_bins_id, "chrom", chrom_id)
        CALL h5_add_dataset_int(grp_bins_id, "start", start_pos)
        CALL h5_add_dataset_int(grp_bins_id, "end", end_pos)

        CALL h5_add_dataset_int(grp_chroms_id, "length", length)
        !CALL h5_add_dataset_str(grp_chroms_id, "name", chrom)  ! TODO add chromosome name

        ! add root attributes
        CALL h5_add_attr_int(file_id, "bin-size", sim_resolution)
        CALL h5_add_attr_int(file_id, "nbins", bins_size)
        CALL h5_add_attr_str(file_id, "bin-type", "fixed")
        !CALL h5_add_attr(file_id, "nchroms", 1)
        CALL h5_add_attr_str(file_id, "format", "HDF5:Cooler3D")
        CALL h5_add_attr_str(file_id, "format-url", "https://gitlab.com/togop/3DPolyS-LE")
        CALL h5_add_attr_str(file_id, "format-version", "1")
        CALL h5_add_attr_str(file_id, "generated", "3DPolyS-LEv2022.7")

        CALL h5gclose_f(grp_pixels_id, error)
        CALL h5gclose_f(grp_bins_id, error)
        CALL h5gclose_f(grp_chroms_id, error)

        ! Close the file.
        CALL h5fclose_f(file_id, error)

        call log%info('in save_hic3d_to_hdf5: Saved ' // trim(filename_h5))

        ! Close FORTRAN interface.
        CALL h5close_f(error)
    end subroutine save_hic3d_to_hdf5

    subroutine save_chip(chip, params, Nmeas, analyse_folder, chrom)
        implicit none
        class (ModelParameters), intent(in) :: params
        real :: chip(Nmeas, params%Nchain)
        integer, intent(in) :: Nmeas
        character(*), intent(in) :: analyse_folder
        character(*), intent(in) :: chrom
        integer, parameter :: sim_resolution = 2000
        integer :: rc, m, i
        character*1 :: tab = char(9)
        integer, dimension(:), allocatable :: start_pos, end_pos

        allocate(start_pos(params%Nchain))
        allocate(end_pos(params%Nchain))
        start_pos =  [0, (i, i = sim_resolution, (params%Nchain-1)*sim_resolution, sim_resolution)]
        end_pos = [sim_resolution, (i, i = 2*sim_resolution, params%Nchain*sim_resolution, sim_resolution)]

        open(11, file = trim(analyse_folder) // 'chip_lef.out', action = 'write', status = 'replace', iostat = rc)
        do m = 1, Nmeas
            write(11, *) chip(m, :)  ! all in one file
            ! output in bedGraph format: chrom, start, end, value
            open(10, file = trim(analyse_folder) // 'chip_lef_' // trim(str0(m)) // '.bedGraph', &
                    action = 'write', status = 'replace', iostat = rc)
                do i = 1, params%Nchain
                    write(10, *) trim(chrom) // tab // trim(str(start_pos(i))) // tab // trim(str(end_pos(i))) &
                            // tab // trim(strff(chip(m, i), '(F8.5)'))
                end do
            close(10)
            call log%info('Saved ' // trim(analyse_folder) // 'chip_lef_' // trim(str0(Nmeas)) // '.bedGraph')
        end do
        call log%info('Saved ' // trim(analyse_folder) // 'chip_lef.out')
        close(11)

    end subroutine save_chip

end module analyse_mod