module logging_mod
    implicit none
    private
    public :: Logger
    public :: LOG_LEVEL, LOG_OFF, LOG_FATAL, LOG_ERROR, LOG_WARN, LOG_INFO, LOG_DEBUG, LOG_TRACE
    integer, public :: global_log_level

    !common/logging/global_log_level

    enum, bind(c)
        enumerator :: LOG_LEVEL = 0
        enumerator :: LOG_OFF = 1
        enumerator :: LOG_FATAL = 2
        enumerator :: LOG_ERROR = 3
        enumerator :: LOG_WARN = 4
        enumerator :: LOG_INFO = 5
        enumerator :: LOG_DEBUG = 6
        enumerator :: LOG_TRACE = 7
    end enum

    type Logger
        private
        character(256), public :: source ! = ''
        integer(kind(LOG_LEVEL)), public :: level ! = LOG_INFO
    contains
        procedure, public :: log, error, warn, info, debug, set_level
    end type Logger

    public :: str, str0, strf, str2loglevel

contains

    character(len = 20) function str(k)
        !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

    character(len = 20) function str0(k)
        !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str0, '(I3.3)') k  ! must be in agreement with plot_hic.py, plot_chi2-min_time_hist.py: f'hic_{i:03}.hdf5'
        str0 = adjustl(str0)
    end function str0

    character(len = 20) function strf(f)
        !   "Convert an float to string."
        real, intent(in) :: f
        write (strf, *) f
        strf = adjustl(strf)
    end function strf

    subroutine set_level(self, level)
        class (Logger), intent(inout) :: self
        integer(kind(LOG_LEVEL)), intent(in) :: level
        self%level = level
    end subroutine set_level

    subroutine log(self, level, message)
        class (Logger), intent(inout) :: self
        integer(kind(LOG_LEVEL)), intent(in) :: level
        character(*), intent(in) :: message
        ! print*, 'my log level: ', self%level, ', message log level: ', level
        if (level <= self%level) then
            ! TODO add log to a file support
            ! TODO add format to output
            print*, '[',trim(level2str(level)),'] ',trim(log_datetime()),' ',trim(self%source),' - ', trim(message)
        end if
    end subroutine log

    subroutine error(self, message)
        class (Logger), intent(inout) :: self
        character(*), intent(in) :: message
        call self%log(LOG_ERROR, message)
    end subroutine error

    subroutine warn(self, message)
        class (Logger), intent(inout) :: self
        character(*), intent(in) :: message
        call self%log(LOG_WARN, message)
    end subroutine warn

    subroutine info(self, message)
        class (Logger), intent(inout) :: self
        character(*), intent(in) :: message
        call self%log(LOG_INFO, message)

    end subroutine info

    subroutine debug(self, message)
        class (Logger), intent(inout) :: self
        character(*), intent(in) :: message
        call self%log(LOG_DEBUG, message)
    end subroutine debug

    function log_datetime()
        character(100) :: log_datetime

        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer,dimension(8) :: values

        call date_and_time(date, time, zone, values)
        !call date_and_time(DATE=date,ZONE=zone)
        !call date_and_time(TIME=time)
        !call date_and_time(VALUES=values)

        write(log_datetime, '(i4.4, "-", i2.2, "-", i2.2, " ", i2.2, ":", i2.2, ":", i2.2," ")') & ! a,
                values(1), values(2), values(3), &
                values(5), values(6), values(7) !, zone
        !write(log_datetime, '(a,":",a,":",a," ")') time(1:2), time(3:4), time(5:6)
        !write(log_datetime, '(a,"/",a,"/",a," ")') date(1:4), date(5:6), date(7:8)
    end function log_datetime

    function level2str(level)  ! log priority
        integer(kind(LOG_LEVEL)), intent(in) :: level
        character(len=5) level2str
        character(len=5), dimension(7), parameter :: levels = ['OFF  ', &
                'FATAL', &
                'ERROR', &
                'WARN ', &
                'INFO ', &
                'DEBUG', &
                'TRACE']

        level2str = levels(level)
    end function level2str

    function str2loglevel(level_str)  ! log level as integer
        character(*), intent(in) :: level_str
        integer(kind(LOG_LEVEL)) :: str2loglevel
        integer :: i
        character(len=5), dimension(7), parameter :: levels = ['OFF  ', &
                'FATAL', &
                'ERROR', &
                'WARN ', &
                'INFO ', &
                'DEBUG', &
                'TRACE']
        str2loglevel = LOG_INFO
        do i = 1,7
            if (level_str == trim(levels(i))) then
                str2loglevel = i
                exit
            end if
        end do
    end function str2loglevel

end module logging_mod