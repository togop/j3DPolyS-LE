!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Short Documentation:
!
! This module provies the Object 'Timer'
! It has two methods:
! - Tic (Timer) starts countint time
! - Tac (Timer) stops clock and prints time in standard output
!
! JJ 2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Timers
    use Kinds

    implicit none
    private
    public :: Timer

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  OBJECT TIMER
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type Timer
        private
        integer(KINT4) :: start, rate=-1
    contains
        procedure, public :: Tic, Tac
    end type Timer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  CREATES AND INITIALISES A TIMER
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Tic(self)
        class (Timer), intent(inout) :: self
        integer(KINT4) :: start, rate

        call system_clock(count_rate=rate)
        call system_clock(start)
        self%start=start
        self%rate=rate
    end subroutine Tic
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  FINISH THE COUNT AND PRINTS THE OUTPUT
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Tac(self, info)
        class (Timer), intent(in) :: self
        character(len=*), intent(in) :: info
        character(100) :: Tac
        integer(KINT4) :: finish

        if(self%rate<0) then
            print*, 'Call to ''Tac'' subroutine must come after call to ''Tic'''
            stop
        endif
        call system_clock(finish)
        write(Tac,'(a,f10.3)') info // ' Elapsed time in seconds:', float(finish-self%start)/self%rate
    end function Tac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Timers