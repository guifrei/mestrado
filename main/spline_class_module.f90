module spline_class_module
    use netlib_module
    implicit none
    private

    type, public :: spline_class
        private
        integer :: field_name
        integer :: nmax
        double precision, dimension(:), allocatable :: t, c
        integer :: nn
        integer :: ord
        double precision, dimension(:), allocatable :: wrkder
    contains
        procedure :: method_name
        procedure :: eval
        procedure :: eval_d
        final :: destructor
    end type spline_class

    interface spline_class
        module procedure new_spline_class
    end interface

contains

    function new_spline_class(vx, vy, nmax) result(r)
        type(spline_class) :: r
        integer, intent(in) :: nmax
        double precision, dimension(:), intent(in) :: vx, vy
        integer :: iopt, lwrk, ier
        double precision, parameter :: a = 0.04
        double precision, dimension(nmax) :: weight
        double precision :: xb, xe
        double precision :: s, fpf
        integer :: nest
        double precision, dimension(:), allocatable :: wrk
        integer, dimension(:), allocatable :: iwrk

        r%ord = 3
        r%nmax = nmax
        nest = nmax+r%ord+1
        iopt = 0
        lwrk = nmax*(r%ord+1)+nest*(7+3*r%ord)
        weight = 1.0
        s = 0.0
        xb = vx(1)
        xe = vx(nmax)

        allocate(r%t(nest), r%c(nest))
        allocate(r%wrkder(nmax))
        allocate(wrk(lwrk), iwrk(nest))

        call curfit(iopt, nmax, vx, vy, weight, xb, xe, r%ord, s, nest, r%nn, r%t, r%c, fpf, wrk, lwrk, iwrk, ier)

        deallocate(iwrk, wrk)
    end function

    subroutine method_name(this)
        class(spline_class), intent(in) :: this
    end subroutine
    
    function eval(this, x) result(y)
        class(spline_class), intent(in) :: this
        double precision, intent(in) :: x
        double precision :: y
        double precision, dimension(1) :: argx, argy
        integer :: ier

        argx(1) = x
        call splev(this%t, this%nn, this%c, this%ord, argx, argy, 1, ier)
        y = argy(1)
    end function

    function eval_d(this, x) result(y)
        class(spline_class), intent(inout) :: this
        double precision, intent(in) :: x
        double precision :: y
        double precision, dimension(1) :: argx, argy
        integer :: ier

        argx(1) = x
        call splder(this%t, this%nn, this%c, this%ord, 1, argx, argy, 1, this%wrkder, ier)
        y = argy(1)
    end function

    subroutine destructor(this)
        type(spline_class), intent(inout) :: this

        deallocate(this%t, this%c, this%wrkder)
    end subroutine

end module spline_class_module
