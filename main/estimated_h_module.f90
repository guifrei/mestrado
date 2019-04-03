module estimated_h_module
    use constants_module
    use netlib_module
    use reciprocity_functions_module
    implicit none

    integer :: nn
    integer, parameter :: ord = 5
    integer, parameter :: nest = tnmax+ord+1
    double precision, dimension((tnmax*(ord+1)+nest*(7+3*ord))) :: wrk
    double precision, dimension(nest) :: t, c
    integer, dimension(nest) :: iwrk

contains
    subroutine spline_interpolation_hest(vx, hy)
        double precision, dimension(tnmax), intent(in) :: vx, hy
        integer :: iopt, lwrk, ier, j
        double precision :: s, fpf
        double precision, dimension(tnmax) :: w

        iopt = 0
        lwrk = tnmax*(ord+1)+nest*(7+3*ord)
        w = 1.0
        s = 0.0

        do j = 1, ord+1
            t(j)=0.0
            t(j+ord+1)=a
        end do

        call curfit(iopt, tnmax, vx, hy, w, 0.0D0, a, ord, s, nest, nn, t, c, fpf, wrk, lwrk, iwrk, ier)

        if (ier > 0) then
            write(*, *)'Error in curfit. Ier = ', ier
            stop
        end if
    end subroutine

    function hest(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        double precision, dimension(1) :: argx, argy
        integer :: ier

        argx(1) = x

        call splev(t, nn, c, ord, argx, argy, 1, ier)
        if (ier /= 0) then
            write(*, *)'Error in curfit. Ier = ', ier
            stop
        end if
        r = argy(1)
    end function
end module estimated_h_module
