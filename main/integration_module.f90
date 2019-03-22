module integration_module
    use iso_c_binding
    use constants_module
    implicit none
contains
    function integrate(f, args, pts) result (r)
        interface
            function f(x, args) result (r)
                import
                double precision, intent(in) :: x
                type(c_ptr), intent(in) :: args
                double precision :: r
            end function
        end interface
        type(c_ptr), intent(in), optional :: args
        double precision, dimension(:), intent(in), optional, target :: pts
        double precision :: r

        double precision, dimension(:), allocatable :: function_pts
        type(c_ptr) :: function_args
        integer :: npts2
        double precision :: abserr
        integer :: neval, ier, lenw, last, limit, cnt
        integer, dimension(:), allocatable :: iwork
        double precision, dimension(:), allocatable :: work
        double precision :: tmp

        if (present(pts)) then
            npts2 = ubound(pts, 1) - lbound(pts, 1) + 1
            allocate(function_pts(npts2))
            function_pts = pts
        else
            npts2 = 2
            allocate(function_pts(npts2))
            function_pts(1) = 0.0D0
            function_pts(2) = a
        end if

        if (present(args)) then
            function_args = args
        else
            function_args = c_null_ptr
        end if

        limit = 1000
        lenw = 4*limit
        allocate(iwork(limit), work(lenw))

        r = 0
        do cnt = 2, npts2
            call dqag(f_aux, function_pts(cnt-1), function_pts(cnt), 1.49D-8, 1.49D-8, 6,&
                tmp, abserr, neval, ier, limit, lenw, last, iwork, work)
            !            call dqags(f_aux, function_pts(cnt-1), function_pts(cnt), 1.49D-8, 1.49D-8,&
            !                tmp, abserr, neval, ier, limit, lenw, last, iwork, work)
            !            call dqng(f_aux, function_pts(cnt-1), function_pts(cnt), 1.49D-8, 1.49D-8, &
            !                tmp, abserr, neval, ier)
            r = r + tmp
        end do
        deallocate(iwork, work)

        !        leniw = 500!3*npts2-2
        !        lenw = leniw*2-npts2
        !        allocate(iwork(leniw), work(lenw))
        !        call dqagp(f_aux, 0.0D0, a, npts2, function_pts, 1.0D-8, 1.0D-8, r, abserr, neval, ier, leniw, lenw, last, iwork, work)
        !        deallocate(iwork, work)

        deallocate(function_pts)
    contains
        function f_aux(x) result(r)
            double precision, intent(in) :: x
            double precision :: r

            r = f(x, args)
        end function
    end function
end module integration_module
