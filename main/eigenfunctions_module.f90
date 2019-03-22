module eigenfunctions_module
    use iso_c_binding
    use constants_module
    implicit none

    integer, private :: i
    double precision, dimension(0: mmax_phi) :: mu = (/(dble(i)*pi/a, i=0, mmax_phi)/)
contains
    function transform(f, m, args, pts) result (r)
        interface
            function f(x, args) result (r)
                import
                double precision, intent(in) :: x
                type(c_ptr), intent(in) :: args
                double precision :: r
            end function
        end interface
        integer, intent(in) :: m
        type(c_ptr), intent(in), optional :: args
        double precision, dimension(:), intent(in), optional, target :: pts
        double precision :: r

        double precision, dimension(:), allocatable :: def_pts
        type(c_ptr) :: def_args
        integer :: npts2, cnt
        double precision :: abserr
        integer :: neval, ier, leniw, maxp1, lenw, last
        integer, dimension(:), allocatable :: iwork
        double precision, dimension(:), allocatable :: work
        double precision :: tmp

        if (present(pts)) then
            npts2 = size(pts)
            allocate(def_pts(npts2))
            def_pts = pts
        else
            npts2 = 2
            allocate(def_pts(npts2))
            def_pts(1) = 0.0
            def_pts(2) = a
        end if

        if (present(args)) then
            def_args = args
        else
            def_args = c_null_ptr
        end if

        leniw = 5000
        maxp1 = 1000
        lenw = leniw*2+maxp1*25
        allocate(iwork(leniw), work(lenw))

        r = 0.0D0
        do cnt = 2, npts2
            call dqawo(f_aux, def_pts(cnt-1), def_pts(cnt), mu(m), &
                1, 1.49D-8, 1.49D-8, tmp, abserr, neval, ier, leniw, maxp1, &
                lenw, last, iwork, work)
            r = r + tmp
        end do

        deallocate(iwork, work)
        deallocate(def_pts)
    contains
        function f_aux(x) result(r)
            double precision, intent(in) :: x
            double precision :: r

            r = f(x, args)
        end function
    end function
end module eigenfunctions_module
