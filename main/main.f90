program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use estimated_h_module
    use iso_c_binding
    use spline_class_module
    implicit none

    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax) :: vx, vy, vy_noerr, tmpy
    integer :: interface_idx, condutance_idx, stdev_idx, nmax, k, j
    double precision :: x, dx, sqrt_rms
    double precision :: h_est
    integer :: kmax
    double precision, dimension(3) :: stdev_values = [0.0D0, 0.1D0, 0.5D0]
    double precision :: stdev, zz
    character(len = 2) :: str_idx, str_cdx, str_stdev
    double precision, dimension(0: mmax_phi) :: werrY
    type(spline_class) :: spline

    wlist(1) = c_funloc(w1)
    wlist(2) = c_funloc(w2)
    wlist(3) = c_funloc(w3)
    wlist(4) = c_funloc(w4)
    wlist(5) = c_funloc(w5)
    wlist(6) = c_funloc(w6)
    wlist(7) = c_funloc(w7)
    wlist(8) = c_funloc(w8)
    wlist(9) = c_funloc(w9)
    wlist(10) = c_funloc(w10)
    wlist(11) = c_funloc(w11)
    wlist(12) = c_funloc(w12)
    wlist(13) = c_funloc(w13)
    wlist(14) = c_funloc(w14)
    wlist(15) = c_funloc(w15)

    dwlist(1) = c_funloc(dw1)
    dwlist(2) = c_funloc(dw2)
    dwlist(3) = c_funloc(dw3)
    dwlist(4) = c_funloc(dw4)
    dwlist(5) = c_funloc(dw5)
    dwlist(6) = c_funloc(dw6)
    dwlist(7) = c_funloc(dw7)
    dwlist(8) = c_funloc(dw8)
    dwlist(9) = c_funloc(dw9)
    dwlist(10) = c_funloc(dw10)
    dwlist(11) = c_funloc(dw11)
    dwlist(12) = c_funloc(dw12)
    dwlist(13) = c_funloc(dw13)
    dwlist(14) = c_funloc(dw14)
    dwlist(15) = c_funloc(dw15)

    hlist(1) = c_funloc(h1)
    hlist(2) = c_funloc(h2)
    hlist(3) = c_funloc(h3)
    hlist(4) = c_funloc(h4)
    hlist(5) = c_funloc(h5)
    hlist(6) = c_funloc(h6)
    hlist(7) = c_funloc(h7)
    hlist(8) = c_funloc(h8)
    hlist(9) = c_funloc(h9)

    call init_uncertainty()

    call c_f_procpointer(wlist(1), w)

    dx = a/1000.0
    open(unit=7,file='../paper/interface.dat')
    open(unit=8,file='../paper/conductance.dat')
    do j = 0, 1000
        x = dx*dble(j)
        zz = 0
        do k = 0, mmax_phi
            zz = zz + werrY(k)*cos(mu(k)*x)
        end do
        write(7, *)x, w1(x), w2(x), w3(x), werr(x), w(x) + spline%eval(x)
        write(8, *)x, h1(x), h2(x), h3(x)
    end do
    close(7)
    close(8)

    do interface_idx = 1, 3
        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)
        do condutance_idx = 1, 3
            call c_f_procpointer(hlist(condutance_idx), h)
            call calculate_temperature_coefficients(h, w, dw, vx, vy_noerr)
            do stdev_idx = 1, 3
                write(str_idx, '(I2.2)') interface_idx
                write(str_cdx, '(I2.2)') condutance_idx
                if (stdev_idx.eq.1) then
                    str_stdev = '00'
                else if (stdev_idx.eq.2) then
                    str_stdev = '01'
                else
                    str_stdev = '05'
                end if

                write(*, *)'Interface = ', interface_idx, ', conductance = ', condutance_idx, ', stdev = ', str_stdev

                stdev = stdev_values(stdev_idx)
                vy = vy_noerr
                call add_error(vy, stdev)
                call least_squares_for_Y(vx, vy, vvY, tnmax)

                call filter_frequencies(vx, vy, vvY, tnmax, 4)!TODO decidir como implementar este filtro

                open(unit=10,file='../paper/difference_interface_' // &
                    str_idx // '_conductance_' // str_cdx // '_stdev_' // str_stdev // '.dat')
                if (stdev.ne.0) then
                    tmpy = 0.0
                    do j = 0, mmax_phi
                        tmpy = tmpy + vvY(j)*cos(mu(j)*vx)
                        sqrt_rms = norm2(tmpy - vy)
                        write(10, *)j, sqrt_rms
                    end do
                end if
                close(10)

                call generate_results(h, w, dw, w, dw)
            end do
        end do
    end do

contains
    subroutine generate_results(h, w, dw, werr, dwerr)
        interface
            function h(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function

            function w(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function

            function dw(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function

            function werr(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function

            function dwerr(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        call calculate_reciprocity_coefficients(werr, dwerr)

        kmax = N
        nmax = kmax

        do j = 0, nmax
            reciprocity_f(j) = calc_reciprocity_f(j)
            reciprocity_g(j) = calc_reciprocity_g(j)
        end do

        open(7, file='../paper/calculated_ctc_interface_' // &
            str_idx // '_conductance_' // str_cdx // '_stdev_' // str_stdev // '.dat')

        dx = a/dble(tmax - 1)
        do j = 1, tmax
            x = dx*(j - 1)
            c_fluxo_calor = 0
            c_delta_temperatura = 0
            do k = 0, nmax
                c_fluxo_calor = c_fluxo_calor + parcela_fluxo_calor(x, k, w, dw)
                c_delta_temperatura = c_delta_temperatura + parcela_delta_temperatura(x, k, w, dw)
            end do
            h_est = c_fluxo_calor/c_delta_temperatura
            write(7, *)x, h_est, h(x)
        end do
        close(7)
    end subroutine

    subroutine filter_frequencies(vx, vy, vvY, nmax, idx)
        integer, intent(in) :: nmax
        double precision, dimension(nmax), intent(in) :: vx
        double precision, dimension(nmax), intent(inout) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        integer, intent(in) :: idx
        integer :: i, j
        double precision, dimension(:, :), allocatable :: mxa
        double precision :: x

        allocate(mxa(nmax, mmax_phi + 1))
        vvY(idx:mmax_phi) = 0.0
        mxa = 0.0
        do i = 1, nmax
            x = vx(i)
            do j = 1, mmax_phi + 1
                mxa(i, j) = cos(mu(j - 1)*x)
            end do
        end do

        do j = 1, mmax_phi + 1
            mxa(:, j) = cos(mu(j - 1)*vx)
        end do

        vy = matmul(mxa, vvY)
        deallocate(mxa)
    end subroutine

    subroutine init_uncertainty()
        integer, parameter :: nmax = 100
        double precision, dimension(nmax) :: vx, rand_u, rand_v, rand_epsilon
        double precision :: stdev
        double precision :: mean
        integer :: i

        stdev = b/50.0
        dx = a/dble(nmax - 1)
        do j = 1, nmax
            vx(j) = dx*(j - 1)
        end do

        call random_seed
        call random_number(rand_u)
        call random_number(rand_v)
        rand_epsilon = cos(2.0*pi*rand_v)*sqrt(-2.0*log(rand_u))
        rand_epsilon = rand_epsilon*stdev

        stdev = 0
        mean = 0
        do i = 1,100
            mean = mean + rand_epsilon(i)
            stdev = stdev + rand_epsilon(i)*rand_epsilon(i)
        end do

        mean = mean/100.0
        stdev = stdev/10000.0

        stdev = sqrt(stdev - mean*mean)

        call least_squares_for_Y(vx, rand_epsilon, werrY, nmax)

        spline = spline_class(vx, rand_epsilon, nmax)
    end subroutine

    function werr(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        integer :: i

        r = w(x)
        do i = 0, mmax_phi
            r = r + werrY(i)*cos(mu(i)*x)
        end do
    end function

    function dwerr(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        integer :: i

        r = dw(x)
        do i = 0, mmax_phi
            r = r - mu(i)*werrY(i)*sin(mu(i)*x)
        end do
    end function
end program main
