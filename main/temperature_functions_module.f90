module temperature_functions_module
    use constants_module
    use eigenfunctions_module
    use netlib_module
    use interfaces_module
    use conductances_module
    implicit none

    double precision, dimension(0: 2*mmax_T+1), target :: vst
contains
    function eta(m, x, w, dw) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x

        interface
            function w(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
            function dw(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        double precision :: r
        double precision :: y, dy, v2, v3

        y = w(x)
        dy = dw(x)
        v2 = sinh(mu(m) * (b - y)) * cos(mu(m) * x) / cosh(mu(m) * b)
        v3 = dy * cosh(mu(m) * (b - y)) * sin(mu(m) * x) / cosh(mu(m) * b)
        r = v2 - v3
    end function


    function sigma(m, x, w, dw) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        double precision :: r

        interface
            function w(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
            function dw(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        double precision :: y, dy, v2, v3

        y = w(x)
        dy = dw(x)
        v2 = cosh(mu(m) * y) * cos(mu(m) * x) / cosh(mu(m) * b)
        v3 = dy * sinh(mu(m) * y) * sin(mu(m) * x) / cosh(mu(m) * b)
        r = v2 + v3
    end function


    function rho(m, x, w, dw, hc) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        interface
            function hc(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        double precision :: r

        interface
            function w(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
            function dw(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        r = sinh(mu(m) * w(x)) * cos(mu(m) * x) * hc(x) * sqrt(1.0 + dw(x) ** 2) / cosh(mu(m) * b)
    end function


    function kappa(m, x, w, dw, hc) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        interface
            function hc(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        double precision :: r

        interface
            function w(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
            function dw(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        r = cosh(mu(m) * (b - w(x))) * cos(mu(m) * x) * hc(x) * sqrt(1.0 + dw(x) ** 2) / cosh(mu(m) * b)
    end function

    function soma_controle_erro(parcela, x, y, v) result(r)
        interface
            function parcela(m, x, y, va) result(r)
                import
                integer, intent(in) :: m
                double precision, intent(in) :: x, y, va
                double precision :: r
            end function
        end interface
        double precision, intent(in) :: x, y
        double precision, dimension(0:), intent(in) :: v
        double precision :: r, partial_r, eps, r_acc
        integer :: m, p
        logical :: keep_1, keep_2, converged

        eps = 0.0
        r_acc = parcela(0, x, y, v(0))
        m = 1
        keep_1 = .true.
        do while (keep_1)
            r_acc = r_acc + parcela(m, x, y, v(m))
            m = m + 1
            keep_1 = m <= mmax_T
        end do

        do while (keep_1)
            p = m
            partial_r = 0.0
            keep_2 = .true.
            do while (keep_2)
                partial_r = partial_r + parcela(m, x, y, v(m))
                m = m + 1
                keep_1 = m <= mmax_T
                keep_2 = keep_1 .and. (m < p + delta_m)
            end do
            r_acc = r_acc + partial_r
            eps = dabs(partial_r/r_acc)
            converged = eps .lt. reltol
            keep_1 = keep_1 .and. (.not. converged)
        end do
        r = r_acc
    end function

    function parcela_t1(j, x, y, vc) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vc
        double precision :: r

        if (j == 0) then
            r = vc / a - q * y / k1
        else
            r = (2.0/a)*vc*(cosh(mu(j)*(b - y))/cosh(mu(j)*b))*cos(mu(j)*x)
        end if
    end function

    function parcela_d_t1_dx(j, x, y, vc) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vc
        double precision :: r

        if (j == 0) then
            r = 0.0
        else
            r = - mu(j)*(2.0/a)*vc*(cosh(mu(j)*(b - y))/cosh(mu(j)*b))*sin(mu(j)*x)
        end if
    end function

    function parcela_d_t1_dy(j, x, y, vc) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vc
        double precision :: r

        if (j == 0) then
            r = - q / k1
        else
            r = - mu(j)*(2.0/a)*vc*(sinh(mu(j)*(b - y))/cosh(mu(j)*b))*cos(mu(j)*x)
        end if
    end function

    function parcela_t2(j, x, y, vb) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vb
        double precision :: r

        if (j == 0) then
            r = vb * y / a
        else
            r = (2.0/a)*vb*(sinh(mu(j)*y)/cosh(mu(j)*b))*cos(mu(j)*x)
        end if
    end function

    function parcela_d_t2_dx(j, x, y, vb) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vb
        double precision :: r

        if (j == 0) then
            r = 0.0
        else
            r = - mu(j)*(2.0/a)*vb*(sinh(mu(j)*y)/cosh(mu(j)*b))*sin(mu(j)*x)
        end if
    end function

    function parcela_d_t2_dy(j, x, y, vb) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x, y, vb
        double precision :: r

        if (j == 0) then
            r = vb / a
        else
            r = mu(j)*(2.0/a)*vb*(cosh(mu(j)*y)/cosh(mu(j)*b))*cos(mu(j)*x)
        end if
    end function

    function t1(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_t1, x, y, vst(1::2))
    end function

    function d_t1_dx(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_d_t1_dx, x, y, vst(1::2))
    end function

    function d_t1_dy(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_d_t1_dy, x, y, vst(1::2))
    end function

    function t2(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_t2, x, y, vst(0::2))
    end function

    function d_t2_dx(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_d_t2_dx, x, y, vst(0::2))
    end function

    function d_t2_dy(x, y) result(r)
        double precision, intent(in) :: x, y
        double precision :: r

        r = soma_controle_erro(parcela_d_t2_dy, x, y, vst(0::2))
    end function

    ! Determinação dos coeficientes via transformação integral
    subroutine calculate_temperature_coefficients(w, dw, hc)
        interface
            function hc(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        integer, target :: j, n
        double precision, dimension(0: 2*mmax_T+1, 0: 2*mmax_T+1) :: mx
        double precision, dimension(0: 2*mmax_T+1) :: vz
        double precision :: x, dx
        double precision, dimension(2*mmax_T+2, 2*mmax_T+2) :: af
        integer, dimension(2*mmax_T+2) :: ipiv
        double precision, dimension(0: 2*mmax_T + 1) :: r, c
        double precision, dimension(1) :: ferr, berr
        double precision, dimension(8*mmax_T+8) :: work
        integer, dimension(2*mmax_T+2) :: iwork
        character(1) :: equed
        integer :: info
        double precision :: rcond

        interface
            function w(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
            function dw(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface

        ! Geracao das temperaturas
!        write(*, *)'Assembling system...'
        !$OMP PARALLEL DO COLLAPSE(2)
        do n = 0, mmax_T
            do j = 0, mmax_T
                mx(n * 2, j * 2) = transform(fa, n, c_loc(j), pts)
                mx(n * 2, j * 2 + 1) = transform(fb, n, c_loc(j), pts)
                mx(n * 2 + 1, j * 2) = transform(fc, n, c_loc(j), pts)
                mx(n * 2 + 1, j * 2 + 1) = transform(fd, n, c_loc(j), pts)
            end do
            vz(n * 2) = transform(fu, n, c_null_ptr, pts)
            vz(n * 2 + 1) = transform(fv, n, c_null_ptr, pts)
        end do
        !$OMP END PARALLEL DO

        equed = 'B'

!        write(*, *)'Invoking LU...'
        call dgesvx('E', 'N', 2*mmax_T+2, 1, mx, 2*mmax_T+2, af, 2*mmax_T+2, ipiv, equed, &
            r, c, vz, 2*mmax_T+2, vst, 2*mmax_T+2, rcond, ferr, berr, work, iwork, info)

    contains
        function fa(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            integer, pointer :: j

            call c_f_pointer(args, j)
            if (j == 0) then
                r = -hc(x) * sqrt(1.0 + dw(x) ** 2)*w(x)
            else
                r = -2.0*rho(j, x, w, dw, hc)
            end if
        end function

        function fb(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            integer, pointer :: j

            call c_f_pointer(args, j)
            if (j == 0) then
                r = hc(x) * sqrt(1.0 + dw(x) ** 2)
            else
                r = 2.0*(k1*mu(j)*eta(j, x, w, dw) + kappa(j, x, w, dw, hc))
            end if
        end function

        function fc(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            integer, pointer :: j

            call c_f_pointer(args, j)
            if (j == 0) then
                r = -(k2 + hc(x) * sqrt(1.0 + dw(x) ** 2)*w(x))
            else
                r = -2.0*(k2*mu(j)*sigma(j, x, w, dw) + rho(j, x, w, dw, hc))
            end if
        end function

        function fd(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            integer, pointer :: j

            call c_f_pointer(args, j)
            if (j == 0) then
                r = hc(x) * sqrt(1.0 + dw(x) ** 2)
            else
                r = 2.0*kappa(j, x, w, dw, hc)
            end if
        end function

        function fu(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r

            r = q*a*(hc(x) * sqrt(1.0 + dw(x) ** 2)*w(x)/k1 - 1.0)
        end function

        function fv(x, args) result(r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r

            r = q*a*hc(x) * sqrt(1.0 + dw(x) ** 2)*w(x)/k1
        end function
    end subroutine
end module temperature_functions_module
