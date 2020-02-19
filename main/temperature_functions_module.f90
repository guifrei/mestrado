module temperature_functions_module
    use constants_module
    use eigenfunctions_module
    use netlib_module
    use interfaces_module
    use conductances_module
    implicit none

    double precision, dimension(0: 2*mmax_T+1), target :: vst
contains
    function eta(m, x, interface_idx) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision :: y, dy, v1, v2, v3

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        y = w(x)
        dy = dw(x)
        v1 = cosh(mu(m) * b)
        v2 = sinh(mu(m) * (b - y)) * cos(mu(m) * x) / v1
        v3 = dy * cosh(mu(m) * (b - y)) * sin(mu(m) * x) / v1
        r = v2 - v3
    end function


    function sigma(m, x, interface_idx) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision :: y, dy, v1, v2, v3

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        y = w(x)
        dy = dw(x)
        v1 = cosh(mu(m) * b)
        v2 = cosh(mu(m) * y) * cos(mu(m) * x) / v1
        v3 = dy * sinh(mu(m) * y) * sin(mu(m) * x) / v1
        r = v2 + v3
    end function


    function rho(m, x, interface_idx, hc) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx
        interface
            function hc(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        r = sinh(mu(m) * w(x)) * cos(mu(m) * x) * hc(x) * sqrt(1.0 + dw(x) ** 2) / cosh(mu(m) * b)
    end function


    function kappa(m, x, interface_idx, hc) result(r)
        integer, intent(in) :: m
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx
        interface
            function hc(x) result(r)
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

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
    subroutine calculate_temperature_coefficients(interface_idx, condutance_idx, hc, vx, vy, write_files)
        integer, intent(in) :: interface_idx, condutance_idx
        interface
            function hc(x) result(r)
                import
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        logical, intent(in), optional :: write_files
        double precision, dimension(0: tnmax - 1), optional, target :: vx, vy
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        integer, target :: j, n
        double precision, dimension(0: 2*mmax_T+1, 0: 2*mmax_T+1) :: mx
        double precision, dimension(0: 2*mmax_T+1) :: vz
        double precision :: x, dx
        character(len = 2) :: str_idx, str_cdx
        double precision, dimension(2*mmax_T+2, 2*mmax_T+2) :: af
        integer, dimension(2*mmax_T+2) :: ipiv
        double precision, dimension(0: 2*mmax_T + 1) :: r, c
        double precision, dimension(1) :: ferr, berr
        double precision, dimension(8*mmax_T+8) :: work
        integer, dimension(2*mmax_T+2) :: iwork
        character(1) :: equed
        integer :: info
        double precision :: rcond
        logical :: opt_write_files
        double precision, dimension(:), pointer :: opt_vx, opt_vy

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        ! Geracao das temperaturas
        !        write(*, *)'Assembling system...'
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

        equed = 'B'

        !        write(*, *)'Invoking LU...'
        call dgesvx('E', 'N', 2*mmax_T+2, 1, mx, 2*mmax_T+2, af, 2*mmax_T+2, ipiv, equed, &
            r, c, vz, 2*mmax_T+2, vst, 2*mmax_T+2, rcond, ferr, berr, work, iwork, info)

        if (present(vy)) then
            opt_vx => vx
            opt_vy => vy
            dx = a/dble(tnmax - 1)
            do n = 0, tnmax - 1
                x = dble(n)*dx
                opt_vx(n) = x
                opt_vy(n) = t1(x, b)
            end do
        end if

        if (present(write_files)) then
            opt_write_files = write_files
        else
            opt_write_files = .false.
        end if

        if (opt_write_files) then
            write(str_idx, '(I2.2)') interface_idx
            write(str_cdx, '(I2.2)') condutance_idx

            open(unit = 1, file = '/home/cx3d/mestrado/' // &
                'data/fortran/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            open(unit = 3, file = '/home/cx3d/mestrado/' // &
                'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            open(unit = 4, file = '/home/cx3d/mestrado/' // &
                'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            dx = a/dble(tnmax - 1)
            do n = 0, tnmax - 1
                x = dble(n)*dx
                if (n == 0) then
                    x = x + dx*0.01
                else if (n == tnmax - 1) then
                    x = x - dx*0.01
                end if
                write(1, *)x, t1(x, b)
                write(3, *)x, t1(x, w(x)) - t2(x, w(x))
                write(4, *)x, -k1*(dw(x)*d_t1_dx(x, w(x)) - d_t1_dy(x, w(x)))/sqrt(1.0 + dw(x)**2)
            end do
            close(4)
            close(2)
            close(1)
        end if

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
                r = -2.0*rho(j, x, interface_idx, hc)
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
                r = 2.0*(k1*mu(j)*eta(j, x, interface_idx) + kappa(j, x, interface_idx, hc))
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
                r = -2.0*(k2*mu(j)*sigma(j, x, interface_idx) + rho(j, x, interface_idx, hc))
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
                r = 2.0*kappa(j, x, interface_idx, hc)
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
