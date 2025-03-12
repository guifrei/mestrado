program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use estimated_h_module
    use tikhonov_module
    implicit none

    double precision, dimension(tnmax, 0: N), target :: m_fluxo_calor, m_delta_temperatura
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax) :: vx, vy
    double precision :: ymax
    character(len = 2) :: str_stdev, str_N
    integer :: nmax, stdev_idx, k, j
    double precision :: desv, x, y, y1, y2, dx, stdev
    double precision :: fluxo_calor_teorico, delta_temperatura_teorico, norm_f, norm_t
    integer :: kmax
    double precision :: start, finish

    dx = a/dble(tnmax - 1)

    open(unit = 1, file = '/home/cx3d/mestrado/data/coordinates.dat')
    do k = 1, tnmax
        x = dble(k - 1)*dx
        if (k == 1) then
            x = x + 0.01*dx
        else
            x = x - 0.01*dx
        end if
        write(1, *)x
    end do
    close(1)

    !===>
    call cpu_time(start)
    call calculate_reciprocity_coefficients(w1, dw1)
    call cpu_time(finish)
    !        write(*, *)'Elapsed time = ', (finish - start), ' s'
    !===>

    ! Salvando o perfil de temperatura calculado no Fortran
    open(unit = 1, file = '/home/cx3d/mestrado/data/fortran/temperaturas_sinteticas.dat')
    do k = 1, tnmax
        x = dble(k - 1)*dx
        y = t1(x, b)
        write(1, *)x, y
    end do
    close(1)

    !    Geracao do arquivo de comparacao de temperaturas medidas
    open(unit = 1, file = '/home/cx3d/mestrado/data/comsol/temperaturas_sinteticas.dat')
    open(unit = 2, file = '/home/cx3d/mestrado/data/fortran/temperaturas_sinteticas.dat')
    open(unit = 3, file = '/home/cx3d/mestrado/data/desvio_relativo.dat')
    do k = 1, tnmax
        read(1, *)x, y1
        read(2, *)x, y2
        desv = 100.0*dabs((y2 - y1)/y1)
        write(3, *)x, desv
    end do
    close(1)
    close(2)
    close(3)

    do stdev_idx = 0, 2
        if (stdev_idx == 0) then
            str_stdev = '00'
            stdev = 0.0
        else if (stdev_idx == 1) then
            str_stdev = '01'
            stdev = 0.1!ymax * 0.1/100.0
        else
            str_stdev = '05'
            stdev = 0.5!ymax * 0.5/100.0
        end if
        !                write(*, *)'        Stdev = ', stdev

        ! Recuperando as temperaturas do COMSOL
        ! e obtendo o valor absoluto maximo
        open(unit = 1, file = '/home/cx3d/mestrado/data/comsol/temperaturas_sinteticas.dat')
        do k = 1, tnmax
            read(1, *)vx(k), vy(k)
            if (k == 1) then
                ymax = dabs(vy(1))
            else if (dabs(vy(k)) > ymax) then
                ymax = dabs(vy(k))
            end if
        end do
        close(1)

        call add_error(vy, stdev)

        m_fluxo_calor = 0.0
        m_delta_temperatura = 0.0


        !call integrate_synthetic_temperatures(vx, vy, tnmax)
        call least_squares_for_Y(vx, vy, stdev_idx)

        open(unit = 1, file = '/home/cx3d/mestrado/data/temperaturas_sinteticas_stdev_'// str_stdev // '.dat')
        do k = 1, tnmax
            write(1, *)vx(k), vy(k)
        end do
        close(1)

        call cpu_time(start)
        do j = 0, N
            reciprocity_f(j) = calc_reciprocity_f(j)
            reciprocity_g(j) = calc_reciprocity_g(j)
        end do
        call cpu_time(finish)
        !                write(*, *)'    Elapsed time = ', (finish - start)*1000.0, ' ms'

        open(unit = 10, file = '/home/cx3d/mestrado/data/erro_rms__stdev_'// str_stdev // '.dat')
        do nmax = 0, N
            write(str_N, '(I2.2)') nmax
            open(unit = 4, file = '/home/cx3d/mestrado/data/fortran/delta_temperatura_stdev_' &
                // str_stdev // '_N_' // str_N // '.dat')
            open(unit = 5, file = '/home/cx3d/mestrado/data/fortran/fluxo_calor_interface_stdev_' &
                // str_stdev // '_N_' // str_N // '.dat')
            open(unit = 14, file = '/home/cx3d/mestrado/data/comsol/delta_temperatura.dat')
            open(unit = 15, file = '/home/cx3d/mestrado/data/comsol/fluxo_calor.dat')

            norm_f = 0.0
            norm_t = 0.0
            do j = 1, tnmax
                x = vx(j)
                if (nmax == 0) then
                    c_fluxo_calor = parcela_fluxo_calor(x, nmax, w1, dw1)
                    c_delta_temperatura = parcela_delta_temperatura(x, nmax, w1, dw1)
                else
                    c_fluxo_calor = m_fluxo_calor(j, nmax - 1) +&
                        parcela_fluxo_calor(x, nmax, w1, dw1)
                    c_delta_temperatura = m_delta_temperatura(j, nmax - 1) +&
                        parcela_delta_temperatura(x, nmax, w1, dw1)
                end if
                m_fluxo_calor(j, nmax) = c_fluxo_calor
                m_delta_temperatura(j, nmax) = c_delta_temperatura
                write(4, *)x, m_delta_temperatura(j, nmax)
                write(5, *)x, m_fluxo_calor(j, nmax)

                read(14, *)x, delta_temperatura_teorico
                read(15, *)x, fluxo_calor_teorico
                norm_t = norm_t + (delta_temperatura_teorico - m_delta_temperatura(j, nmax))**2
                norm_f = norm_f + (fluxo_calor_teorico - m_fluxo_calor(j, nmax))**2
            end do
            norm_t = sqrt(norm_t/tnmax)
            norm_f = sqrt(norm_f/tnmax)
            write(10, *)nmax, norm_t, norm_f
            close(15)
            close(14)
            close(5)
            close(4)
        end do
        close(10)

        !Principio da discrepancia de Morozov
        kmax = N
        if (stdev_idx /= 0) call morozov(stdev, vx, vy, vvY, kmax)

        open(unit = 10, file = '/home/cx3d/mestrado/data/erro_rms_stdev_'// str_stdev // '_morozov.dat')
        open(unit = 4, file = '/home/cx3d/mestrado/data/fortran/delta_temperatura_stdev_' &
            // str_stdev // '_morozov.dat')
        open(unit = 5, file = '/home/cx3d/mestrado/' // &
            'data/fortran/fluxo_calor_stdev_' // str_stdev // '_morozov.dat')
        open(unit = 14, file = '/home/cx3d/mestrado/data/comsol/delta_temperatura.dat')
        open(unit = 15, file = '/home/cx3d/mestrado/data/comsol/fluxo_calor.dat')
        open(unit = 7, file = '/home/cx3d/mestrado/data/estimativa_ctc__stdev_' // str_stdev // '_morozov.dat')

        norm_f = 0.0
        norm_t = 0.0
        do j = 1, tnmax
            x = vx(j)
            c_fluxo_calor = fluxo_calor(x, w1, dw1, kmax)
            c_delta_temperatura = delta_temperatura(x, w1, dw1, kmax)
            write(4, *)x, c_delta_temperatura
            write(5, *)x, c_fluxo_calor
            write(7, *)x, c_fluxo_calor/c_delta_temperatura

            read(14, *)x, delta_temperatura_teorico
            read(15, *)x, fluxo_calor_teorico
            norm_t = norm_t + (delta_temperatura_teorico - c_delta_temperatura)**2
            norm_f = norm_f + (fluxo_calor_teorico - c_fluxo_calor)**2
        end do
        norm_t = sqrt(norm_t/tnmax)
        norm_f = sqrt(norm_f/tnmax)
        write(10, *)nmax, norm_t, norm_f
        close(15)
        close(14)
        close(5)
        close(4)
        close(7)
        close(10)
    end do
end program main
