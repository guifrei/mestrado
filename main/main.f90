program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use algebraic_reconstruction_technique_module
    use estimated_h_module
    implicit none

    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision, dimension(tnmax, 0: N), target :: m_fluxo_calor, m_delta_temperatura
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax, 0: (N + 1)**2 - 1) :: m_condutancia_contato
    double precision, dimension(tnmax) :: vx, vy, tcalc
    double precision :: norm
    character(len = 2) :: str_idx, str_cdx, str_stdev, str_N, str_n_fluxo_calor, str_n_delta_temperatura
    integer :: interface_idx, condutance_idx, nmax, stdev_idx, k, nmax1, nmax2, m, j, n_fluxo_calor, n_delta_temperatura
    double precision :: desv, x, y, y1, y2, dx, stdev, arg, norm_acc, y1prev, y2prev
    integer, target, dimension(2) :: f_args


    !    block
    !        integer, parameter :: nn = 2
    !        integer, parameter :: mm = 3
    !        integer, parameter :: sz = max(nn, mm)
    !        integer, parameter :: minnm = min(nn, mm)
    !        double precision :: aa(sz, mm)
    !        double precision :: ainv(sz, nn)
    !        double precision, dimension(minnm) :: ss, ee
    !        double precision, dimension(sz, nn) :: uu
    !        double precision, dimension(sz, mm) :: vv
    !        double precision, dimension(nn) :: wwork
    !        integer :: irank, ierr
    !
    !        aa(1:2, 1) = [5.0, 3.0]
    !        aa(1:2, 2) = [-2.0, 4.0]
    !        aa(1:2, 3) = [1.0, -1.0]
    !
    !        call mpinv(sz, nn, mm, minnm, aa, ainv, ss, ee, uu, vv, wwork, irank, ierr)
    !        write(*, *)ainv(1, :)
    !        write(*, *)ainv(2, :)
    !        write(*, *)ainv(3, :)
    !        stop
    !    end block

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

    hlist(1) = c_funloc(h1)
    hlist(2) = c_funloc(h2)
    hlist(3) = c_funloc(h3)
    hlist(4) = c_funloc(h4)
    hlist(5) = c_funloc(h5)
    hlist(6) = c_funloc(h6)
    hlist(7) = c_funloc(h7)
    hlist(8) = c_funloc(h8)
    hlist(9) = c_funloc(h9)

    dx = a/dble(tmax - 1)

    ! Geracao dos arquivos de interfaces
    do interface_idx = 1, 3
        write(str_idx, '(I2.2)') interface_idx
        call c_f_procpointer(wlist(interface_idx), w)
        open(unit = 1, file = '/home/cx3d/mestrado/' // &
            'data/interface_' // str_idx // '.dat')
        do k = 0, tmax - 1
            x = dble(k)*dx
            write(1, *)x, w(x)
        end do
        close(1)
    end do

    ! Geracao dos arquivos de condutancias de contato
    do condutance_idx = 1, 3
        write(str_cdx, '(I2.2)') condutance_idx
        call c_f_procpointer(hlist(condutance_idx), h)
        open(unit = 1, file = '/home/cx3d/mestrado/' // &
            'data/conductance_' // str_cdx // '.dat')
        do k = 0, tmax - 1
            x = dble(k)*dx
            write(1, *)x, h(x)
        end do
        close(1)
    end do

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

    do interface_idx = 1, 3
        write(*, *)'Interface = ', interface_idx
        write(str_idx, '(I2.2)') interface_idx
        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)
        call calculate_reciprocity_coefficients(interface_idx)

        open(unit = 1, file='/home/cx3d/mestrado/' // &
            'data/media_beta_gamma_interface_'//str_idx//'.dat')
        do j = 0, N
            f_args(1) = j
            f_args(2) = interface_idx
            write(1, *)j, integrate(f_aux_beta, c_loc(f_args), pts)/a, integrate(f_aux_gamma, c_loc(f_args), pts)/a
        end do
        close(1)

        do condutance_idx = 1, 3
            write(*, *)'    Condutance = ', interface_idx
            write(str_cdx, '(I2.2)') condutance_idx
            call c_f_procpointer(hlist(condutance_idx), h)
            call calculate_temperature_coefficients(interface_idx, condutance_idx, h)

            ! Salvando o perfil de temperatura calculado no Fortran
            open(unit = 1, file = '/home/cx3d/mestrado/' // &
                'data/fortran/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            do k = 1, tnmax
                x = dble(k - 1)*dx
                y = t1(x, b)
                write(1, *)x, y
            end do
            close(1)

            !    Geracao do arquivo de comparacao de temperaturas medidas
            open(unit = 1, file = '/home/cx3d/mestrado/' // &
                'data/comsol/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            open(unit = 2, file = '/home/cx3d/mestrado/' // &
                'data/fortran/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
            open(unit = 3, file = '/home/cx3d/mestrado/' // &
                'data/desvio_relativo_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
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
                    stdev = 0.1
                else
                    str_stdev = '05'
                    stdev = 0.3882
                end if
                !                write(*, *)'        Stdev = ', stdev
                write(*, *)'Interface = ', interface_idx, ', conductance = ', condutance_idx, ', stdev = ', stdev

                ! Recuperando as temperaturas do COMSOL
                open(unit = 1, file = '/home/cx3d/mestrado/' // &
                    'data/fortran/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
                do k = 1, tnmax
                    read(1, *)vx(k), vy(k)
                end do
                close(1)

                call add_error(vy, stdev)
                call regularize(vx, vy, tnmax)

                m_fluxo_calor = 0.0
                m_delta_temperatura = 0.0


                !call integrate_synthetic_temperatures(vx, vy, tnmax)
                call least_squares_for_Y(vx, vy)

                do j = 0, N
                    reciprocity_f(j) = calc_reciprocity_f(j)
                    reciprocity_g(j) = calc_reciprocity_g(j)
                end do

                do nmax = 0, N
                    write(str_N, '(I2.2)') nmax
                    open(unit = 4, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_'//str_cdx// &
                        '_stdev_' // str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 5, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_'//str_cdx// &
                        '_stdev_' // str_stdev // '_N_' // str_N // '.dat')
                    do j = 1, tnmax
                        x = vx(j)
                        if (nmax == 0) then
                            c_fluxo_calor = parcela_fluxo_calor(x, nmax, interface_idx)
                            c_delta_temperatura = parcela_delta_temperatura(x, nmax, interface_idx)
                        else
                            c_fluxo_calor = m_fluxo_calor(j, nmax - 1) +&
                                parcela_fluxo_calor(x, nmax, interface_idx)
                            c_delta_temperatura = m_delta_temperatura(j, nmax - 1) +&
                                parcela_delta_temperatura(x, nmax, interface_idx)
                        end if
                        m_fluxo_calor(j, nmax) = c_fluxo_calor
                        m_delta_temperatura(j, nmax) = c_delta_temperatura
                        write(4, *)x, m_delta_temperatura(j, nmax)
                        write(5, *)x, m_fluxo_calor(j, nmax)
                    end do
                    close(5)
                    close(4)
                end do

                y1 = dabs(reciprocity_f(0))
                y2 = dabs(reciprocity_g(0))
                open(unit=1, file='/home/cx3d/mestrado/' // &
                    'data/metricas_interface_'//str_idx // &
                    '_stdev_' // str_stdev // '.dat')
                do j = 1, N
                    y1prev = y1
                    y2prev = y2
                    y1 = y1 + dabs(reciprocity_f(j))
                    y2 = y2 + dabs(reciprocity_g(j))
                    write(1, *)j, y1/y1prev, y2/y2prev
                end do
                close(1)
            end do
        end do
    end do

    n_delta_temperatura = 6
    n_fluxo_calor = 4
    str_idx = '02'
    str_cdx = '02'
    str_stdev = '05'

    write(str_n_delta_temperatura, '(I2.2)') n_delta_temperatura
    write(str_n_fluxo_calor, '(I2.2)') n_fluxo_calor

    open(unit = 4, file = '/home/cx3d/mestrado/' // &
        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_'//str_cdx// &
        '_stdev_' // str_stdev // '_N_' // str_n_delta_temperatura // '.dat')
    open(unit = 5, file = '/home/cx3d/mestrado/' // &
        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_'//str_cdx// &
        '_stdev_' // str_stdev // '_N_' // str_n_fluxo_calor // '.dat')

    open(unit = 7, file = '/home/cx3d/mestrado/' // &
        'data/estimativa_ctc_interface_'//str_idx//'_conductance_'//str_cdx// &
        '_stdev_' // str_stdev // '.dat')

    do j = 1, tnmax
        read(4, *)x, c_delta_temperatura
        read(5, *)x, c_fluxo_calor
        write(7, *)x, c_fluxo_calor/c_delta_temperatura
    end do

    close(7)
    close(4)
    close(5)

contains

    function f_aux_beta(x, args) result (r)
        double precision, intent(in) :: x
        type(c_ptr), intent(in) :: args
        double precision :: r
        integer, dimension(:), pointer :: ptr

        call c_f_pointer(args, ptr, [2])
        r = fbeta(ptr(1), x, ptr(2))

    end function

    function f_aux_gamma(x, args) result (r)
        double precision, intent(in) :: x
        type(c_ptr), intent(in) :: args
        double precision :: r

        integer, dimension(:), pointer :: ptr

        call c_f_pointer(args, ptr, [2])
        r = fgamma(ptr(1), x, ptr(2))
    end function
end program main
