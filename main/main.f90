!data00 = Import[
!   "/home/cx3d/mestrado/data/temperaturas_sinteticas_interface_03_\
!conductance_01_stdev_00.dat", {"Data", All, 2}];
!
!gaussData01 = GaussianFilter[data01, 20.0, Method -> "Gaussian"];
!ListLinePlot[{data00, gaussData01}]

!plot 'conductance_01.dat', 'estimativa_ctc_interface_01_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_01_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_01_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_03_stdev_05.dat'
!
!
!===============
!
!
!plot 'conductance_01.dat', 'estimativa_ctc_interface_02_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_02_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_02_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_03_stdev_05.dat'
!
!===============
!
!
!plot 'conductance_01.dat', 'estimativa_ctc_interface_03_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_03_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_03_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_03_stdev_05.dat'


program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use estimated_h_module
    use tikhonov_module
    implicit none

    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision, dimension(tnmax, 0: N), target :: m_fluxo_calor, m_delta_temperatura
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax, 0: (N + 1)**2 - 1) :: m_condutancia_contato
    double precision, dimension(tnmax) :: vx, vy, tcalc
    double precision :: norm, ymax
    character(len = 2) :: str_idx, str_cdx, str_stdev, str_N, str_n_fluxo_calor, str_n_delta_temperatura
    integer :: interface_idx, condutance_idx, nmax, stdev_idx, k, nmax1, nmax2, m, j, n_fluxo_calor, n_delta_temperatura
    double precision :: desv, x, y, y1, y2, dx, stdev, arg, norm_acc, y1prev, y2prev
    double precision :: fluxo_calor_teorico, delta_temperatura_teorico, norm_f, norm_t
    integer, target, dimension(2) :: f_args
    integer :: nmax_delta_temperatura, nmax_fluxo_calor
    double precision :: start, finish


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
    !        stop0.38
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

        !===>
        call cpu_time(start)
        call calculate_reciprocity_coefficients(interface_idx)
        call cpu_time(finish)
        !        write(*, *)'Elapsed time = ', (finish - start), ' s'
        !===>

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
                    stdev = 0.1!ymax * 0.1/100.0
                else
                    str_stdev = '05'
                    stdev = 0.5!ymax * 0.5/100.0
                end if
                !                write(*, *)'        Stdev = ', stdev
                write(*, *)'Interface = ', interface_idx, ', conductance = ', condutance_idx, ', stdev = ', stdev

                ! Recuperando as temperaturas do COMSOL
                ! e obtendo o valor absoluto maximo
                open(unit = 1, file = '/home/cx3d/mestrado/' // &
                    'data/comsol/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
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
                call least_squares_for_Y(vx, vy, interface_idx, condutance_idx, stdev_idx)
                if (stdev_idx /= 0) call morozov(stdev, vx, vy, vvY)

                open(unit = 1, file = '/home/cx3d/mestrado/' // &
                    'data/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx // &
                    '_stdev_'// str_stdev // '.dat')
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

                open(unit = 10, file = '/home/cx3d/mestrado/' // &
                    'data/erro_rms_interface_'//str_idx//'_conductance_'//str_cdx // &
                    '_stdev_'// str_stdev // '.dat')
                do nmax = 0, N
                    write(str_N, '(I2.2)') nmax
                    open(unit = 4, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_'//str_cdx// &
                        '_stdev_' // str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 5, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_'//str_cdx// &
                        '_stdev_' // str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 14, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/delta_temperatura_interface_'//str_idx//'_conductance_'//str_cdx// '.dat')
                    open(unit = 15, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/fluxo_calor_interface_'//str_idx//'_conductance_'//str_cdx// '.dat')

                    norm_f = 0.0
                    norm_t = 0.0
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

    !Impressao dos indices otimos

    do interface_idx = 1, 3
        write(str_idx, '(I2.2)') interface_idx
        do condutance_idx = 1, 3
            write(str_cdx, '(I2.2)') condutance_idx
            do stdev_idx = 0, 2
                if (stdev_idx == 0) then
                    str_stdev = '00'
                    stdev = 0.0
                else if (stdev_idx == 1) then
                    str_stdev = '01'
                    stdev = 0.1
                else
                    str_stdev = '05'
                    stdev = 0.5
                end if

                open(unit = 10, file = '/home/cx3d/mestrado/' // &
                    'data/erro_rms_interface_'//str_idx//'_conductance_'//str_cdx // &
                    '_stdev_'// str_stdev // '.dat')
                read(10, *)nmax, y1, y2
                nmax_delta_temperatura = nmax
                nmax_fluxo_calor = nmax
                do j = 1, N
                    read(10, *)nmax, norm_t, norm_f
                    if (norm_t < y1) then
                        nmax_delta_temperatura = nmax
                        y1 = norm_t
                    end if
                    if (norm_f < y2) then
                        nmax_fluxo_calor = nmax
                        y2 = norm_f
                    end if
                end do
                close(10)
                write(*, *)'DT: ', str_idx, ' ', str_cdx, ' ', str_stdev, ' ',nmax_delta_temperatura
                write(*, *)'DQ: ', str_idx, ' ', str_cdx, ' ', str_stdev, ' ',nmax_fluxo_calor

                write(str_n_delta_temperatura, '(I2.2)') nmax_delta_temperatura
                write(str_n_fluxo_calor, '(I2.2)') nmax_fluxo_calor

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
            end do
        end do
    end do
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
