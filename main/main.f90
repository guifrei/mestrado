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
    double precision, dimension(tnmax, szdiv, 0: N), target :: m_fluxo_calor, m_delta_temperatura
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax, 0: (N + 1)**2 - 1) :: m_condutancia_contato
    double precision, dimension(tnmax) :: vx, vy, tcalc
    double precision :: norm
    character(len = 2) :: str_idx, str_cdx, str_stdev, str_N
    integer :: interface_idx, condutance_idx, nmax, stdev_idx, j, k, nmax1, nmax2, m, nlim, ndiv_idx
    double precision :: desv, x, y, y1, y2, dx, stdev, arg, norm_acc, y1prev, y2prev

    double precision, pointer :: tmp


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

    do interface_idx = 2, 2 !2, 2
        write(*, *)'Interface = ', interface_idx
        write(str_idx, '(I2.2)') interface_idx
        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)
        call calculate_reciprocity_coefficients(interface_idx)

        do condutance_idx = 3, 3
            write(*, *)'    Condutance = ', interface_idx
            write(str_cdx, '(I2.2)') condutance_idx
            call c_f_procpointer(hlist(condutance_idx), h)
            call calculate_temperature_coefficients(interface_idx, condutance_idx, h)
            do stdev_idx = 0, 0
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

                !                !    Geracao do arquivo de comparacao de temperaturas medidas
                !                open(unit = 1, file = '/home/cx3d/mestrado/' // &
                !                    'data/comsol/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
                !                open(unit = 2, file = '/home/cx3d/mestrado/' // &
                !                    'data/fortran/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
                !                open(unit = 3, file = '/home/cx3d/mestrado/' // &
                !                    'data/desvio_relativo_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
                !                do k = 1, tnmax
                !                    read(1, *)x, y1
                !                    read(2, *)x, y2
                !                    desv = 100.0*dabs((y2 - y1)/y1)
                !                    write(3, *)x, desv
                !                    ! Armazenando as temperaturas sintéticas teóricas
                !                    vx(k) = x
                !                    vy(k) = y1
                !                end do
                !                close(1)
                !                close(2)
                !                close(3)

                call calculate_temperature_coefficients(interface_idx, condutance_idx, h, .false.)
                do k = 1, tnmax
                    x = dble(k - 1)*dx
                    y = t1(x, b)
                    vx(k) = x
                    vy(k) = y
                end do

                call add_error(vy, stdev)

                ! TODO
                nlim = tnmax/ndiv(szdiv) + mod(tnmax, ndiv(szdiv))
                call integrate_synthetic_temperatures(vx(1::ndiv(szdiv)), vy(1::ndiv(szdiv)), nlim)
                y1 = dabs(reciprocity_f(0))
                y2 = dabs(reciprocity_g(0))
                open(unit=100, file='/home/cx3d/a.txt')
                do j = 1, N
                    y1prev = y1
                    y2prev = y2
                    y1 = y1 + dabs(reciprocity_f(j))
                    y2 = y2 + dabs(reciprocity_g(j))
                    write(100, *)j, y1/y1prev, y2/y2prev
                end do
                close(100)
                ! end TODO

                tmp => m_fluxo_calor(1, 1, 0)

                m_fluxo_calor = 0.0
                m_delta_temperatura = 0.0

                do ndiv_idx = 1, szdiv
                    nlim = tnmax/ndiv(ndiv_idx) + mod(tnmax, ndiv(ndiv_idx))

                    call integrate_synthetic_temperatures(vx(1::ndiv(ndiv_idx)), vy(1::ndiv(ndiv_idx)), nlim)




                    if (N <= nlim) nlim = N

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
                                c_fluxo_calor = m_fluxo_calor(j, ndiv_idx, nmax - 1) +&
                                    parcela_fluxo_calor(x, nmax, interface_idx)
                                c_delta_temperatura = m_delta_temperatura(j, ndiv_idx, nmax - 1) +&
                                    parcela_delta_temperatura(x, nmax, interface_idx)
                            end if
                            m_fluxo_calor(j, ndiv_idx, nmax) = c_fluxo_calor
                            m_delta_temperatura(j, ndiv_idx, nmax) = c_delta_temperatura
                            write(4, *)x, m_delta_temperatura(j, ndiv_idx, nmax)
                            write(5, *)x, m_fluxo_calor(j, ndiv_idx, nmax)
                        end do
                        close(5)
                        close(4)
                    end do
                    write(*, *)'ndiv_idx = ', ndiv_idx
                end do
            end do
        end do
    end do
end program main
