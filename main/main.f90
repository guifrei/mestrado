program main
    use iso_c_binding
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    implicit none

    integer :: interface_idx, condutance_idx, k, j, nmax, stdev_idx, n_delta_temperatura, n_fluxo_calor
    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision :: x, stdev, dx, desv, y1, y2, norma, dnorma, norma_acc, rms, norma_delta_temperatura, norma_fluxo_calor
    character(len = 2) :: str_idx, str_cdx, str_N, str_stdev, str_N_delta_temperatura, str_N_fluxo_calor

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

    do interface_idx = 1, 3
        write(str_idx, '(I2.2)') interface_idx
        do condutance_idx = 1, 3
            write(str_cdx, '(I2.2)') condutance_idx
            call calculate_temperature_coefficients(interface_idx, condutance_idx)
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
                write(*, *)'Interface = ', interface_idx, ', conductance = ', condutance_idx, ', stdev = ', stdev

                call c_f_procpointer(hlist(condutance_idx), h)
                call c_f_procpointer(wlist(interface_idx), w)
                call c_f_procpointer(dwlist(interface_idx), dw)

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

                write(*, *)'Integrating y = b'
                call calculate_integrals_Y(interface_idx, condutance_idx, stdev_idx)
                call calculate_reciprocity_coefficients(interface_idx, condutance_idx)

                do nmax = 0, N
                    !write(*, *)'N = ', nmax
                    write(str_N, '(I2.2)') nmax
                    open(unit = 2, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 3, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    do j = 0, tnmax - 1
                        x = dble(j)*dx
                        if (j == 0) then
                            x = x + dx*0.01
                        else if (j == tnmax - 1) then
                            x = x - dx*0.01
                        end if
                        write(2, *)x, delta_temperatura(x, interface_idx, nmax)
                        write(3, *)x, fluxo_calor(x, interface_idx, nmax)
                    end do
                    close(2)
                    close(3)
                end do

                ! geracao dos arquivos de erro RMS - delta temperatura
                open(unit = 4, file = '/home/cx3d/mestrado/' // &
                    'data/erro_rms_delta_temperatura_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                do nmax = 0, N
                    rms = 0.0
                    write(str_N, '(I2.2)') nmax
                    open(unit = 2, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 3, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/delta_temperatura_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '.dat')
                    do j = 0, tnmax - 1
                        read(2, *)x, y1
                        read(3, *)x, y2
                        rms = rms + (y1 - y2)**2
                    end do
                    rms = sqrt(rms/dble(tnmax))
                    write(4, *)nmax, log10(rms)
                    close(2)
                    close(3)
                end do
                close(4)

                ! geracao dos arquivos de erro RMS - fluxo de calor
                open(unit = 4, file = '/home/cx3d/mestrado/' // &
                    'data/erro_rms_fluxo_calor_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                do nmax = 0, N
                    write(str_N, '(I2.2)') nmax
                    rms = 0.0
                    open(unit = 2, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 3, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/fluxo_calor_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '.dat')
                    do j = 0, tnmax - 1
                        read(2, *)x, y1
                        read(3, *)x, y2
                        rms = rms + (y1 - y2)**2
                    end do
                    rms = sqrt(rms/dble(tnmax))
                    write(4, *)nmax, log10(rms)
                    close(2)
                    close(3)
                end do
                close(4)

                ! Geracao dos arquivos de normas - delta temperatura
                norma_acc = 0.0
                norma_delta_temperatura = 0.0
                open(unit = 4, file = '/home/cx3d/mestrado/' // &
                    'data/norma_delta_temperatura_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                open(unit = 5, file = '/home/cx3d/mestrado/' // &
                    'data/norma_acumulada_delta_temperatura_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                do nmax = 0, N
                    write(str_N, '(I2.2)') nmax
                    norma = 0.0
                    open(unit = 2, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 3, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/delta_temperatura_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '.dat')
                    do j = 0, tnmax - 1
                        read(2, *)x, y1
                        read(3, *)x, y2
                        dnorma = y1 - y2
                        norma = norma + dnorma**2
                    end do
                    norma_acc = norma_acc + norma
                    close(2)
                    close(3)
                    norma = sqrt(norma)
                    write(4, *)nmax, log10(norma)
                    if (nmax == 0) then
                        norma_delta_temperatura = norma
                        n_delta_temperatura = nmax
                    else if (norma < norma_delta_temperatura) then
                        norma_delta_temperatura = norma
                        n_delta_temperatura = nmax
                    end if
                    if (nmax >= 1) write(5, *)nmax, log10(norma_acc)
                end do
                close(4)
                close(5)

                ! Geracao dos arquivos de normas - fluxo de calor
                norma_acc = 0.0
                norma_fluxo_calor = 0.0
                open(unit = 4, file = '/home/cx3d/mestrado/' // &
                    'data/norma_fluxo_calor_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                open(unit = 5, file = '/home/cx3d/mestrado/' // &
                    'data/norma_acumulada_fluxo_calor_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                do nmax = 0, N
                    write(str_N, '(I2.2)') nmax
                    norma = 0.0
                    open(unit = 2, file = '/home/cx3d/mestrado/' // &
                        'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N // '.dat')
                    open(unit = 3, file = '/home/cx3d/mestrado/' // &
                        'data/comsol/fluxo_calor_interface_'//str_idx//'_conductance_' &
                        //str_cdx // '.dat')
                    do j = 0, tnmax - 1
                        read(2, *)x, y1
                        read(3, *)x, y2
                        dnorma = y1 - y2
                        norma = norma + dnorma**2
                    end do
                    close(2)
                    close(3)
                    norma = sqrt(norma)
                    norma_acc = norma_acc + norma
                    write(4, *)nmax, log10(norma)
                    if (nmax == 0) then
                        norma_fluxo_calor = norma
                        n_fluxo_calor = nmax
                    else if (norma < norma_fluxo_calor) then
                        norma_fluxo_calor = norma
                        n_fluxo_calor = nmax
                    end if
                    if (nmax >= 1) write(5, *)nmax, log10(norma_acc)
                end do
                close(4)
                close(5)

                write(*, *)'Generating CTC estimatives'
                write(*, *)'n_delta_temperatura = ', n_delta_temperatura
                write(*, *)'n_fluxo_calor = ', n_fluxo_calor

                ! Geracao das estimativas de condutancia termica de contato
                write(str_N_delta_temperatura, '(I2.2)') n_delta_temperatura
                write(str_N_fluxo_calor, '(I2.2)') n_fluxo_calor
                open(unit = 2, file = '/home/cx3d/mestrado/' // &
                    'data/estimativa_ctc_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '.dat')
                open(unit = 3, file = '/home/cx3d/mestrado/' // &
                    'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N_delta_temperatura // '.dat')
                open(unit = 4, file = '/home/cx3d/mestrado/' // &
                    'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '_N_' // str_N_fluxo_calor // '.dat')

                open(unit = 5, file = '/home/cx3d/mestrado/' // &
                    'data/fortran/delta_temperatura_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '_melhor.dat')
                open(unit = 7, file = '/home/cx3d/mestrado/' // &
                    'data/fortran/fluxo_calor_interface_'//str_idx//'_conductance_' &
                    //str_cdx // '_stdev_'// str_stdev // '_melhor.dat')

                do j = 0, tnmax - 1
                    read(3, *)x, y1
                    read(4, *)x, y2
                    write(2, *)x, y2/y1
                    write(5, *)x, y1
                    write(7, *)x, y2
                end do

                close(7)
                close(5)
                close(4)
                close(3)
                close(2)
            end do
        end do
    end do

end program main
