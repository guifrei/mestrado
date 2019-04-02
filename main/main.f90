program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use algebraic_reconstruction_technique_module
    implicit none

    integer :: interface_idx, condutance_idx, k, j, nmax, stdev_idx, n_delta_temperatura, n_fluxo_calor
    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision :: x, stdev, dx, desv, y1, y2, norma, dnorma, norma_acc, rms, norma_delta_temperatura, norma_fluxo_calor
    character(len = 2) :: str_idx, str_cdx, str_N, str_N_prev, str_stdev, str_N_delta_temperatura, str_N_fluxo_calor
    integer, parameter :: ord = 5
    double precision, dimension(:), allocatable :: t, c
    integer :: nn

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

    ! thorsted m. buzug - computed tomography
!    block
!        integer, parameter :: m = 5
!        integer :: imax = 500
!        double precision, dimension(m, m):: mA
!        double precision, dimension(m):: mX
!        double precision, dimension(m):: mB
!        double precision :: lambda
!
!        mA(1, :) = [5, 3, 2, 0, 4]
!        mA(2, :) = [-1, 4, 6, 2, 0]
!        mA(3, :) = [0, 1, 3, -1, 1]
!        mA(4, :) = [2, 0, 0, -3, 4]
!        mA(5, :) = [5, 2, 0, -1, 0]
!
!        mX = [-1, 2, -3, 4, -5]
!
!
!        mB = [-25, -1, -16, -34, -5]
!        mX = 0
!
!        lambda = 1.8
!
!        call art(m, mA, mB, mX, lambda, imax)
!        write(*, *)mX
!
!
!    end block

!    stop

    block
        double precision, dimension(tnmax) :: vx, vy, hy
        integer :: iopt, nest, lwrk, nmax1, nmax2, ier
        double precision, dimension(tnmax) :: w
        double precision :: r1, r2, tcalc
        double precision :: s, fpf
        integer, dimension(:), allocatable :: iwrk
        double precision, dimension(:), allocatable :: wrk
        double precision :: norm

        do interface_idx = 3, 3
            call calculate_reciprocity_coefficients(interface_idx)
            do condutance_idx = 3, 3
                call c_f_procpointer(hlist(condutance_idx), h)
                do stdev_idx = 2, 2
                    write(str_idx, '(I2.2)') interface_idx
                    write(str_cdx, '(I2.2)') condutance_idx
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

                    call calculate_temperature_coefficients(interface_idx, condutance_idx, h)
                    call calculate_integrals_Y(interface_idx, condutance_idx, stdev_idx)

                    ! Recuperando temperaturas medidas
                    open(unit = 1, file = '/home/cx3d/mestrado/' // &
                        'data/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx // &
                        '_stdev_'// str_stdev // '.dat')
                    do k = 1, tnmax
                        read(1, *)vx(k), vy(k)
                    end do
                    close(1)

                    do nmax1 = 1, 6
                        do nmax2 = 1, 6
                            iopt = 0
                            nest=tnmax+ord+1
                            lwrk = tnmax*(ord+1)+nest*(7+3*ord)
                            w = 1.0

                            allocate(t(nest), c(nest), iwrk(nest))

                            do j = 1, tnmax
                                r1 = delta_temperatura(vx(j), interface_idx, nmax1)
                                r2 = fluxo_calor(vx(j), interface_idx, nmax2)
                                hy(j) = r2/r1
                            end do

                            s = 0.0

                            do j = 1, ord+1
                                t(j)=0.0
                                t(j+ord+1)=a
                            end do

                            allocate(wrk((tnmax*(ord+1)+nest*(7+3*ord))))

                            call curfit(iopt, tnmax, vx, hy, w, 0.0D0, a, ord, s, nest, nn, t, c, fpf, wrk, lwrk, iwrk, ier)

                            if (ier > 0) then
                                write(*, *)'Error in curfit. Ier = ', ier
                                stop
                            end if

                            call calculate_temperature_coefficients(interface_idx, condutance_idx, hest, .false.)

                            open(unit = 2, file = '/home/cx3d/mestrado/hest.dat')
                            do k = 1, tnmax
                                write(2, *)vx(k), hest(vx(k))
                            end do
                            close(2)

                            open(unit = 2, file = '/home/cx3d/mestrado/test.dat')
                            do k = 1, tnmax
                                write(2, *)vx(k), t1(vx(k), b), vy(k)
                            end do
                            close(2)

                            norm = 0.0
                            do k = 1, tnmax
                                tcalc = t1(vx(k), b)
                                norm = norm + (vy(k) - tcalc)**2
                            end do
                            norm = sqrt(norm/dble(tnmax))

                            deallocate(wrk)

                            deallocate(iwrk)
                            deallocate(t, c)
                            write(*, *)nmax1, nmax2, norm
                        end do
                    end do
                end do
            end do
        end do
    end block

    stop

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
        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)
        call calculate_reciprocity_coefficients(interface_idx)
        do condutance_idx = 1, 3
            write(str_cdx, '(I2.2)') condutance_idx
            call c_f_procpointer(hlist(condutance_idx), h)
            call calculate_temperature_coefficients(interface_idx, condutance_idx, h)
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

                ! write(*, *)'Integrating y = b'
                call calculate_integrals_Y(interface_idx, condutance_idx, stdev_idx)

                do nmax = 0, N
                    write(*, *)'N = ', nmax
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
            end do
        end do
    end do

contains
    function hest(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        double precision, dimension(1) :: argx, argy
        integer :: ier

        argx(1) = x

        call splev(t, nn, c, ord, argx, argy, 1, ier)
        if (ier /= 0) then
            write(*, *)'Error in curfit. Ier = ', ier
            stop
        end if
        r = argy(1)
    end function

end program main
