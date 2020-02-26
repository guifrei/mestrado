program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use estimated_h_module
    !use tikhonov_module
    !use morozov_module
    use iso_c_binding
    implicit none

    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax) :: vx, vy, tmpy
    integer :: interface_idx, condutance_idx, stdev_idx, nmax, k, j
    double precision :: x, dx, sqrt_rms
    double precision :: h_est
    integer, dimension(3, 3, 3) :: kmax
    double precision :: lambda
    logical :: success
    double precision, dimension(3) :: stdev = [0.0D0, 0.1D0, 0.5D0]
    character(len = 2) :: str_idx, str_cdx, str_stdev

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

    dx = a/1000.0
    open(unit=7,file='../paper/interface.dat')
    open(unit=8,file='../paper/conductance.dat')
    do j = 0, 1000
        x = dx*dble(j)
        write(7, *)x, w1(x), w2(x), w3(x)
        write(8, *)x, h1(x), h2(x), h3(x)
    end do
    close(7)
    close(8)

    do interface_idx = 1, 3
        do condutance_idx = 1, 3
            do stdev_idx = 1, 3

                write(*, *)'Interface = ', interface_idx, ', conductance = ', condutance_idx, ', stdev = ', stdev_idx
                write(str_idx, '(I2.2)') interface_idx
                write(str_cdx, '(I2.2)') condutance_idx
                if (stdev_idx.eq.1) then
                    str_stdev = '00'
                else if (stdev_idx.eq.2) then
                    str_stdev = '01'
                else
                    str_stdev = '05'
                end if

                call c_f_procpointer(hlist(condutance_idx), h)
                call c_f_procpointer(wlist(interface_idx), w)
                call c_f_procpointer(dwlist(interface_idx), dw)

                call calculate_temperature_coefficients(interface_idx, condutance_idx, h, vx, vy)
                !                stdev = 0.5!maxval(abs(vy))*0.1/100.0
                call add_error(vy, stdev(stdev_idx))
                call least_squares_for_Y(vx, vy, vvY)
                kmax = N

                open(unit=10,file='../paper/difference_interface_' // &
                    str_idx // '_conductance_' // str_cdx // '_stdev_' // str_stdev // '.dat')
                if (stdev(stdev_idx).ne.0) then
                    tmpy = 0.0
                    do j = 0, mmax_phi
                        tmpy = tmpy + vvY(j)*cos(mu(j)*vx)
                        sqrt_rms = norm2(tmpy - vy)
                        write(10, *)j, sqrt_rms !norm2(vvY(i)*cos(mu(i)*vx))
                    end do
                end if
                close(10)

                !    kmax = 6
                !    vy = 0.0
                !    do j = 0, kmax
                !        vy = vy + vvY(j)*cos(mu(j)*vx)
                !    end do

                call calculate_reciprocity_coefficients(interface_idx)

                nmax = kmax(interface_idx, condutance_idx, stdev_idx)

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
                        c_fluxo_calor = c_fluxo_calor + parcela_fluxo_calor(x, k, interface_idx)
                        c_delta_temperatura = c_delta_temperatura + parcela_delta_temperatura(x, k, interface_idx)
                    end do
                    h_est = c_fluxo_calor/c_delta_temperatura
                    write(7, *)x, h_est, h(x)
                end do
                close(7)
            end do
        end do
    end do
end program main
