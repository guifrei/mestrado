program main
    use constants_module
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    !use reciprocity_functions_module
    !use netlib_module
    !use estimated_h_module
    !use tikhonov_module
    implicit none

    double precision, dimension(tnmax, 0: N), target :: m_fluxo_calor, m_delta_temperatura
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax) :: vx, vy
    double precision :: ymax
    character(len = 2) :: str_N
    integer :: nmax, k, j
    double precision :: desv, x, y, y1, y2, dx, stdev
    double precision :: fluxo_calor_teorico, delta_temperatura_teorico, norm_f, norm_t
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

    call calculate_temperature_coefficients(w1, dw1, h1)

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
end program main
