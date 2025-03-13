program main
    use constants_module
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use pytorch_model_module
    implicit none

    double precision :: x, y, dx
    character(len=100) :: line
    character(len=100) :: arg
    integer :: r

    dx = a/dble(tnmax - 1)

    call get_command_argument(1, arg)

    r = init_python()
    r = load_weights()

    if (arg == 'M') then ! temperaturas medidas na superficie
        call calculate_temperature_coefficients(w1, dw1, h3)
    else  if (arg == 'N') then ! temperaturas calculadas via rede neural
        call calculate_temperature_coefficients(w1, dw1, hnn)
    end if

    do
        read(*, '(A)', end=100) line
        read(line, *) x
        y = t1(x, b)
        write(*, *)y
    end do
    100 continue

    call finish_python()

end program main
