module conductances_module
    use iso_c_binding
    use constants_module
    use sigmoid_function_module
    use pytorch_model_module
    implicit none
contains
    function hnn(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = run_pytorch_model(x)
    end function

    function h1(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if ((x < a/4.0).or.(x > 3.0*a/4.0)) then
            r = hmax
        else
            r = 0.0
        end if
    end function

    function h2(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = hmax*sin(pi*x/a)
    end function

    function h3(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        if ((x < a/4.0).or.((x>a/2.0).and.(x<3.0*a/4.0))) then
            r = hmax
        else if ((x>=a/4.0).and.(x<=a/2.0)) then
            r = hmax / 2.0
        else
            r = 0.0
        end if
    end function

    function h4(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if ((x < a / 4.0) .or. ((x > a / 2.0) .and. (x < 3.0 * a / 4.0))) then
            r = hmax
        else
            r = 0.0
        end if
    end function

    function h5(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = dabs(hmax * sin(2.0 * pi * x / a))
    end function

    function h6(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = hmax
    end function

    function h7(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = hmax - (hmax/2.0)*sigmoid(x-a/4.0) + (hmax/2.0)*sigmoid(x-a/2.0) - hmax*sigmoid(x-3.0*a/4.0)
    end function

    function h8(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if ((x < a / 4.0).or.(x > 3.0 * a / 4.0)) then
            r = hmax/2.0
        else
            r = 3.0*hmax/4.0
        end if
    end function

    function h9(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = hmax - hmax*sigmoid(x-a/4.0) + hmax*sigmoid(x-3.0*a/4.0)
    end function
end module conductances_module
