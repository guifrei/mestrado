module interfaces_module
    use iso_c_binding
    use constants_module
    use sigmoid_function_module
    implicit none

    interface
        function w_proc_t(x) result(r)
            import
            double precision, intent(in) :: x
            double precision :: r
        end function

        function dw_proc_t(x) result(r)
            import
            double precision, intent(in) :: x
            double precision :: r
        end function
    end interface

    integer, parameter :: wmax = 15

    type(c_funptr), dimension(wmax) :: wlist
    type(c_funptr), dimension(wmax) :: dwlist

contains
    function w1(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = b/2.0
    end function

    function dw1(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = 0.0
    end function

    function w2(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x <= a/3.0) then
            r = 19.0*b*x**2/(8.0*a**2)-7.0*b*x/(24.0*a)+b/2.0
        else if (x <= 2.0*a/3.0) then
            r = -25.0*b*x**2/(8.0*a**2)+27.0*b*x/(8.0*a)-b/9.0
        else
            r = b*x**2/(8.0*a**2)-23.0*b*x/(24.0*a)+4.0*b/3.0
        end if
    end function

    function dw2(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x <= a/3.0) then
            r = 38.0*b*x/(8.0*a**2)-7.0*b/(24.0*a)
        else if (x <= 2.0*a/3.0) then
            r = -50.0*b*x/(8.0*a**2)+27.0*b/(8.0*a)
        else
            r = 2.0*b*x/(8.0*a**2)-23.0*b/(24.0*a)
        end if
    end function

    function w3(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = b/2.0+(b/20.0)*cos(mp*pi*x/a)
    end function

    function dw3(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = -(b/20.0)*(mp*pi/a)*sin(mp*pi*x/a)
    end function

    function w4(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x <= a / 8.0) then
            r = b/2
        else if (x <= a / 4.0) then
            r = a*a/64 + b/2 - a*x/4+x*x
        else if (x <= 3.0 * a / 4.0) then
            r = -5*a*a/64+b/2+a*x/2-x*x/2
        else if (x <= 7.0 * a / 8.0) then
            r = 49*a*a/64+b/2-7*a*x/4+x*x
        else
            r = b/2
        end if
    end function

    function dw4(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x <= a / 8.0) then
            r = 0
        else if (x <= a / 4.0) then
            r = - a/4+2*x
        else if (x <= 3.0 * a / 4.0) then
            r = a/2-x
        else if (x <= 7.0 * a / 8.0) then
            r = -7*a/4+2*x
        else
            r = 0
        end if
    end function

    function w5(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = (b/4.0)*(sigmoid(x-a/4.0) - sigmoid(x-3.0*a/4.0)) + b/2.0
    end function

    function dw5(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = (b/4.0)*(sigmoid_d(x-a/4.0) - sigmoid_d(x-3.0*a/4.0))
    end function

    function w6(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = b/4.0
    end function

    function dw6(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = 0.0
    end function

    function w7(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = b/2.0 + (b/10.0)*cos(2.0*pi*x/a) + (b/40.0)*sin(4.0*pi*x/a)
    end function

    function dw7(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r =  -(b/10.0)*(2.0*pi/a)*sin(2.0*pi*x/a) + (b/40.0)*(4.0*pi/a)*cos(4.0*pi*x/a)
    end function

    function w8(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = b/2.0 + (b/20.0)*sin(2.0*pi*x/a) - (b/20.0)*cos(8.0*pi*x/a)
    end function

    function dw8(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r =  (b/20.0)*(2.0*pi/a)*cos(2.0*pi*x/a) + (b/20.0)*(8.0*pi/a)*sin(8.0*pi*x/a)
    end function

    function w9(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = 3.0*b/4.0 -(b/2.0)*x/a
    end function

    function dw9(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = -(b/2.0)/a
    end function

    function w10(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/2) then
            r = 3.0*b/4.0 - b*x/(2.0*a)
        else
            r = b/2.0 + b*(x - a/2.0)/(2.0*a)
        end if
    end function

    function dw10(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/2) then
            r = - b/(2.0*a)
        else
            r = b/(2.0*a)
        end if
    end function

    function w11(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/10.0) then
            r = b1-50.0*(b1-b2)*((x/a)**2)
        else if (x < 3.0*a/10.0) then
            r = b2+50.0*(b1-b2)*((1.0/5.0-x/a)**2)
        else if (x < a/2.0) then
            r = b1-50.0*(b1-b2)*((2.0/5.0-x/a)**2)
        else if (x < 7.0*a/10.0) then
            r = b2+50.0*(b1-b2)*((3.0/5.0-x/a)**2)
        else if (x < 9.0*a/10.0) then
            r = b1-50.0*(b1-b2)*((4.0/5.0-x/a)**2)
        else
            r = b2+50.0*(b1-b2)*((1.0-x/a)**2)
        end if

    end function

    function dw11(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/10.0) then
            r = -(100.0/a)*(b1-b2)*(x/a)
        else if (x < 3.0*a/10.0) then
            r = -(100.0/a)*(b1-b2)*(1.0/5.0-x/a)
        else if (x < a/2.0) then
            r = (100.0/a)*(b1-b2)*(2.0/5.0-x/a)
        else if (x < 7.0*a/10.0) then
            r = -(100.0/a)*(b1-b2)*(3.0/5.0-x/a)
        else if (x < 9.0*a/10.0) then
            r = (100.0/a)*(b1-b2)*(4.0/5.0-x/a)
        else
            r = -(100.0/a)*(b1-b2)*(1.0-x/a)
        end if
    end function

    function w12(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/4.0) then
            r = b1-4.0*(b1-b2)*x/a
        else if (x < a/2.0) then
            r = b2+4.0*(b1-b2)*(x-a/4.0)/a
        else if (x < 3.0*a/4.0) then
            r = b1-4.0*(b1-b2)*(x-a/2.0)/a
        else
            r = b2+4.0*(b1-b2)*(x-3.0*a/4.0)/a
        end if
    end function

    function dw12(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x < a/4.0) then
            r = -4.0*(b1-b2)/a
        else if (x < a/2.0) then
            r = 4.0*(b1-b2)/a
        else if (x < 3.0*a/4.0) then
            r = -4.0*(b1-b2)/a
        else
            r = 4.0*(b1-b2)/a
        end if
    end function

    function w13(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = (b/2) + (3*b/4-b/2)*sigmoid(x-a/4)+(b/2-3*b/4)*sigmoid(x-a/2)+ (b/4-b/2)*sigmoid(x-3*a/4)
    end function

    function dw13(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = (3*b/4-b/2)*sigmoid_d(x-a/4)+(b/2-3*b/4)*sigmoid_d(x-a/2)+ (b/4-b/2)*sigmoid_d(x-3*a/4)
    end function

    function w14(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        if (x <= a/4.0) then
            r = b/2.0
        else if (x <= a/2.0) then
            r = 3.0*b/4.0
        else if (x <= 3.0*a/4.0) then
            r = b/2.0
        else
            r = b/4.0
        end if

    end function

    function dw14(x) result(r)
        double precision, intent(in) :: x
        double precision :: r

        r = 0.0
    end function

    function w15(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        double precision :: deltax = a/100.0
        double precision :: deltaw = b/5.0

        if (x <= a/2.0 - deltax) then
            r = b/2.0 - deltaw
        else if (x <= a/2.0 + deltax) then
            r = b/2.0 - deltaw + (deltaw/deltax)*(x - a/2.0 + deltax)
        else
            r = b/2.0 + deltaw
        end if

    end function

    function dw15(x) result(r)
        double precision, intent(in) :: x
        double precision :: r
        double precision :: deltax = a/100.0
        double precision :: deltaw = b/5.0

        if (x <= a/2.0 - deltax) then
            r = 0.0
        else if (x <= a/2.0 + deltax) then
            r = deltaw/deltax
        else
            r = 0.0
        end if
    end function

end module interfaces_module
