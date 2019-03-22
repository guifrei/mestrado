module sigmoid_function_module
    use iso_c_binding
    use constants_module
    implicit none

contains
    function sigmoid(x)
        double precision :: sigmoid
        double precision, intent(in) :: x

        sigmoid = 1.0/(1 + exp(-gamma*x))
    end function

    function sigmoid_d(x)
        double precision :: sigmoid_d, s
        double precision, intent(in) :: x

        s = sigmoid(x)
        sigmoid_d = gamma*s*(1.0 - s)
    end function
end module sigmoid_function_module
