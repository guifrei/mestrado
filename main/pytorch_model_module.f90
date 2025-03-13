module pytorch_model_module
    use iso_c_binding
    implicit none

    interface
        function init_python() bind(C, name='init_python') result(r)
            use iso_c_binding
            integer(c_int) :: r
        end function init_python

        function load_weights() bind(C, name='init_python') result(r)
            use iso_c_binding
            integer(c_int) :: r
        end function load_weights

        function run_pytorch_model(input_value) bind(C, name='run_pytorch_model') result(r)
            use iso_c_binding
            real(c_double), value :: input_value
            real(c_double) :: r
        end function run_pytorch_model

        subroutine finish_python() bind(C, name='finish_python')
        end subroutine finish_python
    end interface

end module pytorch_model_module
