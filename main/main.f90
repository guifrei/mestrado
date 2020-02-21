!data00 = Import[
!   "/home/cx3d/mestrado/data/temperaturas_sinteticas_interface_03_\
!conductance_01_stdev_00.dat", {"Data", All, 2}];
!
!gaussData01 = GaussianFilter[data01, 20.0, Method -> "Gaussian"];
!ListLinePlot[{data00, gaussData01}]

!plot 'conductance_01.dat', 'estimativa_ctc_interface_01_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_01_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_01_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_01_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_01_conductance_03_stdev_05.dat'
!
!
!===============
!
!
!plot 'conductance_01.dat', 'estimativa_ctc_interface_02_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_02_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_02_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_02_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_02_conductance_03_stdev_05.dat'
!
!===============
!
!
!plot 'conductance_01.dat', 'estimativa_ctc_interface_03_conductance_01_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_01_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_01_stdev_05.dat'
!
!plot 'conductance_02.dat', 'estimativa_ctc_interface_03_conductance_02_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_02_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_02_stdev_05.dat'
!
!plot 'conductance_03.dat', 'estimativa_ctc_interface_03_conductance_03_stdev_00.dat' with lines, 'estimativa_ctc_interface_03_conductance_03_stdev_01.dat', 'estimativa_ctc_interface_03_conductance_03_stdev_05.dat'


program main
    use interfaces_module
    use conductances_module
    use temperature_functions_module
    use reciprocity_functions_module
    use netlib_module
    use estimated_h_module
    use tikhonov_module
    use morozov_module
    use iso_c_binding
    implicit none

    procedure(h_proc_t), pointer :: h
    procedure(w_proc_t), pointer :: w
    procedure(dw_proc_t), pointer :: dw
    double precision :: c_fluxo_calor, c_delta_temperatura
    double precision, dimension(tnmax) :: vx, vy
    integer :: interface_idx, condutance_idx, nmax, k, j
    double precision :: x, dx, stdev, sqrt_rms
    double precision :: h_est
    integer :: kmax
    double precision :: lambda
    logical :: success
    double precision, dimension(tnmax) :: tmpvy

    interface generic_spec
        function gsl_rng_alloc(gsl_rng_type) result(gsl_rng) bind(C)
            import
            type(c_ptr), value :: gsl_rng_type
            type(c_ptr) :: gsl_rng
        end function
    end interface generic_spec

!    type(c_ptr), bind(C, name = '') :: gsl_rng_type
    type(c_ptr) ::  gsl_rng

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

    dx = a/dble(tmax - 1)

    interface_idx = 1
    condutance_idx = 1

    call c_f_procpointer(hlist(condutance_idx), h)
    call c_f_procpointer(wlist(interface_idx), w)
    call c_f_procpointer(dwlist(interface_idx), dw)

    call calculate_temperature_coefficients(interface_idx, condutance_idx, h, vx, vy)
    stdev = 0.5
    call add_error(vy, stdev)
    !call least_squares_for_Y(vx, vy, vvY)
    lambda = 0.66
!    call tikhonov2(lambda, vx, vy, vvY)
    call find_root(obj, 0.66D0, 1.0D-4, 1000, lambda, success)
    call tikhonov2(lambda, vx, vy, vvY)
    write(*, *)obj(lambda)

    !Principio da discrepancia de Morozov
    kmax = N
    call morozov(stdev, vx, vy, vvY, kmax)

    call calculate_reciprocity_coefficients(interface_idx)

    do j = 0, kmax
        reciprocity_f(j) = calc_reciprocity_f(j)
        reciprocity_g(j) = calc_reciprocity_g(j)
    end do

    open(7, file = '/tmp/log')
    nmax = kmax
    do j = 1, tnmax
        x = vx(j)
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

contains
    double precision function obj(lambda)
        double precision, intent(in) :: lambda
        double precision, dimension(tnmax) :: tmpvy
        double precision :: sqrt_rms

        call tikhonov2(lambda, vx, vy, vvY)
        tmpvy = 0.0
        do j = 0, N
            tmpvy = tmpvy + vvY(j)*cos(mu(j)*vx)
        end do
        sqrt_rms = norm2(tmpvy - vy)/sqrt(dble(tnmax))
        obj = sqrt_rms - stdev
    end function obj

    subroutine find_root( f, xinit, tol, maxiter, result, success )
        interface
            double precision function f(x)
                double precision, intent(in) :: x
            end function f
        end interface

        double precision, intent(in)     :: xinit
        double precision, intent(in)     :: tol
        integer, intent(in)  :: maxiter
        double precision, intent(out)    :: result
        logical, intent(out) :: success

        double precision                 :: eps = 0.00001
        double precision                 :: fx1
        double precision                 :: fx2
        double precision                 :: fprime
        double precision                 :: x
        double precision                 :: xnew
        integer              :: i

        result  = 0.0
        success = .false.

        x = xinit
        do i = 1,max(1,maxiter)
            fx1    = f(x)
            fx2    = f(x+eps)
            write(*,*) i, fx1, fx2, eps
            fprime = (fx2 - fx1) / eps

            xnew   = x - fx1 / fprime

            if ( abs(xnew-x) <= tol ) then
                success = .true.
                result  = xnew
                exit
            endif

            x = xnew
            write(*,*) i, x
        enddo

    end subroutine find_root
end program main
