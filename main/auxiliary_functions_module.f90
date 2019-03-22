module auxiliary_functions_module
    use iso_c_binding
    use constants_module
    use integration_module
    implicit none

    type :: phi_t
        double precision, dimension(1001) :: x
        double precision, dimension(1001) :: y
        double precision, dimension(1001) :: dy
        double precision, dimension(1001) :: t
        double precision, dimension(1001) :: c
        double precision :: norm
        double precision :: sqrt_norm
        integer :: nest
        integer :: sz
    end type

    type(phi_t), dimension(0: mmax_phi) :: phi
contains
    function kernel(x, interface_idx) result(r)
        use interfaces_module
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx
        double precision :: r
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(dwlist(interface_idx), dw)

        r = sqrt(1 + dw(x)**2)
    end function

    function phi_eval(j, x) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x
        double precision :: r
        double precision, dimension(1) :: xx, yy
        integer :: ier, ord
        if (j == 0) then
            r = 1.0
        else
            ord = 5
            xx(1) = x
            call splev(phi(j)%t, phi(j)%nest, phi(j)%c, ord, xx, yy, 1, ier)
            if (ier /= 0) stop
            r = yy(1)
        end if
    end function

    function phi_eval_norm(j, x) result(r)
        integer, intent(in) :: j
        double precision, intent(in) :: x
        double precision :: r

        r = phi_eval(j, x)/phi(j)%sqrt_norm
    end function

    function transform_imp(f, k, j, interface_idx) result(r)
        use interfaces_module
        integer, intent(in) :: k, j, interface_idx
        interface
            function f(k, x, interface_idx) result(r)
                import
                integer, intent(in) :: k, interface_idx
                double precision, intent(in) :: x
                double precision :: r
            end function
        end interface
        double precision :: r

        r = integrate(f_aux, c_null_ptr, pts)
    contains
        function f_aux(x, args) result (r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            procedure(dw_proc_t), pointer :: dw
            double precision :: weight

            call c_f_procpointer(dwlist(interface_idx), dw)
            weight = kernel(x, interface_idx)
            r = f(k, x, interface_idx)*phi_eval(j, x)
        end function
    end function

    subroutine generate_eigenfunctions(interface_idx)
        use netlib_module
        use interfaces_module

        integer, intent(in) :: interface_idx
        double precision :: twant, tgot, dx, lambda1, lambda2, flambda1, flambda2, dlambda, lambda
        integer :: uflag
        procedure(dw_proc_t), pointer :: dw
        double precision, dimension(64) :: work
        double precision, dimension(2) :: ystart, thres, ygot, ypgot, ymax
        integer :: j, k

        call envirn(outch, mcheps, dwarf)
        call c_f_procpointer(dwlist(interface_idx), dw)

        thres = sqrt(dwarf)
        dx = a/1000.0
        dlambda = (pi/a)/20.0

        lambda1 = 0
        do j = 0, mmax_phi
            if (j /= 0) then
                lambda2 = lambda1 + dlambda
                flambda1 = g(lambda1)
                flambda2 = g(lambda2)
                do while (flambda1*flambda2 > 0)
                    lambda1 = lambda2
                    flambda1 = flambda2
                    lambda2 = lambda2 + dlambda
                    flambda2 = g(lambda2)
                end do

                lambda = zeroin(lambda1, lambda2, g, 1.0D-8)

                ! Compressao dos vetores obtidos da integracao. Por algum motivo a rotina UT repete valores calculados
                ! em alguns momentos. Aqui corrigimos os valores repetidos e ajustamos o tamanho do vetor
                !                k = 2
                !                do while (k .le. phi(j)%sz)
                !                    dif_x = dabs(phi(j)%x(k - 1) - phi(j)%x(k))
                !                    dif_y = dabs(phi(j)%y(k - 1) - phi(j)%y(k))
                !                    !if ((dif_x < 1.0D-8) .and. (dif_y < 1.0D-8)) then
                !                    if (dif_x < 1.0D-12) then
                !                        do i = k, phi(j)%sz
                !                            phi(j)%x(i - 1) = phi(j)%x(i)
                !                            phi(j)%y(i - 1) = phi(j)%y(i)
                !                        end do
                !                        phi(j)%sz = phi(j)%sz - 1
                !                    end if
                !                    k = k + 1
                !                end do

                block
                    integer :: iopt
                    integer :: ord
                    integer :: nn, ier, lwrk
                    double precision :: s, fp
                    double precision, dimension(2000) :: w
                    double precision, dimension(:), allocatable :: wrk
                    integer, dimension(:), allocatable :: iwrk

                    s = 0.0
                    iopt = 0
                    ord = 5
                    phi(j)%nest=phi(j)%sz+ord+1
                    lwrk = phi(j)%sz*(ord+1)+phi(j)%nest*(7+3*ord)
                    w = 1.0

                    allocate(iwrk(phi(j)%nest))

                    do k = 1, ord+1
                        phi(j)%t(k)=0.0
                        phi(j)%t(k+ord+1)=a
                    end do

                    allocate(wrk(lwrk))

                    call curfit(iopt, phi(j)%sz, phi(j)%x, phi(j)%y,&
                        w, 0.0D0, a, ord, s, phi(j)%nest, nn, phi(j)%t, phi(j)%c,fp, &
                        wrk, lwrk, iwrk, ier)
                    if (ier > 0) then
                        write(*, *)'Error in curfit'
                        stop
                    end if

                    deallocate (wrk, iwrk)
                end block
                lambda1 = lambda + dlambda
            end if
            ! Norma
            phi(j)%norm = transform_imp(f_aux, j, j, interface_idx)
            phi(j)%sqrt_norm = sqrt(phi(j)%norm)
        end do
    contains
        function g(lambda0) result (r)
            double precision, intent(in) :: lambda0
            double precision :: r
            logical(c_bool) :: keep, inside

            lambda = lambda0

            ystart(1) = 1.0
            ystart(2) = 0.0

            phi(j)%x = 0.0

            phi(j)%x(1) = 0.0
            phi(j)%y(1) = ystart(1)*cos(ystart(2))
            phi(j)%sz = 1

            call setup(2, 0.0D0, ystart, a, 1.0D-8, thres, 3, 'U', .true., &
                0.0D0, work, 64, .false.)
            twant = dx
            inside = .true.
            do while (inside)
                keep = .true.
                do while (keep)
                    call ut(f, twant, tgot, ygot, ypgot, ymax, work, uflag)
                    if (twant >= tgot) then
                        keep = uflag == 3
                    end if
                end do
                if (dabs(tgot - twant) < 1.0D-8) then
                    phi(j)%sz = phi(j)%sz + 1
                    phi(j)%x(phi(j)%sz) = tgot
                    phi(j)%y(phi(j)%sz) = ygot(1)*cos(ygot(2))
                    inside = dabs(tgot - a) > 1.0D-8
                    twant = tgot + dx
                end if
            end do

            r = ygot(2) - pi*dble(j)
        end function

        subroutine f(t, y, yp)
            double precision, intent(in) :: t
            double precision, dimension(*), intent(in) :: y
            double precision, dimension(*), intent(out) :: yp
            double precision :: weight

            weight = kernel(t, interface_idx)

            yp(1) = ((lambda**2)*weight - 1.0)*sin(y(2))*cos(y(2))*y(1)
            yp(2) = (sin(y(2))**2) + (lambda**2)*weight*(cos(y(2))**2)
        end subroutine

        function f_aux(k, x, interface_idx) result(r)
            integer, intent(in) :: k, interface_idx
            double precision, intent(in) :: x
            double precision :: r

            r = phi_eval(k, x)
        end function
    end subroutine
end module auxiliary_functions_module
