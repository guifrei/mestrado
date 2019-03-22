module reciprocity_functions_module
    use iso_c_binding
    use constants_module
    use eigenfunctions_module
    use interfaces_module
    use netlib_module
    use iso_fortran_env
    use integration_module
    implicit none

    double precision, dimension(0: 2*mmax_F + 1, 0: N), target :: coeffsF
    double precision, dimension(0: mmax_G, 0: N), target :: coeffsG
    double precision, dimension(0: N, 0: mmax_phi), target :: vpsi, vphi
    double precision, dimension(0: mmax_phi) :: integrals_Y
    logical :: factored

    integer, parameter :: ID_F1 = 0
    integer, parameter :: ID_DF1DX = 1
    integer, parameter :: ID_DF1DY = 2
    integer, parameter :: ID_F2 = 3
    integer, parameter :: ID_DF2DX = 4
    integer, parameter :: ID_DF2DY = 5
    integer, parameter :: ID_G1 = 6
    integer, parameter :: ID_DG1DX = 7
    integer, parameter :: ID_DG1DY = 8
contains
    function f_aux_psi(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r

        r = cos(mu(j)*x)

        if (j == 0) then
            r = r/sqrt(a)
        else
            r = r/sqrt(a/2.0)
        end if
    end function

    function f_aux_phi(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r

        r = cos(mu(j)*x)

        if (j == 0) then
            r = r/sqrt(a)
        else
            r = r/sqrt(a/2.0)
        end if
    end function

    function fpsi(j, x) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j
        double precision :: r
        integer :: m

        r = vpsi(j, 0)/a
        do m = 1, mmax_phi
            r = r + (2.0/a)*vpsi(j, m)*cos(mu(m)*x)
        end do
    end function

    function fphi(j, x) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j
        double precision :: r
        integer :: m

        r = vphi(j, 0)/a
        do m = 1, mmax_phi
            r = r + (2.0/a)*vphi(j, m)*cos(mu(m)*x)
        end do
    end function

    subroutine calculate_integrals_Y(interface_idx, condutance_idx, stdev_idx)
        integer, intent(in) :: interface_idx
        integer, intent(in) :: condutance_idx
        integer, intent(in) :: stdev_idx
        double precision :: stdev
        integer :: k, j
        double precision, dimension(tnmax) :: vx
        double precision, dimension(tnmax) :: sample_y
        double precision, dimension(tnmax) :: vy
        character(len = 2) :: str_idx, str_cdx, str_stdev
        integer :: iopt
        integer :: ord
        integer :: nest
        integer :: nn, ier, lwrk
        double precision, dimension(tnmax) :: w
        double precision, dimension(:), allocatable :: wrk
        double precision :: s, fp
        double precision, dimension(:), allocatable :: t, c
        integer, dimension(:), allocatable :: iwrk
        double precision, dimension(tnmax) :: rand_u, rand_v, rand_epsilon

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

        open(unit = 1, file = '/home/cx3d/mestrado/' // &
            'data/comsol/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx//'.dat')
        do k = 1, tnmax
            read(1, *)vx(k), sample_y(k)
        end do
        close(1)

        if (stdev > 0) then
            call random_seed
            call random_number(rand_u)
            call random_number(rand_v)
            rand_epsilon = cos(2.0*pi*rand_v)*sqrt(-2.0*log(rand_u))
            sample_y = sample_y + stdev*rand_epsilon
        end if

        open(unit = 1, file = '/home/cx3d/mestrado/' // &
            'data/temperaturas_sinteticas_interface_'//str_idx//'_conductance_'//str_cdx // &
            '_stdev_'// str_stdev // '.dat')
        do k = 1, tnmax
            write(1, *)vx(k), sample_y(k)
        end do
        close(1)

        iopt = 0
        ord = 5
        nest=tnmax+ord+1
        lwrk = tnmax*(ord+1)+nest*(7+3*ord)
        w = 1.0

        allocate(t(nest), c(nest), iwrk(nest))

        do k = 0, mmax_phi
            do j = 1, tnmax
                vy(j) = sample_y(j)*cos(mu(k)*vx(j))
            end do

            s = 0.0

            do j = 1, ord+1
                t(j)=0.0
                t(j+ord+1)=a
            end do

            allocate(wrk((tnmax*(ord+1)+nest*(7+3*ord))))

            call curfit(iopt,tnmax,vx,vy,w,0.0D0,a,ord,s,nest,nn,t,c,fp, &
                wrk,lwrk,iwrk,ier)
            if (ier > 0) then
                write(*, *)'Error in curfit. Ier = ', ier
                stop
            end if

            deallocate (wrk)
            allocate(wrk(nn))
            integrals_Y(k) = splint(t,nn,c,ord,0.0D0,a,wrk)
            deallocate(wrk)
        end do
        deallocate(t, c, iwrk)
    end subroutine

    function reciprocity_f(j) result(r)
        integer, intent(in) :: j
        double precision :: r
        integer :: m
        double precision, dimension(:), pointer :: vA

        vA(0:) => coeffsF(0::2, j)
        r = -q * vpsi(j, 0) / k1 + (vA(0) - vpsi(j, 0)) * integrals_Y(0) / (a * b)
        do m = 1, mmax_F
            r = r + (2.0/a) * mu(m) * (vA(m) / sinh(mu(m) * b) - vpsi(j, m) / tanh(mu(m) * b)) * integrals_Y(m)
        end do
    end function

    function parcela_delta_temperatura(x, j, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j
        integer, intent(in) :: interface_idx
        double precision :: r

        r = k1*reciprocity_f(j)*fbeta(j, x, interface_idx)
    end function

    function delta_temperatura(x, interface_idx, kmax) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx, kmax
        double precision :: r
        integer :: j

        r = 0
        do j = 0, kmax
            r = r + parcela_delta_temperatura(x, j, interface_idx)
        end do
    end function

    function reciprocity_g(j) result(r)
        integer, intent(in) :: j
        double precision :: r
        integer :: m
        double precision, dimension(:), pointer :: vE

        vE(0:) => coeffsG(0:, j)
        r = -q * vphi(j, 0) / k1 + (vE(0) - vphi(j, 0)) * integrals_Y(0) / (a * b)
        do m = 1, mmax_G
            r = r + (2.0/a) * mu(m) * (vE(m) / sinh(mu(m) * b) - vphi(j, m) / tanh(mu(m) * b)) * integrals_Y(m)
        end do
    end function

    function parcela_fluxo_calor(x, j, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j
        integer, intent(in) :: interface_idx
        double precision :: r

        r = k1*reciprocity_g(j)*fgamma(j, x, interface_idx)
    end function

    function fluxo_calor(x, interface_idx, kmax) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: interface_idx, kmax
        double precision :: r
        integer :: j

        r = 0
        do j = 0, kmax
            r = r + parcela_fluxo_calor(x, j, interface_idx)
        end do
    end function

    subroutine gram_schmidt(interface_idx)
        integer, intent(in) :: interface_idx
        integer :: j, m
        double precision :: prod
        do j = 0, N
            associate (vAj => coeffsF(0::2, j), vDj => coeffsF(1::2, j), vEj => coeffsG(0:, j))
                do m = 0, j - 1
                    associate (vAk => coeffsF(0::2, m), vDk => coeffsF(1::2, m), vEk => coeffsG(0:, m))
                        prod = dot_product_beta(j, m, interface_idx)
                        vAj = vAj - prod*vAk
                        vDj = vDj - prod*vDk
                        vpsi(j, :) = vpsi(j, :) - prod*vpsi(m, :)
                        prod = dot_product_gamma(j, m, interface_idx)
                        vEj = vEj - prod*vEk
                        vphi(j, :) = vphi(j, :) - prod*vphi(m, :)
                    end associate
                end do
                prod = sqrt(dot_product_beta(j, j, interface_idx))
                vAj = vAj/prod
                vDj = vDj/prod
                vpsi(j, :) = vpsi(j, :)/prod
                prod = sqrt(dot_product_gamma(j, j, interface_idx))
                vEj = vEj/prod
                vphi(j, :) = vphi(j, :)/prod
            end associate
        end do
    end subroutine gram_schmidt

    subroutine calculate_reciprocity_coefficients(interface_idx, conductance_idx)
        integer, intent(in) :: interface_idx, conductance_idx
        integer :: j, m
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision, dimension(0: 2*mmax_F + 1, 0: N), target :: tmpcoeffsF
        double precision, dimension(0: mmax_G, 0: N), target :: tmpcoeffsG
        double precision, dimension(0: 2*mmax_F + 1) :: rF, cF
        double precision, dimension(0: mmax_G) :: rG, cG
        integer :: info
        double precision, dimension(2*mmax_F+2, 2*mmax_F+2) :: afF
        double precision, dimension(mmax_G+1, mmax_G+1) :: afG
        integer, dimension(2*mmax_F+2) :: ipivF
        integer, dimension(mmax_G+1) :: ipivG
        double precision :: rcond
        double precision, dimension(N+1) :: ferr, berr
        double precision, dimension(8*mmax_F+8) :: workF
        double precision, dimension(4*mmax_G+4) :: workG
        integer, dimension(2*mmax_F+2) :: iworkF
        integer, dimension(mmax_G+1) :: iworkG
        character(1) :: equed, fact
        double precision, dimension(0: 2*mmax_F + 1, 0: 2*mmax_F + 1) :: mxF
        double precision, dimension(0: mmax_G, 0: mmax_G) :: mxG
        character(len = 2) :: str_idx, str_cdx

        write(str_idx, '(I2.2)') interface_idx
        write(str_cdx, '(I2.2)') conductance_idx

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (factored) then

            !        do j = 0, N
            !            do m = 0, mmax_phi
            !                vpsi(j, m) = ftransform(f_aux_psi, j, m, interface_idx)
            !                vphi(j, m) = ftransform(f_aux_phi, j, m, interface_idx)
            !            end do
            !        end do

            vpsi = 0.0
            vphi = 0.0

            do j = 0, N
                if (j == 0) then
                    vpsi(j, j) = sqrt(a)
                    vphi(j, j) = sqrt(a)
                else
                    vpsi(j, j) = sqrt(a/2.0)
                    vphi(j, j) = sqrt(a/2.0)
                end if
            end do

            write(*, *)'Generating matrices for F'

            mxF = 0.0
            do j = 0, mmax_F
                do m = 0, mmax_F
                    mxF(j*2,m*2) = ftransform(fa, m, j, interface_idx)
                    mxF(j*2,m*2+1) = ftransform(fb, m, j, interface_idx)
                    mxF(j*2+1,m*2) = ftransform(fp, m, j, interface_idx)
                    mxF(j*2+1,m*2+1) = ftransform(fq, m, j, interface_idx)
                end do
                do m = 0, N
                    coeffsF(j*2, m) = ftransform(fc, m, j, interface_idx)
                    coeffsF(j*2+1, m) = ftransform(fr, m, j, interface_idx)
                end do
            end do

            write(*, *)'Generating matrices for G'

            mxG = 0.0
            do j = 0, mmax_G
                do m = 0, mmax_G
                    mxG(j,m) = ftransform(fu, m, j, interface_idx)
                end do
                do m = 0, N
                    coeffsG(j, m) = ftransform(fv, m, j, interface_idx)
                end do
            end do

!            factored = .false.

        end if

        ! Solucao do sistema
        write(*, *)'Solving systems'
        if (factored) then
            fact = 'F'
            equed = 'B'
        else
            fact = 'E'
            equed = 'B'
        end if

        call dgesvx('E', 'N', 2*mmax_F+2, N+1, mxF, 2*mmax_F+2, afF, 2*mmax_F+2, ipivF, equed, &
            rF, cF, coeffsF, 2*mmax_F+2, tmpcoeffsF, 2*mmax_F+2, rcond, ferr, berr, workF, iworkF, info)
        coeffsF = tmpcoeffsF

        call dgesvx('E', 'N', mmax_G+1, N+1, mxG, mmax_G+1, afG, mmax_G+1, ipivG, equed, &
            rG, cG, coeffsG, mmax_G+1, tmpcoeffsG, mmax_G+1, rcond, ferr, berr, workG, iworkG, info)
        coeffsG = tmpcoeffsG

        ! Algoritmo de ortogonalizacao de Gram-Schmidt
        write(*, *)'Gram Schmidt'
        call gram_schmidt(interface_idx)

    !        block
    !            double precision :: x
    !            integer :: r_idx
    !
    !            r_idx = 5
    !            open(unit = 1, file='/home/cx3d/mestrado/data/a.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, F1(r_idx, x, b), fpsi(r_idx, x)
    !            end do
    !            close(1)
    !
    !            open(unit = 1, file='/home/cx3d/mestrado/data/b.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, k1*dF1dn(r_idx, x, interface_idx), k2*dF2dn(r_idx, x, interface_idx)
    !            end do
    !            close(1)
    !
    !            open(unit = 1, file='/home/cx3d/mestrado/data/c.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, F1(r_idx, x, w(x)), F2(r_idx, x, w(x))
    !            end do
    !            close(1)
    !
    !            open(unit = 1, file='/home/cx3d/mestrado/data/d.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, F2(r_idx, x, 0.0_8)
    !            end do
    !            close(1)
    !
    !            open(unit = 1, file='/home/cx3d/mestrado/data/e.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, G1(r_idx, x, b), fphi(r_idx, x)
    !            end do
    !            close(1)
    !
    !            open(unit = 1, file='/home/cx3d/mestrado/data/f.txt')
    !            do j = 0, 100
    !                x = dble(j)*a/100.0
    !                write(1, *)x, dG1dn(r_idx, x, interface_idx)
    !            end do
    !            close(1)
    !        end block
    end subroutine

    function soma_controle_erro(parcela, j, x, y, v, mmax, id) result(r)
        interface
            function parcela(j, m, x, y, va) result(r)
                import
                integer, intent(in) :: j, m
                double precision, intent(in) :: x, y, va
                double precision :: r
            end function
        end interface
        double precision, intent(in) :: x, y
        double precision, dimension(0:), intent(in) :: v
        integer, intent(in) :: id
        integer, intent(in) :: j
        integer, intent(in) :: mmax
        double precision :: r, partial_r, eps, r_acc
        integer :: m, p
        logical :: keep_1, keep_2, converged

        eps = 0.0
        r_acc = parcela(j, 0, x, y, v(0))
        m = 1
        keep_1 = .true.
        do while (keep_1)
            r_acc = r_acc + parcela(j, m, x, y, v(m))
            m = m + 1
            keep_1 = m <= mmax
        end do

        do while (keep_1)
            p = m
            partial_r = 0.0
            keep_2 = .true.
            do while (keep_2)
                partial_r = partial_r + parcela(j, m, x, y, v(m))
                m = m + 1
                keep_1 = m <= mmax
                keep_2 = keep_1 .and. (m < p + delta_m)
            end do
            r_acc = r_acc + partial_r
            eps = dabs(partial_r/r_acc)
            converged = eps .lt. reltol
            keep_1 = keep_1 .and. (.not. converged)
        end do
        r = r_acc
    !        if (.not. converged) then
    !            if (id == ID_F1) then
    !                write(*, *)'F1   ', ': eps = ', eps, ', x = ', x, ', y = ', y, ', j = ', j
    !            else if (id == ID_F2) then
    !                write(*, *)'F2   ', ': eps = ', eps, ', x = ', x, ', y = ', y, ', j = ', j
    !            else if (id == ID_G1) then
    !                if (eps .eq. 0) then
    !                    write(*, *)'G1   ', ': val = ', r, ', x = ', x, ', y = ', y, ', j = ', j
    !                else
    !                    write(*, *)'G1   ', ': eps = ', eps, ', x = ', x, ', y = ', y, ', j = ', j
    !                end if
    !            end if
    !        end if
    end function

    function parcela_F1(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r

        if (m == 0) then
            r = (va*(b - y) + vpsi(j, 0)*y)/(a*b)
        else
            r = (2.0/a)*(va*sinh(mu(m)*(b - y))/sinh(mu(m)*b) +&
                vpsi(j, m)*sinh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function parcela_dF1dx(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r

        if (m == 0) then
            r = 0.0
        else
            r = - mu(m)*(2.0/a)*(va*sinh(mu(m)*(b - y))/sinh(mu(m)*b) +&
                vpsi(j, m)*sinh(mu(m)*y)/sinh(mu(m)*b))*sin(mu(m)*x)
        end if
    end function

    function parcela_dF1dy(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r

        if (m == 0) then
            r = (vpsi(j, 0) - va)/(a*b)
        else
            r = mu(m)*(2.0/a)*(-va*cosh(mu(m)*(b - y))/sinh(mu(m)*b) +&
                vpsi(j, m)*cosh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function parcela_F2(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r

        if (m == 0) then
            r = vd*y/(a*b)
        else
            r = (2.0/a)*(vd*sinh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function parcela_dF2dx(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r

        if (m == 0) then
            r = 0.0
        else
            r = - mu(m)*(2.0/a)*(vd*sinh(mu(m)*y)/sinh(mu(m)*b))*sin(mu(m)*x)
        end if
    end function

    function parcela_dF2dy(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r

        if (m == 0) then
            r = vd/(a*b)
        else
            r = mu(m)*(2.0/a)*(vd*cosh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function parcela_G1(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r

        if (m == 0) then
            r = (ve*(b - y) + vphi(j, 0)*y)/(a*b)
        else
            r = (2.0/a)*(ve*sinh(mu(m)*(b - y))/sinh(mu(m)*b) +&
                vphi(j, m)*sinh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function parcela_dG1dx(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r

        if (m == 0) then
            r = 0.0
        else
            r = - mu(m)*(2.0/a)*(ve*sinh(mu(m)*(b - y))/sinh(mu(m)*b) +&
                vphi(j, m)*sinh(mu(m)*y)/sinh(mu(m)*b))*sin(mu(m)*x)
        end if
    end function

    function parcela_dG1dy(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r

        if (m == 0) then
            r = (vpsi(j, 0) - ve)/(a*b)
        else
            r = mu(m)*(2.0/a)*(-ve*cosh(mu(m)*(b - y))/sinh(mu(m)*b) + &
                vphi(j, m)*cosh(mu(m)*y)/sinh(mu(m)*b))*cos(mu(m)*x)
        end if
    end function

    function F1(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_F1, j, x, y, coeffsF(0::2, j), mmax_F, ID_F1)
    end function

    function dF1dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF1dx, j, x, y, coeffsF(0::2, j), mmax_F, ID_DF1DX)
    end function

    function dF1dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF1dy, j, x, y, coeffsF(0::2, j), mmax_F, ID_DF1DY)
    end function

    function F2(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_F2, j, x, y, coeffsF(1::2, j), mmax_F, ID_F2)
    end function

    function dF2dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF2dx, j, x, y, coeffsF(1::2, j), mmax_F, ID_DF2DX)
    end function

    function dF2dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF2dy, j, x, y, coeffsF(1::2, j), mmax_F, ID_DF2DY)
    end function

    function G1(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_G1, j, x, y, coeffsG(0:, j), mmax_G, ID_G1)
    end function

    function dG1dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dG1dx, j, x, y, coeffsG(0:, j), mmax_G, ID_DG1DX)
    end function

    function dG1dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dG1dy, j, x, y, coeffsG(0:, j), mmax_G, ID_DG1DY)
    end function

    function dF1dn(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        procedure(dw_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        r = (dw(x)*dF1dx(j, x, w(x)) - dF1dy(j, x, w(x)))/sqrt(1 + dw(x)**2)
    end function

    function dF2dn(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        r = (dw(x)*dF2dx(j, x, w(x)) - dF2dy(j, x, w(x)))/sqrt(1 + dw(x)**2)
    end function

    function dG1dn(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        r = (dw(x)*dG1dx(j, x, w(x)) - dG1dy(j, x, w(x)))/sqrt(1 + dw(x)**2)
    end function

    function fbeta(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r

        r = k1*dF1dn(j, x, interface_idx)
    end function

    function fgamma(j, x,interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        procedure(w_proc_t), pointer :: w

        call c_f_procpointer(wlist(interface_idx), w)
        r = G1(j, x, w(x))
    end function

    function dot_product_beta(j1, j2, interface_idx) result(r)
        integer, intent(in) :: j1, j2, interface_idx
        double precision :: r

        r = integrate(f_aux, c_null_ptr, pts)
    contains
        function f_aux(x, args) result (r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            procedure(dw_proc_t), pointer :: dw

            call c_f_procpointer(dwlist(interface_idx), dw)
            r = sqrt(1 + dw(x)**2)*fbeta(j1, x, interface_idx)*fbeta(j2, x, interface_idx)
        end function
    end function

    function dot_product_gamma(j1, j2, interface_idx) result(r)
        integer, intent(in) :: j1, j2, interface_idx
        double precision :: r

        r = integrate(f_aux, c_null_ptr, pts)
    contains
        function f_aux(x, args) result (r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r
            procedure(dw_proc_t), pointer :: dw

            call c_f_procpointer(dwlist(interface_idx), dw)
            r = sqrt(1 + dw(x)**2)*fgamma(j1, x, interface_idx)*fgamma(j2, x, interface_idx)
        end function
    end function

    function ftransform(f, k, j, interface_idx) result(r)
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

        r = transform(f_aux, j, c_null_ptr, pts)
    contains
        function f_aux(x, args) result (r)
            double precision, intent(in) :: x
            type(c_ptr), intent(in) :: args
            double precision :: r

            r = f(k, x, interface_idx)
        end function
    end function

    function fa(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w

        call c_f_procpointer(wlist(interface_idx), w)
        if (j == 0) then
            r = (b - w(x)) / b
        else
            r = 2.0 * (sinh(mu(j) * (b - w(x))) / sinh(mu(j) * b)) * cos(mu(j)*x)
        end if
    end function

    function fb(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w

        call c_f_procpointer(wlist(interface_idx), w)
        if (j == 0) then
            r = - w(x) / b
        else
            r = -2.0 * (sinh(mu(j) * w(x)) / sinh(mu(j) * b)) * cos(mu(j)*x)
        end if
    end function

    function fc(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        integer :: m

        call c_f_procpointer(wlist(interface_idx), w)

        !        if (j == 0) then
        !            r = -a*w(x)/b
        !        else
        !            r = - a* (sinh(mu(j) * w(x)) / sinh(mu(j) * b)) * cos(mu(j)*x)
        !        end if


        r = -vpsi(j, 0)*w(x)/b
        do m = 1, mmax_F
            r = r - 2.0*vpsi(j, m)* (sinh(mu(m) * w(x)) / sinh(mu(m) * b)) * cos(mu(m)*x)
        end do
    end function


    function fp(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -k1 / b
        else
            r = -2.0 * k1 * mu(j)*(-dw(x) * (sinh(mu(j) * (b - w(x))) / sinh(mu(j) * b)) * sin(mu(j)*x) + &
                (cosh(mu(j) * (b - w(x))) / sinh(mu(j) * b)) * cos(mu(j)*x))
        end if
    end function


    function fq(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -k2 / b
        else
            r = - 2.0 * k2 * mu(j) * (dw(x) * (sinh(mu(j) * w(x)) / sinh(mu(j) * b)) * sin(mu(j)*x) + &
                (cosh(mu(j) * w(x)) / sinh(mu(j) * b)) * cos(mu(j)*x))
        end if
    end function

    function fr(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        integer :: m

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        !        if (j == 0) then
        !            r = -a*k1/b
        !        else
        !            r = - a *k1 * mu(j) * (dw(x) * (sinh(mu(j) * w(x)) / sinh(mu(j) * b)) * sin(mu(j)*x) + &
        !                (cosh(mu(j) * w(x)) / sinh(mu(j) * b)) * cos(mu(j)*x))
        !        end if

        r = -vpsi(j, 0) * k1 / b
        do m = 1, mmax_F
            r = r - 2.0 * vpsi(j, m) * k1 * mu(j) * (dw(x) * (sinh(mu(m) * w(x)) / sinh(mu(m) * b)) * sin(mu(m)*x) + &
                (cosh(mu(m) * w(x)) / sinh(mu(m) * b)) * cos(mu(m)*x))
        end do
    end function

    function fu(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -1.0 / b
        else
            r = -2.0 * mu(j)*(-dw(x) * (sinh(mu(j) * (b - w(x))) / sinh(mu(j) * b)) * sin(mu(j)*x) + &
                (cosh(mu(j) * (b - w(x))) / sinh(mu(j) * b)) * cos(mu(j)*x))
        end if
    end function

    function fv(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        integer :: m

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        !        if (j == 0) then
        !            r = -a/b
        !        else
        !            r = - a * mu(j) * (dw(x) * (sinh(mu(j) * w(x)) / sinh(mu(j) * b)) * sin(mu(j)*x) + &
        !                (cosh(mu(j) * w(x)) / sinh(mu(j) * b)) * cos(mu(j)*x))
        !        end if

        r = -vphi(j, 0)/ b
        do m = 1, mmax_G
            r = r - 2.0 * vphi(j, m) * mu(j) * (dw(x) * (sinh(mu(m) * w(x)) / sinh(mu(m) * b)) * sin(mu(m)*x) + &
                (cosh(mu(m) * w(x)) / sinh(mu(m) * b)) * cos(mu(m)*x))
        end do
    end function
end module reciprocity_functions_module
