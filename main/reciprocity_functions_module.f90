module reciprocity_functions_module
    use constants_module
    use eigenfunctions_module
    use interfaces_module
    use netlib_module
    use iso_fortran_env
    use integration_module
    use auxiliary_functions_module
    implicit none

    double precision, dimension(0: 2*mmax_F + 1, 0: N), target :: coeffsF
    double precision, dimension(0: mmax_G, 0: N), target :: coeffsG
    double precision, dimension(0: N, 0: mmax_phi), target :: vpsi, vphi
    double precision, dimension(0: mmax_phi) :: vvY
    double precision, dimension(0: N) :: reciprocity_f, reciprocity_g

contains
    function fpsi(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        integer :: m

        !        r = vpsi(j, 0)/a
        !        do m = 1, mmax_phi
        !            r = r + (2.0/a)*vpsi(j, m)*cos(mu(m)*x)
        !        end do

        !        if (j == 0) then
        !            r = sqrt(1.0/a)
        !        else if (mod(j, 2) == 0) then
        !            r = sqrt(1.0/a)*sin(x*dble(j)*pi/a)
        !        else
        !            r = sqrt(1.0/a)*cos(x*dble(j + 1)*pi/a)
        !        end if

        if (j == 0) then
            r = 1.0
        else
            r = cos(x*dble(j)*pi/a)
        end if

    !        r = phi_eval(j, x)
    end function

    function fphi(j, x, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j, interface_idx
        double precision :: r
        integer :: m

        !        r = vphi(j, 0)/a
        !        do m = 1, mmax_phi
        !            r = r + (2.0/a)*vphi(j, m)*cos(mu(m)*x)
        !        end do

        !                if (j == 0) then
        !                    r = sqrt(1.0/a)
        !                else if (mod(j, 2) == 0) then
        !                    r = sqrt(1.0/a)*sin(x*dble(j)*pi/a)
        !                else
        !                    r = sqrt(1.0/a)*cos(x*dble(j + 1)*pi/a)
        !    end if

        if (j == 0) then
            r = 1.0
        else
            r = cos(x*dble(j)*pi/a)
        end if

    !        r = phi_eval(j, x)
    end function

    subroutine add_error(sample_y, stdev)
        double precision, dimension(tnmax), intent(inout) :: sample_y
        double precision, intent(in) :: stdev
        double precision, dimension(tnmax) :: rand_u, rand_v, rand_epsilon
        integer :: sz, i
        integer, allocatable :: seed(:)

        if (stdev > 0) then
            !            call random_seed(size = sz)
            !            allocate(seed(sz))
            !            call random_seed(get = seed)
            !            do i = 1, sz
            !                write (*, '(A5, I2, A4, I20)') 'seed(', i, ') = ', seed(i)
            !            end do
            !
            !            seed( 1) =           1735964997
            !            seed( 2) =          -1899570570
            !            seed( 3) =          -1916140827
            !            seed( 4) =           -365159221
            !            seed( 5) =          -1707816154
            !            seed( 6) =           1885880775
            !            seed( 7) =           1632974868
            !            seed( 8) =          -1732338435
            !            seed( 9) =          -2144761781
            !            seed(10) =           1558797291
            !            seed(11) =            -76752210
            !            seed(12) =          -1599751705
            !            seed(13) =          -1947760814
            !            seed(14) =           -650042213
            !            seed(15) =           -343692856
            !            seed(16) =            -75608390
            !            seed(17) =           -533548596
            !            seed(18) =             22080825
            !            seed(19) =          -1375472808
            !            seed(20) =           2125428121
            !            seed(21) =          -1221007268
            !            seed(22) =          -2022422129
            !            seed(23) =          -2100777368
            !            seed(24) =           1366608952
            !            seed(25) =            567304476
            !            seed(26) =          -1010752275
            !            seed(27) =           -737203679
            !            seed(28) =          -1255295320
            !            seed(29) =          -1979909672
            !            seed(30) =            998064812
            !            seed(31) =             18609693
            !            seed(32) =            -26961015
            !            seed(33) =                    0
            !
            !            call random_seed(put = seed)
            !            deallocate(seed)

            call random_seed
            rand_epsilon = 0.0
            do i = 1,1000
                call random_number(rand_u)
                call random_number(rand_v)
                rand_epsilon = rand_epsilon + cos(2.0*pi*rand_v)*sqrt(-2.0*log(rand_u))
            end do
            rand_epsilon = rand_epsilon/1000.0
            sample_y = sample_y + stdev*rand_epsilon
        end if
    end subroutine

    subroutine integrate_synthetic_temperatures(vx, sample_y, sz)
        integer, intent(in) :: sz
        double precision, dimension(sz), intent(in) :: vx
        double precision, dimension(sz), intent(in) :: sample_y
        double precision, dimension(sz) :: vy
        integer :: nn, ier
        double precision, dimension(sz) :: w
        double precision :: s, fp
        double precision, dimension(0: mmax_phi) :: integrals_Y
        integer :: k, j

        integer, parameter :: iopt = 0
        integer, parameter :: ord = 5
        integer :: nest
        integer :: lwrk
        integer :: lwrk_int
        double precision, dimension(:), allocatable :: t, c
        integer, dimension(:), allocatable :: iwrk
        double precision, dimension(:), allocatable :: wrk
        double precision, dimension(:), allocatable :: wrk_int

        nest = sz+ord+1
        lwrk = (sz*(ord+1)+nest*(7+3*ord))
        lwrk_int = sz + ord + 1

        w = 1.0

        allocate(t(nest), c(nest), iwrk(nest), wrk(lwrk), wrk_int(lwrk_int))

        do k = 0, mmax_phi
            do j = 1, sz
                vy(j) = sample_y(j)*cos(mu(k)*vx(j))
            end do

            s = 0.0

            do j = 1, ord+1
                t(j)=0.0
                t(j+ord+1)=a
            end do

            call curfit(iopt,sz,vx,vy,w,0.0D0,a,ord,s,nest,nn,t,c,fp, &
                wrk,lwrk,iwrk,ier)
            if (ier > 0) then
                write(*, *)'Error in curfit. Ier = ', ier
                stop
            end if

            integrals_Y(k) = splint(t,nn,c,ord,0.0D0,a,wrk_int)
        end do
        vvY(0) = integrals_Y(0)/a
        vvY(1:mmax_phi) = integrals_Y(1:mmax_phi)*2.0/a

        deallocate(t, c, iwrk, wrk, wrk_int)
    end subroutine

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
        integer :: nn, ier
        double precision, dimension(tnmax) :: w
        double precision :: s, fp
        double precision, dimension(0: mmax_phi) :: integrals_Y

        integer, parameter :: iopt = 0
        integer, parameter :: ord = 5
        integer, parameter :: nest = tnmax+ord+1
        integer, parameter :: lwrk = (tnmax*(ord+1)+nest*(7+3*ord))
        integer, parameter :: lwrk_int = tnmax + ord + 1
        double precision, dimension(nest) :: t, c
        integer, dimension(nest) :: iwrk
        double precision, dimension(lwrk) :: wrk
        double precision, dimension(lwrk_int) :: wrk_int

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

        call add_error(sample_y, stdev)

        w = 1.0

        do k = 0, mmax_phi
            do j = 1, tnmax
                vy(j) = sample_y(j)*cos(mu(k)*vx(j))
            end do

            s = 0.0

            do j = 1, ord+1
                t(j)=0.0
                t(j+ord+1)=a
            end do

            call curfit(iopt,tnmax,vx,vy,w,0.0D0,a,ord,s,nest,nn,t,c,fp, &
                wrk,lwrk,iwrk,ier)
            if (ier > 0) then
                write(*, *)'Error in curfit. Ier = ', ier
                stop
            end if

            integrals_Y(k) = splint(t,nn,c,ord,0.0D0,a,wrk_int)
        end do
        vvY(0) = integrals_Y(0)/a
        vvY(1:nn - 1) = integrals_Y(1:nn - 1)*2.0/a
    end subroutine

    subroutine least_squares_for_Y(vx, vy, vvY, interface_idx, condutance_idx, stdev_idx)
        use tikhonov_module
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(in) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        integer, intent(in), optional :: interface_idx
        integer, intent(in), optional :: condutance_idx
        integer, intent(in), optional :: stdev_idx
        integer, parameter :: mm = tnmax
        integer, parameter :: nn = mmax_phi + 1
        double precision, dimension(mm, nn) :: mxa
        double precision, dimension(mm, 1) :: vb
        integer :: info, i, j
        integer, parameter :: lwork = 64*nn
        double precision, dimension(lwork) :: work
        double precision :: x, y
        double precision, dimension(0: mmax_phi) :: integrals_Y
        character(len = 2) :: str_idx, str_cdx, str_stdev
        double precision, dimension(nn) :: sings
        double precision, dimension(mm) :: uu
        double precision, dimension(mm, nn) :: vt

        mxa = 0.0
        vb = 0.0
        do i = 1, mm
            x = vx(i)
            do j = 1, nn
                mxa(i, j) = cos(mu(j - 1)*x)
            end do
        end do

        vb(1:mm, 1) = vy

        !        call dgesvd('A', 'A', mm, nn, mxa, mm, sings, uu, mm, vt, nn, work, lwork, info)

        call dgels('N', mm, nn, 1, mxa, mm, vb, mm, work, lwork, info)
        integrals_Y(0) = vb(1, 1)*a
        integrals_Y(1:nn-1) = vb(2:nn, 1)*a/2.0
        vvY = vb(1:nn, 1)

        if (present(interface_idx)) write(str_idx, '(I2.2)') interface_idx
        if (present(condutance_idx)) write(str_cdx, '(I2.2)') condutance_idx
        if (present(stdev_idx)) then
            if (stdev_idx == 0) then
                str_stdev = '00'
            else if (stdev_idx == 1) then
                str_stdev = '01'
            else
                str_stdev = '05'
            end if
        end if

        if (present(interface_idx).and.present(condutance_idx).and.present(stdev_idx)) then
            open(unit = 5, file = '/home/cx3d/mestrado/data/coeficientes_interface_'//str_idx//'_conductance_'//str_cdx// &
                '_stdev_' // str_stdev // '.dat')
            do i = 0, nn - 1
                write(5, *)vvY(i)
            end do
            close(5)
        end if
    end subroutine

    function calc_reciprocity_f(j) result(r)
        integer, intent(in) :: j
        double precision :: r
        integer :: m
        double precision, dimension(:), pointer :: vA
        double precision :: arg

        vA(0:) => coeffsF(0::2, j)
        r = -q * vpsi(j, 0) / k1 + (vA(0) - vpsi(j, 0)) * vvY(0) / b
        do m = 1, mmax_F
            arg = mu(m) * b
            r = r + mu(m) * (vA(m) / sinh(arg) - vpsi(j, m) / tanh(arg)) * vvY(m)
        end do
    end function

    function calc_reciprocity_g(j) result(r)
        integer, intent(in) :: j
        double precision :: r
        integer :: m
        double precision, dimension(:), pointer :: vE
        double precision :: arg

        vE(0:) => coeffsG(0:, j)
        r = -q * vphi(j, 0) / k1 + (vE(0) - vphi(j, 0)) * vvY(0) / b
        do m = 1, mmax_G
            arg = mu(m) * b
            r = r + mu(m) * (vE(m) / sinh(arg) - vphi(j, m) / tanh(arg)) * vvY(m)
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

    function parcela_fluxo_calor(x, j, interface_idx) result(r)
        double precision, intent(in) :: x
        integer, intent(in) :: j
        integer, intent(in) :: interface_idx
        double precision :: r

        r = k1
        r = r*reciprocity_g(j)
        r = r*fgamma(j, x, interface_idx)
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

    subroutine calculate_reciprocity_coefficients(interface_idx)
        integer, intent(in) :: interface_idx
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
        character(1) :: equed, fact, trans
        double precision, dimension(0: 2*mmax_F + 1, 0: 2*mmax_F + 1) :: mxF
        double precision, dimension(0: mmax_G, 0: mmax_G) :: mxG
        character(len = 2) :: str_idx

        write(str_idx, '(I2.2)') interface_idx

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)


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

        !        do m = 0, N
        !            do j = 0, mmax_phi
        !                vpsi(m, j) = ftransform(fpsi, m, j, interface_idx)
        !                vphi(m, j) = ftransform(fphi, m, j, interface_idx)
        !            end do
        !        end do

        !write(*, *)'Generating matrices for F'

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

        !write(*, *)'Generating matrices for G'

        mxG = 0.0
        do j = 0, mmax_G
            do m = 0, mmax_G
                mxG(j,m) = ftransform(fu, m, j, interface_idx)
            end do
            do m = 0, N
                coeffsG(j, m) = ftransform(fv, m, j, interface_idx)
            end do
        end do

        ! Solucao do sistema
        !write(*, *)'Solving systems'
        fact = 'E'
        trans = 'N'
        equed = 'B'

        tmpcoeffsF = 0.0
        tmpcoeffsG = 0.0


        !https://www.sintef.no/globalassets/project/evitameeting/2005/lcurve.pdf
        !http://www2.compute.dtu.dk/~pcha/DIP/chap4.pdf
        !http://www2.compute.dtu.dk/~pcha/DIP/chap5.pdf

        call dgesvx(fact, trans, 2*mmax_F+2, N+1, mxF, 2*mmax_F+2, afF, 2*mmax_F+2, ipivF, equed, &
            rF, cF, coeffsF, 2*mmax_F+2, tmpcoeffsF, 2*mmax_F+2, rcond, ferr, berr, workF, iworkF, info)
        coeffsF = tmpcoeffsF

        call dgesvx(fact, trans, mmax_G+1, N+1, mxG, mmax_G+1, afG, mmax_G+1, ipivG, equed, &
            rG, cG, coeffsG, mmax_G+1, tmpcoeffsG, mmax_G+1, rcond, ferr, berr, workG, iworkG, info)
        coeffsG = tmpcoeffsG

        ! Algoritmo de ortogonalizacao de Gram-Schmidt
        !write(*, *)'Gram Schmidt'
        call gram_schmidt(interface_idx)
    end subroutine

    function soma_controle_erro(parcela, j, x, y, v, mmax) result(r)
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
    end function

    function parcela_F1(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = (va*(b - y) + vpsi(j, 0)*y)/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = (2.0/a)*(va*sinh(arg2 - arg1)/sinh(arg2) +&
                vpsi(j, m)*sinh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function parcela_dF1dx(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = 0.0
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = - mu(m)*(2.0/a)*(va*sinh(arg2 - arg1)/sinh(arg2) +&
                vpsi(j, m)*sinh(arg1)/sinh(arg2))*sin(mu(m)*x)
        end if
    end function

    function parcela_dF1dy(j, m, x, y, va) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, va
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = (vpsi(j, 0) - va)/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = mu(m)*(2.0/a)*(-va*cosh(arg2 - arg1)/sinh(arg2) +&
                vpsi(j, m)*cosh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function parcela_F2(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = vd*y/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = (2.0/a)*(vd*sinh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function parcela_dF2dx(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = 0.0
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = - mu(m)*(2.0/a)*(vd*sinh(arg1)/sinh(arg2))*sin(mu(m)*x)
        end if
    end function

    function parcela_dF2dy(j, m, x, y, vd) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, vd
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = vd/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = mu(m)*(2.0/a)*(vd*cosh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function parcela_G1(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = (ve*(b - y) + vphi(j, 0)*y)/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = (2.0/a)*(ve*sinh(arg2 - arg1)/sinh(arg2) +&
                vphi(j, m)*sinh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function parcela_dG1dx(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = 0.0
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = - mu(m)*(2.0/a)*(ve*sinh(arg2 - arg1)/sinh(arg2) +&
                vphi(j, m)*sinh(arg1)/sinh(arg2))*sin(mu(m)*x)
        end if
    end function

    function parcela_dG1dy(j, m, x, y, ve) result(r)
        integer, intent(in) :: j, m
        double precision, intent(in) :: x, y, ve
        double precision :: r
        double precision :: arg1, arg2

        if (m == 0) then
            r = (vpsi(j, 0) - ve)/(a*b)
        else
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = mu(m)*(2.0/a)*(-ve*cosh(arg2 - arg1)/sinh(arg2) + &
                vphi(j, m)*cosh(arg1)/sinh(arg2))*cos(mu(m)*x)
        end if
    end function

    function F1(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_F1, j, x, y, coeffsF(0::2, j), mmax_F)
    end function

    function dF1dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF1dx, j, x, y, coeffsF(0::2, j), mmax_F)
    end function

    function dF1dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF1dy, j, x, y, coeffsF(0::2, j), mmax_F)
    end function

    function F2(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_F2, j, x, y, coeffsF(1::2, j), mmax_F)
    end function

    function dF2dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF2dx, j, x, y, coeffsF(1::2, j), mmax_F)
    end function

    function dF2dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dF2dy, j, x, y, coeffsF(1::2, j), mmax_F)
    end function

    function G1(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_G1, j, x, y, coeffsG(0:, j), mmax_G)
    end function

    function dG1dx(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dG1dx, j, x, y, coeffsG(0:, j), mmax_G)
    end function

    function dG1dy(j, x, y) result(r)
        double precision, intent(in) :: x, y
        integer, intent(in) :: j
        double precision :: r

        r = soma_controle_erro(parcela_dG1dy, j, x, y, coeffsG(0:, j), mmax_G)
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
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        y = w(x)
        if (j == 0) then
            r = (b - y) / b
        else
            arg1 = mu(j)*y
            arg2 = mu(j)*b
            r = 2.0 * (sinh(arg2 - arg1) / sinh(arg2)) * cos(mu(j)*x)
        end if
    end function

    function fb(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        y = w(x)
        if (j == 0) then
            r = - y / b
        else
            arg1 = mu(j)*y
            arg2 = mu(j)*b
            r = -2.0 * (sinh(arg1) / sinh(arg2)) * cos(mu(j)*x)
        end if
    end function

    function fc(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        integer :: m
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)

        y = w(x)
        r = -vpsi(j, 0)*y/b
        do m = 1, mmax_F
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = r - 2.0*vpsi(j, m)* (sinh(arg1) / sinh(arg2)) * cos(mu(m)*x)
        end do
    end function


    function fp(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -k1 / b
        else
            y = w(x)
            arg1 = mu(j)*y
            arg2 = mu(j)*b
            r = -2.0 * k1 * mu(j)*(-dw(x) * (sinh(arg2 - arg1) / sinh(arg2)) * sin(mu(j)*x) + &
                (cosh(arg2 - arg1) / sinh(arg2)) * cos(mu(j)*x))
        end if
    end function


    function fq(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -k2 / b
        else
            y = w(x)
            arg1 = mu(j)*y
            arg2 = mu(j)*b
            r = - 2.0 * k2 * mu(j) * (dw(x) * (sinh(arg1) / sinh(arg2)) * sin(mu(j)*x) + &
                (cosh(arg1) / sinh(arg2)) * cos(mu(j)*x))
        end if
    end function

    function fr(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        integer :: m
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        r = -vpsi(j, 0) * k1 / b
        y = w(x)
        do m = 1, mmax_F
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = r - 2.0 * vpsi(j, m) * k1 * mu(j) * (dw(x) * (sinh(arg1) / sinh(arg2)) * sin(mu(m)*x) + &
                (cosh(arg1) / sinh(arg2)) * cos(mu(m)*x))
        end do
    end function

    function fu(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        if (j == 0) then
            r = -1.0 / b
        else
            y = w(x)
            arg1 = mu(j)*y
            arg2 = mu(j)*b
            r = -2.0 * mu(j)*(-dw(x) * (sinh(arg2 - arg1) / sinh(arg2)) * sin(mu(j)*x) + &
                (cosh(arg2 - arg1) / sinh(arg2)) * cos(mu(j)*x))
        end if
    end function

    function fv(j, x, interface_idx) result(r)
        integer, intent(in) :: j, interface_idx
        double precision, intent(in) :: x
        double precision :: r
        procedure(w_proc_t), pointer :: w
        procedure(dw_proc_t), pointer :: dw
        integer :: m
        double precision :: y, arg1, arg2

        call c_f_procpointer(wlist(interface_idx), w)
        call c_f_procpointer(dwlist(interface_idx), dw)

        y = w(x)
        r = -vphi(j, 0)/ b
        do m = 1, mmax_G
            arg1 = mu(m)*y
            arg2 = mu(m)*b
            r = r - 2.0 * vphi(j, m) * mu(j) * (dw(x) * (sinh(arg1) / sinh(arg2)) * sin(mu(m)*x) + &
                (cosh(arg1) / sinh(arg2)) * cos(mu(m)*x))
        end do
    end function
end module reciprocity_functions_module
