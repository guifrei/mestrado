module tikhonov_module
    use constants_module
    use eigenfunctions_module
    implicit none

contains
    subroutine tikhonov(lambda, vx, vy, vvY)
        !http://www2.compute.dtu.dk/~pcha/DIP/chap4.pdf, chap8
        integer, parameter :: deriv_ord = 0
        double precision, intent(in) :: lambda
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(in) :: vy
        double precision, dimension(0: mmax_phi), intent(out) :: vvY
        integer, parameter :: mm = tnmax
        integer, parameter :: nn = mmax_phi + 1
        double precision, dimension(mm + nn - deriv_ord, nn) :: mxa
        double precision, dimension(mm + nn - deriv_ord, 1) :: vb
        integer :: info, i, j
        integer, parameter :: lwork = 64*nn
        double precision, dimension(lwork) :: work
        double precision :: x, y
        double precision, dimension(0: mmax_phi) :: integrals_Y
        double precision, dimension(nn - deriv_ord, nn) :: mL !operador derivativo

                    !            lambda = 0.0001

        mL = 0.0

        if (deriv_ord == 0) then
            do i = 1, nn - deriv_ord
                mL(i, i) = sqrt(lambda)
            end do
        else if (deriv_ord == 1) then
            do i = 1, nn - deriv_ord
                mL(i, i) = -sqrt(lambda)
                mL(i, i + 1) = sqrt(lambda)
            end do
        else if (deriv_ord == 2) then
            do i = 1, nn - deriv_ord
                mL(i, i) = sqrt(lambda)
                mL(i, i + 1) = -2.0*sqrt(lambda)
                mL(i, i + 2) = sqrt(lambda)
            end do
        end if

        mxa = 0.0
        vb = 0.0
        do i = 1, mm
            x = vx(i)
            do j = 1, nn
                mxa(i, j) = cos(mu(j - 1)*x)
            end do
        end do

        mxa(mm + 1: mm + nn - deriv_ord, :) = mL

        vb(1:mm, 1) = vy
        vb(mm + 1:mm + nn - deriv_ord, 1) = matmul(mL, vvY)

        call dgels('N', mm + nn - deriv_ord, nn, 1, mxa, mm + nn - deriv_ord, vb, mm + nn - deriv_ord, work, lwork, info)
        integrals_Y(0) = vb(1, 1)*a
        integrals_Y(1:nn-1) = vb(2:nn, 1)*a/2.0
        vvY = vb(1:nn, 1)

        open(unit = 5, file = '/home/cx3d/a.txt')
        do i = 1, mm
            x = vx(i)
            y = 0.0
            do j = 1, nn
                y = y + vb(j, 1)*cos(mu(j - 1)*x)
            end do
            write(5, *)x, vy(i), y
        end do
        close(5)
    end subroutine

    subroutine morozov(sigma, vx, vy, vvY, kmax)
        double precision, intent(in) :: sigma
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(in) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        integer, intent(out) :: kmax
        double precision, dimension(tnmax) :: tmpya, tmpyb
        double precision :: sqrt_rms_a, sqrt_rms_b, diff_a, diff_b
        logical :: keep

        keep = .true.
        kmax = 0
        tmpya = approximation_Y(vx, kmax)
        sqrt_rms_a = norm2(tmpya - vy)/sqrt(dble(tnmax))
        diff_a = dabs(sqrt_rms_a - sigma)
        kmax = 1

        do while (keep .and. (kmax <= mmax_phi))
            tmpyb = approximation_Y(vx, kmax)
            sqrt_rms_b = norm2(tmpyb - vy)/sqrt(dble(tnmax))
            diff_b = dabs(sqrt_rms_b - sigma)

            if (sqrt_rms_b <= sigma) then
                keep = .false.
                write(*, *)kmax, sqrt_rms_b
                !zerar os coeficientes deste ponto em diante ("filtrar as frequencias superiores")
                !vvY(kmax + 1:mmax_phi) = 0.0
            end if

            if (keep) then
                tmpya = tmpyb
                sqrt_rms_a = sqrt_rms_b
                diff_a = diff_b
                kmax = kmax + 1
            end if
        end do
    contains
        elemental function approximation_Y(x, nnmax) result(r)
            double precision, intent(in) :: x
            integer, intent(in) :: nnmax
            double precision :: r
            integer :: i

            r = 0.0
            do i = 0, nnmax
                r = r + vvY(i)*cos(mu(i)*x)
            end do
        end function
    end subroutine
end module tikhonov_module
