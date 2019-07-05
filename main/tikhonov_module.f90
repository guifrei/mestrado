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

    subroutine morozov(sigma, vx, vy, vvY)
        double precision, intent(in) :: sigma
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(in) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        double precision, dimension(tnmax) :: tmpy
        double precision :: sqrt_rms
        integer :: nnmax, kk
        logical :: keep

        keep = .true.
        nnmax = 0

        do while (keep .and. (nnmax <= 4))!mmax_phi
            tmpy = approximation_Y(vx)
            sqrt_rms = norm2(tmpy - vy)/sqrt(dble(tnmax))

            if (sqrt_rms <= sigma) then
                keep = .false.
                write(*, *)nnmax, sqrt_rms
                !zerar os coeficientes deste ponto em diante ("filtrar as frequencias superiores")
                vvY(nnmax + 1:mmax_phi) = 0.0
            else
                nnmax = nnmax + 1
            end if
        end do


    contains
        elemental function approximation_Y(x) result(r)
            double precision, intent(in) :: x
            double precision :: r
            integer :: i

            r = 0.0
            do i = 0, nnmax
                r = r + vvY(i)*cos(mu(i)*x)
            end do
        end function
    end subroutine
end module tikhonov_module
