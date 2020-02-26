module tikhonov_module
    use constants_module
    use eigenfunctions_module
    implicit none

contains
    subroutine tikhonov2(lambda, vx, vy, vvY)
        double precision, intent(in) :: lambda
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(in) :: vy
        double precision, dimension(0: mmax_phi), intent(out) :: vvY
        double precision, dimension(tnmax, 0: mmax_phi) :: mA
        double precision, dimension(0: mmax_phi, 0: mmax_phi) :: mAtA
        integer :: i, j

        integer, dimension(mmax_phi + 1) :: ipiv
        double precision, dimension(:), allocatable :: work
        double precision, dimension(1) :: swork
        integer :: lwork
        integer :: info

        do j = 0, mmax_phi
            mA(:, j) = cos(mu(j)*vx)
        end do

        mAtA = matmul(transpose(mA), mA)
        do i = 0, mmax_phi
            mAtA(i, i) = mAtA(i, i) + lambda*lambda
        end do

        vvY = matmul(transpose(mA), vy)
        lwork = -1
        call dsysv('U', mmax_phi + 1, 1, mAtA, mmax_phi + 1, ipiv, vvY, mmax_phi + 1, swork, lwork, info)
        lwork = int(swork(1))
        allocate(work(lwork))
        call dsysv('U', mmax_phi + 1, 1, mAtA, mmax_phi + 1, ipiv, vvY, mmax_phi + 1, work, lwork, info)
    end subroutine

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

        mL = 0.0

        if (deriv_ord == 0) then
            do i = 1, nn - deriv_ord
                mL(i, i) = lambda
            end do
        else if (deriv_ord == 1) then
            do i = 1, nn - deriv_ord
                mL(i, i) = -lambda
                mL(i, i + 1) = lambda
            end do
        else if (deriv_ord == 2) then
            do i = 1, nn - deriv_ord
                mL(i, i) = lambda
                mL(i, i + 1) = -2.0*lambda
                mL(i, i + 2) = lambda
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
        vb(mm + 1:mm + nn - deriv_ord, 1) = 0!matmul(mL, vvY)

        call dgels('N', mm + nn - deriv_ord, nn, 1, mxa, mm + nn - deriv_ord, vb, mm + nn - deriv_ord, work, lwork, info)
        integrals_Y(0) = vb(1, 1)*a
        integrals_Y(1:nn-1) = vb(2:nn, 1)*a/2.0
        vvY = vb(1:nn, 1)
    end subroutine
end module tikhonov_module
