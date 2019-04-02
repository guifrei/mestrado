module algebraic_reconstruction_technique_module
    implicit none

contains

    subroutine kerp(m, mA, mB, mX, alpha, omega, imax)
        integer, intent(in) :: m
        double precision, dimension(m, m), intent(in) :: mA
        double precision, dimension(m), intent(inout) :: mB
        double precision, dimension(m), intent(inout) :: mX
        double precision, intent(in) :: alpha, omega
        integer, intent(in) :: imax
        double precision, dimension(m) :: mY
        integer :: k

        mX = 0.0
        mY = mB

        do k = 0, imax
            mY = capPhi(mY)
            mB = mB - mY
            mX = capF(mB, mX)
            write(*, *)mX
        end do

    contains
        function capF(b, x) result(r)
            double precision, dimension(m), intent(in) :: b, x
            double precision, dimension(m) :: r
            integer :: i

            r = fomega(b, x, M)
            do i = M - 1, 1, -1
                r = fomega(b, r, i)
            end do

        end function

        function capPhi(y) result(r)
            double precision, dimension(m), intent(in) :: y
            double precision, dimension(m) :: r
            integer :: j

            r = fphi(y, M)
            do j = M - 1, 1, -1
                r = fphi(r, j)
            end do

        end function

        function fomega(b, x, i) result(r)
            double precision, dimension(m), intent(in) :: b, x
            integer, intent(in) :: i
            double precision, dimension(m) :: r


            r = (1.0 - omega)*x + omega*f(b, x, i)
        end function

        function fphi(y, j) result(r)
            double precision, dimension(m), intent(in) :: y
            integer, intent(in) :: j
            double precision, dimension(m) :: r


            r = (1.0 - alpha)*y + alpha*phi(y, j)
        end function

        function f(b, x, i) result(r)
            double precision, dimension(m), intent(in) :: b, x
            integer, intent(in) :: i
            double precision, dimension(m) :: r

            r = x - (dot_product(x, mA(i, :)) - b(i))/(norm2(mA(i, :))**2)*mA(i, :)
        end function

        function phi(y, j) result(r)
            double precision, dimension(m), intent(in) :: y
            integer, intent(in) :: j
            double precision, dimension(m) :: r

            r = y - dot_product(y, mA(:, j))/(norm2(mA(:, j))**2)*mA(:, j)
        end function
    end subroutine

    subroutine art(m, mA, mB, mX, lambda, imax)
        integer, intent(in) :: m
        double precision, dimension(m, m), intent(in) :: mA
        double precision, dimension(m), intent(in) :: mB
        double precision, dimension(m), intent(inout) :: mX
        double precision, intent(in) :: lambda
        integer, intent(in) :: imax

        double precision, dimension(m):: rowA
        integer :: i, k
        double precision :: r1, r2, r3

        do k = 1, imax
            do i = 1, m
                rowA = mA(i, :)
                r1 = dot_product(rowA, mX)
                r2 = norm2(rowA)**2
                r3 = lambda*(mB(i) - r1)/r2
                mX = mX + r3*rowA
            !            write(*, *)mX
            end do
        end do

    end subroutine art
end module algebraic_reconstruction_technique_module
