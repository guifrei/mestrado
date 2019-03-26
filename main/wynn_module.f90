module wynn_module
    implicit none

contains

    subroutine ok
        integer :: j, k
        double precision, dimension(0:4) :: s = [4.000, 2.667, 3.467, 2.895, 3.340]
        double precision, dimension(0:4) :: eps_1, eps, eps_new

        do k = 0,4
            if (k == 0) then
                eps_1 = 0
                eps = s
            else
                do j = 0, 3
                    eps_new(j) = eps_1(j + 1) + 1.0/(eps(j + 1) - eps(j))
                end do
                eps = eps_new
            end if
            write(*, *)eps
        end do
    end subroutine
end module wynn_module
