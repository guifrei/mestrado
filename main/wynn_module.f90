module wynn_module
    implicit none


contains

    function wynn(mmax, kmax, s) result(r)
        integer, intent(in) :: mmax
        integer, intent(in) :: kmax !kmax <= mmax + 1
        double precision :: r
        double precision, dimension(1:mmax, 1:kmax) :: eps
        integer :: k, j, cnt
        double precision, dimension(1:mmax), intent(in) :: s

        eps = 0.0

        eps(:,1) = 0.0
        eps(:,2) = s

        cnt = 1
        do k = 3, kmax
            do j = 1, mmax - cnt
                eps(j, k) = eps(j + 1, k - 2) + 1.0/(eps(j + 1, k - 1) - eps(j, k - 1))
            end do
            cnt = cnt + 1
        end do
        r = eps(1, kmax)
    end function
end module wynn_module
