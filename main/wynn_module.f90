module wynn_module
    implicit none


contains
    subroutine ok(nmax, s)
        integer, intent(inout) :: nmax
        double precision, dimension(1:nmax), intent(inout) :: s
        double precision, dimension(1:nmax - 2) :: ts
        integer :: k
        double precision :: r

        nmax = nmax - 2
        do k = 1, nmax
            r = s(k + 2) - ((s(k + 2) - s(k + 1))**2)/((s(k + 2) - s(k + 1)) - (s(k + 1) - s(k)))
            ts(k) = r
        end do

        s = ts
    end subroutine
end module wynn_module
