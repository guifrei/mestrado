module morozov_module
    use constants_module
    use eigenfunctions_module
    implicit none

!https://www.statisticshowto.datasciencecentral.com/wp-content/uploads/2017/06/proof-of-tikhonov.pdf

contains
    subroutine morozov(sigma, vx, vy, vvY, i)
        double precision, intent(in) :: sigma
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(inout) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        integer, intent(out) :: i
        double precision, dimension(tnmax) :: tmpy
        double precision :: sqrt_rms
        logical :: keep

        if (sigma.ne.0) then
            keep = .true.
            tmpy = 0.0
            i = 0

            do while (keep .and. (i <= mmax_phi))
                tmpy = tmpy + vvY(i)*cos(mu(i)*vx)
                sqrt_rms = norm2(tmpy - vy)/sqrt(dble(tnmax))

                if (sqrt_rms <= sigma) then
                    keep = .false.
                    vy = tmpy
                    write(*, *, decimal='comma')i, sqrt_rms
                end if

                if (keep) then
                    i = i + 1
                end if
            end do
        end if
    end subroutine
end module morozov_module
