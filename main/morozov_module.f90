module morozov_module
    use constants_module
    use eigenfunctions_module
    implicit none

!https://www.statisticshowto.datasciencecentral.com/wp-content/uploads/2017/06/proof-of-tikhonov.pdf

contains
    subroutine morozov(sigma, vx, vy, vvY)
        double precision, intent(in) :: sigma
        double precision, dimension(tnmax), intent(in) :: vx
        double precision, dimension(tnmax), intent(inout) :: vy
        double precision, dimension(0: mmax_phi), intent(inout) :: vvY
        integer :: i
        double precision, dimension(tnmax) :: tmpy
        double precision :: sqrt_rms

        if (sigma.ne.0) then
            tmpy = 0.0
            do i = 0, mmax_phi
                tmpy = tmpy + vvY(i)*cos(mu(i)*vx)
                sqrt_rms = norm2(tmpy - vy)
                write(*, '(I3, F20.15)', decimal='comma')i, sqrt_rms !norm2(vvY(i)*cos(mu(i)*vx))
            end do
        end if
    end subroutine
end module morozov_module
