module constants_module
    implicit none

    double precision, parameter :: pi = 4.0*atan(1.0)
    double precision, parameter :: a = 0.04
    double precision, parameter :: b = 0.01
    double precision, parameter :: k1 = 54.0
    double precision, parameter :: k2 = 14.0
    double precision, parameter :: q = -7500.0
    double precision, parameter :: hmax = 400.0
    double precision, parameter :: gamma = 5000.0
    integer, parameter :: mmax_T = 35 !numero total de autofuncoes a serem calculadas para T1 e T2
    integer, parameter :: mmax_F = 20 !numero total de autofuncoes a serem calculadas para F1 e F2
    integer, parameter :: mmax_G = 20 !numero total de autofuncoes a serem calculadas para G1
    integer, parameter :: mmax_phi = max(mmax_T, mmax_F, mmax_G) !maior valor entre mmax_T, mmax_F e mmax_G
    integer, parameter :: delta_m = 60      !incremento no numero de termos
    double precision, parameter :: reltol = 0.0 !1.0D-8
    !    double precision, dimension(*), parameter :: pts = [0.0D0, a/5.0, a/4.0, 2.0*a/5.0, a/2.0, 3.0*a/5.0, 3.0*a/4.0, 4.0*a/5.0, a]
    double precision, dimension(*), parameter :: pts = [0.0D0, a/4.0, a/3.0, a/2.0, 2.0*a/3.0, 3.0*a/4.0, a]
    !        double precision, dimension(*), parameter :: pts = [0.0D0, a]


    double precision, parameter :: mp = 4.0
    double precision, parameter :: b1 = 0.011
    double precision, parameter :: b2 = 0.009

    integer, parameter :: tnmax = 121
    integer, parameter :: N = 20 ! tnmax/ndiv + 1 !min(mmax_F, mmax_G)           !numero total de funcionais de reciprocidade
    integer, parameter :: tmax = 1000          ! para o levantamento do campo de temperaturas
end module constants_module


!http://www.mathnet.or.kr/mathnet/thesis_file/JKMS-48-5-939-952.pdf
!https://www.mat.univie.ac.at/~neum/ms/regtutorial.pdf
!https://ac.els-cdn.com/S1674984715000956/1-s2.0-S1674984715000956-main.pdf?_tid=bec5ff7f-d01d-4b93-8c4e-45c9ff89ffa7&acdnat=1550706622_c2f56ce6e519463af8b8af700a736a2c
!https://web.stanford.edu/~takapoui/preconditioning.pdf

