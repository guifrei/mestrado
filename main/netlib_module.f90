module netlib_module
    implicit none

    double precision :: outch, mcheps, dwarf
    common /rkcom7/  outch, mcheps, dwarf

    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            character(len = 1), intent(in) :: trans
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            integer, intent(in) :: lda
            integer, intent(in) :: ldb
            integer, intent(in) :: lwork
            integer, intent(out) :: info
            double precision, dimension(lda, n), intent(inout) :: a
            double precision, dimension(ldb, nrhs), intent(inout) :: b
            double precision, dimension(lwork), intent(out) :: work
        end subroutine
    end interface

    interface
        subroutine envirn(outch, mcheps, dwarf)
            import
            double precision, intent(out) :: outch, mcheps, dwarf
        end subroutine

        subroutine setup(neq, tstart, ystart, tend, tol, thres, method, task, errass, &
            hstart, work, lenwrk, mesage)
            import
            integer, intent(in) :: neq, method, lenwrk
            double precision, intent(in) :: tstart, tend, tol, hstart
            double precision, dimension(neq), intent(in) :: ystart, thres
            character(len = 1), intent(in) :: task
            logical, intent(in) :: errass, mesage
            double precision, dimension(lenwrk), intent(out) :: work
        end subroutine

        subroutine ut(f, twant, tgot, ygot, ypgot, ymax, work, uflag)
            import
            interface
                subroutine f(t, y, yp)
                    import
                    double precision, intent(in) :: t
                    double precision, dimension(*), intent(in) :: y
                    double precision, dimension(*), intent(out) :: yp
                end subroutine
            end interface
            double precision, intent(in) :: twant
            double precision, intent(out) :: tgot
            double precision, dimension(*), intent(out) :: ygot, ypgot, ymax
            double precision, dimension(*), intent(inout) :: work
            integer, intent(out) :: uflag
        end subroutine

        subroutine ct(f, tnow, ynow, ypnow, work, cflag)
            import
            interface
                subroutine f(t, y, yp)
                    import
                    double precision, intent(in) :: t
                    double precision, dimension(*), intent(in) :: y
                    double precision, dimension(*), intent(out) :: yp
                end subroutine
            end interface
            double precision, intent(in) :: tnow
            double precision, dimension(*), intent(out) :: ynow, ypnow
            double precision, dimension(*), intent(inout) :: work
            integer, intent(out) :: cflag
        end subroutine

        subroutine intrp(twant, reqest, nwant, ywant, ypwant, f, work, wrkint, lenint)
            import
            interface
                subroutine f(t, y, yp)
                    import
                    double precision, intent(in) :: t
                    double precision, dimension(*), intent(in) :: y
                    double precision, dimension(*), intent(out) :: yp
                end subroutine
            end interface
            double precision, intent(in) :: twant
            character(len = 1), intent(in) :: reqest
            integer, intent(in) :: nwant, lenint
            double precision, dimension(nwant), intent(out) :: ywant, ypwant
            double precision, dimension(*), intent(inout) :: work
            double precision, dimension(lenint), intent(out) :: wrkint
        end subroutine
    end interface

    interface
        function zeroin(ax,bx,f,tol) result (r)
            import
            double precision, intent(in) :: ax, bx, tol
            interface
                function f(x) result (r)
                    import
                    double precision, intent(in) :: x
                    double precision :: r
                end function
            end interface
            double precision :: r
        end function

        subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp, &
            wrk,lwrk,iwrk,ier)
            import
            integer, intent(in) :: iopt
            integer, intent(in) :: m
            double precision, dimension(m), intent(in) :: x
            double precision, dimension(m), intent(in) :: y
            double precision, dimension(m), intent(in) :: w
            double precision, intent(in) :: xb, xe
            integer, intent(in) :: k
            double precision, intent(in) :: s
            integer, intent(in) :: nest
            integer, intent(inout) :: n
            double precision, intent(inout), dimension(n) :: t
            double precision, intent(out), dimension(n) :: c
            double precision, intent(out) :: fp
            integer, intent(in) :: lwrk
            double precision, intent(out), dimension(m*(k+1)+nest*(7+3*k)) :: wrk
            integer, intent(inout), dimension(nest) :: iwrk
            integer, intent(out) :: ier
        end subroutine

        subroutine splev(t, n, c, k, x, y, m, ier)
            import
            integer, intent(in) :: n
            integer, intent(in) :: m
            integer, intent(in) :: k
            double precision, intent(in), dimension(n) :: t
            double precision, intent(in), dimension(n) :: c
            double precision, intent(in), dimension(m) :: x
            double precision, intent(out), dimension(m) :: y
            integer, intent(out) :: ier
        end subroutine

        function splint(t,n,c,k,a,b,wrk) result(r)
            import
            integer, intent(in) :: n
            integer, intent(in) :: k
            double precision, intent(in) :: a, b
            double precision, intent(in), dimension(n) :: t
            double precision, intent(in), dimension(n) :: c
            double precision, intent(out), dimension(n) :: wrk
            double precision :: r
        end function

        subroutine dqags(f,a,b,epsabs,epsrel,r,abserr,neval,ier,limit,lenw,last,iwork,work)
            import
            interface
                function f(x) result (r)
                    import
                    double precision, intent(in) :: x
                    double precision :: r
                end function
            end interface
            double precision, intent(in) :: a, b, epsabs,epsrel
            double precision, intent(out) :: r, abserr
            integer, intent(out) :: neval, ier
            integer, intent(in) :: lenw, last, limit
            integer, dimension(limit), intent(out) :: iwork
            double precision, dimension(lenw), intent(out) :: work
        end subroutine

        subroutine dqag(f,a,b,epsabs,epsrel,key,r,abserr,neval,ier, limit,lenw,last,iwork,work)
            import
            interface
                function f(x) result (r)
                    import
                    double precision, intent(in) :: x
                    double precision :: r
                end function
            end interface
            double precision, intent(in) :: a, b, epsabs,epsrel
            integer, intent(in) :: key
            double precision, intent(out) :: r, abserr
            integer, intent(out) :: neval, ier
            integer, intent(in) :: lenw, last, limit
            integer, dimension(limit), intent(out) :: iwork
            double precision, dimension(lenw), intent(out) :: work
        end subroutine

        subroutine dqng(f,a,b,epsabs,epsrel,r,abserr,neval,ier)
            import
            interface
                function f(x) result (r)
                    import
                    double precision, intent(in) :: x
                    double precision :: r
                end function
            end interface
            double precision, intent(in) :: a, b, epsabs,epsrel
            double precision, intent(out) :: r, abserr
            integer, intent(out) :: neval, ier
        end subroutine

        subroutine dqawo(f,a,b,omega,integr,epsabs,epsrel,r,abserr,neval,ier,leniw,maxp1,lenw,last,iwork,work)
            import
            interface
                function f(x) result (r)
                    import
                    double precision, intent(in) :: x
                    double precision :: r
                end function
            end interface
            double precision, intent(in) :: a, b, omega, epsabs,epsrel
            integer, intent(in) :: integr
            double precision, intent(out) :: r, abserr
            integer, intent(out) :: neval, ier
            integer, intent(in) :: leniw, maxp1, lenw, last
            integer, dimension(leniw), intent(out) :: iwork
            double precision, dimension(lenw), intent(out) :: work
        end subroutine

        subroutine dgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, &
            r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info)
            import
            character(1), intent(in) :: fact
            character(1), intent(in) :: trans
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            integer, intent(in) :: lda
            double precision, dimension(lda, n), intent(inout) :: a
            integer, intent(in) :: ldaf
            double precision, dimension(ldaf, n), intent(inout) :: af
            integer, dimension(n), intent(inout) :: ipiv
            character(1), intent(inout) :: equed
            double precision, dimension(n), intent(inout) :: r
            double precision, dimension(n), intent(inout) :: c
            integer, intent(in) :: ldb
            double precision, dimension(ldb, nrhs), intent(inout) :: b
            integer, intent(in) :: ldx
            double precision, dimension(ldx, nrhs), intent(inout) :: x
            double precision, intent(out) :: rcond
            double precision, dimension(nrhs), intent(inout) :: ferr
            double precision, dimension(nrhs), intent(inout) :: berr
            double precision, dimension(4*n), intent(out) :: work
            integer, dimension(n), intent(out) :: iwork
            integer, intent(out) :: info
        end subroutine
    end interface

    interface
        subroutine extrap(init, eps, maxcol, sk, ginit, result, info)
            integer, intent(inout) :: init
            double precision, intent(in) :: eps, sk
            integer, intent(in) :: maxcol
            double precision, dimension(2*(maxcol + 1)), intent(in) :: ginit
            double precision, intent(out) :: result
            integer, intent(out) :: info
        end subroutine
    end interface
end module netlib_module
