      subroutine srange( l, v, i, vmin, vmax)

      implicit NONE

      double precision v(*)

      integer          i, j, k, l

      double precision temp, vmin, vmax

c----------------------------------------------------------------------------                    

      temp  = v(1)
      vmin  = temp
      vmax  = temp

      if (l .eq. 1) return
      if (i .eq. 1) then
        do j = 2, l
          temp = v(j)
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
        end do
      else
        k = 1 + i
        do j = 2, l
          temp = v(k)
          vmin = min(vmin,temp)
          vmax = max(vmax,temp)
          k    = k + i
          end do
      end if

      return
      end

      subroutine lgd1v  ( x, pro, mu, sigsq, n, G, z, hood)

      implicit NONE

      integer            n, G

c     double precision   x(n), pro(G), mu(G), sigsq(G), z(n,G)
      double precision   x(*), pro(*), mu(*), sigsq(*), z(n,*)

      double precision   hood(*)

      integer                 i, k

      double precision        temp, const, tmax, sum, term
      double precision        muk, sigsqk, prok, sigmin, siglog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call srange( G, sigsq, 1, sigmin, temp)  

      if (sigmin .le. zero) then
        do i = 1, n       
          hood(i) = FLMAX
        end do
        return
      end if

      do k = 1, G
        muk    = mu(k)
        sigsqk = sigsq(k)
        siglog = log(sigsqk)
        const  = pi2log + siglog
        do i = 1, n
          temp = abs(x(i) - muk)
          if (temp .ne. zero) then
            term   = two*log(temp) - siglog
            z(i,k) = -(const+exp(term))/two
          else
            z(i,k) = -const/two
          end if 
        end do
      end do

      do i = 1, n
        tmax = -FLMAX
        do k = 1, G
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            z(i,k) = temp
            tmax   = max(tmax,temp)
          end if
        end do
        sum = zero
        do k = 1, G
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood(i) = log(sum)+tmax
      end do

      return
      end

      subroutine lgdvii ( x, pro, mu, sigsq, n, p, G, z, hood)

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G)
      double precision   x(n,*), z(n,*)

c     double precision   pro(G), mu(p,G), sigsq(G)
      double precision   pro(G), mu(p,*), sigsq(*)

c     double precision   hood(max(n,p))
      double precision   hood(*)

      integer                 i, j, k

      double precision        sum, temp, const, tmax
      double precision        prok, sigsqk, sigmin, siglog

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------

      call srange( G, sigsq, 1, sigmin, temp)    

      if (sigmin .le. zero) then
        do i = 1, n       
          hood(i) = FLMAX
        end do
        return
      end if

      do k = 1, G
        sigsqk = sigsq(k)
        siglog = log(sigsqk)
        const  = dble(p)*(pi2log+siglog)
        do i = 1, n
          tmax = zero
          do j = 1, p
            temp    = abs(x(i,j) - mu(j,k))
            tmax    = max( tmax, temp)
            hood(j) = temp
          end do
          sum  = zero
          if (tmax .gt. one) then
            do j = 1, p
              hood(j) = hood(j)/tmax
              sum = sum + hood(j)*hood(j)
            end do
          else
            tmax = one
            do j = 1, p
              sum = sum + hood(j)*hood(j)
            end do
          end if  
          temp   = (log(sum) + two*log(tmax)) - siglog
          z(i,k) = -(const+exp(temp))/two
        end do
      end do

      do i = 1, n
        tmax = -FLMAX
        do k = 1, G
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            z(i,k) = temp
            tmax   = max(tmax,temp)
          end if
        end do
        sum   = zero
        do k = 1, G
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood(i) = log(sum)+tmax
      end do

      return
      end 

      subroutine lgdvvi ( x, pro, mu, scale, shape, n, p, G, z, hood)

      implicit NONE

      integer            n, p, G

c     double precision   x(n,p), z(n,G), hood(max(n,p))       
      double precision   x(n,*), z(n,*), hood(*)

c     double precision   pro(G), mu(p,G), scale(G), shape(p,G)               
      double precision   pro(*), mu(p,*), scale(*), shape(p,*)

      integer                 i, j, k, nz

      double precision        sum, temp, const, tmax
      double precision        smin, smax, prok, scalek

      double precision        zero, one, two
      parameter              (zero = 0.d0, one = 1.d0, two = 2.d0)

      double precision        pi2log
      parameter              (pi2log = 1.837877066409345d0)

      double precision        FLMAX
      parameter              (FLMAX = 1.7976931348623157d308)

      double precision        RTMAX
      parameter              (RTMAX = 1.340780792994260d154)

      double precision        RTMIN
      parameter              (RTMIN = 1.49166814624d-154)

      double precision        SMALOG
      parameter              (SMALOG = -708.d0)

c------------------------------------------------------------------------------   
      call srange( G, scale, 1, smin, smax)

      if (smin .le. zero) then
        do i = 1, n       
          hood(i) = FLMAX
        end do
        return
      end if

      do k = 1, G
        call srange( p, shape(1,k), 1, smin, smax)
        if (smin .le. zero) then
          do i = 1, n       
            hood(i) = FLMAX
          end do
          return
        end if
        temp = log(scale(k))
        do j = 1, p
          shape(j,k) = (temp + log(shape(j,k)))/two
        end do
      end do

      do k = 1, G
        scalek = scale(k)
        const  = dble(p)*(pi2log+log(scalek))
        do i = 1, n
          sum  = zero
          do j = 1, p
            temp = abs(x(i,j) - mu(j,k))
c           temp = temp/shape(j,k)
            if (temp .ne. zero) then
              temp = exp(log(temp) - shape(j,k))
              sum  = sum + temp*temp
            end if
          end do
          z(i,k) = -(const+sum)/two
        end do
      end do

      do i = 1, n
        tmax = -FLMAX
        do k = 1, G
          prok = pro(k)
          if (prok .eq. zero) then
            z(i,k) = zero
          else
            temp   = log(prok) + z(i,k)
            z(i,k) = temp
            tmax   = max(tmax,temp)
          end if
        end do
        sum   = zero
        do k = 1, G
          if (pro(k) .ne. zero) then
            temp =  z(i,k) - tmax
            if (temp .ge. SMALOG) then
              z(i,k) = exp(temp)
              sum    = sum + z(i,k)
            else
              z(i,k) = zero
            end if
          end if
        end do
        hood(i) = log(sum)+tmax
      end do

      return
      end

