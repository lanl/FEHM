      function harmonic_mean(a,b)

      implicit none
      real*8 a,b
      real*8 tol
      real*8 harmonic_mean
      parameter(tol = 1.d-12)

      harmonic_mean = 2.0d0*a*b/(a + b + tol)
      end function harmonic_mean

