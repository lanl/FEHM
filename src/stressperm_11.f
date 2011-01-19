      subroutine stressperm_11(i)
c     s kelkar June 9 2009, simple gangi (1978, Int.J.Rock Mech)
c     model in 2-D
c     assumed model in x-y plane and fracture in x-z plane, so only 
c     x-perm is modified. This formulation is for fracture faces
c     in contact, so str_y = effective stress in y-dir has to be
c     compressive

      use comai
      use combi
      use comdi
      use comsi

      implicit none
      integer i, iispmd
      real*8 gk0, gpmod, gmexp, gn, sigy_eff

      iispmd = ispm(i)    

      if(icnl.ne.0) then
         gk0  =spm1f(iispmd)
         gpmod=spm2f(iispmd)
         gmexp =spm3f(iispmd)
         gn=1./gmexp
         sigy_eff=str_y(i)
         if(sigy_eff.gt.0.d0) then
            pnx(i)=gk0*(1.-(sigy_eff/gpmod)**gmexp)**3.
            e2(i) =  gn*sigy_eff*(gpmod/sigy_eff)**gn
         endif
      endif
      
      return
      
      end
c.....................................................................
