      subroutine stressperm_1(jpt)
c     perm model 1 - not implemented yet (default)
c     not needed unless fully coupled
      use comai
      use combi
      use comdi
      use comsi
      implicit none
      integer jpt

      if(allocated(rlxs)) then
         rlxs(jpt) = 1.0
         rlys(jpt) = 1.0 
         rlzs(jpt) = 1.0
         drlxs(jpt,1) = 0.0
         drlys(jpt,1) = 0.0
         drlzs(jpt,1) = 0.0
         drlxs(jpt,2) = 0.0
         drlys(jpt,2) = 0.0
         drlzs(jpt,2) = 0.0
         drlxs(jpt,3) = 0.0
         drlys(jpt,3) = 0.0
         drlzs(jpt,3) = 0.0
         drlxs(jpt,4) = 0.0
         drlys(jpt,4) = 0.0
         drlzs(jpt,4) = 0.0	  
      endif	 
      
      return
      end
c.....................................................................
