      subroutine allocate_stressperm_3
c     only calculate for model 3, model 5, model 7 and model 8
c     initial setup calcs node neighbor information
      use comai
      use combi
      use comdi
      use comdti
      use comsi
      implicit none
c     
c     
      if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &     .not.allocated(rlxs))then  
         allocate(rlxs(n0))
         allocate(rlys(n0))
         allocate(rlzs(n0))
         allocate(drlxs(n0,4))
         allocate(drlys(n0,4))
         allocate(drlzs(n0,4))	  
         allocate (idum_str1(n0))
         idum_str1 = 0 
      endif
      return
      
      end
c......................................................................
