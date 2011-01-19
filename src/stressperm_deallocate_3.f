      subroutine deallocate_stressperm_3
c     
c     deallocate memory for stress derivatives for fully coupled solution
     
      use comai
      use combi
      use comdi
      use comdti
      use comsi
      implicit none

      if(ipermstr3.ne.0.or.ipermstr5.ne.0.and.
     &     allocated(rlxs))then       
         deallocate(rlxs,rlys,rlzs)
         deallocate(drlxs,drlys,drlzs)
         deallocate(idum_str1)
      endif
      
      return
      
      end
c......................................................................
