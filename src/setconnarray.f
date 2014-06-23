      subroutine setconnarray
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Determine and store node connections at which interface reduction
!D1 factor model is applied.
!D1 
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 20-AUG-99, Programmer: B. Robinson  
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/setconnarray.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:28   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.4.12 Mass transport at interfaces 
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use combi
      use comdti
      use comai
      use comxi
      implicit none

      integer i, ii1, ii2, idg, neqp1, jmi, jmia, jitfc, kitfc
      integer jm, kb, n_levels, ilevel, add_node, add_istrw


      if(idpdp.eq.0) then
         n_levels = 1
      else
         n_levels = 2
      end if


c     Loop for either 1 or 2 continua

      do ilevel = 1, n_levels

         if(ilevel .eq. 1) then
            add_node = 0
            add_istrw = 0
         else
            add_node = neq
            add_istrw = nelm(neq+1)-neq-1
         end if

c     Section for determining the interface reduction
c     array

      do i = 1, neq
         
         ii1=nelm(i)+1
         ii2=nelm(i+1)
         idg=nelmdg(i)-ii1+1
         neqp1=neq+1
         jmi=nelmdg(i)
         jmia=jmi-neqp1
         
c         do jm=jmi+1,ii2
         do jm=ii1,ii2
            kb=nelm(jm)+add_node
            
c     loop over all possible reduction models to search for
c     reduction factor to be applied
            
          reduction: do jitfc = 1, nitfcpairs
            
            if(izonef_itfc(i+add_node).eq.zone_pair(jitfc,1)) then
               if(izonef_itfc(kb).eq.zone_pair(jitfc,2)) then
                  istrw_itfc(jm-neqp1+add_istrw) = jitfc
                  exit reduction
               end if
            elseif(izonef_itfc(i+add_node).eq.zone_pair(jitfc,2)) then
               if(izonef_itfc(kb).eq.zone_pair(jitfc,1)) then
                  istrw_itfc(jm-neqp1+add_istrw) = jitfc
                  exit reduction
               end if
            end if
          end do reduction
          
          colloid: do kitfc = 1, ncoldpairs
            
            if(izonef_itfc(i+add_node).eq.zonec_pair(kitfc,1)) then
               if(izonef_itfc(kb).eq.zonec_pair(kitfc,2)) then
                  istrw_cold(jm-neqp1+add_istrw) = kitfc
                  exit colloid
               end if
            elseif(izonef_itfc(i+add_node).eq.zonec_pair(kitfc,2)) then
               if(izonef_itfc(kb).eq.zonec_pair(kitfc,1)) then
                  istrw_cold(jm-neqp1+add_istrw) = kitfc
                  exit colloid
               end if
            end if
          end do colloid         
      end do
      
      end do

      end do

      return
      end
