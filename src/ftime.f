      real *8 function ftime (iparam1, jparam2, kparam3, fm, concv)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************
!D1 
!D1 PURPOSE
!D1 
!D1 To calculate particle time delay using type curves.
!D1 
!***********************************************************************
!D2 
!D2 REVISION HISTORY
!D2 
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 Initial implementation: 28-Oct-02, Programmer: Z. Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/ftime.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use compfrac
      implicit none
      real*8 concv, cadj, ctime, concv_sc
      integer fm, iparam1, jparam2, kparam3, low, high, midpoint, numv
      
      numv = nump (iparam1, jparam2, kparam3, fm)
      if (numparams .le. 2) then
! Don't scale for Sudicky-Frind curves
         concv_sc = concv
      else
         concv_sc = concv*conc(iparam1, jparam2, kparam3, fm, numv)
      endif
      if (concv_sc .le. conc(iparam1, jparam2, kparam3, fm, 1)) then
         low = 1
         high = 1
         ftime = dtime(iparam1, jparam2, kparam3, fm, 1)
         if(numparams.le.2) then
            ftime = dlog10(ftime)
         end if
         return
      else if (concv_sc .ge. conc(iparam1, jparam2, kparam3, fm, numv)) 
     .        then
         low = numv
         high = numv
         ftime = dtime(iparam1, jparam2, kparam3, fm, numv)
         if(numparams.le.2) then
            ftime = dlog10(ftime)
         end if
         return
      else
         low = 1
         high = numv
 20      midpoint = (low + high) / 2
         if (concv_sc .eq. 
     .        conc(iparam1, jparam2, kparam3, fm, midpoint)) 
     .        then
            low = midpoint
            high = low
            ftime = dtime(iparam1, jparam2, kparam3, fm, 
     .           midpoint)
            if(numparams.le.2) then
               ftime = dlog10(ftime)
            end if
            return
         else if (concv_sc .lt. 
     .           conc(iparam1, jparam2, kparam3, fm, midpoint)) then
            high = midpoint
         else
            low = midpoint
         end if
         if (high - low .eq. 1) goto 30
         goto 20
      endif

 30   continue
      cadj = (concv_sc - conc(iparam1, jparam2, kparam3, fm, low)) / 
     .     (conc(iparam1, jparam2, kparam3, fm, high) - 
     .     conc(iparam1, jparam2, kparam3, fm, low))
      if(numparams.le.2) then
      ctime = dlog10 (dtime(iparam1, jparam2, kparam3, fm, low)) + 
     .     cadj * (dlog10 (dtime(iparam1, jparam2, kparam3, fm, high)) -
     .     dlog10 (dtime(iparam1, jparam2, kparam3, fm, low)))
      ftime = ctime
      else
      ctime = dtime(iparam1, jparam2, kparam3, fm, low) + 
     .     cadj * (dtime(iparam1, jparam2, kparam3, fm, high) -
     .     dtime(iparam1, jparam2, kparam3, fm, low))
      ftime = ctime
      end if

      return
      end

