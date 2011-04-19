      real*8 function interp (ith, cur_node,
     2     par1v, par2v, par3v, fm, concv) 
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To calculate particle time delay using type curves.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 05-MAY-99, Programmer: Zora Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/interp.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:20   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:04   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
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
! Return the interpolated time value

      use compfrac      
      implicit none

      real*8 par1v, par2v, par3v, concv, timev, interp_bilin
      real*8 timec(8), dltime, ftime, interp2, interp4, interp8
      integer fm, ilow, ihigh, jlow, jhigh, klow, khigh
      integer ith, cur_node
      real weight(4)
      integer points(4)

      sigma_low(2) = min (sigma_low(2), par1v)
      sigma_high(2) = max (sigma_high(2), par1v)
      omega_low(2) = min (omega_low(2), par2v)
      omega_high(2) = max (omega_high(2), par2v)
      par3_low(2) = min (par3_low(2), par3v)
      par3_high(2) = max (par3_high(2), par3v)
     
      if(curve_structure.eq.1) then

c     Free format structure - find closest curve and set indices
c     so that the code uses the exact curve w/o interpolation
c     But only do so if the curve has not already been found

         if(itf_curve(ith,cur_node,1).eq.0) then
         
            call find_closest_curve(0,par1v, par2v, par3v, 
     2           ilow, jlow, klow,points,weight)
            itf_curve(ith,cur_node,1) = ilow
            ihigh = ilow
            jhigh = jlow
            khigh = klow
         else
            ilow = itf_curve(ith,cur_node,1)
            ihigh = itf_curve(ith,cur_node,1)
            jlow = 1
            jhigh = 1
            klow = 1
            khigh = 1
         end if
         if(ipartout.eq.1) then
            param_density(itf_curve(ith,cur_node,1),1,1) = 
     2           param_density(itf_curve(ith,cur_node,1),1,1) + 1
         end if
         dltime = ftime(itf_curve(ith,cur_node,1),
     2        1, 1, fm, concv)
         if(numparams.le.2) then
            timev = 10.0 ** dltime 
         else
            timev = dltime
         end if
         interp = timev
         return

      elseif(curve_structure.gt.1) then
         if(itf_curve(ith,cur_node,1).eq.0) then
!            write(6,*) 'ith,cur_node',ith,cur_node
         
            call find_closest_curve(1,par1v, par2v, par3v, 
     2           ilow, jlow, klow,points,weight)
!            write(6,*) 'ith,cur_node after',ith,cur_node
            itf_curve(ith,cur_node,1) = points(1)
            itf_curve(ith,cur_node,2) = points(2)
            itf_curve(ith,cur_node,3) = points(3)
            wt_curve(ith,cur_node,1) = weight(1)
            wt_curve(ith,cur_node,2) = weight(2)
            wt_curve(ith,cur_node,3) = weight(3)
!            write(6,*) 'weight', weight(1), weight(2), weight(3)
!            write(6,*) 'points', points(1), points(2), points(3)
            if(numparams.gt.2) then
               itf_curve(ith,cur_node,4) = points(4)
               wt_curve(ith,cur_node,4) = weight(4)
            end if
         end if
         dltime = 0.
         dltime = dltime + wt_curve(ith,cur_node,1)*
     2        ftime(itf_curve(ith,cur_node,1),
     3        1, 1, fm, concv)
         dltime = dltime + wt_curve(ith,cur_node,2)*
     2        ftime(itf_curve(ith,cur_node,2),
     3        1, 1, fm, concv)
         dltime = dltime + wt_curve(ith,cur_node,3)*
     2        ftime(itf_curve(ith,cur_node,3),
     3        1, 1, fm, concv)
         if(numparams.gt.2) then
            dltime = dltime + wt_curve(ith,cur_node,4)*
     2           ftime(itf_curve(ith,cur_node,4),
     3           1, 1, fm, concv)
         end if
         if(ipartout.eq.1) then
            param_density(itf_curve(ith,cur_node,1),1,1) = 
     2           param_density(itf_curve(ith,cur_node,1),1,1) + 1
         end if

c     Added check to ensure large or small dltime values are not being
c     found by interpolation svd scheme. If they are, revert
c     to nearest neighbor scheme

         if(numparams.le.2) then
            if(abs(dltime).gt.30.) then
               dltime = ftime(itf_curve(ith,cur_node,1),
     3              1, 1, fm, concv)
            end if
            timev = 10.0 ** dltime
            interp = timev
            return
         else
            interp = dltime
            return
         end if


      else
c     Regular structure to curves - do simple interpolation
C     Find the T vs C curves indices
         call indices (ilow, ihigh, par1v, param1, nump1)
         call indices (jlow, jhigh, par2v, param2, nump2)
         if (nump3 .gt. 1) then
            call indices (klow, khigh, par3v, param3, nump3)
         else
            klow = 1
            khigh = 1
         end if
         if(ipartout.eq.1) then
            param_density(ilow,jlow,klow) = 
     2           param_density(ilow,jlow,klow) + 1
         end if
      end if
C Find time values
      if (ilow .eq. ihigh) then
         if (jlow .eq. jhigh) then
            if (klow .eq. khigh) then
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               if(numparams.le.2) then
                  timev = 10.0 ** timec(1)
               else
                  timev = timec(1)
               end if
            else 
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jlow, khigh, fm, concv)
                  dltime = interp2 (param3(klow), param3(khigh), par3v,
     .                 timec(1), timec(2))
               if(numparams.le.2) then
                  timev = 10.0 ** dltime
               else
                  timev = dltime   
               end if
            end if
         else
            if (klow .eq. khigh) then
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jhigh, klow, fm, concv)
               dltime = interp2 (param2(jlow), param2(jhigh), par2v,
     .              timec(1), timec(2))
               if(numparams.le.2) then
                  timev = 10.0 ** dltime
               else
                  timev = dltime
               end if
            else 
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jlow, khigh, fm, concv)
               timec(3) = ftime (ilow, jhigh, klow, fm, concv)
               timec(4) = ftime (ilow, jhigh, khigh, fm, concv)
               dltime = interp4 (param2(jlow), param2(jhigh), par2v,
     .              param3(klow), param3(khigh), par3v, 
     .              timec(1), timec(2), timec(3), timec(4))
               if(numparams.le.2) then
                  timev = 10.0 ** dltime 
               else
                  timev = dltime
               end if
            end if
         end if
      else
         if (jlow .eq. jhigh) then
            if (klow .eq. khigh) then
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ihigh, jlow, klow, fm, concv)
               dltime = interp2 (param1(ilow), param1(ihigh), par1v,
     .              timec(1), timec(2))
               if(numparams.le.2) then
                  timev = 10.0 ** dltime 
               else
                  timev = dltime
               end if

            else
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jlow, khigh, fm, concv)
               timec(3) = ftime (ihigh, jlow, klow, fm, concv)
               timec(4) = ftime (ihigh, jlow, khigh, fm, concv)
               dltime = interp4 (param1(ilow), param1(ihigh), par1v,
     .              param3(klow), param3(khigh), par3v, 
     .              timec(1), timec(2), timec(3), timec(4))

               if(numparams.le.2) then
                  timev = 10.0 ** dltime 
               else
                  timev = dltime
               end if
            end if
         else 
            if (klow .eq. khigh) then
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jhigh, klow, fm, concv)
               timec(3) = ftime (ihigh, jlow, klow, fm, concv)
               timec(4) = ftime (ihigh, jhigh, klow, fm, concv)
               if (curve_structure .eq. 0) then
                  dltime = interp_bilin (param1(ilow), param1(ihigh), 
     .                 par1v, param2(jlow), param2(jhigh), par2v, 
     .                 timec(1), timec(2), timec(3), timec(4))
               else
                  dltime = interp4 (param1(ilow), param1(ihigh), par1v,
     .                 param2(jlow), param2(jhigh), par2v, 
     .                 timec(1), timec(2), timec(3), timec(4))
               end if
               if(numparams.le.2) then
                  timev = 10.0 ** dltime 
               else
                  timev = dltime
               end if

            else
               timec(1) = ftime (ilow, jlow, klow, fm, concv)
               timec(2) = ftime (ilow, jlow, khigh, fm, concv)
               timec(3) = ftime (ilow, jhigh, klow, fm, concv)
               timec(4) = ftime (ilow, jhigh, khigh, fm, concv)
               timec(5) = ftime (ihigh, jlow, klow, fm, concv)
               timec(6) = ftime (ihigh, jlow, khigh, fm, concv)
               timec(7) = ftime (ihigh, jhigh, klow, fm, concv)
               timec(8) = ftime (ihigh, jhigh, khigh, fm, concv)
               dltime = interp8 (param1(ilow), param1(ihigh), par1v,
     .              param2(jlow), param2(jhigh), par2v, 
     .              param3(klow), param3(khigh), par3v, 
     .              timec(1), timec(2), timec(3), timec(4),
     .              timec(5), timec(6), timec(7), timec(8))
               if(numparams.le.2) then
                  timev = 10.0 ** dltime 
               else
                  timev = dltime
               end if

            end if
         end if
      end if
      interp = timev
      return
      end
