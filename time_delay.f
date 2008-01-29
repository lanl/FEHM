      subroutine time_delay (transflag, ith,cur_node,
     2     par1v, par2v, par3v, fm,
     2     ret_factor, rseed, tau_zero, concv, timev)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To calculate particle time delay using type curves. 
CD1 (such as Sudicky and Frind)
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2..10 [10086-STN-2.10-00]
CD2 
CD2 Initial implementation: 05-MAY-99, Programmer: Zora Dash
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/time_delay.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:32:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:00   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
CD2
C**********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 2.3.5 Cell-based particle-tracking module
CD3 2.3.6 Streamline particle-tracking module
CD3 
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************

      use comai
      use compfrac
      implicit none

      real*8  par1v, par2v, par3v, tau_zero, timev, concv, timedelay
      real*8 interp
      real*8 ret_factor
      real ran_sp
      integer fm
      integer rseed
      integer transflag
      real inverf, fact_term, concvr
      integer ith, cur_node

c     if(transflag.eq.1.or.transflag.eq.2) then
! 01-Nov-06 Add flag for colloid diversity
      if( transflag.eq.1.or.transflag.eq.2.or.transflag.ge.11) then
         
         timev = par2v*tau_zero
      else

! Generate random concentration value
         concv = ran_sp (rseed)


         if(par1v.lt.0.) then
c     Error function solution - infinite fracture spacing

            concvr = concv
            timedelay = ret_factor + 
     2           par2v*tau_zero/(par1v**2*inverf(concvr)**2)
            timev = timedelay*tau_zero
         else
            if (.not. pfrac_read) then
               write (ierr, 10) 
               if (iout .ne. 0) write (iout, 10) 
               if (iptty .gt. 0) write (iptty, 10) 
               stop
            else
! Find time delay
               timedelay = interp (ith, cur_node,
     2              par1v, par2v, par3v, fm, concv)
! Compute adjusted time
               timev = tau_zero * (timedelay + ret_factor)
            end if
         end if
      end if
 10   format ("Dispersion type curve data not input")
      end

