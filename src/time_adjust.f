       subroutine time_adjust(mmodel,modmin,time_type,time_cycle,
     * time,days0,day,days,daynew,maxtimes,iptty,iout,l,lchange,
     * node_model,steady_type,isty,n)
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
C***********************************************************************
CD1 
CD1 PURPOSE
CD1 
CD1 To adjust time if necessary for boundary models.  
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2
CD2 Initial implementation: NOV-96, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/time_adjust.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:18   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:12:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:44 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
C***********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 2.5.1 Implement time-step mechanism
CD3 
C***********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4 
C***********************************************************************
c
c adjust time if necessary for boundary models
c
c mmodel - number of models of boundary conditions
c modmax - model which forces time change
c time_type - type of changing time conditions
c time_cycle - cyclic time flow changes
c time   - array of flow rate change times
c days0  - last timestep total days
c days   - current timestep total days
c ichange - change flag
c

      implicit none

      integer maxtimes,mmodel,modmin,time_type(*)
      integer iout,iptty,l,lchange,i,int1,int2,icounter
      integer n,istea,isty
      integer node_model(*),steady_type(*)
      real*8 time_cycle(*)
      real*8 time(maxtimes,*),days0,day,days,daym,daymc,daynew
      real*8 daym0
c
       modmin=0
       int1=0
       int2=0
       lchange=0
       daym=1.d30
       daymc=1.d30
       do i=1,mmodel
c gaz steady state management 10-30-04  (added here 061605)
          if(isty.ne.0) then
             istea = steady_type(i)
          else
             istea = 1
          endif
          if(time_type(i).ne.0.and.istea.gt.0) then
             if(time_type(i).lt.0) then         
c  non cyclic time change
                icounter=abs(time_type(i)-1)
                if(days.ge.time(icounter,i)) then
                   daym0=time(icounter,i)-days0
                   if(daym0.lt.daym) then
                      daym=daym0
                      modmin=i
                   endif
                endif
             else
c  cyclic time change
c first check for repeat of cycle
                icounter=abs(time_type(i)+1)
                if(days.ge.time(icounter,i)) then
                   daym0=time(icounter,i)-days0
                   if(daym0.lt.daym) then
                      daym=daym0
                      modmin=i
                   endif
                endif
             endif
          endif
       enddo
c
c set day =daym, if necessary  
c 
       if(modmin.ne.0) then
          day=daym
          days=days0+day
          daynew=day
          lchange= l+1
       endif
       return
       end
