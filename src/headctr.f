      subroutine headctr(iflg,ii,pres_value,head_value)
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
CD1  PURPOSE
CD1
CD1  This subroutine manages the head input and changes pressures
CD1  into heads and vice versa           
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 Initial implementation: 6-FEB-97, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/headctr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:14   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:07:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:09:48   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:28   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:02:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4  Finite Element Ouput Format
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************
c this routine is called from subroutine input(iflg=0)
c this routine is called from subroutine *****(iflg=1)
c this routine does not call any other routines
c------------------------------------------------------------------

      use comdti
      use comai
      use combi
      use comdi
      use comgi
      use comrxni
      use comfi
      use comii
      use comwt
      implicit none

C     integer iflg,icode
      integer iflg, neq2
      integer mi,mid,ii,iseth
      real*8 cord_val
      real*8 pres_value
      real*8 head_value
      save iseth

c================================================================
      if(ihead.eq.0) return
c================================================================
      neq2 = 2*neq

      if(iflg.eq.0) then
         if(.not.allocated(head)) allocate(head(n0))
         if(grav.le.0.0) then
            iseth = 1
         else
            iseth = 0
         endif
         icons=1
         ihead=1       
      else if(iflg.eq.1) then
c     
c     check for campatibility
c     
         if(ico2.ge.0) then
            if (iout .ne. 0) write(iout,10) 
            if (iptty. ne. 0) write(iptty,10) 
 10         format(/,'**** head macro without air macro :stopping')
            stop
         endif
c     
      else if(iflg.eq.2) then
c     
c     convert from pressure to head(all nodes)
c     
         do mi=1,n
            if(mi.gt.neq2) then
               mid=mi-neq2        
            elseif(mi.gt.neq) then
               mid=mi-neq        
            else
               mid=mi        
            endif
            if(iseth.eq.1) then
               cord_val = 0.0
            else
               cord_val = cord(mid,igrav) 
            endif
c gaz 090119
          if(ichead.eq.0) then                 
            head(mi)=(pho(mi)-pref)/rho1grav + cord_val
     &           -head0
c gaz 112619 ichead uses crl(4,1) still
           else
            head(mi)=(pho(mi)-crl(4,1))/rho1grav + cord_val
     &           -head0
           if(s(mi).le.sat_ich.and.ieos(mi).ne.1) then
            head(mi)= head_id
           endif
          endif            
c     gaz 121306 allow negative heads as necessary
c     head(mi)= max(head(mi),0.0d00)
         enddo

      else if(iflg.eq.3) then
c     
c     convert from head to pressure(all nodes)
c     
         do mi=1,n
            if(mi.gt.neq2) then
               mid=mi-neq2
            elseif(mi.gt.neq) then
               mid=mi-neq
            else
               mid=mi
            endif
            if(iseth.eq.1) then
               cord_val = 0.0
            else
               cord_val = cord(mid,igrav)
            endif
            pho(mi)=rho1grav*(pho(mi)-cord_val+head0)
     &           +pref
         enddo
      else if(iflg.eq.4) then
c     
c     convert from pressure to head(one node)
c     
         mi=ii 
         if(mi.gt.neq2) then
            mid=mi-neq2        
         elseif(mi.gt.neq) then
            mid=mi-neq        
         else
            mid=mi        
         endif
         if(iseth.eq.1) then
            cord_val = 0.0
         else
            cord_val = cord(mid,igrav)
         endif
c gaz 112619 need crl(4,1) for ichead ne 0   
         if(ichead.ne.0) then
         head_value=(pres_value-crl(4,1))/rho1grav +
     &        cord_val-head0
         else
         head_value=(pres_value-pref)/rho1grav +
     &        cord_val-head0  
         endif
c gaz 090119
         if(ichead.ne.0) then
          if(s(mi).le.sat_ich.and.ieos(mi).ne.1) then
           head_value= head_id
          endif
         endif
c         if(irich.eq.0) then
c           head_value= max(head_value,0.0d00)
c         endif 
         
      else if(iflg.eq.5) then
c     
c     convert from head to pressure(one node)
c     this is for input where head resides in pho variable
c     
         mi=ii    
         if(mi.gt.neq2) then
            mid=mi-neq2
         elseif(mi.gt.neq) then
            mid=mi-neq
         else
            mid=mi
         endif
         if(iseth.eq.1) then
            cord_val = 0.0
         else
            cord_val = cord(mid,igrav)
         endif
c gaz 112619  need pref to be pres0 (iconv or ichead via crl(4,1)  
         if(iconv.ne.0.or.ichead.ne.0) then
          pres_value  =rho1grav*(head_value-cord_val+head0)
     &      +crl(4,1)
         else
          pres_value  =rho1grav*(head_value-cord_val+head0)
     &      + pref            
         endif
      else if(iflg.eq.6) then
c     
c     convert from head to pressure(all nodes,pflow)
c     
         if(irich.eq.0) then
            do mi=1,n
               if(mi.gt.neq2) then
                  mid=mi-neq2
               elseif(mi.gt.neq) then
                  mid=mi-neq
               else
                  mid=mi
               endif
               if(iseth.eq.1) then
                  cord_val = 0.0
               else
                  cord_val = cord(mid,igrav) 
               endif
               if(pflow(mi).ne.0.0) then
                  pflow(mi)=rho1grav*(pflow(mi)-cord_val+head0)
     &                 +pref
               endif
            enddo
	 else
            do mi=1,n
               if(mi.gt.neq2) then
                  mid=mi-neq2
               elseif(mi.gt.neq) then
                  mid=mi-neq
               else
                  mid=mi
               endif
               if(iseth.eq.1) then
                  cord_val = 0.0
               else
                  cord_val = cord(mid,igrav)
               endif
               if(pflow(mi).ne.0.0) then
                  pflow(mi)=rho1grav*(pflow(mi)-cord_val+head0)
     &                 +pref
               endif
            enddo
	 endif
      else if(iflg.eq.7) then
c
c adjust head boundary so that it is consistent with height
c
              do mi=1,n
               if(mi.gt.neq2) then
                  mid=mi-neq2
               elseif(mi.gt.neq) then
                  mid=mi-neq
               else
                  mid=mi
               endif
               if(iseth.eq.1) then
                  cord_val = 0.0
               else
                  cord_val = cord(mid,igrav) 
               endif
              enddo
      endif
c     
      return
      end                
