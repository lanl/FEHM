      subroutine loadpcnsk(ja,jb,jc,tmpcnsk,nloc)
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
!D1 This subroutine is used to record the location of the source terms 
!D1 and their strength (pcnsk).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 Initial implementation: 17-APR-1997, Programmer: Chun
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/loadpcnsk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:26   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:56   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:14 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!D4 The location of the source terms are recorded in ptindex and the 
!D4 strengths are recorded in pcnsk. By doing this, we do not need to go
!D4 through from node 1 to node n0 for pcnsk and source term information.
!D4 
!***********************************************************************

      use comai
      use comdti
      use comdi
      use combi
      use compart
      implicit none
      
      integer j,ja,jb,jc,nloc
      real tmpcnsk

      if( ja.lt.0 )then
         ja=-ja
         do j=1,n0
            if(izonef(j).eq.ja)then
               nloc=nloc+1
               pcnsk(nloc)=tmpcnsk
               ptindex(nloc)=j
            endif
         end do
      else
         if(jc.eq.0)then
            ja=1
            jb=n0
            jc=1
         endif
         do j=ja,jb,jc
            nloc=nloc+1
            pcnsk(nloc)=tmpcnsk
            ptindex(nloc)=j
         end do
      endif

      return
      end
