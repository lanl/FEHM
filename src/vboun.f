      subroutine vboun (iz,ndummy)
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
CD1 To set large volumes for effective boundary sinks
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 1-29-09           G. Zyvoloski   N/A     Initial implementation
CD2
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Counter for number of lines of data read
CD5 ivcnd        int         Parameter used for reading
CD5 mid          int         Do loop index parameter
CD5 mi           int         Current node number
CD5 it           int         Current group number
CD5 itp          int         Current conductivity model number
CD5 vc1          real*8      Current thermal conductivity parameter 1
CD5 vc2          real*8      Current thermal conductivity parameter 2
CD5 vc3          real*8      Current thermal conductivity parameter 3
CD5 vc4          real*8      Current thermal conductivity parameter 4
CD5 vc5          real*8      Current thermal conductivity parameter 5
CD5 vc6          real*8      Current thermal conductivity parameter 6
CD5 vc7          real*8      Current thermal conductivity parameter 7
CD5 vc8          real*8      Current thermal conductivity parameter 8
CD5 vc12         real*8      vc1 - vc2
CD5 sqrsat       real*8      Parameter used in calculation
CD5 tmpPor       real*8      Parameter used in calculation
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 Boundary condition options
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
C**********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comdti
      use comai
      use comki
      implicit none
      
      integer iz,ndummy,i,ivcnd,mid,mi,it,itp
  

c read in data
      if(ivboun.ne.0) then
         if(iz.eq.0) then
            if(.not.allocated(ivbounf)) then
               allocate(ivbounf(n0))
               ivbounf = 0
            endif
            narrays = 1
            itype(1) = 4
            default(1) = 1
            macro = "vbou"
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(9), macro, igroup, ireturn,
     2           i4_1=ivbounf(1:n0) )
            
            macroread(9) = .TRUE.
            
         else
c
c change volumes to hold conditions constant
c
            do mid=1,neq
               mi=mid+ndummy
               itp = ivbounf(mi)
               if(itp.gt.0) then
                  sx1(mi)= sx1(mi)*itp
                  ps(mi) = 1.0
               else if(itp.le.0) then 
                  sx1(mi)= sx1(mi)*abs(itp)
               endif
            enddo
           deallocate(ivbounf)
         end if
      end if

      return
      end



