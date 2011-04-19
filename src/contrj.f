      subroutine contrj (inj)
!***********************************************************************
!  Copyright, 1994, 2004,  The  Regents of the University of California.
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
CD1  To write to the contour output file (for use with PATRAN).
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 09-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/contrj.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:00   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:00   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:56   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Wed Jan 17 11:19:58 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.3   Wed Jan 10 10:58:02 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.2   Tue Jan 09 15:26:26 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.1   03/18/94 15:42:56   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:22:34   pvcs
CD2 original version in process of being certified
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
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************
CD5
CD5 INTERFACES
CD5
CD5 Formal Calling Parameters
CD5
CD5   Identifier      Type     Use  Description
CD5
CD5   inj             INT      I    Flag to indicate what to write
CD5                                    to file
CD5
CD5 Interface Tables
CD5
CD5   None
CD5
CD5 Files
CD5
CD5   Name      Use     Description
CD5
CD5   iscon     O       Contour information file
CD5
C***********************************************************************
CD6
CD6 GLOBAL OBJECTS
CD6
CD6 Global Constants
CD6
CD6   None
CD6
CD6 Global Types
CD6
CD6   None
CD6
CD6 Global Variables
CD6
CD6                            COMMON
CD6   Identifier      Type     Block  Description
CD6
CD6 ipqz, neq, iscon, wdd, n, nei, jdate, jtime, verno, cord, ns, nelm,
CD6 days, npn, pnx, pny, pnz, phi, pcp, t, s
CD6
CD6 Global Subprograms
CD6
CD6   Identifier    Type    Description
CD6 
CD6   veloc         N/A     Computes velocity vectors
CD6   concen        N/A     Controls solute operation (called to write
CD6                            data to output file)
CD6
C***********************************************************************
CD7
CD7 LOCAL IDENTIFIERS
CD7
CD7 Local Constants
CD7
CD7   None
CD7
CD7 Local Types
CD7
CD7   None
CD7
CD7 Local variables
CD7
CD7   Identifier      Type     Description
CD7
CD7   i               int      do loop index
CD7   j               int      do loop index
CD7   np2             int      Index for writing arrays
CD7
CD7 Local Subprograms
CD7
CD7   None
CD7
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN contrj
CPS
CPS IF this call is to write coordinate information
CPS
CPS   Write header information
CPS
CPS   FOR each node
CPS      Write coordinate information
CPS   ENDFOR each node
CPS
CPS   FOR each element
CPS     Write information about the element
CPS
CPS     IF there are not 8 points comprising the element
CPS        Write the node numbers for this element
CPS     ELSE there are 8 points
CPS        Write the node numbers for this element
CPS     ENDIF
CPS
CPS   ENDFOR each element
CPS
CPS ELSE variable data are being written
CPS
CPS   veloc - compute velocities
CPS   Write header information
CPS
CPS   FOR each node
CPS     Write velocities, pressures, temperatures, and saturations
CPS   ENDFOR each node
CPS
CPS   concen - write solute concentration information
CPS
CPS ENDIF
CPS
CPS END contrj
CPS 
C***********************************************************************

      use combi
      use comci
      use comdi
      use comdti
      use comai
      implicit none

      integer i,j,np2,inj
      
c**** 2-dimension , 3-dimension output ****
      if (inj.eq.0) then
         ipqz = neq
         write(iscon ,6000)  'node-element'
 6000    format(a12)
         write(iscon ,6010)  25 , 0 , 0 , 1
 6010    format(i2,8i8)
         write(iscon ,6011)  wdd
 6011    format(a80)
         write(iscon ,6020)  26 , 0 , 0 , 1 , n , nei , 0 , 0 , 0
 6020    format(i2,8i8)
         write(iscon ,6021)  jdate , jtime , verno
 6021    format(a8,4x,a8,a30,2x)
         
         do i=1,n
            write(iscon ,6030)  1 , i , 0 , 2
 6030       format(i2,8i8)
            write(iscon ,6031)  cord(i,1) , cord(i,2) , cord(i,3)
 6031       format(3e16.9)
            write(iscon ,6032)  1 , 'G' , 6 , 0 , 0 , 0,0,0,0,0,0
 6032       format(i1,a1,3i8,2x,6i1)
         enddo
         
         do i=1,nei
            write(iscon ,6040)  2 , i , 8 , 2 , 0 , 0
 6040       format(i2,8i8)
            write(iscon ,6041)  ns , 0 , 5 , 0 , 0.0 , 0.0 , 0.0
 6041       format(4i8,3e16.9)
            if (ns.ne.8) then
               write(iscon ,6042)  ( nelm((i-1)*ns+j) , j=1,ns )
            else
               write(iscon ,6042)    nelm((i-1)*ns+1) ,
     *              nelm((i-1)*ns+4) ,
     *              nelm((i-1)*ns+3) ,
     *              nelm((i-1)*ns+2) ,
     *              nelm((i-1)*ns+5) ,
     *              nelm((i-1)*ns+8) ,
     *              nelm((i-1)*ns+7) ,
     *              nelm((i-1)*ns+6)
            endif
 6042       format(10i8)
         enddo
         
      endif
      
c**** variable output ****
c     
      if (inj.ne.0) then
         call  veloc
         write(iscon ,6100)  'calc-data   ' , days
 6100    format(a12,g10.5)
         write(iscon ,6101)  n , n , 0.0 , 0 , 10
 6101    format(2i9,e15.6,2i9)
         write(iscon ,6102)  wdd
 6102    format(a80)
         write(iscon ,6103)  jdate , jtime , verno
 6103    format(a8,4x,a8,a30,2x)
         
         npn=n
         np2=n+n
         
         do i=1,n
            write(iscon,6110) i,pnx(npn+i),pny(npn+i),pnz(npn+i),
     *           pnx(np2+i),pny(np2+i),pnz(np2+i),
     *           phi(i),phi(i)-pcp(i),t(i),s(i)
 6110       format(i8,(5e13.7))
         enddo

         call concen (5,0)

      endif

      return
      end
