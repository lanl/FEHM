      subroutine sice(iz)
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
CD1 To read input, write output, or set variable values for ice
CD1 solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 ?            G. Zyvoloski   N/A     Initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/sice.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:46   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:46   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:48   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:40 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Wed Jun 12 15:21:22 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.4   Mon Jun 03 11:18:36 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.3   Fri May 31 11:02:12 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.2   Fri Feb 02 10:42:34 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.1   03/18/94 16:12:14   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:50   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 iz           int     I       Control parameter
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 inpt                   I    File contains input data
CD3 iout                   I    Main output file
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 ice, inpt, tmelt, wdd1, sii, n,
CD4 ew1, ew2, ew3, ew4, ew5, ew6, ew7, ew8, ew9, ew10, ew11, ew12,
CD4 ev1, ev2, ev3, ev4, ev5, ev6, ev7, ev8, ev9, ev10, ev11, ev12,
CD4 ices, sio, iout, m, nskw, narrays, pointer, itype, default, igroup,
CD4 ireturn, macroread
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier     Type      Description
CD4 
CD4 sther          N/A       Computes thermdynamic parameters for ice
CD4                             solution
CD4 initdata       N/A      Reads and sets input values defined at
CD4                             each node
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
CD5 siin         real*8      Ice saturation at current node
CD5 i            int         Do loop parameter over all nodes
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
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
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
CPS BEGIN sice
CPS 
CPS IF ice solution is required
CPS 
CPS   IF this call is for reading input
CPS   
CPS     Read first line of data
CPS     
CPS     Set parameters for reading ice saturations
CPS       
CPS     initdata - read ice saturations and set at each node
CPS       
CPS     Set flag to denote that the ice macro has been called
CPS       
CPS     IF this is an ice simulation
CPS       Define thermodynamic parameters
CPS       sther - compute thermodynamic parameters
CPS     ENDIF
CPS   
CPS   ELSEIF this call is for initializing variables
CPS   
CPS     IF this is the start of a simulation
CPS   
CPS       FOR each node
CPS     
CPS         Set ice saturation to input value if it is 0
CPS     
CPS         IF ice saturation is 0
CPS           Set flag to indicate no ice
CPS         ELSE
CPS           Set flag to indicate ice
CPS         ENDIF
CPS     
CPS       ENDFOR
CPS   
CPS     ENDIF
CPS     
CPS   ELSEIF this call is to update ice saturations
CPS   
CPS     FOR each node
CPS       Update ice saturation
CPS     ENDFOR
CPS     
CPS   ELSEIF this call is for output
CPS     
CPS     Write node numbers and ice saturations
CPS     
CPS   ENDIF
CPS   
CPS ENDIF
CPS 
CPS END sice
CPS
CPS
C**********************************************************************

      use combi
      use comci
      use comdi
      use comii
      use comgi
      use comfi
      use comei
      use comdti
      use comai
      use comki
      implicit none

      real*8 siin
      integer i,iz

      if ( ice .ne. 0) then
         if ( iz .eq. 0 )  then

c     read in input when ice is present
            read (inpt  ,   *) ice,siin,tmelt

c     initial ice content
            narrays = 1
            itype(1) = 8
            default(1) = 0.
            macro = "ice "
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     2           r8_1 = sii(1:n0) )
            
            macroread(8) = .TRUE.

c     if ice=0 then return and forget ice calculations
            if ( ice .ne. 0) then

c     define linear thermodynamics
               ew1=0.1
               ew2=20.0
               ew3=997.8
               ew4=0.5674
               ew5=-0.1889
               ew6=8.39e-2
               ew7=9.44e-4
               ew8=4.18e-3
               ew9=1.296e-3
               ew10=0.0
               ew11=0.0
               ev1=0.1
               ev2=1.00
               ev3=0.5e-2
               ev4=0.0
               ev5=0.0
               ev6=0.0
               ev7=0.0
               ev8=0.0
               ev9=0.0
               ev10=0.0
               ev11=0.0

c     change to simple thermdynamics
               call sther(2)
            end if

         elseif ( iz .eq. 1 )  then

c     initialize variables
            if ( days .eq. 0.0 )  then
               do 10 i=1,n
                  if ( sii(i) .eq. 0.0) sii(i)=siin
                  if ( sii(i) .eq. 0.0 )  then
                     ices(i)=0
                  else
                     ices(i)=1
                  end if
 10            continue
            end if

         elseif ( iz .eq. 3 )  then

c     update ice saturation
           do i=1,n
               sio(i)=sii(i)
            enddo

         elseif ( iz .eq. 5 )  then

c     write out ice information
            if (iout .ne. 0) then
               write(iout,*)
     &              '********** ice fraction information *********'
               write(iout,*) 'nodes'
               write(iout,50)(nskw(i),i=1,m)
               write(iout,*) 'ice saturation'
               write(iout,60)(sii(nskw(i)),i=1,m)
            end if
 50         format(1x,8i10)
 60         format(1x,8f10.5)
            
         end if
      endif

      return
      end
