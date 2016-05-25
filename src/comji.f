      module comji
!     comji
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global include file for array variables and pointers (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 14-JAN-94    Z. Dash        22      Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comji.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:42   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:36   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:57:12   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:50 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   04/27/95 18:40:24   llt
!D2 added upwind_l, upwind_v, & nbits
!D2 
!D2    Rev 1.1   03/18/94 16:23:18   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:14   pvcs
!D2 original version in process of being certified
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   None
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block it pointers and associated variables *****
!D4   ipit8           POINTER  it     Pointer to variable array it8
!D4   ipit9           POINTER  it     Pointer to variable array it9
!D4   ipit10          POINTER  it     Pointer to variable array it10
!D4   ipit11          POINTER  it     Pointer to variable array it11
!D4   ipit12          POINTER  it     Pointer to variable array it12
!D4   it8             INT      it     Temporary stoarage variable it8
!D4   it9             INT      it     Temporary stoarage variable it9
!D4   it10            INT      it     Temporary stoarage variable it10
!D4   it11            INT      it     Temporary stoarage variable it11
!D4   it12            INT      it     Temporary stoarage variable it12
!D4   
!D4   ***** COMMON Block t pointers and associated variables *****
!D4   ipt1            POINTER  t      Pointer to variable array t1
!D4   ipt2            POINTER  t      Pointer to variable array t2
!D4   ipt3            POINTER  t      Pointer to variable array t3
!D4   ipt4            POINTER  t      Pointer to variable array t4
!D4   ipt5            POINTER  t      Pointer to variable array t5
!D4   ipt5v           POINTER  t      Pointer to variable array t5v
!D4   ipt6            POINTER  t      Pointer to variable array t6
!D4   ipt7            POINTER  t      Pointer to variable array t7
!D4   ipt8            POINTER  t      Pointer to variable array t8
!D4   ipt9            POINTER  t      Pointer to variable array t9
!D4   ipt10           POINTER  t      Pointer to variable array t10
!D4   t1              REAL*8   t      Temporary stoarage variable t1
!D4   t2              REAL*8   t      Temporary stoarage variable t2
!D4   t3              REAL*8   t      Temporary stoarage variable t3
!D4   t4              REAL*8   t      Temporary stoarage variable t4
!D4   t5              REAL*8   t      Temporary stoarage variable t5
!D4   t5v             REAL*8   t      Temporary stoarage variable t5v
!D4   t6              REAL*8   t      Temporary stoarage variable t6
!D4   t7              REAL*8   t      Temporary stoarage variable t7
!D4   t8              REAL*8   t      Temporary stoarage variable t8
!D4   t9              REAL*8   t      Temporary stoarage variable t9
!D4   t10             REAL*8   t      Temporary stoarage variable t10
!D4   
!D4   ***** COMMON Block perm variables *****
!D4   perml           REAL*8   perm    Temporary stoarage variable perml
!D4   permv           REAL*8   perm    Temporary stoarage variable permv
!D4   
!D4   ***** COMMON Block bn_pass pointers and associated variables *****
!D4   ipbn            POINTER  t      Pointer to variable array bn 
!D4   bn              REAL*8   t      Temporary stoarage variable bn
!D4
!D4   ***** COMMON Block upwind pointers and associated variables ****
!D4   ipupwind_l      POINTER  upwind Pointer to variable array upwind_l 
!D4   ipupwind_v      POINTER  upwind Pointer to variable array upwind_r
!D4   upwind_l        REAL*8   upwind Upwind nodes for liquid phase
!D4   upwind_v        REAL*8   upwind Upwind nodes for vapor phase
!D4   nbits           INTEGER  upwind Number of bits
!D4
!D4 Global Subprograms
!D4
!D4   None
!D4
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 N/A
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS
!PS None
!PS
!***********************************************************************

      real*8 perml(3), permv(3) 
      integer nbits

!     ***** Pointers in COMMON Block bn_pass *****
      real*8, allocatable ::  bn(:)

!     ***** Pointers in COMMON Block it *****
      integer, allocatable :: it4(:) 
      integer, allocatable :: it5(:) 
      integer, allocatable :: it6(:) 
      integer, allocatable :: it8(:) 
      integer, allocatable :: it9(:) 
      integer, allocatable :: it10(:) 
      integer, allocatable :: it11(:)     
      integer, allocatable :: it12(:)     
      integer, allocatable :: it13(:)     
      integer, allocatable :: it8a(:) 
      integer, allocatable :: it9a(:) 
      integer, allocatable :: it10a(:) 
      integer, allocatable :: it11a(:)     
      integer, allocatable :: it12a(:)     
      integer, allocatable :: it13a(:)   
      
      integer, allocatable :: its21(:,:)     
      integer, allocatable :: its22(:,:)      
      integer, allocatable :: its31(:,:)   
      integer, allocatable :: its32(:,:)      
      integer, allocatable :: its41(:,:)      
      integer, allocatable :: its42(:,:)            

!     ***** Pointers in COMMON Block t *****
      real*8, allocatable ::  t1(:) 
      real*8, allocatable ::  t2(:) 
      real*8, allocatable ::  t3(:) 
      real*8, allocatable ::  t4(:) 
      real*8, allocatable ::  t5(:) 
      real*8, allocatable ::  t5v(:) 
      real*8, allocatable ::  t6(:) 
      real*8, allocatable ::  t7(:) 
      real*8, allocatable ::  t7a(:) 
      real*8, allocatable ::  t8(:) 
      real*8, allocatable ::  t9(:) 
      real*8, allocatable ::  t10(:) 
      real*8, allocatable ::  t13(:) 
      real*8, allocatable ::  t14(:)
      real*8, allocatable ::  t15(:)
      real*8, allocatable ::  t16(:)
      real*8, allocatable ::  t17(:)
      
      real*8, allocatable :: ts21(:,:)     
      real*8, allocatable :: ts22(:,:)      
      real*8, allocatable :: ts31(:,:)   
      real*8, allocatable :: ts32(:,:)      
      real*8, allocatable :: ts41(:,:)      
      real*8, allocatable :: ts42(:,:) 
	
      real*8, allocatable :: ts51(:)
      real*8, allocatable :: ts52(:) 
      real*8, allocatable :: ts53(:) 	      

!     ***** Pointers in COMMON Block eq2 *****
      real*8, allocatable ::  upwind_l(:)
      real*8, allocatable ::  upwind_v(:)
      
c parameters for anisotropy

      real*8, allocatable ::   dum_ani(:)
      real*8, allocatable ::   dumx_ani(:)
      real*8, allocatable ::   dumy_ani(:)
      real*8, allocatable ::   dumz_ani(:)
      real*8, allocatable ::   dume_ani(:)
      real*8, allocatable ::   dumex_ani(:)
      real*8, allocatable ::   dumey_ani(:)
      real*8, allocatable ::   dumez_ani(:)
      integer, allocatable ::   it4a(:)
      integer, allocatable ::   it5a(:)
      integer, allocatable ::   it6a(:)
      real*8,  allocatable ::   axy_ani(:)
      real*8,  allocatable ::   axyx_ani(:)
      real*8,  allocatable ::   axyy_ani(:)
      real*8,  allocatable ::   axyz_ani(:)   

c parameters for gravity correction
c gaz 051616
      real*8, allocatable ::  grav_wgt(:)
      

      end module comji

