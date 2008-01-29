      module comdti
!     comdti
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
!D1 Global include file for parameters (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 30-SEP-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comdti.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:36 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   06/02/95 10:41:58   llt
!D2 added mdmax (gaz)
!D2 
!D2    Rev 1.1   03/18/94 16:23:04   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:22:00   pvcs
!D2 original version in process of being certified
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3   None
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
!D4   ***** COMMON Block param  *****
!D4   kgmres          INT      param  (based on idof) idof=1, 28*n0+5*n7,
!D4                                     idof=2 or 4, 21*no+5*n7, 
!D4                                     idof=3, 14*n0+5*n7
!D4   ldn             INT      param  Maximum array space needed for jacobian 
!D4                                     array matrix
!D4   ldn1            INT      param  ldn - ldn2
!D4   ldn2            INT      param  ldn / 2 + 1
!D4   lenintg         INT      param  Converts bits to words for allocating 
!D4                                     memory
!D4   lenreal         INT      param  Converts bits to words for allocating 
!D4                                     memory
!D4   maxor           INT      param  Length of gmres working arrays
!D4   n0              INT      param  Maximum number of nodes allowed
!D4   n044            INT      param  n0 * 44
!D4   n2              INT      param  2 * n0, storage parameter
!D4   n3              INT      param  3 * n0, storage parameter
!D4   n4              INT      param  Array storage parameter for 
!D4                                     noncondensible gas solution
!D4   n5              INT      param  Array storage parameter for dual porosity
!D4                                     solution
!D4   n6              INT      param  Array storage for ice solution
!D4   n7              INT      param  Array storage for tracer solution
!D4   n7a             INT      param  Array storage for drc storage     
!D4   n8              INT      param  Array storage for variable  porosity 
!D4                                     solution
!D4   nbd             INT      param  180 * n0 maximum array space for 
!D4                                     incomplete lu decomposition matrix
!D4   ne1             INT      param  Array size of common block /febb/, 
!D4                                     nelmd + nnop + 4n0
!D4   ne2             INT      param  Array size of common /feb/, 
!D4                                     n0 + 9nr + 6nq
!D4   ne3             INT      param  Array size of common block /fbs/, 
!D4                                     3n0 + 336
!D4   ne5             INT      param  Array size of common block /fcc/, 
!D4                                     39n0 + 5n7
!D4   ne6             INT      param  Array size of common block /fdd/, 
!D4                                     36n0
!D4   ne7             INT      param  Array size of common block /fdd1/, 
!D4                                     14n7 + 870
!D4   ne8             INT      param  Array size of common block /ffd2/, 
!D4                                     7n8 + 4
!D4   ne9             INT      param  Array size of common block /fddi/, 
!D4                                     5n0 + 1
!D4   ne10            INT      param  Array size of common block /fhh/, 
!D4                                     14n0
!D4   ne11            INT      param  Array size of common block /co2/, 
!D4                                     32n4 + 1
!D4   ne12            INT      param  Array size of common block /fgg/, 
!D4                                     9nn + 2n3
!D4   ne13            INT      param  Array size of common block /dualp/, 
!D4                                     23n5
!D4   ne14            INT      param  Array size of common block /fice/
!D4                                     2*n6 +1
!D4   ne15            INT      param  Array size of common block /ice/
!D4                                     n6
!D4   nelmd           INT      param  Maximum array space for element 
!D4                                     connectivity array and (later) the nodal
!D4                                     connectivity array
!D4   nelucm          INT      param  nbd / 64
!D4   nn              INT      param  Maximum number of connected elements
!D4   nnop            INT      param  Maximum array space for lu decomposition 
!D4                                     connectivity array
!D4   nq              INT      param  Maximum array space for each finite 
!D4                                     element coefficient array associated
!D4                                     with the stress solution
!D4   nr              INT      param  Maximum space allowed for each finite 
!D4                                     element coefficient array
!D4   mdmax           INT      param  Maximun multiply-connected nodes
!D4   zero_t          REAL*8   param  Value for real parameter zero
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
!D8 Parameter values are input or calculated in the main program (fehm5j) for
!D8 use in dynamic allocation. 
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

      real*8 zero_t
      integer kgmres, ldn, ldn1, ldn2, lenintg, lenreal, maxor
      integer n0, n044, n2, n3, n4, n5, n6, n7, n8, nbd, ne1 
      integer ne2, ne3, ne5, ne6, ne7, ne8, ne9, ne10, ne11 
      integer ne12, ne13, ne14, ne15, nelmd, nelucm, nn, nnop
      integer nq, nr, n7a, mdmax 

      end module comdti
