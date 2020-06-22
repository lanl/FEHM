      subroutine anonp
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
C***********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 Categorize elements and call routines to generate finite element
CD1 coefficients.
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 07-JAN-94    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/anonp.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.13   Fri Apr 26 15:01:48 1996   gaz
CD2 mode mdnode changes as well ae sx_combine calls
CD2 
CD2    Rev 1.12   Thu Mar 21 13:16:52 1996   hend
CD2 Cleaned up format of code
CD2
CD2    Rev 1.11   Wed Feb 07 10:34:12 1996   gaz
CD2 corrections for mdnodes
CD2 
CD2    Rev 1.10   Mon Jan 29 13:10:50 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.9   Mon Jan 29 10:10:14 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   Fri Jan 12 17:44:14 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.7   Tue Jan 09 14:00:48 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.6   11/15/95 16:14:10   gaz
CD2 changes for sx(iw,3) instead of sx(iw,6)
CD2 
CD2    Rev 1.5   08/02/95 15:11:46   gaz
CD2 added coding for multiply defined nodes
CD2 
CD2    Rev 1.4   03/29/95 12:33:34   llt
CD2 checked arrays iplace & idum for references to 0 element
CD2 initialized array istrw
CD2 
CD2    Rev 1.3   03/10/95 10:19:08   llt
CD2 prolog information added - zvd, gaz
CD2 
CD2    Rev 1.2   01/28/95 13:53:56   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.1   03/18/94 15:57:34   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:18   pvcs
CD2 original version in process of being certified
CD2 
c 15-3-94
c got rid of memory calls
c 12/15/94 gaz changed over to nopdum(nei,nemx) and noodum
c 12/22/94 gaz allocate memory for istrw,nelm(change), and sx here
c 12/22/94 gaz defined ldn,nelmd in anonp
c 12/22/94 gaz formed estimate ldn2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   iout                     O    File used for general code output
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   intg            INT      faai   Indicates integration type used
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   istrs           INT      faai   Parameter indicating if the stress
CD4                                     solution is enabled
CD4   istrw           INT      fbb    Starting positions in sx(nr,9) array of
CD4                                     finite element coefficients for each
CD4                                     node
CD4   ldn             INT      param  Maximum array space needed for jacobian
CD4                                     array matrix
CD4   ldn1            INT      param  ldn - ldn2
CD4   ldn2            INT      param  ldn / 2 + 1
CD4   lenintg         INT      param  Converts bits to words for allocating
CD4                                     memory
CD4   lenreal         INT      param  Converts bits to words for allocating
CD4                                     memory
CD4   maxor           INT      param  Length of gmres working arrays
CD4   n0              INT      param  Maximum number of nodes allowed
CD4   n3              INT      param  3 * n0, storage parameter
CD4   n6              INT      param  Array storage for ice solution
CD4   n7              INT      param  Array storage for tracer solution
CD4   n7a             INT      param  Array storage for drc storage
CD4   n8              INT      param  Array storage for variable  porosity
CD4                                     solution
CD4   nbd             INT      param  180 * n0 maximum array space for
CD4                                     incomplete lu decomposition matrix
CD4   nei             INT      faai   Total number of elements in the problem
CD4   nelmd           INT      param  Maximum array space for element
CD4                                     connectivity array and (later) the nodal
CD4                                     connectivity array
CD4   nelucm          INT      param  nbd / 64
CD4   nelm            INT      fbb    Initially information about nodes in each
CD4                                     element, later nodal connectivity
CD4                                     information
CD4   nelmdg          INT      fbb    Contains position of (i,i) element in
CD4                                     connectivity array
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   nop             INT      fbb    Matrix sparsity structure for lu
CD4                                     decomposition
CD4   nnop            INT      param  Maximum array space for lu decomposition
CD4                                     connectivity array
CD4   nq              INT      param  Maximum array space for each finite
CD4                                     element coefficient array associated
CD4                                     with the stress solution
CD4   nr              INT      param  Maximum space allowed for each finite
CD4                                     element coefficient array
CD4   ns              INT      faai   Number of nodes per element
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   determ          N/A      Evaluate determinates.
CD4   gencof          N/A      Generate finite element coefficients (mixed 
CD4                              elements) and perform numerical integration 
CD4                              of the elements.
CD4   mmgetblk        N/A      Allocate space for an array.
CD4   mmrelblk        N/A      Deallocate space for an array.
CD4   storage_derivatives      Allocate memory for derivative arrays and
CD4                   N/A        initialize values, or deallocate   
CD4                              derivative arrays.
CD4   zeroi_out       N/A      Assign value of zero to all elements of an
CD4                              integer array
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   alen1           REAL*8   Intermediate length used in volume calculation
CD5   alen2           REAL*8   Intermediate length used in volume calculation
CD5   alen3           REAL*8   Intermediate length used in volume calculation
CD5   det             REAL*8   Determinate value used in volume calculation
CD5   difvol          REAL*8   Difference in volumes of two elements
CD5   dm              REAL*8   Array of internode coefficients for
CD5                              tetrahedrals
CD5   i               INT      Loop parameter
CD5   i1              INT      Loop parameter
CD5   i2              INT      Loop parameter
CD5   i3              INT      Loop parameter
CD5   i4              INT      Loop parameter
CD5   i5              INT      Loop parameter
CD5   i6              INT      Loop parameter
CD5   ib              INT      Loop parameter
CD5   icode           INT      Error return code.
CD5   iconn           INT      Orthgonality identifier for all elements
CD5   id              INT      Loop parameter
CD5   idiff           INT      Difference in ordering of two elements
CD5   idum            INT      Array of neighbor nodes of a given node
CD5   ie              INT      Element identifer
CD5   ied             INT      Loop parameter
CD5   ig1             INT      Loop parameter
CD5   ij              INT      Loop parameter
CD5   ik              INT      Loop parameter
CD5   in              INT      Loop parameter
CD5   in1             INT      Loop parameter
CD5   in2             INT      Loop parameter
CD5   incon           INT      Counter used in creating connectivity array
CD5   ineluf          INT      Array first elements in unique element list
CD5   inoduf          INT      First node in each unique node set
CD5   int             INT      Specific node in unique node set

CD5   iorth           INT      Orthgonality identifier for each element
CD5   iortho          INT      Specific value of array iorth
CD5   ipiv            INT      Specific value of pivot location array
CD5   is1             INT      Local node in element
CD5   is2             INT      Local node in element
CD5   isx             INT      Counter used in creating istrw array
CD5   iu              INT      Loop parameter
CD5   j               INT      Loop parameter
CD5   je              INT      Loop parameter
CD5   jj              INT      Loop parameter
CD5   k               INT      Local node number in an element
CD5   kb              INT      Neighboring node
CD5   kc              INT      Neighboring node
CD5   kj              INT      Loop parameter
CD5   kk              INT      Loop parameter
CD5   korth           INT      Connectivity of 8 node element accounting for
CD5                              orthogonality
CD5   lei             INT      Loop parameter
CD5   mm              INT      Loop parameter
CD5   n1              INT      Counter used in finding element connectivity
CD5   narr            INT      Array used in comparing nodal ordering within 
CD5                              elements 
CD5   ncon            INT      Connectivity array used when array nelm contains
CD5                              element information
CD5   ncont           INT      Size of ncon array

CD5   ned             INT      Local node in element
CD5   neibr           INT      Neighboring node
CD5   nel             INT      Element identifer
CD5   nele            INT      Element identifer
CD5   nelu            INT      Element identifer
CD5   neluc           INT      Intermediate pointer for unique elements
CD5   nelucd          INT      Specific value of neluc array
CD5   neluf           INT      Array of unique elements
CD5   neqp1           INT      Number of nodes plus 1
CD5   neu             INT      Element identifier
CD5   nf              INT      Number of elements connected to the nodes
CD5   nfd             INT      Specific value of array nf
CD5   nfu             INT      Specific value of array nf
CD5   nj              INT      Loop parameter
CD5   node            INT      Node in an element
CD5   nodeu           INT      A specific unique node
CD5   nodeuf          INT      Unique node group associated with each node
CD5   noo             INT      Local position of node in the elements 
CD5                              connected to it
CD5   nsc             INT
CD5   nsf             INT      Array of element types (nodes per element)
CD5   nsl             INT      Specific value of array nsf
CD5   nslu            INT      Specific value of array nsf
CD5   nt              INT      Specific value of array nodeuf
CD5   numgb           INT      Intermediate result used during sorting of 
CD5                              element nodes
CD5   sumu            REAL*8   Intermediate sum
CD5   tol             REAL*8   Tolerance used in comparing volumes of elements
CD5   tp1             REAL*8   Result of triple product used in volume 
CD5                              calculation
CD5   tp2             REAL*8   Result of triple product used in volume 
CD5                              calculation
CD5   tp3             REAL*8   Result of triple product used in volume 
CD5                              calculation
CD5   tp4             REAL*8   Result of triple product used in volume 
CD5                              calculation
CD5   vlorth          REAL*8   Volume of element assuming orthogonal coordiates
CD5   vole            REAL*8   Volume estimate for elements
CD5   volu            REAL*8   Specific value of array vole
CD5   x1              REAL*8   Argument for triple vector product   
CD5   x2              REAL*8   Argument for triple vector product  
CD5   x3              REAL*8   Argument for triple vector product   
CD5   x12             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x13             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x14             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element  
CD5   x15             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element  
CD5   x23             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x34             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x37             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x45             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x46             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x56             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x58             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x67             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   x78             REAL*8   Difference in x coordinates of local nodes 
CD5                              in an element   
CD5   xie             REAL*8   X coordinate value used in comparing elements
CD5   xnelu           REAL*8   X coordinate value used in comparing elements
CD5   y1              REAL*8   Argument for triple vector product   
CD5   y2              REAL*8   Argument for triple vector product   
CD5   y3              REAL*8   Argument for triple vector product   
CD5   y12             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y13             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y14             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y15             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element  
CD5   y23             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y34             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y37             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y45             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y46             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y56             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y58             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y67             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   y78             REAL*8   Difference in y coordinates of local nodes 
CD5                              in an element   
CD5   yie             REAL*8   Y coordinate value used in comparing elements
CD5   ynelu           REAL*8   Y coordinate value used in comparing elements
CD5   z1              REAL*8   Argument for triple vector product   
CD5   z2              REAL*8   Argument for triple vector product   
CD5   z3              REAL*8    Argument for triple vector product  
CD5   z12             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z13             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z14             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z15             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z23             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z34             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z37             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z45             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z46             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z56             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z58             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z67             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   z78             REAL*8   Difference in z coordinates of local nodes 
CD5                              in an element   
CD5   zero_d          REAL*8   Value of zero assigned to a variable (for
CD5                              mixed precision)
CD5   zie             REAL*8   Z coordinate value used in comparing elements
CD5   znelu           REAL*8   Z coordinate value used in comparing elements
CD5
CD5 Local Subprograms
CD5
CD4   Identifier      Type     Description
CD5   
CD5   al              REAL*8   Calculate the length of a vector.
CD5   tp              REAL*8   Calculate the triple vector product.
CD5   vp              REAL*8   Calculate the vector product.
CD5
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.2 Finite-Element Coefficient Generation
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN anonp
CPS 
CPS   define local functions for vector length, vector product and triple . . .
CPS   . . . vector product
CPS   call storage_derivatives to deallocate derivative storage arrays
CPS   call mmgetblk to allocate required storage for ncon, vole, iorth, . . .
CPS   . . . nf, idum, narr, ineluf, nsf, neluf, nodeuf, inoduf
CPS   call zeroi_out to zero out ncon and iorth
CPS   
CPS   set tolerance
CPS   
CPS   IF stress solution is enabled or Gauss integration is used
CPS      turn on flag to indicate fully connected orthogonal elements
CPS   END IF
CPS   
CPS   set ???
CPS   FOR each node
CPS   
CPS       FOR each node per element
CPS       
CPS           FOR each element
CPS               IF this element includes the current node ???
CPS                  save position of node in element
CPS               END IF
CPS           END FOR
CPS           
CPS       END FOR
CPS       set ???
CPS       
CPS   END FOR
CPS   
CPS   FOR each element
CPS   
CPS       FOR each node in element
CPS           IF ??? node is unique or exists in element
CPS              count node
CPS           END IF
CPS           set ???
CPS       END FOR
CPS       set ???
CPS       
CPS       FOR each node in element
CPS           set ???
CPS           FOR each ???
CPS           EXIT IF ???    
CPS               set ???
CPS           END FOR
CPS           set ???
CPS       END FOR
CPS       
CPS   END FOR
CPS   
CPS   FOR each element (1)
CPS       set ???
CPS       IF this is a 3-D problem
CPS       
CPS          IF ??? 8 (3-D 8 node brick elements)
CPS             FOR each ???
CPS                 set ???
CPS             END FOR
CPS             calculate cords and volume (triple product)
CPS             IF orthogonal elements aren't fully connected
CPS                calculate othogonal volume
CPS                IF volumes are within tolerance
CPS                   set ??? to 1 (orthogonal elements)
CPS                ELSE
CPS                   set ??? to 0 (nonorthogonal elements)
CPS                END IF
CPS             END IF
CPS          END IF
CPS          
CPS          IF ??? 6 (3-D 6 node prism elements)
CPS             FOR each ???
CPS                 set ???
CPS             END FOR
CPS             calculate cords and volume (triple product)
CPS             set ??? to 0 (prisms are always nonorthoganal)
CPS          END IF
CPS          
CPS          IF ??? 4 (3-D 4 node tetrahedral elements)
CPS             FOR each ???
CPS                 set ???
CPS             END FOR
CPS             call determ to calculate volumes
CPS             set ??? to 0 (tetrahedrals are always nonorthoganal)
CPS          END IF
CPS          
CPS       END IF
CPS       
CPS       IF this is a 2-D problem
CPS       
CPS          IF ??? 4 (2-D 4 node quadrilaterals)
CPS             FOR each ???
CPS                 set ???
CPS             END FOR
CPS             calculate cords and volume (triple product)
CPS             IF orthogonal elements aren't fully connected
CPS                calculate othogonal volume
CPS                IF volumes are within tolerance
CPS                   set ??? to 1 (orthogonal elements)
CPS                ELSE
CPS                   set ??? to 0 (nonorthogonal elements)
CPS                END IF
CPS             END IF
CPS          END IF
CPS          
CPS          IF ??? 3 (2-D 3 node triangles)
CPS             FOR each ???
CPS                 set ???
CPS             END FOR
CPS             calculate cords and volume (triple product)
CPS             set ??? to 0 (triangles are always nonorthoganal)
CPS          END IF
CPS          
CPS       END IF
CPS       
CPS       IF ???
CPS          set ???
CPS          FOR each ??? (2)
CPS              set
CPS              IF ??? within tolerance
CPS                 FOR each ???
CPS                     calcualte ???
CPS                     IF ??? not zero
CPS                        goto end of loop (2)
CPS                     END IF
CPS                 END FOR
CPS                 
CPS                 set ???
CPS                 FOR each ???
CPS                     calculate ???
CPS                 END FOR
CPS                 IF ??? within tolerance
CPS                    set ??? and got to end of loop (1)
CPS                 END IF
CPS              END IF
CPS          END FOR (2)
CPS       
CPS       END IF
CPS       
CPS       set ??? 
CPS       
CPS   END FOR (1)
CPS   
CPS   FOR each node
CPS       set ???
CPS       FOR each ???
CPS           set ???
CPS       END FOR
CPS       set ???
CPS   END FOR
CPS   
CPS   FOR each node (3)
CPS       set ???
CPS       IF ???
CPS          FOR each ??? (4)
CPS              set ???
CPS              IF ???
CPS                 set ???
CPS                 FOR each ???
CPS                     calculate ???
CPS                     IF ???
CPS                        exit to end of (4)
CPS                     END IF
CPS                 END FOR
CPS                 set ???
CPS                 FOR each ???
CPS                     calculate ???
CPS                     IF ???
CPS                        exit to end of (4)
CPS                     END IF
CPS                 END FOR
CPS                 set ???
CPS                 exit to end of (3)
CPS              END IF
CPS          END FOR (4)
CPS       END IF
CPS       set ???
CPS   END FOR (3)
CPS   
CPS   set ????
CPS   FOR each node
CPS   
CPS       call zeroi_out to zero out connectivity
CPS       set ???
CPS       
CPS       FOR each ???
CPS       
CPS           set ???
CPS           IF this is a 3-D problem
CPS              FOR each node in element
CPS                  set ???
CPS                  IF this isn't a brick element or this element is . . .
CPS                  . . . nonorthogonal
CPS                     set ???
CPS                  ELSE
CPS                     calculate ???
CPS                  END IF
CPS              END FOR
CPS           ELSE (this is 2-D)
CPS              FOR each node in element
CPS                  set ???
CPS                  IF this isn't a quadrilateral or this element is . . .
CPS                  . . . nonorthogonal
CPS                     set ???
CPS                  ELSE
CPS                     calculate ???
CPS                  END IF
CPS              END FOR
CPS           END IF
CPS           
CPS       END FOR
CPS       
CPS       FOR each node
CPS           IF ??? not zero
CPS              set ???
CPS           END IF
CPS           IF ???
CPS              ???
CPS           END IF
CPS       END FOR
CPS       
CPS       set ???
CPS       
CPS   END FOR
CPS   
CPS   set ???
CPS   FOR each node
CPS       set ???
CPS       IF ???
CPS          set ???
CPS          FOR each ???
CPS              set ???
CPS              FOR each ???
CPS                  set ???
CPS                  IF ???
CPS                     goto (tag) 
CPS                  END IF
CPS              END FOR
CPS              error message but no exit ???
CPS              (tag)
CPS              set ???
CPS          END FOR
CPS          
CPS          FOR each
CPS              set ???
CPS          END FOR
CPS          
CPS       ELSE
CPS       
CPS          set ???
CPS          FOR each ???
CPS              set ???
CPS          END FOR
CPS       
CPS       END IF
CPS       
CPS   END FOR
CPS   
CPS   set ???
CPS   
CPS   call gencof to generate finite element coefficients and perform . . .
CPS   . . . numerical integration of the elements.
CPS   
CPS   call mmrelblk to deallocate storage for ncon, vole, iorth, . . .
CPS   . . . nf, idum, narr, ineluf, nsf, neluf, nodeuf, inoduf
CPS   call storage_derivatives to reallocate derivative storage arrays
CPS 
CPS END anonp
CPS
C***********************************************************************

      use combi
      use comdi
      use comei
      use comsi
      use comdti
      use comai
      implicit none

      real*8 al, alen1, alen2, alen3, det, difvol, dm(4,4), sumu
      real*8 tol, tp, tp1, tp2, tp3, tp4, vlorth, volest, volu
      real*8 vp, xie, xnelu, yie, ynelu, zero_d, zie, znelu 
      real*8 x1, x2, x3, x12, x13, x14, x15, x23, x34, x37, x45, x46 
      real*8 x56, x58, x67, x78, y1, y2, y3, y12, y13, y14, y15, y23 
      real*8 y34, y37, y45, y46, y56, y58, y67, y78, z1, z2, z3, z12
      real*8 z13, z14, z15, z23, z34, z37, z45, z46, z56, z58, z67, z78 
c gaz 051017      
      real*8 vol_tol, vol_tol_chk, vol_chk_default
      integer i, i1, i2, i3, i4, i5, i6, ib, icode, iconn, id
      integer idiff, ie, ied, ig1, ij, ik, in, in1, in2, incon, int
      integer iortho, ipiv, is1, is2, isx, iu, j, je, jj, k, kb, kc, kj
      integer kk, korth(8,8), lei, mm, n1, ncont, ned, neibr, nel, nele
      integer nelu, neluc, nelucd, neqp1, neu, nfd, nfu, node, nodeu
      integer nj, nsc, nsl, nslu, nt, ncoef, numgb(8)
      integer ii, ipivkb, neq_total, n0_save, ipmax, lcnt, lcnt2
      integer ib_min, ib_max
      parameter (vol_chk_default = 1.e-4)

c using storage for a matrix
      real*8, allocatable ::  vole(:)
      integer, allocatable :: ncon(:)
      integer, allocatable :: iorth(:)
      integer, allocatable :: neluf(:)
      integer, allocatable :: nf(:)
      integer, allocatable :: idum(:)
      integer, allocatable :: nodeuf(:)
      integer, allocatable :: inoduf(:)
      integer, allocatable :: narr(:,:)
      integer, allocatable :: nsf(:)
      integer, allocatable :: ineluf(:)
      integer, allocatable :: nopdum(:,:)
      integer, allocatable :: noodum(:,:)
      integer, allocatable :: iplace(:)        

      data korth / 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
     *             0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1,
     *             1, 0, 0, 0, 1, 1, 0 ,1, 0, 1, 0, 0, 1, 1, 1, 0, 
     *             0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1/
      data neluc / 0 /

c     triple vector product
      tp(x1, y1, z1 ,x2, y2, z2, x3, y3, z3) = 
     .     x1 * y2 * z3 + x2 * y3 * z1 + x3 * y1 *z2 - 
     .     x3 * y2 * z1 - x2 * y1 * z3 - x1 * y3 *z2
c     length of vector 
      al(x1, y1, z1) = sqrt(x1 * x1 + y1 * y1 + z1 * z1)
c     vector product
      vp(x1, y1, x2, y2) = x1 * y2 - y1 * x2
      
      zero_d = 0.0
c
c GAZ 02-08-01
c set neq to neq_primary here for gdpm, then set back later
c also set n0 to neq_primary here for gdpm, then set back later
c
      if(gdpm_flag.ne.0 .or. nriver.ne.0) then
         neq_total = neq
         neq = neq_primary
         n0_save= n0
         n0 = neq_primary
      endif
      
      vol_tol_chk = vol_chk_default
      
c     determine if orthogonal elements should be fully connected
      iconn = 0
      if(istrs .ne. 0 .or. intg .gt. 0) iconn = 1
      
      allocate(iplace(n0))

      do i=1,neq
         iplace(i)=0
      enddo
      
c
c calculate # connections for each node
c
      do ie=1,nei
         do j=1,ns
            i=nelm((ie-1)*ns+j)
            if (i .ne. 0) then
               iplace(i)=iplace(i)+1
            endif
         enddo
      enddo

      ipmax=0
c estimate maximum connections
      if(icnl.eq.0) then
       lcnt=8
       lcnt2=27
       ldn2=lcnt2*neq
      else
       lcnt=4
       lcnt2=9
       ldn2=lcnt2*neq
      endif
c add index space for ncon
c note "-2" is based on ns*connected elements vs difference stencils
      ldn2=ldn2+neq+1
      do i=1,neq
         ipmax=max(ipmax,iplace(i))
         if(iplace(i).gt.lcnt) then
          ldn2=ldn2+iplace(i)*(ns-2)+1-lcnt2
         endif
      enddo
      
      allocate(ncon(ldn2))
      ncon=0
      

      nemx=ipmax
      allocate(vole(nei),iorth(nei),nf(n0),idum(n0))
      allocate(narr(nei,8),ineluf(nei),nsf(nei),neluf(nei))
      allocate(nodeuf(n0),inoduf(n0),nopdum(n0,nemx),noodum(n0,nemx))
      iorth=0

      iplace=0
      do ie=1,nei
         do j=1,ns
            i=nelm((ie-1)*ns+j)
            if (i .ne. 0) then
               iplace(i)=iplace(i)+1
               nopdum(i,iplace(i))=ie
               noodum(i,iplace(i))=j 
            endif
         enddo
      enddo
      
c     find number of nodes in each element
      do ij = 1, nei
         nsc = 0
         do ied = 1, ns
            ik = nelm((ij - 1) *ns + ied)
            if(ik .ne. 0) nsc = nsc + 1
            numgb(ied) = ik
            narr(ij, ied) = ied
         end do
         nsf(ij) = nsc
         
c     sort elements into ascending order
         do j = 2, ns
            n1 = numgb(j)
            do i = j - 1, 1, -1
               if(numgb(i) .lt. n1) go to 17
               narr(ij, i + 1) = narr(ij, i)
               numgb(i + 1) = numgb(i)
            end do
            i = 0
 17         narr(ij, i+1) = j
            numgb(i+1) = n1
         end do
      end do

c group elements(volume, area, moment of inertia)
      do ie = 1, nei
         nsl = nsf(ie)
         if(icnl.eq.0) then
c     procedure for 3-d problems
            if(nsl.eq.8) then
               do 99 lei = 1, 8
                  numgb(lei) = nelm((ie-1)*ns+lei)
 99            continue
c     3-d 8 node brick
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)
               z12 = cord(numgb(1), 3)-cord(numgb(2), 3)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)
               z14 = cord(numgb(1), 3)-cord(numgb(4), 3)

               x23 = cord(numgb(2), 1)-cord(numgb(3), 1)
               y23 = cord(numgb(2), 2)-cord(numgb(3), 2)
               z23 = cord(numgb(2), 3)-cord(numgb(3), 3)

               x34 = cord(numgb(3), 1)-cord(numgb(4), 1)
               y34 = cord(numgb(3), 2)-cord(numgb(4), 2)
               z34 = cord(numgb(3), 3)-cord(numgb(4), 3)

               x56 = cord(numgb(5), 1)-cord(numgb(6), 1)
               y56 = cord(numgb(5), 2)-cord(numgb(6), 2)
               z56 = cord(numgb(5), 3)-cord(numgb(6), 3)

               x58 = cord(numgb(5), 1)-cord(numgb(8), 1)
               y58 = cord(numgb(5), 2)-cord(numgb(8), 2)
               z58 = cord(numgb(5), 3)-cord(numgb(8), 3)

               x67 = cord(numgb(6), 1)-cord(numgb(7), 1)
               y67 = cord(numgb(6), 2)-cord(numgb(7), 2)
               z67 = cord(numgb(6), 3)-cord(numgb(7), 3)

               x78 = cord(numgb(7), 1)-cord(numgb(8), 1)
               y78 = cord(numgb(7), 2)-cord(numgb(8), 2)
               z78 = cord(numgb(7), 3)-cord(numgb(8), 3)

               x15 = cord(numgb(1), 1)-cord(numgb(5), 1)
               y15 = cord(numgb(1), 2)-cord(numgb(5), 2)
               z15 = cord(numgb(1), 3)-cord(numgb(5), 3)

               x37 = cord(numgb(3), 1)-cord(numgb(7), 1)
               y37 = cord(numgb(3), 2)-cord(numgb(7), 2)
               z37 = cord(numgb(3), 3)-cord(numgb(7), 3)
c     
c     calculate triple products(volumes)
c     
               tp1=tp(-x15,-y15,-z15,-x12,-y12,-z12,-x14,-y14,-z14)
               tp2=tp(-x37,-y37,-z37,-x34,-y34,-z34,x23,y23,z23)
               tp3=tp(-x15,-y15,-z15,-x56,-y56,-z56,-x58,-y58,-z58)
               tp4=tp(-x37,-y37,-z37,-x78,-y78,-z78,x67,y67,z67)
               volest=(tp1+tp2+tp3+tp4)/4.0
               vole(ie)=volest
c     
c     check orthogonality if nonstress solution
c     
               iorth(ie) = 0
               if(iconn.eq.0) then
c     calculate orthogonal volume
                  alen1 = al(x14, y14, z14)
                  alen2 = al(x12, y12, z12)
                  alen3 = al(x15, y15, z15)
                  vlorth = alen1*alen2*alen3
c     compare with calculated volume
                  vol_tol = (abs(vlorth)-abs(volest))/abs(volest)
                  if(vol_tol.le.vol_tol_chk) then
                     iorth(ie) = 1
                  else
                     iorth(ie) = 0
                  endif
               endif
            endif
c     
            if(nsl.eq.6) then
               do lei = 1, 6
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     procedure for 6 node prisms
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)
               z12 = cord(numgb(1), 3)-cord(numgb(2), 3)

               x13 = cord(numgb(1), 1)-cord(numgb(3), 1)
               y13 = cord(numgb(1), 2)-cord(numgb(3), 2)
               z13 = cord(numgb(1), 3)-cord(numgb(3), 3)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)
               z14 = cord(numgb(1), 3)-cord(numgb(4), 3)

               x45 = cord(numgb(4), 1)-cord(numgb(5), 1)
               y45 = cord(numgb(4), 2)-cord(numgb(5), 2)
               z45 = cord(numgb(4), 3)-cord(numgb(5), 3)

               x46 = cord(numgb(4), 1)-cord(numgb(6), 1)
               y46 = cord(numgb(4), 2)-cord(numgb(6), 2)
               z46 = cord(numgb(4), 3)-cord(numgb(6), 3)

c     calculate triple products(volumes)
               tp1=tp(-x14,-y14,-z14,-x12,-y12,-z12,-x13,-y13,-z13)
               tp2=tp(-x14,-y14,-z14,-x45,-y45,-z45,-x46,-y46,-z46)
               volest = abs((tp1+tp2)/4.0)
               vole(ie) = volest
c     prisms always are non orthogonal
               iorth(ie) = 0
            endif
            if(nsl.eq.3) then
               do  lei = 1, 3
                  ned = nelm((ie-1)*ns+lei)
                  dm(lei, 1) = cord(ned, 1)
                  dm(lei, 2) = cord(ned, 2)
                  dm(lei, 3) = cord(ned, 3)
               enddo
c     procedure for 3 node triangles in 3-d
c     calculate volume
               call area2d_tri(1,dm(1,1),dm(1,2),dm(1,3),
     &              dm(2,1),dm(2,2),dm(2,3),dm(3,1),dm(3,2),dm(3,3),    
     &              volest,x12,x13,x14,y12,y13,y14)                
               vole(ie) = volest
c     triangles always are non orthogonal
               iorth(ie) = 0
            endif
            if(nsl.eq.4) then
               do lei = 1, 4
                  ned = nelm((ie-1)*ns+lei)
                  dm(lei, 1) = 1.0
                  dm(lei, 2) = cord(ned, 1)
                  dm(lei, 3) = cord(ned, 2)
                  dm(lei, 4) = cord(ned, 3)
               enddo
c     procedure for 4 node tetrahedrals
c     calculate volume
               call determ(det, dm, 4)
               volest = 0.166667*det
               vole(ie) = volest
c     tetrahedrals always are non orthogonal
               iorth(ie) = 0
            endif
         endif
c     
c     procedure for 2-d elements
         if(icnl.ne.0) then
c     4-node quadrilaterals
            if(nsl.eq.4) then
               do lei = 1, 4
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)

               x14 = cord(numgb(1), 1)-cord(numgb(4), 1)
               y14 = cord(numgb(1), 2)-cord(numgb(4), 2)

               x23 = cord(numgb(2), 1)-cord(numgb(3), 1)
               y23 = cord(numgb(2), 2)-cord(numgb(3), 2)

               x34 = cord(numgb(3), 1)-cord(numgb(4), 1)
               y34 = cord(numgb(3), 2)-cord(numgb(4), 2)

               tp1 = vp(-x12, -y12, -x14, -y14)
               tp2 = vp(-x34, -y34, x23, y23)
               volest = (tp1+tp2)/2.0
               vole(ie) = volest
               if(iconn.eq.0) then
c     calculate orthogonal volume
                  alen1 = al(x14, y14, zero_d)
                  alen2 = al(x12, y12, zero_d)
                  vlorth = alen1*alen2
c     compare with calculated volume
                  if(abs(vlorth-volest)/volest.le.tol) then
                     iorth(ie) = 1
                  else
                     iorth(ie) = 0
                  endif
               endif
            endif
            if(nsl.eq.3) then
c     procedure for triangles
               do lei = 1, 3
                  numgb(lei) = nelm((ie-1)*ns+lei)
               enddo
c     calculate volume
               x12 = cord(numgb(1), 1)-cord(numgb(2), 1)
               y12 = cord(numgb(1), 2)-cord(numgb(2), 2)

               x13 = cord(numgb(1), 1)-cord(numgb(3), 1)
               y13 = cord(numgb(1), 2)-cord(numgb(3), 2)

               volest = 0.5*(x12*y13-x13*y12)
               vole(ie) = volest
c     
c     triangles are always nonorthogonal
               iorth(ie) = 0
            endif
         endif
c     
c     group elements
c     
c     **************************************************************
c     assume now 4/22/96 that all elements are unique
c     this is because the following code is
c     rotation invariant and we need to identify x,y,z
c     To go back to old delete first c on every line to 7001
c     ************************************************************** 
         go to 7001
cc     ************************************************************** 
c         if(ie.eq.1) then
c            go to 7001
c         endif
c         nelucd = neluc
c         do 1000 iu = 1, nelucd
c            nelu = ineluf(iu)
c            nslu = nsf(nelu)
c            volu = vole(nelu)
c            difvol = abs(volu-volest)+iabs(nsl-nslu)
c            if(difvol/volest.le.tol) then
cc     
cc     check if ordering is consistent with group
c               do 67 mm = 1, ns
c                  idiff = narr(nelu, mm)-narr(ie, mm)
c                  if(idiff.ne.0) go to 1000
c 67            continue
cc     
cc     
cc     check further for match
c               is1 = nelm((nelu-1)*ns+1)
c               is2 = nelm((ie-1)*ns+1)
c               xnelu = cord(is1, 1)
c               ynelu = cord(is1, 2)
c               znelu = cord(is1, 3)
c               xie = cord(is2, 1)
c               yie = cord(is2, 2)
c               zie = cord(is2, 3)
c               sumu = 0.0
c               do 1010 nj = 2, nsl
c                  kj = nelm((nelu-1)*ns+nj)
c                  kk = nelm((ie-1)*ns+nj)
c                  sumu = abs((cord(kj, 1)-xnelu)-(cord(kk, 1)-xie))
c     *                 +abs((cord(kj, 2)-ynelu)-(cord(kk, 2)-yie))
c     *                 +abs((cord(kj, 3)-znelu)-(cord(kk, 3)-zie))+sumu
c                  if(ithic.ne.0) sumu =  abs(thic(kj)-thic(kk))+sumu
c 1010          continue
cc     if the elements still check we have a match
c               if(abs(sumu**3/volest).le.tol) then
cc     old element type
c                  neluf(ie) = iu
c                  go to 100
c               endif
c            endif
c 1000    continue
cc     new element type (ie. we had no match )
 7001    neluc = neluc+1
         neluf(ie) = neluc
         ineluf(neluc) = ie
      end do
c     
c     find sum of element types connected to node
c     
      do i=1,neq
         nfd = 0
         do je = 1,iplace(i)
            nele = nopdum(i,je)
            nfd = nfd+neluf(nele)
         enddo
         nf(i) = nfd
      enddo
c     
c     group nodes
c     
c     find unique node type based on connected elements
c     
      nodeu = 0
      do 1001 i = 1, neq
c     
c     compare to previous nodes
c     
c     **************************************************************
c     assume now 4/22/96 that all elements are unique
c     this is because the following code is
c     rotation invariant and we need to identify x,y,z
c     To go back to old way uncomment first c on every line to 7002
c     ************************************************************** 
         go to 7002
c     ************************************************************** 
c         if(i.eq.1) then
c            go to 7002
c         endif
cc     
c         do 1003 iu = 1, nodeu
c            in = inoduf(iu)
c            nfu = nf(in)
c            if(nfu.eq.nf(i)) then
cc     check for further match
c               i3 = 1          
c               i4 = iplace(in)
c               do 1004 j = i3, i4
c                  neu = max0(nopdum(i,j), 1)
c                  nel = max0(nopdum(in,j), 1)
c                  if(neluf(neu).ne.neluf(nel)) go to 1003
c 1004          continue
c               in1 = 1
c               in2 = iplace(in)
c               do 1005 j = in1, in2
c                  neu = max0(noodum(i,j), 1)
c                  nel = max0(noodum(in,j), 1)
c                  if(neu.ne.nel) go to 1003
c 1005          continue
cc     nodes in and i match so add to set
c               nodeuf(i) = iu
c               go to 1001
c            endif
c 1003    continue
cc     we have found another unique node set
 7002    nodeu = nodeu+1
         nodeuf(i) = nodeu
         inoduf(nodeu) = i
 1001 continue
c     
c     elements and nodes have been grouped
c     orthogonal elements have been noted
c     
c     the connectivity matrix and the coefficient counter
c     can now be constructed
c     
c     
c    
      if (iout .ne. 0) write(iout, 1010)
      if (iptty .ne. 0) write(iptty, 1010)
 1010 format (' >>>>>> Starting connectivity calcs <<<<<< ') 
      ncon(1) = neq+1
      incon = neq+1
c     zero out connectivity
      do  kk = 1, neq
         idum(kk) = 0
      enddo
      do 2000 i = 1, neq
       ib_min = neq
       ib_max = 1
c     connect node to itself
         idum(i)  = i
         nt = nodeuf(i)
         int = inoduf(nt)
c     
c     load ncon
c     
c     loop on elements connected with node i
c     
         i5 = 1
         i6 = iplace(i)
         do jj = i5, i6
            nele = nopdum(i,jj)
c     find position of node i in element nele
            k = noodum(i,jj)
c     find if element is orthogonal
            iortho = iorth(nele)
            nsl = nsf(nele)
c     coding for 3-d elements
            if(icnl.eq.0) then
               do in = 1, ns
                  node = nelm((nele-1)*ns+in)
                  if (node .ne. 0) then
                   ib_min = min(ib_min,node)
                   ib_max = max(ib_max,node)
                     if (nsl.ne.8.or.iortho.eq.0) then
                        idum(node) = node
                     else
                        idum(node) = korth(k, in)*node
                     endif
                  endif
               enddo
c     coding for 2-d elements
            else
               do in = 1, ns
                  node = nelm((nele-1)*ns+in)
                  if (node .ne. 0) then
                   ib_min = min(ib_min,node)
                   ib_max = max(ib_max,node)
                     if (nsl.ne.4.or.iortho.eq.0) then
                        idum(node) = node
                     else
                        idum(node) = korth(k, in)*node
                     endif
                  endif
               enddo
            endif
      end do
c gaz 11-09-2001 comment out mdnode , call after     `
c     
c     call md_nodes to add multiply defined nodes
c     
c        call md_nodes(1,idum,i)                              

         do ib = ib_min,ib_max
            if(idum(ib).ne.0) then
               incon = incon+1
               ncon(incon) = idum(ib)
               idum(ib) = 0
            endif
            if(ib.eq.i) nelmdg(i) = incon
         enddo

c     save position of the end of node i neighbors
c     
         ncon(i+1) = incon
c     
 2000 continue
c     
c     load istrw array
c     
      neqp1 = neq+1
c     call mmgetblk ("istrw", "combi", ipistrw, 
c     &   incon , 1, icode)
      allocate(istrw(incon))
c gaz 11-9-2001 istrw_itfc and istrw_cold now allocated
c in startup (after mdno sorted out)
c     if(idpdp.eq.0) then
c        allocate(istrw_itfc(incon),istrw_cold(incon))
c     else
c        allocate(istrw_itfc(2*incon),istrw_cold(2*incon))
c     end if
c     istrw_itfc = 0
c     istrw_cold = 0

      do i=1,incon
         istrw(i) = 0
      end do
      if(imdnode.eq.0) then
c     code for no multiple defined nodes
         isx = 0
         do id = 1, neq
            nt = nodeuf(id)
            int = inoduf(nt)
            if(id.eq.int) then
                i1 = ncon(id)+1
               i2 = ncon(id+1)
               ipiv = nelmdg(id)
c     previously defined nodes
               do ie = i1, ipiv-1
                  kb = ncon(ie)
c     search for connection
                  i3 = nelmdg(kb)
                  i4 = ncon(kb+1)
                  do ig1 = i3+1, i4
                     kc = ncon(ig1)
                     if(kc.eq.id) go to 2005
                  enddo
                  write(ierr, *) 'failed in anonp:stop'
                  if(iout .ne. 0) write(iout, *) 'failed in anonp:stop'
                  if(iptty.gt.0) write(iptty, *) 'failed in anonp:stop'
 2005             continue
c     found match
                  istrw(ie-neqp1) = istrw(ig1-neqp1)
               enddo
               do ie = ipiv+1, i2
                  kb = ncon(ie)  
                  isx = isx+1
                  istrw(ie-neqp1) = isx
               enddo
            else
c     other member of group
               i1 = ncon(int)+1
               i2 = ncon(int+1)
               i3 = ncon(id)+1
               do ie = i1, i2
                  istrw(i3+ie-i1-neqp1) = istrw(ie-neqp1)
               enddo
            endif
         enddo
      else
c     code for multiple defined nodes
         do id=1,nei
            ineluf(id) = id
            neluf(id) = id   
         enddo
         neluc=nei
         do id=1,neq
            inoduf(id)=id
            nodeuf(id) = id
         enddo
         nodeu=neq
         isx=0
         do id=1,neq
            i1=ncon(id)+1
            i2=ncon(id+1)
            ipiv=nelmdg(id)
            do ii=ipiv+1,i2
               isx=isx+1
               istrw(ii-neqp1)=isx
            enddo
            do ii=i1,ipiv-1
               kb=ncon(ii)
               i3=ncon(kb+1)
               ipivkb=nelmdg(kb)
               do jj=ipivkb+1,i3
                  if(ncon(jj).eq.id) then
                     istrw(ii-neqp1)=istrw(jj-neqp1)
                  endif
               enddo
            enddo
         enddo
      endif
c     
c 
      if (iout .ne. 0) write(iout, 1009)
      if (iptty .ne. 0) write(iptty, 1009)
 1009 format (' >>>>>> Finished connectivity calcs <<<<<< ') 
      if (iout .ne. 0) write(iout, 2009)
      if (iptty .ne. 0) write(iptty, 2009)
 2009 format (' >>>>>> Starting FE coef. calcs <<<<<< ')   
c     printout positions required for geometric coefficients
c     
      neibr = ncon(neq+1)-ncon(neq)
      iad = isx
      iw = isx
      nr=isx
c     
c     determine if connections need to be broken
c     add one more position for possible BC to BC or mdnodes
      nr=isx+1
c     
c     allocate memory for fe coefficients
c    
      if(istrs.ne.0) then	

c     
c     define pointer for stress coefficients 
c     
         i3 = ncon(neqp1)-neqp1
         if(.not.allocated(istrws)) allocate(istrws(i3))
         if(allocated(sxs)) deallocate(sxs)
         if(icnl.eq.0) then
            ncoef = 15
         else
            ncoef = 8
         endif
	 allocate(sxs(i3,ncoef))
c     GAZ 011109 new storage for finv flow simulations	 
	 allocate(sx(nr,3))
	 sxs = 0.0d0
	 sx = 0.0d0 
	 sx1 = 0.0d0
c     
         i3 = 0
         do i = 1,neq
            i1 = ncon(i)+1
            i2 = ncon(i+1)
            do jj = i1,i2
               kb = ncon(jj)
               i3 = i3+1
               istrws(jj-neqp1) = i3
            enddo
	 enddo
      else 
         allocate(sx(nr,3))
	 sx = 0.0d0 
	 sx1 = 0.0d0
      endif
c     
c     call geometric coefficients routine
c     
      nelucm=neluc
      call gencof(nodeu, neluc, nsf, ineluf, inoduf, neluf, ncon
     *     , nodeuf, iplace, noodum, nopdum, n0, isx)
c     
c     fill nelm array with nodal connections
c     
      ncont = ncon(neq+1)
      nelmd = ncont 
c     
c     adjust space for connectivity matrix
c     
      deallocate(idum,narr,ineluf,nsf)
      deallocate(neluf,nodeuf,inoduf,nopdum)
      deallocate(noodum,iplace)
c     
      deallocate(nelm)
      allocate(nelm(ncont))
c     
      do ik = 1, ncont
         nelm(ik) = ncon(ik)
      enddo
c     
c     call md_nodes to break connections for those nodes
c     
      call md_nodes(2,0,0)                            
c     
      deallocate(ncon,vole,iorth,nf)
c     
c     insure isotropy for non stress solutions
c     
      if(isoy.eq.2) then
         call sx_combine(-1)
c     done now in startup.f
c     call sx_combine(1)
         call sx_combine(2)
      else if(isoy.eq.1) then
         call sx_combine(-2)
c     done now in startup.f
c     call sx_combine(1)
         call sx_combine(3)
      endif
c     
c     GAZ 02-08-01
c     setting neq back 
c     
      if(gdpm_flag.ne.0 .or. nriver.ne.0) then
         neq = neq_total
         n0=n0_save
      endif
      if (iout .ne. 0) write(iout, 2010)
      if (iptty .ne. 0) write(iptty, 2010)
 2010 format (' >>>>>> Finished FE coef. calcs <<<<<< ')
      return
      end

