      subroutine setzone(izone, nin, nsl, xz, yz, zz, irad)
c            subroutine setzone(izone, nin, ncord, nsl, xz, yz, zz, irad)
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
CD1 Enter properties using geometric description.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 03-JAN-94    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/setzone.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:36   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.5   Fri Feb 02 10:35:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   08/18/95 10:31:44   llt
CD2 iad already defined, removed for cray
CD2 
CD2    Rev 1.3   08/09/95 15:21:50   gaz
CD2 changed local usage of maxit to maxitz so it is kept separate for newton
CD2 iteation max iterations
CD2 
CD2    Rev 1.2   01/28/95 13:55:44   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 16:03:56   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:36   pvcs
CD2 original version in process of being certified
CD2 
c  12/14/94 gaz changed ncord(i) to ncord(ij) 2 places
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   izone           INT      I    Zone identification number
CD3   ncord           INT      O    Node numbers of nodes in zone
CD3   nin             INT      I/O  Number of nodes in zone
CD3   xz              REAL*8   I    X coordinates defining zone
CD3   yz              REAL*8   I    Y coordinates defining zone
CD3   zz              REAL*8   I    Z coordinates defining zone
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
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
CD4   icnl            INT      faai   Problem dimension
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   w               REAL*8   fbs    Finite element shape functions
CD4   wx              REAL*8   fbs    Derivative of shape functions with
CD4                                     respect to x
CD4   wy              REAL*8   fbs    Derivative of shape functions with
CD4                                     respect to y
CD4   wz              REAL*8   fbs    Derivative of shape functions with
CD4                                     respect to z
CD4
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   Identifier      Type     Description
CD5
CD5
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   a11 - a33       REAL*8
CD5   delx y z        REAL*8
CD5   detja           REAL*8
CD5   eleb            LOGICAL  Does node belong in element
CD5   etad            REAL*8   Local coordinates in a finite element of
CD5                              the numerical integration points

CD5   excid           REAL*8   Local coordinates in a finite element of
CD5                              the numerical integration points
CD5   i               INT      Loop index
CD5   iad             INT      Loop index
CD5   ij              INT      Loop index
CD5   in              INT      Loop index
CD5   jz              INT      Loop index
CD5   maxitz          INT      Number of iterations for local coordinate
CD5                              computation
CD5   nsl             INT      Number of coordinates in element
CD5   rx ry rz        REAL*8
CD5   sa11 - sa33     REAL*8 
CD5   sid             REAL*8   Local coordinates in a finite element ofthe
CD5                              numerical integration points

CD5   sxd
CD5   sy
CD5   sz
CD5   tols            REAL*8
CD5   tolt            REAL*8
CD5   xl              REAL*8   X coordinate
CD5   xmn             REAL*8   Minimum X coordinate value
CD5   xmx             REAL*8   Maximum X coordinate value
CD5   yl              REAL*8   Y coordinate
CD5   ymn             REAL*8   Minimum Y coordinate value
CD5   ymx             REAL*8   Maximum Y coordinate value
CD5   zl              REAL*8   Z coordinate
CD5   zmn             REAL*8   Minimum Z coordinate value
CD5   zmx             REAL*8   Maximum Z coordinate value
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
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
CPS BEGIN  setzone
CPS 
CPS   IF this is a 3D problem
CPS   
CPS      FOR each coordinate position
CPS          determine maximum and minimum X, Y, and Z coordinate in zone
CPS      END FOR
CPS      
CPS      FOR each node
CPS          IF X coordinate is within zone
CPS             IF Y coordinate is within zone
CPS                IF Z coordinate is within zone
CPS                   set node number in zone
CPS                END IF
CPS             END IF
CPS          END IF
CPS      END FOR
CPS      
CPS      FOR each node in the zone
CPS          initialize local coordinates (set to zero)
CPS          FOR each iteration
CPS              call sfn3r to evaluate shape function
CPS              initialize geometric coefficients (set to zero)
CPS              FOR each coordiante position
CPS                  compute geometric coefficients
CPS              END FOR
CPS              save Jacobian and form residuals
CPS              compute corrections to local coordinates 
CPS              IF the square root of the sum of the corrections squared 
CPS               is less than the given tolerance
CPS                 EXIT iterations
CPS              END IF
CPS          END FOR
CPS          
CPS          IF iterations did not converge
CPS             print error message and stop program
CPS          END IF
CPS          
CPS          IF the local coordinates aren't within given tolerance
CPS             set within element flag to false
CPS          END IF
CPS          
CPS          IF the local coordiantes are within the element
CPS             set the zone id
CPS          ELSE
CPS             remove the node from the zone
CPS          END IF
CPS          
CPS      END FOR
CPS      
CPS   ELSE IF this is a 2D problem
CPS   
CPS      FOR each coordinate position
CPS          determine maximum and minimum X and Y coordinate in zone
CPS      END FOR
CPS      
CPS      FOR each node
CPS          IF X coordinate is within zone
CPS             IF Y coordinate is within zone
CPS                set node number in zone
CPS             END IF
CPS          END IF
CPS      END FOR
CPS      
CPS      FOR each node in the zone
CPS          initialize local coordinates (set to zero)
CPS          FOR each iteration
CPS              call sfn2r to evaluate shape function
CPS              initialize geometric coefficients (set to zero)
CPS              FOR each coordiante position
CPS                  compute geometric coefficients
CPS              END FOR
CPS              save Jacobian and form residuals
CPS              compute corrections to local coordinates 
CPS              IF the square root of the sum of the corrections squared 
CPS               is less than the given tolerance
CPS                 EXIT iterations
CPS              END IF
CPS          END FOR
CPS          
CPS          IF iterations did not converge
CPS             print error message and stop program
CPS          END IF
CPS          
CPS          IF the local coordinatesaren't within given tolerance
CPS             set within element flag to false
CPS          END IF
CPS          
CPS          IF the local coordiantes are within the element
CPS             set the zone id
CPS          ELSE
CPS             remove the node from the zone
CPS          END IF
CPS          
CPS       END FOR
CPS       
CPS   END IF
CPS   
CPS END  setzone
CPS
C***********************************************************************

      use combi
      use comdti
      use comai
      implicit none

      logical eleb
c      integer i, ij, in, izone, jz, maxitz, ncord(*), nin, nsl
      integer i, ij, in, izone, jz, maxitz,nin, nsl
      integer i_warn, irad   
      real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8 delx, dely, delz, detja
      real*8 etad, excid, rx, ry, rz, sid, sxd, sy, sz
      real*8 sa11, sa12, sa13, sa21, sa22, sa23, sa31, sa32, sa33
      real*8 tols, tolt, xl, xmn, xmx, yl, ymn, ymx, zl, zmn, zmx
      real*8 xz(*), yz(*), zz(*)

      maxitz = 10
      i_warn = 0
      tols = 1.0e-05
      tolt = 1.0e-05

      if (icnl .eq. 0.and.irad.eq.0)  then
c**** 3-d calculation ****
c**** find coordinates contained in zone ****
c**** determine element type           ****
c**** find extreme coordinates in zone ****
         xmn = 1.0e+20
         xmx = -1.0e+20
         ymn = 1.0e+20
         ymx = -1.0e+20
         zmn = 1.0e+20
         zmx = -1.0e+20
         do in = 1, nsl
            xmn = min(xz(in), xmn)
            xmx = max(xz(in), xmx)
            ymn = min(yz(in), ymn)
            ymx = max(yz(in), ymx)
            zmn = min(zz(in), zmn)
            zmx = max(zz(in), zmx)
         end do
         if(xmn.eq.xmx.or.ymn.eq.ymx.or.
     &      zmn.eq.zmx) then
            write(ierr, 6000) izone
            if (iout .ne. 0) write(iout, 6000)  izone
            if (iptty .gt. 0)  write(iptty, 6000)  izone
 6000       format(/, 1x, 'inconsistent zone coordinates',
     *           ' izone = ', i10, ' please check',
     *           ' icnl in macro CTRL ')
                        
            stop
           
         endif

c**** find nodal coordinates belonging in this zone ****
         do i = 1, neq_primary
            xl = cord(i, 1)
            if (xl .ge. xmn .and. xl .le. xmx)  then
               yl = cord(i, 2)
               if (yl .ge. ymn .and. yl .le. ymx)  then
                  zl = cord(i, 3)
                  if (zl .ge. zmn .and. zl .le. zmx)  then
                     nin = nin + 1
                     ncord(nin) = i
                  end if
               end if
            end if
         end do
c convert coordinates to work in differences
         do in = 1, nsl
          xz(in) = xz(in) -xmn
          yz(in) = yz(in) -ymn
          zz(in) = zz(in) -zmn
         enddo

c**** do calculations on nodal coordinates in this zone ****
         do ij = 1, nin
            i = ncord(ij)
c convert coordinates to work in differences
            xl = cord(i, 1) - xmn
            yl = cord(i, 2) - ymn
            zl = cord(i, 3) - zmn

c**** use newton raphson iteration to find local coordinates ****
c**** set local coordinates initially to zero ****
            sid = 0.0
            etad = 0.0
            excid = 0.0
            
            do iad = 1, maxitz
c**** evaluate shape functions ****
c              sid = max(-1.001d00,min(1.001d00,sid))
c              etad = max(-1.001d00,min(1.001d00,etad))
c              excid = max(-1.001d00,min(1.001d00,excid))
               call  sfn3r   (sid, etad, excid)

c**** form jacobian ****
c**** define geometric coefficients ****
               sxd = 0.0
               sy = 0.0
               sz = 0.0
               sa11 = 0.0
               sa12 = 0.0
               sa22 = 0.0
               sa21 = 0.0
               sa13 = 0.0
               sa23 = 0.0
               sa31 = 0.0
               sa32 = 0.0
               sa33 = 0.0

               do jz = 1, nsl
                  sxd = sxd +w (jz, 1)*xz(jz)
                  sy = sy  +w (jz, 1)*yz(jz)
                  sz = sz  +w (jz, 1)*zz(jz)
                  sa11 = sa11+wx(jz, 1)*xz(jz)
                  sa12 = sa12+wy(jz, 1)*xz(jz)
                  sa13 = sa13+wz(jz, 1)*xz(jz)
                  sa21 = sa21+wx(jz, 1)*yz(jz)
                  sa22 = sa22+wy(jz, 1)*yz(jz)
                  sa23 = sa23+wz(jz, 1)*yz(jz)
                  sa31 = sa31+wx(jz, 1)*zz(jz)
                  sa32 = sa32+wy(jz, 1)*zz(jz)
                  sa33 = sa33+wz(jz, 1)*zz(jz)
               end do

c**** save jacobian information ****
               a11 = sa22*sa33-sa23*sa32
               a21 = -sa21*sa33+sa23*sa31
               a31 = sa21*sa32-sa22*sa31
               a12 = -sa12*sa33+sa13*sa32
               a22 = sa11*sa33-sa13*sa31
               a32 = -sa11*sa32+sa12*sa31
               a13 = sa12*sa23-sa13*sa22
               a23 = -sa11*sa23+sa13*sa21
               a33 = sa11*sa22-sa12*sa21
               detja = sa11*sa22*sa33+sa12*sa23*sa31+
     .              sa13*sa21*sa32-sa13*sa22*sa31-
     .              sa23*sa32*sa11-sa33*sa21*sa12

c**** form residuals ****
               rx = sxd-xl
               ry = sy -yl
               rz = sz -zl

c**** get correction ****
               delx = (a11*rx+a12*ry+a13*rz)/detja
               dely = (a21*rx+a22*ry+a23*rz)/detja
               delz = (a31*rx+a32*ry+a33*rz)/detja
               sid = sid - delx
               etad = etad - dely
               excid = excid - delz

c**** if correction small end iteration ****
               if (sqrt(delx*delx+dely*dely+delz*delz) 
     .              .le. tols) go to 10
            end do

c**** iteration did not converge so stop ****
c write warning to fehmn.err file
            if(i_warn.eq.0) then
               if(iptty.ne.0) then
                  write(iptty,6002) izone
 6002             format('>>>> convergence problem for zone ',i8,
     &                 ' see *.err or *.out files ')
                  i_warn = 1
               endif
            endif
            if (iout .ne. 0) write(iout, 6001)  izone, i
            write(ierr, 6001) izone, i
 6001       format(1x, 'warning: iteration did ',
     &           'not converge, izone ', i6,' node ',i7)
                        

 10         continue

c**** check if local coordinates within the element ****
            eleb = .TRUE.
            if ((abs(sid) .gt. 1.0 + tolt) .or. 
     .           (abs(etad) .gt. 1.0 + tolt) .or.  
     .           (abs(excid) .gt. 1.0 + tolt))  then
               eleb = .FALSE.
            end if

c**** the node i is in zone izone, save information ****`
            if (eleb) then
               izonef(i) = izone
            else
               ncord(ij) = 0
            end if

         end do
      else if (icnl .eq. 0.and.irad.ne.0)  then
c**** 3d geometry wit hradial symmetry ****
c**** find coordinates contained in zone ****
c**** determine element type find extreme coordinates in zone ****
         xmn = 1.0e+20
         xmx = -1.0e+20
         ymn = 1.0e+20
         ymx = -1.0e+20
         do in = 1, nsl
            xmn = min(xz(in), xmn)
            xmx = max(xz(in), xmx)
            ymn = min(yz(in), ymn)
            ymx = max(yz(in), ymx)
         end do

c**** find nodal coordinates belonging in this zone ****
         do i = 1, neq_primary
c calculate radial distance (centered at x = 0)        
            xl = sqrt(cord(i, 1)**2 + cord(i, 2)**2)
            if (xl .ge. xmn .and. xl .le. xmx) then
c y is now the z direction            
               yl = cord(i, 3)
               if (yl .ge. ymn .and. yl .le. ymx) then
                  nin = nin + 1
                  ncord(nin) = i
               end if
            end if
         end do
  
c**** do calculations on nodal coordinates in this zone ****
         do ij = 1, nin
            i = ncord(ij)
            xl = sqrt(cord(i, 1)**2 + cord(i, 2)**2)
            yl = cord(i, 3)

c**** use newton raphson iteration to find local coordinates ****
c**** set local coordinates initially to zero ****
            sid = 0.0
            etad = 0.0
            
            do iad = 1, maxitz
c**** evaluate shape functions ****
               call  sfn2r   (sid, etad)

c**** form jacobian ****
c**** define geometric coefficients ****
               sxd = 0.0
               sy = 0.0
               sa11 = 0.0
               sa12 = 0.0
               sa22 = 0.0
               sa21 = 0.0

               do jz = 1, nsl
                  sxd = sxd +w (jz, 1)*xz(jz)
                  sy = sy  +w (jz, 1)*yz(jz)
                  sa11 = sa11+wx(jz, 1)*xz(jz)
                  sa12 = sa12+wy(jz, 1)*xz(jz)
                  sa21 = sa21+wx(jz, 1)*yz(jz)
                  sa22 = sa22+wy(jz, 1)*yz(jz)
               end do

c**** save jacobian information ****
               a11 = sa22
               a22 = sa11
               a12 = -sa12
               a21 = -sa21
               detja = sa11*sa22-sa12*sa21

c**** form residuals ****
               rx = sxd-xl
               ry = sy -yl

c**** get correction ****
               delx = (a11*rx+a12*ry)/detja
               dely = (a21*rx+a22*ry)/detja
               sid = sid -delx
               etad = etad-dely

c**** if correction small end iteration ****
               if (sqrt(delx*delx+dely*dely) .le. tols) go to 20

            end do

c**** iteration did not converge so stop ****
            write(ierr, 6010)  izone
            if (iout .ne. 0) write(iout, 6010)  izone
            if (iptty .gt. 0)  write(iptty, 6010)  izone
 6010       format(/, 1x, 'zone calcs did not converge',
     *           ' izone = ', i10, ' please check',
     *           ' coordinates ')            
            stop
            
 20         continue

c**** check if local coordinates within the element ****
            eleb = .TRUE.
            if (abs(sid) .gt. 1.0+tolt)  eleb = .FALSE.
            if (abs(etad) .gt. 1.0+tolt)  eleb = .FALSE.

c**** the node i is in zone izone, save information ****
            if (eleb)  then
               izonef(i) = izone
            else
               ncord(ij) = 0
            end if
            
         end do
      else if (icnl .ne. 0)  then
c**** 2-d calculation ****
c**** find coordinates contained in zone ****
c**** determine element type find extreme coordinates in zone ****
         xmn = 1.0e+20
         xmx = -1.0e+20
         ymn = 1.0e+20
         ymx = -1.0e+20
         do in = 1, nsl
            xmn = min(xz(in), xmn)
            xmx = max(xz(in), xmx)
            ymn = min(yz(in), ymn)
            ymx = max(yz(in), ymx)
         end do

c**** find nodal coordinates belonging in this zone ****
         do i = 1, neq_primary
            xl = cord(i, 1)
            if (xl .ge. xmn .and. xl .le. xmx) then
               yl = cord(i, 2)
               if (yl .ge. ymn .and. yl .le. ymx) then
                  nin = nin + 1
                  ncord(nin) = i
               end if
            end if
         end do
  
c**** do calculations on nodal coordinates in this zone ****
         do ij = 1, nin
            i = ncord(ij)
            xl = cord(i, 1)
            yl = cord(i, 2)

c**** use newton raphson iteration to find local coordinates ****
c**** set local coordinates initially to zero ****
            sid = 0.0
            etad = 0.0
            
            do iad = 1, maxitz
c**** evaluate shape functions ****
               call  sfn2r   (sid, etad)

c**** form jacobian ****
c**** define geometric coefficients ****
               sxd = 0.0
               sy = 0.0
               sa11 = 0.0
               sa12 = 0.0
               sa22 = 0.0
               sa21 = 0.0

               do jz = 1, nsl
                  sxd = sxd +w (jz, 1)*xz(jz)
                  sy = sy  +w (jz, 1)*yz(jz)
                  sa11 = sa11+wx(jz, 1)*xz(jz)
                  sa12 = sa12+wy(jz, 1)*xz(jz)
                  sa21 = sa21+wx(jz, 1)*yz(jz)
                  sa22 = sa22+wy(jz, 1)*yz(jz)
               end do

c**** save jacobian information ****
               a11 = sa22
               a22 = sa11
               a12 = -sa12
               a21 = -sa21
               detja = sa11*sa22-sa12*sa21

c**** form residuals ****
               rx = sxd-xl
               ry = sy -yl

c**** get correction ****
               delx = (a11*rx+a12*ry)/detja
               dely = (a21*rx+a22*ry)/detja
               sid = sid -delx
               etad = etad-dely

c**** if correction small end iteration ****
               if (sqrt(delx*delx+dely*dely) .le. tols) go to 25

            end do

c**** iteration did not converge so stop ****
            write(ierr, 6015)  izone
            if (iout .ne. 0) write(iout, 6010)  izone
            if (iptty .gt. 0)  write(iptty, 6010)  izone
 6015       format(/, 1x, 'zone calcs did not converge',
     *           ' izone = ', i10, ' please check',
     *           ' coordinates ')            
            stop
            
 25         continue

c**** check if local coordinates within the element ****
            eleb = .TRUE.
            if (abs(sid) .gt. 1.0+tolt)  eleb = .FALSE.
            if (abs(etad) .gt. 1.0+tolt)  eleb = .FALSE.

c**** the node i is in zone izone, save information ****
            if (eleb)  then
               izonef(i) = izone
            else
               ncord(ij) = 0
            end if
            
         end do
      end if

      return
      end
