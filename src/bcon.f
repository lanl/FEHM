      subroutine bcon (iz)
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
CD1 Manage boundary conditions.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 06-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/bcon.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:52   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Thu Dec 18 15:21:56 1997   gaz
CD2 corrected outflow energy calculation (skipped
CD2 over when apprpriate)
CD2 
CD2    Rev 1.7   Fri Sep 26 15:12:42 1997   llt
CD2 gaz changes
CD2 
CD2    Rev 1.6   Fri Apr 26 15:04:30 1996   gaz
CD2 changes to replenish volume reservoirs for BC s
CD2 
CD2    Rev 1.5   Mon Jan 29 13:15:22 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.4   Mon Jan 29 11:47:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.3   05/30/95 16:44:48   robinson
CD2 Corrected solute mass balance calculation
CD2 
CD2    Rev 1.2   04/03/95 08:43:50   robinson
CD2 Correction to solute mass balance calculation
CD2 
CD2    Rev 1.1   03/18/94 15:53:12   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:26   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   iz              INT      I    Flag to indicate type of boundary
CD3                                   condition modification to be performed
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
CD4   dencj           REAL*8   fcc    Tracer accumulation term
CD4   denei           REAL*8   fcc    Energy accumulation term
CD4   deni            REAL*8   fcc    Mass accumulation term
CD4   denpci          REAL*8   co2    Gas accumulation term at each node
CD4   dtotc           REAL*8   faar   Tracer time step size in seconds
CD4   dtotdm          REAL*8   faar   Last time step size in seconds
CD4   ieos            INT      fddi1  Phase state of fluid at each node
CD4   ka              INT      fbb    Contains boundary type information for
CD4                                     each node
CD4   n               INT      faai   Total number of nodes
CD4   npn             INT      fdd1i  Parameter used in storing tracer data
CD4   nsp             INT      fdd1i  Current species number
CD4   qcout           REAL*8   fdd1   Total produced tracer mass for each
CD4                                     species
CD4   qt              REAL*8   faar   Total outflow for time step
CD4   qtc             REAL*8   co2r   Total mass injected in source term
CD4   qte             REAL*8   faar   Total energy outflow for time step
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   sx1             REAL*8   fbc    Volume associated with each node
CD4   volume          REAL*8   fdd    Volume associated at each node
CD4
CD4 Global Subprograms
CD4
CD4   None
CD4
C***********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   Identifier      Type     Description
CD5
CD5   sx1bc           REAL*8   /1.e12/
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5   cqtout          REAL*8   Total produced tracer mass for each species at
CD5                              the boundary
CD5   cqtin           REAL*8   Total lost tracer mass for each species at
CD5                              the boundary
CD5   cqtrxn          REAL*8   Total produced tracer mass for each species at
CD5                              the boundary (reaction)
CD5   i               INT      Loop index
CD5   j               INT      Loop index
CD5   qtb             REAL*8   Total outflow for time step at the boundary
CD5   qtcb            REAL*8   Total mass injected at the boundary
CD5   qteb            REAL*8   Total energy outflow for time step at the
CD5                              boundary
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
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
CD9 2.3.7 Sources and sinks
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
CPS BEGIN bcon
CPS 
CPS   IF case 1?
CPS   
CPS      FOR each node
CPS          IF the node has fixed boundary conditions
CPS             set values for the node
CPS          END IF
CPS      ENDFOR
CPS      
CPS   ELSE IF case 2?
CPS   
CPS      initialize boundary flow values (mass and energy) to zero
CPS      FOR each node
CPS          IF the node is a boundary node
CPS             compute flow values for the node and add them to the  . . .
CPS             . . . boundary flow value
CPS          ENDIF
CPS      ENDFOR
CPS      adjust the total flow values using the boundary flow values
CPS      
CPS   ELSE IF case 3?
CPS      initialize boundary flow value to zero
CPS      FOR each node
CPS          IF the node is a boundary node
CPS            IF there is solute mass flowing into this node
CPS              Add to running total of mass into the node
CPS            ELSE there is solute mass flowing out of this node
CPS              Add to running total of mass out of the node
CPS            ENDIF
CPS            Add to source/sink reaction term for this node
CPS          ENDIF
CPS      ENDFOR
CPS      Adjust source, sink, and reaction mass balance terms
CPS      
CPS   END IF
CPS 
CPS END bcon
CPS
C***********************************************************************

      use combi
      use comci
      use comdi
C      use comei         
      use comfi
      use comgi
      use davidi
      use comdti
      use comai
      use comsplitts
C      use comflow      
      implicit none

      integer i, iz, j
      real*8 qtb, qtcb, qteb, sx1bc, rxn_term, store_term
      real*8 cqtin, cqtout, aiped_big, bp1_save, bp2_save

      parameter (sx1bc=1.e12)

      if (iz .eq. 1) then
c
c adjust volumes and heat capacities
c
c identify ieos(i)<0 as fixed bc
         if(ice.eq.0) then
            do i = 1, n
               if (ieos(i) .lt. 0) then
                  ka(i) = 0
                  sk(i) = 0
                  sx1(i) = sx1(i) / abs(sx1(i)) * sx1bc
                  ieos(i) = abs(ieos(i))
                  volume(i) = 0.0
               end if
            end do
         end if

      else if (iz .eq. 2) then
c
c remove boundary nodes from mass and energy calculations
c

         qtb  = 0.0
         qteb = 0.0
         qtcb = 0.0
         do i = 1, n
            if (volume(i) .eq. 0.0) then
               qtb = qtb + deni(i) * sx1(i) * dtotdm / sx1bc        
               denh(i) = denh(i) - deni(i)*dtot
               if(ifree.eq.0.and.irdof.ne.13) then
                qteb = qteb + denei(i) * sx1(i) * dtotdm / sx1bc
                deneh(i) = deneh(i) - denei(i)*dtot
               endif
               if(ico2.gt.0) then
                qtcb = qtcb + denpci(i) * sx1(i) * dtotdm / sx1bc
                denpch(i) = denpch(i) - denpci(i)*dtot
               endif
            end if
         end do
         qt   = qt   + qtb * sx1bc
         qte  = qte  + qteb * sx1bc
         qtc  = qtc  + qtcb * sx1bc

      else if (iz .eq. 3) then
c
c remove boundary nodes from tracer calculations
c
         cqtin = 0.0
         cqtout = 0.0
c Corrected solute mass balance calculations
         do j = 1, n
            if (volume(j) .eq. 0.0) then
               i = j + npn
               rxn_term = rc(i)
               store_term = dencj(i) - rxn_term
               if( store_term .lt. 0. ) then
                  cqtin = cqtin + store_term
     2                 * sx1(j) * dtotc / sx1bc
               else
                  cqtout = cqtout + store_term
     2                 * sx1(j) * dtotc / sx1bc
               end if
            end if
         end do
         qcout(nsp) = qcout(nsp) + cqtout
         qcin(nsp) = qcin(nsp) + cqtin

      else if (iz .eq. 4) then
c
c adjust source -sink terms for large aiped      1
c
         if (ifree .eq. 0) then
            aiped_big = 1.e4
         else
            aiped_big = 1.e20
         end if

         if(ice.eq.0.and.ianpe.eq.0) then
c  only perform calcs with non-clathrate simulation
c  only perform calcs with with isotropic permeability
c  need to modify fluxnet for clathrate simulation
            if(ico2.lt.0.and.idpdp.eq.0.and.idualp.eq.0) then
               do i=1,n
c gaz 082123 added more flo3 conditions
                  if(wellim(i).gt.aiped_big.and.ka(i).gt.-2) then
                     bp1_save=bp(i+nrhs(1))
                     bp2_save=bp(i+nrhs(2))
                     bp(i+nrhs(1))=0.0
                     bp(i+nrhs(2))=0.0
                     call flux_net(i)  
                     sk(i)=-bp(i+nrhs(1))
                     if(irdof.ne.13) then
                        qh(i)=-bp(i+nrhs(2))
                     end if
                     bp(i+nrhs(1))=  bp1_save
                     bp(i+nrhs(2))=  bp2_save
                  endif
               enddo
            endif
         end if
      end if

      end
