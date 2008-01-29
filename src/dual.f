      subroutine  dual( inum )
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
CD1 To provide overall control for a dual porosity solution.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-4-93      G. Zyvoloski   N/A     Initial implementation, but
CD2                                        previous non-YMP versions
CD2                                        of FEHM exist, and the
CD2                                        current version may differ
CD2                                        from these
CD2
CD2 $Log:   /pvcs.config/fehm90/src/dual.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:54   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:04   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:38   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Wed Jun 12 15:21:00 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.6   Mon Jun 03 11:17:46 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.5   Fri May 31 10:52:10 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.4   Thu Feb 15 10:07:58 1996   zvd
CD2 Corrected purpose.
CD2 
CD2    Rev 1.3   Mon Jan 29 15:26:16 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   01/28/95 13:54:28   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.1   03/18/94 15:50:02   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:23:20   pvcs
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
CD3 inum         int     I       Flag denoting the reason for calling
CD3                                 the routine
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 None
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4 
CD4 idualp, inpt, wdd1, n, volf1, volf2, rp17f, n0, izonef, apuv1,
CD4 neq, irlp, irlpt, sx1, cord, pnx, pny, pnz, thx, thy, thz, nskw,
CD4 m, m2, nskw2, ico2, icnl, zero_t, narrays, pointer, itype, default,
CD4 igroup,  ireturn,
CD4 macroread
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
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4 
CD4 
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 null            N/A      Determine if there is additional data to
CD4                             read
CD4 varchk          N/A      Initializes EOS parameters
CD4 dualfh          N/A      Compute heat and mass solution for a
CD4                             single phase problem
CD4 dualfa          N/A      Compute heat and mass solution for a
CD4                             two phase problem
CD4 dualex          N/A      Extract heat and mass solution
c no longer in dual, called directly from cnswer 8/11/94
CD4 dualta          N/A      Compute tracer solution
CD4 dualtx          N/A      Extract tracer solution
CD4 initdata        N/A      Reads and sets data input by zone or node
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 None
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop index
CD5 fracm        real*8      Volume fraction of matrix
CD5 frac0        real*8      Volume fraction of primary nodes
CD5 frac1        real*8      Volume fraction of first matrix nodes
CD5 frac2        real*8      Volume fraction of second matrix nodes
CD5 sx1d         real*8      Volume of current node
CD5 icnls        int         Temporary storage to save icnl variable
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
CD7 None
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
CD9 2.4.7 Dual-porosity formulation
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
CPS BEGIN dual
CPS 
CPS IF this is not a dual porosity simulation
CPS   ERROREXIT
CPS ENDIF
CPS 
CPS IF this call is for reading input
CPS 
CPS   Read dual porosity flag
CPS   
CPS   Set parameter values for reading volumes of primary nodes
CPS     
CPS   initdata - read and set values of primary nodes
CPS   
CPS   Set parameter values for reading volumes of secondary nodes
CPS   initdata - read and set values of secondary nodes
CPS   
CPS   Set parameter values for reading area per unit volume parameter
CPS   initdata - read and set values of area per unit volume
CPS   
CPS   Set flag to denote that the dual macro has been called
CPS       
CPS ENDIF this call is for reading input
CPS 
CPS IF this is a call to initialize EOS parameter values
CPS   varchk - initialize EOS parameter values for first set of...
CPS   ... matrix nodes
CPS   varchk - initialize EOS parameter values for second set of...
CPS   ... matrix nodes
CPS ENDIF
CPS 
CPS IF this is a call to modify fracture volumes
CPS 
CPS   IF there was insufficient space allocated for dual porosity
CPS     Write error message
CPS     EXIT the program
CPS   ENDIF
CPS   
CPS   FOR each node
CPS   
CPS     IF van Ganuchten relative permeabilities are used
CPS       Calculate volume fractions from fracture porosity
CPS     ENDIF
CPS     
CPS     IF the primary node volume fraction is less than zero
CPS       Write error message
CPS       EXIT the program
CPS     ENDIF
CPS     
CPS     IF the volume fraction of the second matrix node is less...
CPS     ... than zero
CPS       Set volume fraction of first matrix node to correct value...
CPS       ... such that second node is zero
CPS     ENDIF
CPS   
CPS     Compute volume associated with primary and matrix nodes
CPS     Set parameter to eliminated gravity in matrix nodes
CPS     
CPS   ENDFOR each node
CPS   
CPS   FOR each node
CPS   
CPS     Initialize parameters so that permeabilities are...
CPS     ... proportional to volume fractions
CPS     
CPS     IF permeabilities are not to be proportional to volume fraction
CPS       Set parameters so that permeabilities are the input values
CPS       IF there is no volume in the second set of matrix nodes
CPS         Set parameter for permeability of second set of nodes
CPS       ELSE
CPS         Set parameter to zero out permeability of second matrix node
CPS       ENDIF
CPS     ENDIF
CPS     
CPS     Set permeability and conductivity values for matrix nodes
CPS     Modify permeability and conductivity values for primary nodes
CPS     
CPS   ENDFOR each node
CPS   
CPS ENDIF
CPS 
CPS IF this call is to add matrix node numbers to output node numbers
CPS   Compute total number of nodes
CPS     
CPS   FOR each node for which output is written
CPS     Compute node number for matrix nodes, store in array
CPS   ENDFOR
CPS     
CPS   FOR each node for which output is written to dp file
CPS     Compute node number for matrix nodes, store in array
CPS   ENDFOR
CPS     
CPS ENDIF
CPS 
CPS IF this is a call to obtain a heat and mass transfer solution
CPS 
CPS   IF this is a single phase problem
CPS     dualfh - compute heat and mass solution
CPS   ELSEIF this is a two phase solution
CPS     dualfa - compute heat and mass solution
CPS   ENDIF
CPS   
CPS ENDIF
CPS 
CPS IF this is a call to extract the solution for dual porosity nodes
CPS   dualex - extract solution for dual porosity node
CPS ENDIF
CPS 
CPS IF this is a call to obtain a tracer solution
CPS   dualta - compute the tracer solution
CPS ENDIF
CPS 
CPS IF this is a call to extract the tracer solution for dual...
CPS ... porosity nodes
CPS   dualtx - extract tracersolution for dual porosity node
CPS ENDIF
CPS 
CPS END dual
CPS 
C**********************************************************************
changes to version FEHM5.1J
c 17-sept-91
c put zone input in for dual porosity
c 1-oct-91
c read in idualp in input for dual porosity
c 1-oct-91
c corrected mistake introduced sept-17,effect volf1=volf2
c 1-oct-91
c put in capability for volf2=0.(1-zone  dual porosity)
c multiplied flow parameters bt fracture volume
c 10-oct-91
c fixed error in saving sx12s in dualta
c 17-oct-91
c printout of dual por level nodes automatically
c 28-Nov-91
c started changing to the symmetry format
c 2-dec-91
c changed to symmetry format in dualfa(air only)
c 3-dec-91
c changed def of volf1 to refer to fracture
c 5-dec-91
c added dualfh for heat transfer
c 16-dec-91
c corrected dist01 to dist12 in dualfa and dualfh
c deleted r2e-a33ee lines were duplicated in dualfa,dualha
c added accumulation terms to matrix eqs(oops!)
c defined sx1d in dualfa
c removed zeroing of transfer terms from daulfa dualfh
c zeroed t6 and t7 in dualfa dualfh
c changed frac2.le.0.0 to frac2.le.tmch
c corrected iniyial volume calcs in dual
c made perms for second matrix node idll values
c zeroed a11mp etc(oops!)
c put changes of air model into heat/mass
c being careful with accumulation terms
c 17-dec-91
c changed over tracer dual solution
c 19-mar-92
c defined swi in dualta
c 26-may-92
c made t9=fid always
c 11-june-92
c changed codeing subroutine dual so volume in matrix nodes can be zero
c checked for enough storage for dual porosity in dual
c 17-june-92
c commented out 2nd def of coef1(duplicate)
c 22-july-92
c put stop statement in dual if n>n0
c 31-july-92
c commented out sx4d,sx4h definitions,no numerical change
c  20-nov-92
c if frac porosity available from rlperm input,use for frac volume sub dual
c 1-feb-93
c changed nb to nb0 in dualfh
c 23-feb-93 llt
c changed equivalences to pointers
C***********************************************************************

      use comhi
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comki
      implicit none

      integer inum,i,icnls
      real*8 fracm,frac0,frac1,frac2,sx1d

      if ( idualp .eq. 0)  go to 999
      if ( inum .eq. 0 )  then
c     
c     read in dual porosity parameter idualp
c     
         read(inpt,*) idualp
c     
c**** read input for dual porosity   ****
c**** read volume fractions of nodes ****
c     
         narrays = 1
         itype(1) = 8
         default(1) = 0.001
         macro = "dual"
         igroup = 2
         
         call initdata2( inpt, ischk, neq, narrays,
     &        itype, default, macroread(3), macro, igroup, ireturn,
     2        r8_1=volf1(1:neq) )
         
         default(1) = 0.5
         igroup = 3
         call initdata2( inpt, ischk, neq, narrays,
     &        itype, default, macroread(3), macro, igroup, ireturn,
     2        r8_1=volf2(1:neq) )
         
         default(1) = 5.
         igroup = 4
         call initdata2( inpt, ischk, neq, narrays,
     &        itype, default, macroread(3), macro, igroup, ireturn,
     2        r8_1=apuv1(1:neq) )
         
         macroread(3) = .TRUE.
      
      endif

      if ( inum .eq. 1 )  then
         call  varchk  ( 0,neq     )
         call  varchk  ( 0,neq+neq )
      endif

      if ( inum .eq. 2 )  then
c**** modify fracture volume at nodes for dual porosity calcs ****
c**** set cord(i,3) for dual porosity nodes ****
c     
c     check for enough storage
         if (n.gt.n0) then
            write(ierr, 100) n, n0
            if (iout .ne. 0) write(iout, 100) n, n0
            if(iatty.ne.0) write(iatty, 100) n, n0
            stop
 100        format ('***** n=', i8,' > n0= ', i8, 'stopping ****')
         endif
         do i=1,neq
c     check if van Genutchen rl perms were selected
            if (rlp_flag .eq. 1) then
               if(irlpt(irlp(i)).eq.4) then
c     calculate volume factor from fracture porosity
                  volf1(i)=rp17f(irlp(i))
                  fracm=1.0-volf1(i)
                  volf2(i)=0.1*fracm
               endif
            end if
            frac0=volf1(i)
            frac1=volf2(i)
            frac2=1.0-frac0-frac1
            if(frac0.le.0.0) then
c     
c     bad parameters,stop calculation
c     
               write(ierr, 200) 
               write(ierr, 210)                
               if (iout .ne. 0) then
                  write(iout, 200) 
                  write(iout, 210) 
               end if
               if(iatty.ne.0) then
                  write(iatty,*)
                  write(iatty,*)
               endif
 200           format ('**** check fracture volumes,stopping **** ')
 210           format ('**** check equivalent continuum VGs **** ')
               stop

            endif
            if(frac2.lt.0.0) then
               frac1=1.0-frac0
               volf2(i)=frac1
            endif
            sx1d=sx1(i)
            sx1 (i        ) =  sx1d*frac0
            sx1 (i+neq    ) =  sx1d*frac1
            sx1 (i+neq*2  ) =  sx1d*( 1.0-frac0-frac1 )
            cord(i+neq  ,3) =  1.
            cord(i+neq*2,3) =  1.
         enddo
c     
c**** modify matrix permeability and conductivity to be ****
c**** maximum of directional values ****
c     
         do i=1,neq
            frac0=volf1(i)
            frac1=volf2(i)
            frac2=1.0-frac0-frac1
c     
c     if idualp=2 then set permeabilities to k*fraction
c     
            if(idualp.ne.2) then
               frac0=1.0
               frac1=1.0
               if(frac2.gt.0.0) then
                  frac2=1.0
               else
                  frac2=0.0
               endif
            endif
            pnx(i+neq  ) =  frac1*max( pnx(i+neq  ),pny(i+neq  ),
     *           pnz(i+neq  ),zero_t )
            thx(i+neq  ) =  frac1*max( thx(i+neq  ),thy(i+neq  ),
     *           thz(i+neq  ),zero_t )
            pnx(i+neq*2) =  frac2*max( pnx(i+neq*2),pny(i+neq*2),
     *           pnz(i+neq*2),zero_t )
            thx(i+neq*2) =  frac2*max( thx(i+neq*2),thy(i+neq*2),
     *           thz(i+neq*2),zero_t )
            pny(i+neq  ) =  zero_t
            pny(i+neq*2) =  zero_t
            pnz(i+neq  ) =  zero_t
            pnz(i+neq*2) =  zero_t
            thy(i+neq  ) =  zero_t
            thy(i+neq*2) =  zero_t
            thz(i+neq  ) =  zero_t
            thz(i+neq*2) =  zero_t
c     
c     modify fracture nodes
c     
            pny(i)=pny(i)*frac0
            pnz(i)=pnz(i)*frac0
            pnx(i)=pnx(i)*frac0
            thx(i)=thx(i)*frac0
            thy(i)=thy(i)*frac0
            thz(i)=thz(i)*frac0
            thy(i)=thy(i)*frac0
            thz(i)=thz(i)*frac0
         enddo

      endif

      if (inum.eq.3) then

         n=neq*3
c     
c     add nodes to nskw and nskw2
c     
         do i=1,m
            nskw(m+i)=nskw(i)+neq
            nskw(m+m+i)=nskw(i)+neq+neq
         enddo
         m=3*m
         do i=1,m2
            nskw2(m2+i)=nskw2(i)+neq
            nskw2(m2+m2+i)=nskw2(i)+neq+neq
         enddo
         m2=3*m2

      endif

      if ( inum .eq. 10 )  then
c     
c**** factor in dual porosity contribution(heat and mass) ****
c     
c         icnls  =  icnl
c         icnl   =  1
         if(ico2.eq.0) then
            call dualfh
         else if(ico2.lt.0) then
            call  dualfa
         endif
c         icnl   =  icnls

      endif
c     
c**** extract solution for dual porosity nodes(heat and mass) ****
c     
      if ( inum .eq. 11 )  then

         call  dualex

      endif
c     changed 8/11/94 to incorporate coupling
c     dual(12) is no longer called.  Dualta is called from cnswer
c     
c**** extract solution for dual porosity nodes(tracer) ****
c     
      if ( inum .eq. 13 )  then

         call  dualtx

      endif

 999  continue
      
      return
      end


