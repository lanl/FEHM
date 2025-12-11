      subroutine  wrtcon(iz)
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
CD1 To write tracer simulation information to output files.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-17-93     G. Zyvoloski   00022   Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/wrtcon.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Mon Jun 10 14:59:36 1996   hend
CD2 Added trac output per time
CD2 
CD2    Rev 1.7   Fri May 24 11:56:16 1996   hend
CD2 Removed excess column & kg-->mol
CD2 
CD2    Rev 1.6   Wed May 08 14:16:42 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.5   Fri Feb 02 14:30:40 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   09/29/95 16:12:36   llt
CD2 added small number on dividing -- was dividing by zero
CD2 
CD2    Rev 1.3   08/16/95 16:27:02   robinson
CD2 Revised to write either node information or mass balance information
CD2 
CD2    Rev 1.2   01/28/95 14:21:14   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.1   03/18/94 16:16:06   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:34   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier   Type     Description
CD3
CD3 iz           INTEGER  Flag to determine which output to write
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                Description
CD3
CD3 tape number iatty    File number associated with screen output
CD3 tape number iout     File number used for output if desired by user
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4 Global Types
CD4
CD4 None
CD4
CD4 Global Variables
CD4
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 dtotc, ntty, nsp, nspeci, npn, npt, iout, iadd, iatty, m, nskw, neq,
CD4 rc, icns, ico2, sk, s, qh, an, anl, anv, bp, rcss, qcrxn
CD4 
CD4 Global Subprograms
CD4
CD4 Name    Type     Description
CD4 
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 idp          int         Flag used in writing to file for dual or
CD5                          dpdp solution
CD5 idq          int         Flag used in writing to file for dual or
CD5                          dpdp solution
CD5 i            int         Do loop index
CD5 mdd          int         Temporary storage for node number
CD5 rcd          real*8      Source/sink and reaction term at current node
CD5 rcdss        real*8      Source/sink term at current node
CD5 md           int         Temporary storage for tracer unknown number
CD5 rcc          real*8      Temporary storage of source/sink
CD5 srmim        real*8      Source/sink flow rate at current node
CD5 
CD5 Local Subprograms
CD5 
CD5 None
CD5 
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following functions are carried out in this routine:
CD6 
CD6   Convert the time step size from seconds to days;
CD6   
CD6   Then, if the information is to be written at all for this time:
CD6   
CD6   Loop through each solute, writing the species number, time step
CD6   and a heading;
CD6   
CD6     Then, in an inner loop over all nodes for which information is
CD6     written, write the concentration and residual information;
CD6     
CD6   Then, to complete the information for this time step, write the
CD6   overall mass balance information;
CD6   
CD6   Return to calling routine.
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
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN wrtcon
CPS 
CPS Compute time step in days
CPS 
CPS IF information is to be written this time step
CPS 
CPS   FOR each solute
CPS   
CPS     Write header, time in days
CPS     Write solute number, time step, number of iterations
CPS     
CPS     IF information is to be written to file iatty
CPS       Write solute number, time step, number of iterations
CPS     ENDIF
CPS     
CPS     IF we are writing nodal information
CPS       FOR each node at which we are writing information
CPS     
CPS         IF we are writing information for matrix level 1
CPS           Write heading
CPS         ELSEIF we are writing information for matrix level 2
CPS           Write heading
CPS         ENDIF
CPS       
CPS         IF the current solute is a liquid phase solute
CPS           IF this is an air-water system
CPS             Assign fluid source term
CPS           ELSE this is not an air water system
CPS             Compute fluid source term
CPS           ENDIF
CPS         ELSE  the solute is a gas phase solute
CPS           IF this is an air-water system
CPS             Assign fluid source term
CPS           ELSE this is not an air-water system
CPS             Compute fluid source term
CPS           ENDIF
CPS         ENDIF
CPS       
CPS         IF there is a finite fluid mass input or output
CPS           Compute solute source or sink term
CPS         ENDIF
CPS       
CPS         Write node number, concentrations, source/sink values, ...
CPS         ... and residuals for this node
CPS     
CPS         IF information is to be written to file iatty
CPS           Write node number, concentrations, source/sink values,...
CPS           ... and residuals for this node
CPS         ENDIF
CPS 
CPS       ENDFOR each node at which we are writing information
CPS
CPS     ENDIF we are writing nodal information
CPS
CPS     Compute global mass balance
CPS   
CPS     Write overall mass balance information
CPS     
CPS     IF information is to be written to file iatty
CPS       Write overall mass balance information
CPS     ENDIF
CPS 
CPS   ENDFOR each solute
CPS 
CPS ENDIF information is to be written this time step
CPS
CPS END wrtcon
CPS 
C**********************************************************************

      use comai
      use combi
      use comdi
      use comdti
      use comgi
      use comrxni, only : scl, scv
      use comxi
      use comchem, only: species

      implicit none

      integer idp,idq,i,mdd,md,iz,sehi,sehindex, jjj, kkk, lll
      integer nts_trac
      real*8 rcd,rcdss,rcc,srmim,cbal,sehratein,sehrateout
      
c recalculate tracer timesteps
      nts_trac = dtotdm/dtotc
      
c**** write out output to appropriate device ****  
      dtotc  =  dtotc/86400.0

c**** print output for each tracer ****
c changing print interval for tracer information, only print at 
c HM timestep
c for Don Neeper  PHS  4-13-05
c-----------------------------------------------------------------------
c      if(ntty.eq.2 .or. iz .eq. 0 ) then
      if((ntty.eq.2).AND.(iz.eq.1)) then
c         if( iz .eq. 0 ) then
         if (iout .ne. 0) then
            write(iout,31) 
            write(iout,6050) days
            write(iout  ,6001)  nts_trac, dtotc ,  iadd(1), iaddt(1)
         end if
 31      format(/,19x,'*************************')
         if ( iatty .gt. 0 )  then
            write(iatty  ,6001)  nts_trac, dtotc ,  iadd(1), iaddt(1)
         endif
c         end if
         do nsp=1,nspeci
            npn=npt(nsp)
            if (iz.ne.1) then
               if (iout .ne. 0) write(iout,6000) species(nsp)
               if (iatty.gt.0) then
                  write(iatty,6000) species(nsp) 
               end if
            else if( iz. eq. 1 ) then
c If there are nodes to output
               if (m .gt. 0) then
                  if (iout .ne. 0) then
                     write(iout,777) 
                     write(iout,778) species(nsp)
c     write(iout,6010)
                     if (iadsfl(nsp,itrc(nskw(1))).ne.0.or.
     &                    iadsfv(nsp,itrc(nskw(1))).ne.0) then
                        write(iout,779)
                        write(iout,6012)
                     else
                        write(iout,679) 
                        write(iout,6610)
                     end if
                  end if
                  if ( iatty .gt. 0 )  then
                     write(iatty,777)
                     write(iatty,6000) species(nsp)
c     write(iatty,6010)
                     if (iadsfl(nsp,itrc(nskw(1))).ne.0.or.
     &                    iadsfv(nsp,itrc(nskw(1))).ne.0) then 
                        write(iatty,779) 
                        write(iatty,6012)
                     else
                        write(iatty,679) 
                        write(iatty,6610)
                     end if
                  end if
 777              format(/,20x,'Nodal Information (Tracer)')
                  idp=0
                  idq=0
                  do i=1,m
                     mdd    =  nskw(i)
c     write headings for matrix levels
                     if(mdd.gt.neq.and.mdd.le.neq+neq.and.idp.eq.0) then
                        if (iout .ne. 0) write(iout,*) 
     &                       ' matrix level = 1'
                        if(iatty.ne.0) write(iatty,*) 
     &                       ' matrix level = 1'
                        idp=1
                     else if(mdd.gt.neq+neq.and.idq.eq.0) then
                        if (iout .ne. 0) write(iout,*) 
     &                       ' matrix level = 2'
                        if(iatty.ne.0) write(iatty,*) 
     &                       ' matrix level = 2'
                        idq=1
                     endif
                     md     =  mdd+npn
                     rcd    =  rc(md)
                     rcdss = rcss(md)
                     rcc    =  0.0
                     if(icns(nsp).gt.0) then
                        if(ico2.lt.0) then
                           srmim=sk(mdd)
                        else
                           srmim=sk(mdd)*s(mdd)
                        endif
                     else
c     gas phase tracer
                        if( ico2 .lt. 0 ) then
                           srmim=qh(mdd)
                        else
                           srmim=sk(mdd)*(1.0-s(mdd))
                        endif
                     endif
                     if ( abs( srmim ) .gt. 1.0d-15 )  then
                        rcc  =  rc(md)/(sk(mdd)+1e-30)
                     end if
                  
c                  if (iout .ne. 0) 
c     &                 write(iout  ,6011)  mdd, an(md), anl(md),
c     &                 anv(md), rcdss, bp(mdd)
c                  if ( iatty .gt. 0 )
c     &                 write(iatty ,6011)  mdd, an(md), anl(md), 
c     &                 anv(md), rcdss, bp(mdd)
                     if (iout .ne. 0) then
                        if (iadsfl(nsp,itrc(mdd)).ne.0.or.
     &                       iadsfv(nsp,itrc(mdd)).ne.0) then
                           write(iout  ,6013)  mdd , an(md) , anl(md),
     &                          anv(md), scl(md), scv(md), rcdss,bp(mdd)
                        else
                           write(iout  ,6611)  mdd, an(md), anl(md),
     &                          anv(md), rcdss, sinkint(md), bp(mdd)
                        end if
                     end if
                     if ( iatty .gt. 0 ) then
                        if (iadsfl(nsp,itrc(mdd)).ne.0.or.
     &                       iadsfv(nsp,itrc(mdd)).ne.0) then
                           write(iatty  ,6013)  mdd , an(md) , anl(md),
     &                          anv(md), scl(md), scv(md), rcdss,bp(mdd)
                        else
                           write(iatty ,6611)  mdd, an(md), anl(md), 
     &                          anv(md), rcdss, sinkint(md), bp(mdd)
                        end if
                     end if
                  end do
                  if ( iout .ne. 0 ) write(iout, *)
                  if ( iatty .gt. 0 ) write(iatty, *)
               end if
            end if

c            if( iz .eq. 0 ) then
            if( iz .eq. 1 ) then
               cbal  =  cm(nsp) - cm0(nsp) + qcout(nsp) + qcin(nsp) +
     2              qcrxn(nsp)
               sehratein=0.
               sehrateout=0.
               do sehi=1,n0
                  sehindex=sehi+npt(nsp)
                  if (rcss(sehindex).le.0.) then
                     sehratein=sehratein-rcss(sehindex)
                  else
                     sehrateout=sehrateout+rcss(sehindex)
                  endif
               enddo
               
               if (iout .ne. 0) then
                  if (m .eq. 0) write (iout, 6000) species(nsp)
                  write(iout,6020)  cm0(nsp),cm(nsp),abs(qcin(nsp)),
     &                 sehratein,qcout(nsp),sehrateout,-qcrxn(nsp),cbal
               end if
               if ( iatty .gt. 0 ) then
                  if (m .eq. 0) write (iatty, 6000) nsp
                  write(iatty,6020) cm0(nsp),cm(nsp),abs(qcin(nsp)),
     &                 sehratein,qcout(nsp),sehrateout,-qcrxn(nsp),cbal
               end if
            end if
         end do
      endif

      return
 6000 format(/,1x,'Solute output information, species number ',a15)
 778  format(1x,'Solute output information, species number ',a15)
 6001 format(1x,'Num of solute timesteps ', i6,' Avg tstep = ',
     & 1p,g13.6,' SAI Iter = ',i6,' Tot SAI iter ', i8)
 6010 format(4x,'Node',5x,'an',9x,'anl',7x,'anv',
     *     8x,'mol/s',5x,'residual')
 6012 format(4x,'Node',5x,'an',9x,'anl',7x,'anv',7x,'scl',7x,'scv',
     *     8x,'mol/s',5x,'residual')
 6610 format(4x,'Node',6x,'an',10x,'anl',10x,'anv',
     *     9x,'mol/s',19x,'residual')
 779  format(42x,'src/sink',3x,'equation')
 679  format(49x,'src/sink',5x,'sinkint',6x,'equation')
c6011 format(1x,i5,1x,g10.3,1x,g10.3,1x,g10.3,1x,f9.2,1x,g10.3,1x,g10.3)
 6011 format(1x,i7,1x,g10.3,1x,g10.3,1x,g10.3,1x,g10.3,1x,g10.3)
 6611 format(1x,i7,6(1x,g12.5))
 6013 format(1x,i7, 7(1x,g12.5))
 6020 format(3x,'initial mass =',20x,e14.6,' mol',/,3x,
     &     'current mass = ',19x,e14.6,' mol',/,3x,
     &     'total injected mass = ',12x,e14.6,' mol (',
     &     e13.6,' mol/s)',/,3x,
     &     'total produced mass = ',12x,e14.6,' mol (',
     &     e13.6,' mol/s)',/,3x,'total ',
     &     'mass produced by reaction = ',e14.6,' mol',/,3x,'net mass ',
     &     'balance = ',15x,e14.6,' mol')
 6050 format(1x,'Solute information at time = ',1p,g14.6,' days')
      end

