      subroutine  cnswer(tol_value, sia_iter)
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
CD1 To compute the new solute concentrations of the current chemical
CD1 at the current time.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 ?            G. Zyvoloski   N/A     Initial implementation -
CD2                                     however, previous non-YMP
CD2                                     versions of FEHM exist, and the
CD2                                     current version may differ from
CD2                                     these
CD2 04-06-93      B. Robinson   22      Revised to implement the
CD2                                     chemical reaction models
CD2                                     between solutes
CD2                                     
CD2 07-15-96      H. Viswanathan N/A    Revised for new chemical model
CD2
CD2  $Log:   /pvcs.config/fehm90/src/cnswer.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:10   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:25:54   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:20   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Fri May 24 09:51:44 1996   hend
CD2 Updated trac for mdnodes
CD2 
CD2    Rev 1.8   Thu Mar 21 13:17:18 1996   hend
CD2 Fixed for new coneq1 calling sequence
CD2 
CD2    Rev 1.7   Mon Jan 29 13:30:52 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   08/16/95 16:24:24   robinson
CD2 Corrected iteration flag for the -maxit option
CD2 
CD2    Rev 1.5   08/07/95 11:04:46   awolf
CD2 Fixed for dpdp problems. neq changed to n0 for
CD2 matrix setup
CD2 
CD2    Rev 1.4   01/28/95 14:19:56   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.3   03/28/94 16:35:46   robinson
CD2 Removed unneeded array.
CD2
CD2   Rev 1.3   11/08/93 17:51:12   llt
CD2 Changed include files to .h files.
CD2
CD2   Rev 1.2   11/05/93 16:34:08   llt
CD2 Added $Log$ keyword
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type     Use   Description
CD3 
CD3 tol_value    int      I/O   Flag indicating the status of the
CD3                             solution for the current solute
CD3 iter_counter int      I/O   Current iteration number
CD3              int      I     Flag denoting if full solution or
CD3                             nspeci x nspeci solution is being
CD3                             performed
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
CD4 Global Constants
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
CD4 nsp, iadd, maxit, iout, fdum, iatty, idpdp, neq, an, bp, npn, epc,
CD4 tmch, npt, gamma
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier      Type     Description
CD4 
CD4 react           N/A      Compute reaction terms and add to 
CD4                          Jacobian and residual for the main
CD4                          or fracture nodes
CD4 thermc          N/A      Computes solute source/sink and reaction
CD4                          terms
CD4 coneq1          N/A      Computes advection and dispersion terms
CD4                             for a given node
CD4 dualta          N/A      Computes terms specific to dual porosity
CD4                             concentration values
CD4 gencon          N/A      Computes new residuals, calculates new
CD4                          concentration values
CD4 gentdp          N/A      Same as gencon for dual porosity/dual
CD4                          permeability simulation
CD4 zeror_out       N/A      Sets all array values to zero (real array)
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
CD5 id           int         Do-loop index parameter
CD5 nsp          int         Species number for variable
CD5                          we are computing rate contribution
CD5                          for (need better description)
CD5 spec_num     int         Current species number
CD5 deriv_spec_num int       Species number for variable we are
CD5                          taking the concentration with respect
CD5                          to (need better description)
CD5 icouple      int         Species number for variable
CD5                          we are computing rate contribution
CD5                          for (need better description)
CD5 jcouple      int         Species number for variable we are
CD5                          taking the concentration with respect
CD5                          to (need better description)
CD5 igrp       int         Do loop index indicating the current
CD5                          group number
CD5 matnum       int         Index denoting the current submatrix 
CD5                          number
CD5 neqp1        int         Number of equations plus 1
CD5 nsizea       int         Size of A matrix
CD5 c_node       int         The position of the current species
CD5                          in the concentration array
CD5 mi           int         The postion of current species in
CD5                          concentration array
CD5 icnls        int         Temporary storage to save icnl variable
CD5 irxn         int         Current reaction number
CD5                        
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 The following is a description of the functions are performed by
CD6 this routine:
CD6 
CD6   The concentrations for the current chemical are computed
CD6   iteratively until convergence is achieved in an outer LOOP-
CD6   ENDLOOP.
CD6   
CD6 Inside this LOOP:  
CD6 
CD6   The number of iterations is checked to see if the maximum has
CD6   been reached, and if it has, warning messages are written.
CD6   
CD6   thermc is called to compute the source/sink and chemical
CD6   reaction terms of the mass balance.
CD6
CD6   coneq1 is called to calculate the advection and 
CD6   dispersion terms for each node for each species in 
CD6   the current grouping
CD6
CD6   Within another loop over all nodes, the code computes the
CD6   chemical reaction terms by:
CD6   
CD6     Checking the for reaction and computing the reaction rate term
CD6     based on the model chosen.  If the multiple reaction model was
CD6     selected, a call to mult_rxn is used.
CD6     
CD6   
CD6   The convergence criterion is initialized to 0.
CD6   
CD6   An IF block is used to determine if a conventional or dual
CD6   porosity/dual permeability (dpdp) solution is being performed,
CD6   and calls either gencon or gentdp to obtain the new solution at
CD6   the current iteration.
CD6   
CD6   The loop is exited if the tolerance parameter returned from
CD6   gencon indicates that new concentrations were not calculated
CD6   this iteration.
CD6   
CD6   If an iteration was taken, the counter is incremented and new
CD6   concentrations are computed.
CD6   
CD6   The loop is exited if we achieved convergence, after setting the
CD6   value of the tolerance parameter to indicate convergence.
CD6   
CD6
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
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
CD9 2.3.4 Solute-transport equations
CD9 2.5.2 Solve nonlinear equation set at each time step
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
CPS BEGIN cnswer
CPS 
CPS Set number of SIA iterations performed to 0
CPS 
CPS LOOP for solving for total aqueous concentration at 
CPS ... each position for current component
CPS 
CPS   IF the maximum number of iterations has been reached
CPS     Write warning message to output file
CPS     
CPS     IF another output file for this information exists
CPS       Write warning message to this output file
CPS     ENDIF
CPS     Set parameter to indicate no convergence
CPS     ERROREXIT
CPS   ENDIF
c changed  9/20/95 to add component loop
CPS
CPS 
CPS   
CPS   FOR each unknown in a total aqueous component
CPS     Set residual to 0
CPS   ENDFOR
CPS 
CPS   Compute number of elements in A matrix array
CPS   FOR each element in a matrix array for the current aqueous component
CPS     Set element in A matrix array to 0
CPS   ENDFOR
CPS 
CPS 
CPS   thermc - compute solute mass terms for this specie
CPS   FOR each unknown in the total aqueous component
CPS     call coneq1 - compute advection and dispersion terms for ... 
CPS     ... this unknown for the main or fracture nodes
CPS   ENDFOR each unknown
CPS   call react - compute reaction terms and add to Jacobian ...
CPS   ... and residual for the main or fracture nodes
CPS   IF a dual porosity simulation is performed
CPS     thermc - compute solute mass and reaction terms ... 
CPS     ... for this specie for the 1st set of matrix nodes
CPS     react - compute reaction terms and add to Jacobian ...
CPS     ... and residual for the 1st set of matrix nodes
CPS     thermc - compute solute mass and reaction terms ... 
CPS     ... for this specie for the 2nd set of matrix nodes
CPS     react - compute reaction terms and add to Jacobian ...
CPS     ... and residual for the 2nd set of matrix nodes
CPS         dualta - compute dual porosity terms
CPS   ENDIF
CPS   
CPS   Set convergence criterion to 0
CPS   
CPS   IF dpdp solution is not being performed
CPS     gencon - compute new residuals, solve for new ...
CPS     ... concentrations for each grouping
CPS   ELSE a dpdp solution is being performed
CPS     gentdp - compute new residuals, solve for new ...
CPS     ... concentrations for each grouping
CPS   ENDIF
CPS   FOR each node
CPS     Compute new total aqueous concentration value
CPS   ENDFOR each node
CPS ENDFOR endloop over an total aqueous species component
CPS   
CPS EXITIF the solution is converged over all total aqueous species
CPS   
CPS Add one to the number of SIA iterations 
CPS   
CPS Perform another SIA iteration by looping over all components again   
CPS   
CPS 
CPS ENDLOOP
CPS 
CPS END cnswer
CPS 
C**********************************************************************

      use comai
      use combi
      use comchem
      use comci
      use comcouple
      use comdi
      use comdti
      use comei
      use comfi
      use comgi
      use comrxni
      implicit none

      integer neqp1
      integer nsizea
      integer mi
      integer icnls
      integer sia_iter
      integer isolute
      integer tol_value
      integer iter_flag
      real*8 dc_sia
      real*8 sia_old_imm
      real*8 oldimm
      real*8 rc_molkghr
      real*8 newimm
      real*8 h_const
      real*8 dvap_conc
      integer igrp
      integer aq_compnum
      integer d_aq_compnum
      integer ic,im,in,iv
      integer matnum, spec_num, ic2
      integer matnum_dp,spec_num_dp
      integer itemp
      integer c_node
      integer num_left_to_solve
      integer reset_flag

      save reset_flag

      sia_iter = 0
      c_node = 0
      do isolute = 1, nspeci
         do in = 1,n0
            mi=in+(isolute-1)*n0
            anlo(mi)=an(mi)
         enddo
      enddo
      ndconv=0
      tol_value = 0
c     start SIA loop
 1000 continue
      iter_flag = 0
      if ( sia_iter .gt. abs(iaccmx) )   then
         tol_value = 1
         reset_flag = 1
         if (iout .ne. 0) write(iout, 6000) 
         if ( iptty .ne. 0 ) write(iptty ,6000)  
         goto 9000
      end if
      neqp1=neq+1
      if(rxn_flag.eq.1)then
         if(iskip.eq.0.or.reset_flag.eq.1)then
            ndconv=0
         else
            num_left_to_solve = 0
            do in = 1,n0
               if(ndconv(in).lt.(ncpnt+nimm+nvap))then
                  num_left_to_solve = num_left_to_solve+1
                  ndconv(in)=0
               endif
            enddo
            if(iptty.ne.0)
     2           write(iptty,*)'# nodes chemistry on ',num_left_to_solve
         endif   

         call chemod(dtotc/3600,tol_value)

         if (tol_value .eq. 2) then
            if (iout .ne. 0) write(iout, 6020) 
            if ( iptty .ne. 0 ) write(iptty, 6020)  
            goto 9000
         elseif (tol_value .eq. 3) then
            goto 9000 
         endif
      else
c zvd 08-Oct-08 : check to ensure all species have converged at node
         do in = 1,n0
            if(ndconv(in).lt.(ncpnt+nimm+nvap))then
               ndconv(in) = 0
            end if
         end do
      endif

      do igrp = 1, ngroups
         matnum = 1
         spec_num = 1
         matnum_dp = 1
         spec_num_dp = 1
         bp = 0 
         nsizea=(nelm(neqp1)-neqp1)*n_couple_species(igrp)**2
c     if(idpdp.ne.0)nsizea=4*nsizea
         a=0
         do aq_compnum = 1, n_couple_species(igrp)
            ic = pos(igrp,aq_compnum)
            nsp = pcpnt(ic)
            npn = npt(nsp)
            do d_aq_compnum = 1, n_couple_species(igrp)
               ic2 = pos(igrp,d_aq_compnum)
               if(ic.eq.ic2)then
                  call  thermc  (0)
                  if(idpdp.eq.0) then
c     gaz 2-13-03        call coneq1(0, matnum,spec_num)
                     if(ianpe.ne.0) then
                        call coneq1_ani(0, matnum,spec_num)
                     else if (gdpm_flag.ne.0) then
                        call coneq1_gdpm(0, matnum,spec_num)
                     else
                        call coneq1(0, matnum,spec_num)
                     endif
                     if (imdnode.ne.0) 
     2                    call coneq1mdnode(matnum,spec_num)
                  else
                     call gentdp(0,igrp,matnum_dp,spec_num_dp
     2                    ,sia_iter,ic,n_couple_species(igrp))
                  endif
                  if(rxn_flag.eq.1)then
                     if(idpdp.eq.0)then
                        call react(1, matnum, spec_num, 
     2                       ic,ic2,igrp)
                     else
                        call react(1, matnum_dp, spec_num_dp, 
     2                       ic,ic2,igrp)
                        call react(2,matnum_dp,spec_num_dp,
     2                       ic,ic2,igrp)
                     endif    
                  endif
                  if ( idualp .ne. 0) then
                     call thermc(neq)
                     if(rxn_flag.eq.1)call react(2,matnum,spec_num,
     2                    ic,ic2,igrp)
                     call thermc(neq+neq)
                     if(rxn_flag.eq.1)call react(3,matnum,spec_num,
     2                    ic,ic2,igrp)
c     icnls = icnl
c     icnl = 1
                     call  dualta(spec_num,matnum)
c     icnl = icnls
                  end if     
               else
                  if(rxn_flag.eq.1)then
                     if(idpdp.eq.0)then
                        call react(1, matnum, spec_num,
     2                       ic,ic2,igrp)
                     else
                        call react(1, matnum_dp, spec_num_dp,
     2                       ic,ic2,igrp)
                        call react(2, matnum_dp, spec_num_dp,
     2                       ic,ic2,igrp)
                     endif  
                  endif
               endif
               matnum=matnum+1
               matnum_dp=matnum_dp+2
            enddo
            spec_num = spec_num + 1
            matnum_dp=matnum_dp+n_couple_species(igrp)*2
            spec_num_dp = spec_num_dp+2
         enddo
         if(idpdp.eq.0) then
            call  gencon(igrp,sia_iter,n_couple_species(igrp))
         else
            call gentdp(1,igrp,matnum,aq_compnum,sia_iter,1,
     2           n_couple_species(igrp))
         endif

         do itemp=1,n_couple_species(igrp)
            ic = pos(igrp,itemp)
            do in=1,n0
               if(ndconv(in).ne.(ncpnt+nimm+nvap))then
                  mi = in+(itemp-1)*n0
                  c_node = in + (pcpnt(ic)-1)*n0
                  an(c_node) =  an(c_node) - bp(mi)
                  if(icns(nsp) .eq. -2 ) then
                     if(henry_model(nsp).eq.1)then
                        h_const= a_henry(nsp)*
     2                       exp(dh_henry(nsp)/
     3                       (gas_const*0.001)*(1/298.16-1/
     4                       (t(in)+temp_conv)))
                     elseif(henry_model(nsp).eq.2)then
                        h_const = 10**(hawwa(nsp,1)+
     2                       hawwa(nsp,2)*(t(in)+
     3                       temp_conv)+hawwa(nsp,3)
     4                       /(t(in)+temp_conv)+hawwa(nsp
     5                       ,4)*dlog10(t(in)+temp_conv)+hawwa
     6                       (nsp,5)/(t(in)+temp_conv)**2)
                        h_const= (101325*rolf(in)*1e-3)/
     2                       (h_const*1e6*mw_water)
                     elseif(henry_model(nsp).eq.3)then
                        h_const= (phi(in)-pci(in))*dh_henry(nsp)
                     endif
                     dvap_conc = ( mw_water * h_const ) /
     2                    ( phi(in) * avgmolwt(in) )
                     anl(c_node)=an(c_node)
                     anv(c_node)= dvap_conc* an(c_node)
                  endif
                  if (an(c_node) < 0.d0 .and.
     2                 neg_conc_flag.eq.0) then
                     if(rxn_flag.ne.0) then
                        if (.not.(ifxconc(ic).ge.2 .or.
     2                       neg_conc_possible(ic).eq.1)) then
                           tol_value = 2
                           if (iout .ne. 0) write(iout, 6010) in, 
     &                          cpntnam(ic)
                           if (iptty .ne. 0) write(iptty, 6010) in, 
     &                          cpntnam(ic)
                           return
                        endif
                     else
                        tol_value = 2
                        if (iout .ne. 0) write(iout, 6010) in, 
     &                       cpntnam(ic)
                        if (iptty .ne. 0) write(iptty, 6010) in, 
     &                       cpntnam(ic)
                        return
                     end if
                  endif

                  if(abs(bp(mi)).gt.max(epc*an(mi),epc))then
                     iter_flag=1
                  elseif(abs(an(c_node)).gt.epc*epc)then
                     if(abs(bp(mi)/an(c_node)).gt.0.01)then
                        iter_flag=1
                     else
                        ndconv(in)=ndconv(in)+1
                     endif
                  else
                     ndconv(in)=ndconv(in)+1
                  endif
               endif
            enddo
         enddo
      enddo

      do iv=1,nvap
         nsp = pvap(iv)
         npn = npt(nsp)
         bp=0
         nsizea=(nelm(neqp1)-neqp1)
         a=0
         call thermc (0)
         if(idpdp.eq.0) then
c     gaz 2-13-03 call coneq1(0,1,1)
            if(ianpe.ne.0) then
               call coneq1_ani(0,1,1)
            else if (gdpm_flag.ne.0) then
               call coneq1_gdpm(0,1,1)
            else
               call coneq1(0,1,1)
            endif
            if (imdnode.ne.0)
     2           call coneq1mdnode(1,1)
         else
            call gentdp(0,1,1,1,sia_iter,1,1)
         endif
         if(rxn_flag.eq.1)then
            if(idpdp.eq.0)then
               call react(1,1,1,1,iv,1)
            else
               call react(1,1,1,1,iv,1)
               call react(2,1,1,1,iv,1)
            endif
         endif
         if ( idualp .ne. 0) then
            call thermc(neq)
            if(rxn_flag.eq.1)call react(2,1,1,
     2           1,iv,1)
            call thermc(neq+neq)
            if(rxn_flag.eq.1)call react(3,1,1,
     2           1,iv,1)
c     icnls = icnl
c     icnl = 1
            call  dualta(spec_num,matnum)
c     icnl = icnls
         end if     
         if(idpdp.eq.0) then
            call  gencon(1,sia_iter,1)
         else
            call gentdp(1,1,1,1,sia_iter,1,1)
         endif
         
         do in=1,n0
            if(ndconv(in).ne.(ncpnt+nimm+nvap))then
               mi = in+(pvap(iv)-1)*n0
               an(mi) =  an(mi) - bp(in)
               if(abs(bp(in)).gt.max(epc*an(mi),epc))then
                  iter_flag=1
               elseif(abs(an(mi)).gt.epc*epc)then
                  if(abs(bp(in)/an(mi)).gt.0.01)then
                     iter_flag=1
                  else
                     ndconv(in)=ndconv(in)+1
                  endif
               else
                  ndconv(in)=ndconv(in)+1
               endif
            endif
         enddo
      enddo

      do im = 1, nimm
         nsp = pimm(im)
         call thermc(0)
         do in=1,n0
            if(ndconv(in).ne.(ncpnt+nimm+nvap))then
               if(idpdp.ne.0)call thermc(neq)
               mi = in + (pimm(im)-1)*n0
               sia_old_imm = an(mi)
               oldimm = anlo(mi)
               rc(mi)=-rc(mi)
               rc_molkghr = -(rc(mi)*3600)/(denr(in)*sx1(in))
               an(mi)=((dtotc/3600)*(rc_molkghr
     2              -drdcimm(im,in)*sia_old_imm)+oldimm)/
     3              (1-(dtotc/3600)*drdcimm(im,in))
               dc_sia = abs(an(mi)-sia_old_imm) 
               if(abs(dc_sia).gt.max(epc*an(mi),epc))then
                  iter_flag = 1
               elseif(abs(an(mi)).gt.epc*epc)then
                  if(dc_sia/an(mi).gt.0.01)then
                     iter_flag=1
                  else
                     ndconv(in)=ndconv(in)+1
                  endif
               else
                  ndconv(in)=ndconv(in)+1
               endif
            endif
         enddo
      enddo
c     if( iter_flag.eq.0.and.sia_iter.gt.1 ) then
      if( iter_flag.eq.0.and.sia_iter.gt.0 ) then
         sia_iter = sia_iter+1
         if(iptty.ne.0)write(iptty,*)'Sia Iteration Number ',sia_iter
         goto 3000
      end if
      sia_iter = sia_iter+1
      if(iptty.ne.0)write(iptty,*)'Sia Iteration Number ',sia_iter
      go  to  1000
 3000 continue
      reset_flag = 0
      do im = 1,nimm
         do in = 1,n0
            mi = in+(pimm(im)-1)*n0
            newimm = an(mi)
            if(newimm.le.conc_min)then
               an(mi)=0.d0
c     No mineral present at node
               pd_flag(im,in)=1
c     Added 12/30/99 by geh
            else
               pd_flag(im,in)=0
c     end add - geh
            endif
         enddo
      enddo
 9000 continue     
      return
 6000 format(/,1x,'SIA iterations gt maxit  ')
 6010 format(/,1x,'Negative concentrations encountered at node: ', 
     2     i7, ', Component: ', a12, '; Cut time step.')
 6020 format(/,1x,'Negative concentrations encountered in CHEMOD', 
     2     ', cut time step.')
 6030 format(/,1x,'Negative mineral concentrations encountered at',
     2     '  node: ', i7, ', Mineral: ', i2)
 6040 format(/,1x,'Negative h+ concentrations encountered at',
     2     '  node: ', i7)
      end
