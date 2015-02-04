      subroutine allocmem
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
CD1 Allocate memory to dynamic variable arrays.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 18-OCT-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/allocmem.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:38   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:25:20   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:10   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:46 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.21   Wed May 29 14:25:30 1996   hend
CD2 Added variable diffusion with water content
CD2 
CD2    Rev 1.20   Mon Mar 04 15:57:26 1996   hend
CD2 Removed uneccessary calculations from coneq1 and added trac input option
CD2 
CD2    Rev 1.19   Wed Feb 07 12:10:36 1996   gaz
CD2 took out 2 comments (no changes)
CD2 
CD2    Rev 1.18   Mon Jan 29 13:10:46 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.17   Mon Jan 29 10:07:56 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.16   Fri Jan 12 17:43:32 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.15   12/13/95 10:25:46   robinson
CD2 Incorporated new zeolite hydration module
CD2 
CD2    Rev 1.14   12/13/95 08:28:32   gaz
CD2 variable dimension for rlp parameters
CD2 
CD2    Rev 1.13   08/02/95 11:48:04   awolf
CD2 No longer allocate a_axy here (done in startup)
CD2 If num of particles is 0, do not allocates those arrays (LLT)
CD2 
CD2    Rev 1.12   06/02/95 08:43:06   llt
CD2 increased allocation for iflx1 and iflx2 to 500
CD2 
CD2    Rev 1.11   06/02/95 08:36:42   llt
CD2 increased allocation for flx12l and flx12v to 500
CD2 
CD2    Rev 1.10   04/27/95 18:25:30   llt
CD2 allocated memory for upwind_l & upwind_v
CD2 
CD2    Rev 1.9   04/10/95 11:43:32   robinson
CD2 No longer allocate array p here, it is done in set_ptrk
CD2 
CD2    Rev 1.8   03/16/95 09:07:46   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.7   03/10/95 10:24:38   llt
CD2 changed allocation of it8-it12 to integers
CD2
CD2    Rev 1.6   01/28/95 13:53:48   llt
CD2 water balance equation was modified
CD2
CD2    Rev 1.5   07/28/94 10:06:56   robinson
CD2 Put allocation of denci in this routine 
CD2 
CD2    Rev 1.4   03/28/94 16:34:50   robinson
CD2 Removed unneeded array.
CD2 
CD2    Rev 1.3   03/23/94 14:42:00   robinson
CD2  
CD2    Rev 1.2   03/18/94 15:53:08   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   03/02/94 13:30:46   llt
CD2 no change
CD2 
CD2    Rev 1.0   01/20/94 10:21:14   pvcs
CD2 original version in process of being certified
CD2 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 None
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
CD4   The definitions for the parameters and dynamic variables are found in
CD4   common files comai.h combi.h, comci.h, comdi.h, comdti.h, comei.h,
CD4   comfi.h comgi.h, comhi.h, comii.h, comrxnb.h, comrxni.h, davidi.h
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   mmgetblk        N/A      Allocate memory to an array
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
CD5   icode           INT      ?Error flag?
CD5   nsize           INT      Length of volf1 and apuv1 arrays
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
CD8 This routine is a general purpose routine necessary with dynamic
CD8 memory allocation.
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 Not Applicable.  See special comments. 
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
CPS BEGIN allocmem
CPS 
CPS   call allocate to allocate memory for each dynamic array variable
CPS   defined in the common files (combi.h, comci.h, comdi.h,
CPS   comei.h, comfi.h, comgi.h, comhi.h, comii.h, comrxnb.h,
CPS   comrxni.h, davidi.h)
CPS 
CPS END allocmem
CPS
C***********************************************************************

      use comflow
      use combi
      use comchem
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comhi
      use comii
      use comji
      use comrxnb
      use comrxni
      use davidi
      use comcouple
      use compart
      use comzeoli
      use comdti
      use comai
      use comsptr
      use comsk
      use comnire
	use trxnvars
c RJP 12/4/06 added following
      use comriv
c RJP 04/17/06 added following
      use comco2

      implicit none

      integer icode
      integer nsize
      integer totcom, totcomalloc, totnum, totvnum
      integer total_max_colloids,total_max_daughters

      if(gdpm_flag.ne.0) then
         allocate(igdpm(neq_primary))
         igdpm = 0
      end if
c*** allocate memory to all arrays in compart ***
      if(ptrak) then 
         allocate(npt(nspeci+1))
         allocate(num_particles(nspeci),ori(nspeci))
c bhl_12/8/05
         allocate(num_part_mem(nspeci))
         allocate(cfraction(nspeci))
         num_part_mem=0
         cfraction=0.0d0
c bhl_12/8/05
         if(nspeci.eq.1) then
            num_particles(1) = max_particles
         end if
         allocate(obj(nspeci),kfact(nspeci))
         allocate(trak_type(nspeci))
         allocate(timeleft(max_particles,nspeci))
         timeleft = 0
         allocate(box(max_particles,nspeci))
         box = 0
         allocate(start_time(max_particles,nspeci))
         start_time = 0
         allocate(frac_done(max_particles,nspeci))
         frac_done = 0
         allocate(theta(max_particles,nspeci))
         theta = 0
         allocate(flow_ot(n0))
         allocate(flow_ottot(neq))
         allocate(Rf(n0,nspeci))
         Rf = 0
         allocate(mass(n0))
         allocate(pconc(n0))
         allocate(kd(maxlayers,nspeci))
         allocate(rd_frac(maxlayers,nspeci))
         allocate(kcoll(maxlayers,nspeci))
         allocate(rcoll(maxlayers,nspeci))
         allocate(fcoll(maxlayers,nspeci))
         allocate(flag_log(nspeci))
         flag_log = 0
         allocate(flag_method(nspeci))
         flag_method = 1
c         allocate(matrix_por(n0))
c         allocate(aperture(n0))
         allocate(matrix_por(maxlayers))
         allocate(aperture(maxlayers))
         allocate(gamma_afm(maxlayers))
         allocate(sresidual(maxlayers))
         allocate(dispflag(maxlayers,nspeci))
         allocate(diffflag(maxlayers,nspeci))
         allocate(sud_frind(maxlayers,nspeci))
         sud_frind = 0
         allocate(num_exit(n0,nspeci))
         num_exit = 0
         if(ncolsizes.ne.0) then
            allocate(sizes(nspeci))
            allocate(probsize(0:ncolsizes,ncolspec))
            allocate(porsize(0:ncolsizes,ncolspec))
            allocate(partsize(max_particles,ncolspec))
            sizes = 0
            probsize = 0.
            porsize = 0.
            partsize = 0.
         end if
cHari 01-Nov-06 colloid diversity model  ..........................
         total_max_colloids = max(1,total_colloids)
         allocate (flag_diversity(nspeci))
         allocate (divs(nspeci))
         allocate (divs_d(nspeci))
         allocate (probdiv(0:maxprobdivs,total_max_colloids))
         allocate (nprobdivs(total_colloids))
         allocate (rcdiv(0:maxprobdivs, total_max_colloids))
c         allocate (rcoll_div(max_particles, total_irrev))
         allocate (flag_col_irrev(nspeci))
         allocate (flag_col_daughter(nspeci))
         allocate (mean_rcoll_reve(total_rev))
         allocate (reves(total_max_colloids))
         allocate (irrevs(total_max_colloids))
         allocate (tprpflag(nspeci))
c s kelkar march 24 07
         allocate (k_rev(nspeci))
         allocate (r_min(nspeci))
         allocate (r_max(nspeci))
         allocate (slope_kf(nspeci))
         allocate (cint_kf(nspeci))
cHari 01-Nov-06 initialize colloid diversity variables
         flag_diversity = .false.
         flag_col_irrev = .false.
         flag_col_daughter = 0
         reves = 0
         irrevs = 0
         tprpflag = 0
         nprobdivs = 0
c s kelkar march 24 07
         divs = 0
         divs_d = 0
         k_rev=0.
         r_min=0.
         r_max=0.
         slope_kf=0.
         cint_kf=0.

      endif
    
c***allocate memory for all arrays in comsptr***groemer 5/1/98***
c*** allocate memory for arrays in comsk*** zvd 08/04/2003
      if(sptrak) then
         allocate(kd(nsize_tprp,1))
         allocate(diffmfl(1,nsize_tprp))
         allocate(rd_frac(nsize_tprp,1))
         allocate(matrix_por(nsize_tprp))
         allocate(aperture(nsize_tprp))
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         allocate(secondary(nsize_tprp))
         secondary = 0.
         allocate(dispersivity1(nsize_tprp))
         allocate(dispersivity2(nsize_tprp))
         allocate(dispersivity3(nsize_tprp))
         allocate(dispersivity4(nsize_tprp))
         allocate(dispersivity6(nsize_tprp))
         allocate(dispersivity7(nsize_tprp))
         allocate(dispersivity8(nsize_tprp))
         allocate(vratio(nsize_tprp))
         allocate(divdwt(nsize_tprp))
         divdwt = 1.
         allocate(gradd(neq,3))
         allocate(tprpflag(nsize_tprp))
         if(nzbtc.gt.0) then
            allocate(zbtc(nzbtc))
            allocate(izonebtc(nzbtc,num_part))
            allocate(totalpart(nzbtc))
            if(div_flag) then
               allocate(totalpart_ret(nzbtc),totlast_ret(nzbtc))
               totalpart_ret = 0
               totlast_ret = 0
            end if
            allocate(totlast(nzbtc,2))
            allocate(ttbtc(nzbtc,num_part))
            izonebtc = 0
            zbtc = 0
            totalpart = 0
            totlast = 0
            ttbtc = 0
         end if
c..... zvd 08/04/2003
         if (.not. omr_flag) then
            omr_nodes = neq
            omr_flag = .true.
         end if
         if (omr_nodes .eq. 0) then
            omr_flag = .false.
         else
            allocate (iboulist(n0,7))
            allocate(node_omr(omr_nodes))
         end if
         allocate (ipc_save(num_part))

c..... skelkar  july23 01..........
         if(nplum.gt.0) then
            allocate(plum(nplum))
            allocate(izoneplum(nplum,num_part))
            allocate(plumpart(nplum))
            izoneplum = 0
            plum = 0
            plumpart = 0
         end if
c..................................

         allocate(omega_partial(n0))
         allocate(sigma_partial(n0))
         allocate(ddx(0:n0))
         allocate(ddy(0:n0))
         allocate(ddz(0:n0))
         allocate(ggg(0:n0,-3:3))
         allocate(corn(0:n0,3))
         allocate(ddxv(num_part))
         allocate(vx1(num_part))
         allocate(axv(num_part))
         allocate(x1(num_part))
         allocate(x2(num_part))
         allocate(x3(num_part))
         allocate(vx1bv(num_part))
         allocate(vx2bv(num_part))
         allocate(dtx(num_part))
         allocate(x_new(num_part))
         allocate(ddyv(num_part))
         allocate(vy1(num_part))
         allocate(ayv(num_part))
         allocate(y1(num_part))
         allocate(y2(num_part))
         allocate(y3(num_part))
         allocate(vy1bv(num_part))
         allocate(vy2bv(num_part))
         allocate(dty(num_part))
         allocate(y_new(num_part))
         allocate(ddzv(num_part))
         allocate(vz1(num_part))
         allocate(azv(num_part))
         allocate(z1(num_part))
         allocate(z2(num_part))
         allocate(z3(num_part))
         allocate(vz1bv(num_part))
         allocate(vz2bv(num_part))
         allocate(dtz(num_part))
         allocate(z_new(num_part))
         allocate(tt(num_part))
         allocate(dt(num_part))
         allocate(x61(num_part))
         allocate(axyzm(num_part))
         allocate(ttt(num_part))
         allocate(ttp1(num_part))
         allocate(ioutt(num_part))
         allocate(istop(num_part))
         allocate(ijkv(num_part))
         allocate(exit_node(num_part))
         allocate(ic_new(num_part))
         allocate(irray(0:n0,-3:3))
         allocate(ijkvs(num_part))
         allocate(divd_omr(n0,3))
         allocate(dpordx_omr(n0,3))
         allocate(xo(num_part))
         allocate(yo(num_part))
         allocate(zo(num_part))
         allocate(ttpo(num_part))
         allocate(lastnode(num_part))
         allocate(pstart_out(num_part))
         allocate(ijkvss(num_part))
         allocate(vx(num_part))
         allocate(vy(num_part))
         allocate(vz(num_part))
         allocate(edtmax(num_part))
         allocate(npt(2))
c... april 12 04 s kelkar 3domr
         allocate(oldnode(num_part))
         allocate(oldnode_rand(num_part))
         allocate(oldnode2(num_part))
         allocate(count_steps(num_part))
! ZVD - 10 August 05, for sptr restarts
         allocate(part_id(num_part,2))
! ZVD - 14 October 05, use ttp1 instead of time_sptr
!         allocate(time_sptr(num_part))

         npt(1) = 0
         npt(2) = n0
cHari 01-Nov-06 colloid diversity model for sptr
         allocate (divs(nspeci))
         allocate (probdiv(0:maxprobdivs,1))
         allocate (nprobdivs(1))
         allocate (rcdiv(0:maxprobdivs, 1))
c s kelkar march 24 07
         allocate (k_rev(nspeci))
         allocate (r_min(nspeci))
         allocate (r_max(nspeci))
         allocate (slope_kf(nspeci))
         allocate (cint_kf(nspeci))
c         allocate (rcoll_div(num_part, 1))
c         allocate (ret_weight(num_part))
cHari 01-Nov-06 initialize colloid diversity variables
         nprobdivs = 0
         divs = 0
         rcdiv = 0.0
         probdiv = 0.0
         k_rev=0.
         r_min=0.
         r_max=0.
         slope_kf=0.
         cint_kf=0.
c         ret_weight = 1./num_part
c..................................................................

      endif 

c***  allocate memory to all arrays in combi ***
      allocate(izonef(n0))
      izonef = 0
      if(interface_flag.ne.0) then
         allocate(izonef_itfc(n0))
         izonef_itfc = 0
      end if
      allocate(zone_pair(max(1,nitfcpairs),2),
     2     red_factor(0:max(1,nitfcpairs+1)))
      allocate(zonec_pair(max(1,ncoldpairs),2),
     2     ftn_factor(0:max(1,ncoldpairs)))
      allocate(filter_flag(nspeci))
      red_factor = 1.
      ftn_factor = 1.
      filter_flag = 0
      zone_pair = 0
      zonec_pair = 0
      if(nitfcitfc.gt.0) then
         allocate(itfcsize(max(1,ncoldpairs)))
         itfcsize = 0
         allocate(itfcporsize(0:nitfcsizes,nitfcitfc))
         allocate(itfcprobsize(0:nitfcsizes,nitfcitfc))
         itfcporsize = 0.
         itfcprobsize = 0.
      end if
      allocate(ka(n0))
      allocate(nar(n0))
      allocate(nelm(nelmd))
      allocate(nelmdg(n0))
c     ***** COMMON Block fbc *****
      allocate(sx1(n0))
c      if(compute_flow) allocate(sxs(nq,6))
c     ***** COMMON Block fbs *****men
      allocate(cord(n0,3))
      allocate(dp(6))
      allocate(dr(8))
      allocate(eta(8))
      allocate(exci(8))
      allocate(si(8))
      allocate(w(8,8))
      allocate(wp(6,6))
      allocate(wr(8,8))
      allocate(wx(8,8))
      allocate(wxp(6,6))
      allocate(wxr(8,8))
      allocate(wy(8,8))
      allocate(wyp(6,6))
      allocate(wyr(8,8))    
      allocate(wz(8,8))
      allocate(wzp(6,6))
      allocate(wzr(8,8))
      allocate(xd(8))
      allocate(xt(8))
      allocate(yd(8))
      allocate(yt(8))
      allocate(zd(8))
      allocate(zt(8))
** allocate memory to arrays in common block comchem
c     ***** COMMON Block place *****
      allocate(pcpnt(ncpnt),pimm(nimm),pvap(nvap))
c     ***** COMMON Block print_flag *****
      allocate(cpntprt(ncpnt),cplxprt(101:ncplx+100))
      allocate(immprt(nimm),vapprt(nvap))
c     ***** COMMON Block chem_name *****
      allocate(cpntnam(ncpnt),cplxnam(101:ncplx+100))
      allocate(immnam(nimm),vapnam(nvap))
      if (iccen .eq. 1) then
         allocate(ndconv(n0))
      else
         allocate(ndconv(1))
      end if
c     ***** COMMON Block rxn_switch
      if(rxn_flag.ne.0)then
c     ***** COMMON Block idntrxn *****
         allocate(naqsp(numrxn),nimsp(numrxn),nivsp(numrxn))
         allocate(irxnic(numrxn,(ncpnt+ncplx)))
         allocate (irxnim(numrxn,nimm),irxniv(numrxn,nvap))
c     ***** COMMON Block rdsp *****
         allocate (ckeq(101:ncplx+100),heq(101:ncplx+100,5))
         allocate (spstoic(101:ncplx+100,ncpnt))
         allocate (cpntgs(ncpnt),idrxn(numrxn),ifxconc(ncpnt))
         allocate(neg_conc_possible(ncpnt))
c     ***** COMMON Block initchem *****
         allocate (cpntsv(ncpnt,n0),calcderiv(ncpnt))
c     ***** COMMON Block couple_index (not allocated dynamically *****
c     ***** COMMON Block rxnpara *****
         allocate(ckeqlb(numrxn),ckmtrn(numrxn),simmmx(numrxn))
c     ***** COMMON Block irrpara *****
         allocate(kfor(numrxn),krev(numrxn))
         allocate(sticirrv(numrxn,(ncpnt+ncplx)))
         allocate(stimirrv(numrxn,nimm))
         allocate(stivirrv(numrxn,nvap))
c     ***** COMMON Block biopara *****
         allocate(ckc(numrxn),cka(numrxn))
         allocate(decay(numrxn),biofac(numrxn),hfac(numrxn))
         allocate(carbfac(numrxn),ammfac(numrxn))
         allocate(phthresh(numrxn),qm(numrxn),yield(numrxn))
         allocate(xminit(numrxn),nbiofrm(numrxn))
         allocate(icbio(numrxn,(ncpnt+ncplx)))
c     ***** COMMON Block precdis *****
         allocate(sarea(numrxn),pdstic(numrxn,ncpnt),pd_flag(nimm,n0))
	   allocate(mw_mineral(nimm), rho_mineral(nimm))
	   allocate(ps_delta_rxn(n0),ps_delta_rxn_s(n0))
c     ***** COMMON Block concs *****
         allocate(totaq(ncpnt),cpnt(ncpnt),cplx(101:ncplx+100))
c     ***** COMMON Block rates *****
         allocate(rrcplx(101:ncplx+100),rrcpnt(ncpnt),rrimm(nimm))
         allocate(rrvap(nvap))
c     ***** COMMON Block deriv_rates *****
         allocate(drdctaq(n7a),drdcimm(nimm,n0),drdcvap(nvap,n0))
         allocate(drcpnt(ncpnt,ncpnt),drimm(nimm),drvap(nvap))
         allocate(drtaq(ncpnt))
         allocate(drcplx(101:ncplx+100))
c     ***** COMMON Block jacobian *****
         allocate(xjac(ncpnt,ncpnt),dxct(ncpnt,ncpnt))
         allocate(resid(ncpnt),ipiv(ncpnt))
c     ***** COMMON Block scaling *****
         allocate(ipiv2(2*ncpnt),nterms(2*ncpnt))
         allocate(sclmtx(2*ncpnt,2*ncpnt),sclfactr(2*ncpnt))
c     ***** COMMON Block rxn_switch ****
         allocate(rxnon(numrxn,n0))
c     ***** COMMON Block temperature ****
         allocate(temp_model(101:ncplx+100))
         allocate(temp_model_kin(numrxn))
         allocate(tcoeff(numrxn,3))
      endif
      if (isalt.ne.0) then
       	   if(.not.allocated(ps_delta_rxn))
     &         allocate(ps_delta_rxn(n0),ps_delta_rxn_s(n0))
        ps_delta_rxn = 0.0d0 
        ps_delta_rxn_s = 0.0d0  
      endif
c***  allocate memory to all arrays in comci ***
c      call storage_derivatives(0)
      allocate(rolf(n0))
      allocate(dil(n0))
      if(irdof.ne.13) then
        allocate(rovf(n0))
        allocate(div(n0))
      else
        allocate(rovf(1))
        allocate(div(1))
      endif	
      allocate(denci(n7))
c***  allocate memory to all arrays in comdi ***
c     ***** COMMON Block fdd *****
      if(irdof.ne.13 .or. ifree .ne. 0) then
        allocate(pcp(n0))
        allocate(dpcef(n0))
        pcp = 0.
        dpcef = 0.
      else
        allocate(pcp(1))
        allocate(dpcef(1))
      end if     
      if(irdof.ne.13) then
        allocate(deneh(n0))
        allocate(denej(n0))
      else
        allocate(deneh(1))
        allocate(denej(1))
      endif	
      allocate(eflow(n0))
      allocate(esk(n0))
      allocate(pflow(n0))
      allocate(phi(n0))
      allocate(pho(n0))
      if (idoff .ne. -1) then
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            allocate(pnx(n3))
            allocate(pny(n3))
            allocate(pnz(n3))
            allocate(pnxi(n3))
            allocate(pnyi(n3))
            allocate(pnzi(n3))
         else
            allocate(pnx(n2))
            allocate(pny(n2))
            allocate(pnz(n2))
            allocate(pnxi(n2))
            allocate(pnyi(n2))
            allocate(pnzi(n2))
         end if
      else
         allocate(pnx(1))
         allocate(pny(1))
         allocate(pnz(1))
         allocate(pnxi(1))
         allocate(pnyi(1))
         allocate(pnzi(1))
      endif
      allocate(ps(n0))
      if (iccen .eq. 1 .or. sptrak) then
         allocate(ps_trac(n0))
      else
         allocate(ps_trac(1))
      end if
      allocate(qh(n0))
      if (irdof .ne. 13 .or. ifree .ne. 0 ) then
         allocate(s(n0))
         allocate(so(n0))
      else
         allocate(s(1))
         allocate(so(1))
      end if
      allocate(sk(n0))
      allocate(t(n0))
      allocate(cpr(n0))
      allocate(denr(n0))
      allocate(denh(n0))
      allocate(denj(n0))
      if(ico2.ge.0.or.ice.ne.0) then
c gaz 10-19-2001
         allocate(thx(n0))
         allocate(thy(n0))
         allocate(thz(n0))
         allocate(qflux(n0))
         allocate(qflxm(n0))
      else
         allocate(thx(1))
         allocate(thy(1))
         allocate(thz(1))
         allocate(qflux(1))
         allocate(qflxm(1))
      endif
      allocate(dvas(n0))
      allocate(to(n0))
      allocate(vf(n0))
      allocate(volume(n0))
      allocate(wellim(n0))
c     ***** COMMON Block fddi *****
      allocate(nskw(n0))
      allocate(nskw2(n0))
c     ***** COMMON Block fdd1 *****
      totcom = nvap+nimm+ncpnt
      totcomalloc = max(totcom,1)
      totnum = max (numd, numsorp, 1)
      allocate(a1adfl(totcomalloc,totnum),a1adfv(totcomalloc,totnum),
     2     a2adfl(totcomalloc,totnum))
      allocate(a2adfv(totcomalloc,totnum),an(n7))
      allocate(anl(n7),anlo(n7),anv(n7),betadfl(totcomalloc,totnum))
      allocate(betadfv(totcomalloc,totnum),cm(totcomalloc),
     2     cm0(totcomalloc),cnsk(n7),pcnsk(n7))
      allocate(henry_model(totcomalloc),hawwa(totcomalloc,5))
      allocate(a_henry(totcomalloc),dh_henry(totcomalloc))
      allocate(conc_read(totcomalloc))
      henry_model = 0
      pcnsk = 0
      conc_read = .false.
      allocate(cp1f(nrlp),cp2f(nrlp),cp3f(nrlp),cp4f(nrlp))

      if(ptrak) then
c         allocate(diffmfl(nspeci,n0))
c         allocate(tclx(1,n0),tcly(1,n0),tclz(1,n0))
         allocate(diffmfl(nspeci,maxlayers))
         allocate(tclx(1,maxlayers),tcly(1,maxlayers),tclz(1,maxlayers))
      elseif(.not.sptrak) then
         allocate(diffmfl(totcomalloc,totnum))
         allocate(tclx(totcomalloc,totnum),tcly(totcomalloc,totnum),
     2        tclz(totcomalloc,totnum))
      end if

      allocate(dench(n7),dencj(n7),diffmfv(totcomalloc,totnum))
c     allocate(dit(300),fc(n7),flx12l(500),flx12v(500))
c     allocate(fc(n7),flx12l(500),flx12v(500))
      allocate(flx12l(nflxt),flx12v(nflxt))
      allocate(iflx1(nflxt),iflx2(nflxt))
      if(.not.allocated(nflxc)) then
        allocate(nflxc(iflxc))
        iflxc = 0
      endif
      allocate(fc(n7))
      allocate(qcin(totcomalloc),qcout(totcomalloc),qcrxn(totcomalloc),
     2     rc(n7))
      allocate(rcss(n7),rp1f(nrlp),rp2f(nrlp),rp3f(nrlp))
      allocate(rp4f(nrlp),rp5f(nrlp),rp6f(nrlp),rp7f(nrlp))
      allocate(rp8f(nrlp),rp9f(nrlp),rp10f(nrlp),rp11f(nrlp))
      allocate(rp12f(nrlp),rp13f(nrlp),rp14f(nrlp),rp15f(nrlp))
      allocate(rp16f(nrlp),rp17f(nrlp),rp18f(nrlp),rp19f(nrlp))
      allocate(rp20f(nrlp),rp21f(nrlp),rp22f(nrlp),rp23f(nrlp))
      allocate(rp24f(nrlp))
      allocate(t1sk(n7),t2sk(n7))
      allocate(tcvx(totcomalloc,totnum),tcvy(totcomalloc,totnum),
     2     tcvz(totcomalloc,totnum))
      allocate(sehdiff(totnum),sehdiffv(totnum))
! Conca diffusion model no longer requires thetaresid
!      allocate(thetaresid(totcomalloc,totnum),thetaresidv(totcomalloc,totnum)) 
      totvnum = max (numvcon, 1)
      allocate(vc1f(totvnum),vc2f(totvnum),vc3f(totvnum),vc4f(totvnum),
     2    vc5f(totvnum),vc6f(totvnum),vc7f(totvnum),vc8f(totvnum))
c     ***** COMMON Block fdd1i *****
      allocate(iadd(totcomalloc),iaddt(totcomalloc),
     2     iadsfl(totcomalloc,totnum), iadsfv(totcomalloc,totnum))
      allocate(icapt(nrlp),icns(totcomalloc))
c     allocate(irlpt(nrlp),itc(100),ivcn(n0),ivcon(totvnum))
      allocate(irlpt(nrlp),ivcn(n0),ivcon(totvnum))
      allocate (mflagl(totcomalloc,totnum),mflagv(totcomalloc,totnum))
c      if(.not.allocated(npt)) allocate(npt(max(1,totcomalloc)))
      if(.not.allocated(npt)) allocate(npt(totcomalloc))
c     ***** COMMON Block fdd2 *****
      allocate(phini(n8),psini(n8),tini(n8))
      if(.not.compute_flow) then
         n8 = 1
         n6 = 1
      end if
      if(irlp_fac.eq.1) then
         allocate(rlp_fac(2*n0))
      endif                     
****   TENMA   ****
      if(abs(iporos).eq.1) then
         allocate(amgang(n8),dporp(n8),dport(n8))
      else if(abs(iporos).eq.2) then
         allocate(amgang(n8),dporp(n8),dport(n8),pgangi(n8))
         allocate(strgan(n8),alphae(n8))
      else if(abs(iporos).eq.4) then
         allocate(nskw3(n8)) 
         allocate(amgang(n8),dporp(n8),dport(n8),pgangi(n8),
     &        wgangi(n8),sgangi(n8),agangi(n8),tenma_ww(n8))
      else if(iporos.eq.5.or.iporos.eq.7) then
c        parameters for temperature-dependent porosity
         allocate(dporp(n8),dport(n8))
         allocate(porTemp1(n8),porTemp2(n8),porTemp3(n8),porTemp4(n8))
      else if(iporos.eq.-5) then
         allocate(amgang(n8),dporp(n8),dport(n8),wgangi(n8))
      else if(iporos.eq.6.or.iporos.eq.7) then
      	 allocate(amgang(n8),dporp(n0),dport(n0))
      else
         allocate(dporp(n0),dport(n0))
      endif
****   TENMA   ****
c     ***** COMMON Block fice *****
c  RJP added icarb flag below
      if ((ice .eq. 1) .or. (icarb.eq.1)) then
         allocate(sii(n6),sio(n6))
c     ***** COMMON Block iice *****
         allocate(ices(n6))
      else
         allocate(sii(1),sio(1))
         allocate(ices(1))
      end if
c     ***** COMMON Block fddi1 *****
      allocate(ieos(n0),ieos_ch(n0),iieos(n0),iporf(n0))
      if (rlp_flag .eq. 1 .or. nrlp .ne. 0) then
         allocate(icap(n0), irlp(n0))
      else
         allocate(icap(1), irlp(1))
      end if
      allocate(itrc(n7),itrcdsp(n7))
c***  allocate memory to all arrays in comei ***
c     ***** COMMON Block fee ***** 
c     ***** COMMON Block fff ***** 
c     ***** COMMON Block fhh ***** 
      allocate(iirb(n0),irb(n0),nopt(n0),npvt(n0))
c     ***** COMMON Block fii ***** 
      allocate(c(maxor),g(maxor),h(n0,maxor),ss(maxor),y(maxor)) 
c***  allocate memory to all arrays in comfi ***
c     ***** COMMON Block co2 ***** 
      allocate (qc(n4))
      if (irdof .ne. 13 .or. ifree .ne. 0) then
         allocate(pci(n4),pcio(n4))
         if(.not.compute_flow) then
            n4 = 1
         end if
      else
         allocate(pci(1),pcio(1))
         n4 = 1
      end if
      allocate(cnlf(n4),cnvf(n4),dcc(n4),dce(n4))
      allocate(dclcf(n4),dclef(n4),dclf(n4),dcp(n4))
      allocate(dcqc(n4),dcqh(n4),dcvcf(n4),dcvef(n4))
      allocate(dcvf(n4),dec(n4),delcf(n4),denpch(n4))
      allocate(denpci(n4),denpcj(n4),deqc(n4),devcf(n4))
      allocate(dglc(n4),dgvc(n4),dilc(n4),divc(n4))
      allocate(dmc(n4),dqc(n4),dtpac(n4),dqpc(n4),qng(n4))
      allocate(eskc(n4))
      allocate(sici(n4))
c***  allocate memory to all arrays in comgi ***
c     ***** COMMON Block fgg ***** 
      allocate(dfee(nn),dfep(nn),dfepc(nn),dfme(nn))
      allocate(dfmp(nn),dfmpc(nn),dfpce(nn),dfpcp(nn))
      allocate(dfpcpc(nn))
c***  allocate memory to all arrays in comhi ***
c***  ***** COMMON Block dualp ***** 
      if(.not.compute_flow) then
         n5 = 1
      end if
      allocate(a21eef(n5),a21epf(n5),a21mef(n5),a21mpf(n5))
      allocate(a32eef(n5),a32epf(n5),a32mef(n5),a32mpf(n5))
      allocate(r3ef(n5),r3mf(n5),rb2ef(n5),rb2mf(n5))
      allocate(tb11(n5),tb12(n5),tb21(n5),tb22(n5))
      allocate(volf2(n5),wb11(n5),wb12(n5),wb21(n5))
      allocate(wb22(n5))
      if( idpdp .ne. 0 ) then
         nsize = neq
      else
         nsize = n5
      end if
      allocate(apuv1(nsize),volf1(nsize)) 
c***  allocate memory to all arrays in comii ***
c     ***** COMMON Block coeff ***** 
      allocate(cel(20,4),cev(20,4),crl(20,4),crv(20,4),cvl(20,4))
      allocate(cvv(20,4))
c     ***** COMMON Block coeff1 ***** 
      allocate(pmax(3),pmin(3),tmax(3),tmin(3)) 
c***  allocate memory to all arrays in comrxnb ***
c     ***** COMMON Block comrxnb ***** 
      allocate(avgmolwt(n0))
      if (iccen .eq. 1) then
         allocate(akc(n7),danl(n7),danv(n7))
         allocate(dsccl(n7),dsccv(n7),scl(n7),scv(n7))
      else
         allocate(akc(1),danl(1),danv(1))
         allocate(dsccl(1),dsccv(1),scl(1),scv(1))
      end if
      if(compute_flow .or. iccen .eq. 1) then
c      allocate(stored_derivative(n7),stored_residual(n7))
c***  allocate memory to all arrays in comrxni ***
c     ***** COMMON Block stor ***** 
c     ***** COMMON Block dispersion ***** 
         if (iccen .eq. 1) then
            allocate(displx(n0),disply(n0),displz(n0))
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               allocate(dispvx(n0),dispvy(n0),dispvz(n0))
            end if
         end if
      end if
c***  allocate memory to all arrays in davidi ***
c     ***** COMMON Block david2 ***** 
c  gaz 2-25-03 needed idof+1, idof*idof +1 for ice
      allocate(nb(37),nmat(37),nmatb(36),nrhs(7),nrhsb(7))
      allocate(nb_sol(mdof_sol**2),nmat_sol(mdof_sol**2))
      allocate(nrhs_sol(mdof_sol))
      allocate(nvar_select(max_variables))
c     open t1-t13 space
      allocate(t1(nn),t2(nn),t3(nn),t4(nn))
      allocate(t5(nn),t5v(nn),t6(nn),t7(nn))
c RJP 04/08/07 added t13 for diffusion terms.
      allocate(t8(nn),t9(nn),t10(nn),t13(nn),t14(nn),t17(nn))
      allocate(it4(nn),it5(nn),it6(nn),it8(nn),it8a(nn))
      allocate(it9(nn),it10(nn),it11(nn),it12(nn),it13(nn))
      allocate(it9a(nn),it10a(nn),it11a(nn),it12a(nn),it13a(nn))
c     allocate memotry to upwind arrays ***
      if(compute_flow .or. iccen .eq. 1) then
         allocate(upwind_l(4*n0))
         if (irdof .ne. 13 .or. ifree .ne. 0) then
            allocate(upwind_v(4*n0))
         else
            allocate(upwind_v(1))
         end if
      else
         allocate(upwind_l(1),upwind_v(1))
      end if
      if( izeolites .eq. 1 ) then
         allocate(kzeol(n0),fwater(n0),fwater_old(n0))
      end if

c---------- PHS 4/13/05 added sink_integ to integrate the source 
c                 sink at a single node over a run in comdi
      allocate(sinkint(n7))
      if (n7 .ne. 0) sinkint = 0

      end
