      subroutine  data
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
CD1 Initialize scalar variables, zero all arrays, and load thermodynamic
CD1 coefficients.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 06-DEC-93    Z. Dash        22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/data.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:01:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:52   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:02   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:10   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:24 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.19   Tue Jun 04 09:49:08 1996   hend
CD2 Fixed size of sehdiff and sehdiffv
CD2 
CD2    Rev 1.18   Wed May 29 14:35:06 1996   hend
CD2 Added variable diffusion with water content
CD2 
CD2    Rev 1.17   Thu Feb 15 09:37:52 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.16   Wed Feb 07 10:53:12 1996   gaz
CD2 changed iad_up to 1000 added icons=1000
CD2 
CD2    Rev 1.15   Mon Jan 29 14:16:54 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.14   Fri Jan 12 16:49:46 1996   gaz
CD2 initialized irlp to 1 ,not 0
CD2 
CD2    Rev 1.13   Tue Jan 09 14:01:54 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.12   12/13/95 08:31:50   gaz
CD2 changed nbits to 256 and upwind_l(and v) to 4*n0
CD2 
CD2    Rev 1.11   08/16/95 16:27:44   robinson
CD2 Changed name of variable to set print out interval
CD2 
CD2    Rev 1.10   06/02/95 08:44:42   llt
CD2 increased initialization size of flux variables to 500
CD2 
CD2    Rev 1.9   04/27/95 18:24:54   llt
CD2 initialized iad_up & nbits
CD2 
CD2    Rev 1.8   04/25/95 09:23:06   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.7   03/24/95 09:41:48   llt
CD2 initialized iupk to 0
CD2 
CD2    Rev 1.6   01/28/95 13:54:10   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.5   03/28/94 16:37:32   robinson
CD2 Removed unneeded array.
CD2 
CD2    Rev 1.4   03/23/94 14:39:56   robinson
CD2 Additional cleanup of memory management
CD2 
CD2    Rev 1.3   03/18/94 15:42:58   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.2   03/09/94 10:44:28   llt
CD2 corrected initialization of arrays in common block fdd1i (from 
CD2 zeror_out to zeroi_out)
CD2
CD2    Rev 1.1   03/08/94 15:41:52   llt
CD2 no changes
CD2 
CD2    Rev 1.0   01/20/94 10:22:46   pvcs
CD2 original version in process of being certified
CD2 
c 17-mar-94 gaz don't zero out b
c 12/14/94 gaz zeroed out ithic , imdnode
c 12/16/94 gaz don't zero out nop(not alloc till slvesu)
c 12/22/94 gaz don't zero istrw,sx,sxs,nelm
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
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
CD4
CD4   The definitions for the parameters, scalar, and dynamic variables are
CD4   found in common files comai.h, combi.h, comci.h, comdi.h, comdti.h,
CD4   comei.h, comfi.h, comgi.h, comhi.h, comii.h, comji.h, comrxnb.h, 
CD4   comrxni.h, davidi.h
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
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
CD5   i               INT      Loop index
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
CD8 This routine initializes variables to default values.
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9 2.8   Provide Multiple Realization Option
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
CPS BEGIN  data
CPS 
CPS   assign a value to variables in common 
CPS   
CPS   assign thermodynamic coefficients
CPS   
CPS   FOR each macro flag
CPS       initialize macro read flag to false
CPS   END FOR
CPS   
CPS END  data
CPS
C***********************************************************************

      use comai
      use combi
      use comchem
      use comci
      use comdi
      use comdti
      use comei
      use comfi
      use comflow
      use comgi
      use comhi
      use comii
      use comji
      use compart
      use comrtd
      use comrxni
c gaz 050421      
      use comsi, only : ihms
      use comsplitts
      use comwt, only : head_id, irich
      use comzone
      use davidi
c gaz 123020      
      use com_exphase
c gaz 041222      
      use commass_AWH, only : imass_phase
      implicit none

      integer i

c initialize scalar integers in comzone
      node_azones = 0
      ozflag = 0
      pflag = 0
      tflag = 0
      hflag = 0
      eflag = 0
c initialize number of saved zones      
      izone_sv_cnt = 0
c gaz 061323 num_sv_zones initialized in scanin
c      num_sv_zones = 0
c initialize number of flux countour files       
      icflux = 0
      icconc = 0
c initialize soilvision output to false        
      sv_combine = .false.
c initialize character and boolean variables in comai

      accm = '    '
      hist_flag = .FALSE.
      out_zones = .FALSE.
      boun_out  = .TRUE.
      connect_out = .FALSE.

c initialize scalar integers in comai

      iab    = 0
      iac    = 0
      iaccmx = 0
      iad    = 0
      iad_up = 1000
      iadif  = 0
      iamm   = 0
      iamx   = 500
      icapp  = 0
      ice    = 0
      icf    = 0
! Set in scanin
!      icgts  = 0
!      iccen  = 0
      ichng  = 0
      isplitts  = 0
      icont = 0
      igdpm_rate = 0
c      ico2   = 0
c ico2 set in scanin.f
      icontr = 0
      ics    = 0
      idof   = 2
      ifinsh = 0
      iflag  = 0
      igrav  = 0
      ihf    = 0
      ihs    = 0
      iprtout    = 0
      ilt    = 0
      impf   = 0
      intg   = -1
      ipest = 0
      iwt_uz = 0
! Read in scanin
!      ihead = 0
      ichead = 0
      ivar = 0
c gaz 052322 added ivar_switch
      ivar_switch = 0
c gaz 061422      
      nact_elem_geo = 0
      iexrlp = 0
      iboun = 0
      iflux_ts = 0
      inobr = 0
      iporos = 0
      ipqz   = 0
      ipsat  = 0
      ireord = 0
      irpd   = 0
      interface_flag = 0
      isubbc = 0
      iter   = 0
      itert  = 0
      itotal = 0
      itotals = 0
      itsat  = 0
      iupk   = 0
      ivapl  = 0
c 23-Feb-12 Default value for thermal conductivity models
      ivcond = 0
      ivrock = 0
      ivf    = 1
      ivfcal = 0
      ithic  = 0
      imdnode= 0
      isubmodel = 0
      iw     = 0
      iwelb  = 0
      l      = 0
      lda    = 0
      m      = 0
      m2     = 0
      m3     = 0
      maxit  = 0   
      mink   = 0
      maxsolve = 0
      icons  = 1000
      istea_pest = 0
      mlz    = 1
      iaprf = 0
      if(irun.eq.1) n      = 0
      nbnd   = 0
      ncntr  = 10000000
      if(irun.eq.1) nei    = 0
      neigh  = 0
      if(irun.eq.1) nemx   = 0
      nflx   = 0
c gaz 013014 nflxz = 0 moved to scanin (where it can change)
c      nflxz   = 0
      if(irun.eq.1) ni     = 0
      nhist  = 1
      nicg   = 1
      nsave  = 0
      nstep  = 0
      ntty   = 0
c zero iden_vis and ideng_vis
      iden_vis = 0
      ideng_vis = 0
      bin_flag = 0
      app_flag = 0
      flux_flag = 'no fluxes  '
      maxmix_flag = 0
      sat_ich = 0
c gaz 111520 initialize new variable iter_intrvl  
      iter_intrvl = 0
      header_flag = 'new'
C *** integer in comrxni.h
C *** integers in davidi.h 
! Read in scanin
!      iback  = 0
!      icoupl = 0 
!      irdof  = 0  ! Needs to be kept from scanin
!      islord = 0
      itest  = 0
      ncolsizes = 0
      ncolspec = 0
      igrad = 0
      iconv = 0
      ngrad = 0
      nconv = 0
      isteady = 0
      ipini = 0
c gaz 052323 initialized in scanin
c      jswitch = 0
      joff = 0
C from comwt
      irich = 0
      ifdm_elem = 0
      ibcfar = 0
      nr_stop = 0
      
c gaz 082622 nr_completed indicates status of NR update      
      nr_completed = 0
c gaz 050223 initialized  imass_chk
      imass_chk = 0     
c initialize scalar reals in comai
      aiaa   =   1.0
      am0    =   0.0
      amass  =   0.0
      ame    =   0.0
      an0    =   0.0
      asteam =   0.0
      astmo  =   0.0
      aw     =   0.0
      awc    =   0.0
      awt    =   0.0
      ay     =   0.0
      ayc    =   0.0
      contim =   1.0d+30
      contim_rip =   1.0d+30
      day    =   0.0
      daycf  =   0.0
      daycm  =   0.0
      daycmm =   0.0
      daycmx =   0.0
      daycs  =   0.0
      dayhf  =   0.0
      dayhs  =   0.0
      daymax =  30.0
      daymin =   1.0d-05
      daynew =   0.0
      days   =   0.0
      dayscn =   0.0
      daysi  =   0.0
      daysp  =   0.0
      delat  =   0.0
      delpt  =   0.0
      delst  =   0.0
      deltt  =   0.0
      den    =   0.0
      dene   =   0.0
      depcng =   0.0
      dife   =   0.0
      difm   =   0.0
      ditnd  =   0.0
      dnwgt  =   0.0
      dnwgta =   0.0
      dtot   =   0.0
      dtotc  =   0.0
      dtotdm =   0.0
c gaz 121319 updated 100723   (1.d-5, overwritten if Richard's Eq)  
      pchng = 1.d-3
c gaz 010924 for h2 amanzi problem
c      schng = 1.d-3
       schng = 1.d-4
c gaz 100318
      initdata_pad = 0
C new variables in comsplitts
      dtot_split = 0.0
      dtot_next = 0.0
      explicit_factor = 0.0
      emiss  =   0.0
      epc    =   0.0
      epe    =   0.0
      eps    =   0.0
      f0     =   0.0
      fdum   =   0.0
      fdum1  =   0.0
      fimp   =   0.0
      g1     =   1.0d-06
      g2     =   1.0d-06
      g3     =   1.0d-03
      grad2  =   0.0
      gradnt =   0.0
      grav   =   0.0
      histime =  1.d+30
C from comwt
      head_id =  0.0
      overf  =   0.0
      pein   =   0.0
C added phi_inc to comai so it only needs to be computed in one place
      phi_inc =  0.0
c gaz 110919 initialized perf,tref      
      pref = 0.101325
      tref = 20.0
      pow    =   0.0
      qt     =   0.0
      qte    =   0.0
      qtot   =   0.0
      qtote  =   0.0
      qtotei =   0.0
      qtoti  =   0.0
      quad   =   0.0
      rnmax  =   1.0d+11
      str    =   1.0
      strd   =   1.0
      sv     =   0.0
      teoutf =   0.0
      tims   =   0.0
      tin    =   0.0
      tin0   =   0.0
      tin1   =   0.0
      tmch   =   1.0d-09
      time_ch=   0.0
      tort   =   0.0
      toutfl =   0.0
      upwgt  =   1.0
      upwgta =   1.0
      vlmax  =   0.0
      vvmax  =   0.0
      vtot   =   0.0
      tol_phase = 1.e-4
      strd_iter = 0.99
      p_stop = 0.0
      pa_stop = 0.0
      t_stop = 0.0
      s_stop = 0.0 
      s2_stop = 0.0
      co2f_stop = 0.0
      h_stop = 0.0  
c gaz 032722 initialize more variables used in nr_stop_ctr.f  
      eqwm_stop =0.0
      eqen_stop = 0.0
      eqnm_stop = 0.0
      node_water_err_max = 0
      node_energy_err_max = 0 
      node_ngas_err_max = 0
c gaz 112721 added   critical pt parameters - used in h2o EOS tables   
      pcrit_h2o_true=22.064d0
      tcrit_h2o_true=373.946
      pcrit_h2o=22.064d0
      tcrit_h2o=373.946
c gaz 112021 added tolerance parameters for EOS sampling      
      p_tol = 1.d-19
      pc_tol = 1.d-19
      t_tol = 1.d-19      
c gaz 120421 for fluid_prop_control.f      
      eval_test_h2o = 0      
c logical variables in comai
      wflux_flag = .false.
      vflux_flag = .false.
      eflux_flag = .false.

c zero out include combi
      if(irun.eq.1) then
         ka=0
         nar=0
         nelm=0
         nelmdg=0
         cord=0
         dp=0
         dr=0
         eta=0
         exci=0
         si=0
         sx1=0
         w=0
         wp=0
         wr=0
         wx=0
         wxp=0
         wxr=0
         wy=0
         wyp=0
         wyr=0
         wz=0
         wzp=0
         wzr=0
         xd=0
         xt=0
         yd=0
         yt=0
         zd=0
         zt=0
      end if
c zero out include comci
      akc=0
      danl=0
      danv=0
c Don't do here because they are not allocated yet
c zvd - Initalize dil, div, rolf, rovf, denci - now allocated in allocmem
c      ddvac=0
c      ddvae=0
c      ddvap=0
c      deef=0
c      delef=0
c      delf=0
      denci=0
c      denei=0
c      deni=0
c      denvac=0
c      denvae=0
c      denvap=0
c      depf=0
c      deqh=0
c      devef=0
c      devf=0
c      dgle=0
c      dglp=0
c      dgve=0
c      dgvp=0
      dil=0
c      dile=0
c      dilp=0
      div=0
c      dive=0
c      divp=0
c      dmef=0
c      dmpf=0
c      dpcef=0
c      dq=0
c      dqh=0
c      dqt=0
c      drc=0
c      dstm=0
c      dtpa=0
c      dtpae=0
c      dva=0
c      enlf=0
c      enva=0
c      envf=0
      rolf=0
      rovf=0
c zero out include comchem
      pcpnt=0
      pimm=0
      pvap=0
      cpntprt=0
      cplxprt=0
      immprt=0
      vapprt=0
      cpntnam=' '
      cplxnam=' '
      immnam=' '
      vapnam=' ' 
      ndconv=0
      cflxz = 0
c gaz 040606 to hopefully correct usage
      neg_conc_flag= 1
      if(rxn_flag.ne.0)then
         naqsp=0
         nimsp=0
         nivsp=0
         irxnic=0
         irxnim=0
         irxniv=0
         ckeq=0
         heq=0
         spstoic=0
         cpntgs=0
         idrxn=0
         ifxconc=0
         cpntsv=0
         calcderiv=.false.
         cden = .false.
         ckeqlb=0
         ckmtrn=0
         simmmx=0
         kfor=0
         krev=0
         sticirrv=0
         stimirrv=0
         stivirrv=0
         ckc=0
         cka=0
         decay=0
         biofac=0
         hfac=0
         carbfac=0
         ammfac=0
         phthresh=0
         qm=0
         yield=0
         sarea=0
         pdstic=0
         pd_flag=0
         totaq=0
         cpnt=0
         cplx=0
         rrcplx=0
         rrcpnt=0
         rrimm=0
         rrvap=0
         drdctaq=0
         drdcimm=0
         drdcvap=0
         drcpnt=0
         drtaq=0
         drcplx=0
         drimm=0
         drvap=0
         xjac=0
         dxct=0
         ipiv=0
         resid=0
         ipiv2=0
         nterms=0
         sclmtx=0
         sclfactr=0
         rxnon=0
         temp_model = ' '
         temp_model_kin = ' '
         tcoeff = 0
	   ps_delta_rxn = 0.0d0
         mw_mineral = 0.0d0

      endif
         
c zero out include comdi
      iadsfl=0
      iadsfv=0
      icap=0
      ices=0
      ieos=0
c gaz 072120 zero ieos_sc
      ieos_sc=0
      iieos=0
      iporf=0
      irlp=1
      itrc=0
      itrcdsp=0
c gaz 113020 (initialize AWH henry's law to old model)
      imod_sol = 0
      alpha0 = 1.6111d-4
c gaz 123020 intialize i_ex_update = 0 (no explicit update)
      i_ex_update = 0
c zvd 17-Aug-09 initialize mmodel here
      mmodel=0     
      nskw=0
      nskw2=0
      a1adfl=0
      a1adfv=0
      a2adfl=0
      a2adfv=0
      if (allocated (amgang)) amgang=0
      an=0
      anl=0
      anlo=0
      anv=0
      betadfl=0
      betadfv=0
      cm=0
      cm0=0
      cnsk=0
      cp1f=0
      cp2f=0
      cp3f=0
      cp4f=0
      cpr=0
      dench=0
      dencj=0
      deneh=0
      denej=0
      denh=0
      denj=0
      denr=0
      diffmfl=0
      diffmfv=0
c     dit=0
      dporp=0
      dport=0
      eflow=0
      esk=0
      fc=0
      flx12l=0
      flx12v=0
      head0=0
      iadd=0
      iaddt=0
      icapt=0
      icns=0
      iflx1=0
      iflx2=0
      irlpt=0
c     itc=0
      ivcn=0
      ivcon=0
      npt=0
      i_subsid = 0
      pcp=0
      pflow=0
      if (allocated (pgangi)) pgangi=0
      phi=0
      phini=0
      pho=0
      pnx=0
      pny=0
      pnz=0
      ps=0
      ps_trac=0
      psini=0
      qcin=0
      qcout=0
      qcrxn=0
      qflux=0
      qflxm=0
      qh=0
      rc=0
      rcss=0
      rol0 = 997.
      rp10f=0
      rp11f=0
      rp12f=0
      rp13f=0
      rp14f=0
      rp15f=0
      rp16f=0
      rp17f=0
      rp18f=0
      rp19f=0
      rp1f=0
      rp20f=0
      rp21f=0
      rp22f=0
      rp23f=0
      rp24f=0
      rp2f=0
      rp3f=0
      rp4f=0
      rp5f=0
      rp6f=0
      rp7f=0
      rp8f=0
      rp9f=0
      s=1.
      sii=1.
      sio=1.
      sk=0.
c gaz 120120      
      sk0=0.
      so=0.
      t=0.
      t1sk=0.
      t2sk=0.
      if (allocated (tclx)) tclx=0
      if (allocated (tcly)) tcly=0
      if (allocated (tclz)) tclz=0
      if (allocated (tcvx)) tcvx=0
      if (allocated (tcvy)) tcvy=0
      if (allocated (tcvz)) tcvz=0
      sehdiff=0
      sehdiffv=0
!      thetaresid=0
!      thetaresidv=0
      thx=0
      thy=0
      thz=0
      tini=0
      to=0
      vc1f=0
      vc2f=0
      vc3f=0
      vc4f=0
      vc5f=0
      vc6f=0
      vc7f=0
      vc8f=0
      if(irun.eq.1) vf=0
      if(irun.eq.1) volume=0
      wellim=0
      npn=0
      nsp=0
      phydro=0.0
      sigini=0.0
      thexp=0.0
      tmelt=0.0
      young=0.0

c zero out include comei
      if(irun.eq.1) then
         iirb=0
         irb=0
         nopt=0
         npvt=0
      end if
c gaz 050421 zero out ihms
      ihms = 0
      c=0
      g=0
      h=0
      ss=0
      y=0

c zero out include comfi
      cnlf=0
      cnvf=0
      dcc=0
      dce=0
      dclcf=0
      dclef=0
      dclf=0
      dcp=0
      dcqc=0
      dcqh=0
      dcvcf=0
      dcvef=0
      dcvf=0
      dec=0
      delcf=0
      denpch=0
      denpci=0
      denpcj=0
      deqc=0
      dqpc=0
      devcf=0
      dglc=0
      dgvc=0
      dilc=0
      divc=0
      dmc=0
      dqc=0
      dtpac=0
      eskc=0
      pci=0
      pcio=0
      qc=0
      sici=0
      acner=0.0
      amc=0.0
      difc=0.0
      qtc=0.0
      qtotc=0.0
c gaz 040621      
      qtotin = 0.0

c zero out include comgi
      dfee=0
      dfep=0
      dfepc=0
      dfme=0
      dfmp=0
      dfmpc=0
      dfpce=0
      dfpcp=0
      dfpcpc=0

c zero out include comhi
      a21eef=0
      a21epf=0
      a21mef=0
      a21mpf=0
      a32eef=0
      a32epf=0
      a32mef=0
      a32mpf=0
      apuv1=0
      r3ef=0
      r3mf=0
      rb2ef=0
      rb2mf=0
      tb11=0
      tb12=0
      tb21=0
      tb22=0
      volf1=0
      volf2=0
      wb11=0
      wb12=0
      wb21=0
      wb22=0

c zero out include comji
      it8 = 0
      it9 = 0
      it10 = 0
      it11 = 0
      it12 = 0
      t1 = 0.
      t2 = 0.
      t3 = 0.
      t4 = 0.
      t5 = 0.
      t5v = 0.
      t6 = 0.
      t7 = 0.
      t8 = 0.
      t9 = 0.
      t10 = 0.

c zero out include comii
      cel=0
      cev=0
      crl=0
      crv=0
      cvl=0
      cvv=0
      pmax=0
      pmin=0
      tmax=0
      tmin=0
      ev1= 0.0
      ev2= 0.0
      ev3= 0.0
      ev4= 0.0
      ev5= 0.0
      ev6= 0.0
      ev7= 0.0
      ev8= 0.0
      ev9= 0.0
      ev10= 0.0
      ev11 = 0.0
      ew1 = 0.0
      ew2 = 0.0
      ew3 = 0.0
      ew4 = 0.0
      ew5 = 0.0
      ew6 = 0.0
      ew7 = 0.0
      ew8 = 0.0
      ew9 = 0.0
      ew10 = 0.0
      ew11 = 0.0
      psa0 = 0.0
      psb0 = 0.0
      psta1 = 0.0
      psta2 = 0.0
      psta3 = 0.0
      psta4 = 0.0
      pstb1 = 0.0
      pstb2 = 0.0
      pstb3 = 0.0
      pstb4 = 0.0
      tsa0 = 0.0
      tsb0 = 0.0
      tspa1 = 0.0
      tspa2 = 0.0
      tspa3 = 0.0
      tspa4 = 0.0
      tspb1 = 0.0
      tspb2 = 0.0
      tspb3 = 0.0
      tspb4 = 0.0

c**** coefficients for saturation equation ****
c
c**** date created 01/14/91, max error(percent) 0.350 ****
c**** tmin  10.0 tmax 340.0 pmin 0.00123 pmax 14.59410 ****
c
      tsa0   = -0.25048121d-05
      tspa1  =  0.45249584d-02
      tspa2  =  0.33551528d+00
      tspa3  =  0.10000000d+01
      tspa4  =  0.12254786d+00
      tsb0   =  0.20889841d-06
      tspb1  =  0.11587544d-03
      tspb2  =  0.31934455d-02
      tspb3  =  0.45538151d-02
      tspb4  =  0.23756593d-03
c
c**** coefficients for saturation equation ****
c
c**** date created 01/14/91, max error(percent) 0.300 ****
c**** tmin  10.0 tmax 340.0 pmin 0.00123 pmax 14.59410 ****
c
      psa0   =  0.71725602d-03
      psta1  =  0.22607516d-04
      psta2  =  0.26178556d-05
      psta3  = -0.10516335d-07
      psta4  =  0.63167028d-09
      psb0   =  0.10000000d+01
      pstb1  = -0.22460012d-02
      pstb2  =  0.30234492d-05
      pstb3  = -0.32466525d-09
      pstb4  =  0.0

c**** enthalpy of liquid ****
c
c**** date created 10/31/88, max error(percent) 0.060 ****
c
      pmin(1) =    0.001
      pmax(1) =  110.000
      tmin(1) =   15.000
      tmax(1) =  360.000

      cel(01, 1) =  0.25623465d-03
      cel(02, 1) =  0.10184405d-02
      cel(05, 1) =  0.41769866d-02
      cel(03, 1) =  0.22554970d-04
      cel(06, 1) = -0.21244879d-04
      cel(08, 1) =  0.89557885d-04
      cel(04, 1) =  0.34836663d-07
      cel(07, 1) =  0.25493516d-07
      cel(10, 1) = -0.21720560d-06
      cel(09, 1) =  0.10855046d-06
      cel(11, 1) =  0.10000000d+01
      cel(12, 1) =  0.23513278d-01
      cel(15, 1) = -0.50770309d-02
      cel(13, 1) =  0.48716386d-04
      cel(16, 1) =  0.57780287d-05
      cel(18, 1) = -0.58981537d-04
      cel(14, 1) = -0.19935046d-08
      cel(17, 1) =  0.90972916d-09
      cel(20, 1) =  0.45872518d-08
      cel(19, 1) = -0.12990752d-07
c
c**** density of liquid ****
c
c**** date created 10/31/88, max error(percent) 0.055 ****
c
      pmin(1) =    0.001
      pmax(1) =  110.000
      tmin(1) =   15.000
      tmax(1) =  360.000

      crl(01, 1) =  0.10000000d+01
      crl(02, 1) =  0.17472599d-01
      crl(05, 1) =  0.49564109d-02
      crl(03, 1) = -0.20443098d-04
      crl(06, 1) = -0.40757664d-04
      crl(08, 1) =  0.50330978d-04
      crl(04, 1) = -0.17442012d-06
      crl(07, 1) =  0.50676664d-07
      crl(10, 1) = -0.18383009d-06
      crl(09, 1) =  0.33914814d-06
      crl(11, 1) =  0.10009476d-02
      crl(12, 1) =  0.16812589d-04
      crl(15, 1) =  0.48841156d-05
      crl(13, 1) = -0.24582622d-07
      crl(16, 1) = -0.32967985d-07
      crl(18, 1) =  0.53249055d-07
      crl(14, 1) = -0.17014984d-09
      crl(17, 1) =  0.28619380d-10
      crl(20, 1) = -0.12221899d-09
      crl(19, 1) =  0.30456698d-09
c
c**** viscosity of liquid ****
c
c**** date created 11/07/88, max error(percent) 0.340 ****
c
      pmin(1) =    0.001
      pmax(1) =  110.000
      tmin(1) =    0.0001
      tmax(1) =  360.000

      cvl(01, 1) =  0.17409149d-02
      cvl(02, 1) =  0.18894882d-04
      cvl(05, 1) = -0.31534914d-05
      cvl(03, 1) = -0.66439332d-07
      cvl(06, 1) =  0.11120716d-07
      cvl(08, 1) =  0.28006861d-07
      cvl(04, 1) = -0.23122388d-09
      cvl(07, 1) = -0.48576020d-10
      cvl(10, 1) =  0.47180171d-10
      cvl(09, 1) =  0.23225035d-09
      cvl(11, 1) =  0.10000000d+01
      cvl(12, 1) =  0.10523153d-01
      cvl(15, 1) =  0.29869141d-01
      cvl(13, 1) = -0.22658391d-05
      cvl(16, 1) =  0.21844248d-03
      cvl(18, 1) =  0.41690362d-03
      cvl(14, 1) = -0.31796607d-06
      cvl(17, 1) = -0.87658855d-06
      cvl(20, 1) =  0.22144660d-05
      cvl(19, 1) = -0.25147022d-05
c
c**** enthalpy of vapour ****
c
c**** date created 10/31/88, max error(percent) 0.003 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =   15.000
      tmax(2) =  360.000

      cev(01, 1) =  0.31290881d+00
      cev(02, 1) = -0.10000000d+01
      cev(05, 1) =  0.11319298d-01
      cev(03, 1) =  0.25748596d-01
      cev(06, 1) =  0.20966376d-04
      cev(08, 1) =  0.19206133d-02
      cev(04, 1) =  0.38846142d-03
      cev(07, 1) =  0.74228083d-08
      cev(10, 1) =  0.59104245d-07
      cev(09, 1) = -0.10372453d-03
      cev(11, 1) =  0.12511319d+00
      cev(12, 1) = -0.36061317d+00
      cev(15, 1) =  0.44331611d-02
      cev(13, 1) =  0.58668929d-02
      cev(16, 1) =  0.50902084d-05
      cev(18, 1) =  0.90918809d-03
      cev(14, 1) =  0.99059715d-04
      cev(17, 1) = -0.10812602d-08
      cev(20, 1) = -0.36454880d-06
      cev(19, 1) = -0.26960555d-04
c
c**** density of vapour ****
c
c**** date created 10/31/88, max error(percent) 0.062 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =   15.000
      tmax(2) =  360.000

      crv(01, 1) =  0.15089524d-05
      crv(02, 1) =  0.10000000d+01
      crv(05, 1) =  0.40111210d-07
      crv(03, 1) = -0.10000000d+01
      crv(06, 1) =  0.25625316d-10
      crv(08, 1) =  0.43379623d-01
      crv(04, 1) = -0.16676705d-02
      crv(07, 1) = -0.40479650d-12
      crv(10, 1) = -0.94755043d-04
      crv(09, 1) =  0.24991800d-02
      crv(11, 1) =  0.12636224d+00
      crv(12, 1) = -0.30463489d+00
      crv(15, 1) =  0.59318010d-02
      crv(13, 1) =  0.27981880d-02
      crv(16, 1) =  0.80972509d-05
      crv(18, 1) =  0.53046787d-03
      crv(14, 1) =  0.51132337d-05
      crv(17, 1) = -0.43798358d-07
      crv(20, 1) =  0.48444919d-06
      crv(19, 1) = -0.84916607d-05
c
c**** viscosity of vapour ****
c
c**** date created 11/07/88, max error(percent) 0.160 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =   15.000
      tmax(2) =  360.000

      cvv(01, 1) = -0.13920783d-03
      cvv(02, 1) =  0.98434337d-02
      cvv(05, 1) =  0.27105772d-04
      cvv(03, 1) = -0.51504232d-03
      cvv(06, 1) =  0.84981906d-05
      cvv(08, 1) = -0.25524682d-03
      cvv(04, 1) =  0.62554603d-04
      cvv(07, 1) =  0.34539757d-07
      cvv(10, 1) =  0.12316788d-05
      cvv(11, 1) =  0.10000000d+01
      cvv(12, 1) =  0.10000000d+01
      cvv(15, 1) =  0.10000000d+01
      cvv(13, 1) = -0.10000000d+01
      cvv(16, 1) =  0.10000000d+01
      cvv(18, 1) =  0.10000000d+01
      cvv(14, 1) = -0.10000000d+01
      cvv(17, 1) = -0.22934622d-03
      cvv(20, 1) =  0.25834551d-01
c
c**** enthalpy of liquid ****
c
c**** date created 07/14/89, max error(percent) 0.216 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      cel(01, 2) = -0.28892688d-04
      cel(02, 2) =  0.10155128d-02
      cel(05, 2) =  0.42068671d-02
      cel(03, 2) =  0.38182267d-04
      cel(06, 2) = -0.26722745d-04
      cel(08, 2) =  0.14983417d-03
      cel(04, 2) =  0.29406408d-06
      cel(07, 2) =  0.39965615d-07
      cel(10, 2) = -0.44963038d-06
      cel(09, 2) =  0.11199162d-05
      cel(11, 2) =  0.10000000d+01
      cel(12, 2) =  0.38028489d-01
      cel(15, 2) = -0.62817403d-02
      cel(13, 2) =  0.32800006d-03
      cel(16, 2) =  0.87410801d-05
      cel(18, 2) = -0.11452490d-03
      cel(14, 2) =  0.38164755d-07
      cel(17, 2) =  0.18991534d-08
      cel(20, 2) =  0.19903338d-08
      cel(19, 2) = -0.11341777d-06
c
c**** enthalpy of vapour ****
c
c**** date created 07/14/89, max error(percent) 0.005 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      cev(01, 2) =  0.49023415d+00
      cev(02, 2) = -0.10000000d+01
      cev(05, 2) =  0.86459576d-02
      cev(03, 2) =  0.24474474d-01
      cev(06, 2) =  0.38256791d-04
      cev(08, 2) =  0.16237610d-02
      cev(04, 2) =  0.23476073d-03
      cev(07, 2) =  0.19190905d-07
      cev(10, 2) = -0.78086106d-06
      cev(09, 2) = -0.74126396d-04
      cev(11, 2) =  0.19602927d+00
      cev(12, 2) = -0.35954866d+00
      cev(15, 2) =  0.33114850d-02
      cev(13, 2) =  0.54884993d-02
      cev(16, 2) =  0.12829588d-04
      cev(18, 2) =  0.78784157d-03
      cev(14, 2) =  0.58496026d-04
      cev(17, 2) = -0.20053974d-08
      cev(20, 2) = -0.52896691d-06
      cev(19, 2) = -0.18512345d-04
c
c**** density of liquid ****
c
c**** date created 07/14/89, max error(percent) 0.337 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      crl(01, 2) =  0.10000000d+01
      crl(02, 2) = -0.50430220d-01
      crl(05, 2) =  0.11719555d-01
      crl(03, 2) =  0.12147449d-02
      crl(06, 2) = -0.10272834d-03
      crl(08, 2) =  0.74802254d-03
      crl(04, 2) = -0.29566543d-04
      crl(07, 2) =  0.16483547d-06
      crl(10, 2) = -0.16978281d-05
      crl(09, 2) =  0.17552861d-05
      crl(11, 2) =  0.10020170d-02
      crl(12, 2) = -0.52711077d-04
      crl(15, 2) =  0.11718816d-04
      crl(13, 2) =  0.14548166d-05
      crl(16, 2) = -0.93182060d-07
      crl(18, 2) =  0.72200359d-06
      crl(14, 2) = -0.36472636d-07
      crl(17, 2) =  0.12768238d-09
      crl(20, 2) = -0.14167944d-08
      crl(19, 2) =  0.18887078d-08
c
c**** density of vapour ****
c
c**** date created 07/14/89, max error(percent) 0.072 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      crv(01, 2) =  0.13299942d-04
      crv(02, 2) =  0.10000000d+01
      crv(05, 2) = -0.32791354d-06
      crv(03, 2) = -0.10000000d+01
      crv(06, 2) =  0.21636240d-08
      crv(08, 2) =  0.40896880d-01
      crv(04, 2) = -0.56746973d-02
      crv(07, 2) = -0.38485869d-11
      crv(10, 2) = -0.94741649d-04
      crv(09, 2) =  0.27696827d-02
      crv(11, 2) =  0.12789230d+00
      crv(12, 2) = -0.28996744d+00
      crv(15, 2) =  0.55690966d-02
      crv(13, 2) =  0.26873883d-02
      crv(16, 2) =  0.72603809d-05
      crv(18, 2) =  0.49878874d-03
      crv(14, 2) =  0.33783903d-04
      crv(17, 2) = -0.44323127d-07
      crv(20, 2) =  0.72041771d-06
      crv(19, 2) = -0.13186635d-04
c
c**** viscosity of liquid ****
c
c**** date created 07/14/89, max error(percent) 2.320 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      cvl(01, 2) =  0.17395487d-02
      cvl(02, 2) =  0.18724784d-04
      cvl(05, 2) =  0.64958498d-04
      cvl(03, 2) = -0.15981722d-06
      cvl(06, 2) = -0.31524525d-06
      cvl(08, 2) =  0.17392236d-05
      cvl(04, 2) =  0.15081123d-08
      cvl(07, 2) =  0.16965774d-09
      cvl(10, 2) = -0.17278880d-08
      cvl(09, 2) =  0.27900189d-08
      cvl(11, 2) =  0.10000000d+01
      cvl(12, 2) =  0.12149818d-01
      cvl(15, 2) =  0.56642913d-01
      cvl(13, 2) = -0.82786317d-04
      cvl(16, 2) =  0.20657804d-02
      cvl(18, 2) =  0.11471862d-02
      cvl(14, 2) =  0.90820763d-06
      cvl(17, 2) = -0.86925068d-05
      cvl(20, 2) =  0.47391172d-04
c
c**** viscosity of vapour ****
c
c**** date created 07/14/89, max error(percent) 0.203 ****
c
      pmin(2) =    0.001
      pmax(2) =   20.000
      tmin(2) =    0.500
      tmax(2) =  360.000

      cvv(01, 2) = -0.67484241d-04
      cvv(02, 2) =  0.36800173d-02
      cvv(05, 2) =  0.21890632d-04
      cvv(03, 2) =  0.10553076d-02
      cvv(06, 2) =  0.86065590d-05
      cvv(08, 2) = -0.18856602d-03
      cvv(04, 2) =  0.75936247d-04
      cvv(07, 2) =  0.33613345d-07
      cvv(10, 2) =  0.64108570d-06
      cvv(09, 2) = -0.44264826d-05
      cvv(11, 2) =  0.10000000d+01
      cvv(12, 2) = -0.10000000d+01
      cvv(15, 2) =  0.10000000d+01
      cvv(13, 2) = -0.10000000d+01
      cvv(16, 2) =  0.10000000d+01
      cvv(18, 2) = -0.10000000d+01
      cvv(14, 2) = -0.16348067d+00
      cvv(17, 2) = -0.25908581d-03
      cvv(20, 2) =  0.11181278d-01

c  Initialize values for macro flag parameters

      do i = 1, nmacros
         macroread(i) = .FALSE.
      end do

c  Initialize nbits -- number of bits used with setbit

      nbits = 256
      upwind_l=0
      upwind_v=0

c     Initialize restarting flag for particle tracking to false
      restarting = .false.

      end
