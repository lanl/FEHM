      subroutine chemod(dt,tol_value)
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
CD1 To compute the reaction terms
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2 
CD2 Initial implementation: ?, Programmer: Hari Viswanathan
CD2 new reaction method using chemical components
CD2
CD2 Subrxn 8 fixed to handle precipitation-dissolution reactions
CD2 by Glenn Hammond 4/99
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/chemod.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:06   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:08   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:10 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type   Use         Description
CD3 
**********************************************************************
CD4 Global Subprograms
CD4
CD4 Name            Type     Description
CD4
CD4 speciate        N/A      perform equilibrium speciation 
CD4                          of aqueous reactions
CD4 derivxu         N/A      calculate jacobian for speciation
CD4 subrxn1         N/A      kinetic subroutine for linear sorption
CD4 subrxn2         N/A      kinetic subroutine for langmuir sorption
CD4 subrxn3         N/A      kinetic subroutine for general reaction
CD4 subrxn4         N/A      kinetic subroutine for biodeg (monod 
CD4                          kinetics reaction
CD4 subrxn5         N/A      radioactive decay reaction
CD4 subrxn6         N/A      Henry's law reaction
CD4 subrxn7         N/A      kinetic precipitation dissolution reaction
CD4 subrxn8         N/A      kinetic precipitation dissolution reaction 
CD4                          with rates based on free-ion concentrations
CD4
*****************************************************************
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5
*****************************************************************
CD6
CD6  REQUIREMENTS TRACEABILITY
CD6
CD6  2.3.4 Solute-transport equations
CD6  2.4.5 Adsorbing solutes
CD6  2.4.6 Multiple, interacting solutes
CD6
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS AND REFERENCES
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************

      use comai
      use combi
      use comchem
      use comci
      use comdi
      use comrxni
      use comdti
	use comco2
      use davidi, only : irdof
      implicit none
      integer info, tol_value
      integer ic,ic2,ic3,im,iv,ix,in,mi
      integer irxn,i,j,isp,itemp,itemp2,matnum,mat
      integer isolute
      real*8 dt,dum
      real*8 danl_subst,danv_subst,stemp,disco2_flow


      do in=1,n0
         if(ndconv(in).eq.0)then
            do isolute = 1,nspeci
               mi=in+(isolute-1)*n0
               rc(mi)=0
            enddo
            do ic = 1,ncpnt
               do itemp = 1,nderivs(ic)
                  ic2=drcpos(ic,itemp)
                  matnum=(ic2+ncpnt*(ic-1))
                  mat = in + (matpos(matnum)-1)*n0
                  drdctaq(mat)=0
               enddo
            enddo
            do im = 1,nimm
               drdcimm(im,in) = 0
            enddo
            do iv = 1,nvap
               drdcvap(iv,in) = 0
            enddo

c.......Proceed if complexation is involved.
            if(ncplx.ge.1) then

               do ic=1,ncpnt
                  mi = in+(pcpnt(ic)-1)*n0
                  totaq(ic)=an(mi)
               enddo

c .... Handle dissolved CO2 case
c icarb = 1 for co2 flow, should add check to make sure co2 solubility is turned on
c               if(icarb.eq.1.and.co2_couple.eq.1)then
!               write(*,*) 'new code'
c loop over aqueus components
c                 do ic = 1,ncpnt
c	              mi = in+(pcpnt(ic)-1)*n0
c check that we want couple co2 flow solubility to the trac macro
c		          if (cpntnam(ic).eq.'H2CO3')then
c convert mass fraction to moles/kg water 
c	                  disco2_flow = (yc(in)/(1-yc(in)))/0.044
c if the flow solution solubility is above a specified tolerance we need to adjust the
c total h2co3
c                    write(iout,*) 'inside ',in,disco2_flow,totaq(ic)
c	                  if(disco2_flow.gt.1e-5)then
c		                 totaq(ic) = disco2_flow+cpnt(ic)
c                           totaq(ic) = max(disco2_flow,totaq(ic))
c                          an(mi)=totaq(ic)
c no need to respeciate since we are just changing the total carbonate conc
c	                     ifxconc(ic)=1
c	                  endif
c	                  call speciate(in,info,tol_value)
c	                  ifxconc(ic)=2
c                       dum = cpnt(ic)
c                        Do i=101,ncplx+100
c                          dum = dum + cplx(i)*spstoic(i,ic)
c                       Enddo
c                        totaq(ic) = dum
c
c		          endif 
c	            enddo
c	          endif


c.........Speciation
               call speciate(in,info,tol_value)

               if (tol_value > 0) return

c........Save each aq. component concentration for Node IN
               do ic=1,ncpnt
c Added Dec 30, 1999 GEH
                  if (cpnt(ic) <= 0.d0.and.neg_conc_flag .eq. 0) then
c                    negative free-ion concentrations not allowed
                     if (iout .ne. 0) write(iout,*) 
     1                    'WARNING: Negative concentration, ',
     2                    in, ic, cpnt(ic), ' Save = ', 
     3                    cpntsv(ic,in)
                     if (iptty .ne. 0)  write(iptty,*)
     1                   'WARNING: Negative concentration, ',
     2                   in, ic, cpnt(ic), ' Save = ', 
     3                   cpntsv(ic,in)
	               tol_value = 2
	               return
	            endif
c End add
                  cpntsv(ic,in)=cpnt(ic)
               enddo

c.........If speciation succeeded, generate DXCT
               If(info.eq.0) then
                  Call derivxu
               else
                  do ic=1,ncpnt
                     do i=1,ncpnt
                        dxct(ic,i)=0.d0
                     enddo
                  enddo
               Endif

c.......For non-complexation system
            else
               
               do ic=1,ncpnt
                  do i=1,ncpnt
                     dxct(ic,i)=0.d0
                  enddo
                  dxct(ic,ic)=1.d0
                  mi = in+(pcpnt(ic)-1)*n0
                  cpnt(ic)=an(mi)
               enddo
               
            endif

            do irxn=1,numrxn
               
               do i=1,ncpnt
                  rrcpnt(i)=0.d0
                  do j= 1, ncpnt
                     drcpnt(i,j)=0.d0
                  enddo
               enddo
               
               do i=101,ncplx+100
                  rrcplx(i)=0.d0
                  drcplx(i)=0.d0
               enddo
               
               do i=1,nimm
                  rrimm(i)=0.d0
                  drimm(i)=0.d0
               enddo
               
               do i=1,nvap
                  rrvap(i)=0.d0
                  drvap(i)=0.d0
               enddo
               
c.........Call user-defined subroutines
               
               if(rxnon(irxn,in).eq.1)then
                  goto (210,220,230,240,250,260,270,280) idrxn(irxn)
               else
                  goto 299
               endif
               
 210           call subrxn1(dt,in,irxn)
               goto 299
               
 220           call subrxn2(dt,in,irxn)
               goto 299
               
 230           call subrxn3(dt,in,irxn)
               goto 299
               
 240           call subrxn4(dt,in,irxn)
               goto 299
               
 250           call subrxn5(dt,in,irxn)
               goto 299

 260           call subrxn6(dt,in,irxn)
               goto 299

 270           call subrxn7(dt,in,irxn)
               goto 299

 280           call subrxn8(dt,in,irxn)
               goto 299
               
 299           continue
               if (irdof .ne. 13) then
cHari 3/26/08
                  stemp = min(strac_max,s(in))
                  danl_subst = max(rolf(in)*stemp,rtol)
               else
                  danl_subst = max(rolf(in),rtol)
               endif
               do isp=1,naqsp(irxn)
                  ic=irxnic(irxn,isp)
                  if(ic.gt.100) then
                     ix = ic
                     do i = 1,NCPNT
                        mi=in+(pcpnt(i)-1)*n0
                        rc(mi)=rc(mi)-(ps_trac(in)*danl_subst*
     2                         sx1(in)*spstoic(ix,i)*
     3                         rrcplx(ix)/3600)
                     enddo
                     do ic2=1,NCPNT
                        dum = 0.d0
                        do j=1,NCPNT
                           if(cpnt(j).gt.0.d0) then
                              dum=dum+(spstoic(ix,j)*CPLX(ix)/
     2                             CPNT(j))*dxct(j,ic2)
                           endif
                        enddo
                        do itemp = 1, nderivs(ic2)
                           ic = drcpos(ic2,itemp)
                           matnum = (ic2+ncpnt*(ic-1))
                           mat = in + (matpos(matnum)-1)*n0
                           drdctaq(mat)=drdctaq(mat)-(ps_trac(in)*
     2                          danl_subst*sx1(in)*spstoic(ix,ic)*
     3                          drcplx(ix)*dum/3600)
                        enddo
                        
                     enddo
                  else                  
                     mi=in+(pcpnt(ic)-1)*n0
                     if (irdof .ne. 13) then

cHari 3/26/08
                        stemp = min(strac_max,s(in))
                        danl_subst = max(rolf(in)*stemp,rtol)
                     else
                        danl_subst = max(rolf(in),rtol)
                     end if
                     rc(mi)=rc(mi)-(ps_trac(in)*danl_subst*
     2                    sx1(in)*rrcpnt(ic)/3600)
                     do itemp= 1, nderivs(ic)
                        ic2 = drcpos(ic,itemp)
                        matnum = (ic2+ncpnt*(ic-1))
                        mat = in + (matpos(matnum)-1)*n0 
                        do itemp2 = 1,nderivs(ic)
                           ic3 = drcpos(ic,itemp2)
                           drdctaq(mat)=drdctaq(mat)-(ps_trac(in)*
     2                          danl_subst*sx1(in)*drcpnt(ic,ic3)*
     3                          dxct(ic3,ic2)/3600)
                        enddo
                     enddo
                  endif
               enddo
               
               do isp=1,nimsp(irxn)
                  im=irxnim(irxn,isp)
                  mi=in+(pimm(im)-1)*n0
                  drdcimm(im,in)=drdcimm(im,in)+drimm(im)
                  rc(mi)=rc(mi)+rrimm(im)*denr(in)*sx1(in)/3600
               enddo
               
               if (nivsp(irxn) .gt. 0) then
                  do isp=1,nivsp(irxn)
                     iv=irxniv(irxn,isp)
                     mi=in+(pvap(iv)-1)*n0
                     drdcvap(iv,in)=drdcvap(iv,in)+drvap(iv)
cHari 3/26/08
                     stemp = min(strac_max,s(in))
                     danv_subst = max( rovf(in)*(1-stemp),rtol)
                     rc(mi)=rc(mi)-(ps_trac(in)*danv_subst*sx1(in)*
     2                    rrvap(iv)/3600)
                     drdcvap(iv,in)=drdcvap(iv,in)-(ps_trac(in)*
     2                    danv_subst*sx1(in)*drvap(iv)/3600)

                  enddo
               end if

c.......End of loop w.r.t. IRXN
            enddo

         endif

c.....End of loop w.r.t. IN
      enddo

      return
      end
c **********************************************************************
c <<< Subroutine SPECIATE >>>

      subroutine speciate(in,info,tol_value)
      use comai, only : ierr
      use comchem
      use comdi
      use comdti
	use comco2
      implicit none
      integer itrmxnr
      parameter(itrmxnr=100)

      integer in,info,job,nsolve,nreset,numitr
      integer itr,i,j,k
	integer tol_value

      real*8 weight,tkeq(200)
      real*8 err,dum,denom
      real*8 gas_const
      parameter (gas_const=8.314)
      real*8 temp_conv,stemp, disco2_flow
      parameter (temp_conv=273.16)
      character*1 trans


c====================================================================
c       NCPNT = number of components in system
c       TOTAQ = vector of total analytical component concentrations
c       CPNT = vector of component concentrations 
c       CPLX = vector of complex concentrations
c       CKEQ = vector of equilibrium constants for complex formation
c       CPLXSPSTOICH = matrix of complex stoichiometries
c       IFXCONC = vector signaling fixed component concentrations
c               0 = Not fixed
c               1 = Fixed at TOTAQ concentration
c       RESID = vector of mass balances on components
c       XJAC = Jacobian Matrix
c       RSDMAX = error tolerance for convergence criterion
c       MAXNR = maximum number of iterations allowed
c====================================================================

      job=1
      info=0
      TRANS='N'
      nreset=0
      numitr=2
      weight=1.d0

c====================================================================
c   Set initial guesses for component concentrations.
c   Note: MINTEQ uses 1.0e-9 moles/L as a standard initial guess.
c         In the Co/EDTA/Fe problem, Dan Pastor used an initial 
c         guess of 1.0e-10 M for Co and EDTA and an initial guess of
c         1.0e-15 M for Fe. 
c====================================================================
      do i = 101,ncplx+100
         if(temp_model(i).eq.'l')then 
            tkeq(i)=heq(i,1)+heq(i,2)*t(in)+heq(i,3)*t(in)**2
            tkeq(i)=10**tkeq(i)
	elseif (temp_model(i) .eq. 't') then
		tkeq(i) = heq(i, 1) + heq(i, 2) * t(in) + heq(i, 3) / t(in)
     2		+ heq(i, 4) * log10(t(in)) + heq(i, 5) / t(in) ** 2
		tkeq(i) = 10 ** tkeq(i)
         else
            tkeq(i)=ckeq(i)*exp((heq(i,1)/gas_const*
     2           (1/(25+temp_conv)-1/(t(in)+temp_conv))))
         endif
      enddo
      NSOLVE = 0

      Do j=1,ncpnt
         If((totaq(j).le.1.d-90).and.(ifxconc(j).eq.0)) then
            totaq(j) = 0.d0
            cpnt(j) = 0.d0
            NSOLVE = NSOLVE + 1
	   Else If (ifxconc(j).eq.1) then

            cpnt(j)=totaq(j)
c	   Else if (ifxconc(j).eq.1.and.icarb.eq.1.and.ico2dis(in).eq.1) then
c         Else if (ifxconc(j).eq.0.and.icarb.eq.1) then
CHari convert mass fraction co2 (yc) to moles co2/kg water
CHari if mass fraction is x and mole fraction y they are 
CHari related by: x = y/(1-y)*1/MW

            

         Else 
            cpnt(j)=cpntsv(j,in)
            if((ifxconc(j).eq.0).and.((cpnt(j).gt.totaq(j)).or.
     &           (cpnt(j).le.0.d0))) 
     &           cpnt(j)=totaq(j)
         Endif
         resid(j)=0.d0
      Enddo

      If(NSOLVE.eq.ncpnt) then
         info = 1
         Do i=101,ncplx+100
            cplx(i) = 0.d0
         Enddo
         return
      Endif

 120  Continue

      call sucaprx(tkeq,numitr)

 125  Continue

c====================================================================
c   Main Loop
c====================================================================
      Do 110 itr=1,itrmxnr
c====================================================================
c   Use component concentrations to calculate complex concentrations.
c====================================================================
         Do i=101,ncplx+100
            cplx(i)=tkeq(i)
            Do j=1,ncpnt
               if(spstoic(i,j).ne.0)
     &              cplx(i)=cplx(i)*(cpnt(j))**spstoic(i,j)
            Enddo
         Enddo

c====================================================================
c                        Calculate residuals.
c====================================================================
         Do j=1,ncpnt
            If (((ifxconc(j).eq.0).and.(totaq(j).gt.0.d0)).or.
     &           (ifxconc(j).ge.2)) then
               resid(j)=0.d0
               Do i=101,ncplx+100
                  resid(j)=resid(j)+spstoic(i,j)*cplx(i)
               Enddo
               resid(j)=resid(j)+cpnt(j)-totaq(j)
            Else
               resid(j)=0.d0
            Endif
         Enddo

c====================================================================
c                      Calculate relative error.
c====================================================================
         ERR=0.d0
         Do j=1,ncpnt
            DENOM=dabs(resid(j)+totaq(j))
            If(dabs(totaq(j)).lt.DENOM) DENOM=totaq(j)
            If (DENOM.eq.0.d0) DENOM = 1.d0
            DUM=dabs(resid(j)/DENOM)
            If (DUM.ge.ERR) ERR=DUM
         Enddo

c====================================================================
c   Compute Jacobian Matrix for the equation: 
c       DCPNT[n] = CPNT[n-1]-CPNT[n] = (XJAC[n]^-1)*(RESID[n-1])
c       CPNT[n] = CPNT[n-1]-(XJAC[n]^-1)*(RESID[n-1])
c====================================================================
         Do j=1,ncpnt
            Do k=1,ncpnt
               xjac(j,k)=0.d0
            Enddo
         Enddo

         Do j=1,ncpnt
            If((ifxconc(j).ne.1).and.(cpnt(j).gt.0.d0)) then
               Do k=1,ncpnt
                  If (cpnt(k).gt.0.d0) then
                     Do i=101,ncplx+100
                        xjac(j,k) = xjac(j,k) +
     &                       spstoic(i,j)*spstoic(i,k)*cplx(i)
                     Enddo
                     xjac(j,k) = xjac(j,k)/cpnt(k)
                  Endif
               Enddo
            Endif
         Enddo

         Do k=1,ncpnt
            xjac(k,k)=xjac(k,k)+1.d0
         Enddo

c====================================================================
c   Call subroutine DGETRF to get LU factorization of Jacobian.
c====================================================================
         call DGETRF(ncpnt,ncpnt,xjac,ncpnt,ipiv,info)

c====================================================================
c        Check whether residuals satisfy error tolerance.
c     This is done here so that the most recent Jacobian is 
c     calculated and passed on to DERIVXU.
c====================================================================
         If(ERR.le.rsdmax) then
            return
         Elseif (itr.ge.itrmxnr) then
            If(nreset.eq.0) then
               Do j=1,ncpnt
                  If(ifxconc(j).eq.1) then
                     cpnt(j)=totaq(j)
                  Else If ((totaq(j).le.0.d0).and.(ifxconc(j).eq.0)) 
     &                    then
                     cpnt(j)=0.d0
                  Else
                     cpnt(j)=cpntgs(j)
                  Endif
               Enddo
               nreset=1
               info=0
               GO TO 120
            Else If(nreset.eq.1) then
               weight=0.5d0
               nreset=2
               info=0
               GO TO 125
            Else If(nreset.eq.2) then
               Do j=1,ncpnt
                  If(ifxconc(j).eq.1) then
                     cpnt(j)=totaq(j)
                  Else If ((totaq(j).le.0.d0).and.(ifxconc(j).eq.0)) 
     &                    then
                     cpnt(j)=0.d0
                  Else
                     cpnt(j)=cpntgs(j)
                  Endif
               Enddo
               call sucaprx(tkeq,numitr)
               call sclespec(tkeq,in,info,tol_value)
               return
            Endif
            write (ierr, 195)
 195        format (5X,'Newton-Raphson iteration limit exceeded ',
     &           'in speciation subroutine!')
c Added Dec 30, 1999 by GEH
c            stop
           tol_value = 3
	     return
         Endif

c====================================================================
c   Call matrix solver DGETRS to solve for DPCNT in the equation:
c       (RESID[n-1]) = (XJAC[n])*(DCPNT[n])
c
c   Upon output from DGETRS, RESID = DCPNT.
c====================================================================
         If(info.ne.0) then
            If(nreset.eq.0) then
               Do j=1,ncpnt
                  If(ifxconc(j).eq.1) then
                     cpnt(j)=totaq(j)
                  Else If ((totaq(j).le.0.d0).and.(ifxconc(j).lt.2)) 
     &                    then
                     cpnt(j)=0.d0
                  Else
                     cpnt(j)=cpntgs(j)
                  Endif
               Enddo
               nreset=1
               info=0
               GO TO 120
            Endif
            Write(ierr, 8001)
            stop
         Endif

         call DGETRS(TRANS,ncpnt,job,xjac,ncpnt,ipiv,RESID,
     &        ncpnt,info)
 
         Do j=1,ncpnt
            If((ifxconc(j).ne.1).and.(cpnt(j).gt.0.d0)) then
               DUM=cpnt(j)-resid(j)*weight
               If(DUM.le.0.d0) then
                  cpnt(j)=cpnt(j)/10.d0
               Else If((DUM.gt.totaq(j)).and.(ifxconc(j).eq.0)) then
                  cpnt(j)=cpnt(j)/10.d0
               Else
                  cpnt(j)=DUM
               Endif
            Endif
         Enddo
 110  continue

 8001 format(5X,'Speciation Jacobian matrix is singular!')

      RETURN
      END

c********************************************************************
c <<< Subroutine SUCAPRX >>>
c====================================================================
c  In this subroutine, the "continuous fraction" successive 
c  approximation method is used to estimate the uncomplexed 
c  concentrations of all components (before speciating with
c  the Newton-Raphson method).
c====================================================================

      subroutine sucaprx(tkeq,numitr)

      use comchem
      use comdti

      implicit none

      integer numitr,i,j,n
      real*8  denom,dum,err
      real*8  tkeq(200)

c====================================================================
c   Main Loop
c====================================================================
      Do n=1,numitr

         Do i=101,ncplx+100
            cplx(i)=tkeq(i)
            Do j=1,ncpnt
               if(spstoic(i,j).ne.0)
     &              cplx(i)=cplx(i)*(cpnt(j))**spstoic(i,j)
            Enddo
         Enddo

         Do j=1,ncpnt
            If ((ifxconc(j).eq.0).and.(totaq(j).gt.0.d0)) then
               resid(j)=cpnt(j)
               Do i=101,ncplx+100
                  resid(j)=resid(j)+spstoic(i,j)*cplx(i)
               Enddo
            Endif
         Enddo

         err=0.0d0
         Do j=1,ncpnt
            If (ifxconc(j).eq.0) then
               denom=resid(j)
               If (totaq(j).lt.denom) denom=totaq(j)
               If (denom.eq.0.d0) denom = 1.d0
               dum=dabs((resid(j)-totaq(j))/denom)
               If (dum.ge.err) err=dum
            Endif
         Enddo
         If (err.le.rsdmax) return

         Do j=1,ncpnt
            If((ifxconc(j).eq.0).and.(cpnt(j).gt.0.d0)) then
               dum = (totaq(j)/resid(j))*cpnt(j)
               If(dum.le.0.0) then
                  cpnt(j)=cpnt(j)/10.0d0
               Else If(dum.gt.totaq(j)) then
                  cpnt(j)=cpnt(j)/10.0d0
               Else
                  cpnt(j)=dum
               Endif
            Endif
         Enddo

      Enddo

      RETURN
      END

C **********************************************************************
c <<< Subroutine SCLESPEC >>>

      subroutine sclespec(tkeq,in,info,tol_value)

      use comai, only : ierr
      use comchem
      use comdti

      implicit none

      integer itrmxnr
      parameter(itrmxnr=100)

      integer in,info,numitr
      integer i,j,jj,k,kk,l,ll,itr,job,nreset,ncpnt2
	integer tol_value

      real*8 tkeq(200)
      real*8 denom,dum,err,weight

      character*1 trans


      job=1
      info=0
      TRANS='N'
      nreset=0
      weight=1.d0

      ncpnt2 = ncpnt*2

      Do i=101,ncplx+100
         cplx(i)=tkeq(i)
         Do j=1,ncpnt
            if(spstoic(i,j).ne.0)
     &           cplx(i)=cplx(i)*(cpnt(j))**spstoic(i,j)
         Enddo
      Enddo

      Do j=1,ncpnt
         If (((ifxconc(j).eq.0).and.(totaq(j).gt.0.d0)).or.
     &        (ifxconc(j).ge.2)) then
            dum = 0.d0
            Do i=101,ncplx+100
               dum=dum+spstoic(i,j)*cplx(i)
            Enddo
            resid(j)=dum+cpnt(j)-totaq(j)
         Else
            resid(j)=0.d0
         Endif
      Enddo

c====================================================================
c                  SCALING ALGORITHM
c====================================================================
      Do k=1,ncpnt
         If((ifxconc(k).eq.0).and.(totaq(k).le.0.d0)) then
            nterms(k) = 1.d0
         Else
            nterms(k) = 2.d0
            Do i=101,ncplx+100
               If ((dabs(spstoic(i,k)).gt.0.d0).and.(cplx(i).gt.0.d0))
     &              nterms(k) = nterms(k) + 1.d0
            Enddo
         Endif
      Enddo

      Do l=1,ncpnt2
         Do k=1,ncpnt2
            sclmtx(l,k) = 0.d0
         Enddo
         sclfactr(l) = 0.d0
      Enddo

      Do l=1,ncpnt
         ll = l + ncpnt
         If (((ifxconc(l).eq.0).and.(totaq(l).le.0.d0)).or.
     &        (ifxconc(l).eq.1)) then
            sclmtx(l,l) = 1.d0
            sclmtx(ll,ll) = 1.d0
         Else
            sclmtx(l,l) = nterms(l)
            sclmtx(l,ll) = 1.d0
            sclmtx(ll,ll) = 1.d0
            sclmtx(ll,l) = 1.d0
            Do k=1,ncpnt
               Do i=101,ncplx+100
                  If((cplx(i).gt.0.d0).and.(dabs(spstoic(i,k)).gt.0.d0))
     &                 then
                     kk = k + ncpnt
                     If((dabs(spstoic(i,l))).gt.0.d0) sclmtx(l,kk) = 
     &                    sclmtx(l,kk) + spstoic(i,k)
                     sclmtx(ll,k) = sclmtx(ll,k) + spstoic(i,l)
                     Do j=1,ncpnt
                        jj = j + ncpnt
                        sclmtx(ll,jj) = sclmtx(ll,jj) + 
     &                       spstoic(i,l)*spstoic(i,j)
                     Enddo
                  Endif
               Enddo
            Enddo
         Endif
      Enddo

      Do l=ncpnt2+1,ncpnt2
         sclmtx(l,l) = 1.d0
      Enddo

      Do j=1,ncpnt
         jj = ncpnt + j
         If (((ifxconc(j).eq.0).and.(totaq(j).gt.0.d0)).or.
     &        (ifxconc(j).ge.2)) then
            sclfactr(j) = dlog10(dabs(totaq(j)))
            Do i=101,ncplx+100
               If (cplx(i).gt.0.d0) then
                  If (dabs(spstoic(i,j)).gt.0.d0) sclfactr(j) = 
     &                 sclfactr(j) + dlog(dabs(spstoic(i,j)*tkeq(i)))
                  Do k=1,ncpnt
                     dum = dabs(spstoic(i,k)*tkeq(i))
                     If (dum.gt.0.d0) sclfactr(jj) = sclfactr(jj) + 
     &                    spstoic(i,j)*dlog(dum)
                  Enddo
               Endif
            Enddo
         Endif
      Enddo

      call DGETRF(ncpnt2,ncpnt2,sclmtx,ncpnt2,ipiv2,info)

      If(info.ne.0) then
         Write(6,8003)
         stop
      Endif

      call DGETRS(TRANS,ncpnt2,job,sclmtx,ncpnt2,ipiv2,sclfactr,
     &     ncpnt2,info)

 120  Continue

c====================================================================
c   Main Loop
c====================================================================
      Do 110 itr=1,itrmxnr

         Do i=101,ncplx+100
            cplx(i)=tkeq(i)
            Do j=1,ncpnt
               if(spstoic(i,j).ne.0)
     &              cplx(i)=cplx(i)*(cpnt(j))**spstoic(i,j)
            Enddo
         Enddo

         Do j=1,ncpnt
            Do k=1,ncpnt
               xjac(j,k)=0.d0
            Enddo
         Enddo

         Do j=1,ncpnt
            If((ifxconc(j).ne.1).and.(cpnt(j).gt.0.d0)) then
               Do k=1,ncpnt
                  If (cpnt(k).gt.0.d0) then
                     Do i=101,ncplx+100
                        xjac(j,k) = xjac(j,k) +
     &                       spstoic(i,j)*spstoic(i,k)*cplx(i)
                     Enddo
                     xjac(j,k) = xjac(j,k)*(10.d0**(sclfactr(j)+
     &                    sclfactr(ncpnt+k)))/cpnt(k)
                  Endif
               Enddo
            Endif
            xjac(j,j)= xjac(j,j) + 10.d0**(sclfactr(j)+
     &           sclfactr(ncpnt+j))
         Enddo

c====================================================================
c        Check whether residuals satisfy error tolerance.
c     This is done here so that the most recent Jacobian is 
c     calculated and passed on to DERIVXU.
c====================================================================

         Do j=1,ncpnt
            If (((ifxconc(j).eq.0).and.(totaq(j).gt.0.d0)).or.
     &           (ifxconc(j).ge.2)) then
               resid(j)=0.d0
               Do i=101,ncplx+100
                  resid(j)=resid(j)+spstoic(i,j)*cplx(i)
               Enddo
               resid(j)=resid(j)+cpnt(j)-totaq(j)
            Else
               resid(j)=0.d0
            Endif
         Enddo

         err=0.d0
         Do j=1,ncpnt
            denom=dabs(resid(j)+totaq(j))
            If(dabs(totaq(j)).lt.denom) denom=totaq(j)
            If (denom.eq.0.d0) denom = 1.d0
            dum=dabs(resid(j)/denom)
            If (dum.ge.err) err=dum
         Enddo

         If(err.le.rsdmax) then
            Do j=1,ncpnt
               Do k=1,ncpnt
                  xjac(j,k)= xjac(j,k)*10**(-(sclfactr(j)+
     &                 sclfactr(ncpnt+k)))
                  call DGETRF(ncpnt,ncpnt,xjac,ncpnt,ipiv,info)
               Enddo
            Enddo
            return
         Elseif (itr.ge.itrmxnr) then
            If (nreset.eq.0) then
               weight = 0.5d0
               nreset = 1
               GO TO 120
            Else If (nreset.eq.1) then
               Do j=1,ncpnt
                  If(ifxconc(j).eq.1) then
                     cpnt(j)=totaq(j)
                  Else If ((totaq(j).le.0.d0).and.(ifxconc(j).eq.0)) 
     &                    then
                     cpnt(j)=0.d0
                  Else
                     cpnt(j)=cpntgs(j)
                  Endif
               Enddo
               weight = 1.0d0
               nreset = 2
               info=0
               GO TO 120
            Else If (nreset.eq.2) then
               weight = 0.5d0
               nreset = 3
               GO TO 120
            Else If (nreset.eq.3) then
               numitr = 10
               call sucaprx(tkeq,numitr)
               nreset = 4
               GO TO 120               
            Else
               write (ierr, 8002)
 8002          format (5X,'Newton-Raphson iteration limit exceeded ',
     &              'in scaled speciation subroutine!')
               Write(ierr, 8005) in
               Write(ierr, *) (totaq(j),j=1,ncpnt)
               Write(ierr, *) (cpnt(j),j=1,ncpnt)
               Write(ierr, *) (sclfactr(k),k=1,ncpnt2)
c Added Dec 30, 1999 by GEH
c               stop
               tol_value = 3
	         return
            Endif
         Endif

         call DGETRF(ncpnt,ncpnt,xjac,ncpnt,ipiv,info)

         Do j=1,ncpnt
            resid(j)=resid(j)*10.d0**sclfactr(j)
         Enddo

c====================================================================
c   Call matrix solver DGETRS to solve for DPCNT in the equation:
c       (RESID[n-1]) = (XJAC[n])*(DCPNT[n])
c
c   Upon output from DGETRS, RESID = DCPNT.
c====================================================================
         If(info.ne.0) then
            If (nreset.lt.2) then
               Do j=1,ncpnt
                  If(ifxconc(j).eq.1) then
                     cpnt(j)=totaq(j)
                  Else If ((totaq(j).le.0.d0).and.(ifxconc(j).eq.0)) 
     &                    then
                     cpnt(j)=0.d0
                  Else
                     cpnt(j)=cpntgs(j)
                  Endif
               Enddo
               weight = 1.0d0
               nreset = 2
               info=0
               GO TO 120
            Else
               Write(ierr, 8001)
               stop
            Endif
         Endif

         call DGETRS(TRANS,ncpnt,job,xjac,ncpnt,ipiv,resid,
     &        ncpnt,info)
 
         Do j=1,ncpnt
            If((ifxconc(j).ne.1).and.(cpnt(j).gt.0.d0)) then
               dum=cpnt(j)-resid(j)*10.d0**(sclfactr(ncpnt+j))*weight
               If(dum.le.0.d0) then
                  cpnt(j)=cpnt(j)/10.d0
                  err = 1.d0
               Else If((dum.gt.totaq(j)).and.(ifxconc(j).eq.0)) then
                  cpnt(j)=cpnt(j)/10.d0
                  err = 1.d0
               Else
                  cpnt(j)=dum
               Endif
            Endif
         Enddo
 110  continue

 8001 format(5X,'Scaled Speciation Jacobian matrix is singular!')
 8003 format(5X,'Speciation scaling matrix is singular!')
 8005 format(/,5x,'Failure at node ',I3)

      RETURN
      END
c **********************************************************************
c <<< Subroutine varph >>>

      subroutine varph(in)
      use comchem
      use comdti
      implicit none
      integer i,in,info,j
c Added Dec 30, 1999 by GEH
	integer dum_int
	integer num_ifx2, ifxconc2_id(20)
c End add

      real*8 dum


c====================================================================
c   First, set the initial guesses for the free component 
c   concentrations and determine the component number of H+ (idh).
c====================================================================

c Altered to allow for other negative total concentrations - Dec 30,
c 1999 by GEH

      num_ifx2 = 0

      Do j=1,ncpnt
         If (ifxconc(j).eq.2) then
	      num_ifx2 = num_ifx2 + 1
	      ifxconc2_id(num_ifx2) = j
            dum = totaq(j)
	      totaq(j)=10.0d0**(-dum)
            ifxconc(j) = 1
         Endif
      Enddo

c====================================================================
c   Call solver.
c====================================================================

      call speciate(in,info,dum_int)

c====================================================================
c   Reset ifxconc for H+.
c====================================================================
c  Not necessary anymore - 12/30/99 - geh; it is performed below
c      ifxconc(idh) = 2

c====================================================================
c   Determine the total aqueous concentration of H+ (TOTH).
c====================================================================
c Altered to allow for other negative total concentrations - Dec 30,
c 1999 by GEH

      do j=1,num_ifx2

        ifxconc(ifxconc2_id(j)) = 2

        dum = cpnt(ifxconc2_id(j))
        Do i=101,ncplx+100
           dum = dum + cplx(i)*spstoic(i,ifxconc2_id(j))
        Enddo

        totaq(ifxconc2_id(j)) = dum
      enddo

c      dum = cpnt(idh)
c      Do i=101,ncplx+100
c         dum = dum + cplx(i)*spstoic(i,idh)
c      Enddo
c
c      totaq(idh) = dum

      return
      end
c **********************************************************************
c <<< Subroutine derivxu >>>

      subroutine derivxu 
      use comchem
      use comdti
      implicit none
      integer info,job,j,k


      CHARACTER*1  TRANS



      job = 1
      TRANS='T'
      info = 0

      Do j=1,NCPNT
         Do k=1,NCPNT
            DXCT(j,k) = 0.d0
         Enddo
         RESID(j) = 0.d0
      Enddo

      Do j=1,NCPNT
         If (CALCDERIV(j)) then
            RESID(j) = 1.d0
            call DGETRS(TRANS,NCPNT,job,XJAC,ncpnt,IPIV,
     &           RESID,ncpnt,INFO)
            Do k=1,NCPNT
               DXCT(j,k) = RESID(k)
               RESID(k) = 0.d0
            Enddo
         Endif
      Enddo

      return
      end

c **********************************************************************
c <<< Subroutine chckderiv >>>

c User-defined subroutine for generating 'calcderiv'

      subroutine chckderiv
      use comchem
      use comdti
      implicit none


c ----Users can write codes below----

      integer i,j,k,isp

      Do j=1,NCPNT
         CALCDERIV(j) = .FALSE.
      Enddo

      Do i=1,NUMRXN
         do isp=1,naqsp(i)
            k = irxnic(i,isp)
            If(k.gt.100) then
               Do j=1,NCPNT
                  If(spstoic(k,j).ne.0.d0) CALCDERIV(j) = .TRUE.
               Enddo
            Elseif(k.gt.0) then
               CALCDERIV(k) = .true.
            Endif
         enddo
      Enddo

      RETURN
      END

c **********************************************************************
c <<< Subroutine subrxn1 >>>

c     User-defined subroutine #1 (Nodal reaction)

      subroutine subrxn1(dt,in,irxn)
c ----Inherent statements
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn

      real*8 den,dt,por,danl_subst,stemp

c  ...Inherent common blocks

c ----Users can write codes below----

c     Type of reaction:
c        Linear rate-limiting adsorption
      

c  ...User-introduced variables
      integer ic,im,mi
      real*8 oimm,ckeqlb1,ckmtrn1

c  ...Substitutions
      ic=irxnic(irxn,1)
      im=irxnim(irxn,1)
      mi = in+(pimm(im)-1)*n0
      oimm=anlo(mi)
      if(temp_model_kin(irxn).ne.'l')then
         ckeqlb1=ckeqlb(irxn)
      else
         ckeqlb1=tcoeff(irxn,1)+t(in)*tcoeff(irxn,2)+t(in)**2*
     2        tcoeff(irxn,3)
         ckeqlb1=10**ckeqlb1
      endif
      ckmtrn1=ckmtrn(irxn)
cHari 3/26/08
      if (irdof .ne. 13) then
         stemp = min(s(in),strac_max)
         danl_subst = max(rolf(in)*stemp,rtol)
      else
         danl_subst = max(rolf(in),rtol)
      end if
      por = ps_trac(in)
      den = denr(in)

c  ...Calculation for complexes
      if(ic.gt.100) then  

C    ...Compute reaction rate and dr/dCPLX

         rrcplx(ic) = -(den*ckmtrn1) * (ckeqlb1*cplx(ic) -
     &        oimm) / (danl_subst*por*dt*ckmtrn1 + 
     &        den*ckeqlb1)

         drcplx(ic) = - (den*ckmtrn1*ckeqlb1) / 
     &    ( danl_subst*por*dt*ckmtrn1 + den*ckeqlb1) 

c
         rrimm(im)=-(por*danl_subst)/den*rrcplx(ic)
         drimm(im)=0.d0

c  ...Calculation for components
      else

         rrcpnt(ic) = -(den*ckmtrn1) * (ckeqlb1*cpnt(ic) -
     &        oimm) / (danl_subst*por*dt*ckmtrn1 + 
     &        den*ckeqlb1)
        
         drcpnt(ic,ic) = -(den*ckmtrn1*ckeqlb1) / 
     &        ( danl_subst*por*dt* ckmtrn1 + den*ckeqlb1)

         rrimm(im)=-(por*danl_subst)/den*rrcpnt(ic)
         drimm(im)=0.d0

      endif

      return
      end   

c **********************************************************************
c <<< Subroutine subrxn2 >>>

c     User-defined subroutine #2 (Nodal reaction)

      subroutine subrxn2(dt,in,irxn)
c ----Inherent statements
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn,mi
      real*8 den,dt,stemp

c  ...Inherent common blocks

c ----Users can write codes below----

c     Type of reaction:
c        Langmuir adsorption

c  ...User-introduced variables
      integer ic,im
      real*8 oimm,ckeqlb1,ckmtrn1,simmmx1
      real*8 danl_subst,por

      ic=irxnic(irxn,1)
      im=irxnim(irxn,1)
      mi = in+(pimm(im)-1)*n0
      oimm=anlo(mi)
      if(temp_model_kin(irxn).ne.'l')then
         ckeqlb1=ckeqlb(irxn)
      else
         ckeqlb1=tcoeff(irxn,1)+t(in)*tcoeff(irxn,2)+t(in)**2*
     2        tcoeff(irxn,3)
         ckeqlb1=10**ckeqlb1
      endif
      ckmtrn1=ckmtrn(irxn)
      simmmx1=simmmx(irxn)
cHari 3/26/08
      if (irdof .ne. 13) then
         stemp = min(strac_max,s(in))
         danl_subst = max(rolf(in)*stemp,rtol)
      else
         danl_subst = max(rolf(in),rtol)
      end if
      por = ps_trac(in)
      den = denr(in)

      if(IC.gt.100) then  

         rrcplx(ic) = - den*(ckeqlb1*CPLX(IC)*simmmx1 - 
     &        OIMM*(ckeqlb1*CPLX(IC) + 1.d0))/
     &        (danl_subst*POR*dt*(ckeqlb1*CPLX(IC) + 1.d0 +
     &        1.d0/(ckmtrn1*dt))) 

         drcplx(ic) = - DEN*(ckeqlb1*simmmx1*(1.d0+
     &        1.d0/(ckmtrn1*dt)) - 
     &        (1.d0/(ckmtrn1*dt))*ckeqlb1*OIMM)
     &        /(danl_subst*POR*dt*(ckeqlb1*CPLX(IC) + 1.d0
     &        + 1.d0/(ckmtrn1*dt))**2) 

         rrimm(im)=-(por*danl_subst)/den*rrcplx(ic)
         drimm(im)=0.d0

      else

         rrcpnt(ic) = - DEN*(ckeqlb1*CPNT(IC)*simmmx1 - 
     &        OIMM*(ckeqlb1*CPNT(IC) + 1.d0))/
     &        (danl_subst*POR*dt*(ckeqlb1*CPNT(IC) + 1.d0
     &        + 1.d0/(ckmtrn1*dt))) 

         drcpnt(ic,ic) = - DEN*(ckeqlb1*simmmx1*(1.d0+ 
     &        1.d0/(ckmtrn1*dt)) - 
     &        (1.d0/(ckmtrn1*dt))*ckeqlb1*OIMM)
     &        /(danl_subst*POR*dt*(ckeqlb1*CPNT(IC) + 1.d0
     &        + 1.d0/(ckmtrn1*dt))**2) 

         rrimm(im)=-(por*danl_subst)/den*rrcpnt(ic)
         drimm(im)=0.d0

      endif

      return
      end

c **********************************************************************
c <<< Subroutine subrxn3 >>>

c     User-defined subroutine #3 (Nodal reaction)

      subroutine subrxn3(dt,in,irxn)
c ----Inherent statements

c  ...Common maximum array dimensions are written in 'dimpara.f'
c       < contain: maxbrk,maxcplx,maxcpnt,maximm,maxnn,maxny,
c       <          maxprt,maxrxn,maxspcs,maxtbc
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn

      real*8 den,dt,por,stemp
      
c  ...Inherent common blocks

c ----Users can write codes below----

c     Type of reaction:
c        Generic (Reversible or Irreversible) Reaction

c  ...User-introduced variables
      integer ic,im,j,k,l,m,ifor,ifor2,irev, ic2, ic3
      integer mi,irev2,nzerof,nzeror
      real*8 newimm,forrate,revrate,rate
      real*8 danl_subst,wr_conv
      por = ps_trac(in)
      den = denr(in)
cHari 3/26/08
      if (irdof .ne. 13) then
         stemp = min(strac_max,s(in))
         danl_subst = max(rolf(in)*stemp,rtol)
      else
         danl_subst = max(rolf(in),rtol)
      end if
      wr_conv = danl_subst*por/den 
      nzerof = 0
      nzeror = 0
      ifor = 0
      irev = 0
      ifor2 = 0
      irev2 = 0
      forrate = kfor(irxn)
      revrate = krev(irxn)
      If(revrate.le.0.d0) then
         revrate = 0.d0
         nzeror = 2
      Endif

c ====================================================================
c  Loop over aqueous species involved in the reaction.
c     If the stoichiometry is postitive (>= zero), the species is a
c     reactant and contributes to the forward rate.
c     If the stoichiometry is negative (< zero), the species is a
c     product and contributes to the reverse rate.
c
c ====================================================================
      do j=1,naqsp(irxn)
         ic = irxnic(irxn,j)
         If (sticirrv(irxn,j).ge.0.d0) then
            If (ic.gt.100) then
               If (cplx(ic).gt.0.d0) then
                  forrate = forrate*cplx(ic)**sticirrv(irxn,j)
               Else
                  nzerof = nzerof + 1
                  ifor = ic
               Endif
            Else
               If (cpnt(ic).gt.0.d0) then
                  forrate = forrate*cpnt(ic)**sticirrv(irxn,j)
               Else
                  nzerof = nzerof + 1
                  ifor = ic
               Endif
            Endif
         Else If (revrate.gt.0.d0) then
            If (ic.gt.100) then
               If (cplx(ic).gt.0.d0) then
                  revrate = revrate*cplx(ic)**(-sticirrv(irxn,j))
               Else
                  nzeror = nzeror + 1
                  irev = ic
               Endif
            Else
               If (cpnt(ic).gt.0.d0) then
                  revrate = revrate*cpnt(ic)**(-sticirrv(irxn,j))
               Else
                  nzeror = nzeror + 1
                  irev = ic
               Endif
            Endif
         Endif
      enddo

c ====================================================================
c  Loop over immobile species involved in the reaction.
c ====================================================================
      do m=1,nimsp(irxn)
         im = irxnim(irxn,m)
         mi = in+(pimm(im)-1)*n0
         newimm=an(mi)
         If (stimirrv(irxn,m).ge.0.d0) then
            If (newimm.gt.0.d0) then
               forrate = forrate*newimm**stimirrv(irxn,m)
            Else
               nzerof = nzerof + 1
               ifor2 = im
            Endif
         Else If (revrate.gt.0.d0) then
            If (newimm.gt.0.d0) then
               revrate = revrate*
     &              newimm**(-stimirrv(irxn,m))
            Else
               nzeror = nzeror + 1
               irev2 = im
            Endif
         endif
      enddo

c ====================================================================
c  Calculate the overall rate of production of the products 
c     (or equivalently, the rate of consumption of the reactants)
c                     rate = revrate - forrate
c
c  NOTE: We assume that the rate constants (kfor and krev) are
c        chosen to give rate units of [moles/kg water/hr]
c        (i.e. units of aqueous concentration per time).
c        The exception is the case when only immobile species are 
c        involved in the reaction.  In this case, we assume that
c        the rate units are [moles/kg rock/hr].  We correct for
c        this case below.
c ====================================================================
      rate = 0.d0
      If (naqsp(irxn).eq.0) then
         forrate = forrate/wr_conv
         revrate = revrate/wr_conv
      endif
      If (nzerof.eq.0) rate = rate-forrate
      If (nzeror.eq.0) rate = rate+revrate

c ====================================================================
c  Calculate the rates of change in concentration for each of the 
c  species involved in the reaction.
c ====================================================================
      If ((nzerof.eq.0).or.(nzeror.eq.0)) then
         do k=1,naqsp(irxn)
            ic = irxnic(irxn,k)
            If(ic.gt.100) then
               If (sticirrv(irxn,k).eq.0.d0) then
                  rrcplx(ic) = rate
               Else
                  rrcplx(ic) = rate*sticirrv(irxn,k)
               Endif
            Else
               If (sticirrv(irxn,k).eq.0.d0) then
                  rrcpnt(ic) = rate
               Else
                  rrcpnt(ic) = rate*sticirrv(irxn,k)
               Endif
            Endif
         enddo
         
         do l=1,nimsp(irxn)
            im = irxnim(irxn,l)
            If (stimirrv(irxn,l).eq.0.d0) then
               rrimm(im)=rate*wr_conv
            Else
               rrimm(im)=rate*stimirrv(irxn,l)*wr_conv
            Endif
         enddo

      Endif
         
c ====================================================================
c  Calculate the the derivatives of the rates of change in 
c  concentration for each of the species involved in the reaction.
c  The derivatives of the forward and reverse rates are calculated
c  separately.
c ====================================================================
      If (nzerof.eq.1) then
         if(sticirrv(irxn,ifor).eq.0.0)drcpnt(ifor,ifor) = 1e8
         If((ifor.gt.0).and.(sticirrv(irxn,ifor).eq.1.0))then
            If (ifor.gt.100) then
               drcplx(ifor) = -forrate
            Else if (ifor.gt.0) then
               drcpnt(ifor,ifor) = -forrate
            endif
         Else if((ifor2.gt.0).and.(stimirrv(irxn,ifor2).eq.1.d0)) then
            drimm(ifor2) = -forrate*wr_conv
         Endif
      Endif

      If (nzeror.eq.1) then
         If ((irev.gt.100).and.(sticirrv(irxn,irev).eq.-1.d0)) then
            drcplx(irev) = revrate
         Else if ((irev.gt.0).and.(sticirrv(irxn,irev).eq.-1.d0)) then
            drcpnt(irev,irev) = -revrate
         Else if (stimirrv(irxn,irev2).eq.-1.d0) then
            drimm(irev2) = revrate*wr_conv
         Endif
      Endif

      If (nzerof.eq.0) then
         do k=1,naqsp(irxn)
            If (sticirrv(irxn,k).gt.0.d0) then
               ic = irxnic(irxn,k)
               If(ic.gt.100) then
                  drcplx(ic) = -forrate*sticirrv(irxn,k)**2/cplx(ic)

c Note that cross derivatives for complexes are not calculated yet!
c Need to add this!

               Else
                  do l = 1, naqsp(irxn)
                     if(sticirrv(irxn,l).gt.0.d0)then
                        ic2 = irxnic(irxn,l)
                        if(ic2.gt.100)then

c check whether complexes are being handled correctly!

                           do ic3 = 1, ncpnt 
                              drcpnt(ic,ic3) = -spstoic(ic2,ic3)*
     2                             forrate*sticirrv(irxn,k)*
     3                             sticirrv(irxn,l)/cpnt(ic3)
                           enddo
                        else
                           drcpnt(ic,ic2) = -forrate*sticirrv(irxn,k)*
     2                          sticirrv(irxn,l)/cpnt(ic2)
                        endif
                     else
                        ic2 = irxnic(irxn,l)
                        if(ic2.gt.100)then

c check whether complexes are being handled correctly!

                           do ic3 = 1, ncpnt 
                              drcpnt(ic,ic3) = -spstoic(ic2,ic3)*
     2                             revrate*sticirrv(irxn,k)*
     3                             sticirrv(irxn,l)/cpnt(ic3)
                           enddo
                        else
                           drcpnt(ic,ic2) = -revrate*sticirrv(irxn,k)*
     2                          sticirrv(irxn,l)/cpnt(ic2)
                        
                        endif
                     endif
                  enddo
               Endif
            Endif
         enddo

         do l=1,nimsp(irxn)
            im = irxnim(irxn,l)
            mi = in+(pimm(im)-1)*n0
            newimm=an(mi)
            If (stimirrv(irxn,l).gt.0.d0) drimm(im)=-forrate*
     &           stimirrv(irxn,l)**2/newimm*wr_conv
         enddo
      Endif

      If (nzeror.eq.0) then
         do k=1,naqsp(irxn)
            If (sticirrv(irxn,k).lt.0.d0) then
               ic = irxnic(irxn,k)
               If(ic.gt.100) then
                  drcplx(ic) = revrate*sticirrv(irxn,k)**2/cplx(ic)

c Note that cross derivatives for complexes are not calculated yet!
c Need to add this!

               Else
                  do l = 1, naqsp(irxn)
                     if(sticirrv(irxn,l).lt.0.d0)then
                        ic2 = irxnic(irxn,l)

c check whether complexes are being handled correctly!

                        if(ic2.gt.100)then
                           do ic3 = 1, ncpnt 
                              drcpnt(ic,ic3)= -spstoic(ic2,ic3)*revrate*
     2                             sticirrv(irxn,k)*sticirrv(irxn,l)/
     3                             cpnt(ic3)
                           enddo
                        else
                           drcpnt(ic,ic2) = -revrate*sticirrv(irxn,k)*
     2                          sticirrv(irxn,l)/cpnt(ic2)
                        endif
                     else
                        ic2 = irxnic(irxn,l)
                        if(ic2.gt.100)then

c check whether complexes are being handled correctly!

                           do ic3 = 1, ncpnt 
                              drcpnt(ic,ic3) = -spstoic(ic2,ic3)*
     2                             forrate*sticirrv(irxn,k)*
     3                             sticirrv(irxn,l)/cpnt(ic3)
                           enddo
                        else
                           drcpnt(ic,ic2) = -forrate*sticirrv(irxn,k)*
     2                          sticirrv(irxn,l)/cpnt(ic2)
                        endif

                     endif
                  enddo
               Endif
            Endif
         enddo

         do l=1,nimsp(irxn)
            im = irxnim(irxn,l)
            mi = in+(pimm(im)-1)*n0
            newimm=an(mi)
            If (stimirrv(irxn,l).lt.0.d0) drimm(im) = revrate*
     &              stimirrv(irxn,l)**2/newimm*wr_conv
         enddo 
      Endif
c      do mi = 1,naqsp(irxn)
c         do im = 1, naqsp(irxn)
c            write(*,*)mi, im, drcpnt(mi,im)
c         enddo
c      enddo
c      pause
      return
      end
c **********************************************************************
c <<< Subroutine subrxn4 >>>

c     User-defined subroutine #4 (Nodal reaction)

      subroutine subrxn4(dt,in,irxn)
c ----Inherent statements
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      implicit none

c  ...Declarations of variables
      integer in,irxn
      integer mi

      real*8 dt, stemp

c  ...Inherent common blocks

c ----Users can write codes below----

c     Type of reaction:
c        Biological degradation

c  ...User-introduced variables
      integer ia,ib,ic,id,ih,j,k
      real*8 ckc1,cka1,decay1,fac1,fac2,hmax,qm1,rate,yield1
      real*8 newtaqid,newtaqia,newimmib
      integer ico3, inh4
      real*8 fac3, fac4, phmax, xminit1, ph, deriv, substrc
      real*8 dum1,dum2
c  ...Electron donor
      id=irxnic(irxn,1)
      mi = in+(pcpnt(id)-1)*n0
      newtaqid=an(mi)
c  ...Electron acceptor
      ia=irxnic(irxn,2)
      mi = in+(pcpnt(ia)-1)*n0
      newtaqia=an(mi)
c  ...H+
      ih=irxnic(irxn,3)
c  ... CO3-2
      ico3=irxnic(irxn,4)
c  ... NH4+
      inh4= irxnic(irxn,5)
c  ...Microbes
      ib=irxnim(irxn,1)
      mi = in+(pimm(ib)-1)*n0
      newimmib = an(mi)

c  ...Substitutions
      ckc1=ckc(irxn)
      cka1=cka(irxn)
      decay1=decay(irxn)
      fac1=biofac(irxn)
      fac2=hfac(irxn)
      fac3=carbfac(irxn)
      fac4=ammfac(irxn)
      phmax=phthresh(irxn)
      qm1=qm(irxn)
      yield1=yield(irxn)
      xminit1=xminit(irxn)
      if(ih.ne.0)then
         ph= -log10(cpnt(ih))
      else
         ph=0
      endif
      if (ph.gt.phmax) return
      if(nbiofrm(irxn).eq.0)then

         rate=qm1*newimmib*newtaqid*newtaqia
     &        /(ckc1+newtaqid)/(cka1+newtaqia)
         rrcpnt(id)=-rate
         rrcpnt(ia)=-fac1*rate
         rrcpnt(ih)=-fac2*rate
         rrcpnt(ico3)=-fac3*rate
         rrcpnt(inh4)=-fac4*rate
         rrimm(ib)=yield1*rate-decay1*newimmib+decay1*xminit1
         deriv=qm1*newimmib*ckc1*newtaqia
     &        /(ckc1+newtaqid)**2/(cka1+newtaqia)
         drtaq(id)=-deriv

         deriv=qm1*newimmib*newtaqid*cka1
     &        /(ckc1+newtaqid)/(cka1+newtaqia)**2
         drtaq(ia)=-fac1*deriv

         deriv=qm1*newtaqid*newtaqia
     &        /(ckc1+newtaqid)/(cka1+newtaqia)
         drimm(ib)=yield1*deriv-decay1

      else
         substrc = 0.d0
         do j = 1,nbiofrm(irxn)
            If(icbio(irxn,j).gt.100) then
               substrc = substrc + spstoic(icbio(irxn,j),id)*
     &              cplx(icbio(irxn,j))
            else
               substrc = substrc + cpnt(icbio(irxn,j))
            endif
         enddo

         rate=qm1*newimmib*substrc*newtaqia
     &        /(ckc1+substrc)/(cka1+newtaqia)
         rrcpnt(id)=-rate
         rrcpnt(ia)=-fac1*rate
         rrcpnt(ih)=-fac2*rate
         rrcpnt(ico3)=-fac3*rate
         rrcpnt(inh4)=-fac4*rate
         rrimm(ib)=yield1*rate-decay1*newimmib+decay1*xminit1

         deriv = qm1*newimmib*ckc1*newtaqia
     &        /(ckc1+substrc)**2/(cka1+newtaqia)
         do j = 1,nbiofrm(irxn)
            ic = icbio(irxn,j)
            If(ic.gt.100) then
               dum1 = spstoic(ic,id)*cplx(ic)*deriv
               do k=1,ncpnt
                  If ((cpnt(k).gt.0.d0).and.(spstoic(ic,k).ne.0.d0)) 
     &                 then
                     dum2 = spstoic(ic,k)*dum1/cpnt(k)
c                     drcpnt(id,k) = drcpnt(id,k) - dum2
c                     drcpnt(ia,k) = drcpnt(ia,k) - fac1*dum2
c                     drcpnt(ih,k) = drcpnt(ih,k) - fac2*dum2
c                     drcpnt(ico3,k) = drcpnt(ico3,k) - fac3*dum2
c                     drcpnt(inh4,k) = drcpnt(inh4,k) - fac4*dum2
                  endif
               enddo
            else
c               drcpnt(id,id)= drcpnt(id,id) - deriv
c               drcpnt(ia,id) = drcpnt(ia,id) - fac1*deriv
c               drcpnt(ih,id) = drcpnt(ih,id) - fac2*deriv
c               drcpnt(ico3,id) = drcpnt(ico3,id) - fac3*deriv
c               drcpnt(inh4,id) = drcpnt(inh4,id) - fac4*deriv
            endif
         enddo

         deriv=qm1*newimmib*substrc*cka1
     &        /(ckc1+substrc)/(cka1+newtaqia)**2
         drtaq(ia) = -fac1*deriv

         deriv=qm1*substrc*newtaqia
     &        /(ckc1+substrc)/(cka1+newtaqia)
         drimm(ib) = yield1*deriv - decay1

      endif

      return
      end
c **********************************************************************
c <<< Subroutine subrxn5 >>>

c     User-defined subroutine #5 (Nodal reaction)
      subroutine subrxn5(dt,in,irxn)
c ----Inherent statements
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn,mi

      real*8 dt, stemp

c  ...Inherent common blocks

c ----Users can write codes below----

c     Type of reaction:
c        Radioactive Decay
      

c  ...User-introduced variables
      integer par,dau  
      real*8 newtaqpar
      real*8 newimmpar
      real*8 newvappar
      real*8 h_correct
c     New variables for enhancement to include equilibrium Kd
c     with decay
      real*8 convfact, sorbconc, conc_subst, betam1l, betam1v
      if(naqsp(irxn).ne.0)then
         par = irxnic(irxn,1)
         if (naqsp(irxn) .eq. 1) then
            dau = 0
         else
            dau = irxnic(irxn,2)
         endif
         mi = in+(pcpnt(par)-1)*n0

c     Routine now handles decay when there is sorption
c     through the code's original Kd approach
c     convfact is the equivalent aqueous concentration for
c     the portion that is sorbed to the rock
c     The total amount of the species now decays, and the 
c     portion of the derivative associated with the sorbed part
c     is included in drcpnt

cHari 3/26/08
         if (irdof .ne. 13) then
            stemp = min(strac_max,s(in))
            convfact = rolf(in)*stemp*ps_trac(in)
         else
            convfact = rolf(in)*ps_trac(in)
         end if
         scl(mi) = 0.
         scv(mi) = 0.
         dsccl(mi) = 0.
         dsccv(mi) = 0.
         h_correct = 1.0 
         
c     Update scl and dsccl if the species is sorbing
c     Because thermc and solstore have not yet been called
c     for this iteration

CPS   IF species is sorbing        
         if(iadsfl(pcpnt(par),itrc(mi)).ne.0)then
CPS   calculate sorption terms
            conc_subst = max(an(mi),1.d-90)
            scl(mi) = (denr(in) * a1adfl(pcpnt(par),itrc(mi)) *
     2           conc_subst
     2           **betadfl(pcpnt(par),itrc(mi)) / (1.0 +
     3           a2adfl(pcpnt(par),itrc(mi)) *
     4           conc_subst**betadfl(pcpnt(par),itrc(mi)) ))
            betam1l = betadfl(pcpnt(par),itrc(mi)) - 1
            conc_subst = max(an(mi),1.d-90)
            dsccl(mi) = denr(in) * a1adfl(pcpnt(par),itrc(mi)) *
     2           betadfl(pcpnt(par),itrc(mi)) *
     3           conc_subst**betam1l / ( 1.0 + 
     2           a2adfl(pcpnt(par),itrc(mi)) *
     4           conc_subst**betadfl(pcpnt(par),itrc(mi)))**2
         endif
c======================================================
c= = = = ADDING Henrys Fix to decay in vapor component
c= = = = 11/21/2003  P. Stauffer
c======================================================
CPS   IF species is Henrys

cHari 3/26/08
 
         if (abs(icns(par)).eq.2) then
           stemp = min(strac_max,s(in))               
           h_correct = ( anv(mi)*rovf(in)*(1-stemp) 
     2                  +anl(mi)*rolf(in)*(stemp) )
     3                  /(anl(mi)*rolf(in)*(stemp))

         end if

         sorbconc = scl(mi)/max(1.d-20,convfact)
         newtaqpar=an(mi)*h_correct + sorbconc
         rrcpnt(par)= -ckmtrn(irxn)*newtaqpar
         if(dau.ne.0) rrcpnt(dau)=  ckmtrn(irxn)*newtaqpar
         drcpnt(par,par)= -ckmtrn(irxn)*(1.+dsccl(mi)/
     2        max(1.d-20,convfact))
      endif
      if(nimsp(irxn).ne.0)then
         par = irxnim(irxn,1)
         if (nimsp(irxn) .eq. 1) then
            dau = 0
         else
            dau = irxnim(irxn,2)
         endif
         mi = in+(pimm(par)-1)*n0
         newimmpar=an(mi)
         rrimm(par)= -ckmtrn(irxn)*newimmpar
         if(dau.ne.0)rrimm(dau)=  -ckmtrn(irxn)*newimmpar            
         drimm(par)= -ckmtrn(irxn)
      endif
      if(nivsp(irxn).ne.0)then
         par = irxniv(irxn,1)
         if (nivsp(irxn) .eq. 1) then
            dau = 0
         else
            dau = irxniv(irxn,2)
         endif
         mi = in+(pvap(par)-1)*n0

c     Decay is corrected for sorbed portion, see above for details
CPS   IF species is sorbing        
         if(iadsfv(pvap(par),itrc(mi)).ne.0)then
CPS   calculate sorption terms
            conc_subst = max(an(mi),1.d-90)
            scv(mi) = (denr(in) * a1adfv(pvap(par),itrc(mi)) * 
     2           conc_subst
     2           **betadfv(pvap(par),itrc(mi)) / (1.0 + 
     2           a2adfv(pvap(par),itrc(mi))
     3           * conc_subst**betadfv(pvap(par),itrc(mi)) ))
            betam1v = betadfv(pvap(par),itrc(mi)) - 1
            conc_subst = max(an(mi),1.d-90)
            dsccv(mi) = denr(in)*a1adfv(pvap(par),itrc(mi))*
     2           betadfv(pvap(par),itrc(mi))*
     3           conc_subst**betam1v/ ( 1.0 +
     4           a2adfv(pvap(par),itrc(mi)) * 
     5           conc_subst**betadfl(pvap(par),itrc(mi)))**2
CPS   ENDIF
         endif

cHari 3/26/08
         stemp = min(strac_max,s(in))
         convfact = rovf(in)*(1.-stemp)*ps_trac(in)
         sorbconc = scv(mi)/max(1.d-20,convfact)
         newvappar=an(mi) + sorbconc
         rrvap(par)= -ckmtrn(irxn)*newvappar
         if(dau.ne.0)rrvap(dau)=  ckmtrn(irxn)*newvappar            
         drvap(par)= -ckmtrn(irxn)*(1.+dsccv(mi)/max(1.d-20,convfact))
      endif
      return
      end
c *********************************************************************
c <<< Subroutine subrxn6 >>>

c     User-defined subroutine #6 (Nodal reaction)

      subroutine subrxn6(dt,in,irxn)
      use comchem
      use comrxni
      use comci
      use comdi
      use comdti
      implicit none

c  ...Declarations of variables
      integer in,irxn

      real*8 dt, stemp
c  ...Inherent common blocks
c local variables
      integer ic,iv,mi
      real*8 ovap,ckmtrn1,ckeqlb1,conv
c ----Users can write codes below----

c     Type of reaction:

c        --- Henry's law reaction  ---
      ic=irxnic(irxn,1)
      iv=irxniv(irxn,1)
      mi=in+(pvap(iv)-1)*n0
      ovap=anlo(mi)
      ckeqlb1=ckeqlb(irxn)
      ckmtrn1=ckmtrn(irxn)

cHari 3/26/08
      stemp = min(strac_max,s(in))
      conv = (rovf(in)*(1-stemp))/(mw_air*stemp*rolf(in)*phi(in)*
     2     (1000/101.325))
      if(ic.gt.100)then
         rrcplx(ic) = -((cplx(ic)-ckeqlb1*ovap)/(1/ckmtrn1+(dt/conv)*
     2        ckeqlb1))
         drcplx(ic) = -1/(1/ckmtrn1+dt*ckeqlb1/conv)
         rrvap(iv) = -rrcplx(ic)/conv
         drvap(iv) = 0.0
      else
         rrcpnt(ic) = -((cpnt(ic)-ckeqlb1*ovap)/(1/ckmtrn1+(dt/conv)*
     2        ckeqlb1))
         drcpnt(ic,ic) = -1/(1/ckmtrn1+dt*ckeqlb1/conv)
         rrvap(iv) = -rrcpnt(ic)/conv
         drvap(iv) = 0.0
      endif
      return
      end


c *********************************************************************
c <<< Subroutine subrxn7 >>>

c     User-defined subroutine #7 (Nodal reaction)

      subroutine subrxn7(dt,in,irxn)

      use comchem
      use comrxni
      use combi
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn

      real*8 dt, stemp
c  ...Inherent common blocks
c local variables
      integer i,ic,mi,im,ic2,j
      real*8 aqconc(3)
      real*8 newimm
      real*8 qequil
      real*8 danl_subst
      real*8 mol_kgh2ohr_conv, mol_kgshr_conv
      real*8 den,por,proddr,ckeqlb1,eqcheck
c ----Users can write codes below----
      if (irdof .ne. 13) then

cHari 3/26/08
         stemp = min(strac_max,s(in))
         danl_subst = max(rolf(in)*stemp,rtol)
      else
         danl_subst = max(rolf(in),rtol)
      end if
      por = ps_trac(in)
      den = denr(in)
      if(temp_model_kin(irxn).ne.'l')then
         ckeqlb1=ckeqlb(irxn)
      else
         ckeqlb1=tcoeff(irxn,1)+t(in)*tcoeff(irxn,2)+t(in)**2*
     2        tcoeff(irxn,3)
         ckeqlb1=10**ckeqlb1
      endif


c     Type of reaction:

c        --- Precipitation/dissolution reaction for the
c            total aqueous concentration ---
      Qequil = 1
      do i=1,naqsp(irxn)
         ic=irxnic(irxn,i)
         mi=in+(pcpnt(ic)-1)*n0
         aqconc(i)=an(mi)
c  '**pdstic(irxn,i)' added 30 dec 99 -- geh
         Qequil = Qequil*aqconc(i)**pdstic(irxn,i)
      enddo
      im=irxnim(irxn,1)
      mi=in+(pimm(im)-1)*n0
      newimm=an(mi)
      if(pd_flag(im,in).eq.0.or.Qequil/ckeqlb1.gt.1)then
         pd_flag(im,in)=0
         do i = 1,naqsp(irxn)
            ic=irxnic(irxn,i)
            eqcheck = Qequil/ckeqlb1
            mol_kgh2ohr_conv = sx1(in)*ps_trac(in)*
     2           danl_subst/3600
c  multiplication by sx1(in) added 30 dec 99 -- geh
            rrcpnt(ic)=(sx1(in)*sarea(irxn)*pdstic(irxn,i)*
     2           ckmtrn(irxn)*(1-eqcheck))/mol_kgh2ohr_conv
            proddr = 1
            do j = 1,naqsp(irxn)
               ic2 = irxnic(irxn,j)
               if(ic.eq.ic2)then
                  proddr = proddr*pdstic(irxn,j)*aqconc(j)**
     2                 (pdstic(irxn,j)-1)
               else
                  proddr = proddr*aqconc(j)**pdstic(irxn,j)
               endif
            enddo
c  multiplication by sx1(in) added 30 dec 99 -- geh
            drcpnt(ic,ic)=(-sx1(in)*sarea(irxn)*pdstic(irxn,i)*
     2           proddr*ckmtrn(irxn)/ckeqlb1)/(mol_kgh2ohr_conv*
     3           dxct(ic,ic))
         enddo
         mol_kgshr_conv = sx1(in)*denr(in)/3600 
c  multiplication by sx1(in) added 30 dec 99 -- geh
         rrimm(im) = (-sx1(in)*sarea(irxn)*pdstim*ckmtrn(irxn)*
     2        (1-eqcheck))/mol_kgshr_conv
         drimm(im) = 0
      endif
      return
      end
c *********************************************************************
c <<< Subroutine subrxn8 >>>

c     User-defined subroutine #8 (Nodal reaction)

      subroutine subrxn8(dt,in,irxn)

      use comchem
      use comrxni
      use combi
      use comci
      use comdi
      use comdti
      use davidi, only : irdof
      implicit none

c  ...Declarations of variables
      integer in,irxn

      real*8 dt, stemp
c  ...Inherent common blocks
c local variables
      integer i,ic,mi,im,ic2,j
      real*8 aqconc(3)
      real*8 newimm
      real*8 qequil
      real*8 danl_subst
      real*8 mol_kgh2ohr_conv, mol_kgshr_conv
      real*8 den,por,prod,ckeqlb1,eqcheck
c ----Users can write codes below----
      if (irdof .ne. 13) then

cHari 3/26/08
         stemp = min(strac_max,s(in))
         danl_subst = max(rolf(in)*stemp,rtol)
      else
         danl_subst = max(rolf(in),rtol)
      end if
      por = ps_trac(in)
      den = denr(in)
      if(temp_model_kin(irxn).ne.'l')then
         ckeqlb1=ckeqlb(irxn)
      else
         ckeqlb1=tcoeff(irxn,1)+t(in)*tcoeff(irxn,2)+t(in)**2*
     2        tcoeff(irxn,3)
         ckeqlb1=10**ckeqlb1
      endif


c     Type of reaction:

c        --- Precipitation/dissolution reaction for complexes
c            and free ion concentrations ---

c Code altered Dec 30, 1999 by GEH

      Qequil = 0.d0
      do i=1,naqsp(irxn)
         ic=irxnic(irxn,i)
         if(ic.gt.100)then
            Qequil = Qequil + pdstic(irxn,i)*log(cplx(ic))
         else
            Qequil = Qequil + pdstic(irxn,i)*log(cpnt(ic))
         endif
      enddo
    
      Qequil = exp(Qequil)

c Compute Q/K
      eqcheck = Qequil/ckeqlb1  
c End of coding change    

      im=irxnim(irxn,1)
      mi=in+(pimm(im)-1)*n0
      newimm=an(mi)

      if(pd_flag(im,in).eq.0.or.eqcheck.gt.1)then
         
         mol_kgh2ohr_conv = sx1(in)*ps_trac(in)*danl_subst/3600

         do i = 1,naqsp(irxn)
            ic=irxnic(irxn,i)
c  results in mols/kgh20/hr
            prod = sx1(in)*sarea(irxn)*pdstic(irxn,i)*
     2             ckmtrn(irxn)/mol_kgh2ohr_conv
            if(ic.gt.100)then
               rrcplx(ic) = prod * (1 - eqcheck)
               drcplx(ic) = -prod*Qequil/ckeqlb1/cplx(ic)*
     2                       pdstic(irxn,i)
            else
               rrcpnt(ic) = prod * (1-eqcheck)

               do j = 1,naqsp(irxn)
                  ic2 = irxnic(irxn,j)
                  if(ic2.gt.100)then
                     drcpnt(ic,ic2)=-prod*Qequil/ckeqlb1/
     2                               cplx(ic2)*pdstic(irxn,j)
                  else
                     drcpnt(ic,ic2)=-prod*Qequil/ckeqlb1/
     2                               cpnt(ic2)*pdstic(irxn,j)
                  endif
               enddo   
            endif
         enddo

         mol_kgshr_conv = sx1(in)*denr(in)/3600 

         prod = sx1(in)*sarea(irxn)*pdstim*ckmtrn(irxn)/
     2          mol_kgshr_conv

         rrimm(im) = -prod*(1-eqcheck)
        
         drimm(im) = 0.0

      endif
      return
      end
c *********************************************************************
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) REAL*8 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END

      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) REAL*8 array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END

      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 20, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or real*8.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END

      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END

c
c
c-------solver---
c
c 
c
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) REAL*8 array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) REAL*8 array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END


      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END

c==========================================================================
c   BLAS SUBROUTINES 
c==========================================================================
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end

************************************************************************
*
*     File of the REAL*8 Level-3 BLAS.
*     ==========================================
*
*     SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
*    $                   BETA, C, LDC )
*
*     SUBROUTINE DSYMM ( SIDE,   UPLO,   M, N,    ALPHA, A, LDA, B, LDB,
*    $                   BETA, C, LDC )
*
*     SUBROUTINE DSYRK ( UPLO,   TRANS,     N, K, ALPHA, A, LDA,
*    $                   BETA, C, LDC )
*
*     SUBROUTINE DSYR2K( UPLO,   TRANS,     N, K, ALPHA, A, LDA, B, LDB,
*    $                   BETA, C, LDC )
*
*     SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
*    $                   B, LDB )
*
*     SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
*    $                   B, LDB )
*
*     See:
*
*        Dongarra J. J.,   Du Croz J. J.,   Duff I.  and   Hammarling S.
*        A set of  Level 3  Basic Linear Algebra Subprograms.  Technical
*        Memorandum No.88 (Revision 1), Mathematics and Computer Science
*        Division,  Argonne National Laboratory, 9700 South Cass Avenue,
*        Argonne, Illinois 60439.
*
*
************************************************************************
*
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL*8   ALPHA, BETA
*     .. Array Arguments ..
      REAL*8   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL*8.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL*8 array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL*8 array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL*8.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL*8 array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL*8   TEMP
*     .. Parameters ..
      REAL*8   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END

c***********************************************************************
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      real*8 da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL*8   ALPHA
*     .. Array Arguments ..
      REAL*8   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL*8.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - REAL*8 array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL*8 array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL*8   TEMP
*     .. Parameters ..
      REAL*8   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END

c***********************************************************************
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      real*8 dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end

************************************************************************
*
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      REAL*8   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      REAL*8   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL*8.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL*8 array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL*8 array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL*8 array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL*8   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      REAL*8   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER
*
      END

c *********************************************************************
c *********************************************************************
c==================================================================
c   LAPACK SUBROUTINES
c==================================================================
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL*8   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) REAL*8 array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL*8   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END






