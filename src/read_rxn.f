      subroutine read_rxn
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
CD1 To read the data specified in the rxn macro.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 12-08-92     B. Robinson    22      Initial Implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/read_rxn.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:06   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:12   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:42   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Fri Feb 02 10:21:00 1996   hend
CD2 Added to Requirements Traceability
CD2 
CD2    Rev 1.7   Thu Feb 01 15:47:02 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   Wed Jan 17 14:13:20 1996   hend
CD2 Added use of parser for input lines to function on sgi
CD2
CD2    Rev 1.5   08/07/95 11:57:46   awolf
CD2 Hari's changes on read(inpt and commented otu terms gamma_check, rate_factor, and round_tol
CD2 
CD2    Rev 1.4   04/25/95 10:06:18   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.3   01/28/95 14:20:46   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.2   06/20/94 11:14:28   zvd
CD2 Added ierr unit number for error output.
CD2
CD2    Rev 1.1   03/18/94 16:15:52   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:10   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 None
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name              Use     Description
CD3 
CD3 fehm input file    I      The fehm input file containing all data
CD3                           needed for the run
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4 
CD4 
c added 8/1/94
CD4 ngroups, group
CD4 
CD4 Global Constants
CD4
CD4 Identifier           Type     Description
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4
CD4                         COMMON
CD4 Identifier      Type    Block      Description
CD4
CD4 
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4
CD4 Identifier      Type         Description
CD4 
CD4 null1            LOGICAL      Checks for data on the next line

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
CD%
CD5 Identifier       Type        Description
CD5
CD5                                 information is present
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
CD9 2.4.5 Adsorbing solutes
CD9 2.4.6 Multiple, interacting solutes
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9 
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHM User's manual and Robinson's memo
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN read_rxn
CPS 
CPS 
CPS END read_rxn
CPS 
C**********************************************************************

      use comai
      use combi
      use comchem
      use comco2, only : icarb,carbon_tracer
      use comcouple
      use comdi
      use comdti
      use comrxni
      implicit none

* local variables
      integer logkeq
      integer ic,im,iv,ix,ii,idum
      character*80 input_msg
      integer msg(11)
      integer nwds
      real*8 xmsg(11)
      integer imsg(11)
      character*32 cmsg(11)
      character*4 por_change
      character*10 co2_flag
      logical null1
      
c ----Users can write codes below----

c  ...User-defined variables and common blocks
 

      integer i,j 
      integer igrp
      integer itype
      integer junk
      integer a1,a2,a3
      integer num_temps
      real*8 lkeq(20),eqtemp(20)
c ----------------------------
c Chemical species information
c ----------------------------

c Read in coupling information
      read(inpt,*)
      read(inpt,*)ncplx,numrxn
      read(inpt,*)
      read(inpt,*)ngroups
      do igrp = 1, ngroups
         read(inpt,*)
      enddo

      co2_couple = 0
      read(inpt, '(a10)') co2_flag
      if (co2_flag .eq. 'co2_couple') then
         if (icarb .eq. 1) then
            co2_couple = 1
            if (iout .ne. 0) write(iout,*) 'CO2 coupling with trac ON'
            if (iptty .ne. 0) write(iptty,*) 'CO2 coupling with trac ON'
         else
            write (ierr, *) '**** WARNING **** CO2 coupling not active',
     &           ' for non-CO2 problem'
         end if
      else
         backspace inpt
      end if

c Identification numbers and names of components
      read(inpt,*)
      do ic=1,ncpnt
         read(inpt,*) junk,cpntnam(ic),ifxconc(ic),cpntprt(ic),
     2        cpntgs(ic)
     	species(ic)=cpntnam(ic)
 	 if(cpntnam(ic).eq.'H2CO3') then
 		carbon_tracer=ic
           if (iout .ne. 0) write(iout,*) 'coupling with species# ',ic
           if (iptty .ne. 0) write(iptty,*) 'coupling with species# ',ic
 	 endif
 	if(cpntnam(ic).eq.'H'.or.cpntnam(ic).eq.'H+') then
 		ph_species=ic
 	endif

      enddo
      
c Identification numbers and names of complexes
      read(inpt,*)
      do ix=101,ncplx+100
         read(inpt,*) junk,cplxnam(ix),cplxprt(ix)
      enddo

c Identification numbers and names of immobile species
      read(inpt,*)
      do im=1,nimm
         read(inpt,*) junk,immnam(im),immprt(im)
         species(ncpnt+im)=immnam(im)
      enddo
c Identification numbers and names of vapor species
      read(inpt,*)
      do iv=1,nvap
         read(inpt,*) junk,vapnam(iv),vapprt(iv)
      enddo
C Arrays indicating which components/complex
      ii=0
      do ic = 1,ncpnt
         if(cpntprt(ic).eq.0)then
            ii = ii+1
            cpntprt(ii)=ic
         endif
      enddo
      ncpntprt=ii
      ii=100
      do ix=101,ncplx+100
         if(cplxprt(ix).eq.0)then
            ii= ii+1
            cplxprt(ii)=ix
         endif
      enddo
      ncplxprt=ii-100
      ii=0
      do im=1,nimm
         if(immprt(im).eq.0)then
            ii=ii+1
            immprt(ii)=im
         endif
      enddo
      nimmprt=ii
      ii=0
      do iv=1,nvap
         if(vapprt(iv).eq.0)then
            ii=ii+1
            vapprt(ii)=iv
         endif
      enddo
      nvapprt=ii

c -------------------------------
c Computation control information
c -------------------------------
c Decide whether to skip nodes 
      read(inpt,*)
      read(inpt,*) iskip


c Stopping criteria for speciation
c ..For Newton-Raphson
c     Parse the line to detect if the flag for handling negative
c     concentrations is present
      read(inpt,*)
      read (inpt,'(a)') input_msg
      call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)

c     Assign rsdmax
      if(msg(1).eq.1) then
         rsdmax = imsg(1)
      else
         rsdmax = xmsg(1)
      end if

c     Assign neg_conc_flag
      if(nwds.gt.1) then
         neg_conc_flag = imsg(2)
      else
         neg_conc_flag = 0
      end if

c -------------------------------------------
c Speciation and reaction control information
c -------------------------------------------
c Stability constants
      read(inpt,*)
      read(inpt,*)
      If(ncplx.ge.1) read(inpt,*) logkeq
      read(inpt,*)
      do ix=101,ncplx+100
         read(inpt,'(a1)')temp_model(ix)(1:1)
         backspace inpt
         if(temp_model(ix).eq.'l')then
            read(inpt,'(a)')input_msg
            call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
            num_temps=imsg(2)
            read(inpt,*)(eqtemp(i),i=1,num_temps)
            read(inpt,*)(lkeq(i),i=1,num_temps)
            do i = 1,num_temps
               if(logkeq.eq.0)lkeq(i)=dlog10(lkeq(i))
            enddo
            call lstsq(num_temps,lkeq,eqtemp)
            do i = 1,3
               heq(ix,i)=lkeq(i)
            enddo
            ckeq(ix)=heq(ix,1)+heq(ix,2)*25.0+heq(ix,3)*25.0**2
            ckeq(ix)=10**ckeq(ix)
         else
            read(inpt,'(a)')input_msg
            call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
            if(msg(1).eq.3)then
               if(msg(2).eq.1)then
                  ckeq(ix)=imsg(2)
               else
                  ckeq(ix)=xmsg(2)
               endif
               if(msg(3).eq.1)then
                  heq(ix,1)=imsg(3)
               else
                  heq(ix,1)=xmsg(3)
               endif
            else
               if(msg(1).eq.1)then
                  ckeq(ix)=imsg(1)
               else
                  ckeq(ix)=xmsg(1)
               endif
               if(msg(2).eq.1)then
                  heq(ix,1)=imsg(2)
               else
                  heq(ix,1)=xmsg(2)
               endif
            endif
            If(logkeq.eq.1) ckeq(ix)=10.d0**(ckeq(ix))
         endif
      enddo
               

c Stoichiometric coefficients
      read(inpt,*)
      do ix=101,ncplx+100
         read(inpt,*) (spstoic(ix,ic), ic=1,ncpnt)
      enddo

c Determine which components have negative stoichiometries, and
c therefore, can have ligitimate negative total concentrations 
c - geh apr20,1999

      do ic=1,ncpnt
         neg_conc_possible(ic) = 0
      enddo

      do ix=101,ncplx+100
         do ic=1,ncpnt
           if (spstoic(ix,ic) < 0.) neg_conc_possible(ic) = 1
         enddo
      enddo

c  ...If no reaction involved, no need for reading chemical parameters
      if(numrxn.lt.1) goto 999
      do j = 1, numrxn
         read(inpt,*)
         read(inpt,*)idrxn(j)
         read(inpt,*)
c  ...Read in where the reaction takes place
 10      read(inpt,'(a80)') wdd1
         if (null1(wdd1)) goto 20
         read(wdd1,*) a1,a2,a3
         if (a1.lt.0) then
c we are dealing with zone format
            a1=abs(a1)
            do i=1,n0
               if (izonef(i).eq.a1) then
                  rxnon(j,i)=1
               endif
            enddo
         else
c nodes are in block format
            if(a1.eq.1.and.a2.eq.0)then
               do i = 1,n0
                  rxnon(j,i)=1
               enddo
            else
               do i=a1,a2,a3
                  rxnon(j,i)=1
               enddo
            endif
         endif
         goto 10 
******************************************************************
c If linear kinetic reaction then read in the appropriate data
 20      if(idrxn(j).eq.1)then
c only 1 liquid and 1 solid can participate in a linear kinetic reaction
            naqsp(j)=1
            nimsp(j)=1
            nivsp(j)=0
c specify which free ion conc/complex and solid component participate 
            read(inpt,*)
            read(inpt,*)irxnic(j,1),irxnim(j,1)
c read in the distribution coefficient for rate-limiting 
c adsorption (kg water/kg rock) (add del H)
            read(inpt,*)
            read(inpt,'(a1)')temp_model_kin(j)(1:1)
            backspace inpt
            if(temp_model_kin(j).ne.'l')then
               read(inpt,*)ckeqlb(j)
            else
               read(inpt,'(a)')input_msg
               call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
               num_temps=imsg(2)
               read(inpt,*)(eqtemp(i),i=1,num_temps)
               read(inpt,*)(lkeq(i),i=1,num_temps)
               do i = 1,num_temps
                  lkeq(i)=log10(lkeq(i))
               enddo
               call lstsq(num_temps,lkeq,eqtemp)
               do i = 1,3
                  tcoeff(j,i)=lkeq(i)
               enddo
               ckeqlb(j)=tcoeff(j,1)+tcoeff(j,2)*25+tcoeff(j,3)*25**2
            endif
c read in the mass transfer coefficient for rate-limiting 
c adsorption (1/hr)(add Ea)
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
*******************************************************************
c Elseif Langmuir kinetic reaction then read in the appropriate data
         elseif(idrxn(j).eq.2)then
c only 1 liquid and 1 solid can participate in a Langmuir kinetic reaction
            naqsp(j)=1
            nimsp(j)=1
            nivsp(j)=0
c specify which free ion conc/complex and solid component participate 
            read(inpt,*)
            read(inpt,*)irxnic(j,1),irxnim(j,1)
c read in the distribution coefficient for rate-limiting 
c adsorption (kg water/moles) (add del H)
            read(inpt,*)
            read(inpt,'(a1)')temp_model_kin(j)(1:1)
            backspace inpt
            if(temp_model_kin(j).ne.'l')then
               read(inpt,*)ckeqlb(j)
            else
               read(inpt,'(a)')input_msg
               call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
               num_temps=imsg(2)
               read(inpt,*)(eqtemp(i),i=1,num_temps)
               read(inpt,*)(lkeq(i),i=1,num_temps)
               do i = 1,num_temps
                  lkeq(i)=log10(lkeq(i))
               enddo
               call lstsq(num_temps,lkeq,eqtemp)
               do i = 1,3
                  tcoeff(j,i)=lkeq(i)
               enddo
               ckeqlb(j)=tcoeff(j,1)+tcoeff(j,2)*25+tcoeff(j,3)*25**2
            endif
c read in the mass transfer coefficient for rate-limiting 
c adsorption (1/hr)(add Ea)
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
c read in the maximum sorption capacity (moles/kg rock)
            read(inpt,*)
            read(inpt,*)simmmx(j)
******************************************************************
c Elseif Generic (Reversible or Irreversible) Reaction
         elseif(idrxn(j).eq.3)then
c read in the number of solid, liquid and vapor species for rxn
            read(inpt,*)
            read(inpt,*)nimsp(j), naqsp(j), nivsp(j)
c read Forward and Reverse rate constants for rxn 
            read(inpt,*)
            read(inpt,*)kfor(j),krev(j)
c read in stoichiometry
            if(nimsp(j).gt.0)then
               read(inpt,*)
               read(inpt,*)(irxnim(j,i),i=1,nimsp(j))
               read(inpt,*)
               read(inpt,*)(stimirrv(j,i),i=1,nimsp(j))
            endif
            if(naqsp(j).gt.0)then
               read(inpt,*)
               read(inpt,*)(irxnic(j,i),i=1,naqsp(j))
               read(inpt,*)
               read(inpt,*)(sticirrv(j,i),i=1,naqsp(j))
            endif
            if(nivsp(j).gt.0)then
               read(inpt,*)
               read(inpt,*)(irxniv(j,i),i=1,nivsp(j))
               read(inpt,*)
               read(inpt,*)(stivirrv(j,i),i=1,nivsp(j))
            endif            
******************************************************************
C ELSEIF biodegradation reaction then read in appropriate data
         elseif(idrxn(j).eq.4)then
            read(inpt,*)
            read(inpt,*)naqsp(j),nimsp(j)
            nivsp(j)=0
            read(inpt,*)
            read(inpt,*)(irxnic(j,i),i=1,naqsp(j))
            read(inpt,*)
            read(inpt,*)(irxnim(j,i),i=1,nimsp(j))
            read(inpt,*)
            read(inpt,*)ckc(j)
            read(inpt,*)
            read(inpt,*)cka(j)
            read(inpt,*)
            read(inpt,*)decay(j)
            read(inpt,*)
            read(inpt,*)biofac(j)
            read(inpt,*)
            read(inpt,*)hfac(j)
            read(inpt,*)
            read(inpt,*)carbfac(j)
            read(inpt,*)
            read(inpt,*)ammfac(j)
            read(inpt,*)
            read(inpt,*)phthresh(j)
            read(inpt,*)
            read(inpt,*)qm(j)
            read(inpt,*)
            read(inpt,*)yield(j)
            read(inpt,*)
            read(inpt,*)xminit(j)
            read(inpt,*)
            read(inpt,*)nbiofrm(j)
            read(inpt,*)
            read(inpt,*)(icbio(j,i),i=1,nbiofrm(j))
******************************************************************
c If radioactive decay then read in the appropriate data
         elseif(idrxn(j).eq.5)then
c read in the half life (years)
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
            ckmtrn(j) = log(2.0)/(ckmtrn(j)*(365.25*24))
c read in reaction type (solid,liquid.gas) 
            read(inpt,*)
            read(inpt,*)itype
c if solid reaction
            if(itype.eq.0)then
               nimsp(j) = 2
c read in parent and daughter product for a solid species
               read(inpt,*)
               read(inpt,*) irxnim(j,1),idum
               if(idum.eq.0) then
                  nimsp(j)=1
               else
                  irxnim(j,2) = idum
               endif
	
            elseif(itype.eq.1)then
c read in parent and daughter product for a liquid species
               naqsp(j) = 2
               read(inpt,*)
               read(inpt,*) irxnic(j,1),idum
               if(idum.eq.0) then
                  naqsp(j)=1
               else
                  irxnic(j,2) = idum
               endif

            elseif(itype.eq.-1)then
c read in parent and daughter product for a vapor species
               nivsp(j) = 2
               read(inpt,*)
               read(inpt,*) irxniv(j,1),idum
               if(idum.eq.0) then
                  nivsp(j)=1
               else
                  irxniv(j,2) = idum
               endif
            endif
******************************************************************
c If Henry's law reaction then read in the appropriate data
         elseif(idrxn(j).eq.6)then
c only 1 liquid and 1 vapor can participate in a Henry's law reaction
            naqsp(j)=1
            nimsp(j)=0
            nivsp(j)=1
c specify which free ion conc/complex and vapor component participate
            read(inpt,*)
            read(inpt,*)irxnic(j,1),irxniv(j,1)
c read in the Henry's law constant (moles/atm*L)
            read(inpt,*)
            read(inpt,*)ckeqlb(j)
c read in the reaction rate constant
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
******************************************************************
c If precipiation/dissolution reaction then read in the appropriate data
         elseif(idrxn(j).eq.7)then
c only 1 immobile species can dissolve/precipitate
            nimsp(j)=1
c specify which immobile species participates in precip/diss
            read(inpt,*)
            read(inpt,*)irxnim(j,1)
c specify how many total aqueous species participate
            read(inpt,*)
            read(inpt,*)naqsp(j)
c specify which total aqueous components participate
            read(inpt,*)
            read(inpt,*)(irxnic(j,i),i=1,naqsp(j))
c specify the stoichiometry for the immobile components
            read(inpt,*)
            read(inpt,*)pdstim
c specify the stoichiometry for the aqueous components
            read(inpt,*)
            read(inpt,*)(pdstic(j,i),i=1,naqsp(j))
c read in the solubility product
            read(inpt,*)
            read(inpt,'(a1)')temp_model_kin(j)(1:1)
            backspace inpt
            if(temp_model_kin(j).ne.'l')then
               read(inpt,*)ckeqlb(j)
            else
               read(inpt,'(a)')input_msg
               call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
               num_temps=imsg(2)
               read(inpt,*)(eqtemp(i),i=1,num_temps)
               read(inpt,*)(lkeq(i),i=1,num_temps)
               do i = 1,num_temps
                  lkeq(i)=log10(lkeq(i))
               enddo
               call lstsq(num_temps,lkeq,eqtemp)
               do i = 1,3
                  tcoeff(j,i)=lkeq(i)
               enddo
               ckeqlb(j)=tcoeff(j,1)+tcoeff(j,2)*25+tcoeff(j,3)*25**2
            endif
c read in the rate of reaction in moles/(m^2*sec)
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
c read in the surface area of the mineral in m^2
            read(inpt,*)
            read(inpt,*)sarea(j)
******************************************************************
c If precipiation/dissolution reaction then read in the appropriate data
         elseif(idrxn(j).eq.8)then
c only 1 immobile species can dissolve/precipitate
            nimsp(j)=1
c specify which immobile species participates in precip/diss
            read(inpt,*)
            read(inpt,*)irxnim(j,1)
c specify how many total aqueous species participate
            read(inpt,*)
            read(inpt,*)naqsp(j)
c specify which free ion concs/aqueous complexes participate
            read(inpt,*)
            read(inpt,*)(irxnic(j,i),i=1,naqsp(j))
c specify the stoichiometry for the immobile components
            read(inpt,*)
            read(inpt,*)pdstim
c specify the stoichiometry for the free ion concs/aqueous complexes
            read(inpt,*)
            read(inpt,*)(pdstic(j,i),i=1,naqsp(j))
c read in the solubility product
            read(inpt,*)
            read(inpt,'(a1)')temp_model_kin(j)(1:1)
            backspace inpt
            if(temp_model_kin(j).ne.'l')then
               read(inpt,*)ckeqlb(j)
            else
               read(inpt,'(a)')input_msg
               call parse_string(input_msg,imsg,msg,xmsg,cmsg,nwds)
               num_temps=imsg(2)
               read(inpt,*)(eqtemp(i),i=1,num_temps)
               read(inpt,*)(lkeq(i),i=1,num_temps)
               do i = 1,num_temps
                  lkeq(i)=log10(lkeq(i))
               enddo
               call lstsq(num_temps,lkeq,eqtemp)
               do i = 1,3
                  tcoeff(j,i)=lkeq(i)
               enddo
               ckeqlb(j)=tcoeff(j,1)+tcoeff(j,2)*25+tcoeff(j,3)*25**2
            endif
c read in the molecular weight and density of the mineral for the optional 
c porosity change due to precip/diss reaction
	      read(inpt,'(a1)')por_change
		  if(por_change.eq.'p')then
	         read(inpt,*)
	         read(inpt,*) mw_mineral(irxnim(j,1)), 
     &                rho_mineral(irxnim(j,1))
	      else
	         backspace(inpt)
		  endif
c read in the rate of reaction in moles/(m^2*sec)
            read(inpt,*)
            read(inpt,*)ckmtrn(j)
c read in the surface area of the mineral in m^2
            read(inpt,*)
            read(inpt,*)sarea(j)
		  
         endif
      enddo
 999  return
      end

