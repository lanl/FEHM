      subroutine cappr(iz,ndummy) 
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
CD1 Calculate capillary pressure.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 05-20-92     G. Zyvoloski   00022   Initial implementation.
CD2                                        However, previous non-YMP
CD2                                        versions of FEHM exist, and
CD2                                        the current version differs
CD2                                        from these in minor ways.  
CD2
CD2 $Log:   /pvcs.config/fehm90/src/cappr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:58   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:04   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.6   Wed Feb 14 11:35:26 1996   zvd
CD2 Updated purpose to reflect input now done with rlp
CD2 
CD2    Rev 1.5   Mon Jan 29 13:16:58 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.4   12/13/95 08:30:04   gaz
CD2 cap pressure model inputted in macro rlp
CD2 
CD2    Rev 1.3   01/23/95 10:58:26   robinson
CD2 BAR: Fixed bug in computation of linear cap. pressure model
CD2 
CD2    Rev 1.2   09/02/94 12:14:28   robinson
CD2 Fixed bug in reading of parameters - BAR
CD2 
CD2    Rev 1.1   03/18/94 16:11:24   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:21:32   pvcs
CD2 original version in process of being certified
CD2 
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3 Identifier              Type     Use  Description
CD3
CD3 iz                      int       I   Control parameter used to
CD3                                          determine the reason for
CD3                                          this call
CD3 ndummy                  int       I   Parameter to set the correct
CD3                                          node number for dual
CD3                                          porosity or dpdp
CD3
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name           Description
CD3
CD3 input file     Contains all input data from FEHMN macros in ASCII
CD3                form
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
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
CD4 icapp, inpt, wdd1, icapt, cp1f, cp2f, cp3f, narrays, pointer, itype,
CD4 default, igroup, ireturn,
CD4 macroread
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4 
CD4 
CD4 
CD4 Global Subprograms
CD4 
CD4 Name      Type     Description
CD4 
CD4 initdata  N/A      Reads and sets input values defined at each node
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
CD5 Local Variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Do loop parameter
CD5 mid          int         Do loop index
CD5 mi           int         Current node number
CD5 it           int         Current region number of capillary
CD5                             pressure model
CD5 itp          int         Current capillary pressure model number
CD5 itperm       int         Current model number for relative
CD5                             permeability model
CD5 itpperm      int         Current relative permeability model number
CD5 cp1          real*8      Capillary pressure model parameter
CD5 cp3          real*8      Capillary pressure model parameter
CD5 im           int         Do loop parameter
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
CD9 2.4.4 Relative-permeability and capillary-pressure functions
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
CPS BEGIN cappr
CPS 
CPS IF capillary pressures are used in this simulation
CPS 
CPS   IF the data are to be read in on this pass
CPS 
CPS     LOOP to read in first line of data
CPS     
CPS       null - Check for null line or zero
CPS       IF there is more data to read
CPS         Read in capillary pressure parameters
CPS         IF the capillary pressure parameters were set correctly
CPS           Compute slope of linear capillary pressure model
CPS         ELSE the parameters are incorrect
CPS           Set slope to 0 to zero out capillary pressure at that node
CPS         ENDIF
CPS       ELSE
CPS         IF no data was present
CPS           Set parameter to neglect capillary pressures
CPS         ENDIF
CPS     EXITIF no more data in this group
CPS       ENDIF
CPS       
CPS     ENDLOOP
CPS     
CPS     Set parameters for reading capillary pressure flag
CPS     
CPS     initdata - read flag and set at each node
CPS     
CPS     Set flag to denote that the cap macro has been called
CPS       
CPS   ELSE parameter values are being set
CPS   
CPS     FOR each node
CPS       Set integer parameters
CPS       IF a van Ganuchten relative permeability model is used
CPS         Set parameter to skip capillary pressure calculation
CPS       ELSE
CPS         Set capillary pressure model flag
CPS       ENDIF
CPS       IF a linear model is assumed
CPS         Compute capillary pressure and derivative
CPS         IF the computed value is less than 0
CPS           Set value and derivative to 0
CPS         ENDIF
CPS       ENDIF
CPS     ENDFOR
CPS   
CPS     FOR each node
CPS       IF this is a compressed liquid
CPS         Set capillary pressure and derivative to 0
CPS       ENDIF
CPS     ENDFOR
CPS     
CPS   ENDIF the data are to be read in on this pass
CPS   
CPS ENDIF capillary pressures are used in this simulation
CPS 
CPS END cappr
CPS 
C**********************************************************************
c
c subroutine to calculate the capillary pressure
c
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comdti
      use comai
      use comki
      use comwt

      implicit none

      integer i,iz,ndummy,mid,it,mi,itp,itperm,itpperm,im, iz_temp
	real*8 hp, sl, dhp,rp1,rp2,rp3,rp4,sucut,slcut
      real*8 cp1,cp3,cap_max
      integer jj, ireg
	real*8 s_corr,ds_corrp,s_vg,ds_vgp
	real*8 ph_dum,hmin,hmax,pl,del_h,del_s,dr_nr
	real*8 ac3, ac4, bc3, bc4, hlcut, hucut, dslh
	real*8 smcut,smcutm,scutm,alamda,alpi,hcut,alpha,beta
c PJJohnson code changes - variable definitions
      real*8 pjni, pjSri, pjCpmaxi, pjn, pjlmax, pjsat, pjm, pjb
      real*8 pjSr, pjCpmax, pjconstant
c end PJJohnson variable changes      
	parameter (scutm = 1.d-03)	
	
c read in data
      if(icapp.eq.0) return
         if(iz.gt.0) then
c     
c     load nodal capillary pressures
c     
            do mid=1,neq
               mi=mid+ndummy
                  it=icap(mi)
               itperm = irlp(mi)
               itpperm = irlpt(itperm)
              if(it.eq.0.or.(itpperm .ge. 3 .and. itpperm .le. 20)) then
                  itp = itpperm 
              else
                  itp=icapt(it)
              end if
               if(itp.eq.1) then
c     
c     linear forsythe(1988) model
c     
c * * * * * * * * * * PJJOHNSON CODE MODIFICATIONS, SUMMER 2017 * * * * * * *  
c added four flags for por-dependent capillary functions
c these are tied to the linear rlp model - user assigns as normal except Slmax is given a flag, -111, -333, -666, or -999
c which activates the f(sat, n) one
c pjn, pjni, pjsr, pjSri, pjcpmax, pjcpmaxi are variables used for that
c ALL PJJOHNSON variables used for calcs throughout FEHM are noted by prefix pj                   
c                                  

                   if(irlpt(it).eq.-666) then
c set residual saturation and max capillary pressure
c                      open(666,file='Cp_out.txt')

                       pjSri= rp1f(it)
                       pjCpmaxi=cp1f(it)
                       
c read porosity from global variables                       
                       pjn=ps(mi)
                      if(pjn.ge. 0.9 .or. pjn.le.1e-6) then
                          pcp(mi) = 0
                          dpcef(mi) = 0
                      else
                          
c flag for initial porosity: if 0, then read it from global
c otherwise, user can specify value to use (e.g. for same material with different properties)
                      if(rp5f(it) .gt. 0) then
                          pjni = rp7f(it)
                      else
                      
                       pjni = psini(mi)
                      endif
                      
c flag for changing saturation above which Pc = 0; if 0, change it, otherwise stays same                       
                      if(cp3f(it).ne.0) then
                          pjlmax = 1
                      else
                          pjlmax = 1-pjn
                      endif
c linear residual saturation calculation                      
                       pjSr=(pjSri/(1-pjni))*(1-pjn)
                       
c extrapolate max capillary pressure - note that this is at sat = Sr, not sat = 0 (different from rlp 1)                       
c                       pjCpmax=(pjCpmaxi/pjsri)*(pjsr)
                       pjCpmax=pjCpmaxi*(1-pjn)/(1-pjni)
c read saturation                       
                       pjsat=s(mi)
c                         pjb = pjCpmax + (pjCpmax/(pjlmax-pjSr)*pjSr)
c calculate Pc
                       if(pjsat.ge.pjlmax) then
                           pcp(mi)=0
                       else if(pjsat.le.pjSr) then
                          pcp(mi)=pjCpmax
                       else
                         pjm = -1*(pjCpmax / (pjlmax - pjSr))
c define derivative                         
                         pcp(mi) = pjCpmax*(pjlmax-pjsat)/(pjlmax-pjSr)
c
c      
                          dpcef(mi) = pjm
c                          
                       end if
                       end if


                  else if(irlpt(it).eq. -333) then
c Leverett function   

c if user enters a number for porosity entry in rlperm.f, use that
c if they enter 0, then use the global variable
c this allows them to have multiple different units be compared to same standard
                      if(rp7f(it) .gt. 0) then
                          pjni = rp7f(it)
                      else
                      
                       pjni = psini(mi)

                      endif
c take initial maximum capillary pressure and residual saturation from input                       
                       pjCpmaxi = cp1f(it)

                       
                       pjSri= rp1f(it)
c
c read current porosity from global variables                       
                       pjn=ps(mi)
c Permeability comes from salt macro, but is sometimes fed 0 (e.g. timestep 1)
c so address that if need be                       
                       if(pjk(mi) .eq. 0) then
                           pjk(mi) = pjki(it)
                       endif

c Calculate leverett capillary pressure                       
                       pjCpmax = pjCpmaxi*(sqrt(pjki(it)/pjni) /
     &                 (sqrt(pjk(mi)/pjn)))
                       
                          pjlmax = 1

c handle residual saturation based on user preference
c if they use 0, hold Sr constant; otherwise, vary with porosity
                      if(cp3f(it).ne.0) then
                       pjSr=(pjSri/(1-pjni))*(1-pjn)
                      
                      else
                          pjSr = pjSri
                      endif
c truncate Sr in case it tries to go >1 or to .99999 with function divergence issues                      
                       if(pjSr .ge. 0.99) then
                           pjSr = 0.99
                       endif
c read saturation                       
                       pjsat=s(mi)
c calculate Pc
c max sat is currently 1 but if modified, this line handles that                       
                       if(pjsat.ge.pjlmax) then
                           pcp(mi)=0
                           dpcef(mi)=0
c if Sl < Sr, use maximum value                           
                       else if(pjsat.le.pjSr) then
                          pcp(mi)=pjCpmax
                          dpcef(mi)=0
c otherwise fit saturation function                          
                       else
                         pjm = -1*(pjCpmax / (pjlmax - pjSr))
                         pcp(mi) = pjCpmax*(pjlmax-pjsat)/(pjlmax-pjSr)
                         
c I added a truncation here to hold Pc to no more than an order of magnitude higher than the starting value
c to avoid crashes                         
                         if(pcp(mi) .gt. (10*pjCpmaxi)) then
                             pcp(mi) = 10*pjCpmaxi
c failsafe in case of model crash                             
                         else if (pcp(mi).lt.0) then
                             pcp(mi) = 0
                         else
                             pcp(mi) = pcp(mi)
                         end if

                         pjm = -1*(pjCpmax / (pjlmax - pjSr))
                         dpcef(mi) = pjm
                       end if
c following comment lines are output files that can be activated to check calcs                       
c                      write(333,*) pjni, pjn
c                      write(666,*) pjki(it), pjk(mi)
c                      write(999,*) pjCpmaxi, pjCpmax
c                      write(111,*) pjSri, pjSr
c                      write(777,*) pjsat, pcp(mi)
c gaz debug 011219
c                       end if
                       

   

                  else
                  pjn = 0
                  pjsat = 0
c *********END PJJOHNSON CHANGES **********
c 
c gaz debug 011219                  
                  endif
                  cp1=cp1f(it)
                  cp3=cp3f(it)
                  if(ieos(mi).eq.2) then
                     pcp(mi)=cp1*(cp3-s(mi))
                     dpcef(mi)=-cp1
                     if(pcp(mi).lt.0.0) then
                        pcp(mi)=0.0
                        dpcef(mi)=0.0
                     endif
                   else if(ieos(mi).eq.3) then
                     pcp(mi)=cp1*cp3
                     dpcef(mi)=0.0  
                   else if(ieos(mi).eq.1) then
                     pcp(mi)=0.0        
                     dpcef(mi)=0.0  
                   endif

               end if
               if(itp.eq.3) then
c     
c     van Genucten (1980) model (calculated in rlperm)
c     
               endif
             continue
            enddo
 
         else if(iz.lt.0) then
c
c inverse models (031707 gaz)
c     
c     load nodal capillary pressures
c   
         

         do mid=1,neq
              mi=mid+ndummy
c inverse calculation valid only for wtsi nodes
              iz_temp = izone_free_nodes(mi)
	        if(iz_temp .gt. 1) then  
               it=icap(mi)
               itperm = irlp(mi)
               itpperm = irlpt(itperm)
               if(itpperm .ge. 3 .and. itpperm .le. 20) then
                  itp = itpperm 
               else
                  itp=icapt(it)
               end if    
	   if(itp.eq.1) then          
c     
c     linear forsyth(1988) model
c     inverse solved explicitly
c
                  cp1=cp1f(it)
                  cp3=cp3f(it)
c  "+" removes negative cap pressure
                  hmin = head12(mi,1)+cp1*cp3 
		        hmax = head12(mi,2)
	            del_h = hmax-hmin
	            pl = phi(mi)
	             if(irich.eq.0) then
c richard's equation  with wtsi correction
                    if(iz_temp.eq.2) then
                      rlxyf(mi) = (pl-(hmin-cp1*cp3))/
     &				 (del_h+cp1) + rlptol   
			        drlxyf(mi) = 1.d0/(del_h+cp1) 
				    s(mi) = rlxyf(mi) - rlptol    
				  else
				   	rlxyf(mi) = 0.0
					s(mi) = 0.0
				    drlxyf(mi) = 1.d0/(del_h+cp1) 
				  endif	 
	            else
c richard's equation only
                   call vgcap_inv_calc(2,mi)   
	            endif
 
          else if(itp.eq.3) then
c     
c     van Genucten (1980) model
c  
c s = cap_inv  explicit calc 
c hmin and hmax are "pressures" 
c 
c richard's equation  with wtsi correction
	      if(irich.eq.0) then
c get richards equation terms 
	        call vgcap_inv_calc(2,mi)   
c wtsi terms  (still needs work gaz-112007)
             	  hmin = head12(mi,1) 
		      hmax = head12(mi,2)
	          del_h = hmax-hmin
	          pl = phi(mi) 
                s_corr = (pl-hmin)/del_h  
                ds_corrp = 1.d0/del_h
	          if(s_corr.lt.0.0d0) then
	            s_corr = 0.0d0
	            ds_corrp = 0.0d0
	          else if(s_corr.gt.1.0d0) then
	           s_corr = 1.0d0
	           ds_corrp = 0.0d0
 	          endif
c
c  gaz dslh is not defined 
	           s_vg = s(mi) 
                 ds_vgp = drlxyf(mi)			   
c 			   	
c Now the combination
c			   
			   	rlxyf(mi) = s_corr*(1.d0-s_vg) + s_vg + rlptol
                  s(mi) = rlxyf(mi) - rlptol
                  drlxyf(mi) = ds_corrp*(1.d0-s_vg) 
     &				+ (1.d0-s_corr)*ds_vgp   
	       
	       else 
c
c richard's equation only (test solution)
c
                 call vgcap_inv_calc(2,mi)   	          	 			                          
             endif
        endif
c end itp block
        endif
c end wtsi block
       enddo   
      endif

     
      return
      
      


      end
