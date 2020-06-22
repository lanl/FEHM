      subroutine vgcap_inv_calc(iflg,ndummy)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract with the U.S. Department of Energy (DOE). 
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
CD1 To compute the saturation from capillary pressure and derivatives for the 
CD1 van Genuchten model. The inverse of subroutine vg_cap.f .
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2 4-04-07      george Zyvoloski initial implementation
CD2
CD2 $Log:   /pvcs.config/fehm90/src/vgcap_inv.f_a  $                                     
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 sl           real*8   I/O    Liquid saturation on input,
CD3                                  normalized liquid saturation on
CD3                                  output
CD3 slr          real*8   I      Residual liquid saturation
CD3 smr          real*8   I      Maximum liquid saturation
CD3 alpha        real*8   I      Van Ganuchten parameter
CD3 beta         real*8   I      Van Ganuchten parameter
CD3 ac1          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac2          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 ac4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 smcut        real*8   I      Lower cutoff saturation below which a
CD3                                  spline fit is used for capillary
CD3                                  pressure versus saturation
CD3 sucut        real*8   I      Upper cutoff saturation above which a
CD3                                  linear fit is used for capillary
CD3                                  pressure versus saturation
CD3 bc3          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 bc4          real*8   I      Parameter in spline fit for capillary
CD3                                  pressure at low saturations
CD3 hp           real*8   O      Capillary pressure
CD3 dhp          real*8   O      Derivative of capillary pressure with
CD3                                  respect to saturation
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3 
CD3 NONE
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4 
CD4 NONE
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 NONE
CD4 
CD4 Global Subprograms
CD4 
CD4 NONE
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5 
CD5 NONE
CD5 
CD5 Local Types
CD5
CD5 NONE
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 alamda       real*8      Exponent in correlation
CD5 ds           real*8      Reciprocal of demoninator in normalized
CD5                              saturation expression
CD5 denom        real*8      Denominator in normalized saturation
CD5                              expression
CD5 star         real*8      Normalized saturation
CD5 alpi         real*8      Exponent in correlation
CD5 hp           real*8      Capillary pressure before unit conversion
CD5 dhp          real*8      Derivative of hp with respect to saturation
CD5 termstar1    real*8      Term used in capillary pressure calculation
CD5 termstar2    real*8      Term used in capillary pressure calculation
CD5 termb1       real*8      Term used in capillary pressure calculation
CD5 termb2       real*8      Term used in capillary pressure calculation
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 ASSUMPTIONS AND LIMITATIONS
CD6 
CD6 N/A
CD6
C**********************************************************************
CD7
CD7 SPECIAL COMMENTS
CD7 
CD7  Requirements from SDN: 10086-RD-2.20-00
CD7    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD7    FEHM Application Version 2.20
CD7
C**********************************************************************
CD8
CD8 REQUIREMENTS TRACEABILITY
CD8 
CD8 2.4.4 Relative-permeability and capillary-pressure functions
CD8
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See GZSOLVE SRS, MMS, and SDD for documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN vgcap_inv
CPS 
CPS 
CPS END vgcap_inv
CPS
C**********************************************************************
      use comhi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
	use comwt
	use comii
      use comki
      implicit none

      real*8 scutm,hmin,darcyf,tol_l,tol_u,su_cut
      parameter(scutm = 1.d-03)
      parameter(hmin = 1.d-18)
      parameter(darcyf = 1.d12)
      parameter(tol_l  = 1.d-5)
      parameter(tol_u  = 1.d-5)
      parameter(su_cut = 0.99d00)
 

      integer iz,ndummy,i,irlpd,mi,ieosd,it,ir,j,num_models,ireg
	integer iflg,inode, inr, inr_max, iz_temp
      real*8 alpha,beta,alamda,alpi,smcut,slcut,fac,ds,dhp,rl,rv
      real*8 drls,drvs,rp1,rp2,rp3,rp4,denom,star,hp,rl1,rv1
      real*8 drls1,drvs1,akf,akm,porf,permb,sl,sl_init
      real*8 smcutm,smcutf,alpham,alamdam,facf
      real*8 rpa1, rpa2, rpa3, rpa4, rpa5, cp1, cp3
	real*8 dresids,resid,dels,resid_max ,deriv_tol, approx_deriv
	real*8 ac3,ac4,bc3,bc4,hucut  
      real*8 sucut,hlcut,pl,dslh
	parameter (inr_max= 30, deriv_tol = 1.d-7)    
      logical null1,ex
      resid_max = max(0.0000001d0,rlptol)
      if(iflg.eq.1) then
C calculate the inverse (S) for VG cap models (all nodes)
 
         do i = 1,n
           mi = i+ndummy
	     iz_temp = izone_free_nodes(mi)
           if(iz_temp.gt.1) then
            pl = phi(mi)
            if (rlp_flag .eq. 0) then
            ieosd = ieos(mi)
               it = 0
            else
               it = irlp(mi)
            end if
            if(it.eq.0) then
               irpd=0
            else
               irpd = irlpt(it)
            endif
            if(irpd.ge.3.and.irpd.le.8) then 
               ieosd = 2
            endif 
	         if(irpd.eq.1) then
c     
c     linear forsyth(1988) model
c     inverse solved explicitly
c
                 if(iz_temp.eq.2) then
                  cp1=cp1f(it)
                  cp3=cp3f(it)
	             pl = pl - pref
				 rlxyf(mi) = (pl+cp1*cp3)/cp1 + rlptol  
			     drlxyf(mi) = 1.d0/cp1 
				 s(mi) = rlxyf(mi) - rlptol
			   else 
	            cp1=cp1f(it)
	             rlxyf(mi) = rlptol
				 s(mi) = 0.0
				 drlxyf(mi) = 1.d0/cp1 
			   endif	       
               else if(irpd.eq.3) then
c VG model
c richard's equation only
		      it = irlp(mi)               
	          rp1 = rp1f(it)
                rp2 = rp2f(it)
                rp3 = rp3f(it)
                rp4 = rp4f(it)
	          ac3 = rp9f(it)
	          ac4 = rp10f(it)   
	          bc3 = cp1f(it)
	          bc4 = cp2f(it)
c linear fit from slcut to sl=0.0
c first calculate star
                smcut=(rp6f(it)-rp1)/(rp2-rp1)
                smcutm = max(smcut,scutm)  
c now calculate sl
	          slcut = smcutm*(rp2-rp1) + rp1
	          alpha = rp3
	          beta = rp4
	          alamda = 1.0-1.0/beta	           
                alpi = 1.0/alamda
                hucut = 1.0/alpha*(1.0/smcutm**alpi-1.0)**(1.0-alamda)
c cutoff for near sucut = 1 region
                sucut = (0.99d0*(rp2-rp1) + rp1)
	          hlcut = -bc3*(1.-sucut)
	          pl = phi(mi)-pref
	          hp = pl/h_to_p      
	             call vgcap_inv(1, sl, rp1, rp2, rp3, rp4, 
     &              slcut,sucut,hlcut, hucut, ac3, ac4, bc3, bc4,
     &              hp, dslh)	
                 if(sl.lt.0.0d0) then
	            sl = 0.0d0
 
	          else if(sl.gt.1.d0) then
	           sl = 1.0d0
	           dslh = 0.0d0	         
			   endif            
	           rlxyf(mi) = sl + rlptol
                 drlxyf(mi) = -dslh/h_to_p
			   s(mi) = sl 	 
			 endif
			 if(iz_temp.eq.3) then
	          rlxyf(mi) =  rlptol
	          s(mi) = 0.0
			 endif	  
	  	else
	     sl = 1.0
           rlxyf(mi) = sl
           drlxyf(mi) = 0.0d0
		 s(mi) = sl 
		endif	 
        enddo	   
c end loop on nodes	  					       
      else if(iflg.eq.2) then   
c	
c single node calculation 
c

           mi = ndummy
           if(izone_free_nodes(mi).gt.1) then
            pl = phi(mi)
            if (rlp_flag .eq. 0) then
            ieosd = ieos(mi)
               it = 0
            else
               it = irlp(mi)
            end if
            if(it.eq.0) then
               irpd=0
            else
               irpd = irlpt(it)
            endif
            if(irpd.ge.3.and.irpd.le.8) then 
               ieosd = 2
            endif 
	         if(irpd.eq.1) then
c     
c     linear forsyth(1988) model
c     inverse solved explicitly
c
                   cp1=cp1f(it)
                   cp3=cp3f(it)
	             pl = pl - pref
				 rlxyf(mi) = (pl+cp1*cp3)/cp1 + rlptol  
			     drlxyf(mi) = 1.d0/cp1 
				 s(mi) = rlxyf(mi) - rlptol	       	       
               else if(irpd.eq.3) then
c VG model
c richard's equation only
		      it = irlp(mi)               
	          rp1 = rp1f(it)
                rp2 = rp2f(it)
                rp3 = rp3f(it)
                rp4 = rp4f(it)
	          ac3 = rp9f(it)
	          ac4 = rp10f(it)   
	          bc3 = cp1f(it)
	          bc4 = cp2f(it)
c linear fit from slcut to sl=0.0
c first calculate star
                smcut=(rp6f(it)-rp1)/(rp2-rp1)
                smcutm = max(smcut,scutm)  
c now calculate sl
	          slcut = smcutm*(rp2-rp1) + rp1
	          alpha = rp3
	          beta = rp4
	          alamda = 1.0-1.0/beta	           
                alpi = 1.0/alamda
                hucut = 1.0/alpha*(1.0/smcutm**alpi-1.0)**(1.0-alamda)
c cutoff for near sucut = 1 region
                sucut = (0.99d0*(rp2-rp1) + rp1)
	          hlcut = -bc3*(1.-sucut)
	          pl = phi(mi)-pref
	          hp = pl/h_to_p      
	             call vgcap_inv(1, sl, rp1, rp2, rp3, rp4, 
     &              slcut,sucut,hlcut, hucut, ac3, ac4, bc3, bc4,
     &              hp, dslh)	
                 if(sl.lt.0.0d0) then
	            sl = 0.0d0
c	            dslh = 0.0d0	 
	          else if(sl.gt.1.d0) then
	           sl = 1.0d0
	           dslh = 0.0d0	         
			   endif            
	           rlxyf(mi) = sl + rlptol
                 drlxyf(mi) = -dslh/h_to_p
			   s(mi) = sl 	 
			 endif	  
	  	else
	     sl = 1.0
           rlxyf(mi) = sl
           drlxyf(mi) = 0.0d0
		 s(mi) = sl 
		endif	 
      
	   
      endif     
      return
      end
