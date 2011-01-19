      subroutine stress_mech_props(iflg,model_flag,ndummy)
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
CD1 To calculate coefficients and derivatives for mech properties
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 01-20-2007     G. Zyvoloski   00022   Initial implementation.
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
CD9 2.3.2 Heat- and mass-transfer equations
CD9 2.3.3 Noncondensible gas flow equations
CD9 2.3.7 Sources and sinks
CD9 2.4.1 Pressure- and temperature-dependent water properties
CD9 2.4.2 Properties of air and air/water vapor mixtures
CD9       Stress and displacement calculations
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
CPS 
C**********************************************************************
c

c
      use davidi
      use comai
      use comii
      use comgi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comwt
      use comsi
      use comfem, only  : ifem
      implicit none
      
      integer ndummy,mid,mi,ieosd,kq,iflg, i, model_flag
      real*8 pl,tl,dtin
c s kelkar may 2010, tmeperature dependant properties       
      real*8 ddsdde(6,6)

      integer                      :: itmp, iModel
      real*8,  dimension(6)        :: gp_stress, gp_strain
      real*8,  dimension(6)        :: gp_strain_mech
      real*8,  dimension(6, 6)     :: DSai
      real*8 bulk_tol, bulk_mod

      parameter(bulk_tol=1.d-12)

c     
c     misc. constants
      
c     
      dtin=1.0/dtot
c     
c     generate mechanical properties as function of 
c     s kelkar, April 20, 2010    
      if(isNonlinear.eq.1) then
c     stiffness moduli as function of temperature only
         do mid=1,neq
            i=mid+ndummy 
            pl = phi(i) 
            tl = t(i)
            if(model_flag.eq.1) then
               elastic_mod(i) = e_ini(i) + dEdt(i)*(tl - tini(i))
               poisson(i) = poisson_ini(i) + dNuedt(i)*(tl - tini(i))
               if(istrs.ne.2) then
                  e1(i) = elastic_mod(i)*(1.0d0-poisson(i))/
     &                 (1.d0+poisson(i))/(1.0d0-2.0d0*poisson(i))
                  e2(i) = e1(i)*poisson(i)/(1.0d0-poisson(i))
                  e3(i) = e1(i)*(1.0d0-2.0d0*poisson(i))/
     &                 2.0d0/(1.0d0-poisson(i))
               else
                  e1(i) = elastic_mod(i)/(1.d0-poisson(i)*poisson(i))
                  e2(i) = e1(i)*poisson(i)
                  e3(i) = e1(i)*(1.0d0-poisson(i))/2.0d0
               endif
            elseif(model_flag.eq.91) then
c     s kelkar Oct 2010, table lookup 
               call young_temp_table(i)
            endif
            
         enddo      
      elseif(iPlastic.eq.1) then
         do mid=1,neq
            i=mid+ndummy 
            itmp = modelNumber(i)
            iModel = plasticModel(itmp)
            if(ifem.eq.1) then
c
            endif
         enddo
      endif 
c     
      return
      end
c..............................................................

      subroutine young_temp_table(i)
c table lookup for young's modulus as function of temperature

      use comai, only:istrs
      use comdi, only: t
      use comsi

      implicit none
      integer i,itable
      real*8 youngt, fact,tempi, poisst

      tempi = t(i)

      if(tempi.le.e_temp91(1,1)) then
         youngt=e_temp91(1,2)
         poisst=e_temp91(1,3)
      else
         do itable=2,nentries_young
            if(tempi.lt.e_temp91(itable,1)) then
               fact=(e_temp91(itable,2)-e_temp91(itable-1,2))
     &              /(e_temp91(itable,1)-e_temp91(itable-1,1))
               youngt=(tempi-e_temp91(itable-1,1))*fact
     &              +e_temp91(itable-1,2)
               fact=(e_temp91(itable,3)-e_temp91(itable-1,3))
     &              /(e_temp91(itable,1)-e_temp91(itable-1,1))
               poisst=(tempi-e_temp91(itable-1,1))*fact
     &              +e_temp91(itable-1,3)
               goto 9193
            endif
         enddo
         youngt=e_temp91(nentries_young,2)
         poisst=e_temp91(nentries_young,3)
         
 9193    continue
      endif

         write(97,*)i,youngt,poisst

      elastic_mod(i) = youngt
      poisson(i) = poisst

      if(istrs.ne.2) then
         e1(i) = youngt*(1.0d0-poisst)/
     &        (1.d0+poisst)/(1.0d0-2.0d0*poisst)
         e2(i) = e1(i)*poisst/(1.0d0-poisst)
         e3(i) = e1(i)*(1.0d0-2.0d0*poisst)/
     &        2.0d0/(1.0d0-poisst)
      else
         e1(i) = youngt/(1.d0-poisst*poisst)
         e2(i) = e1(i)*poisst
         e3(i) = e1(i)*(1.0d0-poisst)/2.0d0
      endif
      
      return

      end
c..............................................................
