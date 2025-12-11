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
      
      integer ndummy,mid,mi,ieosd,kq,iflg, i, model_flag, i_tab
      real*8 pl,tl,dtin
c s kelkar may 2010, tmeperature dependant properties       
      real*8 ddsdde(6,6)

      integer                      :: itmp, iModel
      real*8,  dimension(6)        :: gp_stress, gp_strain
      real*8 bulk_tol, bulk_mod

      parameter(bulk_tol=1.d-12)

c     
c     misc. constants
c 
c     
      dtin=1.0/dtot
c     
c     generate mechanical properties as function of 
c     s kelkar, April 20, 2010    
      if(isNonlinear.eq.1) then
c     stiffness moduli and poisson ratio as function of temperature only
         do mid=1,neq
            i=mid+ndummy 
            pl = phi(i) 
            tl = t(i)
            i_tab = iy_tab(i)
            if(istr_non_model(i_tab).eq.1) then
               elastic_mod(i) = e_ini(i_tab) + 
     &             dEdt(i_tab)*(tl - t_non_ref(i_tab))
               poisson(i) = poisson_ini(i_tab) + 
     &             dNuedt(i_tab)*(tl - t_non_ref(i_tab))
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
            else
c     s kelkar Oct 2010, table lookup 
               call young_temp_table(1,i,i_tab)
            endif            
         enddo   
       endif  
       if(isbiotNonLin.eq.1) then
c     stiffness thermal expansion and biot term as function of temperature only  
           do mid=1,neq
            i=mid+ndummy 
            pl = phi(i) 
            tl = t(i)
            i_tab = iy_tab_biot(i)
            if(istr_non_model_biot(i_tab).eq.1.or. 
     &          istr_non_model_biot(i_tab).eq.91) then
              call biot_temp_table(1,i,i_tab)
            endif
           enddo
       endif
       if(iPlastic.eq.1) then
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
c calculate elastic constants
c 
      call elastic_constants(1)
c      
      return
      end
c..............................................................

      subroutine young_temp_table(iflg,j,i)
c
c table lookup for young's modulus as function of temperature
c gaz 042216 modified for multiple tables
c iflg = 0, read input. = 1, evaluate table
c j = node number
c i = table number
c 
      use comai
      use comdi, only: t
      use comsi

      implicit none
      integer i,j,itable,ifile,iflg,i91,j91, idum
      integer open_file
      real*8 youngt, fact,tempi, poisst
      character*80 young_temp_file
      character*100 temp_junk
      if(iflg.eq.0) then
c read input

               read(inpt,*) idum, young_temp_file
               ifile = open_file( young_temp_file, 'old')
c read title
               read(ifile,'(a)') temp_junk
c check size of table
               do i91=1,nentries_young_max
                  read(ifile,*,end=91913) temp_junk
               enddo
               write(iptty,*)'error in stress input. Too many entries'
               write(iptty,*)' in the E vs Temperature file. STOP'
               write(ierr,*)'error in stress input. Too many entries'
               write(ierr,*)' in the E vs Temperature file. STOP'
               stop
91913          nentries_young = i91-1
c
c  this is where loop on tables should end
c
               rewind (ifile)
               read(ifile,'(a)') temp_junk
               do i91=1,nentries_young
                  read(ifile,*)(e_temp91(i91,j91,i),j91=1,3)
               enddo
             close (ifile)
             
      else
c evaluate table

      tempi = t(j)

      if(tempi.le.e_temp91(1,1,i)) then
         youngt=e_temp91(1,2,i)
         poisst=e_temp91(1,3,i)
      else
         do itable=2,nentries_young
            if(tempi.lt.e_temp91(itable,1,i)) then
               fact=(e_temp91(itable,2,i)-e_temp91(itable-1,2,i))
     &              /(e_temp91(itable,1,i)-e_temp91(itable-1,1,i))
               youngt=(tempi-e_temp91(itable-1,1,i))*fact
     &              +e_temp91(itable-1,2,i)
               fact=(e_temp91(itable,3,i)-e_temp91(itable-1,3,i))
     &              /(e_temp91(itable,1,i)-e_temp91(itable-1,1,i))
               poisst=(tempi-e_temp91(itable-1,1,i))*fact
     &              +e_temp91(itable-1,3,i)
               goto 9193
            endif
         enddo
         youngt=e_temp91(nentries_young,2,i)
         poisst=e_temp91(nentries_young,3,i)
        
 9193    continue
      endif
c gaz after testing, comment out next lines
c      write
c     & (97,'(t1,i6,t10,f10.3,t20,f12.4,t40,f12.4)')j,t(j),youngt,poisst

      elastic_mod(j) = youngt
      poisson(j) = poisst


      endif
      return

      end
c..............................................................

      subroutine biot_temp_table(iflg,j,i)
c
c table lookup for bulk() and alp() as function of temperature
c gaz  modified for multiple tables
c iflg = 0, read input. = 1, evaluate table
c j = node number
c i = table number
c 
      use comai
      use comdi, only: t
      use comsi

      implicit none
      integer i,j,itable,ifile,iflg,i91,j91, idum
      integer open_file
      real*8 bulkt, bulk_tol, bulk_mod, fact, tempi, alpt
      parameter(bulk_tol=1.d-12)
      character*80 biot_temp_file
      character*100 temp_junk
      if(iflg.eq.0) then
c read input (T,alp,bulk)

               read(inpt,*) idum, biot_temp_file
               ifile = open_file(biot_temp_file, 'old')
c read title
               read(ifile,'(a)') temp_junk
c check size of table
               do i91=1,nentries_biot_max
                  read(ifile,*,end=91913) temp_junk
               enddo
               write(iptty,*)'error in biot input. Too many entries'
               write(iptty,*)' in the biot vs Temperature file. STOP'
               write(ierr,*)'error in biot input. Too many entries'
               write(ierr,*)' in the biot vs Temperature file. STOP'
               stop
91913          nentries_biot = i91-1
c
c  this is where loop on tables should end
c
               rewind (ifile)
               read(ifile,'(a)') temp_junk
               do i91=1, nentries_biot
                  read(ifile,*)(biot_temp91(i91,j91,i),j91=1,3)
               enddo
             close (ifile)
             
      else
c evaluate table

      tempi = t(j)

      if(tempi.le.biot_temp91(1,1,i)) then
         alpt=biot_temp91(1,2,i)
         bulkt=biot_temp91(1,3,i)
      else
         do itable=2,nentries_biot
            if(tempi.lt.biot_temp91(itable,1,i)) then
               fact=(biot_temp91(itable,2,i)-biot_temp91(itable-1,2,i))
     &              /(biot_temp91(itable,1,i)-biot_temp91(itable-1,1,i))
               alpt=(tempi-biot_temp91(itable-1,1,i))*fact
     &              +biot_temp91(itable-1,2,i)
               fact=(biot_temp91(itable,3,i)-biot_temp91(itable-1,3,i))
     &              /(biot_temp91(itable,1,i)-biot_temp91(itable-1,1,i))
               bulkt=(tempi-biot_temp91(itable-1,1,i))*fact
     &              +biot_temp91(itable-1,3,i)
               goto 9193
            endif
         enddo
         alpt=biot_temp91(nentries_biot,2,i)
         bulkt=biot_temp91(nentries_biot,3,i)   
 9193    continue
      endif
      alp0(j) = alpt
      bulk0(j) = bulkt
      endif
      return

      end
c..............................................................
      subroutine elastic_constants(iflg)
c
c calculate elastic constants
c
      use comai
      use comdi, only: t
      use comdti
      use comsi

      implicit none
      integer i,j,itable,ifile,iflg,i91,j91, idum
      integer open_file
      real*8 bulkt, bulk_tol, bulk_mod, fact, tempi, alpt
      real*8 young_p , young_t, pois_p, pois_t, pois_sq
      real*8 fac1,fac2,fac3, ezzi,ezzkb,ezzbar,efacxy,efacz 
      parameter(bulk_tol=1.d-12)
      if(iflg.eq.1) then
c..............................................................
c     
c     linear isotropic or anisotropic (at present plain strain and 3D)
c     plain stress has different combinations
c     
         do i = 1,n0
c     change from volumetric to linear coef. of thermal expansion
c             
            alp(i) = alp0(i)/3.0
            if(istrs.ne.2) then
c     plain strain and 3-D                   
c     s kelkar 12/6/09 axisymmetric anisotropy
c     in the notation used in the notes
c     e1=c11, e2=c12=c21, e3=c66=Gp, e4=c13, and ezz=c33
c     these goto isotropic limit when Ep=Et and Nue-p=Nue-t
               if(stress_anisotropy_in) then
                  young_p= elastic_mod(i)
                  young_t= elastic_mod_t(i)
                  pois_p= poisson(i)
                  pois_t= poisson_t(i)
                  pois_sq= pois_t*pois_t
                  fac1= young_p/young_t
                  fac2= 1.0d0- pois_p -2.0d0*fac1*pois_sq
                  fac3= fac2*(1.0d0+pois_p)
                  e1(i)= young_p*(1.0d0-fac1*pois_sq)/fac3
                  e2(i)= young_p*(pois_p+fac1*pois_sq)/fac3
                  e3(i)= 0.5d0*young_p/(1.d0+pois_p)
                  e4(i)= young_p*pois_t/fac2
                  ezz(i)= young_t*(1.0d0-pois_p)/fac2
c     for thermal expansion, we input Alpha which is a small number, but
c     for pore pressure 
c     we want to be able to input number such that 0<=beta_p<=1 and also
c     have the temperature and pore pressure terms look similalr in the 
c     balance equations.Hence the term beta_p/3Hp is saved, not  beta_p
c     See Keita's notes dated 2/25/2010, Here bulk_mod is
c     defined as bulk_mod=Hp=(C11+C12+C13)/3 and biot=beta_p/3Hp. Then
c     Ks=Hp/(1-beta_p)=bulk_mod/(1-3*bulk_mod*biot).
c     later beta_t calculated from
c     =1-Ht/Ks where Ht=(2C13+C33)/3
                  bulk_mod=(e1(i)+e2(i)+e4(i))/3.
                  if(bulk_mod.gt.bulk_tol) then
                     bulk(i) = bulk0(i)/(3.0*bulk_mod)
                  else
                     bulk(i) = bulk_tol
                  endif
c..................................................
               elseif(stress_anisotropy_use) then
c     calculate the Biot term 
                  bulk_mod = elastic_mod(i)/(3.
     &                 *(1.0d0-2.0d0*poisson(i)))
c     bulk will be biot/(3K)
                  bulk_mod=(e1(i)+e2(i)+e4(i))/3.
                  if(bulk_mod.gt.bulk_tol) then
                     bulk(i) = bulk0(i)/(3.0*bulk_mod)
                  else
                     bulk(i) = bulk_tol
                  endif
               else
                  e1(i) = elastic_mod(i)*(1.0d0-poisson(i))/
     &                 (1.d0+poisson(i))/(1.0d0-2.0d0*poisson(i))
                  e2(i) = e1(i)*poisson(i)/(1.0d0-poisson(i))
                  e3(i) = e1(i)*(1.0d0-2.0d0*poisson(i))/
     &                 2.0d0/(1.0d0-poisson(i))
c     calculate the Biot term 
                  bulk_mod = elastic_mod(i)/(3.
     &                 *(1.0d0-2.0d0*poisson(i)))
c     bulk will be biot/(3K)
                  if(bulk_mod.gt.bulk_tol) then
                     bulk(i) = bulk0(i)/(3.0*bulk_mod)
                  else
                     bulk(i) = bulk_tol
                  endif
               endif
            else
c     plain strain
               e1(i) = elastic_mod(i)/(1.d0-poisson(i)*poisson(i))
               e2(i) = e1(i)*poisson(i)
               e3(i) = e1(i)*(1.0d0-poisson(i))/2.0d0
            endif                    
         enddo
      endif

      return
      end