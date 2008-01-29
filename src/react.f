      subroutine react(i_node_set, temp1,temp2, 
     2     icouple, jcouple,igrp)
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
CD1 To add reaction terms to the Jacobian and the
CD1 residual for each node.
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                     
CD2 Date        Programmer       Number   Comments
CD2
CD2 ?           Hari Viswanathan  N/A     Initial  implementation -
CD2                                       broken out from coneq1
CD2                                       in order to
CD2                                       to allow for the calculation
CD2                                       of the Jacobian cross  
CD2                                       derivative terms
CD2
CD2 $Log:   /pvcs.config/fehm90/src/react.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:46   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:12   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:06:38   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:48 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.4   Thu Feb 01 15:41:42 1996   hend
CD2 Added Requirements Traceability
CD2 
CD2    Rev 1.3   08/30/95 17:02:28   robinson
CD2 Fixed source term calculation
CD2 
CD2    Rev 1.2   08/30/95 13:28:50   robinson
CD2 Fixed problem with dpdp transport with reaction
CD2 
CD2    Rev 1.1   04/25/95 10:04:44   llt
CD2 retrieved lost log history information
CD2
CD2    Rev 1.0   01/28/95 14:05:46   llt
CD2 revised reactive transport module
CD2
C**********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.4 Solute-transport equations
CD3  2.4.6 Multiple, interacting solutes
CD3
C**********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use combi
      use comchem
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comrxni
      use davidi
      use comcouple
      use comdti
      use comai
      implicit none

      integer i_node_set
      integer i
      integer mi
      integer mim
      integer matnum
      integer spec_num
      integer icouple
      integer jcouple
      integer ndummy
      integer iz
      integer jmi
      integer neqp1
      integer jmia
      integer nsizea_single
      integer mat
      integer dmatnum
      integer igrp
      integer temp1
      integer temp2
      if(i_node_set.eq.1)then
         matnum = temp1
         spec_num = temp2
      else
         matnum=temp1+n_couple_species(igrp)*2+1
         spec_num = temp2+1
      endif
      ndummy = (i_node_set-1) * neq
      neqp1 = neq + 1
      nsizea_single=(nelm(neqp1)-neqp1)
      dmatnum = matpos(jcouple+ncpnt*(icouple-1))
      if(icns(nsp).eq.1.or.abs(icns(nsp)).eq.2)then
         do i = 1, neq
            mi=i+ndummy+npn
            mim=mi-npn
            iz = i
            jmi=nelmdg(i)
            jmia=jmi-neqp1
            mat = i+ndummy+(dmatnum-1)*n0
            if(icouple.eq.jcouple)then
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num)) 
     2              + rc(mi)
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))
     2              +drdctaq(mat)
            else
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))
     2              +drdctaq(mat)
            endif   
         end do
      else
         do i = 1, neq
            mi=i+ndummy+npn
            mim=mi-npn
            iz = i
            jmi=nelmdg(i)
            jmia=jmi-neqp1
            if(icouple.eq.jcouple)then
               bp(iz+nrhs_sol(spec_num))=bp(iz+nrhs_sol(spec_num)) 
     2              + rc(mi)
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))
     2              +drdcvap(1,mim)
            else
               a(jmia+nmat_sol(matnum))=a(jmia+nmat_sol(matnum))
     2              +drdcvap(icouple,mim)
            endif   
         end do
      endif
      end


