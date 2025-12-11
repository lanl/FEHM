      subroutine  veloc
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine calculates fluid velocities.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 11-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/veloc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:50   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:02   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:18   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:06 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Mon Jun 10 13:05:58 1996   hend
CD2 Fixed dpdp mdnodes
CD2 
CD2    Rev 1.10   Fri May 24 09:56:14 1996   hend
CD2 Updated for mdnodes
CD2 
CD2    Rev 1.9   Fri Feb 16 13:22:26 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.8   Thu Jan 11 12:17:52 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.7   09/29/95 16:10:14   llt
CD2 added small number on dividing -- was dividing by zero.
CD2 
CD2    Rev 1.6   09/11/95 17:30:14   awolf
CD2 seh fixed in conjunction with rev 1.4 of flxo.f
CD2 
CD2    Rev 1.5   04/05/95 12:41:50   robinson
CD2 Calculated velocity from the larger of inflow and outflow
CD2 
CD2    Rev 1.4   04/03/95 08:48:56   robinson
CD2 Corrected velocities for 2D problems, variable thickness
CD2 
CD2    Rev 1.3   03/10/95 11:10:20   llt
CD2 only calculate velocity from outflow - gaz
CD2 
CD2    Rev 1.2   06/15/94 10:36:06   robinson
CD2 Corrected bug in velocity calc. when saturation goes to 0
CD2 
CD2    Rev 1.1   03/18/94 15:43:22   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:29:16   pvcs
CD2 original version in process of being certified
CD2 
c 1/6/95 gaz only calculate velocity from outflow
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  Not Applicable.  See Special Comments.
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4  
CD4  This is a general utility routine used in the code as needed.
CD4  
c Henderson --> 11 Sep 1995
c This code currently computes a darcy velocity. If you want the 
c velocity particles will travel at, for instance,
c 
c FRACTURES: (Computed Velocity) / (porosity*sat*fracture volume frac)
c
c MATRIX: (Computed Velocity) / (porosity*sat)
CD4
C***********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comji
      use comflow
      use comdti
      use davidi
      implicit none
 
      integer i, ij, iq, jj, kb, neqp1, iwd, iw_max
      real*8 pnxl_in, pnyl_in, pnzl_in, pnxv_in, pnyv_in, pnzv_in
      real*8 pnxl_out, pnyl_out, pnzl_out, pnxv_out, pnyv_out, pnzv_out
      real*8 area_t, axy, axyf, dili, dilkb, dis, divi, divkb 
      real*8 fid, fid1
      real*8 vld_in, vld_out, vvd_in, vvd_out, vxy, vxyf
      real(8) :: vld = 0., vvd = 0.
      real*8 xdis_cos, ydis_cos, zdis_cos
      real*8 tol, dis_tol
      parameter(iw_max = 10000000)
      parameter(tol=1.d-30, dis_tol=1.d-20)
c gaz 11-09-2001 added dis_tol so a velocity can be calculated at
c coincident mdnode nodes

c ************
c     return for anisotropy or heat-only
c
      if(ianpe.ne.0 .or. idoff .eq. -1) return
c     
      vlmax  =  0.0
      vvmax  =  0.0
      neqp1  =  neq+1

c     TAKE CARE OF NODES 1 TO NEQ (ONLY FRACTURES FOR DPDP)
      do i=1,neq
c     zero out velocities
         pnxl_in = 0.
         pnyl_in = 0.
         pnzl_in = 0.
         pnxv_in = 0.
         pnyv_in = 0.
         pnzv_in = 0.
         pnxl_out = 0.
         pnyl_out = 0.
         pnzl_out = 0.
         pnxv_out = 0.
         pnyv_out = 0.
         pnzv_out = 0.
c     calculate lower diagonal geometric coefficient
c     gaz 1-25-03 used abs(iw)      
         iq=0
         do ij=nelm(i)+1,nelmdg(i)-1
            kb=nelm(ij)
            do jj=nelmdg(kb)+1,nelm(kb+1)
               if(nelm(jj).eq.i) then
                  iq=iq+1 
                  it9(iq)=ij
                  it10(iq)=istrw(jj-neqp1)
                  iw=abs(it10(iq))
                  if (sx(iw,isox)+sx(iw,isoy)+
     &                 sx(iw,isoz) .eq.0.) iq=iq-1
               endif
            enddo
         enddo
         do jj=nelmdg(i)+1,nelm(i+1)     
            iq=iq+1 
            it9(iq)=jj
            it10(iq)=istrw(jj-neqp1)
            iw=abs(it10(iq))
            if(iw.eq.0.or.iw.gt.iw_max)then
               iq=iq-1
            else 
               if (sx(iw,isox)+sx(iw,isoy)+
     &              sx(iw,isoz) .eq.0.) iq=iq-1
            endif
         enddo
         do ij=1,iq                   
            jj=it9(ij)
            kb=nelm(jj)
            iwd=it10(ij)
            iw = abs(iwd)
c     use internode fluxes already stored
            axy=a_axy(jj-neqp1)
            if (irdof .ne. 13 .and. jswitch .eq. 0) then
               vxy=a_vxy(jj-neqp1)
            end if
c     velocities
            if(icnl.eq.0) then
               dis=sqrt((cord(i,1)-cord(kb,1))**2
     &              + (cord(i,2)-cord(kb,2))**2
     &              + (cord(i,3)-cord(kb,3))**2)+dis_tol
               area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
               if(iwd.lt.0) area_t=area_t*sx_mult
            else
               dis=sqrt((cord(i,1)-cord(kb,1))**2
     &              + (cord(i,2)-cord(kb,2))**2)+dis_tol
c     area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
               area_t=-(sx(iw,isox)+sx(iw,isoy))*dis
     &              *0.5*(cord(i,3)+cord(kb,3))
               if(iwd.lt.0) area_t=area_t*sx_mult
            endif
            area_t=area_t+tol
            if(axy.gt.0.0) then
               fid=upwgt
               fid1=1.0-fid
            else if(axy.lt.0.0) then
               fid1=upwgt
               fid=1.0-fid1
            else
               fid=0.5
               fid1=0.5
            endif
            dili=dil(i)
            dilkb=dil(kb)
            axyf=max(tol,(fid*dilkb+fid1*dili))
            if(abs(axy).gt.tol) then
               axy=axy/(axyf+tol  )/area_t
     &              *(fid1*dil(i)/(rolf(i)+tol)+
     &              fid*dil(kb)/(rolf(kb)+tol))
               if(axy.gt.0.0) then
                  if(icnl.ne.0) then
                     xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                     ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                     pnxl_in=pnxl_in+xdis_cos*axy
                     pnyl_in=pnyl_in+ydis_cos*axy
                  else
                     xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                     ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                     zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                     pnxl_in=pnxl_in+xdis_cos*axy
                     pnyl_in=pnyl_in+ydis_cos*axy
                     pnzl_in=pnzl_in+zdis_cos*axy
                  endif
               else
                  if(icnl.ne.0) then
                     xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                     ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                     pnxl_out=pnxl_out+xdis_cos*axy
                     pnyl_out=pnyl_out+ydis_cos*axy
                  else
                     xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                     ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                     zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                     pnxl_out=pnxl_out+xdis_cos*axy
                     pnyl_out=pnyl_out+ydis_cos*axy
                     pnzl_out=pnzl_out+zdis_cos*axy
                  endif
               endif
            endif
            if (irdof .ne. 13 .and. jswitch .eq. 0) then
               if(vxy.gt.0.0) then
                  fid=upwgt
                  fid1=1.0-fid
               else if(vxy.lt.0.0) then
                  fid1=upwgt
                  fid=1.0-fid1
               else
                  fid=0.5
                  fid1=0.5
               endif
               divi=div(i)
               divkb=div(kb)
               vxyf=max(tol,(fid*divkb+fid1*divi))
               if(abs(vxy).gt.tol) then
                  vxy=vxy/(vxyf+tol)/area_t
     &                 *(fid1*div(i)/(rovf(i)+tol)+fid*
     &                 div(kb)/(rovf(kb)+tol))
                  if(vxy.gt.0.0) then
                     if(icnl.ne.0) then
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        pnxv_in=pnxv_in+xdis_cos*vxy
                        pnyv_in=pnyv_in+ydis_cos*vxy
                     else
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                        pnxv_in=pnxv_in+xdis_cos*vxy
                        pnyv_in=pnyv_in+ydis_cos*vxy
                        pnzv_in=pnzv_in+zdis_cos*vxy
                     endif
                  else
                     if(icnl.ne.0) then
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        pnxv_out=pnxv_out+xdis_cos*vxy
                        pnyv_out=pnyv_out+ydis_cos*vxy
                     else
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                        pnxv_out=pnxv_out+xdis_cos*vxy
                        pnyv_out=pnyv_out+ydis_cos*vxy
                        pnzv_out=pnzv_out+zdis_cos*vxy
                     endif
                  endif
               endif
            end if
         enddo
c     calculate magnitude of vapor velocity
         if (irdof .ne. 13 .and. jswitch .eq. 0) then
            vvd_in = sqrt(pnxv_in**2 + pnyv_in**2 + pnzv_in**2)
            vvd_out = sqrt(pnxv_out**2 + pnyv_out**2 + pnzv_out**2)
            if( vvd_in .gt. vvd_out ) then
               vvd = vvd_in
               pnx(n*2+i) = pnxv_in
               pny(n*2+i) = pnyv_in
               pnz(n*2+i) = pnzv_in
            else
               vvd = vvd_out
               pnx(n*2+i) = pnxv_out
               pny(n*2+i) = pnyv_out
               pnz(n*2+i) = pnzv_out
            end if
         end if
c     calculate magnitude of liquid velocity
c     
         vld_in = sqrt(pnxl_in**2 + pnyl_in**2 + pnzl_in**2)
         vld_out = sqrt(pnxl_out**2 + pnyl_out**2 + pnzl_out**2)
         if( vld_in .gt. vld_out ) then
            vld = vld_in
            pnx(n+i) = pnxl_in
            pny(n+i) = pnyl_in
            pnz(n+i) = pnzl_in
         else
            vld = vld_out
            pnx(n+i) = pnxl_out
            pny(n+i) = pnyl_out
            pnz(n+i) = pnzl_out
         end if
c     find maximum velocities
c     
         vlmax  =  max( vlmax,abs(vld) )
         vvmax  =  max( vvmax,abs(vvd) )
      enddo

c     TAKE CARE OF NEQ+1 TO N0 IF THEY EXIST (MATRIX FOR DPDP) 
      if (idpdp.ne.0) then
         do i=1,neq
c     zero out velocities
            pnxl_in = 0.
            pnyl_in = 0.
            pnzl_in = 0.
            pnxv_in = 0.
            pnyv_in = 0.
            pnzv_in = 0.
            pnxl_out = 0.
            pnyl_out = 0.
            pnzl_out = 0.
            pnxv_out = 0.
            pnyv_out = 0.
            pnzv_out = 0.
c     calculate lower diagonal geometric coefficient
            iq=0
            do ij=nelm(i)+1,nelmdg(i)-1
               kb=nelm(ij)
               do jj=nelmdg(kb)+1,nelm(kb+1)
                  if(nelm(jj).eq.i) then
                     iq=iq+1 
                     it9(iq)=ij
                     it10(iq)=istrw(jj-neqp1)
                     if (sx(it10(iq),isox)+sx(it10(iq),isoy)+
     &                    sx(it10(iq),isoz) .eq.0.) iq=iq-1
                  endif
               enddo
            enddo
            do jj=nelmdg(i)+1,nelm(i+1)     
               iq=iq+1 
               it9(iq)=jj
               it10(iq)=istrw(jj-neqp1)
               if (sx(it10(iq),isox)+sx(it10(iq),isoy)+
     &              sx(it10(iq),isoz) .eq.0.) iq=iq-1
            enddo
            do ij=1,iq                   
               jj=it9(ij)
               kb=nelm(jj)
               iw=it10(ij)                
c     use internode fluxes already stored
               axy=a_axy(jj-neqp1+nelm(neq+1)-neq-1)
               if (irdof .ne. 13 .and. jswitch .eq. 0) then
                  vxy=a_vxy(jj-neqp1+nelm(neq+1)-neq-1)
               end if
c     velocities
               if(icnl.eq.0) then
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2
     &                 + (cord(i,3)-cord(kb,3))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
               else
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy))*dis
     &                 *0.5*(cord(i,3)+cord(kb,3))
               endif
               area_t=area_t+tol
               if(axy.gt.0.0) then
                  fid=upwgt
                  fid1=1.0-fid
               else if(axy.lt.0.0) then
                  fid1=upwgt
                  fid=1.0-fid1
               else
                  fid=0.5
                  fid1=0.5
               endif
               dili=dil(i+neq)
               dilkb=dil(kb+neq)
               axyf=max(tol,(fid*dilkb+fid1*dili))
               if(abs(axy).gt.tol) then
                  axy=axy/axyf/area_t
     &                 *(fid1*dil(i+neq)/(rolf(i+neq)+tol)+
     &                 fid*dil(kb+neq)/(rolf(kb+neq)+tol))
                  if(axy.gt.0.0) then
                     if(icnl.ne.0) then
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        pnxl_in=pnxl_in+xdis_cos*axy
                        pnyl_in=pnyl_in+ydis_cos*axy
                     else
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                        pnxl_in=pnxl_in+xdis_cos*axy
                        pnyl_in=pnyl_in+ydis_cos*axy
                        pnzl_in=pnzl_in+zdis_cos*axy
                     endif
                  else
                     if(icnl.ne.0) then
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        pnxl_out=pnxl_out+xdis_cos*axy
                        pnyl_out=pnyl_out+ydis_cos*axy
                     else
                        xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                        ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                        zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                        pnxl_out=pnxl_out+xdis_cos*axy
                        pnyl_out=pnyl_out+ydis_cos*axy
                        pnzl_out=pnzl_out+zdis_cos*axy
                     endif
                  endif
               endif
               if (irdof .ne. 13 .and. jswitch .eq. 0) then
                  if(vxy.gt.0.0) then
                     fid=upwgt
                     fid1=1.0-fid
                  else if(vxy.lt.0.0) then
                     fid1=upwgt
                     fid=1.0-fid1
                  else
                     fid=0.5
                     fid1=0.5
                  endif
                  divi=div(i+neq)
                  divkb=div(kb+neq)
                  vxyf=max(tol,(fid*divkb+fid1*divi))
                  if(abs(vxy).gt.tol) then
                     vxy=vxy/vxyf/area_t
     &                    *(fid1*div(i+neq)/(rovf(i+neq)+tol)+
     &                    fid*div(kb+neq)/(rovf(kb+neq)+tol))
                     if(vxy.gt.0.0) then
                        if(icnl.ne.0) then
                           xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                           ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                           pnxv_in=pnxv_in+xdis_cos*vxy
                           pnyv_in=pnyv_in+ydis_cos*vxy
                        else
                           xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                           ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                           zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                           pnxv_in=pnxv_in+xdis_cos*vxy
                           pnyv_in=pnyv_in+ydis_cos*vxy
                           pnzv_in=pnzv_in+zdis_cos*vxy
                        endif
                     else
                        if(icnl.ne.0) then
                           xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                           ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                           pnxv_out=pnxv_out+xdis_cos*vxy
                           pnyv_out=pnyv_out+ydis_cos*vxy
                        else
                           xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                           ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                           zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                           pnxv_out=pnxv_out+xdis_cos*vxy
                           pnyv_out=pnyv_out+ydis_cos*vxy
                           pnzv_out=pnzv_out+zdis_cos*vxy
                        endif
                     endif
                  endif
               endif
            enddo
c     calculate magnitude of vapor velocity
            if (irdof .ne. 13 .and. jswitch .eq. 0) then
               vvd_in = sqrt(pnxv_in**2 + pnyv_in**2 + pnzv_in**2)
               vvd_out =sqrt(pnxv_out**2 + pnyv_out**2 + pnzv_out**2)
               if( vvd_in .gt. vvd_out ) then
                  vvd = vvd_in
                  pnx(n*2+i+neq) = pnxv_in
                  pny(n*2+i+neq) = pnyv_in
                  pnz(n*2+i+neq) = pnzv_in
               else
                  vvd = vvd_out
                  pnx(n*2+i+neq) = pnxv_out
                  pny(n*2+i+neq) = pnyv_out
                  pnz(n*2+i+neq) = pnzv_out
               end if
            end if
c     calculate magnitude of liquid velocity
c     
            vld_in = sqrt(pnxl_in**2 + pnyl_in**2 + pnzl_in**2)
            vld_out = sqrt(pnxl_out**2 + pnyl_out**2 + pnzl_out**2)
c     if( vld_in .gt. vld_out ) then
            if( vld_in .eq. -1.d40  ) then
               vld = vld_in
               pnx(n+i+neq) = pnxl_in
               pny(n+i+neq) = pnyl_in
               pnz(n+i+neq) = pnzl_in
            else
               vld = vld_out
               pnx(n+i+neq) = pnxl_out
               pny(n+i+neq) = pnyl_out
               pnz(n+i+neq) = pnzl_out
            end if
c     find maximum velocities
c     
            vlmax  =  max( vlmax,abs(vld) )
            vvmax  =  max( vvmax,abs(vvd) )
         enddo
      endif

      return
      end
