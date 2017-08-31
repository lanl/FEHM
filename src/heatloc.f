      subroutine  heatloc
!***********************************************************************
! Copyright 2014. Los Alamos National Security, LLC.  This material was 
! produced under U.S. Government contract DE-AC52-06NA25396 for Los 
! Alamos National Laboratory (LANL), which is operated by Los Alamos 
! National Security, LLC for the U.S. Department of Energy. The U.S. 
! Government has rights to use, reproduce, and distribute this software.
! Neither the U.S. Government nor Los Alamos National Security, LLC or 
! persons acting on their behalf, make any warranty, express or implied,
! or assumes any liability for the accuracy, completeness, or usefulness
! of the software, any information pertaining to the software, or 
! represents that its use would not infringe privately owned rights.
!
! The software being licensed may be Export Controlled.   It may not be 
! distributed or used by individuals or entities prohibited from having 
! access to the software package, pursuant to United States export 
! control laws and regulations.  An export control review and 
! determination must be completed before LANS will provide access to the
! identified Software.
!***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine calculates heat fluxes.
CD1   patterned after veloc.f
CD1 currently (7/7/14) doing liquid water oinly, single porosity only
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 7-July-2014 S. Kelkar
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
      real*8 adv_x_out,adv_y_out,adv_z_out
      real*8 cond_x_out,cond_y_out,cond_z_out
      real*8 e_cond, e_adv
      real*8 e_cond_area, e_adv_area
      real*8 area_t, dis
      real*8 fid, fid1
      real*8 vld_in, vld_out
      real(8) :: vld = 0., vvd = 0.
      real*8 xdis_cos, ydis_cos, zdis_cos
      real*8 tol, dis_tol
      parameter(iw_max = 10000000)
      parameter(tol=1.d-30, dis_tol=1.d-20)
c gaz 11-09-2001 added dis_tol so a velocity can be calculated at
c coincident mdnode nodes

c ************
c     return for anisotropy
c
      if(ianpe.ne.0) return
c     
      vlmax  =  0.0
      neqp1  =  neq+1

c     TAKE CARE OF NODES 1 TO NEQ (ONLY FRACTURES FOR DPDP)
      do i=1,neq
c     zero out heat fluxes
         adv_x_out = 0.
         adv_y_out = 0.
         adv_z_out = 0.
         cond_x_out = 0.
         cond_y_out = 0.
         cond_z_out = 0.
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
c     s kelkar 7/7/14 for heat flux
            if (idoff .ne. -1) e_adv  = e_axy_adv(jj-neqp1)
            e_cond = e_axy_cond(jj-neqp1)
c................
            
c     outgoing advective heat flux
            if(e_adv.lt.0.0) then
               if(icnl.ne.0) then
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy))*dis
     &                 *0.5*(cord(i,3)+cord(kb,3))
                  if(iwd.lt.0) area_t=area_t*sx_mult
                  area_t=area_t+tol
                  e_adv_area = e_adv/area_t
                  xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                  ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                  adv_x_out=adv_x_out+xdis_cos*e_adv_area
                  adv_y_out=adv_y_out+ydis_cos*e_adv_area
               else
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2
     &                 + (cord(i,3)-cord(kb,3))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
                  if(iwd.lt.0) area_t=area_t*sx_mult
                  area_t=area_t+tol
                  e_adv_area = e_adv/area_t
                  xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                  ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                  zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                  adv_x_out=adv_x_out+xdis_cos*e_adv_area
                  adv_y_out=adv_y_out+ydis_cos*e_adv_area
                  adv_z_out=adv_z_out+zdis_cos*e_adv_area
               endif
            endif
c     outgoing conductive heat flux
            if(e_cond.lt.0.0) then
               if(icnl.ne.0) then
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy))*dis
     &                 *0.5*(cord(i,3)+cord(kb,3))
                  if(iwd.lt.0) area_t=area_t*sx_mult
                  area_t=area_t+tol
                  e_cond_area = e_cond/area_t
                  xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                  ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                  cond_x_out=cond_x_out+xdis_cos*e_cond_area
                  cond_y_out=cond_y_out+ydis_cos*e_cond_area
               else
                  dis=sqrt((cord(i,1)-cord(kb,1))**2
     &                 + (cord(i,2)-cord(kb,2))**2
     &                 + (cord(i,3)-cord(kb,3))**2)+dis_tol
                  area_t=-(sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz))*dis
                  if(iwd.lt.0) area_t=area_t*sx_mult
                  area_t=area_t+tol
                  e_cond_area = e_cond/area_t
                  xdis_cos=(cord(kb,1)-cord(i,1))/dis 
                  ydis_cos=(cord(kb,2)-cord(i,2))/dis 
                  zdis_cos=(cord(kb,3)-cord(i,3))/dis 
                  cond_x_out=cond_x_out+xdis_cos*e_cond_area
                  cond_y_out=cond_y_out+ydis_cos*e_cond_area
                  cond_z_out=cond_z_out+zdis_cos*e_cond_area
               endif
            endif
         enddo
c
         if (idoff .ne. -1) then
            e_adv_nodal(i,1) = adv_x_out
            e_adv_nodal(i,2) = adv_y_out
            e_adv_nodal(i,3) = adv_z_out
         end if
         e_cond_nodal(i,1) = cond_x_out
         e_cond_nodal(i,2) = cond_y_out
         e_cond_nodal(i,3) = cond_z_out
            
      enddo
      
      return
      end subroutine heatloc
