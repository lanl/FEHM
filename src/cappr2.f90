subroutine cap_pressure(cap_select,sw,it,iphase,mi,cp,dpcp)
!*************************************************************************
! Copyright  2015.   Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
! National  Laboratory  (LANL),  which is operated by  Los Alamos National
! Security, LLC  for the U. S. Department of Energy.  The U. S. Government
! has rights to use, reproduce, and distribute this software.  Neither the
! U. S. Government nor Los Alamos National Security, LLC or persons acting
! on their behalf,  make any warranty,  express or implied, or assumes any
! liability for the accuracy, completeness, or usefulness of the software,
! any information pertaining to the software,  or  represents that its use
! would not infringe privately owned rights.

! The  software  being licensed  may  be Export Controlled.  It may not be
! distributed  or  used by individuals  or entities prohibited from having
! access to the software package, pursuant to United States export control
! laws and regulations. An export control review and determination must be
! completed before LANS will provide access to the identified Software.
!*************************************************************************

  ! cap_select is the type of model (linear, exponential, etc.)
  ! it is the group
  ! iphase refers to the line in the rlpm model
  !  				
  use comrlp , only: cap_param, rlp_fparam,vg1,vg2,vg3,vg4,cp1,cp2  &
       ,cp1f,cp2f,cap_fparam,max_cp,vg1f,vg2f,vg3f,vg4f
  use comco2, only: icarb
  use comai, only: idpdp ,idualp ,neq
  use comdi, only: ieos
  implicit none
  real*8 sw,cp,dpcp,dpcp1,hp,dhp,hp1,dhp1,su_cut,slr,smr,ccp1,ccp2,ccp3
  integer cap_select,it,j,iphase,ireg,mi,itbl,k,mm,i
  ireg=1;su_cut = 0.99d00
  j=(iphase - 1) * max_cp
!  write(*,*) 'line 20 ',sw,cap_param(it,j+2)
 
  select case (cap_select)
  case (0)
     ! No capillary model is defined, do nothing
  case (1)
     !     Linear
     call linear(sw, cap_param(it, j + 1),   &
          cap_param(it, j + 2), cp, dpcp)   

  case(2)				
     !     linear forsythe(1988) model  cp=f(S)
     ccp1=cap_param(it, j + 1);ccp2=cap_param(it, j + 2);ccp3=cap_param(it, j + 3)
                  ccp1=cap_param(it, j + 1)
                  ccp3=cap_param(it, j + 2)
                  if(ieos(mi).eq.2) then   ! two phase
                    cp=ccp1*(ccp2-sw)
                     dpcp=-ccp1
                     if(cp.lt.0.0) then
                        cp=0.0
                        dpcp=0.0
                     endif
                   else if(ieos(mi).eq.3) then  ! single phase gas.  this is only to improve convergence
                     cp=ccp1*ccp3
                     dpcp=0.0  
                   else if(ieos(mi).eq.1) then  ! single phase water
                     cp=0.0        
                     dpcp=0.0  
                  endif

   case (3)
     !     Exponential
     call exponential(sw, cap_param(it, j + 1),    &
          cap_param(it, j + 2), cap_param(it, j + 3),     &
          cap_param(it, j+4) ,cp, dpcp,hp1,dhp1)
  
  case(4)
     ! Brooks-Corey
     call brooks_corey_cap(sw, cap_param(it, j + 1),  &
          cap_param(it, j + 2), cap_param(it, j + 3),  &
          cap_param(it, j + 4), cap_param(it, j + 5),  &
          cap_param(it, j + 6), cp, dpcp)
  
  case(5,6)
     !     Van Genuchten.   

     if(cap_fparam(it, j + 1) .lt. 0.) then   ! does not have fractures
        if (cap_select .eq. 6) then  ! simple linear apprx at low and high S
           call vgcap_ek(sw, cap_param(it, j + 1), cap_param(it, j + 2),  cap_param(it, j + 3),   &
                cap_param(it, j + 4), cap_param(it, j + 5),   &
                cap_param(it, j + 6), cp, dpcp)

        else  ! cubic spline
 
           call vgcap(sw, cap_param(it, j + 1),    &
                cap_param(it, j + 2), cap_param(it, j + 3),  &
                cap_param(it, j + 4), vg1(it, iphase),   &
                vg2(it, iphase), vg3(it, iphase),   &
                vg4(it, iphase), cap_param(it, j + 6),  &
                su_cut, cp1(it, iphase), cp2(it, iphase),  &
                cp, dpcp, ireg)
        endif

     else   ! either not a carbon problem or has fractures
        !     check for fractures, only for air/water case

        if( idpdp .ne. 0 .or. idualp .ne. 0 ) then  ! this is a fracture problem
           if(mi.le.neq) then   ! fractures

              call vgcap(sw, cap_fparam(it, j + 1),  &
                   cap_fparam(it, j + 2),        &
                   cap_fparam(it, j + 3),    &
                   cap_fparam(it, j + 4), vg1f(it, iphase),    & 
                   vg2f(it, iphase), vg3f(it, iphase),     &
                   vg4f(it, iphase), cap_fparam(it, j + 6),      &
                   su_cut, cp1f(it, iphase),           &
                   cp2f(it, iphase),       &
                   cp, dpcp, ireg)
           else  ! matrix
              call vgcap(sw,  cap_param(it, j + 1),     &
                   cap_param(it, j + 2),               &
                   cap_param(it, j + 3),              &
                   cap_param(it, j + 4), vg1(it, iphase),   & 
                   vg2(it, iphase), vg3(it, iphase),      &
                   vg4(it, iphase), cap_param(it, j + 6),   &
                   su_cut, cp1(it, iphase),        &
                   cp2(it, iphase),           &
                   cp, dpcp ,ireg)
                   
           end if
 
!           hp =  cp 
!           dhp = dpcp
  
        else
           !     we're just using matrix cap pressure here
           call vgcap(sw,  cap_param(it, j + 1),       &
                cap_param(it, j + 2), cap_param(it, j + 3),   &
                cap_param(it, j + 4), vg1(it, iphase),       &
                vg2(it, iphase), vg3(it, iphase),      &
                vg4(it, iphase), cap_param(it, j + 6),     &
                su_cut, cp1(it, iphase),          &
                cp2(it, iphase), cp, dpcp, ireg)
           !
           if (cap_select .eq. 7) then
              call vgcap(sw, cap_fparam(it, j + 1),    &
                   cap_fparam(it, j + 2),       &
                   cap_fparam(it, j + 3),       &
                   cap_fparam(it, j + 4), vg1f(it, iphase),    & 
                   vg2f(it, iphase), vg3f(it, iphase),      &
                   vg4f(it, iphase), cap_fparam(it, j + 6),      &
                   su_cut, cp1f(it, iphase),        &
                   cp2f(it, iphase),          &
                   hp1, dhp1, ireg)
 

           else
!              hp1 = hp
!              dhp1 = dhp
           end if
        end if
     end if
     !     conversion from (1/m) to MPa
     ! not sure if hp1, dhp1 need this conversion
        cp = 9.8e-3 * cp
        dpcp=9.8e-3 * dpcp 

  case (9)
     !     values from table, called with wetting phase saturation, table number
     itbl = iphase
     call rlp_cap_table (sw, itbl, 4, cp, dpcp)
     continue
!                  
  end select 
  return
end subroutine cap_pressure

