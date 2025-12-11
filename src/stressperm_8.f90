subroutine stressperm_8(i)
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.      
!***********************************************************************


      use comai
      use combi
      use comdi
      use comsi
      use comdti, only : n0

      implicit none
      integer i, iispmd
      integer kbx1,kbx2,kby1,kby2,kbz1,kbz2,idir

      real*8 ksx, ksy, ksz
      real*8  ipmd8_dsx,ipmd8_dsy,ipmd8_dsz
      real*8  ipmd8_dsxy,ipmd8_dsyz,ipmd8_dsxz
      real*8  ipmd8_tns,ipmd8_clb,cohes,sh_cohe
      real*8  ipmd8_nmx,ipmd8_nmy,ipmd8_nmz
      real*8  ipmd8_shx,ipmd8_shy,ipmd8_shz
      real*8  esxi, esyi, eszi, phid, phib,pi
      real*8  dsx12,dsy12,dsz12,dsxy12,dsxz12,dsyz12
      real*8  L11,L12,L13,L21,L22,L23,L31,L32,L33
      real*8  frac_tol, disz,disy,disx
      real*8 frac_bx, frac_by, frac_bz
      parameter (pi=3.1415926535)

!     **************failure criteria + stress-perm (for intact rock)
         
      iispmd = ispm(i)
      
      if(icnl.ne.0) then 
         
         !     ******** 2D is not implemented yet
         
      endif
      
      !     ******** 3D 
      !     
      if(icnl.eq.0) then 
         if (.not. allocated(es_f_x0)) then
            allocate(es_f_x0(n0,3))
            allocate(es_f_y0(n0,3))
            allocate(es_f_z0(n0,3))           
            allocate(s_f_xy0(n0,3))  
            allocate(s_f_xz0(n0,3))
            allocate(s_f_yz0(n0,3))
            allocate(frc_zen(n0,3))
            allocate(frc_azm(n0,3))
         endif
         frac_bx = spm1f(iispmd)
         frac_by = spm2f(iispmd)
         frac_bz = spm3f(iispmd)
         knx_stressperm = spm4f(iispmd)
         ksx = spm5f(iispmd)
         phid = spm6f(iispmd)
         ipmd8_tns = spm7f(iispmd)
         phib = spm8f(iispmd)
         cohes = spm9f(iispmd)
         
         ipmd8_clb = (1+sin(phib/180*pi))/(1-sin(phib/180*pi))
         sh_cohe=2.*cohes*cos(phib/180*pi)/(1-sin(phib/180*pi))
         
         
         perx_m=1.e4
         pery_m=1.e4
         perz_m=1.e4
         !     Fracture toughness is input (can be a funciton of stress)	    
         kny_stressperm = knx_stressperm
         knz_stressperm = knx_stressperm
         ksy = ksx
         ksz = ksx
         
         ipmd8_nmx = 0.
         ipmd8_nmy = 0.
         ipmd8_nmz = 0.
         ipmd8_shx = 0.
         ipmd8_shy = 0.
         
         !     Distance	    
         kbx1 = ipermx(i,1)
         kbx2 = ipermx(i,2)
         kby1 = ipermy(i,1)
         kby2 = ipermy(i,2)
         kbz1 = ipermz(i,1)
         kbz2 = ipermz(i,2)
         
         disx = (cord(kbx2,1)-cord(kbx1,1))/2
         disy = (cord(kby2,2)-cord(kby1,2))/2
         disz = (cord(kbz2,3)-cord(kbz1,3))/2
         
         !     frac_flg = 0 (no fracture)
         !     = 1 (fra! - xz plane, vertical   -> kx,kz)
         !     = 2 (frac - yz plane, vertical   -> ky,kz)
         !     = 3 (frac - xy plane, horizontal -> kx,ky)
         !     = 4 (frac - xz and xy plane 1&3  -> kx*,ky,kz)
         !     = 5 (frac - yz and xy plane 2&3  -> kx,ky*,kz)
         !     = 6 (frac - xz and yz plane 1&2  -> kx,ky,kz*)
         !     = 7 (frac - all planes           -> kx*,ky*,kz*)	   
         
         !**********************************check rock failure
         esxi = str_x(i)-phi(i)
         esyi = str_y(i)-phi(i)
         eszi = str_z(i)-phi(i)
         
         !**********************************Frac on xz-plane, frac_flg = 1
         if(esyi.lt.ipmd8_tns.or.  &
              ((esyi.gt.esxi).and.(esyi.gt.eszi).and. &
              ((esyi.gt.sh_cohe+eszi*ipmd8_clb).or. &
              (esyi.gt.sh_cohe+esxi*ipmd8_clb)))) then
            !     if(esyi.lt.ipmd8_tns) then
            !     Save stresses when xz-plane frac, normal to y(=2), 
            ! is initiated for futher dilation calc            
            if(frac_flg(i).eq.0.or.frac_flg(i).eq.3.or.frac_flg(i) &  
                 .eq.2.or.frac_flg(i).eq.5) then
               if(esyi.lt.ipmd8_tns) then
                  es_f_y0(i,2) = ipmd8_tns
                  frc_zen(i,2) = 0.
                  frc_azm(i,2) = 0.
               elseif(esyi.gt.sh_cohe+eszi*ipmd8_clb) then
                  es_f_y0(i,2) = esyi
                  frc_zen(i,2) = cohes
                  frc_azm(i,2) = 0.
               elseif(esyi.gt.sh_cohe+esxi*ipmd8_clb) then
                  es_f_y0(i,2) = esyi
                  frc_zen(i,2) = 0.
                  frc_azm(i,2) = cohes
               endif
               es_f_x0(i,2) = esxi
               es_f_z0(i,2) = eszi
               s_f_xy0(i,2) = str_xy(i)
               s_f_yz0(i,2) = str_yz(i)
               s_f_xz0(i,2) = str_xz(i)
               !     save frac directions in frac_flg           
               if(frac_flg(i).eq.0) then
                  frac_flg(i) = 1
               elseif(frac_flg(i).eq.3) then
                  frac_flg(i) = 4
               elseif(frac_flg(i).eq.2) then
                  frac_flg(i) = 6
               elseif(frac_flg(i).eq.5) then
                  frac_flg(i) = 7
               endif
            endif
         endif
            
!**********************************Frac on yz-plane, frac_flg = 2         
          if(esxi.lt.ipmd8_tns.or. &
         ((esxi.gt.esyi).and.(esxi.gt.eszi).and. &
          ((esxi.gt.sh_cohe+eszi*ipmd8_clb) &
          .or.(esxi.gt.sh_cohe+esyi*ipmd8_clb)))) then
!          if(esxi.lt.ipmd8_tns) then
! Save stresses when yz-plane frac, normal to x(=1), 
! is initiated for futher dilation calc            
            if(frac_flg(i).eq.0.or.frac_flg(i).eq.3.or.frac_flg(i).eq.1 &
             .or.frac_flg(i).eq.4) then
              if(esxi.lt.ipmd8_tns) then
                es_f_x0(i,1) = ipmd8_tns
                frc_zen(i,1) = 0.
                frc_azm(i,2) = 0.
              elseif(esxi.gt.sh_cohe+eszi*ipmd8_clb) then
                es_f_x0(i,1) = esxi
                frc_zen(i,1) = cohes
                frc_azm(i,1) = 0.
              elseif(esxi.gt.sh_cohe+esyi*ipmd8_clb) then
                es_f_x0(i,1) = esxi
                frc_zen(i,1) = 0.
                frc_azm(i,1) = cohes
              endif               
              es_f_y0(i,1) = esyi
              es_f_z0(i,1) = eszi
              s_f_xy0(i,1) = str_xy(i)
              s_f_yz0(i,1) = str_yz(i)
              s_f_xz0(i,1) = str_xz(i)
! save frac directions in frac_flg            
              if(frac_flg(i).eq.0) then
                frac_flg(i) = 2
              elseif(frac_flg(i).eq.3) then
                frac_flg(i) = 5
              elseif(frac_flg(i).eq.1) then
                frac_flg(i) = 6
              elseif(frac_flg(i).eq.4) then
                frac_flg(i) = 7
              endif
            endif
          endif
          
!********************************** Frac on xy-plane, frac_flg = 3        
          if(eszi.lt.ipmd8_tns.or. &
         ((eszi.gt.esyi).and.(eszi.gt.esxi).and. &
          ((eszi.gt.sh_cohe+esyi*ipmd8_clb) &
          .or.(eszi.gt.sh_cohe+esxi*ipmd8_clb)))) then
!          if(eszi.lt.ipmd8_tns) then
! Save stresses when xy-plane frac, normal to z(=3), 
! is initiated for futher dilation calc 
             if(frac_flg(i).eq.0.or.frac_flg(i).eq.1.or.frac_flg(i) &
                 .eq.2.or.frac_flg(i).eq.6) then            
                if(eszi.lt.ipmd8_tns) then
                   es_f_z0(i,3) = ipmd8_tns
                   frc_zen(i,3) = 0.
                   frc_azm(i,3) = 0.
                elseif(eszi.gt.sh_cohe+esyi*ipmd8_clb) then
                   es_f_z0(i,3) = eszi
                   frc_zen(i,3) = cohes
                   frc_azm(i,3) = 90.
                elseif(eszi.gt.sh_cohe+esxi*ipmd8_clb) then
                   es_f_z0(i,3) = eszi
                   frc_zen(i,3) = cohes
                   frc_azm(i,3) = 0.
                endif               
                es_f_x0(i,3) = esxi
                es_f_y0(i,3) = esyi
                s_f_xy0(i,3) = str_xy(i)
                s_f_yz0(i,3) = str_yz(i)
                s_f_xz0(i,3) = str_xz(i)
                
                if(frac_flg(i).eq.0) then
                   frac_flg(i) = 3
                elseif(frac_flg(i).eq.1) then
                   frac_flg(i) = 4
                elseif(frac_flg(i).eq.2) then
                   frac_flg(i) = 5
                elseif(frac_flg(i).eq.6) then
                   frac_flg(i) = 7
                endif
             endif         
          endif 
          
!     frac_flg = 1, xz-plane -> kx, kz          
          if(frac_flg(i).eq.1) then
             dsx12  = es_f_x0(i,2) - esxi
             dsy12  = es_f_y0(i,2) - esyi
             dsz12  = es_f_z0(i,2) - eszi
             dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
             dsyz12 = abs(s_f_yz0(i,2)-str_yz(i))
             dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
             L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
             L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
             L13 = -sin(frc_zen(i,2)/180*pi)
             L21 = -sin(frc_azm(i,2)/180*pi)
             L22 = cos(frc_azm(i,2)/180*pi)
             L23 = 0.
             L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
             L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
             L33 = cos(frc_zen(i,2)/180*pi)
             
           ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                 + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
           
           ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
               + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
           
           ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
               + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
           
           ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
               + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
               + (L11*L23+L13*L21)*dsxz12
           
           ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
               + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
               + (L21*L33+L23*L31)*dsxz12
           
           ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
               + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
               + (L11*L33+L13*L31)*dsxz12  
           
                ipmd8_nmy = ipmd8_dsy/kny_stressperm + disy &
               *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i) &
               +alp(i)*(t(i)-tini(i)))
                
                if(esyi.lt.0) then
                   ipmd8_shy = 0.
                else
                   ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy &
                       +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmy = max(0.d0,ipmd8_nmy)
                ipmd8_shy = max(0.d0,ipmd8_shy)
                
                pnx(i) = pnx0(i)+ &
                    1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                pnz(i) = pnz0(i)+ &
                    1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                
!     frac_flg = 2, yz-plane -> ky, kz            
             elseif(frac_flg(i).eq.2) then
                dsx12  = es_f_x0(i,1) - esxi
                dsy12  = es_f_y0(i,1) - esyi
                dsz12  = es_f_z0(i,1) - eszi
                dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
                dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))	
                L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L13 = -sin(frc_zen(i,1)/180*pi)
                L21 = -sin(frc_azm(i,1)/180*pi)
                L22 = cos(frc_azm(i,1)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L33 = cos(frc_zen(i,1)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                    + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                    + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                    + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                    + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                    + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                    + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                    + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                    + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                    + (L11*L33+L13*L31)*dsxz12     
                
                ipmd8_nmx = ipmd8_dsx/knx_stressperm + disx &
                    *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i) &
                    +alp(i)*(t(i)-tini(i)))
                
                if(esxi.lt.0) then
                   ipmd8_shx = 0.
                else
                   ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy &
                       +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmx = max(0.d0,ipmd8_nmx)
                ipmd8_shx = max(0.d0,ipmd8_shx)
                
                pny(i) = pny0(i)+ &
                    1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                pnz(i) = pnz0(i)+ &
                    1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                
                
!     frac_flg = 3, xy-plane -> kx, ky     
             elseif(frac_flg(i).eq.3) then
                dsx12  = es_f_x0(i,3) - esxi
                dsy12  = es_f_y0(i,3) - esyi
                dsz12  = es_f_z0(i,3) - eszi         
                dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
                dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
                
                L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L13 = -sin(frc_zen(i,3)/180*pi)
                L21 = -sin(frc_azm(i,3)/180*pi)
                L22 = cos(frc_azm(i,3)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                    + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                    + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                    + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                    + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                    + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                    + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                    + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                    + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                    + (L11*L33+L13*L31)*dsxz12  
                L33 = cos(frc_zen(i,3)/180*pi)          
                
                ipmd8_nmz = ipmd8_dsz/knz_stressperm + disz &
                    *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i) &
                    +alp(i)*(t(i)-tini(i)))
                
                if(eszi.lt.0) then
                   ipmd8_shz = 0.
                else
                   ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz &
                       +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmz = max(0.d0,ipmd8_nmz)
                ipmd8_shz = max(0.d0,ipmd8_shz)
                
                pnx(i) = pnx0(i)+ &
                    1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                pny(i) = pnz0(i)+ &
                    1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                
!     frac_flg = 4, xz and xy plane -> kx*, ky, kz   
             elseif(frac_flg(i).eq.4) then
!     
!     frac on xz
!     
                dsx12  = es_f_x0(i,2) - esxi
                dsy12  = es_f_y0(i,2) - esyi
                dsz12  = es_f_z0(i,2) - eszi
                dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
                dsyz12 = abs(s_f_yz0(i,2)-str_yz(i))     
                dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
                L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L13 = -sin(frc_zen(i,2)/180*pi)
                L21 = -sin(frc_azm(i,2)/180*pi)
                L22 = cos(frc_azm(i,2)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L33 = cos(frc_zen(i,2)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12  
                
                ipmd8_nmy = ipmd8_dsy/kny_stressperm + disy &
                     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esyi.lt.0) then
                   ipmd8_shy = 0.
                else
                   ipmd8_shy = ((disx/e3(i)+1/ksy)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
!     
!     frac on xy plane
!     
                dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
                dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
                dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
                dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
                dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
                L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L13 = -sin(frc_zen(i,3)/180*pi)
                L21 = -sin(frc_azm(i,3)/180*pi)
                L22 = cos(frc_azm(i,3)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L33 = cos(frc_zen(i,3)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 & 
                   + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12   
               	
                ipmd8_nmz = ipmd8_dsz/knz_stressperm + disz &
                     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(eszi.lt.0) then
                   ipmd8_shz = 0.
                else
                   ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz &
                        +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmy = max(0.d0,ipmd8_nmy)
                ipmd8_shy = max(0.d0,ipmd8_shy)
                ipmd8_nmz = max(0.d0,ipmd8_nmz)
                ipmd8_shz = max(0.d0,ipmd8_shz)
                
                pnx(i) = pnx0(i)+ &
                     1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy) &
                     +1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                pny(i) = pny0(i)+ &
                     1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                pnz(i) = pnz0(i)+ &
                     1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                
!     frac_flg = 5, yz and xy plane -> kx, ky*, kz   
             elseif(frac_flg(i).eq.5) then
!     
!     frac on yz plane
!     
                dsx12  = ipmd8_tns - esxi
                dsy12  = es_f_y0(i,1) - esyi
                dsz12  = es_f_z0(i,1) - eszi
                dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
                dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
                L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L13 = -sin(frc_zen(i,1)/180*pi)
                L21 = -sin(frc_azm(i,1)/180*pi)
                L22 = cos(frc_azm(i,1)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L33 = cos(frc_zen(i,1)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12                         
                
                ipmd8_nmx = ipmd8_dsx/knx_stressperm + disx &
                     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esxi.lt.0) then 
                   ipmd8_shx = 0.
                else
                   ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksy)*ipmd8_dsxz)*tan(phid/180*pi)
                endif
                
!     
!     frac on xy plane
!     
                dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
                dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
                dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
                dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
                dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
                L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L13 = -sin(frc_zen(i,3)/180*pi)
                L21 = -sin(frc_azm(i,3)/180*pi)
                L22 = cos(frc_azm(i,3)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L33 = cos(frc_zen(i,3)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12  
            	
                ipmd8_nmz = ipmd8_dsz/knz_stressperm + disz &
                     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(eszi.lt.0) then
                   ipmd8_shz = 0.
                else
                   ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz &
                        +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)	
                endif
                
                ipmd8_nmx = max(0.d0,ipmd8_nmx)
                ipmd8_shx = max(0.d0,ipmd8_shx)
                ipmd8_nmz = max(0.d0,ipmd8_nmz)
                ipmd8_shz = max(0.d0,ipmd8_shz)
                
                pnx(i) = pnx0(i)+ &
                     1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                pny(i) = pny0(i)+ &
                     1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz) &
                     +1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                pnz(i) = pnz0(i)+ &
                     1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                
!     frac_flg = 6, xz and yz plane -> kx, ky, kz*   
             elseif(frac_flg(i).eq.6) then 
!     
!     frac on yz plane
!     
                dsx12  = es_f_x0(i,1) - esxi
                dsy12  = es_f_y0(i,1) - esyi
                dsz12  = es_f_z0(i,1) - eszi
                dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
                dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))	
                dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
                L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L13 = -sin(frc_zen(i,1)/180*pi)
                L21 = -sin(frc_azm(i,1)/180*pi)
                L22 = cos(frc_azm(i,1)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L33 = cos(frc_zen(i,1)/180*pi)   
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12    
                
                ipmd8_nmx = ipmd8_dsx/knx_stressperm + disx &
                     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esxi.lt.0) then
                   ipmd8_shx = 0.
                else
                   ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
                endif
                
!     
!     frac on xy plane
!     
                dsx12  = es_f_x0(i,2) - (str_x(i)-phi(i))
                dsy12  = es_f_y0(i,2) - (str_y(i)-phi(i))
                dsz12  = es_f_z0(i,2) - (str_z(i)-phi(i))
                dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
                dsyz12 = abs(s_f_yz0(i,2)-str_yz(i)) 
                dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
                L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L13 = -sin(frc_zen(i,2)/180*pi)
                L21 = -sin(frc_azm(i,2)/180*pi)
                L22 = cos(frc_azm(i,2)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L33 = cos(frc_zen(i,2)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12
                
                ipmd8_nmy = ipmd8_dsy/kny_stressperm + disy &
                     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esyi.lt.0) then
                   ipmd8_shy = 0.
                else
                   ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmx = max(0.d0,ipmd8_nmx)
                ipmd8_shx = max(0.d0,ipmd8_shx)
                ipmd8_nmy = max(0.d0,ipmd8_nmy)
                ipmd8_shy = max(0.d0,ipmd8_shy)
                
                pnx(i) = pnx0(i)+ &
                     1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                pny(i) = pny0(i)+ &
                     1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                pnz(i) = pnz0(i)+ &
                     1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx) &
                     +1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy)
                
!     frac_flg = 7, all plane -> kx*, ky*, kz*   
             elseif(frac_flg(i).eq.7) then  
!     
!     frac on yz plane
!     
                dsx12  = es_f_x0(i,1) - esxi
                dsy12  = es_f_y0(i,1) - esyi
                dsz12  = es_f_z0(i,1) - eszi
                dsxy12 = abs(s_f_xy0(i,1)-str_xy(i))
                dsxz12 = abs(s_f_xz0(i,1)-str_xz(i))	
                dsyz12 = abs(s_f_yz0(i,1)-str_yz(i))
                L11 = cos(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L12 = cos(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L13 = -sin(frc_zen(i,1)/180*pi)
                L21 = -sin(frc_azm(i,1)/180*pi)
                L22 = cos(frc_azm(i,1)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,1)/180*pi)*cos(frc_azm(i,1)/180*pi)
                L32 = sin(frc_zen(i,1)/180*pi)*sin(frc_azm(i,1)/180*pi)
                L33 = cos(frc_zen(i,1)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12  
                
                ipmd8_nmx = ipmd8_dsx/knx_stressperm + disx &
                     *((ipmd8_dsx-poisson(i)*(ipmd8_dsy+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esxi.lt.0) then
                   ipmd8_shx = 0.
                else
                   ipmd8_shx = ((disy/e3(i)+1/ksy)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksz)*ipmd8_dsxz)*tan(phid/180*pi)
                endif
                
!     
!     frac on xz plane
!     
                
                dsx12  = es_f_x0(i,2) - (str_x(i)-phi(i))
                dsy12  = es_f_y0(i,2) - (str_y(i)-phi(i))
                dsz12  = es_f_z0(i,2) - (str_z(i)-phi(i))
                dsxy12 = abs(s_f_xy0(i,2)-str_xy(i))
                dsyz12 = abs(s_f_yz0(i,2)-str_yz(i)) 
                dsxz12 = abs(s_f_xz0(i,2)-str_xz(i))
                L11 = cos(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L12 = cos(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L13 = -sin(frc_zen(i,2)/180*pi)
                L21 = -sin(frc_azm(i,2)/180*pi)
                L22 = cos(frc_azm(i,2)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,2)/180*pi)*cos(frc_azm(i,2)/180*pi)
                L32 = sin(frc_zen(i,2)/180*pi)*sin(frc_azm(i,2)/180*pi)
                L33 = cos(frc_zen(i,2)/180*pi)
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12  
                
                ipmd8_nmy = ipmd8_dsy/kny_stressperm + disy &
                     *((ipmd8_dsy-poisson(i)*(ipmd8_dsx+ipmd8_dsz))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(esyi.lt.0) then
                   ipmd8_shy = 0.
                else
                   ipmd8_shy = ((disx/e3(i)+1/ksx)*ipmd8_dsxy &
                        +(disz/e3(i)+1/ksz)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
!     
!     frac on xy plane
!     
                dsx12  = es_f_x0(i,3) - (str_x(i)-phi(i))
                dsy12  = es_f_y0(i,3) - (str_y(i)-phi(i))
                dsz12  = es_f_z0(i,3) - (str_z(i)-phi(i))         
                dsxz12 = abs(s_f_xz0(i,3)-str_xz(i))
                dsyz12 = abs(s_f_yz0(i,3)-str_yz(i))
                dsxy12 = abs(s_f_xy0(i,3)-str_xy(i))
                L11 = cos(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L12 = cos(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L13 = -sin(frc_zen(i,3)/180*pi)
                L21 = -sin(frc_azm(i,3)/180*pi)
                L22 = cos(frc_azm(i,3)/180*pi)
                L23 = 0.
                L31 = sin(frc_zen(i,3)/180*pi)*cos(frc_azm(i,3)/180*pi)
                L32 = sin(frc_zen(i,3)/180*pi)*sin(frc_azm(i,3)/180*pi)
                L33 = cos(frc_zen(i,3)/180*pi)   
                
                ipmd8_dsx = L11*L11*dsx12 + L12*L12*dsy12 + L13*L13*dsz12 &
                     + 2.*L11*L12*dsxy12 + 2.*L11*L13*dsxz12 + 2.*L12*L13*dsyz12
                
                ipmd8_dsy = L21*L21*dsx12 + L22*L22*dsy12 + L23*L23*dsz12 &
                     + 2.*L21*L22*dsxy12 + 2.*L21*L23*dsxz12 + 2.*L22*L23*dsyz12
                
                ipmd8_dsz = L31*L31*dsx12 + L32*L32*dsy12 + L33*L33*dsz12 &
                     + 2.*L31*L32*dsxy12 + 2.*L31*L33*dsxz12 + 2.*L32*L33*dsyz12
                
                ipmd8_dsxy = L11*L21*dsx12 + L12*L22*dsy12 + L13*L23*dsz12 &
                     + (L11*L22+L12*L21)*dsxy12 + (L12*L23+L13*L22)*dsyz12 &
                     + (L11*L23+L13*L21)*dsxz12
                
                ipmd8_dsyz = L21*L31*dsx12 + L22*L32*dsy12 + L23*L33*dsz12 &
                     + (L21*L32+L22*L31)*dsxy12 + (L22*L33+L23*L32)*dsyz12 &
                     + (L21*L33+L23*L31)*dsxz12
                
                ipmd8_dsxz = L11*L31*dsx12 + L12*L32*dsy12 + L13*L33*dsz12 &
                     + (L11*L32+L12*L31)*dsxy12 + (L12*L33+L13*L32)*dsyz12 &
                     + (L11*L33+L13*L31)*dsxz12  
                
                ipmd8_nmz = ipmd8_dsz/knz_stressperm + disz &
                     *((ipmd8_dsz-poisson(i)*(ipmd8_dsx+ipmd8_dsy))/elastic_mod(i) &
                     +alp(i)*(t(i)-tini(i)))
                
                if(eszi.lt.0) then
                   ipmd8_shz = 0.
                else
                   ipmd8_shz = ((disx/e3(i)+1/ksx)*ipmd8_dsxz &
                        +(disy/e3(i)+1/ksy)*ipmd8_dsyz)*tan(phid/180*pi)
                endif
                
                ipmd8_nmx = max(0.d0,ipmd8_nmx)
                ipmd8_shx = max(0.d0,ipmd8_shx)
                ipmd8_nmy = max(0.d0,ipmd8_nmy)
                ipmd8_shy = max(0.d0,ipmd8_shy)
                ipmd8_nmz = max(0.d0,ipmd8_nmz)
                ipmd8_shz = max(0.d0,ipmd8_shz)                
                
                pnx(i) = pnx0(i)+ &
                     1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy) &
                     +1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz)
                pny(i) = pny0(i)+ &
                     1.e6*(frac_bz+ipmd8_nmz+ipmd8_shz)**3/(12*disz) &
                     +1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx)
                pnz(i) = pnz0(i)+ &
                     1.e6*(frac_bx+ipmd8_nmx+ipmd8_shx)**3/(12*disx) &
                     +1.e6*(frac_by+ipmd8_nmy+ipmd8_shy)**3/(12*disy) 
                
                
             endif          
             pnx(i) = min(perx_m*pnx0(i),pnx(i))
             pny(i) = min(pery_m*pny0(i),pny(i))
             pnz(i) = min(perz_m*pnz0(i),pnz(i))
             
!     pnx(i) = max(pnx0(i),pnx(i))
!     pny(i) = max(pny0(i),pny(i))  
!     pnz(i) = max(pnz0(i),pnz(i))
             
             
          endif 
          
          return

          end
!......................................................................
