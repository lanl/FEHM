      subroutine stress_perm(iflg,ndummy)         
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

c     
c     directional_perm_information
c     max and min directional quantities
c     
      
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comsi    
      use comfem, only: permFactor, ifem
      implicit none 

      integer iflg, ndummy, i, iispmd,ipermx_max, ipermy_max
      integer ipermz_max,ipr_str

      real*8 permchng_tol, frac_tol, dis_tol 
      real*8 permx_max,permy_max,permz_max,sxx_min,syy_min
      real*8 szz_min,permsum,permsum0,permchng
      real*8 coorxz_max,cooryz_max,coorzz_max

      parameter (permchng_tol = 0.01,frac_tol=0.00001)
      parameter (dis_tol = 1.d-8)
c     
c     
      if(istrs.eq.0) return
      if(ipermstr.eq.0) return
c     
      if(iflg.eq.0) then

         call allocate_stressperm
c     
c     
c     only calculate for model 3, model 5, model 7 and model 8
c     initial setup calcs node neighbor information
         if(ipermstr3.ne.0.or.ipermstr5.ne.0.or.ipermstr7.ne.0.
     &      .or.ipermstr8.ne.0) then
            call setup_stressperm
         endif
c    
      else if (iflg.eq.-1) then
         call allocate_stressperm_3  
      else if (iflg.eq.-2) then
         call deallocate_stressperm_3  

c     
c     allocate memory for stress derivatives for fully coupled solution
c     just before call to generate equations
c    
      else if (iflg.eq.-3) then
c      only for mode 31 (ice to check for initial ice)
         if(.not.allocated(ispm)) return
c     
         do i = 1,n0
            iispmd = ispm(i)    
            ispmd = ispmt(iispmd) 
            if(ispmd.eq.31) then
c new routine for changing mechanical properties with ice formation          
               call stressperm_31(i)                  
            endif         
         enddo
      else if (iflg.eq.1) then
c     
c     general loop to calculate permeability models and derivatives
c     
c     skip loop if perm macro is read
         if(.not.allocated(ispm)) return

c s kelkar 4/9/12 for now, if using fem, and ihms>0,  
c must use edge-based permfactors. This needs to be generalized
         if(ifem.eq.1.and.ihms.gt.0) then
            write(ierr,*)'From stress_perm. Currently when using ifem=1'
            write(ierr,*)'and ihms>0, nodal permmodels are not allowed' 
            return
         endif
         
c     CV approach, node based perm models 
         do i = 1,n0
            iispmd = ispm(i)    
            ispmd = ispmt(iispmd) 
            if(ispmd.eq.1.and.ihms.eq.1) then
               call stressperm_1(i)
            else if(ispmd.eq.2) then
               call stressperm_2(i)
            else if(ispmd.eq.3) then
               call stressperm_3(i)
            else if(ispmd.eq.4) then
               call stressperm_4(i)
            else if(ispmd.eq.5) then
               call stressperm_5(i)
            else if(ispmd.eq.6) then
               call stressperm_6(i)
            else if(ispmd.eq.7) then
               call stressperm_7(i)
            else if(ispmd.eq.8) then
               call stressperm_8(i)
            else if(ispmd.eq.222) then
               call stressperm_222(i)
            else if(ispmd.eq.21) then
               call stressperm_21(i)
            else if(ispmd.eq.22) then
               call stressperm_22(i)
            else if(ispmd.eq.23) then
               call stressperm_23(i)
            else if(ispmd.eq.24) then
               call stressperm_24(i)
            else if(ispmd.eq.25) then
               call stressperm_25(i)
            else if(ispmd.eq.31) then
c     new routine for changing mechanical properties with ice formation          
               call stressperm_31(i)               
            elseif(ispmd.eq.11) then
               call stressperm_11(i)
c     s kelkar June 9 2009, simple gangi (1978, Int.J.Rock Mech)
            elseif(ispmd.eq.91) then
c     s kelkar June 201, table lookup ....................................
               call stressperm_91(i)
            endif
         enddo
            
      if(spm1f(iispmd).lt.0.0.and.ifem.eq.1) call compute_average_stress

c
c     check for damage zone (permeability changes)
c     and maximum allowable changes
c     
         if(ipermstr2.ne.0.or.ipermstr5.ne.0.or.ipermstr7.ne.0
     &    .or.ipermstr6.ne.0.or.ipermstr8.ne.0)then
            permx_max = 0.0
            permy_max = 0.0
            permz_max = 0.0
            ipermx_max = 1
            ipermy_max = 1
            ipermz_max = 1
            sxx_min = 0.0
            syy_min = 0.0
            szz_min = 0.0
            ipr_str = 0 
            do i = 1,n0
               if(icnl.eq.0) then
                  permsum = pnx(i)+pny(i)+pnz(i)
                  permsum0 = pnx0(i)+pny0(i)+pnz0(i)
               else
                  pnx(i) = min(perx_m*pnx0(i),pnx(i))
                  pny(i) = min(pery_m*pny0(i),pny(i))        
                  permsum = pnx(i)+pny(i)
                  permsum0 = pnx0(i)+pny0(i)       
               endif
               permchng = (permsum - permsum0)/permsum0
               if(permchng.gt.permchng_tol) then
c     count damaged zone nodes
                  ipr_str = ipr_str +1
c     find minimum stresses in the damaged zone           
                  if(str_x(i).lt.sxx_min) then
                     sxx_min  = str_x(i) 
                     ipermx_max = i
                  endif
                  if(str_y(i).lt.syy_min) then
                     syy_min  = str_y(i) 
                     ipermy_max = i
                  endif  
                  if(icnl.eq.0.and.str_z(i).lt.szz_min) then
                     szz_min  = str_z(i) 
                     ipermz_max = i
                  endif                 
               endif
            enddo
c     
c     write out information for the damaged zones    
c     
            if(icnl.eq.0) then
               coorxz_max = cord(ipermx_max,3)
               cooryz_max = cord(ipermy_max,3)
               coorzz_max = cord(ipermz_max,3)
            else
               coorxz_max = 0.0
               cooryz_max = 0.0
               coorzz_max = 0.0        
            endif
            if(iout.ne.0) write(iout,99) l,days 
            if(iptty.ne.0) write(iptty,99) l,days  
            if(iout.ne.0) write(iout,100) ipr_str 
            if(iptty.ne.0) write(iptty,100) ipr_str
            if(iout.ne.0) write(iout,101) 
            if(iptty.ne.0) write(iptty,101)   
            if(iout.ne.0)
     &          write(iout,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max
*******************************
            
c     &  write(iout,201)ipermx_max,str_y(ipermx_max)-phi(ipermx_max),
c     &        frac_flg(ipermx_max)
c*******************************
            if(iptty.ne.0)
     &         write(iptty,102)ipermx_max,sxx_min,pnx(ipermx_max)*1.e-6,
     &           cord(ipermx_max,1),cord(ipermx_max,2),coorxz_max 
            if(iout.ne.0)
     &         write(iout,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            
c     &  write(iout,106)ipermy_max,
c     &       estr_x0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       estr_y0(ipermx_max) - (str_x(ipermx_max)-pho(ipermx_max)),
c     &       pny(ipermx_max)*1.e-6,
c     &       cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
c*******************************
            if(iptty.ne.0)
     &         write(iptty,103)ipermy_max,syy_min,pny(ipermy_max)*1.e-6,
     &           cord(ipermy_max,1),cord(ipermy_max,2),cooryz_max  
            if(icnl.eq.0) then
               if(iout.ne.0)
     &          write(iout,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &              cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max  
               if(iptty.ne.0)
     &         write(iptty,104)ipermz_max,szz_min,pnz(ipermz_max)*1.e-6,
     &              cord(ipermz_max,1),cord(ipermz_max,2),coorzz_max    
            endif     
         endif         
      endif
 99   format(/,1x'Time step ',i6,' Days'1x,f9.2)   
 100  format(1x,'Number of damaged gridblocks (gt 0.01 k/k0 ) ', 1x,i6)
 101  format(1x,'Largest tensile stresses in damaged zone')
 102  format(1x,'Node ',i6,1x,'Sxx ',g12.4,1x,'Kxx ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 103  format(1x,'Node ',i6,1x,'Syy ',g12.4,1x,'Kyy ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 104  format(1x,'Node ',i6,1x,'Szz ',g12.4,1x,'Kzz ',
     &     g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)            
      
 105  format(1x,'Node ',i6,1x,'Sxx ',g12.4,1x,'Kxo',g12.4,1x
     &     'Kxx ',g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)
 106  format(1x,'Node ',i6,1x,'Syy ',g12.4,1x,'Kyo ',g12.4,1x
     &     'Kyy ',g12.4,1x,'x ',g10.4,' y ',g10.4,' z ',g10.4)	 
 201  format(1x,'Node ',i6,1x,'eSyy',g12.4,'frac_flg ',i6)   	 
      
      return 
      end	 
      
c......................................................................
