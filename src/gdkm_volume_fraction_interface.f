      subroutine gdkm_volume_fraction_interface(iflg)   
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

c  gaz 080816
c   caculates volume ratios for gdkm model
c   applies volume ratios to flow and diffusion parameters
c   selectively to achieve various gdkm conceptualizations of fracture
c   orientation
c 
c
      use comai
      use combi
      use comco2
      use comdi
      use comdti
      implicit none
      
      integer iflg,i,ncon_size,i1,i2,kb,ik,i_pri,i_sec
      integer jk,kc,j,jj,n_loop
      integer neqp1,imodel, i_dir_gdkm
c gaz 041718      
      integer i_pri_kb, i_dir_gdkm_kb
c gaz 120722
      integer num_fracture_i, num_fracture_kb
      real*8 vfrac,vfrac2,vfrac_sec,tol_dir,tol_dis
      real*8 cordxa,cordya,cordza,cordxb,cordyb,cordzb
      real*8 dl2,dx2,dy2,dz2
      real*8 sx1_primary, sx1_total, red_factor_old
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 disx1, disy1, disz1, disx2, disy2, disz2
      real*8 length_total, length_pri, dis_p
      parameter(tol_dir = 1.d-12,tol_dis = 1.d-6)
      if(gdkm_flag.eq.0) return
      
      if(iflg.eq.-1) then
c allocate memory and other initialization tasks
       if (.not. allocated(gdkm_volume_fraction)) then
        allocate(gdkm_volume_fraction(n0))
        gdkm_volume_fraction = 1.0d00
       endif
      else if(iflg.eq.0) then
c add previous required space to gdkm required space       
       ncon_size=nelm(neq+1) + nitfcpairs+1
       if (allocated(istrw_itfc)) deallocate (istrw_itfc)
        allocate(istrw_itfc(ncon_size))
        istrw_itfc = 0 
c
        if (allocated (red_factor)) deallocate (red_factor)
         allocate(red_factor(0:ncon_size))
         red_factor = 1.0 
      else if(iflg.eq.1) then
c calculate cell lengths for GDPM and GDKM calculations    
       if(.not.allocated(dzrg)) allocate(dzrg(n))
       if(.not.allocated(dyrg)) allocate(dyrg(n))
       if(.not.allocated(dxrg)) allocate(dxrg(n))   
       if(gdkm_flag.ne.0) then   
        n_loop = neq_primary
       else if(idpdp.ne.0) then
        n_loop = neq_primary
       else
        n_loop = n
       endif
       do i = 1, n_loop
          cord1 = cord(i,1)
	    cord2 = cord(i,2)
	    cord3 = 0.0
	    if(icnl.eq.0) then
           cord3 = cord(i,3)
          endif
c find lengths of all connecting primary gridblocks
          disx1 = 0.0
          disy1 = 0.0
          disz1 = 0.0      
          disx2 = 1.e20
          disy2 = 1.e20
          disz2 = 1.e20
          i1 = nelm(i)+1
          i2 = nelm(i+1)
          do jj = i1,i2
           kb = nelm(jj)
           cord1j = cord(kb,1)
           cord2j = cord(kb,2)
           cord3j = 0.0
           if(icnl.eq.0) then
            cord3j = cord(kb,3)
           endif 
              disx1=max(cord1j-cord1,disx1)
              disx2=min(cord1j-cord1,disx2)
              disy1=max(cord2j-cord2,disy1)
              disy2=min(cord2j-cord2,disy2) 
              disz1=max(cord3j-cord3,disz1)
              disz2=min(cord3j-cord3,disz2)                           
          enddo  
            if(ivf.eq.-1) then
               if(disx1.eq.0.0.and.disx2.eq.0.0)
     &             disx1 = abs(cord1-x_orig)*2.
               if(disy1.eq.0.0.and.disy2.eq.0.0) 
     &             disy1 = abs(cord2-y_orig)*2.
               if(disz1.eq.0.0.and.disz1.eq.0.0)
     &             disz1 = abs(cord3-z_orig)*2.
               dzrg(i) = max(disz1,abs(disz2))
               dyrg(i) = max(disy1,abs(disy2))
               dxrg(i) = max(disx1,abs(disx2))
               continue               
            else
               dzrg(i) = abs(disz1-disz2)/2.	
               dyrg(i) = abs(disy1-disy2)/2.    
               dxrg(i) = abs(disx1-disx2)/2.                        
            endif      
       enddo  
 
      else if(iflg.eq.2) then
c transport (tracer) related parameters done in coneqi.f   
c not used
      else if(iflg.eq.3) then 
c calculate and store volume fractions     
       do i = 1,neq_primary
        imodel= igdpm(i)
c first gdpm node is given last connection of primary node
c check this logic
c gaz 120722 check for negative number of fractures        
c        if(gdkm_dir(igdpm(i)).eq.4.and.gdpm_x(igdpm(i),1).lt.0) then
c          gdkm_dir(igdpm(i)) = 0
c          num_fracture_i = abs(gdpm_x(igdpm(i),1))
c          length_pri = sx1(i)**0.333333/((num_fracture_i+1)*2)
c          gdpm_x(igdpm(i),1) = length_pri
c        endif
        kb = nelm(nelm(i+1))
        if(imodel.ne.0) then
          sx1_primary = sx1(i)
	    vfrac = vfrac_primary(imodel)    
c save gdkm volume fraction for node 
           gdkm_volume_fraction(i) = vfrac           
           sx1_total = sx1_primary/vfrac
           do jk = 1, ngdpm_layers(imodel)
            kc = kb+jk-1
            vfrac2 = sx1(kc)/sx1_total 
c save gdkm volume fraction for node 
            gdkm_volume_fraction(kc) = vfrac2
           enddo
        endif
       enddo 
      else if(iflg.eq.4) then   
c calculate interface factors 
c  gaz debug next line       
       ik_gdkm_red = sx(1,1)
       ik_gdkm_red = 0     
       neqp1 = neq +1
       ik = nitfcpairs+1
       do i = 1,neq
        i1 = nelmdg(i) + 1
        i2 = nelm(i+1)
        if(i.le.neq_gdkm) then
         i_dir_gdkm = gdkm_dir(igdpm(i))  
         num_fracture_i = gdpm_x(igdpm(i),1)
        else
         i_pri = nelm(nelm(i)+1)
         num_fracture_i = gdpm_x(igdpm(i_pri),1)
         i_dir_gdkm = gdkm_dir(igdpm(i_pri)) 
        endif
        cordxa=cord(i,1)
        cordya=cord(i,2)
        cordza=cord(i,3)
c identify secondary node associated with primary node i
c if i_sec.lt.neq_gdkm, then not a gdkm node        
        i_pri = nelm(nelmdg(i))
        i_sec = nelm(nelm(i+1))
c insure that the full area is used for gdkm models   
c neq_gdkm is the original primary nodes        
        if(i.le.neq_gdkm.and.i_sec.gt.neq_gdkm) then 
         vfrac_sec = 2.0 + i
        else
         vfrac_sec = 0.0
        endif
        vfrac = gdkm_volume_fraction(i)        
        do j = i1,i2
         kb = nelm(j)
         vfrac2 = gdkm_volume_fraction(kb)
c gaz 041718 defined gdkm model for node kb         
         if(kb.le.neq_gdkm) then
c gaz debug 110722 intersection if fracs different directions             
          i_dir_gdkm = max(i_dir_gdkm,gdkm_dir(igdpm(kb))) 
          num_fracture_kb = gdpm_x(igdpm(kb),1)
         else
          i_pri_kb = nelm(nelm(kb)+1)
          num_fracture_kb = gdpm_x(igdpm(i_pri_kb),1)
          i_dir_gdkm = max(i_dir_gdkm,gdkm_dir(igdpm(i_pri_kb)))  
         endif
         cordxb=cord(kb,1)
 	   cordyb=cord(kb,2)
         cordzb=cord(kb,3)
         dx2 = (cordxb-cordxa)**2
         dy2 = (cordyb-cordya)**2
         dz2 = (cordzb-cordza)**2
         if(icnl.eq.0) then
          dl2 = dx2+dy2+dz2
         else
          dl2 = dx2+dy2
         endif
         if(gdkm_flag.eq.1.and.i_dir_gdkm.eq.1) then
c  fracture plane is orthogonal to x axis               
          ik = ik + 1
          istrw_itfc(j-neqp1) = ik
c check for primary-secondary connection 
          if(vfrac.eq.1.0.and.vfrac2.eq.1.0) then
c non gdkm node to non gdkm node connection
            red_factor(ik) = 1.0  
c vfrac_sec represents the sum i + 2 (so the primary node can be extracted)               
          else if(kb.eq.i_sec.and.vfrac_sec.gt.2.) then 
           red_factor(ik) = vfrac_sec
c check for primary-primary connection in x direction where
c node i is gdkm and kb is non-gdkm node 
c gaz 091016 the next line commented out because orthogonality not required           
          else if(dl2.gt.tol_dis.and.abs(dx2-dl2).le.tol_dir) then
c          else if(dl2.gt.tol_dis) then    
c correction for centered distance for gdkm node  to cell edge          
             length_total = 0.5*(dxrg(i) + dxrg(kb))
             if(vfrac.lt.1.0) then
c node i is gdkm node                 
              length_pri = vfrac   
              dis_p =   length_pri*dxrg(i)/((num_fracture_i+1)*2)
              red_factor(ik) = length_total/(dis_p + 0.5*dxrg(kb))
c gaz 1015016 debug  gaz debug 041718           
c              red_factor(ik)=1.0
             else
c node kb is gdkm node
              length_pri = vfrac2   
              dis_p =   length_pri*dxrg(kb)/((num_fracture_kb+1)*2)
              red_factor(ik) = length_total/(0.5*dxrg(i)+dis_p)
c gaz 1015016 debug  gaz debug 041718            
c              red_factor(ik)=1.0              
             endif
          else 
c whats left is primary-primary in non-gdkm direction or 2nd-2nd in non-gdkm direction 
c gdkm direction : 1 = x, 2 = y, 3 = z              
           red_factor(ik) = 0.5*(vfrac + vfrac2)
           continue
          endif                         
         elseif(gdkm_flag.eq.1.and.i_dir_gdkm.eq.2) then
c  fracture plane is orthogonal to y axis    
          ik = ik + 1
          istrw_itfc(j-neqp1) = ik
c check for primary-secondary connection 
c vfrac_sec represents the sum i + 2 (so the primary node can be extracted)          
          if(kb.eq.i_sec.and.vfrac_sec.gt.2.) then  
           red_factor(ik) = vfrac_sec
c check for primary-primary connection in x direction           
          else if(dl2.gt.tol_dis.and.abs(dy2-dl2).le.tol_dir) then
c          else if(dl2.gt.tol_dis) then   
c correction for centered distance for gdkm node  to cell edge         
             length_total = 0.5*(dyrg(i) + dyrg(kb))
             if(vfrac.lt.1.0) then
c node i is gdkm node                 
              length_pri = vfrac   
              dis_p =   length_pri*dyrg(i)/((num_fracture_i+1)*2)
              red_factor(ik) = length_total/(dis_p + 0.5*dyrg(kb))
c gaz 1015016 debug             
c              red_factor(ik)=1.0                
             else
c node kb is gdkm node
              length_pri = vfrac2   
              dis_p =   length_pri*dyrg(kb)/((num_fracture_kb+1)*2)
              red_factor(ik) = length_total/(0.5*dyrg(i)+dis_p)
c gaz 1015016 debug             
              red_factor(ik)=1.0                
             endif
          else 
c whats left is primary-primary in non-z direction or 2nd-2nd in non-x direction          
           red_factor(ik) = 0.5*(vfrac + vfrac2)
          endif       
         elseif(gdkm_flag.eq.1.and.i_dir_gdkm.eq.3) then
c  fracture plane is orthogonal to z axis 
          ik = ik + 1
          istrw_itfc(j-neqp1) = ik
c check for primary-secondary connection          
          if(kb.eq.i_sec.and.vfrac_sec.gt.2.) then  
           red_factor(ik) = vfrac_sec
c check for primary-primary connection in z direction           
          else if(dl2.gt.tol_dis.and.abs(dz2-dl2).le.tol_dir) then
c          else if(dl2.gt.tol_dis) then   
c correction for centered distance for gdkm node  to cell edge         
             length_total = 0.5*(dzrg(i) + dzrg(kb))
             if(vfrac.lt.1.0) then
c node i is gdkm node                 
              length_pri = vfrac   
              dis_p =   length_pri*dzrg(i)/((num_fracture_i+1)*2)
              red_factor(ik) = length_total/(dis_p + 0.5*dzrg(kb))
c gaz 1015016 debug             
              red_factor(ik)=1.0                
             else
c node kb is gdkm node
              length_pri = vfrac2   
              dis_p =   length_pri*dzrg(kb)/((num_fracture_kb+1)*2)
              red_factor(ik) = length_total/(0.5*dzrg(i)+dis_p)
c gaz 1015016 debug             
              red_factor(ik)=1.0                
             endif
          else 
c whats left is primary-primary in non-z direction or 2nd-2nd in non-z direction          
           red_factor(ik) = 0.5*(vfrac + vfrac2)
          endif
         else
c  fracture plane is "general or non directional"  
           ik = ik + 1
           istrw_itfc(j-neqp1) = ik
c gaz 020217           
           if(kb.eq.i_sec.and.vfrac_sec.gt.2.) then 
            red_factor(ik) = vfrac_sec
           else
            red_factor(ik) = 0.5*(vfrac + vfrac2)
           endif
         endif
        enddo
        enddo
        ik_gdkm_red = ik
        call gdkm_dir_check(1)
      else if(iflg.eq.5) then    
 
      endif
      
      return
      end

      subroutine gdkm_dir_check(iflg)
c gaz 111322
c check and print gdkm coefficients for directional models
c 
c
      use comai
      use combi
      use comco2
      use comdi
      use comdti
      implicit none
      
      integer iflg,i,ncon_size,i1,i2,kb,ik,i_pri,i_sec
      integer jk,kc,j,jj,n_loop, iriver
      integer neqp1,imodel, i_dir_gdkm    
      integer i_pri_kb, i_dir_gdkm_kb
      integer iq, icd, jm, kz, jmi, ii2, ii1,idg 
      integer it11(20),it8(20)
      real*8 vfrac,vfrac2,vfrac_sec,tol_dir,tol_dis
      real*8 cordxa,cordya,cordza,cordxb,cordyb,cordzb
      real*8 dl2,dx2,dy2,dz2
      real*8 sx1_primary, sx1_total, red_factor_old
      real*8 cord1,cord2,cord3,cord1j,cord2j,cord3j
      real*8 disx1, disy1, disz1, disx2, disy2, disz2
      real*8 length_total, length_pri, dis_p
      real*8 reduction_factor
      parameter(tol_dir = 1.d-12,tol_dis = 1.d-6)    
c      
      if(iflg.eq.0) then
c initialization if necessary          
      else if(iflg.eq.1) then
c modify interface factors
      neqp1 = neq +1
      if(i.gt.neq.and.idualp.eq.0) then
      icd=neq
      else
      icd=0
      endif
      do i = 1,neq
      iq = 0
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      jmi=nelmdg(i-icd)
      do jm=jmi+1,ii2
       kb=nelm(jm)+icd
       iq=iq+1
       it8(iq)=kb
       it11(iq)=jm-neqp1
      enddo          
      if(i.le.neq_gdkm) then
         i_dir_gdkm = gdkm_dir(igdpm(i))  
      else
         i_pri = nelm(nelm(i)+1)
         i_dir_gdkm = gdkm_dir(igdpm(i_pri)) 
      endif           
      do jm=1,iq
       kb=it8(jm)
       kz=kb-icd
        Vfrac2 = gdkm_volume_fraction(kb)
         if(kb.le.neq_gdkm) then          
          i_dir_gdkm_kb = gdkm_dir(igdpm(kb)) 
         else
          i_pri_kb = nelm(nelm(kb)+1)
          i_dir_gdkm_kb= gdkm_dir(igdpm(i_pri_kb))  
         endif
         if(i_dir_gdkm.ne.0.and.i_dir_gdkm_kb.ne.0) then
c  apply  volume fraction 
          if(i_dir_gdkm.ne.i_dir_gdkm_kb) then
           vfrac2 = gdkm_volume_fraction(kb)
           vfrac  = gdkm_volume_fraction(i)
           red_factor(istrw_itfc(it11(jm))) = 0.5*(vfrac+vfrac2)
          endif
         endif
         reduction_factor = red_factor(istrw_itfc(it11(jm)))
        enddo
       enddo
      else
      endif
      return
      end