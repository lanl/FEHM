c gaz 280325 subroutines for non-darcy flow
      subroutine nd_flow_vel(iflg,icd,i1,iq,axyd,velij,aij,kij,
     &                      fid,dvelpi,dvelpj,dvelei,dvelej,i,j,jm)
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comwellphys
      use comai
      use com_nondarcy
      implicit none

      integer iflg,inr,i,j,i1,iq,iz,kz
      integer jm,kb,icd,neighc, ik, ik_max
      real*8 r_velij,dr_velij_dv,betaij,velij,fid,fid1
      real*8 axyd,axyf,dlapi,dlapj,dili,dilpi,dilj,dilpj,rol_i
      real*8 rol_j,drolp_i,drolp_j,dvelpi,dvelpj,dvelei,dvelej
      real*8 axyij,dilij,pxy,kij,axy_nd,aij,disij
      real*8 dis, phi_i, phi_j
      real*8 dis2, delx2, dely2, delz2
      real*8 sx4d, pxyi, term1, term2
      real*8 a_nd, b_nd, c_nd, devli, devlj
c gaz 270125 moved v_tol to com_nondarcy
      real*8 r_vel,dr_velv,fac_nd
      real*8  velij0,coef0,coef1,coef2
c gaz 180425 test variables
      real*8 b_nd_0, a_nd_0
      real*8  dcndpi,dcndpj,dandpi,dandpj,dbndpi,dbndpj  
      real*8  dcndei,dcndej,dandei,dandej,dbndei,dbndej  
      real*8 s_i,s_j, s_avg
      real*8 velij_2,d_velij_2, dvel, der, tol_and
      integer i_vel_deriv_test,inorm_aij
      parameter(i_vel_deriv_test = 0, inorm_aij =0)
      if(iflg.eq.0) then
c darcy law liquid phase
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
        phi_j = phi(j)
        phi_i = phi(i)
c  pressure is in Pa 
        if(irdof.ne.13) then     
         phi_j = phi_j - pcp(j)
         phi_i = phi_i - pcp(i)
        endif
         phi_grad = -(1.d6*(phi_j-phi_i)/dis)         
         velij = phi_grad*(kij/muij) + g_term/muij 
         continue
      else if(iflg.eq.-1) then  
c darcy derivatives (Pa/m)
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
c
        dvelpi = 1.d6/dis*(kij/muij) + phi_grad*(-kij/muij**2)*dmuijpi+
     &     dg_termpi/muij
        dvelpj =-1.d6/dis*(kij/muij) + phi_grad*(-kij/muij**2)*dmuijpj+
     &     dg_termpkb/muij
       if(irdof.ne.13) then 
c note muij only depends on total pressure
c phi_i = phi_i - pcp(i)
        dvelei = -1.d6*dpcef(i)/dis*(kij/muij) 
        dvelej = 1.d6*dpcef(j)/dis*(kij/muij) 
       endif
      continue
      else if(iflg.eq.-2) then
c darcy law gas phase
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
        phi_j = phi(j)
        phi_i = phi(i)
c  pressure is in Pa  
c  add     
        phi_grad = -(1.d6*(phi_j-phi_i)/dis)          
        velij = phi_grad*(kij/muvij) + g_term/muvij  
        continue
      else if(iflg.eq.-3) then  
c darcy derivatives (Pa/m)
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
      dvelpi = 1.d6/dis*(kij/muvij) + phi_grad*(-kij/muvij**2)*dmuvijpi+
     &     dg_termpi/muvij
      dvelpj =-1.d6/dis*(kij/muvij) + phi_grad*(-kij/muvij**2)*dmuvijpj+
     &     dg_termpkb/muvij
      continue
      else if(iflg.eq.1) then

c      this is now called after variable update 
    

        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
        phi_j = phi(j)
        phi_i = phi(i)
c gaz 111324     
         if(irdof.ne.13) then
          phi_j = phi(j)-pcp(j)
          phi_i = phi(i)-pcp(i)
         endif
c  pressure is in Pa   
c gaz 090125 (da.mo.yr)     
c        phi_grad = -(1.d6*(phi_j-phi_i)/dis)    
c gaz 300325 
c        rolij = 0.5d0*(rolf(j)+rolf(i))
        betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 

        c_nd = axyd 
        b_nd = muij*aij
        a_nd = (kij*aij)*rolij*betaij
c        a_nd = 0.0d0
c gaz 092224 use NR  
        ik_max = 20
c gaz 030425 velij estimated  in calling routine
c        velij = c_nd/(b_nd+vel_tol_min)
c gaz 290325

         v_tol = vel_tol_min        
c
c         dvel =  max(abs(velij*1.d-3),vel_tol_min)
c gaz 050525 setting tols to best vapor tols
c gaz 080625
        if(inorm_aij.eq.1)  then
          c_nd = c_nd*aij
          b_nd = b_nd*aij
          a_nd = a_nd*aij
          velij = velij*aij
        endif
        dvel =  max(abs(velij*1.d-8),vel_tol_min)
c
        do ik = 1, ik_max
c note  velij**2 to velij*abs(velij)
        velij_2 = velij*abs(velij)
c gaz 060425 changed to numerical derivative
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der
          r_vel = -c_nd+b_nd*velij+a_nd*velij_2
          dr_velv = b_nd + a_nd*d_velij_2
          if(abs(r_vel).lt.dvel.and.ik.gt.1) then
           velij = velij-r_vel/(dr_velv+v_tol)
           go to 99
          else
           velij = velij-r_vel/(dr_velv+v_tol)
          endif
        enddo
c gaz debug 121024
       if(ik.ge.ik_max) then
        write(ierr,444) 'liquid', l,i,j,iad,r_vel, velij
444     format(a6,1x,'ts ',i7,' i ',i7,' j ',i7,' iad ',i4,
     &    ' resid ', g14.7,' velij ', g14.7,/)
       write(iout,444) 'liquid', l,i,j,iad,r_vel, velij
       continue
       endif
99     continue
c gaz  280325 fid calculated after call
c99     fid = 0.5d0
c       if(velij.gt.0.0) then
c         fid= dnwgt
c       else if(velij.lt.0.0) then
c         fid=upwgt
c       endif
        continue
      else if(iflg.eq.2) then
c     
c liquid phase calculations
c velocity derivatives
c

              delx2=(cord(j,1)-cord(i,1))**2
              dely2=(cord(j,2)-cord(i,2))**2
              delz2=(cord(j,3)-cord(i,3))**2            
              dis2=delx2+dely2+delz2
              dis = sqrt(dis2)
              if(irdof.ne.13) then
               phi_j = phi(j)-pcp(j)
               phi_i = phi(i)-pcp(i)
              endif
c  pressure is in now Pa includes cap pressure      
               fid1 = 1.d0-fid
               betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
               b_nd = muij*aij
               a_nd = (kij*aij)*rolij*betaij
c gaz 100125
               dbndpi = aij*dmuijpi
               dbndpj = aij*dmuijpj
               a_nd = (kij*aij)*rolij*betaij
c gaz 100125   
               dandpi =  (kij*aij)*betaij*drolijpi
               dandpj =  (kij*aij)*betaij*drolijpj
               c_nd =  axyd
               dcndpi= daxydpi
               dcndpj= daxydpkb


c residual eq: 0.0 = r_vel = -c_nd+b_nd*velij+a_nd*velij_2
c dvelij/dpi:

c gaz mod 270125

             dvel =  max(abs(velij*1.d-3),vel_tol_min)
             velij_2 = velij*abs(velij)
c             d_velij_2 = 2.0d0*velij + vel_tol_min 
c gaz 080525 tigher tol
c             der = 1.d-3*velij +1.d-8
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der 
c   -dcndei+velij(  not -dcndei-velij(    
             dvelpi = (-dcndpi+velij*(dbndpi+dandpi*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelpj = (-dcndpj+velij*(dbndpj+dandpj*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)   

c gaz add d/ds for cap pressure
            if(irdof.ne.13) then
              s_i = s(i)
              s_j = s(j)
              dcndei = daxydei
              dcndej = daxydekb
              dbndei = aij*dmuijei
              dbndej = aij*dmuijej
c gaz  060425 added derivative of kij=kij*xrl wrt S
c d(kij)/sj = kij/(xrl_nd-vel_tol_min)*drlef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drolijei + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(i)*fid1
c              dandej = (kij*aij)*betaij*drolijej 
              dandej = (kij*aij)*betaij*drolijej + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(j)*fid

c   -dcndei+velij(  not -dcndei-velij(         
             dvelei = (-dcndei+velij*(dbndei+dandei*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelej = (-dcndej+velij*(dbndej+dandej*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)                
            endif
            continue
      else if(iflg.eq.5) then
c     
c liquid phase calculations
c velocity derivatives
c thermal WH
c

              delx2=(cord(j,1)-cord(i,1))**2
              dely2=(cord(j,2)-cord(i,2))**2
              delz2=(cord(j,3)-cord(i,3))**2            
              dis2=delx2+dely2+delz2
              dis = sqrt(dis2)
              if(irdof.ne.13) then
               phi_j = phi(j)-pcp(j)
               phi_i = phi(i)-pcp(i)
              endif
c  pressure is in now Pa includes cap pressure      
               fid1 = 1.d0-fid
               betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
               b_nd = muij*aij
               a_nd = (kij*aij)*rolij*betaij
c gaz 100125
               dbndpi = aij*dmuijpi
               dbndpj = aij*dmuijpj
               a_nd = (kij*aij)*rolij*betaij
c gaz 100125   
               dandpi =  (kij*aij)*betaij*drolijpi
               dandpj =  (kij*aij)*betaij*drolijpj
               c_nd =  axyd
               dcndpi= daxydpi
               dcndpj= daxydpkb


c residual eq: 0.0 = r_vel = -c_nd+b_nd*velij+a_nd*velij_2
c dvelij/dpi:

c gaz mod 270125

             dvel =  max(abs(velij*1.d-3),vel_tol_min)
             velij_2 = velij*abs(velij)
c             d_velij_2 = 2.0d0*velij + vel_tol_min 
c gaz 080525 tigher tol
c             der = 1.d-3*velij +1.d-8
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der 
c   -dcndei+velij(  not -dcndei-velij(    
             dvelpi = (-dcndpi+velij*(dbndpi+dandpi*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelpj = (-dcndpj+velij*(dbndpj+dandpj*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)   

c gaz add d/ds for cap pressure
            if(irdof.ne.13) then
              s_i = s(i)
              s_j = s(j)
              dcndei = daxydei
              dcndej = daxydekb
              dbndei = aij*dmuijei
              dbndej = aij*dmuijej
c gaz  060425 added derivative of kij=kij*xrl wrt S
c d(kij)/sj = kij/(xrl_nd-vel_tol_min)*drlef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drolijei + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(i)*fid1
c              dandej = (kij*aij)*betaij*drolijej 
              dandej = (kij*aij)*betaij*drolijej + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(j)*fid

c   -dcndei+velij(  not -dcndei-velij(         
             dvelei = (-dcndei+velij*(dbndei+dandei*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelej = (-dcndej+velij*(dbndej+dandej*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)                
            endif
            continue
      else if(iflg.eq.3) then
c gaz 310325 check saturation
c       if(xrv_nd.le.vel_tol_min) then        
c        return  
c        continue
c       endif
c calculate gas velocity
c      this is now called after variable update     
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
        phi_j = phi(j)
        phi_i = phi(i)
c  pressure is in Pa        
c        phi_grad = -(1.d6*(phi_j-phi_i)/dis)     
c gaz 103024    
        betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
c c_nd pressure gradient
c        c_nd = phi_grad
c        b_nd = kij/muvij
c        a_nd = (kij/muvij)*rovij*betaij
c axyd is vxyd on pass through
        c_nd = axyd 
        if(i_vel_deriv_test.eq.0) then
         b_nd = muvij*aij
         a_nd = (kij*aij)*rovij*betaij
        else
c gaz 180425 rovij = 2, muvij = 1.d-5
         muvij = 1.d-5
         rovij = 2.0d0
         b_nd = muvij*aij
         a_nd = (kij*aij)*rovij*betaij
        endif
c gaz mod 270125
        ik_max = 20
c gaz 010425 velij initial comes from geneq1_w_nondarcy
c        velij = c_nd/(b_nd+vel_tol_min)*(1.d0-s_avg)

        v_tol = vel_tol_min        
c
        dvel =  max(abs(velij*1.d-8),vel_tol_min)
        velij_2 = velij*abs(velij)
        r_vel = -c_nd+b_nd*velij+a_nd*velij_2
        do ik = 1, ik_max
        velij_2 = velij*abs(velij)
c gaz 060425 changed to numerical derivative
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der
          r_vel = -c_nd+b_nd*velij+a_nd*velij_2
          dr_velv = b_nd + a_nd*d_velij_2
          if(abs(r_vel).lt.dvel.and.ik.gt.1) then
           velij = velij-r_vel/(dr_velv+v_tol)
           go to 199
          else
           velij = velij-r_vel/(dr_velv+v_tol)
          endif
        enddo
        continue
       if(ik.ge.ik_max) then
        write(ierr,444) 'vapor ', l,i,j,iad,r_vel, velij
        write(iout,444) 'vapor ', l,i,j,iad,r_vel, velij 
        continue
       endif
199    continue
c gaz 280325 fid calculated in calling subroutine
c199    fid=0.5d0
c       if(velij.gt.0.0) then
c         fid= dnwgt
c       else if(velij.lt.0.0) then
c         fid=upwgt
c       endif

        continue
      else if(iflg.eq.4) then
c     
c gas phase calculations
c velocity derivatives
c
              delx2=(cord(j,1)-cord(i,1))**2
              dely2=(cord(j,2)-cord(i,2))**2
              delz2=(cord(j,3)-cord(i,3))**2            
              dis2=delx2+dely2+delz2
              dis = sqrt(dis2)
               fid1 = 1.d0-fid
               betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
c gaz 180425               
        c_nd = axyd 
        if(i_vel_deriv_test.eq.0) then
         b_nd = muvij*aij
         a_nd = (kij*aij)*rovij*betaij
        else
c gaz 180425 rovij = 2, muvij = 1.d-5
         muvij = 1.d-5
         rovij = 2.0d0
         b_nd = 1.d-5*aij
         a_nd = (kij*aij)*rovij*betaij
         dmuvijpi = 0.0d0
         dmuvijpj = 0.0d0
         dmuvijei = 0.0d0
         dmuvijej = 0.0d0
         drovijpi = 0.0d0
         drovijpj = 0.0d0
        endif
               dbndpi = -kij*(1.0d0/muvij)**2*dmuvijpi
               dbndpj = -kij*(1.0d0/muvij)**2*dmuvijpj
               dbndei = -kij*(1.0d0/muvij)**2*dmuvijei
               dbndej = -kij*(1.0d0/muvij)**2*dmuvijej
c gaz 100125  
            if(i_vel_deriv_test.eq.0) then 
               dandpi =  (kij*aij)*betaij*drovijpi 
               dandpj =  (kij*aij)*betaij*drovijpj
c gaz  060425 added derivative of kij=kij*xrv wrt S
c d(kij)/sj = kij/(rvf_nd(j)-vel_tol_min)*drvef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drovijei + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(i)*fid1
              dandej = (kij*aij)*betaij*drovijej + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(j)*fid
           else
              dandpi = 0.0d0
              dandpj = 0.0d0
              dandei = 0.0d0
              dandej = 0.0d0
           endif
c axyd is vxyd for velocity
               c_nd =  axyd
               dcndpi= daxydpi
               dcndpj= daxydpkb               
               dcndei= daxydei
               dcndej= daxydekb
c residual eq: 0.0 = -c_nd*b_nd+velij+a_nd*velij**2
c r_vel = -(c_nd*b_nd)+g_term/muvij+velij+a_nd*velij**2
c dvelij/dpi:
c gaz 121324 velij_2
             v_tol = vel_tol_min        
c
             dvel =  max(abs(velij*1.d-8),vel_tol_min)
             velij_2 = velij*abs(velij)
c       
c             d_velij_2 = 2.d0*velij  +  vel_tol_min
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der               
             dvelpi = (-dcndpi+velij*(dbndpi+dandpi*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelpj = (-dcndpj+velij*(dbndpj+dandpj*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)   
             dvelei = (-dcndei+velij*(dbndei+dandei*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelej = (-dcndej+velij*(dbndej+dandej*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)               
            continue 
      else if(iflg.eq.6) then
c     
c gas phase calculations
c velocity derivatives
c thermal WH
c
              delx2=(cord(j,1)-cord(i,1))**2
              dely2=(cord(j,2)-cord(i,2))**2
              delz2=(cord(j,3)-cord(i,3))**2            
              dis2=delx2+dely2+delz2
              dis = sqrt(dis2)
               fid1 = 1.d0-fid
               betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
c gaz 180425               
        c_nd = axyd 
        if(i_vel_deriv_test.eq.0) then
         b_nd = muvij*aij
         a_nd = (kij*aij)*rovij*betaij
        else
c gaz 180425 rovij = 2, muvij = 1.d-5
         muvij = 1.d-5
         rovij = 2.0d0
         b_nd = 1.d-5*aij
         a_nd = (kij*aij)*rovij*betaij
         dmuvijpi = 0.0d0
         dmuvijpj = 0.0d0
         dmuvijei = 0.0d0
         dmuvijej = 0.0d0
         drovijpi = 0.0d0
         drovijpj = 0.0d0
        endif
               dbndpi = -kij*(1.0d0/muvij)**2*dmuvijpi
               dbndpj = -kij*(1.0d0/muvij)**2*dmuvijpj
               dbndei = -kij*(1.0d0/muvij)**2*dmuvijei
               dbndej = -kij*(1.0d0/muvij)**2*dmuvijej
c gaz 100125  
            if(i_vel_deriv_test.eq.0) then 
               dandpi =  (kij*aij)*betaij*drovijpi 
               dandpj =  (kij*aij)*betaij*drovijpj
c gaz  060425 added derivative of kij=kij*xrv wrt S
c d(kij)/sj = kij/(rvf_nd(j)-vel_tol_min)*drvef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drovijei + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(i)*fid1
              dandej = (kij*aij)*betaij*drovijej + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(j)*fid
           else
              dandpi = 0.0d0
              dandpj = 0.0d0
              dandei = 0.0d0
              dandej = 0.0d0
           endif
c axyd is vxyd for velocity
               c_nd =  axyd
               dcndpi= daxydpi
               dcndpj= daxydpkb               
               dcndei= daxydei
               dcndej= daxydekb
c residual eq: 0.0 = -c_nd*b_nd+velij+a_nd*velij**2
c r_vel = -(c_nd*b_nd)+g_term/muvij+velij+a_nd*velij**2
c dvelij/dpi:
c gaz 121324 velij_2
             v_tol = vel_tol_min        
c
             dvel =  max(abs(velij*1.d-8),vel_tol_min)
             velij_2 = velij*abs(velij)
c       
c             d_velij_2 = 2.d0*velij  +  vel_tol_min
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der               
             dvelpi = (-dcndpi+velij*(dbndpi+dandpi*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelpj = (-dcndpj+velij*(dbndpj+dandpj*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)   
             dvelei = (-dcndei+velij*(dbndei+dandei*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelej = (-dcndej+velij*(dbndej+dandej*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)               
            continue 
      endif

      return
      end                                                  
c
      subroutine nd_props(iflg,icd,i1,iq,i,j,jm,fid)
c gaz 102724
c properties for ND vel calc
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comwellphys
      use comai
      use com_nondarcy
      use com_prop_data, only : den_h2o, enth_h2o, visc_h2o, humid_h2o,
     & psat_h2o, den_ngas, enth_ngas, visc_ngas, xnl_ngas, ieval_flag  
  
      implicit none

      integer iflg,inr,i,j,i1,iq,iz,kz
      integer jm,kb,icd,neighc, ik, ik_max
      real*8      rol_i,   rol_j     
      real*8  	   drolp_i, drolp_j	 	 
      real*8  	   drolt_i, drolt_j	 	 
      real*8  	   rov_i,   rov_j	 	 
      real*8  	   ros_i,   ros_j     
      real*8  	   drovp_i, drovp_j	 	 
      real*8  	   drovt_i, drovt_j   
      real*8      enl_i,   enl_j	 	 
      real*8      dhlp_i,  dhlp_j	 	 
      real*8      dhlt_i,  dhlt_j	 	 
      real*8      env_i,   env_j	 	 
      real*8	     ens_i,   ens_j     
      real*8      dhvp_i,  dhvp_j	 	 
      real*8      dhvt_i,  dhvt_j	 	 
      real*8      xvisl_i, xvisl_j   
      real*8 	    dvislp_i,dvislp_j  
      real*8 	    dvislt_i, dvislt_j 
      real*8 	    xvisv_i, xvisv_j   
      real*8 	    vis_i,   vis_j     
      real*8 	    dvisvp_i,dvisvp_j  
      real*8 	    dvisvt_i,dvisvt_j  
c 
      real*8  fid, fid1     
      real*8 xvisl_temp
c gaz 050425
      real*8 pv,dtdp,dtps,dpct,dpsats
c gaz 280325 added testing
      integer itest_nd
      itest_nd = 0
c
      if(iflg.eq.0) then
       fid1 = 1.d0-fid
       if(ico2.lt.0.and.irdof.eq.13) then     
c fully saturated isothermal flow  
c gaz 110324         
        rol_i     =  den_h2o(i,1) 
        rol_j     =  den_h2o(j,1) 
        drolp_i	  =  den_h2o(i,2) 
        drolp_j	  =  den_h2o(j,2) 
        drolijpi = drolp_i*fid1
        drolijpj = drolp_j*fid
        xvisl_i = visc_h2o(i,1) 
        xvisl_j = visc_h2o(j,1)
        dvislp_i = visc_h2o(i,2) 
        dvislp_j = visc_h2o(i,2) 
        muij = xvisl_j*fid + xvisl_i*fid1
        dmuijpi =  dvislp_i*fid1
        dmuijpj =  dvislp_j*fid
        rolij = rol_j*fid + rol_i*fid1

       else if(ico2.lt.0) then
c isothermal 2 phase
c i       
        rol_i     =  den_h2o(i,1) 
        drolp_i	  =  den_h2o(i,2) 
        rov_i	  =  den_h2o(i,4) 
        ros_i     =  rov_i
        drovp_i	  =  den_h2o(i,5) 
        xvisl_i   =  visc_h2o(i,1)
        dvislp_i  =  visc_h2o(i,2)
        xvisv_i   =  visc_h2o(i,4)
        vis_i     =  xvisv_i
        dvisvp_i  =  visc_h2o(i,5)  
c j
        rol_j     =  den_h2o(j,1) 
        drolp_j	  =  den_h2o(j,2) 
        rov_j	  =  den_h2o(j,4) 
        ros_j     =  rov_j
        drovp_j	  =  den_h2o(j,5) 
        xvisl_j   =  visc_h2o(j,1)
        dvislp_j  =  visc_h2o(j,2)
        xvisv_j   =  visc_h2o(j,4)
        vis_j     =  xvisv_j
        dvisvp_j  =  visc_h2o(j,5)  
        muij = xvisl_j*fid + xvisl_i*fid1
        dmuijpi =  dvislp_i*fid1
        dmuijpj =  dvislp_j*fid
        rolij = rol_j*fid + rol_i*fid1
        drolijpi =  drolp_i*fid1
        drolijpj =  drolp_j*fid
        muvij = xvisv_j*fid + xvisv_i*fid1
        dmuvijpi =  dvisvp_i*fid1
        dmuvijpj =  dvisvp_j*fid
        rovij = rov_j*fid + rov_i*fid1
        drovijpi =  drovp_i*fid1
        drovijpj =  drovp_j*fid
       else if(ico2.eq.0) then
c WH
c i       
        rol_i     =  den_h2o(i,1) 
        drolp_i	  =  den_h2o(i,2) 
        drolt_i	  =  den_h2o(i,3) 
        rov_i	  =  den_h2o(i,4) 
        ros_i     =  rov_i
        drovp_i	  =  den_h2o(i,5) 
        drovt_i   =  den_h2o(i,6) 
        enl_i	  =  enth_h2o(i,1)
        dhlp_i	  =  enth_h2o(i,2)
        dhlt_i	  =  enth_h2o(i,3)
        env_i	  =  enth_h2o(i,4)
        ens_i     =  env_i
        dhvp_i	  =  enth_h2o(i,5)
        dhvt_i	  =  enth_h2o(i,6)
        xvisl_i   =  visc_h2o(i,1)
        dvislp_i  =  visc_h2o(i,2)
        dvislt_i  =  visc_h2o(i,3)
        xvisv_i   =  visc_h2o(i,4)
        vis_i     =  xvisv_i
        dvisvp_i  =  visc_h2o(i,5)  
        dvisvt_i  =  visc_h2o(i,6)
c gaz 050425    modify derivatives for 2-phase
         if(ieos(i).eq.2) then
          pv         = psat_h2o(i,1)  
          dtdp       = psat_h2o(i,2)
          dpct       = psat_h2o(i,3)
c gaz 100721 (dpsats needed for vapor pressure lowering)       
          dpsats     = psat_h2o(i,4)      
          dtps       = dtdp       
            drolp_i=drolp_i+drolt_i*dtps
            drovp_i=drovp_i+drovt_i*dtps            
            dhlp_i=dhlp_i+dhlt_i*dtps
            dhvp_i=dhvp_i+dhvt_i*dtps
            dvisvp_i=dvisvp_i+dvisvt_i*dtps
            dvislp_i=dvislp_i+dvislt_i*dtps
            drolt_i=0.0
            drovt_i=0.0
            dhlt_i=0.0
            dhvt_i=0.0
            dvisvt_i=0.0
            dvislt_i=0.0            
         endif
c j
        rol_j     =  den_h2o(j,1) 
        drolp_j	  =  den_h2o(j,2) 
        drolt_j	  =  den_h2o(j,3) 
        rov_j	  =  den_h2o(j,4) 
        ros_j     =  rov_j
        drovp_j	  =  den_h2o(j,5) 
        drovt_j   =  den_h2o(j,6) 
        enl_j	  =  enth_h2o(j,1)
        dhlp_j	  =  enth_h2o(j,2)
        dhlt_j	  =  enth_h2o(j,3)
        env_j	  =  enth_h2o(j,4)
        ens_j     =  env_j
        dhvp_j	  =  enth_h2o(j,5)
        dhvt_j	  =  enth_h2o(j,6)
        xvisl_j   =  visc_h2o(j,1)
        dvislp_j  =  visc_h2o(j,2)
        dvislt_j  =  visc_h2o(j,3)
        xvisv_j   =  visc_h2o(j,4)
        vis_j     =  xvisv_j
        dvisvp_j  =  visc_h2o(j,5)  
        dvisvt_j  =  visc_h2o(j,6)
        muij = xvisl_j*fid + xvisl_i*fid1
        dmuijpi =  dvislp_i*fid1
        dmuijpj =  dvislp_j*fid
        dmuijei =  dvislt_i*fid1
        dmuijej =  dvislt_j*fid  
c gaz 050425    modify derivatives for 2-phase

         if(ieos(j).eq.2) then
          pv         = psat_h2o(j,1)  
          dtdp       = psat_h2o(j,2)
          dpct       = psat_h2o(j,3)
c gaz 100721 (dpsats needed for vapor pressure lowering)       
          dpsats     = psat_h2o(j,4)      
          dtps       = dtdp       
            drolp_j=drolp_j+drolt_j*dtps
            drovp_j=drovp_j+drovt_j*dtps            
            dhlp_j=dhlp_j+dhlt_j*dtps
            dhvp_j=dhvp_j+dhvt_j*dtps
            dvisvp_j=dvisvp_j+dvisvt_j*dtps
            dvislp_j=dvislp_j+dvislt_j*dtps
            drolt_j=0.0
            drovt_j=0.0
            dhlt_j=0.0
            dhvt_j=0.0
            dvisvt_j=0.0
            dvislt_j=0.0            
         endif 
c gaz 120724        
        muvij = xvisv_j*fid + xvisv_i*fid1
        dmuvijpi =  dvisvp_i*fid1
        dmuvijpj =  dvisvp_j*fid
        dmuvijei =  dvisvt_i*fid1
        dmuvijej =  dvisvt_j*fid         
        rolij = rol_j*fid + rol_i*fid1
        drolijpi =  drolp_i*fid1
        drolijpj =  drolp_j*fid
        drolijei =  drolt_i*fid1
        drolijej =  drolt_j*fid        
        enlij = enl_j*fid + enl_i*fid1
        denlijpi =  dhlp_i*fid1
        denlijpj =  dhlp_j*fid
        denlijei =  dhlt_i*fid1
        denlijej =  dhlt_j*fid
        rovij = rov_j*fid + rov_i*fid1
        drovijpi =  drovp_i*fid1
        drovijpj =  drovp_j*fid
        drovijei =  drovt_i*fid1
        drovijej =  drovt_j*fid        
        envij = env_j*fid + env_i*fid1
        denvijpi =  dhvp_i*fid1
        denvijpj =  dhvp_j*fid
        denvijei =  dhvt_i*fid1
        denvijej =  dhvt_j*fid
       if(itest_nd.eq.1) then
         muij = 1.d-4
        dmuijpi =  0.0d0
        dmuijpj =  0.0d0
        dmuijei =  0.0d0
        dmuijej =  0.0d0  
         muvij = 1.d-5
        dmuvijpi =  0.0d0
        dmuvijpj =  0.0d0
        dmuvijei =  0.0d0
        dmuvijej =  0.0d0  
         rolij = 800.d0
        drolijpi =  0.0d0
        drolijpj =  0.0d0
        drolijei =  0.0d0
        drolijej =  0.0d0    
         enlij = 1.d0
        denlijpi =  0.0d0
        denlijpj =  0.0d0
        denlijei =  0.0d0
        denlijej =  0.0d0          
       endif
      else if(ico2.gt.0) then
c AWH   Not completed yet
      endif                  
      else if(iflg.eq.1) then
      endif
      return 
      end

 
      