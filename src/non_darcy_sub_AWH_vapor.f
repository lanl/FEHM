      subroutine nd_flow_AWH_vapor(iflg,jm,neighc,iz,i,kz,kb,
     &                         vxyd,vxy_nd,fid)
c gaz 180925 subroutines with non-darcy AWH flow
      use comai
      use combi
      use com_nondarcy
      use com_nondarcy_AWH, only : vel_nd,sx4d,kij_tol,
     &    kij,aij,dvapi_nd,dvapkb_nd,dvaei_nd,dvaekb_nd,
     &    daxydci,daxydckb,axyd_nd,pxyi,pxy,s_i,s_j,
     &    dcndpi,dcndpj,dandpi,dandpj,dbndpi,dbndpj,
     &    dcndei,dcndej,dandei,dandej,dbndei,dbndej,
     &    dvevpi,dvevpj,dvevei,dvevej,dvevci,dvevcj,
     &    dandci,dandcj,dcndci,dcndcj,dbndci,dbndcj,
     &    dvepi_nd, dvepkb_nd, dvepi_nd, dvepkb_nd,
     &    vxyf,vexyf,vcxyf,dcndcj,vexy_nd,vcxy_nd,
     &    dveekb_nd,dveei_nd,
     &    dvcpkb_nd,dvcpi_nd,dvcekb_nd,dvcei_nd, 
     &    envi,envkb,cnvi,cnvkb,
     &    divpi,divpkb,divei,divekb,divci,divckb,
     &    devi,devkb,devei,devekb,devci,devckb,
     &    dcvi,dcvkb,dcvei,dcvekb,dcvci,dcvckb,
     &    dveckb_nd,dveci_nd,dvcckb_nd,dvcci_nd,
     &    dvackb_nd,dvaci_nd,
     &    divi,divkb,
     &    vxyd_nd,dvelpi,dvelpj,dvelei,dvelej,
     &    dvelci,dvelcj
      use comci
      use comdi
      use davidi
      use comji
      use comfi
      implicit none 
      integer iflg
      integer i,kb,jm,j,ik,ik_max,kz,iz
      integer neighc
      real*8 velij_2,velij,r_vel,phi_j,phi_i,b_nd,a_nd
      real*8 dvel,dr_velv,dpvti,dis2,dis,betaij,fid,fid1
      real*8 der,delx2,dely2,delz2,d_velij_2,c_nd,vxyd,vxy_nd
      fid1 =1.d0-fid
      if(iflg.eq.0) then  
      else if(iflg.eq.1) then
c generate terms required for eq assembly
c axyd units m**2*(area/dis)*Mpa
c initialize grav_wgt
             if(jm.eq.1) grav_wgt = 0.5d0
             if(rovf(i).le.0.0.or.rovf(kb).le.0.0) grav_wgt(jm) = 1.0d0
c gaz 161125
c               vxyd=pxyh+grav_wgt(jm)*sx4d*(rovf(i)+rovf(kb))
c     &        *(cord(kz,igrav)-cord(iz,igrav))
               g_term = grav_wgt(jm)*sx4d*(rovf(i)+rovf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
               dg_termpi = grav_wgt(jm)*sx4d*dgvp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
               dg_termpkb = grav_wgt(jm)*sx4d*dgvp(kb)*
     &                (cord(kz,igrav)-cord(iz,igrav))
               daxydpi = -pxy+dg_termpi
               daxydpkb = pxy+dg_termpkb
               
               daxydei=grav_wgt(jm)*sx4d*dgve(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
               daxydekb=grav_wgt(jm)*sx4d*dgve(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav)) 
c gaz 011025 added deriv wrt ngas variable
               daxydci=grav_wgt(jm)*sx4d*dgvc(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
               daxydckb=grav_wgt(jm)*sx4d*dgvc(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
      else if(iflg.eq.2) then
c calculate velocity in AWH applications      
c this is now called after variable update 
        j = kb
        aij = abs(t5_nd(neighc)) 
        xrv_nd = fid*rvf_nd(kb)+fid1*rvf_nd(i)
        kij = t15(neighc)*1.d-6*xrv_nd
        vel_nd = vxyd/(aij*muvij+kij_tol)
        delx2=(cord(j,1)-cord(i,1))**2
        dely2=(cord(j,2)-cord(i,2))**2
        delz2=(cord(j,3)-cord(i,3))**2            
        dis2=delx2+dely2+delz2
        dis = sqrt(dis2)
        phi_j = phi(j)
        phi_i = phi(i)
c gaz gas phase uses total pressure    
c         if(irdof.ne.13) then
c          phi_j = phi(j)-pcp(j)
c          phi_i = phi(i)-pcp(i)
c         endif
c  pressure is in Pa   
c gaz 090125 (da.mo.yr)     
c        phi_grad = -(1.d6*(phi_j-phi_i)/dis)    
c gaz 300325 
c        rovij = 0.5d0*(rovf(j)+rovf(i))
        betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 

        c_nd = vxyd 
        b_nd = muvij*aij
        a_nd = (kij*aij)*rovij*betaij
c use estimate (darcy vel) of velij 
        velij = vel_nd
c gaz 092224 use NR  
        ik_max = 20
c gaz 030425 velij estimated  in calling routine
c        velij = c_nd/(b_nd+vel_tol_min)
c gaz 290325

         v_tol = vel_tol_min        
c
c         dvel =  max(abs(velij*1.d-3),vel_tol_min)
c gaz 050525 setting tols to best vapor tols
        dvel =  max(abs(velij*1.d-8),vel_tol_min)
c
        do ik = 1, ik_max
c note  velij**2 to velij*abs(velij)
        velij_2 = velij*abs(velij)
c gaz 060425 changed to numerical derivative
        der = 1.d-9*abs(velij) + 1.d-14
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
        write(ierr,444) 'gas', l,i,j,iad,r_vel, velij
444     format(a6,1x,'ts ',i7,' i ',i7,' j ',i7,' iad ',i4,
     &    ' resid ', g14.7,' velij ', g14.7,/)
       write(iout,444) 'gas', l,i,j,iad,r_vel, velij
       continue
       endif
99     continue
      continue               
       vxyd_nd = velij*aij*muvij
      else if(iflg.eq.3) then
c     
c liquid phase calculations
c velocity derivatives
c thermal WH
c gaz 100225 modified for AWH
             j = kb
             aij = abs(t5_nd(neighc)) 
             xrv_nd = fid*rvf_nd(kb)+fid1*rvf_nd(i)
             kij = t15(neighc)*1.d-6*xrv_nd
              delx2=(cord(j,1)-cord(i,1))**2
              dely2=(cord(j,2)-cord(i,2))**2
              delz2=(cord(j,3)-cord(i,3))**2            
              dis2=delx2+dely2+delz2
              dis = sqrt(dis2)
c              if(irdof.ne.13) then
c               phi_j = phi(j)-pcp(j)
c               phi_i = phi(i)-pcp(i)
c              endif
c  pressure is in now Pa includes cap pressure      
               fid1 = 1.d0-fid
               betaij = 0.5d0*(nd_beta(j) + nd_beta(i)) 
               b_nd = muvij*aij
               a_nd = (kij*aij)*rovij*betaij
c gaz 100125
               dbndpi = aij*dmuvijpi
               dbndpj = aij*dmuvijpj
               a_nd = (kij*aij)*rovij*betaij
c gaz 100125   
               dandpi =  (kij*aij)*betaij*drovijpi
               dandpj =  (kij*aij)*betaij*drovijpj
               c_nd =  vxyd
               dcndpi= daxydpi
               dcndpj= daxydpkb
c gaz 021025 
               velij = vxyd_nd/(aij*muvij)
c residual eq: 0.0 = r_vel = -c_nd+b_nd*velij+a_nd*velij_2
c dvelij/dpi:

c gaz mod 270125

             dvel =  max(abs(velij*1.d-3),vel_tol_min)
             velij_2 = velij*abs(velij)
c gaz 080525 tigher tol
c             der = 1.d-3*velij +1.d-8
        der = 1.d-8*abs(velij) + 1.d-14
        d_velij_2 = ((velij+der)*abs(velij+der)-velij*abs(velij))/
     &              der 
  
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
              dcndci = daxydci
              dcndcj = daxydckb
              dbndei = aij*dmuvijei
              dbndej = aij*dmuvijej
              dbndci = aij*dmuvijci
              dbndcj = aij*dmuvijcj
c gaz  060425 added derivative of kij=kij*xrv wrt S
c d(kij)/sj = kij/(xrv_nd-vel_tol_min)*drvef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drovijei + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(i)*fid1
              dandej = (kij*aij)*betaij*drovijej + betaij*rovij*
     &                aij*(kij/(xrv_nd-vel_tol_min))*drvef_nd(j)*fid
              dandci =  (kij*aij)*betaij*drovijci
              dandcj =  (kij*aij)*betaij*drovijcj
      
             dvelei = (-dcndei+velij*(dbndei+dandei*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelej = (-dcndej+velij*(dbndej+dandej*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)    
             dvelci = (-dcndci+velij*(dbndci+dandci*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)  
             dvelcj = (-dcndcj+velij*(dbndcj+dandcj*abs(velij)))
     &                /(-b_nd-a_nd*d_velij_2+vel_tol_min)                           
            endif
            continue

       t13(neighc) = dvelpi*aij*muvij + vel_nd*aij*dmuvijpi
       t14(neighc) = dvelpj*aij*muvij + vel_nd*aij*dmuvijpj
       t18(neighc) = dvelei*aij*muvij + vel_nd*aij*dmuvijei
       t19(neighc) = dvelej*aij*muvij + vel_nd*aij*dmuvijej  	        
       t20a(neighc) = dvelci*aij*muvij + vel_nd*aij*dmuvijci       
       t21a(neighc) = dvelcj*aij*muvij + vel_nd*aij*dmuvijcj
c  gaz 191025
            envkb=envf(kb)
            cnvkb=cnvf(kb)
            divpkb=divp(kb)
            divekb=dive(kb)
            divckb=divc(kb)
            devkb=devf(kb)
            devekb=devef(kb)
            devckb=devcf(kb)
            dcvkb=dcvf(kb)
            dcvekb=dcvef(kb)
            dcvckb=dcvcf(kb)
            divkb=div(kb)
            divei=dive(i)
            divekb=dive(kb)
            vexyf=(fid*divkb*envkb+fid1*divi*envi)
            vxyf=(fid*divkb*(1.0-cnvkb)+fid1*divi*(1.0-cnvi))
            vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
            vxy_nd = vxyd_nd*vxyf 
            vexy_nd= vxyd_nd*vexyf
            vcxy_nd= vxyd_nd*vcxyf

          dvapi_nd  = t13(neighc)*vxyf+vxyd_nd*fid1*
     &              (divpi*(1.0-cnvi)-divi*dcvi)
          dvapkb_nd = t14(neighc)*vxyf+vxyd_nd*fid*
     &               (divpkb*(1.0-cnvkb)-divkb*dcvkb)
          dvaei_nd  = t18(neighc)*vxyf+vxyd_nd*fid1*
     &              (divei*(1.0-cnvi)-divi*dcvei)
          dvaekb_nd = t19(neighc)*vxyf+vxyd_nd*fid*
     &               (divekb*(1.0-cnvkb)-divkb*dcvekb)
          dvaci_nd  = t20a(neighc)*vxyf+vxyd_nd*fid1*
     &              (divci*(1.0-cnvi)-divi*dcvci)
          dvackb_nd = t21a(neighc)*vxyf+vxyd_nd*fid*
     &               (divckb*(1.0-cnvkb)-divkb*dcvckb)
          dvepi_nd  = t13(neighc)*vexyf+vxyd_nd*fid1*
     &              (divpi*envi+divi*devi)
          dvepkb_nd = t14(neighc)*vexyf+vxyd_nd*fid*
     &               (divpkb*envkb+divkb*devkb)  
          dveei_nd  = t18(neighc)*vexyf+vxyd_nd*fid1*
     &              (divei*envi+divi*devei)
          dveekb_nd = t19(neighc)*vexyf+vxyd_nd*fid*
     &               (divekb*envkb+divkb*devekb) 
          dveci_nd  = t20a(neighc)*vexyf+vxyd_nd*fid1*
     &               (devci*divi+envi*divci)
          dveckb_nd = t21a(neighc)*vexyf+vxyd_nd*fid*
     &               (devckb*divkb+envkb*divckb)  
         dvcpi_nd   = t13(neighc)*vcxyf+vxyd_nd*fid1*
     &             (dcvi*divi+cnvi*divpi)
         dvcpkb_nd  = t14(neighc)*vcxyf+vxyd_nd*fid*
     &              (dcvkb*divkb+cnvkb*divpkb)
         dvcei_nd   = t18(neighc)*vcxyf+vxyd_nd*fid1*
     &             (dcvei*divi+cnvi*divei)
         dvcekb_nd  = t19(neighc)*vcxyf+vxyd_nd*fid*
     &              (dcvekb*divkb+cnvkb*divekb)
         dvcci_nd   = t20a(neighc)*vcxyf+vxyd_nd*fid1*
     &             (dcvci*divi+cnvi*divci)
         dvcckb_nd  = t21a(neighc)*vcxyf+vxyd_nd*fid*               
     &              (dcvckb*divkb+cnvkb*divckb)
c gaz 191025
      t11c(neighc)  = dvapi_nd 			  
	t12c(neighc)  =	dvapkb_nd			  
	t13c(neighc)  =	dvaei_nd 		        
	t14c(neighc)  =	dvaekb_nd		        
	t15c(neighc)  =	dvaci_nd 		        
	t16c(neighc)  =	dvackb_nd		               
	t21c(neighc)  =	dvepi_nd                                     		        
	t22c(neighc)  =	dvepkb_nd		        
	t23c(neighc)  =	dveei_nd 	              
	t24c(neighc)  =	dveekb_nd	              
	t25c(neighc)  =	dveci_nd 		        
	t26c(neighc)  =	dveckb_nd	                
	t31c(neighc)  =	dvcpi_nd 		        
	t32c(neighc)  =	dvcpkb_nd 		               
	t33c(neighc)  =	dvcei_nd  		        
	t34c(neighc)  =	dvcekb_nd 		        
	t35c(neighc)  =	dvcci_nd  			  
	t36c(neighc)  =	dvcckb_nd 
      continue
      else if(iflg.eq.4) then
c get nd_derivatives

            j = kb
            envi = envf(i)
            envkb=envf(kb)
            cnvi = cnvf(i)
            cnvkb=cnvf(kb)
            divi=div(i)
            divkb=div(kb)
            vxyd_nd = t8_nd(neighc)
            vxyf=(fid*divkb*(1.0-cnvkb)+fid1*divi*(1.0-cnvi))
            vexyf=(fid*divkb*envkb+fid1*divi*envi)
            vcxyf=(fid*divkb*cnvkb+fid1*divi*cnvi)
            vxy_nd = vxyd_nd*vxyf 
            vexy_nd= vxyd_nd*vexyf
            vcxy_nd= vxyd_nd*vcxyf

      dvapi_nd  = t11c(neighc)
 	dvapkb_nd =	t12c(neighc)
 	dvaei_nd  =	t13c(neighc) 
 	dvaekb_nd =	t14c(neighc) 
 	dvaci_nd  =	t15c(neighc) 
 	dvackb_nd =	t16c(neighc) 
 	dvepi_nd  = t21c(neighc) 
 	dvepkb_nd =	t22c(neighc) 
 	dveei_nd  =	t23c(neighc) 
 	dveekb_nd =	t24c(neighc) 
 	dveci_nd  =	t25c(neighc)       
 	dveckb_nd =	t26c(neighc)       
 	dvcpi_nd  =	t31c(neighc)   
 	dvcpkb_nd =	t32c(neighc) 
 	dvcei_nd  =	t33c(neighc) 
 	dvcekb_nd =	t34c(neighc) 
 	dvcci_nd  =	t35c(neighc)
 	dvcckb_nd =	t36c(neighc)
      continue
      endif
      return
      end