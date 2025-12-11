      subroutine nd_flow_AWH_liquid(iflg,jm,neighc,iz,i,kz,kb,axyd,
     &                              axy_nd,fid)
c gaz 180925 subroutines with non-darcy AWH flow
      use comai
      use combi
      use com_nondarcy
      use com_nondarcy_AWH, only : vel_nd,sx4d,kij_tol,
     &    kij,aij,dlapi_nd,dlapkb_nd,dlaei_nd,dlaekb_nd,
     &    daxydci,daxydckb,axyd_nd,pxyi,pxy,s_i,s_j,
     &    dcndpi,dcndpj,dandpi,dandpj,dbndpi,dbndpj,
     &    dcndei,dcndej,dandei,dandej,dbndei,dbndej,
     &    dvelpi,dvelpj,dvelei,dvelej,dvelci,dvelcj,
     &    dandci,dandcj,dcndci,dcndcj,dbndci,dbndcj,
     &    dlepi_nd, dlepkb_nd, dvepi_nd, dvepkb_nd,
     &    axyf,aexyf,acxyf,dcndcj,aexy_nd,acxy_nd,
     &    dleekb_nd,dleei_nd,
     &    dlcpkb_nd,dlcpi_nd,dlcekb_nd,dlcei_nd, 
     &    enli,enlkb,cnli,cnlkb,
     &    dilpi,dilpkb,dilei,dilekb,dilci,dilckb,
     &    deli,delkb,delei,delekb,delci,delckb,
     &    dcli,dclkb,dclei,dclekb,dclci,dclckb,
     &    dleckb_nd,dleci_nd,dlcckb_nd,dlcci_nd,
     &    dlackb_nd,dlaci_nd,
     &    dili,dilkb
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
      real*8 der,delx2,dely2,delz2,d_velij_2,c_nd,axyd,axy_nd
      fid1 =1.d0-fid
      if(iflg.eq.0) then  
      else if(iflg.eq.1) then
c generate terms required for eq assembly
c axyd units m**2*(area/dis)*Mpa
c initialize grav_wgt
             if(jm.eq.1) grav_wgt = 0.5d0
             if(rolf(i).le.0.0.or.rolf(kb).le.0.0) grav_wgt(jm) = 1.0d0
               axyd=pxyi+grav_wgt(jm)*sx4d*(rolf(i)+rolf(kb))
     &        *(cord(kz,igrav)-cord(iz,igrav))
               g_term = grav_wgt(jm)*sx4d*(rolf(i)+rolf(kb))
     &              *(cord(kz,igrav)-cord(iz,igrav))
               dg_termpi = grav_wgt(jm)*sx4d*dglp(i)*
     &                 (cord(kz,igrav)-cord(iz,igrav))
               dg_termpkb = grav_wgt(jm)*sx4d*dglp(kb)*
     &                (cord(kz,igrav)-cord(iz,igrav))
               daxydpi = -pxy+dg_termpi
               daxydpkb = pxy+dg_termpkb
               
               daxydei=pxy*dpcef(i)+grav_wgt(jm)*sx4d*dgle(i)*
     &                    (cord(kz,igrav)-cord(iz,igrav))
               daxydekb=-pxy*dpcef(kb)+grav_wgt(jm)*sx4d*dgle(kb)
     &                    *(cord(kz,igrav)-cord(iz,igrav)) 
c gaz 011025 added deriv wrt ngas variable
               daxydci=grav_wgt(jm)*sx4d*dglc(i)*
     &           (cord(kz,igrav)-cord(iz,igrav))
               daxydckb=grav_wgt(jm)*sx4d*dglc(kb)*
     &            (cord(kz,igrav)-cord(iz,igrav))
      else if(iflg.eq.2) then
c calculate velocity in AWH applications      
c this is now called after variable update 
        j = kb
        aij = abs(t5_nd(neighc)) 
        xrl_nd = fid*rlf_nd(kb)+fid1*rlf_nd(i)
        kij = t15(neighc)*1.d-6*xrl_nd
        vel_nd = axyd/(aij*muij+kij_tol)
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
        write(ierr,444) 'liquid', l,i,j,iad,r_vel, velij
444     format(a6,1x,'ts ',i7,' i ',i7,' j ',i7,' iad ',i4,
     &    ' resid ', g14.7,' velij ', g14.7,/)
       write(iout,444) 'liquid', l,i,j,iad,r_vel, velij
       continue
       endif
99     continue
      continue               
       axyd_nd = velij*aij*muij
      else if(iflg.eq.3) then
c     
c liquid phase calculations
c velocity derivatives
c thermal WH
c gaz 100225 modified for AWH
             j = kb
             aij = abs(t5_nd(neighc)) 
             xrl_nd = fid*rlf_nd(kb)+fid1*rlf_nd(i)
             kij = t15(neighc)*1.d-6*xrl_nd
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
c gaz 021025 
               velij = axyd_nd/(aij*muij)
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
              dbndei = aij*dmuijei
              dbndej = aij*dmuijej
              dbndci = aij*dmuijci
              dbndcj = aij*dmuijcj
c gaz  060425 added derivative of kij=kij*xrl wrt S
c d(kij)/sj = kij/(xrl_nd-vel_tol_min)*drlef_nd(i)*fid1
              dandei = (kij*aij)*betaij*drolijei + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(i)*fid1
              dandej = (kij*aij)*betaij*drolijej + betaij*rolij*
     &                aij*(kij/(xrl_nd-vel_tol_min))*drlef_nd(j)*fid
              dandci =  (kij*aij)*betaij*drolijci
              dandcj =  (kij*aij)*betaij*drolijcj
      
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

       t13(neighc) = dvelpi*aij*muij + vel_nd*aij*dmuijpi
       t14(neighc) = dvelpj*aij*muij + vel_nd*aij*dmuijpj
       t18(neighc) = dvelei*aij*muij + vel_nd*aij*dmuijei
       t19(neighc) = dvelej*aij*muij + vel_nd*aij*dmuijej  	        
      	t20a(neighc) = dvelci*aij*muij + vel_nd*aij*dmuijci       
       t21a(neighc) = dvelcj*aij*muij + vel_nd*aij*dmuijcj
c  gaz 191025
            enlkb=enlf(kb)
            cnlkb=cnlf(kb)
            dilpkb=dilp(kb)
            dilekb=dile(kb)
            dilckb=dilc(kb)
            delkb=delf(kb)
            delekb=delef(kb)
            delckb=delcf(kb)
            dclkb=dclf(kb)
            dclekb=dclef(kb)
            dclckb=dclcf(kb)
            dilkb=dil(kb)
            dilei=dile(i)
            dilekb=dile(kb)
            aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            axyf=(fid*dilkb*(1.0-cnlkb)+fid1*dili*(1.0-cnli))
            acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
            axy_nd = axyd_nd*axyf 
            aexy_nd= axyd_nd*aexyf
            acxy_nd= axyd_nd*acxyf

          dlapi_nd  = t13(neighc)*axyf+axyd_nd*fid1*
     &              (dilpi*(1.0-cnli)-dili*dcli)
          dlapkb_nd = t14(neighc)*axyf+axyd_nd*fid*
     &               (dilpkb*(1.0-cnlkb)-dilkb*dclkb)
          dlaei_nd  = t18(neighc)*axyf+axyd_nd*fid1*
     &              (dilei*(1.0-cnli)-dili*dclei)
          dlaekb_nd = t19(neighc)*axyf+axyd_nd*fid*
     &               (dilekb*(1.0-cnlkb)-dilkb*dclekb)
          dlaci_nd  = t20a(neighc)*axyf+axyd_nd*fid1*
     &              (dilci*(1.0-cnli)-dili*dclci)
          dlackb_nd = t21a(neighc)*axyf+axyd_nd*fid*
     &               (dilckb*(1.0-cnlkb)-dilkb*dclckb)
          dlepi_nd  = t13(neighc)*aexyf+axyd_nd*fid1*
     &              (dilpi*enli+dili*deli)
          dlepkb_nd = t14(neighc)*aexyf+axyd_nd*fid*
     &               (dilpkb*enlkb+dilkb*delkb)  
          dleei_nd  = t18(neighc)*aexyf+axyd_nd*fid1*
     &              (dilei*enli+dili*delei)
          dleekb_nd = t19(neighc)*aexyf+axyd_nd*fid*
     &               (dilekb*enlkb+dilkb*delekb) 
          dleci_nd  = t20a(neighc)*aexyf+axyd_nd*fid1*
     &               (delci*dili+enli*dilci)
          dleckb_nd = t21a(neighc)*aexyf+axyd_nd*fid*
     &               (delckb*dilkb+enlkb*dilckb)  
         dlcpi_nd   = t13(neighc)*acxyf+axyd_nd*fid1*
     &             (dcli*dili+cnli*dilpi)
         dlcpkb_nd  = t14(neighc)*acxyf+axyd_nd*fid*
     &              (dclkb*dilkb+cnlkb*dilpkb)
         dlcei_nd   = t18(neighc)*acxyf+axyd_nd*fid1*
     &             (dclei*dili+cnli*dilei)
         dlcekb_nd  = t19(neighc)*acxyf+axyd_nd*fid*
     &              (dclekb*dilkb+cnlkb*dilekb)
         dlcci_nd   = t20a(neighc)*acxyf+axyd_nd*fid1*
     &             (dclci*dili+cnli*dilci)
         dlcckb_nd  = t21a(neighc)*acxyf+axyd_nd*fid*               
     &              (dclckb*dilkb+cnlkb*dilckb)
c gaz 191025
      t11c(neighc)  = dlapi_nd 			  
	t12c(neighc)  =	dlapkb_nd			  
	t13c(neighc)  =	dlaei_nd 		        
	t14c(neighc)  =	dlaekb_nd		        
	t15c(neighc)  =	dlaci_nd 		        
	t16c(neighc)  =	dlackb_nd		               
	t21c(neighc)  =	dlepi_nd                                     		        
	t22c(neighc)  =	dlepkb_nd		        
	t23c(neighc)  =	dleei_nd 	              
	t24c(neighc)  =	dleekb_nd	              
	t25c(neighc)  =	dleci_nd 		        
	t26c(neighc)  =	dleckb_nd	                
	t31c(neighc)  =	dlcpi_nd 		        
	t32c(neighc)  =	dlcpkb_nd 		               
	t33c(neighc)  =	dlcei_nd  		        
	t34c(neighc)  =	dlcekb_nd 		        
	t35c(neighc)  =	dlcci_nd  			  
	t36c(neighc)  =	dlcckb_nd 
      continue
      else if(iflg.eq.4) then
c get nd_derivatives

            j = kb
            enli = enlf(i)
            enlkb=enlf(kb)
            cnli = cnlf(i)
            cnlkb=cnlf(kb)
            dili=dil(i)
            dilkb=dil(kb)
            axyd_nd = t8_nd(neighc)
            axyf=(fid*dilkb*(1.0-cnlkb)+fid1*dili*(1.0-cnli))
            aexyf=(fid*dilkb*enlkb+fid1*dili*enli)
            acxyf=(fid*dilkb*cnlkb+fid1*dili*cnli)
            axy_nd = axyd_nd*axyf 
            aexy_nd= axyd_nd*aexyf
            acxy_nd= axyd_nd*acxyf

      dlapi_nd  = t11c(neighc)
 	dlapkb_nd =	t12c(neighc)
 	dlaei_nd  =	t13c(neighc) 
 	dlaekb_nd =	t14c(neighc) 
 	dlaci_nd  =	t15c(neighc) 
 	dlackb_nd =	t16c(neighc) 
 	dlepi_nd  = t21c(neighc) 
 	dlepkb_nd =	t22c(neighc) 
 	dleei_nd  =	t23c(neighc) 
 	dleekb_nd =	t24c(neighc) 
 	dleci_nd  =	t25c(neighc)       
 	dleckb_nd =	t26c(neighc)       
 	dlcpi_nd  =	t31c(neighc)   
 	dlcpkb_nd =	t32c(neighc) 
 	dlcei_nd  =	t33c(neighc) 
 	dlcekb_nd =	t34c(neighc) 
 	dlcci_nd  =	t35c(neighc)
 	dlcckb_nd =	t36c(neighc)
      continue
      endif
      return
      end