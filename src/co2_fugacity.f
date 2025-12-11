      subroutine co2_fugacity(var1,var2,var3,prop,der1,der2)
!***********************************************************************
! Copyright 2011. Los Alamos National Security, LLC.  This material was
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

c      use comclath
      implicit none
      real*8 denc,tc,rg,pc,rho3, mol, mco2
      real*8 t,p,rho,dddt,dddp,ent,dhdt,dhdp,visc,dvdt,dvdp
      real*8 del,tau,dddp1,dddt1,fir,fg,dfgdp,fg1,dfgdt,dfgdd
      real*8 vpartial, rho_h2o, Cs, xco2, rho_wco2,dhdd
      real*8 fiot,fiott,fird,firt,firdt,firdd,firtt,ftau,fdel
      real*8 n(42),ti(40),c(40),d(40),a(8),phic(8)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 var1,var2,prop,der1,der2,var3,var4
      real*8 var5,dum2,der3,der4
      real*8 a1,a2,a3,a4,a5,t1,t2,t3,t4,t5,sum,sum2,t11,t21  
      real*8 a6,x,x2,x3,x4,x5,fn,dfn
      integer i, i1, iflg,iphase
      common/crit/ denc,tc,rg,pc
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params2/a,phic
      common/params3/aco2,bco2,capa,capb,capc,capd
c     equation 5.2 pg 1091
c     units are P : MPa
c     T : degree celsius
c     h(enthalpy) : MJ/Kg
c     dhdt : MJ/kg/C
c     dhdp : MJ/Kg/MPa
c     rho(density) : kg/m3
c     dddt : kg/m3/C
c     dddp : kg/m3/MPa
c     visc(viscosity) : Pa.s
c     dvdt: Pa.s/C
c     dvdp : Pa.s/Mpa
      include 'params_eosco2.h'
c calculate the dissolved CO2 concentration based on Duan (2003) eqn.
c need to calculate fugacity first based on the density.
	del = var3/denc
	tau = tc/var2

      call phir(fir,del,tau)
	call dphirddel(fird,del,tau)
	call dphirdddel(firdd,del,tau)
	call dphirdtau(firt,del,tau)
	call dphirddeldtau(firdt,del,tau)

	ftau = (del*del*firdt)-(1000.d0*var1/(denc*rg*tc))
        fdel = (2.d0*del*fird)+(del*del*firdd)+1.
        dddt = (denc*tc/(var2*var2))*(ftau/fdel)/denc

      dddp = 1000.d0/(rg*var2*(1.d0+(del*del*firdd)+(2.d0*del*fird)))

      fg1 = dexp(fir+(del*fird)-dlog(1.d0+(del*fird)))
      dfgdp = (fird*(1.d0+(2.d0*del*fird)+(del*del*firdd))/
     &     (1.d0+(del*fird)))*dddp/denc
      dfgdp = fg1*dfgdp
      
      dfgdt =fg1*(firt+(del*firt*fird)+(del*del*firdt*fird))/
     & (1.d0+(del*fird))
      dfgdd =fg1*(fird+(2.d0*del*fird*fird)+(del*del*fird*firdd))/
     & (1.d0+(del*fird))
      dfgdt = (dfgdt*(-tc/(var2*var2)))+(dfgdd*dddt)

	prop = fg1
	der1 = dfgdp*0.1d0
	der2 = dfgdt      
      end


      subroutine dphiodtau(dr,del2,tau2)
      real*8 a(8),phic(8),del2,tau2,dr,derti_helm
      real*8 ideal_helm
      common/params2/a,phic

      ideal_helm = dlog(del2)+a(1)+(a(2)*tau2)+(a(3)*dlog(tau2))
      do i = 4, 8
         ideal_helm = ideal_helm+(a(i)*dlog(1.d0-dexp(-phic(i)*tau2)))
      enddo

c     Table 34 derivative of ideal helmholtz wrt tau

      derti_helm = a(2) + (a(3)/tau2)
c      write(*,*) derti_helm
      do i = 4, 8
       derti_helm = derti_helm+(a(i)*phic(i)*((1.d0/(1.d0-dexp(-phic(i)*
     & tau2)))-1.d0))
      enddo
c      write(*,*) derti_helm
      dr = derti_helm

      end

      subroutine dphiodtautau(dr,del2,tau2)
      real*8 a(8),phic(8),del2,tau2,dr,dihelm_dtautau
      common/params2/a,phic

c     Table 34 double derivative of ideal helmholtz wrt tau

      dihelm_dtautau = -a(3)/(tau2*tau2)
c      write(*,*) dihelm_dtautau
      do i = 4, 8
       dihelm_dtautau = dihelm_dtautau-(a(i)*(phic(i)**2.)*
     & ((1.d0-dexp(-phic(i)*tau2))**(-2.d0))*(dexp(-phic(i)*tau2)))
      enddo
c      write(*,*) dihelm_dtautau
      dr = dihelm_dtautau

      end

      subroutine phir(r_helm,del2,tau2)

      implicit none
      
      real*8 n(42),ti(40),del2,tau2,r_helm,psi1
      real*8 c(40),d(40),psi,capdel,capdel1
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)     
      integer i
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     equation 5.3 for residual part of Helmholtz function

      r_helm = 0.0d0

      do i = 1, 7
         r_helm = r_helm+(n(i)*(del2**d(i))*(tau2**ti(i)))
      enddo

      do i = 8, 34
         r_helm = r_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*dexp
     &        (-del2**c(i))) 
      enddo

      do i = 35, 39
         r_helm = r_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2))-(beta(i)*
     & ((tau2-gamma(i))**2))))
      enddo

      do i = 40, 42
         psi1=psi(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)

         r_helm=r_helm+(n(i)*(capdel1**bco2(i))*del2*psi1)
      enddo
      end

      subroutine dphirddel(dr,del2,tau2)

      implicit none
      integer i
      real*8 n(42),ti(40),del2,tau2,dr,derdr_helm
      real*8 c(40),d(40)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 psi1,psi,capdel1,capdel,dsidd,dpsiddel,ddelbdd,ddelbiddel
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*d(i)*(del2**(d(i)-1.d0))*
     & (tau2**ti(i)))
      enddo
      do i = 8, 34
         derdr_helm=derdr_helm+(n(i)*(del2**(d(i)-1.d0))*(tau2**ti(i))
     &        *dexp(-del2**c(i))*(d(i)-(c(i)*(del2**c(i))))) 
      enddo
      do i = 35, 39
         derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2.0d0))-(beta(i)*
     & ((tau2-gamma(i))**2.0d0)))*
     & ((d(i)/del2)-(2.0d0*alpha(i)*(del2-ipsilon(i)))))
      enddo
      do i = 40, 42
         psi1=psi(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)
         dsidd=dpsiddel(i,del2,tau2)
         ddelbdd=ddelbiddel(i,del2,tau2)
         derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i))*(psi1+
     &        (del2*dsidd)))+(del2*psi1*ddelbdd)))
      enddo
      dr = derdr_helm
      
      end

      subroutine dphirdddel(dpdd,del2,tau2)
      implicit none
      integer i
      real*8 n(42),ti(40),del2,tau2,derdr_helm
      real*8 c(40),d(40),dpdd,psi1,psi,capdel1,capdel
      real*8 dsidd,dpsiddel,d2sidd,d2psiddel2,ddelbdd,ddelbiddel
      real*8 d2delbdd,d2delbiddel2
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)

      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*d(i)*(d(i)-1.d0)*
     &        (del2**(d(i)-2.d0))*(tau2**ti(i)))
      enddo
      do i = 8, 34
         derdr_helm=derdr_helm+(n(i)*(del2**(d(i)-2))*(tau2**ti(i))*dexp
     &        (-del2**c(i))*(((d(i)-(c(i)*(del2**c(i))))*
     &        (d(i)-1.d0-(c(i)*(del2**c(i)))))-((c(i)**2.d0)*
     &        (del2**c(i)))))
      enddo
      do i = 35, 39
         derdr_helm = derdr_helm+(n(i)*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2.0d0))-(beta(i)*
     & ((tau2-gamma(i))**2.0d0)))*
     & ((-2.d0*alpha(i)*(del2**d(i)))+(4.d0*(alpha(i)**2.d0)*
     & (del2**d(i))*((del2-ipsilon(i))**2.d0))-(4.d0*d(i)*
     & alpha(i)*(del2**(d(i)-1))*(del2-ipsilon(i)))+
     & (d(i)*(d(i)-1.d0)*(del2**(d(i)-2.d0)))))
      enddo
      do i = 40, 42
         psi1=psi(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)
         dsidd=dpsiddel(i,del2,tau2)
         d2sidd=d2psiddel2(i,del2,tau2)
         ddelbdd=ddelbiddel(i,del2,tau2)
         d2delbdd=d2delbiddel2(i,del2,tau2)
         derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i))*((2.d0*
     & dsidd)+(del2*d2sidd)))+(2.d0*ddelbdd*(psi1+(del2*dsidd)))
     & +(d2delbdd*del2*psi1)))
      enddo
      dpdd = derdr_helm
      
      end
      
      subroutine dphirdtau(dpdtau,del2,tau2)

      implicit none
      integer i
      real*8 n(42),ti(40),del2,tau2,derdr_helm
      real*8 c(40),d(40),dpdtau
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 psi1,psi,dbidt,ddelbidtau,capdel1,capdel,dsidt
      real*8 dpsidtau
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)*
     & (tau2**(ti(i)-1))*(del2**d(i)))
      enddo
      do i = 8, 34
         derdr_helm=derdr_helm+(n(i)*ti(i)*(del2**d(i))*
     &        (tau2**(ti(i)-1))*dexp(-del2**c(i)))
      enddo
      do i = 35, 39
         derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2.0d0))-(beta(i)*
     & ((tau2-gamma(i))**2.0d0)))*
     &        ((ti(i)/tau2)-(2.0d0*beta(i)*(tau2-gamma(i)))))
      enddo
      do i = 40, 42
         psi1=psi(i,del2,tau2)
         dbidt=ddelbidtau(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)
         dsidt=dpsidtau(i,del2,tau2)
         derdr_helm=derdr_helm+(n(i)*del2*((dbidt*psi1)+((capdel1**
     &        bco2(i))*dsidt)))
      enddo
      dpdtau = derdr_helm
      
      end

      subroutine dphirdtautau(dpdtt,del2,tau2)

      implicit none
      integer i
      real*8 n(42),ti(40),del2,tau2,derdr_helm
      real*8 c(40),d(40),dpdtt
      real*8 psi,psi1,d2delbidtau2,ddelbidtau,dpsidtau,capdel
      real*8 d2psidtau2,d2bidtt,capdel1,dbidt,dsidt,d2sidtt
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)*(ti(i)-1.d0)*
     &        (tau2**(ti(i)-2.d0))*(del2**d(i)))
      enddo
      do i = 8, 34
         derdr_helm=derdr_helm+(n(i)*ti(i)*(ti(i)-1.d0)*(del2**d(i))*
     &        (tau2**(ti(i)-2))*dexp(-del2**c(i)))
      enddo
      do i = 35, 39
         derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2.0d0))-(beta(i)*
     & ((tau2-gamma(i))**2.0d0)))*(
     & (((ti(i)/tau2)-(2.0d0*beta(i)*(tau2-gamma(i))))**2.d0)
     & -(ti(i)/(tau2*tau2))-(2.d0*beta(i))))
      enddo
      do i = 40, 42
         psi1=psi(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)
         d2bidtt=d2delbidtau2(i,del2,tau2)
         dbidt=ddelbidtau(i,del2,tau2)
         dsidt=dpsidtau(i,del2,tau2)
         d2sidtt=d2psidtau2(i,del2,tau2)
         derdr_helm=derdr_helm+(n(i)*del2*((d2bidtt*psi1)+
     & (2.d0*dbidt*dsidt)+((capdel1**bco2(i))*d2sidtt)))

      enddo
      dpdtt = derdr_helm
      
      end

      subroutine dphirddeldtau(dpddt,del2,tau2)

      implicit none
      integer i
      real*8 n(42),ti(40),del2,tau2,derdr_helm
      real*8 c(40),d(40),dpddt
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 psi,capdel,d2psiddeltau,ddelbiddel,dpsidtau
      real*8 ddelbidtau,dpsiddel,d2delbiddeltau
      real*8 psi1,capdel1,dsidt,dsiddt,dbidd,dbidt,dsidd,d2biddt
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon
      common/params3/aco2,bco2,capa,capb,capc,capd

c     table 36, derivative of residual helmholtz wrt delta

      derdr_helm = 0.0d0
      do i = 1, 7
        derdr_helm=derdr_helm+(n(i)*ti(i)*d(i)*
     &        (tau2**(ti(i)-1))*(del2**(d(i)-1)))
      enddo
      do i = 8, 34
         derdr_helm=derdr_helm+(n(i)*ti(i)*(del2**(d(i)-1))*
     &        (tau2**(ti(i)-1))*dexp(-del2**c(i))*
     & (d(i)-(c(i)*(del2**c(i)))))
      enddo
      do i = 35, 39
         derdr_helm = derdr_helm+(n(i)*(del2**d(i))*(tau2**ti(i))*
     & dexp((-alpha(i)*((del2-ipsilon(i))**2.0d0))-(beta(i)*
     & ((tau2-gamma(i))**2.0d0)))*
     & ((ti(i)/tau2)-(2.0d0*beta(i)*(tau2-gamma(i))))*
     & ((d(i)/del2)-(2.d0*alpha(i)*(del2-ipsilon(i)))))
      enddo
      do i = 40, 42
         psi1=psi(i,del2,tau2)
         capdel1=capdel(i,del2,tau2)
         dsidt=dpsidtau(i,del2,tau2)
         dsiddt=d2psiddeltau(i,del2,tau2)
         dbidd=ddelbiddel(i,del2,tau2)
         dbidt=ddelbidtau(i,del2,tau2)
         dsidd=dpsiddel(i,del2,tau2)
         d2biddt=d2delbiddeltau(i,del2,tau2)
         derdr_helm=derdr_helm+(n(i)*(((capdel1**bco2(i))*(dsidt+
     & (del2*dsiddt)))+(del2*dbidd*dsidt)+(dbidt
     & *(psi1+(del2*dsidd)))+(d2biddt*del2*psi1)))

      enddo
      dpddt = derdr_helm
      
      end

      double precision function psi(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi=-(capc(i)*((del2-1.d0)**2.d0))-(capd(i)*((tau2-1.d0)**2.d0))
      psi=dexp(psi)

      end

      double precision function dpsiddel(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,psi1,psi
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi1=psi(i,del2,tau2)
      dpsiddel = -2.d0*capc(i)*(del2-1.d0)*psi1

      end

      double precision function d2psiddel2(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,psi1,psi
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi1=psi(i,del2,tau2)
      d2psiddel2 = ((2.d0*capc(i)*((del2-1.d0)**2.d0))-1.d0)*2.d0*
     & capc(i)*psi1

      end

      double precision function dpsidtau(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,psi1,psi
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi1=psi(i,del2,tau2)
      dpsidtau = -2.d0*capd(i)*(tau2-1.d0)*psi1

      end

      double precision function d2psidtau2(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,psi1,psi
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi1=psi(i,del2,tau2)
      d2psidtau2 = ((2.d0*capd(i)*(tau2-1.d0)**2.d0)-1.d0)*2.d0*
     & capd(i)*psi1

      end

      double precision function d2psiddeltau(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,psi1,psi
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      psi1=psi(i,del2,tau2)
      d2psiddeltau = 4.d0*capc(i)*capd(i)*(del2-1.d0)*(tau2-1.d0)*psi1

      end

      double precision function theta(i,del2,tau2)
      implicit none
      real*8 n(42),ti(40),c(40),d(40)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon

      theta=(1.d0-tau2)+(capa(i)*(((del2-1.d0)**2.d0)**(1.d0/(2.d0*
     &     beta(i)))))

      end

      double precision function capdel(i,del2,tau2)
      implicit none
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 theta1,theta,del2,tau2
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd

      theta1=theta(i,del2,tau2)
      capdel=(theta1*theta1)+(capb(i)*(((del2-1.d0)**2.d0)**
     &        aco2(i)))
      end

      double precision function dcapdelddel(i,del2,tau2)
      implicit none
      real*8 del2,tau2,theta,theta1
      real*8 n(42),ti(40),c(40),d(40)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon

      theta1=theta(i,del2,tau2)
      dcapdelddel=(del2-1.d0)*((capa(i)*theta1*(2.d0/beta(i))
     & *((del2-1.d0)**2.d0
     & )**((1/(2.d0*beta(i)))-1.d0))+(2.d0*capb(i)*aco2(i)*
     & (((del2-1.d0)**2.d0)**(aco2(i)-1.d0))))

      end

      double precision function d2capdelddel2(i,del2,tau2)
      implicit none
      real*8 tmp1,del2,tau2,theta,theta1
      real*8 n(42),ti(40),c(40),d(40)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 dcapdelddel,ddd
      integer i
      common/params3/aco2,bco2,capa,capb,capc,capd
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon

      theta1=theta(i,del2,tau2)
      tmp1=4.d0*capb(i)*aco2(i)*(aco2(i)-1.d0)*(((del2-1.d0)
     & **2.d0)**(aco2(i)-2.d0))

      tmp1=tmp1+(2.d0*capa(i)*capa(i)*((1.d0/beta(i))**2.d0)*
     & ((((del2-1.d0)**2.d0)**((1.d0/(2.d0*beta(i)))-1.d0))**2.d0))

      tmp1=tmp1+capa(i)*theta1*(4.d0/beta(i))*((1.d0/(2.d0*beta(i)))
     & -1.d0)*(((del2-1.d0)**2.d0)**((1.d0/(2.d0*beta(i)))-2.d0))

      tmp1=tmp1*((del2-1.d0)**2.d0)

      ddd=dcapdelddel(i,del2,tau2)
      d2capdelddel2=tmp1+((1.d0/(del2-1.d0))*ddd)

      end

      double precision function ddelbiddel(i,del2,tau2)
      implicit none
      integer i
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 del2,tau2,capdel1,capdel,ddd,dcapdelddel
      common/params3/aco2,bco2,capa,capb,capc,capd
      
      capdel1=capdel(i,del2,tau2)
      ddd=dcapdelddel(i,del2,tau2)
      ddelbiddel=bco2(i)*(capdel1**(bco2(i)-1.d0))*ddd

      end

      double precision function d2delbiddel2(i,del2,tau2)
      implicit none
      integer i
      real*8 del2,tau2,ddd1,ddd2,capdel1
      real*8 capdel,dcapdelddel,d2capdelddel2
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params3/aco2,bco2,capa,capb,capc,capd

      ddd1=dcapdelddel(i,del2,tau2)
      ddd2=d2capdelddel2(i,del2,tau2)
      capdel1=capdel(i,del2,tau2)
      d2delbiddel2=bco2(i)*((ddd2*(capdel1**(bco2(i)-1.d0)))
     & +((bco2(i)-1.d0)*(capdel1**(bco2(i)-2.d0))*((ddd1)**2.d0)))

      end

      double precision function ddelbidtau(i,del2,tau2)
      implicit none
      integer i
      real*8 del2,tau2,theta,theta1,capdel,capdel1
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params3/aco2,bco2,capa,capb,capc,capd

      theta1=theta(i,del2,tau2)
      capdel1=capdel(i,del2,tau2)

      ddelbidtau=-2.d0*theta1*bco2(i)*(capdel1**(bco2(i)-1.d0))

      end

      double precision function d2delbidtau2(i,del2,tau2)
      implicit none
      real*8 del2,tau2
      integer i
      real*8 capdel,capdel1,theta,theta1
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      common/params3/aco2,bco2,capa,capb,capc,capd

      capdel1=capdel(i,del2,tau2)
      theta1=theta(i,del2,tau2)
      d2delbidtau2=(2.d0*bco2(i)*(capdel1**(bco2(i)-1.d0)))+(4.d0*
     & (theta1**2.d0)*bco2(i)*(bco2(i)-1.d0)*(capdel1**(bco2(i)-2.d0)))

      end

      double precision function d2delbiddeltau(i,del2,tau2)
      implicit none
      real*8 tmp3,del2,tau2
      integer i
      real*8 n(42),ti(40),c(40),d(40)
      real*8 alpha(39),beta(42),gamma(39),ipsilon(39)
      real*8 aco2(42),bco2(42),capa(42),capb(42),capc(42),capd(42)
      real*8 capdel,capdel1,dcapdelddel,ddd,theta,theta1
      common/params3/aco2,bco2,capa,capb,capc,capd
      common/params/ n,c,d,ti,alpha,beta,gamma,ipsilon

      capdel1=capdel(i,del2,tau2)
      theta1=theta(i,del2,tau2)
      ddd=dcapdelddel(i,del2,tau2)
      tmp3=-capa(i)*bco2(i)*(2.d0/beta(i))*(capdel1**(bco2(i)-1.d0))*
     & (del2-1.d0)*(((del2-1.d0)**2.d0)**((1.d0/(2.d0*beta(i)))-1.d0))

      tmp3=tmp3-(2.d0*theta1*bco2(i)*(bco2(i)-1.d0)*(capdel1**
     & (bco2(i)-2.d0))*ddd)
    
      d2delbiddeltau=tmp3

      end

