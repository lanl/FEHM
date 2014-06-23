      subroutine dispersion_node(current_node, d)
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
!***********************************************************************
!D1 
!D1 PURPOSE
!D1 
!D1 Find dispersivity values needed for POD basis functions.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 Initial implementation: 09-SEP-04, Programmer: B. Robinson
!D2
!D2 $Log:   /pvcs.config/fehm90/src/dispersion_node.f_a  $
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!***********************************************************************

      use comai, only : ierr
      use comdi
      use comsptr
      use compart

      implicit none
      integer current_model, current_node
      real*8 d(3),ux,uy,uz,uxux,uyuy,uzuz,uxuz,uxuy,uyuz,uu,u
      real*8 dm,al,ath,atv,a1,a2,a3,a4,asx,asy,asz,sig,sig2
      real*8 ap,bp,cp,dp,term,eig1,eig2,eig3,root2,beta,sbeta
      real*8 alh,alv,at
      
      current_model = itrc(current_node)
      d = 0.
      if(tprpflag(current_model).eq.2.or.
     $     tprpflag(current_model).eq.4) then
         
         
c     initialize random walk displacement to zero.
c     unless isotropic_flag = 0,1,2,4,5,6;  on return dx=0
c     and x2-in = x2-out
         
         dm=dispersivity1(current_model)
         
         
         ux=(ggg(current_node,1)-ggg(current_node,-1))*.5
         uy=(ggg(current_node,2)-ggg(current_node,-2))*.5
         uz=(ggg(current_node,3)-ggg(current_node,-3))*.5
         
         uxux=ux*ux
         uyuy=uy*uy
         uzuz=uz*uz
         uxuz=ux*uz
         uxuy=ux*uy
         uyuz=uy*uz
         
         uu=uxux+uyuy+uzuz
         u = sqrt(uu)
         if(u.gt.1.d-30) then
            
            
c.......................................
            
c     pre-dec-2000 way of doing the random tensor
            
            if(abs(itensor).eq.5) then
               
c     
               al =dispersivity1(current_model)
               ath=dispersivity2(current_model) 
               atv=dispersivity3(current_model)
               dm =dispersivity4(current_model)
               
               d(1)=al*u+dm
               d(2)=(ath*(uxux+uyuy)+atv*uzuz)/u+dm
               d(3)=atv*u+dm
               
               
            elseif(abs(itensor).eq.1) then
c........................
               
c     generalized form of the dispersion tensor, axisymmetric 
c     formulation, 
c     assumed that axis of symmetry oblique to verticle and to the flow
c     Lichtner et al 99, eq 41-47 and my notes dated 10/25/99
c     here a1,a2,a3,a4 are teh 4 generatized dispersivities from eq 41
c     and asx, asy,asz are the direction cosins of the axis of symmetry 
c     denoted by
c     Lambda in eq41
c     vxlambda is an eigen vector of D, the other two are constructed 
c     from linear
c     combination of v and lambda
c     the columns of U (denoted here by r(i,j)) are made up of the
c     eigenvectors. 
c     vas = unit vector in velocity direction. dot. axis os symmetry
               
               a1=dispersivity2(current_model)
               a2=dispersivity3(current_model)
               a3=dispersivity4(current_model)
               a4=dispersivity6(current_model)
               asx=dispersivity7(current_model)
               asy=dispersivity8(current_model)
               asz=sqrt(1.-asx*asx-asy*asy)
               sig=(asx*ux+asy*uy+asz*uz)/u
               sig2=sig*sig
               if((1.-sig2).lt.1.e-20) then
c     v parallel to axis of symmetry, requires a different formulation
                  write(ierr,*)'v parallel to axis of symmetry, '
                  write(ierr,*)' requires a different formulation. stop'
                  stop
               else
                  ap=a1+a2+sig2*a3+sig*a4
                  bp=sig*a3+0.5*a4
                  cp=(1.-sig2)*bp
                  dp=a1+(1.-sig2)*a3
                  
                  term=(ap-dp)**2.+4.*bp*cp
                  if(term.lt.0.) then
                     write(ierr,*)'error, term.le.0 in random_walk.'
                     write(ierr,*)'check dispersiv coeff in sptr macro.'
                     write(ierr,*)'STOP'
                     stop
                  endif
                  term=sqrt(term)
                  
                  eig1=0.5*((ap+dp)+term)
                  eig2=0.5*((ap+dp)-term)
                  
                  root2=-cp/(ap-eig2)
c     11/18/99 s kelkar NOTE: refere to eq 37-40 in lochtner et al
c     writing xsi-2=(root2*v-hat+omega-hat)/norm but
c     xsi-1=(v-hat+inv.root1*omega-hat)/norm
c     in order to facilitate the limit of the axissymetrc case
c     with alpha-3 = alpha-4 =0.
c     also note that the equations bellow are cast in terms of lamda-hat 
c     instead of omega-hat
                  
                  eig3=a1
                  
                  d(1)=eig1*u
                  d(2)=eig2*u
                  d(3)=eig3*u
               end if   
            elseif(abs(itensor).eq.2) then
                  
c..............................
               
c     Burnett and Frind tensor, axisymmetric 
c     formulation, 
c     assumed that axis of symmetry is along z axis.
c     Lichtner et al 99, eq 59-62 and 76-78
               
               
               beta=(uxux+uyuy)
               
               if(beta.ge.1.e-34) then
                  sbeta=sqrt(beta)
                  sig=uz/u
                  sig2=sig*sig
                  
               endif
               
               al  =dispersivity2(current_model)
               ath =dispersivity3(current_model) 
               atv =dispersivity4(current_model)
               
               a1=ath+sig2*(atv-ath)
               
               
               d(1)=(al*u+dm)
               d(3)=(atv*u+dm)
               d(2)=(a1*u+dm)
               
            elseif(abs(itensor).eq.3) then
               
c................................................
               
c     modified Burnett and friend, Lichtner et al eq 57-61,
c     formulation, Lamda assumed to be along z axis 
c     assumed that axis of symmetry oblique to verticle and to the flow
c     here a1,a2,a3,a4 are teh 4 generatized dispersivities from eq 41
c     vxlambda is an eigen vector of D, the other two are constructed 
c     from linear
c     combination of v and lambda
c     the columns of U (denoted here by r(i,j)) are made up of the
c     eigenvectors. 
c     vas = unit vector in velocity direction. dot. axis os symmetry
               
               alh=dispersivity2(current_model)
               alv=dispersivity3(current_model)
               ath=dispersivity4(current_model)
               atv=dispersivity6(current_model)
               
               asx=0.
               asy=0.
               asz=1.
               
               beta=(uxux+uyuy)
               sig=(asx*ux+asy*uy+asz*uz)/u
               sig2=sig*sig
               
               if((1.-sig2).lt.1.e-20) then
c     v parallel to axis of symmetry, requires a different formulation
                  
                  d(1)=(alv*u+dm)
                  d(2)=(ath*u+dm)
                  d(3)=(ath*u+dm)
                  
               else
                  
                  al=alh+sig2*(alv-alh)
                  at=atv+sig2*(ath-atv)
                  
                  sbeta=sqrt(beta)
                  
                  d(1)=(al*u+dm)
                  d(2)=(at*u+dm)
                  d(3)=(ath*u+dm)
                  
               endif
               
            elseif(abs(itensor).eq.4) then
               
c..............................
               
c     isotropic formulation of the tensor as in thompson's report
               
               al=dispersivity2(current_model)
               at=dispersivity3(current_model) 
               
               
c     Looks wrong - BAR 9-9-2004
               d(1) = (al*u+dm)
               d(2) = (al*u+dm)
               d(3) = (al*u+dm)
               
            end if
c     End of itensor options, this else is for
c     no velocity (i.e. diffusion only)
         else
c     if velocity is zero, only molecular dispersion
            
            
            d(1) = dm
            d(2) = dm
            d(3) = dm
            
         endif
         
      end if
      
      return
      
      end

