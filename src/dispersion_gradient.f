      subroutine dispersion_gradient
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To calculate the dispersion gradient.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: 25-AUG-99, Programmer: S. Kelkar
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/dispersion_gradient.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:35:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:27:08   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!D4 grad.D term from Tompson and Gelhar, WRR, Oct 90, eq. 8
!D4
!**********************************************************************

c 10/31/00 s kelkar
c grad.D term for the general dispersion tensor

      use comai
      use comdi
      use comsptr
      implicit none


      integer i,j,k,current_model, current_flag

      real*8 dvd(3),dvxd(3),dvyd(3),dvzd(3),v(3),vv
      real*8 a1,a2,a3,a4,alamda(3),dvjdxi(3,3)
      real*8 al, ath, atv, cost, alh, alv, at
      real*8 dvidxi,dvjdxivi,dvdxili,dvdxivi,dvjdxili
      real*8 term1,term2,term3,term4,term5,term6

      do i=1,neq
         
         current_model = itrc(i)
         current_flag = tprpflag(current_model)
         if(current_flag.eq.2.or.current_flag.eq.4) then

            call velocity_derivatives(i,dvd,dvxd,dvyd,dvzd,v,vv)
c NOTE v is the normalized velicity, ie v(1)=Vx/vv etc.
            
c     select coeffecients depending on the type of dispersion tensor
            
            if(itensor.eq.1) then
c     general tensor
               a1=dispersivity2(current_model)
               a2=dispersivity3(current_model)
               a3=dispersivity4(current_model)
               a4=dispersivity6(current_model)
               alamda(1)=dispersivity7(current_model)
               alamda(2)=dispersivity8(current_model)
               term1=alamda(1)*alamda(1)+alamda(2)*alamda(2)
               if(term1.gt.1.0) then
                  write(ierr,*)'error in grad_general. sum of direction'
                  write(ierr,*)'cosines must be le 1. Check sptr input.'
                  write(ierr,*)' STOP'
                  stop
               endif
               alamda(3)=sqrt(1.-term1)
               
            elseif (itensor.eq.2) then
c     Burnett & Frind tensor
               alamda(1)=0.
               alamda(2)=0.
               alamda(3)=1.
               al  =dispersivity2(current_model)
               ath =dispersivity3(current_model) 
               atv =dispersivity4(current_model)
c NOTE v is the normalized velicity, ie v(1)=Vx/vv etc.
               cost=v(3)
               a1=ath+cost*cost*(atv-ath)
               a2=al-ath
               a3=atv-ath
               a4=-2.*cost*a3
               
            elseif(itensor.eq.3) then
c     Modified B&F tensor
               alamda(1)=0.
               alamda(2)=0.
               alamda(3)=1.
               alh =dispersivity2(current_model)
               alv =dispersivity3(current_model)
               ath =dispersivity4(current_model) 
               atv =dispersivity6(current_model)
c NOTE v is the normalized velicity, ie v(1)=Vx/vv etc.
               cost=v(3)
               a1=ath
               a2=alh-ath+cost*cost*(alv-alh+atv-ath)
               a3=atv-ath
               a4=-2.*cost*a3
               
            elseif(itensor.eq.4) then
c     Thompson tensor
               alamda(1)=0.
               alamda(2)=0.
               alamda(3)=0.
               al =dispersivity2(current_model)
               at =dispersivity3(current_model)
               a1=at
               a2=al-at
               a3=0.
               a4=0.
               
            elseif(itensor.eq.5) then
c     old fehm V2.10
               alamda(1)=0.
               alamda(2)=0.
               alamda(3)=0.
               al =dispersivity1(current_model)
               at =dispersivity2(current_model)
               a1=at
               a2=al-at
               a3=0.
               a4=0.
               
            endif
            
            do j=1,3
               dvjdxi(1,j)=dvxd(j)
               dvjdxi(2,j)=dvyd(j)
               dvjdxi(3,j)=dvzd(j)
            enddo
            
            dvidxi=0.
            dvjdxivi=0.
            dvdxili=0.
            dvdxivi=0.
            dvjdxili=0.
            
            do j=1,3
               dvidxi=dvidxi+dvjdxi(j,j)
               dvdxili=dvdxili+dvd(j)*alamda(j)
               dvdxivi=dvdxivi+dvd(j)*v(j)
            enddo
            
            do j=1,3
               
               dvjdxivi=0.
               dvjdxili=0.
               do k=1,3
                  dvjdxivi=dvjdxivi+dvjdxi(j,k)*v(k)
                  dvjdxili=dvjdxili+dvjdxi(j,k)*alamda(k)
               enddo
               
c NOTE v is the normalized velicity, ie v(1)=Vx/vv etc.
               term1=a1*dvd(j)
               term2=(a2*v(j)+0.5*a4*alamda(j))*dvidxi
               term3=a2*dvjdxivi
               term4=a3*alamda(j)*dvdxili
               term5=-a2*v(j)*dvdxivi
               term6=0.5*a4*dvjdxili
               
               gradd(i,j)=term1+term2+term3+term4+term5+term6
               
            enddo
            
         end if         
      enddo
      
      return
      
      end
