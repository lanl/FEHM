      subroutine wellimped_ctr(iflg)
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
!!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To create a wellbore resistence term that represents a "peaceman" well 
!D1 and other type subgrid scale resistance terns
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 9-Jan-02, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/submodel_bc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai

      implicit none

      integer, allocatable :: kq_dum(:), icount(:)
      integer i,j,ii,ij,jj,kb,i1,i2,neqp1,max_wellmod
      integer izone1,izone2,iflg,ibnd,iroot,idsubm 
      integer k,open_file
c gaz 110819 removed tref, pref (now global)       
      real*8 aiped, pnxa, pnya 
      real*8 head_value, parm1, parm2, parm3
      real*8 delx,dely,delz,xi,yi,zi 
      real*8 rw,rthick,r0,rw_tol
      real*8 peaceman_term,aperm,term1,term2
      real*8 pi2
      logical null1
      character*4 keyword
      parameter(max_wellmod = 50,pi2 = 6.283185,rw_tol = 1.e-6)
c

         if(iflg.eq.0) then

            if (.not. allocated(izonewel1)) then
c 
c this should be the fist time called
c
c 
               allocate (izonewel1(n0))
               allocate (itypwel(max_wellmod))         
               allocate (keywel1(max_wellmod),keywel2(max_wellmod))
               allocate (keywel3(max_wellmod),keywel4(max_wellmod))
               allocate (parmwel1(max_wellmod),parmwel2(max_wellmod))
               allocate (parmwel3(max_wellmod))
               allocate (icountw(max_wellmod))
               itypwel = 0
               izonewel1 = 0
               keywel1 = ''
	         keywel2 = ''
	         keywel3 = ''
               keywel4 = ''
               isubwd = 0
             endif
c model 1            
c  parmwel1 = peaceman's rw wellbore radius
c  parmwel2 = peaceman's h reservoir thickness        
c  parmwel3 = print parameter
c  
c          
c        read input parameters
c
            do
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) exit      
               isubwd = isubwd + 1
               read(wdd1,*) itypwel(isubwd)
             if(itypwel(isubwd).eq.1) then  
               read(wdd1,*) itypwel(isubwd),keywel1(isubwd),
     &              parmwel1(isubwd),parmwel2(isubwd),parmwel3(isubwd)
             else if(itypwel(isubwd).eq.2) then 
               read(wdd1,*) itypwel(isubwd),keywel1(isubwd),
     &              parmwel1(isubwd),parmwel2(isubwd),parmwel3(isubwd)
             endif
            end do
            
            igroup = 1
            narrays = 1
            itype(1) = 4
            igroup = 1
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           i4_1 = izonewel1(1:n0)) 

c
c  adjust izonewel1 for multiple calls 
c


         else if(iflg.eq.1) then
c
c  adjust impedance or connection to represnt a variety of 
c  near-wellbore analytic solutions
c
c  first check for applied impedance model that is not associated with a source or sink
c gaz 061210 
c          do i = 1,n0
c           if(ka(i).eq.0) izonewel1(i) = 0
c          enddo
c          
          do i = 1,n0
c ij is the model number (in numerical sequence)          
           ij = izonewel1(i)
c jj is the ijth model type (1 = peaceman, 2 = horizontal,etc) 
           jj = 0
           if(ij.ne.0) then        
            jj = itypwel(ij)
           endif
	     if(jj.eq.1) then
c
c Peaceman model 
c
             delx = 0.0
             dely = 0.0
             delz = 0.0
c find gridblock dimensions
             xi = cord(i,1)
             yi = cord(i,2)
             zi = cord(i,3)
             i1 = nelm(i)+1
             i2 = nelm(i+1)
             do k = i1,i2
              kb = nelm(k)
              delx = max(abs(xi-cord(kb,1)),delx)
              dely = max(abs(yi-cord(kb,2)),dely)
              delz = max(abs(zi-cord(kb,3)),delz)
             enddo
c delx,dely,delz represent the max diff between neighbors (upper bound of gridblock size)              
             if(keywel1(jj).eq.'pwxy') then          
               rw = parmwel1(ij)    
               if(rw.lt.rw_tol) then
                if(iout.ne.0) write(iout,*) 
     &       'wellbore radius too small (wellimped_ctr.f),stopping'
                if(iptty.ne.0) write(iptty,*) 
     &       'wellbore radius too small (wellimped_ctr.f),stopping'
                stop
               endif        
               rthick = parmwel2(ij)
   	         pnxa = 1.e-6*pnx(i)
   	         pnya = 1.e-6*pny(i)
   	         term1 = (pnya/pnxa)**0.25 + (pnxa/pnya)**0.25
   	         term2 = sqrt(sqrt(pnya/pnxa)*delx*delx + 
     &                sqrt(pnxa/pnya)*dely*dely) 
    	         r0 = 0.28*term2/term1	         
   	         aperm = sqrt(pnxa*pnya)
   	         peaceman_term = 1.e06*pi2*aperm*rthick/log(r0/rw)
c
c  note: this well impedance still needs to be divided by dil(i)
c  	         
   	         wellim(i) = peaceman_term
   	       endif
   	    
	     else if(jj.eq.2) then
c
c horizontal well
c	     
	     endif 
          enddo       
      else if(iflg.eq.2) then
c
c pre-factor well solution
c        
      else if(iflg.eq.3) then  
c
c back substitute well solution
c 
      else if(iflg.eq.4) then  
c
c printout well solution
c  
         if(l.eq.1) then
               welifile = 
     &           open_file('wellimpef.welli', 'unknown')   
               write(welifile,*)
     &          'this file contains the wellbore pressure,flowrate,time'
         endif      
         do i = 1,n0
c ij is the model number (in numerical sequence)          
           ij = izonewel1(i)
c jj is the ijth model type (1 = peaceman, 2 = horizontal,etc)    
           jj = 0
           if(ij.ne.0) then     
            jj = itypwel(ij)
           endif
           if(jj.eq.1) then
            if(parmwel3(ij).ne.0.0) then
             write(welifile,100) i,ij,jj,days,phi(i),pflow(i),sk(i)
            endif
           else
           endif
         enddo
                                 
      endif
100   format(1x,i9,1x,i3,1x,i3,1p,3(1x,3g14.5))      
      return 
      end
