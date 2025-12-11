       subroutine stream_tube(i,iflg,icnl,nd,neq,ncon,a_axy,a_vxy,cord)
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
!D1 To load array for stream tube model.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/stream_tube.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:32   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:46:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
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
!D4 iflg - flag to implement procedures(dpdp,vapor,etc)
!D4 nd - first dimension of cord and ggg arrays
!D4 neq - number of grid points
!D4 ncon - connectivity array
!D4 a_axy - mass fluxes for liquid phase
!D4 a_vxy - mass fluxes for vapor phase
!D4 cord - coordinate array
!D4 ggg - array (output) for stream tube calculations
!D4 
!***********************************************************************

      use comsptr
      use comsk

      implicit none
      integer iflg,neq,nd,ncon(*)
      integer i,i1,i2,jj,kb,neqp1,ipos
      integer i3,i33, upper_limit, icnl
      integer connect_flag
      real*8 a_axy(*),a_vxy(*),cord(nd,*)
      real*8 tol_c,axy,x33

      real*8 epsilon

      parameter(tol_c=1.d-20)
c
      epsilon=1.e-10

      if(iflg.eq.1) then
c single porosity

         if(icnl.eq.0) then
            upper_limit = 3
         else
            upper_limit = 2
         end if
         neqp1=neq+1
         
         do i3=-3,3
            ggg(i,i3)=0.0d00
         enddo
         
c.......................................................
c     feb 21, 02 s kelkar, commenting out the next if statement because
c     want to include pumping nodes, thus ggg is calculated for
c     all 3 cases 
c     irray(i,o)    < 0   for pumping nodes
c     =0    for OMR nodes
c     =i >0    for regualr nodes
c     if(irray(i,0).ge.0) then
c.... s kelkar 11/7/01 changed the above if from (gt.0) to
c     (ge.0) because want to do load gg array for OMR nodes as well
c.......................................................
         
         i1=ncon(i)+1
         i2=ncon(i+1)
         do jj=i1,i2
            ipos= jj-neqp1
            axy=a_axy(ipos)
            kb=ncon(jj)
            
c     Do loop filters out any connections that aren't 
c     oriented only in the x, y, or z directions
            
            connect_flag = 0
            do i3=1,upper_limit
               x33=cord(kb,i3)-cord(i,i3)
               if(x33.gt. tol_c ) connect_flag = connect_flag + 1
               if(x33.lt.-tol_c ) connect_flag = connect_flag + 1
            end do
            if(connect_flag.lt.2) then
               i33=0
               do i3=1,upper_limit
                  x33=cord(kb,i3)-cord(i,i3)
                  if(x33.gt. tol_c ) i33= i3
                  if(x33.lt.-tol_c ) i33=-i3
               enddo
               ggg(i,i33)=axy
            end if
            
         enddo
c......................................................
c     feb 21, 02 s kelkar, commenting out the next if statement because
c     want to include pumping nodes, thus ggg is calculated for
c     all 3 cases 
c     end if
c..........................................................
         
      else if(iflg.eq.2) then
c     double porosity(not implemented)
      endif
      return
      end
