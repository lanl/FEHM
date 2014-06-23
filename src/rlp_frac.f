      subroutine rlp_frac(iflg,i,sl,star,rp1,rp2,
     &                    rl,drls,rv,drvs,rp24)
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
!D1 To calculate fracture - matrix interaction terms.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2
!D2 $Log:   /pvcs.config/fehm90/src/rlp_frac.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:14:14   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of rlperm 08-Feb-02
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.4.4 Relative-permeability and capillary-pressure functions
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

      use comdti
      use comai
      use combi
      use comci
      use comdi
      implicit none

      integer iflg
      integer icode
      integer i
      real*8 sl
      real*8 star
      real*8 rv 
      real*8 drvs 
      real*8 rl 
      real*8 drls 
      real*8 rp1
      real*8 rp2
      real*8 rp24

      if(iflg.eq.0) then
c error condition       
         if(idpdp.eq.0.and.idualp.eq.0) then
            write(ierr,100)
            if (iptty .ne. 0) write(iptty,100)
            if (iout .ne. 0) write(iout,100)
            stop
         endif
         if(rp24.gt.0.0.and.rp24.lt.1.0) then
            write(ierr,100)
            if (iptty .ne. 0) write(iptty,100)
            if (iout .ne. 0) write(iout,100)
 100        format(/,'*************************************'
     1           ,'f-m terms but no dpdp : stopping'
     2           ,'*************************************')
            stop
         endif
c allocate memory
       if(.not. allocated(fmvf)) allocate(fmvf(n0))
       if(.not. allocated(dfmvf)) allocate(dfmvf(n0))
       if(.not. allocated(fmlf)) allocate(fmlf(n0))
       if(.not. allocated(dfmlf)) allocate(dfmlf(n0))
      else
c form f-m functions
         if(rp24.le.0.0) then
c rl perm is f-m function  
            fmlf(i)=rl
            dfmlf(i)=drls
            fmvf(i)=rv
            dfmvf(i)=drvs
         else
c     S to a power
c     dstar= 1.0/(rp2-rp1)
c     fmlf(i)=(star)**rp24
c     dfmlf(i)=rp24*star**(rp24-1.0)
            if(rp24.eq.1.) then
               fmlf(i)=sl
               dfmlf(i)=1.
               fmvf(i)=1.-sl
               dfmvf(i)=-1.
            else
               fmlf(i)=sl**rp24
               dfmlf(i)=rp24*sl**(rp24-1.0)
               sv=1.-sl
               fmvf(i)=sv**rp24
               dfmvf(i)=-rp24*sv**(rp24-1.0)
            end if
         endif
      endif
      return
      end
