      subroutine simplify_ncon
     &(iflg,ncon,nelmdg,idum,istrw,idum1,neq,idof,ka,ps,icount)      
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
!D1 To simplify connectivity based on porosity for non heat conduction 
!D! problems.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/simplify_ncon.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:56   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 09:15:54   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of startup 08-Feb-02
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 None
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 None
!D4 
!**********************************************************************
c gaz 4-23-2001 made porosity condition strictly lt 0
c gaz 5-27-2002 icount on exit is number of removed nodes

      implicit none

      integer iflg, idof, neq, neqp1
      integer i,jj,i1,i2,nmat,icount
      integer ka(*),ncon(*),idum(*),nelmdg(*), idum1(*), istrw(*)
      real*8 ps(*)
     
      icount = 0
      if(idof.le.1) return
      if(iflg.eq.0) then    
       neqp1=neq+1
       nmat=ncon(neqp1)
       do i=1,nmat
        idum(i)=ncon(i)
       enddo
       do i=1,nmat-neqp1
        idum1(i)=istrw(i)
       enddo
       icount=neqp1
       ncon(1)=neqp1
       do i=1,neq
        if(ps(i).ge.0.0) then
         i1=idum(i)+1
         i2=idum(i+1)
         do jj=i1,i2
         if(ps(idum(jj)).ge.0.0) then
           icount=icount+1
           ncon(icount)=idum(jj)
           istrw(icount-neqp1)=idum1(jj-neqp1)
           if(idum(jj).eq.i) nelmdg(i)=icount
         endif
         enddo
         ncon(i+1)=icount
        else
         icount=icount+1
         ncon(i+1)=icount
         nelmdg(i)=icount
         ncon(icount)=i   
        endif
       enddo
       icount = 0
       do i =1,neq
        if(ps(i).lt.0.0) then
          icount = icount+1
          ps(i) = 0.0
        endif
       enddo
      else if(iflg.eq.1) then    
c
      endif
      return
      end
