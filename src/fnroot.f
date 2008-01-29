      subroutine fnroot(root,ncon,mask,nlvl,xls,ls)
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
!D1 Subroutine used in RCM algorthm.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/fnroot.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:02   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:05:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2 Split subroutine out of renum 08-Feb-02
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

      implicit none

      integer ncon(*),ls(*),mask(*),xls(*)
      integer ccsize,j,jstrt,k,kstop,kstrt,mindeg,nabor
      integer ndeg,nlvl,node,nunlvl,root
c
      call rootls(root,ncon,mask,nlvl,xls,ls)
      ccsize=xls(nlvl+1) 
      if(nlvl.eq.1.or.nlvl.eq.ccsize) return
c
c gax 1-18-2002
c
      nunlvl= 1000000
100   jstrt=xls(nlvl)
      mindeg=ccsize
      if(ccsize.ne.jstrt) then
       do j=jstrt,ccsize
        node=ls(j)
        ndeg=0
        kstrt=ncon(node)+1
        kstop=ncon(node+1) 
         do k=kstrt,kstop
          nabor=ncon(k)
          if(mask(nabor).gt.0) ndeg=ndeg+1
         enddo
        if(ndeg.lt.mindeg) then
         root=node
         mindeg=ndeg
        endif
       enddo
      endif
c
      call rootls(root,ncon,mask,nlvl,xls,ls)
      if(nunlvl.le.nlvl) return
        nlvl = nunlvl
        if(nlvl.lt.ccsize) go to 100
        return
      end
