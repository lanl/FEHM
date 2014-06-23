      subroutine genrcm(neqns,ncon,perm,mask,xls)
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
!D1 General reversed Cuthill-Mckee algorthm.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/genrcm.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:06:30   pvcs
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

      integer ncon(*),mask(*),perm(*),xls(*)
      integer ccsize,i,neqns, nlvl
      integer num,root
c
      do 100 i=1,neqns
       mask(i)=1
100   continue
      num=1
      do 200 i=1,neqns
       if(mask(i).eq.0) go to 200
        root=i
        call fnroot(root,ncon,mask,nlvl,xls,perm(num))
        call rcm(root,ncon,mask,perm(num),ccsize,xls)
        num=num+ccsize
        if(num.gt.neqns) return
200   continue
      end
