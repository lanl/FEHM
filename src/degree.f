      subroutine degree (root,ncon,mask,deg,ccsize,ls)
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
!D1 To compute the degrees of a node.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/degree.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:06   pvcs
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

      integer ncon(*),deg(*),ls(*),mask(*)
      integer ccsize,i,ideg,j,jstop
      integer jstrt,lbegin,lvlend,lvsize,nbr
      integer node,root
c
      ls(1)=root
      ncon(root)=-ncon(root)
      lvlend=0
      ccsize=1
100   lbegin=lvlend +1
      lvlend=ccsize
      do 400 i=lbegin,lvlend
         node = ls(i)
         jstrt= -ncon(node)+1
         jstop= iabs(ncon(node+1)) 
         ideg = 0
         if(jstop.lt.jstrt) go to 300
         do 200 j=jstrt,jstop
            nbr=ncon(j)
            if(mask(nbr).eq.0) go to 200
            ideg=ideg+1
            if(ncon(nbr).lt.0) go to 200
            ncon(nbr)=-ncon(nbr)
            ccsize=ccsize+1
            ls(ccsize)=nbr
 200     continue
 300     deg(node)=ideg
 400  continue
      lvsize= ccsize-lvlend
      if(lvsize.gt.0) go to 100
      do 500 i =1,ccsize
         node=ls(i)
         ncon(node)=-ncon(node)
 500  continue
      return
      end
