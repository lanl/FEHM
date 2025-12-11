      subroutine rcm(root,ncon,mask,perm,ccsize,deg)
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
!D1 RCM ordering.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/rcm.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:44   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:02   pvcs
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

      integer ncon(*),deg(*),mask(*),perm(*)
      integer ccsize,fnbr,i,j,jstop
      integer jstrt,k,l,lbegin,lnbr,lperm
      integer lvlend,nbr,node,root
c
      call degree(root,ncon,mask,deg,ccsize,perm)
      mask(root)=0
      if(ccsize.le.1) return
      lvlend=0
      lnbr=1
100   lbegin=lvlend+1
      lvlend=lnbr
      do 600 i=lbegin,lvlend
       node=perm(i)
       jstrt=ncon(node)+1
       jstop=ncon(node+1)
       fnbr=lnbr+1
       do 200 j=jstrt,jstop
        nbr=ncon(j)
        if(mask(nbr).eq.0) go to 200
         lnbr=lnbr+1
         mask(nbr)=0
         perm(lnbr)=nbr
200    continue
       if(fnbr.ge.lnbr) go to 600
        k=fnbr
300     l=k
        k=k+1
        nbr=perm(k)
400     if(l.lt.fnbr) go to 500
         lperm=perm(l)
         if(deg(lperm).le.deg(nbr)) go to 500
          perm(l+1)=lperm
          l=l-1
          go to 400
500      perm(l+1)=nbr
         if(k.lt.lnbr) go to 300
600    continue
       if(lnbr.gt.lvlend) go to 100
       k=ccsize/2
       l=ccsize
       do 700 i=1,k
        lperm=perm(l)
        perm(l)=perm(i)
        perm(i)=lperm
        l=l-1
700    continue
       return
       end
