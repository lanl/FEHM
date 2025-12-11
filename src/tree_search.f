      subroutine tree_search(inp1,level_max,ierr, iptty,
     $     xp,yp,zp,flag_box,nout_face,save_out_face)
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
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
!D1 To find which control volume a point belongs to.
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.30
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/tree_search.f_a  $
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
c s kelkar Jan 20 05
c starting with the node inp1, do a search on neighbours to find
c which control volume the point xp,yp,zp belongs to. the search
c progresses down successive levels of neighbours of neighbours
      
      use comsk, only : omr_flag

      implicit none

      integer inp1,i,n_oldlist,n_newlist,n1,n2,index
      integer n_level,level_max, maxdim,flag_box,n_total
      integer, allocatable :: oldlist(:)
      integer, allocatable :: newlist(:)
      integer ierr,iptty
      integer nout,iout_face(3),save_out_face(6), nout_face

      real*8 xp,yp,zp

c first do a level 0 search- is the particle in the cc of inp1?
      flag_box=0
      call check_inside_box(inp1,xp,yp,zp,flag_box,nout,
     $     iout_face,nout_face,save_out_face)
c flag_box > 0 :particle inside the box,flag_box=inp1
c flag_box < 0 :particle outside the model domain
c flag_box = 0 :particle outside the box, but may be within model domain

      if(flag_box.eq.0) then
c     continue search
         maxdim=1000000
         allocate(oldlist(maxdim))
         allocate(newlist(maxdim))
         
         n_oldlist=0
         n_newlist=1
         oldlist(1)=inp1
         flag_box=0
         
         do n_level=1,level_max
            
c     going to compress_list, n_oldlist is the #nodes 2 levels, back
c     and n_newlist is the # of nodes 1 level back. On return,
c     n_oldlist is the #nodes 1 level back, and n_newlist is the 
c     #nodes at the current level
            call compress_list(maxdim,n_oldlist,oldlist,n_newlist,
     1           newlist,n_level)
            
            do n1=1,n_newlist
               i=newlist(n1)
               flag_box=0
               call check_inside_box(i,xp,yp,zp,flag_box,nout,
     $              iout_face,nout_face,save_out_face)
c flag_box > 0 :particle inside the box,flag_box=inp1
c flag_box < 0 :particle outside the model domain
c flag_box = 0 :particle outside the box, but may be within model domain
               if(flag_box.ne.0) then
                  goto 99999
               endif
            enddo
            
c     augement oldlist with newlist. Keeping the previous oldlist
c     to go 2 levels back in compress_list
            n_total=n_oldlist+n_newlist
            if(n_total.gt.maxdim) goto 99998
            do n1=1,n_newlist
               oldlist(n1+n_oldlist)=newlist(n1)
            enddo
            
         enddo
c box not found. Set flag_box and return
         flag_box=0
         
99998    continue
c         write(ierr,*)'error in tree_search.maxdim,n_total= ',
c     1        maxdim,n_total
c         write(ierr,*)' inp1,level_max,xp,yp,zp=',inp1,
c     $        level_max,xp,yp,zp
c         write(ierr,*)'STOP'
c         if (iptty .ne. 0) then
c            write(iptty,*)'error in tree_search.maxdim,n_total= ',
c     1           maxdim,n_total
c            write(ierr,*)' inp1,level_max,xp,yp,zp=',inp1,
c     $           level_max,xp,yp,zp
c         write(iptty,*)'STOP'
c         end if
c         stop

99999    deallocate(oldlist)
         deallocate(newlist)

      endif

      return

      end

c..................................................................

      subroutine compress_list(maxdim,n_oldlist,oldlist,n_newlist,
     1     newlist,n_level)

c s kelkar Jan 21-05
c make a list of all the nodes connected to all the nodes in oldlist
c remove repeating node #s from newlist, and also remove the node #s
c that already appeared in oldlist. use a sort routine to help
c going to compress_list, n_oldlist is the #nodes 2 levels, back
c and n_newlist is the # of nodes 1 level back. On return,
c n_oldlist is the #nodes 1 level back, and n_newlist is the 
c #nodes at the current level
c on entry, oldlist(1:n_oldlist) is the list 2 levels back
c and oldlist(1+n_oldlist:n_total) is the list 1 level back
c on return, oldlist(1:n_oldlist) is the list 1 level back
c and newlist(1:n_newlist) is the nodes at current level

      use combi

      implicit none

      integer n_oldlist,n_newlist,newn1,n_level
      integer index,index2,n1,n2,i,ii1,ii2,kb,n_total,maxdim
      integer oldlist(maxdim),newlist(maxdim)


      if(n_level.eq.1) then 

         index=0
         i=oldlist(1)
         ii1=nelm(i)+1
         ii2=nelm(i+1)
         do n2=ii1,ii2
            kb=nelm(n2)
            if(kb.ne.i) then
               index=index+1
               newlist(index)=kb
            endif
         enddo
         index2=index

      else

         n_total=n_oldlist+n_newlist
         index=0
         do n1=1+n_oldlist,n_total
            i=oldlist(n1)
            ii1=nelm(i)+1
            ii2=nelm(i+1)
            do n2=ii1,ii2
               kb=nelm(n2)
               if(kb.ne.i) then
                  index=index+1
                  newlist(index)=kb
               endif
            enddo
         enddo
         
c     sort  newlist in ascending order
         call hpsorti(index,newlist)
         
         index2=1
         do n1=2,index
            newn1=newlist(n1)
            if(newn1.ne.newlist(index2)) then
               
c     oldlist has two ordered segments, these are treated seperately
               
c     check node list from previous level
               if(newn1.gt.oldlist(n_total)) then
c     newn1 is not in list for the previous level, goto checking 
c     the list from 2 levels back
                  goto 99997
               else
                  do n2=1+n_oldlist,n_total
                     if(newn1.eq.oldlist(n2)) then
c     newn1 appeared in previous level. skip and continue with next  
c     element of newlist
                        goto 99999
                     elseif(newn1.lt.oldlist(n2)) then
c     newn1 is not in list for the previous level, goto checking 
c     the list from 2 levels back
                        goto 99997
                     endif
                  enddo
               endif
               
c     check node list from 2 levels back
99997          if(newn1.gt.oldlist(n_oldlist)) then
c     newn1 is new, end check, augment newlist and continue
                  goto 99998
               else
                  do n2=1,n_oldlist
                     if(newn1.eq.oldlist(n2)) then
c     newn1 already appeared. skip and continue with next element of 
c     newlist
                        goto 99999
                     elseif(newn1.lt.oldlist(n2)) then
c     newn1 is new, end check, augment newlist and continue with next
c     element of newlist
                        goto 99998
                     endif
                  enddo
               endif
99998          index2=index2+1
               newlist(index2)=newn1
            endif
99999       continue
         enddo
         
      endif

c     reorganize oldlist so that 1 to n_newlist contain node#s from the
c     previous level
      do n1=1,n_newlist
         oldlist(n1)=oldlist(n_oldlist+n1)
      enddo
c     update n_oldlist and n_newlist
      n_oldlist=n_newlist
      n_newlist=index2
      
      return

      end

c..................................................................

      SUBROUTINE hpsorti(n,ia)
C
C
C#######################################################################
C
C     PURPOSE -
C
C     THIS IS BASED ON THE NUMERICAL RECIPES ROUTINE FOR SORTING
C     AN ARRAY USING THE 'HEAP SORT'.
C
C
C     INPUT ARGUMENTS -
C
C        n - NO. OF ELEMENTS TO SORT INTO ASCENDING ORDER.
C        ia - INTEGER ARRAY TO BE SORTED.
C
C
C
C     OUTPUT ARGUMENTS -
C
C        ia - SORTED INTEGER ARRAY.
C
C
C     CHANGE HISTORY -
C
C        $Log:   /pvcs.config/utilities/src/hpsorti.f_a  $
CPVCS    
CPVCS       Rev 1.0   22 Aug 2000 10:21:48   dcg
CPVCS    Initial revision.
CPVCS    
CPVCS       Rev 1.0   Tue Aug 24 21:58:34 1999   kuprat
CPVCS    Initial revision.
C
C#######################################################################
C
 
      implicit none
 
      INTEGER n
      INTEGER ia(n)
      INTEGER i,ir,j,l
      INTEGER iia
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          iia=ia(l)
        else
          iia=ia(ir)
          ia(ir)=ia(1)
          ir=ir-1
          if(ir.eq.1)then
            ia(1)=iia
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ia(j).lt.ia(j+1))j=j+1
          endif
          if(iia.lt.ia(j))then
            ia(i)=ia(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ia(i)=iia
      goto 10

      END

c..................................................................
      
      subroutine check_inside_box(inp1,xp,yp,zp,flag_box,nout,
     $     iout_face,nout_face,save_out_face)

c s kelkar Jan 21, 05
c check if the point xp,yp,zp is within the brick shape around the 
c node i. If so, flag_box=inp1, if not, then
c if the particle exited model domain, flag_box=-inp1
c if not sure that particle exited model domain, then flag_box=0
      use comsk
      use comsptr

      implicit none
      
      integer flag_box,inp1,nout,iout_face(3),ibou_flag,i
      integer save_out_face(6), nout_face

      real*8 xp,yp,zp,x(3),ddd(3), corner(-3:3)

      flag_box=0
      x(1)=xp
      x(2)=yp
      x(3)=zp
      ddd(1)=ddx(inp1)
      ddd(2)=ddy(inp1)
      ddd(3)=ddz(inp1)
      do i=1,3
         corner(-i)=corn(inp1,i)
         corner(+i)=corn(inp1,i)+ddd(i)
      enddo
     
c     check if particle within box limits
      nout=0
      do i=1,3
         iout_face(i)=0
      enddo
      do i=1,3
         if(corner(-i).gt.x(i)) then
            nout=nout+1
            iout_face(nout)=-i
         endif
         if(corner(+i).lt.x(i)) then
            nout=nout+1
            iout_face(nout)=+i
         endif
      enddo

         flag_box=inp1
         nout_face=0
         if(nout.ge.1) then
c if particle not in box, 
            flag_box=0
            do i=1,nout
               call included_in_ibou(inp1,iout_face(i),ibou_flag)
               if(ibou_flag.eq.1) then
                  flag_box = -inp1
                  nout_face=nout_face+1
                  save_out_face(nout_face)=iout_face(i)
               endif
            enddo
         endif
    
      return

      end

c...................................................................

      subroutine included_in_ibou(inp1,ifc,ibou_flag)
c find if ic is in the ibou_list of node inp1

      use comsptr
      use comsk
      
      integer inp1,ifc,ibou_flag,ibou,i

      ibou_flag=0

      if(omr_flag) then
         ibou=iboulist(inp1,7)
         do i=1,ibou
            if(ifc.eq.iboulist(inp1,i)) then
               ibou_flag=+1
               goto 9999
            endif
         enddo
      else
         if(irray(inp1,ifc).eq.0) ibou_flag=+1
      endif

 9999 continue

      return

      end

c...................................................................
