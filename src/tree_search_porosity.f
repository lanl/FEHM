      subroutine tree_search_porosity(inp1,level_max,flag_box)
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
c s kelkar oct 20 08
cmodified from tree_search, to check for positive porosity neighbour
      
      use comsk, only : omr_flag
      use comdi, only : ps

      implicit none

      integer inp1,i,n_oldlist,n_newlist,n1,index
      integer n_level,level_max, maxdim,flag_box,n_total
      integer, allocatable :: oldlist(:)
      integer, allocatable :: newlist(:)

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
               if(ps(i).gt.0.) then
                  flag_box=i
                  goto 99999
               endif
            enddo
            
c     augement oldlist with newlist. Keeping the previous oldlist
c     to go 2 levels back in compress_list
            n_total=n_oldlist+n_newlist
            if(n_total.gt.maxdim) goto 99999
            do n1=1,n_newlist
               oldlist(n1+n_oldlist)=newlist(n1)
            enddo
            
         enddo
c positive porosity not found. Set flag_box and return
         flag_box=0
         
99999    deallocate(oldlist)
         deallocate(newlist)

      return

      end

c..................................................................
