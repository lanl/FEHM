      subroutine cascade_sat(iflg)
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
!D1 To distribute saturations over 1.0.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Dec-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/cascade_sat.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:42:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2
!D2    Rev 2.4   29 Jan 2003 08:55:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 ?
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
c
c subroutine to distribute saturations over 1.0
c
c should work for steady state,may work for transient
c
      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
c
      implicit none
c  
      integer iflg,i,j,i1,i2,ii,kb,addnode,inneq
      integer indexa_axy,add_node0,idummy,addnode0
      real*8 sat_tol,sat_dif,fluid_dif,sumfout,sumfrac
      logical matrix_node
      parameter(sat_tol=1.d-6)
      addnode0=nelm(neq+1)-neq-1
c
      if(iflg.eq.-1) then
c allocate memory
       if(.not.allocated(flux_ts)) allocate(flux_ts(n0))
      else if(iflg.eq.0) then
c zero out timestep flux term
       flux_ts=0.0d00
      elseif(iflg.eq.1) then
       if(iad.gt.2) return
       do i=1,n0
        if(s(i).gt.1.0+sat_tol) then
c
         sat_dif=s(i)-1.0d00
         fluid_dif=sat_dif*rolf(i)*ps(i)*sx1(i)
         s(i)=1.0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
         if(i.gt.neq) then
            inneq = i-neq
            matrix_node = .true.
            addnode = addnode0                   
            idummy = neq
         else
            matrix_node = .false.
            inneq = i
            addnode = 0
            idummy = 0
         end if
         i1=nelm(i)+1
         i2=nelm(i+1)
         do ii=i1,i2
          indexa_axy = ii-neq-1+addnode
c     add to sum if it is flow out of the node
          kb=nelm(ii)
          if(a_axy(indexa_axy).gt.0.0) then
           sumfout = sumfout + a_axy(indexa_axy)
          endif
         enddo
         do ii=i1,i2
          indexa_axy = ii-neq-1+addnode
          kb=nelm(ii)
          if(a_axy(indexa_axy).gt.0.0) then
           sumfrac = a_axy(indexa_axy)/sumfout
           flux_ts(kb)=-sumfrac*fluid_dif/dtot
          endif
         enddo
        endif
       enddo
      elseif(iflg.eq.2) then
c load flux corrections into rhs
       do i=1,n
        bp(i+nrhs(1))= bp(i+nrhs(1))+flux_ts(kb)
       enddo
      endif
      return
      end
