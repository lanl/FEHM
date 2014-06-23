      subroutine inflogh(iflg)
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
!D1 Read in sources and sinks for seepage faces and drainage.
!D1 Read in sources and sinks for different inputs
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: 2005, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inflogh.f_a  $
!D2 
!***********************************************************************
!D3
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************

      use comai
      use combi
      use comdi
      use comdti
      use comki
      implicit none

      integer i,j,icode, max_gh,iflg
      real*8 rhodvis, aiptmp, wellim_tol
      real*8, allocatable :: aiped(:)
      real*8, allocatable::  sktmp(:)
      real*8, allocatable ::  esktmp(:)
      real*8, allocatable ::  distmp(:)
      integer, allocatable ::  katmp(:)
      integer, allocatable ::  idir_tmp(:)
c     density divided by viscosity at 0.1 Mpa and 20 C
      parameter (rhodvis = 995360.835533628,max_gh = 100000)
      parameter (wellim_tol= 1.e-35)

      macro = 'flgh'

      if(iflg.eq.0) then
         
         if(.not.allocated(pflow_gh)) then
            ngh = 0
            allocate (wellim_gh(max_gh), pflow_gh(max_gh),
     &           node_gh(max_gh))
            allocate (idir_gh(max_gh))
            wellim_gh = 0.0
            node_gh = 0
            pflow_gh = 0.0
            idir_gh = 0
         endif
         allocate(aiped(n0),esktmp(n0),sktmp(n0),distmp(n0),katmp(n0))
         allocate(idir_tmp(n0))      
c**** read flow data ****
         narrays = 5
         itype(1) = 8
         itype(2) = 8
         itype(3) = 8
         itype(4) = 8
         itype(5) = 4
         default(1) = 0.
         default(2) = 0.
         default(3) = 0.
         default(4) = 0.
         default(5) = 0
         igroup = 1

c     generalized head condition
c     only applicable to saturated only conditions
c     esk carries the permebility
c     wellim(i) [aiped] has the area divided by distance
c     distmp has the  distance
c     this formulation allows accumulation
         
         call initdata2 (inpt, ischk, n0, narrays, itype, 
     *        default, macroread(12), macro, igroup, ireturn,
     *        r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0),
     *        r8_4=distmp(1:n0),i4_1=idir_tmp(1:n0)) 
         
         do i = 1, n0
            if (sktmp(i) .ne. default(1) .or. esktmp(i) .ne. default(2)
     *           .or. aiped(i) .ne. default(3))then
               aiptmp = aiped(i)/distmp(i)*esktmp(i)*rhodvis*1.d06
               ngh = ngh + 1
               if(ngh.gt.max_gh) then
                  if (iptty .ne. 0) write(iptty,110) max_gh
                  write(ierr,110) max_gh   
               endif 
               node_gh(ngh) = i
               idir_gh(ngh) = idir_tmp(i)
               pflow_gh(ngh) = sktmp(i)
               wellim_gh(ngh) = aiptmp
            end if

         end do      
         deallocate(aiped,esktmp,sktmp,katmp,distmp,idir_tmp)

      else  if(iflg.eq.1) then

         allocate(sktmp(n0))
         allocate(aiped(n0))
         sktmp = 0.0
         aiped = 0.0
         do j = 1,ngh
            i = node_gh(j)
            esk(i)=1.0
            ka(i) = -1	  
            sktmp(i) = pflow_gh(j)*wellim_gh(j)+sktmp(i)
            aiped(i) = wellim_gh(j)+aiped(i)
         enddo
         do j = 1,ngh
            i = node_gh(j)
            pflow(i) = sktmp(i)/aiped(i)
            call headctr(5,i,pflow(i),pflow(i))
            wellim(i) = aiped(i)
         enddo

         deallocate(sktmp,aiped)

      endif
 110  format(1x,' Error: Number of Generalized head BCs', 
     &     ' exceed maximum allowed(', i8,')')
      
      end
