      subroutine evaporation(iflag)
!***********************************************************************
! Copyright 2008 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Output concentration flux for a zone.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.30
!D2 Initial implementation: 14-Sep-07,  Programmer: Phil Stauffer
!D2 Made into subrotuine  : 03-Sep-08, Programmer: Z. Dash 
!D2
!***********************************************************************

      use comai , only : days, inpt, neq
      use combi , only : sx1
      use comdi, only : s, sk
      use comevap

      implicit none
      integer i, iflag, open_file, ievap, pflag, pnode
      real*8  s_dum, fyr
      character*1000 evap_file

      save ievap

      if (iflag .eq. 0) then
c Read evaporation function parameters and evaporation node list
         read (inpt, '(a1000)') evap_file
         ievap = open_file(evap_file, 'old')

c Read past first header line: Truncation Saturation for function
         read(ievap,*)
         read(ievap,*) evap_trunc
c Read past second header line: Four sets of A + B*sat + Csat^2  time1 time2
         read(ievap,*)
c Read evaporation equation parameters
         read(ievap,*) par(1,1), par(1,2), par(1,3), tyr(1), tyr(2)
         read(ievap,*) par(2,1), par(2,2), par(2,3), tyr(3), tyr(4)
         read(ievap,*) par(3,1), par(3,2), par(3,3), tyr(5), tyr(6)
         read(ievap,*) par(4,1), par(4,2), par(4,3), tyr(7), tyr(8)
c Read past third header line:  Number of Evaporation Nodes
         read(ievap,*)
         read(ievap,*) num_evap
         allocate(evap_node(num_evap),area_node(num_evap))
         allocate(evap_flag(neq))
         evap_node = 0
         area_node = 0.
         evap_flag = .false.
c Read past fourth header line: Evaporation nodes, Area factor
         read(ievap,*) 
         read(ievap,*) (evap_node(i), area_node(i), i = 1,num_evap)
         close(ievap)
         pflag = 1
         do i = 1, neq
            if(i .eq. evap_node(pflag)) then
               evap_flag(evap_node(pflag)) = .true.
               pflag = pflag + 1
               if(pflag .gt. num_evap) exit
            end if
         end do
         pflag = 1
      else if (iflag .eq. 1) then
         ievap = open_file('radial_nodal_volumes.txt','unknown')
         do i = 1, neq
            write(ievap,*) i, sx1(i)
         end do
         close (ievap)
         ievap = open_file('time_vs_evap.txt','unknown')
      else if (iflag .eq. 2) then
        fyr = dmod(days,365.d0)/365.d0

        do i = 1, num_evap
	  pnode = evap_node(i)
          s_dum = s(pnode)
          if(s_dum .gt. 0.00001) then
            if(s_dum .gt. evap_trunc) s_dum =  evap_trunc               

            if ((fyr .ge. tyr(1)) .and. (fyr .le. tyr(2))) then
              sk(pnode) = par(1,1) + par(1,2)*s_dum 
     x                             + par(1,3)*s_dum*s_dum  

            else if ((fyr.ge.tyr(3)) .and. (fyr .le. tyr(4))) then
              sk(pnode) = par(2,1) + par(2,2)*s_dum
     x                             + par(2,3)*s_dum*s_dum

            else if ((fyr .ge. tyr(5)).and. (fyr.le. tyr(6))) then
              sk(pnode) = par(3,1) + par(3,2)*s_dum
     x                             + par(3,3)*s_dum*s_dum

            else if ((fyr .ge. tyr(7)) .and. (fyr.le.tyr(8))) then
              sk(pnode) = par(4,1) + par(4,2)*s_dum
     x                             + par(4,3)*s_dum*s_dum

            end if

            sk(pnode) = sk(pnode)*area_node(i)                 
         else
	    sk(pnode) = 1.e-19
         end if
      end do

      write(ievap, 666) days,sk(evap_node(1)),s(evap_node(1)),
     2     sk(evap_node(num_evap)),s(evap_node(num_evap))

      end if
         
 666  format(f12.2,4e15.3)

      end subroutine evaporation
