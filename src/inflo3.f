      subroutine inflo3
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Read in sources and sinks for seepage faces and drainage.
!D1 
!***********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: ?, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/inflo3.f_a  $
!D2
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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

      integer i,icode
      real*8, allocatable :: aiped(:)
      real*8, allocatable::  sktmp(:)
      real*8, allocatable ::  esktmp(:)
      integer, allocatable ::  katmp(:)

      macro = 'flo3'
      
      allocate(aiped(n0),esktmp(n0),sktmp(n0),katmp(n0))
	      
c**** read flow data ****
      narrays = 4
      itype(1) = 8
      itype(2) = 8
      itype(3) = 8
      itype(4) = 4
      default(1) = 0.
      default(2) = 0.
      default(3) = 0.
      default(4) = 0
      igroup = 1
      
      call initdata2 (inpt, ischk, n0, narrays, itype, 
     *     default, macroread(12), macro, igroup, ireturn,
     *     r8_1=sktmp(1:n0),r8_2=esktmp(1:n0),r8_3=aiped(1:n0),
     *     i4_1=katmp(1:n0)) 
      
      do i = 1, n0
         if (sktmp(i) .ne. default(1) .or. esktmp(i) .ne. default(2)
     *        .or. aiped(i) .ne. default(3).or.
     *         katmp(i) .ne. default(4)) then
            esk(i) = esktmp(i)
            if(esk(i).gt.0.and.igrav.ne.0.and.ico2.ge.0) then
             esk(i) = esktmp(i)-grav*cord(i,igrav)
            endif            
            if (katmp(i).eq.1) then
               sk(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = 0.0d00
               ka(i) = 1
            else if (katmp(i).eq.-1) then
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = aiped(i) * 1.0e+06
               ka(i) = -1             
            else if (katmp(i).eq.-2) then
c outflow only                        
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = aiped(i) * 1.0e+06
               ka(i) = -2               
            else if (katmp(i).eq.-23.or.katmp(i).eq.-24) then
c free drainage                
               esk(i) = esktmp(i)
               wellim(i) = sktmp(i)
	         pflow(i) = 0.0
               ka(i) = katmp(i)
            else if (katmp(i).eq.-3) then
c seepage face  - 2 phase and wtsi  
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -3
            else if (katmp(i).eq.-4) then
c seepage face  average pressure
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -4
            else if (katmp(i).eq.-5) then
c seepage face  (only flow if cell is full) wtsi only
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -5
            else if (katmp(i).eq.-6) then
c set ponding condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -6
            else if (katmp(i).eq.-7) then 
c set air only pressure with no water source
c gaz 081323 allow for outflow only 
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
c               wellim(i) = abs(aiped(i)) * 1.0e+06
               wellim(i) = aiped(i)*1.0e+06
               ka(i) = -7
            else if (katmp(i).eq.-9) then 
c set air  pressure, set water saturation
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -9               
            else if (katmp(i).eq.-10) then 
c set air  pressure, specified flow
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -10
            else if (katmp(i).eq.-11) then 
c set air  pressure, set water saturation
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
               wellim(i) = aiped(i)
               ka(i) = -11
             else if (katmp(i).eq.-12) then 
c set air  pressure, set water saturation
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
               wellim(i) = aiped(i)
               ka(i) = -12
            else if (katmp(i).eq.-13.or.katmp(i).eq.-14.or.
     &               katmp(i).eq.-15.or.katmp(i).eq.-17) then 
c set air  pressure, set water saturation
               pflow(i) = sktmp(i) 
               esk(i) = esktmp(i)
               wellim(i) = aiped(i)
               ka(i) = katmp(i)   
             else if (katmp(i).eq.-16) then 
c set air  flow, set water flow
               sk(i) = sktmp(i) 
               esk(i) = esktmp(i)
               ka(i) = -16                                        
            else if (katmp(i).eq.-8) then
c set ponding condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -8
            else if (katmp(i).eq.-22) then
c set special outflow only condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -22
            else if (katmp(i).eq.-21) then
c gaz 111519 
c set special air only flow when sat < 1. condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -21   
            else if (katmp(i).eq.-20.or.katmp(i).eq.-25) then 
c gaz 121219 
c set special air only flow when sat < 1. condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = katmp(i)    
            else if (katmp(i).eq.-18) then 
c gaz 082023 outflows for air and water
c linear relperms
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = katmp(i)                           
            else if (katmp(i).eq.-19) then
c gaz 121219 
c set special air only flow when sat < 1. condition
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -19               
            else if (katmp(i).eq.-101) then
c gaz 051120 this BC seekes to adjust bounday pressure to maintain a positive (or small) flow
c different than other models as it adjust boundary pressure - not not gridblock pressure
c pflow() is the initial boundary pressure
c esk() is a tolerance
c wellim() is the usual impedance              
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -101   
            else if (katmp(i).eq.-201) then
c gaz 051720 this BC seekes to adjust bounday pressure based on changing densities (cden related)
c will explicitly change boundary pressures after a transport simulation is finished,
c pflow() is the initial boundary pressure
c esk() is a tolerance
c wellim() is the usual impedance              
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -201 
           else if (katmp(i).eq.-202) then     
c gaz 112923 test for dissolved species   
c specified pressure and saturation   
               pflow(i) = sktmp(i)
               esk(i) = esktmp(i)
               wellim(i) = abs(aiped(i)) * 1.0e+06
               ka(i) = -202  
           end if
       endif
      end do
      
      deallocate(aiped,esktmp,sktmp,katmp)                                 
           
      end
