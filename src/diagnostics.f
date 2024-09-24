      subroutine diagnostics(iflg)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  Check for input consistency
CD1  To print out the largest residuals
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/diagnostics.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:50   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:02:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:08:18   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:10   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:59:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:40:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.13   Mon Mar 31 08:33:30 1997   gaz
CD2 minor changes
CD2 
CD2    Rev 1.12   Tue May 14 14:32:16 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.11   Fri Mar 01 14:37:18 1996   gaz
CD2 changes to improve output
CD2 
CD2    Rev 1.10   Wed Feb 07 10:54:54 1996   gaz
CD2 cleaned up output
CD2 
CD2    Rev 1.9   Wed Jan 17 09:58:56 1996   zvd
CD2 Minor updates to prolog
CD2 
CD2    Rev 1.8   Wed Jan 10 11:01:32 1996   hend
CD2 Corrected ECD Number
CD2 
CD2    Rev 1.7   Tue Jan 09 15:51:14 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.6   12/19/95 14:20:40   gaz
CD2 minor fix for 1dof problems
CD2 
CD2    Rev 1.5   11/15/95 10:24:34   gaz
CD2 corrected minor bug
CD2 
CD2    Rev 1.4   09/20/95 16:13:14   gaz
CD2 more information outputted
CD2 
CD2    Rev 1.3   09/14/95 16:14:40   gaz
CD2 corrected error to allow 4 dof output
CD2 
CD2    Rev 1.2   09/11/95 16:58:30   gaz
CD2 added 4dof to routine, made correction to 3dof
CD2 
CD2    Rev 1.1   08/09/95 15:53:32   zvd
CD2 Corrected write to iptty,900 when unassigned.
CD2 
CD2    Rev 1.0   06/02/95 08:51:30   llt
CD2 original version
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use davidi
      use comgi
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comco2
      implicit none

      integer ibp,ibp1,ibp2,ibp3,ibp4,ibp5,ibp6
      integer iflg,i,idiagnos
      real*8  s_upper,s_lower
      real*8  bpmax,bp1max,bp2max,bp3max,bp4max,bp5max,bp6max
      real*8  bpmin_min, bpmax_max
      parameter(s_upper=1.0000,s_lower=0.0000)
      parameter(bpmin_min=-1.e10, bpmax_max=1.e10)
      save ibp,ibp1,ibp2,ibp3,ibp4,ibp5,ibp6
      save bpmax,bp1max,bp2max,bp3max,bp4max,bp5max,bp6max
      if(.not.allocated(bp)) return
      if(iflg.eq.0) then
         if(ihead.ne.0) then
            if(igrav.eq.0) then
               if (iout .ne. 0) write(iout, 100) 
               if (ischk .ne. 0) write(ischk, 100) 
            endif
         endif
 100     format ('>>>> gravity not set for head problem: stopping <<<<')
         if(ivf.eq.-1) then
            if(icnl.ne.0) then
               if (iout .ne. 0) write(iout, 200) 
               if (ischk .ne. 0) write(ischk, 200) 
            endif
 200        format ('>>>> dimension (icnl) not set to 3 for FDM: ',
     &           'stopping <<<<')
         endif
      else if(iflg.eq.-1) then
c
c     find largerst residuals  
c
         if(idof.eq.1 .or. idof_co2.eq.1) then
            ibp=0
            ibp1=0
            bpmax=-1.0
            bp1max=-1.0
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.bp1max) then
                  bp1max=abs(bp(i+nrhs(1)))
                  ibp1=i
               endif
               if(abs(bp(i+nrhs(1))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(1)))
                  ibp=i
               endif
            enddo

         else if(idof.eq.2 .or. idof_co2.eq.2) then
            ibp=0
            ibp1=1
            ibp2=1
            bpmax=-1.0
            bp1max=-1.0
            bp2max=-1.0
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.bp1max) then
                  bp1max=abs(bp(i+nrhs(1)))
                  ibp1=i
               endif
               if(abs(bp(i+nrhs(1))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(1)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bp2max) then
                  bp2max=abs(bp(i+nrhs(2)))
                  ibp2=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(2)))
                  ibp=i
               endif
            enddo

         else if(idof.eq.3 .or. idof_co2.eq.3) then
            ibp=0
            ibp1=0
            ibp2=0
            ibp3=0
            bpmax=-1.0
            bp1max=-1.0
            bp2max=-1.0
            bp3max=-1.0
            do i=1,neq
               if(abs(bp(i+nrhs(1))).gt.bp1max) then
                  bp1max=abs(bp(i+nrhs(1)))
                  ibp1=i
               endif
               if(abs(bp(i+nrhs(1))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(1)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bp2max) then
                  bp2max=abs(bp(i+nrhs(2)))
                  ibp2=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(2)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bp3max) then
                  bp3max=abs(bp(i+nrhs(3)))
                  ibp3=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(3)))
                  ibp=i
               endif
            enddo
            
            if (ibp.eq.0.or.ibp1.eq.0.or.ibp2.eq.0.or.
     &          ibp3.eq.0) then
             if (iptty .ne. 0) 
     &        write(iptty,*)"Nans for all residuals: stopping"
             if (ntty .ne. 0) 
     &        write(ntty,*)"Nans for all residuals: stopping"  
              write(ierr,*)"Nans for all residuals: stopping"   
              stop
            endif         


         else if(idof.eq.4 .or. idof_co2.eq.4) then
            ibp=0
            ibp1=0
            ibp2=0
            ibp3=0
            ibp4=0
            bpmax=-1.0
            bp1max=-1.0
            bp2max=-1.0
            bp3max=-1.0
            bp4max=-1.0
            do i=1,neq
               if(bp(i+nrhs(1)).lt.bpmin_min) bp(i+nrhs(1))=bpmin_min
               if(bp(i+nrhs(2)).lt.bpmin_min) bp(i+nrhs(2))=bpmin_min
               if(bp(i+nrhs(3)).lt.bpmin_min) bp(i+nrhs(3))=bpmin_min
               if(bp(i+nrhs(4)).lt.bpmin_min) bp(i+nrhs(4))=bpmin_min
c     
               if(bp(i+nrhs(1)).gt.bpmax_max) bp(i+nrhs(1))=bpmax_max
               if(bp(i+nrhs(2)).gt.bpmax_max) bp(i+nrhs(2))=bpmax_max
               if(bp(i+nrhs(3)).gt.bpmax_max) bp(i+nrhs(3))=bpmax_max
               if(bp(i+nrhs(4)).gt.bpmax_max) bp(i+nrhs(4))=bpmax_max
               if(abs(bp(i+nrhs(1))).gt.bp1max) then
                  bp1max=abs(bp(i+nrhs(1)))
                  ibp1=i
               endif
               if(abs(bp(i+nrhs(1))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(1)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bp2max) then
                  bp2max=abs(bp(i+nrhs(2)))
                  ibp2=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(2)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bp3max) then
                  bp3max=abs(bp(i+nrhs(3)))
                  ibp3=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(3)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(4))).gt.bp4max) then
                  bp4max=abs(bp(i+nrhs(4)))
                  ibp4=i
               endif
               if(abs(bp(i+nrhs(4))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(4)))
                  ibp=i
               endif
            enddo


         else if(idof.eq.6 .or. idof_co2.eq.6) then
            ibp=0 
            ibp1=0
            ibp2=0
            ibp3=0
            ibp4=0
            ibp5=0
            ibp6=0
            bpmax=-1.0
            bp1max=-1.0
            bp2max=-1.0
            bp3max=-1.0
            bp4max=-1.0
            bp5max=-1.0
            bp6max=-1.0
            do i=1,neq
               if(bp(i+nrhs(1)).lt.bpmin_min) bp(i+nrhs(1))=bpmin_min
               if(bp(i+nrhs(2)).lt.bpmin_min) bp(i+nrhs(2))=bpmin_min
               if(bp(i+nrhs(3)).lt.bpmin_min) bp(i+nrhs(3))=bpmin_min
               if(bp(i+nrhs(4)).lt.bpmin_min) bp(i+nrhs(4))=bpmin_min
               if(bp(i+nrhs(5)).lt.bpmin_min) bp(i+nrhs(5))=bpmin_min
               if(bp(i+nrhs(6)).lt.bpmin_min) bp(i+nrhs(6))=bpmin_min
c
               if(bp(i+nrhs(1)).gt.bpmax_max) bp(i+nrhs(1))=bpmax_max
               if(bp(i+nrhs(2)).gt.bpmax_max) bp(i+nrhs(2))=bpmax_max
               if(bp(i+nrhs(3)).gt.bpmax_max) bp(i+nrhs(3))=bpmax_max
               if(bp(i+nrhs(4)).gt.bpmax_max) bp(i+nrhs(4))=bpmax_max
               if(bp(i+nrhs(5)).gt.bpmax_max) bp(i+nrhs(5))=bpmax_max
               if(bp(i+nrhs(6)).gt.bpmax_max) bp(i+nrhs(6))=bpmax_max
               if(abs(bp(i+nrhs(1))).gt.bp1max) then
                  bp1max=abs(bp(i+nrhs(1)))
                  ibp1=i
               endif
               if(abs(bp(i+nrhs(1))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(1)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bp2max) then
                  bp2max=abs(bp(i+nrhs(2)))
                  ibp2=i
               endif
               if(abs(bp(i+nrhs(2))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(2)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bp3max) then
                  bp3max=abs(bp(i+nrhs(3)))
                  ibp3=i
               endif
               if(abs(bp(i+nrhs(3))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(3)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(4))).gt.bp4max) then
                  bp4max=abs(bp(i+nrhs(4)))
                  ibp4=i
               endif
               if(abs(bp(i+nrhs(4))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(4)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(5))).gt.bp5max) then
                  bp5max=abs(bp(i+nrhs(5)))
                  ibp5=i
               endif
               if(abs(bp(i+nrhs(5))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(5)))
                  ibp=i
               endif
               if(abs(bp(i+nrhs(6))).gt.bp6max) then
                  bp6max=abs(bp(i+nrhs(6)))
                  ibp6=i
               endif
               if(abs(bp(i+nrhs(6))).gt.bpmax) then
                  bpmax=abs(bp(i+nrhs(6)))
                  ibp=i
               endif
            enddo

          
         endif


      else if(iflg.eq.1) then
c     print out worst residuals  
         if(iad.le.abs(maxit)) then
            if (ntty .eq. 2) then
               if (iout .ne. 0) write(iout, 300) 
            endif
            if (iptty.ne.0 ) write(iptty, 300) 
 300        format(20x,'Largest Residuals')
         else
            if (iout .ne. 0) then 
               write(iout, *) '   '
               write(iout, 400) l
            end if
            if (iptty .ne. 0) then
               write(iptty, *) '   '
               write(iptty, 400) l
            endif
         endif
 400     format ('#### largest N-R corrections, timestep ',i8, ' ####')
         if(idof.eq.1 .or. idof_co2.eq.1) then
c gaz 080819  added zone ino to output 
            if (iptty .ne. 0) then
               write(iptty,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $              cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
            endif
            if (ntty .eq. 2) then
             if (iout .ne. 0) write(iout,900) 'EQ1',bp1max,izonef(ibp1),
     $              ibp1,cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
            endif
         else if(idof.eq.2 .or. idof_co2.eq.2) then

            if (iptty .ne. 0) then
               write(iptty,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $           cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
               write(iptty,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $           cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
            endif
            if (ntty .eq. 2) then
               if (iout .ne. 0) then
                  write(iout,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iout,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
               endif
            end if
         else if(idof.eq.3 .or. idof_co2.eq.3) then

            if (ibp.eq.0.or.ibp1.eq.0.or.ibp2.eq.0.or.
     &          ibp3.eq.0) then
             if (iptty .ne. 0) 
     &        write(iptty,*)"Nans for all residuals: stopping"
             if (ntty .ne. 0) 
     &        write(ntty,*)"Nans for all residuals: stopping"  
              write(ierr,*)"Nans for all residuals: stopping"   
              stop
            endif         
            if (iptty .ne. 0) then
                  write(iptty,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iptty,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iptty,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
            endif
            if (ntty .eq. 2) then
               if (iout .ne. 0) then
                  write(iout,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iout,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iout,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
               end if
            endif
         else if(idof.eq.4 .or. idof_co2.eq.4) then

            if (iptty .ne. 0) then
                  write(iptty,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iptty,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iptty,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
                  write(iptty,900) 'EQ4',bp4max,ibp4,izonef(ibp4),
     $            cord(ibp4,1),cord(ibp4,2),cord(ibp4,3)
            endif
            if (ntty .eq. 2) then
               if (iout .ne. 0) then
                  write(iout,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iout,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iout,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
                  write(iout,900) 'EQ4',bp4max,ibp4,izonef(ibp4),
     $            cord(ibp4,1),cord(ibp4,2),cord(ibp4,3)
               end if
            endif
         else if(idof.eq.6 .or. idof_co2.eq.6) then

            if (iptty .ne. 0) then
                  write(iptty,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iptty,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iptty,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
                  write(iptty,900) 'EQ4',bp4max,ibp4,izonef(ibp4),
     $            cord(ibp4,1),cord(ibp4,2),cord(ibp4,3)
                  write(iptty,900) 'EQ5',bp5max,ibp5,izonef(ibp5),
     $            cord(ibp5,1),cord(ibp5,2),cord(ibp5,3)
                  write(iptty,900) 'EQ6',bp6max,ibp6,izonef(ibp6),
     $            cord(ibp6,1),cord(ibp6,2),cord(ibp6,3)
            endif
            if (ntty .eq. 2) then
               if (iout .ne. 0) then
                  write(iout,900) 'EQ1',bp1max,ibp1,izonef(ibp1),
     $               cord(ibp1,1),cord(ibp1,2),cord(ibp1,3)
                  write(iout,900) 'EQ2',bp2max,ibp2,izonef(ibp2),
     $               cord(ibp2,1),cord(ibp2,2),cord(ibp2,3)
                  write(iout,900) 'EQ3',bp3max,ibp3,izonef(ibp3),
     $            cord(ibp3,1),cord(ibp3,2),cord(ibp3,3)
                  write(iout,900) 'EQ4',bp4max,ibp4,izonef(ibp4),
     $            cord(ibp4,1),cord(ibp4,2),cord(ibp4,3)
                  write(iout,900) 'EQ5',bp5max,ibp5,izonef(ibp5),
     $            cord(ibp5,1),cord(ibp5,2),cord(ibp5,3)
                  write(iout,900) 'EQ6',bp6max,ibp6,izonef(ibp6),
     $            cord(ibp6,1),cord(ibp6,2),cord(ibp6,3)
               endif
            end if
         endif
c 
c add largest residual to printout list
c
c     if(idualp.ne.0) then
c       nskw(m)=ibp
c       nskw(m+m)=ibp+neq
c       nskw(m+m+m)=ibp+neq+neq
c       if(m2.ne.0) then
c        nskw2(m2)=ibp
c        nskw2(m2+m2)=ibp+neq
c        nskw2(m2+m2+m2)=ibp+neq+neq
c       endif                
c     else if(idpdp.ne.0) then
c       nskw(m)=ibp
c       nskw(m+m)=ibp+neq
c       if(m2.ne.0) then
c        nskw2(m2)=ibp
c        nskw2(m2+m2)=ibp+neq
c       endif                
c     else
c       nskw(m)=ibp
c       if(m2.ne.0) then
c        nskw2(m2)=ibp
c       endif                
c     endif 
      else if(iflg.eq.2) then
         do i=1,n
            if (irdof .ne. 13 .or. ifree .ne. 0) then
c gaz 110123
              if(idoff.ne.-1) then
               if(s(i).gt.s_upper.and.strd.lt.1.0) then
                  s(i)=1.0
               else if(s(i).lt.s_lower.and.strd.lt.1.0) then
                  s(i)=0.0
               endif
              endif
            end if
         enddo
      endif
 900  format(1x,a3,1x,'R=',g12.4,1x,'node=',i7,1x,'zone=',i6,
     &     1x,'x=',g10.4,1x,'y=',g10.4,1x,'z=',g10.4)
      return
      end
