      subroutine set_ptrk(resflag, lstep, cur_part)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine sets up the starting nodes and times for particle 
CD1  tracking, activated with macro ptrk.  It is called only once, the 
CD1  first time the particle tracker is invoked.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2 $Log:   /pvcs.config/fehm90/src/set_ptrk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:30   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:24   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.8   Wed Mar 27 11:35:48 1996   hend
CD2 Updated for ptrk output options
CD2 
CD2    Rev 1.7   Tue Mar 05 16:20:20 1996   hend
CD2 Added initialization for p array
CD2 
CD2    Rev 1.6   Fri Jan 12 17:56:26 1996   llt
CD2 changed mmgetblk agruments
CD2 
CD2    Rev 1.5   Thu Jan 11 10:25:56 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.4   08/07/95 11:39:40   awolf
CD2 Fixed for DPDP (n0 instead of neq). Also changed for 
CD2 breakthrough output
CD2 
CD2    Rev 1.3   05/08/95 12:49:10   robinson
CD2 Calls plotc1 to write coordinate info to .trc file
CD2 
CD2    Rev 1.2   03/15/95 17:05:34   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.1   02/02/95 15:22:50   llt
CD2 added pvcs log info
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.5 Cell-based particle-tracking module
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

      use combi
      use comci
      use comdi
      use comflow
      use compart
      use comdti
      use comai
      use davidi, only : irdof
      
      implicit none
      integer i,i2,cur_part,positive,num,split,nsizep,icode,resflag
      integer lstep
      real*8 inc, inflow, denom, sdum

c Set the initial node position and starting times for every particle. The
c manner in which these are set is dependent on the type of problem being 
c run. Pcnsk must either be all positive or all negative. Positive 
c indicates injecting according to source flow, while negative indicates
c evenly distributing the particles over the nodes. Pcnsk gives the 
c relative factor to distribute -- i.e. a 2 indicates twice as many 
c particles in this node or set of nodes. 
      if( resflag .eq. 0 ) then
c make sure pcnsk is either always positive or always negative
         positive=3
         do i=1,n0
            if (pcnsk(i).gt.0.) then
               if (positive.eq.-1) goto 14
               positive=1
            else if (pcnsk(i).lt.0.) then
               if (positive.eq.1) goto 14
               positive=-1
            endif
         enddo
         goto 15
         
 14      write(ierr,*) 'ERROR: Pcnsk in ptrk must be either always '
         write(ierr,*) 'positive or always negative. '
         write(ierr,*) 'Code aborted in set_ptrk.f'
         if (iatty.ne.0) then
            write(iatty,*) 'ERROR: Pcnsk in ptrk must be either always '
            write(iatty,*) 'positive or always negative. '
            write(iatty,*) 'Code aborted in set_ptrk.f'
         endif
         stop
         
 15      continue
         
c     Take care of positive cnsk (apportion by source rate)
         
         if (positive.eq.1) then
            inflow=0
c     Liquid only and isothermal air-water
            if ((trak_type(1).eq.1).and.(ico2.lt.0)) then
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(sk(i).lt.0.)) then
                     inflow=inflow+abs(sk(i))*pcnsk(i)
                  endif
               enddo
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(sk(i).lt.0.)) then
                     num=int(num_particles(1)*abs(sk(i))*
     2                    pcnsk(i)/inflow)
                     if(num.gt.0) then
                        inc=(t2sk(i)-t1sk(i))/num
                     end if
                     do i2=1,num
                        box(cur_part,1)=i
                        start_time(cur_part,1)=86400.*(t1sk(i)+
     2                       inc*(i2-1))
                        cur_part=cur_part+1
                     enddo
                  endif
               enddo
c     else if vapor only and isothermal air-water
            else if ((trak_type(1).eq.2).and.(ico2.lt.0)) then
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(qh(i).lt.0.)) then
                     inflow=inflow+abs(qh(i))*pcnsk(i)
                  endif
               enddo
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(qh(i).lt.0.)) then
                     num=int(num_particles(1)*abs(qh(i))*
     2                    pcnsk(i)/inflow)
                     if(num.gt.0) then
                        inc=(t2sk(i)-t1sk(i))/num
                     end if
                     do i2=1,num
                        box(cur_part,1)=i
                        start_time(cur_part,1)=86400.*
     2                       (t1sk(i)+inc*(i2-1))
                        cur_part=cur_part+1
                     enddo
                  endif
               enddo
c     else if liquid only and not isothermal air-water
            else if ((trak_type(1).eq.1).and.(ico2.ge.0)) then
               do i=1,n0
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     sdum = s(i)
                  else
                     sdum = 1.
                  end if
                  if ((pcnsk(i).gt.0.).and.(sk(i)*sdum.lt.0.)) then
                     inflow=inflow+abs(sk(i))*pcnsk(i)
                  endif
               enddo
               do i=1,n0
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     sdum = s(i)
                  else
                     sdum = 1.
                  end if
                  if ((pcnsk(i).gt.0.).and.(sk(i)*sdum.lt.0.)) then
                     num=int(num_particles(1)*abs(sk(i))*
     2                    sdum*pcnsk(i)/inflow)
                     if(num.gt.0) then
                        inc=(t2sk(i)-t1sk(i))/num
                     end if
                     do i2=1,num
                        box(cur_part,1)=i
                        start_time(cur_part,1)=86400.*
     2                       (t1sk(i)+inc*(i2-1))
                        cur_part=cur_part+1
                     enddo
                  endif
               enddo
c     else if vapor only and not isothermal air-water
            else if ((trak_type(1).eq.2).and.(ico2.ge.0)) then
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(sk(i)*(1-s(i)).lt.0.)) then
                     inflow=inflow+abs(sk(i))*pcnsk(i)
                  endif
               enddo
               do i=1,n0
                  if ((pcnsk(i).gt.0.).and.(sk(i)*(1-s(i)).lt.0.)) then
                     num=int(num_particles(1)*abs(sk(i))*
     2                    (1-s(i))*pcnsk(i)/inflow)
                     if(num.gt.0) then
                        inc=(t2sk(i)-t1sk(i))/num
                     end if
                     do i2=1,num
                        box(cur_part,1)=i
                        start_time(cur_part,1)=86400.*(t1sk(i)+
     2                       inc*(i2-1))
                        cur_part=cur_part+1
                     enddo
                  endif
               enddo
               
            endif    

c     Take care of negative pcnsk
         else
            
            split=0
            do i=1,n0
               split=split-pcnsk(i)
            enddo
            do i=1,n0
               num=int(num_particles(1)*(-pcnsk(i))/split)
               if(num.gt.0) then
                 inc=(t2sk(i)-t1sk(i))/num
                 do i2=1,num
                   box(cur_part,1)=i
                   start_time(cur_part,1)=86400.*(t1sk(i)+inc*(i2-1))
                   cur_part=cur_part+1
                 enddo
	       endif
            end do
         endif
         
         num_particles(1)=cur_part-1
c     nsport needs to be set here because the number of particles
c     is not computed until this point
         nsport(2) = num_particles(1)

c     Option 5 has particle concentrations

         if(pout.eq.5) then
            do i = 1, num_particles(1)
               partconc(i) = nodepconc(box(i,1))
            end do
         end if


         do i2 = 1, num_particles(1)
            frac_done(i2,1) = 0.
         end do
      end if
      if (idpdp.eq.0) then
         nsizep = nelm(neq+1) - neq - 1
      else
         nsizep = (nelm(neq+1)-neq-1)*2
      endif
c      call mmgetblk ("p","compart",ipp,
c     +     nsizep,2,icode)
	if(.not. allocated(p)) then
           allocate(p(nsizep))
        end if
        p = -5

      if (pout.eq. -7) then
         write(istrc,100)
      else if (pout.lt.0) then
      write(istrc,*) 'Node Number, Time A Particle Left This Node (sec)'
      else
         call plotc1(0,0)
      endif
 100  format (2x, 'Particle #', 4x, 'Node Number', 5x, 'Zone', 2x, 
     &     'Time in node', 3x, 'Current time (sec)')

      if (pout.eq.0) then
         if (trak_type(1).eq.1) then
            do i2=1,num_particles(1)
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  denom = (rolf(box(i2,1))*ps(box(i2,1))
     2                 *s(box(i2,1))*sx1(box(i2,1)))
                  if(denom .ne. 0.) then
                     an(box(i2,1))=an(box(i2,1))+1/denom
                  else
                     an(box(i2,1)) = 0
                  end if
                  anv(box(i2,1))=an(box(i2,1))
               else
                  sdum = 1.
                  denom = (rolf(box(i2,1))*ps(box(i2,1))
     2                 *sdum*sx1(box(i2,1)))
                  if(denom .ne. 0.) then
                     an(box(i2,1))=an(box(i2,1))+1/denom
                  else
                     an(box(i2,1)) = 0
                  end if
               end if
            enddo
         else
            do i2=1,num_particles(1)
               denom = (rovf(box(i2,1))*ps(box(i2,1))
     2              *(1-s(box(i2,1)))*sx1(box(i2,1)))
               if(denom .ne. 0.) then
                  an(box(i2,1))=an(box(i2,1))+1/denom
               else
                  an(box(i2,1)) = 0.
               end if
               anv(box(i2,1))=an(box(i2,1))
            enddo
         endif
      endif
      
      return 
      end





