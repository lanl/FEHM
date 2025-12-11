      subroutine move_particle(ipart, spacing_subst, 
     2     local_pos, local_spacing,
     3     too_many_jumps, coord_flag)
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
!D1 Perform random walk shift of particle for dispersion between cells 
!D1 in the streamline particle tracking algorithm.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.10 [10086-STN-2.10-00]
!D2 
!D2 Initial implementation: ?, Programmer: ?
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/move_particle.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:10:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:04   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:31:08   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:14   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
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

      use comsptr
      use comdti
      use comai
      use combi
      use comdi
      implicit none

      integer maxjumps
      parameter(maxjumps = 10)
      integer maxjumps1
      parameter(maxjumps1 = maxjumps +1)
      integer nodewrite(maxjumps1)
      integer ipart, current_node, current_model, transflag
      integer ijump, coord_flag
      logical done
      real*8 spacing_subst(0:*)
      real*8 nondim, local_pos, local_spacing
      logical too_many_jumps

      too_many_jumps = .false.
      current_node = ijkv(ipart)
      if(current_node.ne.0) then
         current_model = itrc(ijkv(ipart))
         transflag = tprpflag(current_model)
      else
         transflag = 0
      end if
      if(transflag.eq.2.or.transflag.eq.4.or.transflag.eq.13.or.
     &     transflag.eq.14) then
         if(coord_flag.eq.1) ijkvs(ipart) = ijkv(ipart)
         ic_new(ipart) = 0
c     Only move the particle if it hasn't already left the system
c     BAR - 12-4-98
            
         if(ijkv(ipart).ne.0) then
c            if(irray(ijkv(ipart),0).gt.0) then
               done = .false.
               ijump = 0
               nodewrite = 0
c     do while loop handles the case of particles that move
c     more than one grid cell in the x direction
               do while(.not.done)
                  ijump = ijump + 1
                  nodewrite(ijump) = ijkv(ipart)
                  if(ijump.gt.maxjumps) then
                     call write_error(coord_flag)
                     too_many_jumps = .true.
                     done = .true.
                     goto 1000
c                     stop
                  end if
                  nondim = local_pos/local_spacing
                  if(nondim.ge.1.) then
                     if(irray(ijkv(ipart), coord_flag).ne.0) then
c     Move particle forward to adjacent cell in x direction
                        ic_new(ipart)=ic_new(ipart)+1
                        ijkv(ipart)=irray(ijkv(ipart), coord_flag)
                        local_pos=local_pos-local_spacing
                        local_spacing = spacing_subst(ijkv(ipart))
                     else
c     Reflect particle back in x direction
                        ic_new(ipart)=ic_new(ipart)+1
                        local_pos=2*local_spacing-local_pos
                     end if
                  elseif(nondim.le.0.) then
                     if(irray(ijkv(ipart),-coord_flag).ne.0) then
c     Move particle backward to adjacent cell in x direction
                        ic_new(ipart)=ic_new(ipart)+1
                        ijkv(ipart)=irray(ijkv(ipart),-coord_flag)
                        local_pos=local_pos+
     2                       spacing_subst(ijkv(ipart))
                        local_spacing = spacing_subst(ijkv(ipart))
                     else
c     Reflect particle forward in x direction
                        ic_new(ipart)=ic_new(ipart)+1
                        local_pos=-local_pos
                     end if
                  else
                     done = .true.
                  end if
c                  if(irray(ijkv(ipart),0).lt.0) done=.true.
                  
c     Internal subroutine to determine if particle has previously
c     or just now reached a btc zone
                  
                  call seek_btczone(ipart, ijkv(ipart))
 1000             continue
                  
               end do
               
c     Node may have changed, need to change ddyv and ddzv
               if(coord_flag.eq.1) then
                  ddyv(ipart) = ddy(ijkv(ipart))
                  ddzv(ipart) = ddz(ijkv(ipart))
               elseif(coord_flag.eq.2) then
                  ddxv(ipart) = ddx(ijkv(ipart))
                  ddzv(ipart) = ddz(ijkv(ipart))
               elseif(coord_flag.eq.3) then
                  ddxv(ipart) = ddx(ijkv(ipart))
                  ddyv(ipart) = ddy(ijkv(ipart))
               end if
               
c            end if  
         end if
      end if




      
      return

      contains

      subroutine write_error(idirection)
      implicit none
      integer idirection

      if(iptty.gt.0) then
c         write(iptty,*) 'Fatal Error in ptrac3'
c         write(iptty,*) 'Particle has attempted to jump'
c         write(iptty,*) 'more than ',maxjumps,'cells'
c         write(iptty,*) 'in this random walk time step.'
c         write(iptty,*) 'Reduce particle tracking time step.'
c         write(iptty,*) 'Particle number', ipart
         if(idirection.eq.1) then
c            write(iptty,*) 
c     2'Warning - too many jumps in x for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iptty,*) 'Fatal direction for jump: x'
         elseif(idirection.eq.2) then
c            write(iptty,*) 
c     2'Warning - too many jumps in y for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iptty,*) 'Fatal direction for jump: y'
         else
c            write(iptty,*) 
c     2'Warning - too many jumps in z for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iptty,*) 'Fatal direction for jump: z'
         end if
c         write(iptty,*) 'Final cell number:',ijkv(ipart)
c         write(iptty,*) 'Final relative x value:',x2(ipart)
c         write(iptty,*) 'Final relative y value:',y2(ipart)
c         write(iptty,*) 'Final relative z value:',z2(ipart)
c         write(iptty,*) 'Node number sequence for jumps:'
c         write(iptty,2100) nodewrite
      end if
      if(iout.gt.0) then
c         write(iout,*) 'Fatal Error in ptrac3'
c         write(iout,*) 'Particle has attempted to jump'
c         write(iout,*) 'more than ',maxjumps,'cells'
c         write(iout,*) 'in this random walk time step.'
c         write(iout,*) 'Reduce particle tracking time step.'
c         write(iout,*) 'Particle number', ipart
         if(idirection.eq.1) then
c            write(iout,*) 
c     2'Warning - too many jumps in x for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iout,*) 'Fatal direction for jump: x'
         elseif(idirection.eq.2) then
c            write(iout,*) 
c     2'Warning - too many jumps in y for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iout,*) 'Fatal direction for jump: y'
         else
c            write(iout,*) 
c     2'Warning - too many jumps in z for particle', ipart,
c     3'cell number',ijkv(ipart)
c            write(iout,*) 'Fatal direction for jump: z'
         end if
c         write(iout,*) 'Final cell number:',ijkv(ipart)
c         write(iout,*) 'Final relative x value:',x2(ipart)
c         write(iout,*) 'Final relative y value:',y2(ipart)
c         write(iout,*) 'Final relative z value:',z2(ipart)
c         write(iout,*) 'Node number sequence for jumps:'
c         write(iout,2100) nodewrite
      end if
 2100 format(15i8)
      return
      end subroutine write_error



      end subroutine move_particle

