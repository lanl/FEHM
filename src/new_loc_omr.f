      subroutine new_loc_omr(itime,istep,np1,edt,x_ts,y_ts,z_ts,
     $     time_step_flag,vx_zheng,vy_zheng,vz_zheng,rw,inp_adv)
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
!D2 $Log:   /pvcs.config/fehm90/src/new_loc_omr.f_a  $
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

      use comai, only : ierr, iptty
      use comdi, only : itrc
      use combi, only : cord
      use comsk
      use comsptr

      implicit none

      logical done,inside

      integer inp1,icapture,ijkv_in,walk_counter,upcoming_node
      integer ijkvsave,ijkvssave,current_node,ijkv_last
      integer current_model,np1,time_step_flag,itime,istep
      integer max_level,flag_box,inp_adv
      integer iii,ipc
      integer save_out_face(6), nout_face

      real*8 dx,dy,dz,x,xx,y,yy,z,zz,xpg,ypg,zpg
      real*8 vx_zheng,vy_zheng,vz_zheng,rw,v_x,v_y,v_z,vcurrent
      real*8 x2save,y2save,z2save,ddxvsave,ddyvsave,ddzvsave
      real*8 x_ts,y_ts,z_ts,edt,vnew,vdivd
      real*8 xc,yc,zc,epsilonv

      epsilonv=+1.e-24
      
      x2save = x2(np1)
      y2save = y2(np1)
      z2save = z2(np1)
      ijkvsave = ijkv(np1)
      ijkv_in=ijkv(np1)
      ijkvssave = ijkvs(np1)
      ddxvsave = ddxv(np1)
      ddyvsave = ddyv(np1)
      ddzvsave = ddzv(np1)

      if(ijkvsave.gt.0) then
         
         done = .false.
c     Cell-based velocities
c     v_x = 0.5*(ggg(ijkvsave,+1)-ggg(ijkvsave,-1))
c     v_y = 0.5*(ggg(ijkvsave,+2)-ggg(ijkvsave,-2))
c     v_z = 0.5*(ggg(ijkvsave,+3)-ggg(ijkvsave,-3))
c     vcurrent = sqrt(v_x**2+v_y**2+v_z**2)
c     Local velocities
         v_x = (ggg(ijkvsave,+1)+ggg(ijkvsave,-1))/ddxv(np1)
         v_x=v_x*x2(np1)-ggg(ijkvsave,-1)
         v_y = (ggg(ijkvsave,+2)+ggg(ijkvsave,-2))/ddyv(np1)
         v_y=v_y*y2(np1)-ggg(ijkvsave,-2)
         v_z = (ggg(ijkvsave,+3)+ggg(ijkvsave,-3))/ddzv(np1)
         v_z=v_z*z2(np1)-ggg(ijkvsave,-3)
         vcurrent = sqrt(v_x**2+v_y**2+v_z**2)
         
c     norm of velocity based on the random walk Div(D) term
c         
c         vdivd=0.
c         do iii=1,3
c            vdivd=vdivd+divd_omr(ijkvsave,iii)**2.
c         enddo
c         vdivd=sqrt(vdivd)
c     do not do random walk if Div(D) is greater than flow velocity
c         if(vdivd.gt.vcurrent) done=.true.
c         
         walk_counter = 0
         
         do while(.not.done)
            
            walk_counter = walk_counter + 1
            
            x2(np1) = x2save
            y2(np1) = y2save
            z2(np1) = z2save
            ijkv(np1) = ijkvsave
            ijkvs(np1) = ijkvssave
            ddxv(np1) = ddxvsave
            ddyv(np1) = ddyvsave
            ddzv(np1) = ddzvsave
            current_node = ijkvs(np1)
            upcoming_node= ijkv(np1)
            
c     Perform random walk for this particle
            
            if(upcoming_node.ne.0) then
c     if(upcoming_node.ne.0.and.
c     2        irray(upcoming_node,0).gt.0) then
               current_model = itrc(current_node)
               
               call random_walk_omr(1,np1,current_node,
     2              current_model ,edt,x_ts,y_ts,z_ts,itime,istep
     3              ,time_step_flag)
               
c...  check if near a capture node, 1/22/01 s kelkar 
               inp1=ijkv(np1)
               if(irray(inp1,0).eq.-inp1.or.irray(inp1,0).
     1              lt.-10000000) then
                  dx=ddx(inp1)
                  dy=ddy(inp1)
                  
                  x=x1(np1)-0.5*dx
                  y=y1(np1)-0.5*dy
                  xx=x2(np1)-0.5*dx
                  yy=y2(np1)-0.5*dy
                  
                  call check_well_capture(inp1,x,y,xx,yy,icapture,
     1                 vx_zheng,vy_zheng,vz_zheng,rw)               
                  if(icapture.eq.1) then
c  particle has entered the extraction well, remove it from the system
                     exit_node(np1)=-ijkv(np1)
                     call capture_particle(np1)
                     goto 9999                  
                  endif
               endif
               
c.................
               
c     s kelkar 3/29/05 random walk on OMR grids. use directional_search
               
               inside = .false.
               ijkv_last=ijkv(np1)
               max_level=10

               call directional_search(max_level,ijkv_last,np1,
     1              x1(np1),y1(np1),z1(np1),x2(np1),y2(np1),z2(np1),
     2              xc,yc,zc,ipc,flag_box,nout_face,save_out_face)

               if(flag_box.eq.ijkv_last) then
c still in same cc, no updating necessary
               done = .true.               
               elseif(flag_box.gt.0) then
c      the particle has entered a new cc, 
                  x2(np1)=x2(np1)+(corn(ijkv_last,1)-corn(flag_box,1))
                  y2(np1)=y2(np1)+(corn(ijkv_last,2)-corn(flag_box,2))
                  z2(np1)=z2(np1)+(corn(ijkv_last,3)-corn(flag_box,3))
                  oldnode_rand(np1)=ijkv(np1)
                  ijkv(np1)=flag_box
                  ijkv_last=flag_box
                  ipc_save(np1)=0
                  inside = .true.
               elseif(flag_box.lt.0) then
c     particle outside model domain, reflect back
                  call reflect_particle(inp1,np1,x2(np1),y2(np1),
     $                 z2(np1), nout_face,save_out_face)

                  call directional_search(max_level,ijkv_last,np1,
     1                 x1(np1),y1(np1),z1(np1),x2(np1),y2(np1),z2(np1),
     2                 xc,yc,zc,ipc,flag_box,nout_face,save_out_face)

                  if(flag_box.eq.ijkv_last) then
c     still in same cc, no updating necessary

                  elseif(flag_box.gt.0) then
                   x2(np1)=x2(np1)+(corn(ijkv_last,1)-corn(flag_box,1))
                   y2(np1)=y2(np1)+(corn(ijkv_last,2)-corn(flag_box,2))
                   z2(np1)=z2(np1)+(corn(ijkv_last,3)-corn(flag_box,3))
                     oldnode_rand(np1)=ijkv(np1)
                     ijkv(np1)=flag_box
                     ijkv_last=flag_box
                     ipc_save(np1)=0
                     inside = .true.
                  else
c                     if (iptty .ne. 0) then
c                        write(iptty,*)'reflected particle outside ',
c     $                       'domain'
c                        write(iptty,*)'inp1,np1,x1,y1,z1=',
c     $                       ijkv(np1),np1,x1(np1),y1(np1),z1(np1)
c                        write(iptty,*)' will try another random step.'
c                     end if
c                     write(ierr,*)'reflected particle outside domain'
c                     write(ierr,*)'inp1,np1,x1,y1,z1=',
c     $                    ijkv(np1),np1,x1(np1),y1(np1),z1(np1)
c                     write(ierr,*)' will try another random step.'
                  endif
               elseif(flag_box.eq.0) then   
                  if (iptty .ne. 0) then
c     write(iptty,*)'in new_loc_omr, max_level=',max_level
                     write(iptty,*)'Cant find the CC, ',
     $                    'another random step.'
                     write(iptty,*)'inp1,np1,x2,y2,z2=',
     $                    ijkv(np1),np1,x2(np1),y2(np1),z2(np1)
                  end if
c     write(iptty,*)'STOP'
c     write(ierr,*)'in new_loc_omr, max_level=',max_level
                  write(ierr,*)'Cant find the CC, another random step.' 
                  write(ierr,*)'inp1,np1,x2,y2,z2=',
     $                 ijkv(np1),np1,x2(np1),y2(np1),z2(np1)
c     write(ierr,*)'STOP'
c     stop
               endif            
               
c...  check if near a capture node, 1/22/01 s kelkar 
               if(inside) then
                  inp1=ijkv(np1)
                  if(irray(inp1,0).eq.-inp1.or.irray(inp1,0).
     1              lt.-10000000) then
                     icapture=1
c  particle has entered the extraction well, remove it from the system
                     exit_node(np1)=-ijkv(np1)
                     call capture_particle(np1)
                     goto 9999
                  endif
                  
c     Now determine if the particle has artifically jumped into 
c     slow moving fluid
c     First, skip this section if the tree-search routine
c     returned too many jumps for the particle, because we want to
c     try over
 2000             continue
                  v_x = (ggg(ijkv(np1),+1)+ggg(ijkv(np1),-1))/
     2                 ddxv(np1)
                  v_x=v_x*x2(np1)-ggg(ijkv(np1),-1)
                  v_y = (ggg(ijkv(np1),+2)+ggg(ijkv(np1),-2))/
     2                 ddyv(np1)
                  v_y=v_y*y2(np1)-ggg(ijkv(np1),-2)
                  v_z = (ggg(ijkv(np1),+3)+ggg(ijkv(np1),-3))/
     2                 ddzv(np1)
                  v_z=v_z*z2(np1)-ggg(ijkv(np1),-3)
                  vnew = sqrt(v_x**2+v_y**2+v_z**2)
                  if(vcurrent.gt.epsilonv.and.vnew/vcurrent.ge.
     2                 abs(vratio(current_model))) then
c     Velocity in new cell is ok, we're done
                     done = .true.
                  end if
c     else
c     done = .true.
               end if
            end if
            
            if(.not.done.and.walk_counter.ge.10) then
c     Can't get a valid jump, leave particle where it is and proceed
               x2(np1) = x2save
               y2(np1) = y2save
               z2(np1) = z2save
               ijkv(np1) = ijkvsave
               ijkvs(np1) = ijkvssave
               done = .true.
            end if
            
         end do
         
 9999    continue
         
c     if  particle has made a random jump on the
c     previous call to new_part _loc_struct. In that case 
c     reset oldnode and oldnode2 to ijkv(np1) to avoid confusing
c     this with a bouncing particle
         If(ijkv_in.ne.ijkv(np1)) then
            oldnode(np1)=ijkv(np1)
            oldnode2(np1)=ijkv(np1)
         endif

      endif
         
      return
      
      end subroutine new_loc_omr

c...................................................


      subroutine reflect_particle(inp1,np1,x,y,z,
     $         nout_face,save_out_face)

      use comsptr

      implicit none

      integer inp1,np1,nout_face,save_out_face(6),i,ifc

      real*8 x,y,z,d

      do i=1,nout_face
         ifc=save_out_face(i)
         if(ifc.eq.-1) then
            x=-x
         elseif(ifc.eq.+1) then
            x=2.*ddx(inp1)-x
         elseif(ifc.eq.-2) then
            y=-y
         elseif(ifc.eq.+2) then
            y=2.*ddy(inp1)-y
         elseif(ifc.eq.-3) then
            z=-z
         elseif(ifc.eq.+3) then
            z=2.*ddz(inp1)-z
         endif
      enddo

      return

      end

c..........................................................
