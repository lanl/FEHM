      subroutine directional_search(max_level,inp1,np1,
     1     xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,ipc,newnode,
     2     nout_face,save_out_face)  
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
!D1 To do a directional search to find the node a particle has moved
!D1 to.
!D1
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
c s kelkar July 14 05. 
c given inp1, xp1,yp1,zp1, Do a directional search to find the 
c newnode where xp2,yp2,zp2 lies and list all the brick shapes
c that are crossed in between. 

      use comai, only : ierr, iptty 
      use comsptr
      use comsk
      
      implicit none
      
      integer i,inp1,np1,newnode,ipc,nodep,itemp
      integer nplanes, iplane(3), j,max_level
      integer save_out_face(6), nout_face

      real*8 epsilon,xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc,dum,temp,xpg,ypg,zpg
      real*8 xpin,ypin,zpin,xpend,ypend,zpend

      max_level=10
      newnode=inp1
      nodep=inp1
      dum=0.

      if(xp1.eq.xp2) then
         if(yp1.eq.yp2) then
            if(zp1.eq.zp2) then
!               write(iptty,*)
!               write(iptty,*)'WARNING in DIRECTIONAL SEARCH'
!               write(iptty,*)'start and end point identical'               
!               write(ierr,*)
!               write(ierr,*)'WARNING in DIRECTIONAL SEARCH'
!               write(ierr,*)'start and end point identical'               
               newnode=inp1
               nout_face=0
               save_out_face=0
               goto 99999
            endif
         endif
      endif

      xpin=xp1
      xpend=xp2
      ypin=yp1
      ypend=yp2
      zpin=zp1
      zpend=zp2

      do j=1,max_level
c check if the end-point is in the cc of ip
c         call find_exit_planes(nodep,np1,xpend,ypend,zpend,nplanes,
c     1        iplane)
         xpg=xpend+corn(nodep,1)
         ypg=ypend+corn(nodep,2)
         zpg=zpend+corn(nodep,3)
         call check_inside_box(nodep,xpg,ypg,zpg,newnode,nplanes,
     $     iplane,nout_face,save_out_face)
         if(newnode.eq.nodep.or.newnode.lt.0) then
c     if newnode=nodep particle still in the same cc. if 
c     newnode<0 particle outside the model domain. Done.
            goto 99999
         elseif(newnode.eq.0) then
c     find the intersection (xc) of line (xp1)-(xp2) on the 
c     exit face. Exit face is calculated as that with shortest 
c     travel distance out of the cc of ip in the direction
c     (xp1)-(xp2)
c     also find new control volume for the point xc,yc,zc.
            call find_exit_point2(nodep,np1,nplanes,iplane,xpin,ypin,
     1           zpin,xpend,ypend,zpend,ipc,xc,yc,zc,newnode)
c     if newnode=0 particle has exited the system, 
            if(newnode.eq.0) then
               newnode=-nodep
               goto 99999
            endif
c     check if particle has crossed a breakthrough zone
c temporarily change ijkv(np1) because of the way it is used
c in seek_btczone, then swithc it back to old value
            itemp=ijkv(np1)
            ijkv(np1)=newnode
            call seek_btczone(np1, newnode)
            ijkv(np1)=itemp
c update (xpin),(xpend) and (xc) with reference to corn(newnode,-)
            temp=corn(nodep,1)-corn(newnode,1)
            xpend=xpend+temp
            xpin=xpin+temp
            xc=xc+temp
            temp=corn(nodep,2)-corn(newnode,2)
            ypend=ypend+temp
            ypin=ypin+temp
            yc=yc+temp
            temp=corn(nodep,3)-corn(newnode,3)
            zpend=zpend+temp
            zpin=zpin+temp
            zc=zc+temp
            nodep=newnode
         endif

      enddo
c particle jumped too far, try another random step
      newnode=0

99999 continue
         
      return

      end

c..............................................................


      subroutine find_exit_planes(inp1,np1,
     $      xp2,yp2,zp2,nplanes,iplane)

      use comsptr

      implicit none

      integer inp1,np1
      integer iplane(3),nplanes

      real*8 x2dim,y2dim,z2dim
      real*8 xp2,yp2,zp2

      x2dim = xp2/ddxv(np1)
      y2dim = yp2/ddyv(np1)
      z2dim = zp2/ddzv(np1)
      
c     ... s kelkar  march 29 04........................................
c     find how many control voulme faces are crossed by the line segme
      nplanes = 0
      iplane = 0
      if(x2dim.gt.1.) then
         nplanes=nplanes+1
         iplane(nplanes) = +1
      elseif(x2dim.lt.0.) then
         nplanes=nplanes+1
         iplane(nplanes) = -1
      endif
      if(y2dim.gt.1.) then
         nplanes=nplanes+1
         iplane(nplanes) = +2
      elseif(y2dim.lt.0.) then
         nplanes=nplanes+1
         iplane(nplanes) = -2
      endif
      if(z2dim.gt.1.) then
         nplanes=nplanes+1
         iplane(nplanes) = +3
      elseif(z2dim.lt.0.) then
         nplanes=nplanes+1
         iplane(nplanes) = -3
      endif

      return

      end subroutine find_exit_planes

c...............................................................

      subroutine find_exit_point2(inp1,np1,nplanes,iplane,xp1,yp1,zp1,
     1     xp2,yp2,zp2,ipc,xc,yc,zc,newnode)
      
c     find exit point(xc,yc,zc with reference to corn(inp1,)) 
c     and exit plane(ipc), and the new neighbour
c     xc is taken to be the intersection of the STRAIGHT line xp1-xp2
c     with the plane ipc
      
      use comsptr
      use comsk
      
      implicit none
      
      integer inp1,np1,icnl,newnode
      integer iplane(3),nplanes,ipc
      
      real*8 xp1,yp1,zp1,xp2,yp2,zp2
      real*8 xc,yc,zc
      
      newnode=0
      ipc=0
      
c     if nplanes > 0 then particle has entered a new cc
      if(nplanes.eq.1) then
         call exit_face1(inp1,np1,iplane(1),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     3        ipc,xc,yc,zc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.2) then
         call exit_face2(inp1,np1,iplane(1),iplane(2),
     1        xp1,yp1,zp1,xp2,yp2,zp2,
     3        ipc,xc,yc,zc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      elseif(nplanes.eq.3) then
         call exit_face3(inp1,np1,iplane(1),iplane(2),
     1        iplane(3), xp1,yp1,zp1,xp2,yp2,zp2,
     3        ipc,xc,yc,zc)
         call new_neighbour(inp1,np1,ipc,xc,yc,zc,newnode)
      endif
            
      return
      
      end subroutine find_exit_point2

c...........................................................
      

      subroutine exit_face1(i,np1,ip1,
     1     xp1,yp1,zp1,xp2,yp2,zp2,ipc,xc,yc,zc)
c     s kelkar march 29 04
c     find exit point accross the plane given by ip
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip,ip1,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1
      
      call face_intersect_line(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc1,yc1,zc1)
      
      xc=xc1
      yc=yc1
      zc=zc1
      ipc=ip1
      
      return
      
      end subroutine exit_face1
      
c...........................................................
      
      subroutine exit_face2(i,np1,ip1,ip2,
     1     xp1,yp1,zp1,xp2,yp2,zp2,ipc,xc,yc,zc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1, and ip2 and taking the closest of the two
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,d1,d2,dc

      call face_intersect_line(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc1,yc1,zc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call face_intersect_line(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc2,yc2,zc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
      endif
      
      return
      
      end subroutine exit_face2
      
c............................................................
      
      subroutine exit_face3(i,np1,ip1,ip2,ip3,
     1     xp1,yp1,zp1,xp2,yp2,zp2,ipc,xc,yc,zc)
      
c     s kelkar march 29 04
c     find  point of exit from the control volume of i
c     by finding intersection of line xp1-xp2 accross the planes 
c     given by ip1,ip2 and ip3 and taking the closest of the three
      
      use comsptr
      use comsk
      
      implicit none
      
      integer i,ip1,ip2,ip3,np1,ipc
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc,alamda,amue,anue,ap
      real*8 xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,d1,d2,d3,dc
      
      call face_intersect_line(i,np1,ip1,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc1,yc1,zc1)
      d1=sqrt((xp1-xc1)**2.+(yp1-yc1)**2.+(zp1-zc1)**2.)
      
      call face_intersect_line(i,np1,ip2,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc2,yc2,zc2)
      d2=sqrt((xp1-xc2)**2.+(yp1-yc2)**2.+(zp1-zc2)**2.)
      
      call face_intersect_line (i,np1,ip3,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc3,yc3,zc3)
      d3=sqrt((xp1-xc3)**2.+(yp1-yc3)**2.+(zp1-zc3)**2.)
      
      if(d1.le.d2) then
         xc=xc1
         yc=yc1
         zc=zc1
         ipc=ip1
         dc=d1
      else
         xc=xc2
         yc=yc2
         zc=zc2
         ipc=ip2
         dc=d2
      endif
      if(dc.gt.d3) then
         xc=xc3
         yc=yc3
         zc=zc3
         ipc=ip3
         dc=d3
      endif
      
      return
      
      end subroutine exit_face3
c......................................................................

      subroutine face_intersect_line (i,np1,ipc,xp1,yp1,zp1,xp2,yp2,zp2,
     1     xc,yc,zc)

c find the intersection of the straight line (xp1,yp1,zp1)-
c (xp2,yp2,zp2) with the plane of the face ipc
c xp1,xp2,xc etc are all refered to corn(current ndoe,-)
c s kelkar july 14,05

      use comai, only : ierr, iptty 
      use comsptr
      use comsk
      
      implicit none
      
      integer i,np1,ipc,ipc_sign
      real*8 xp1,yp1,zp1,xp2,yp2,zp2,xc,yc,zc
      real*8 epsilon,epsilon2
      real*8 dxp,dyp,dzp,dxc,dyc,dzc

      epsilon=+1.e-22
      epsilon2=+1.e-8
c epsilon2 is used as a fraction of original cell-size added
c to move the point slightly within the new cell

      if(abs(ipc).eq.+1) then
         if(ipc.eq.+1) then
            xc=(1.+epsilon2)*ddx(i)
         else
            xc=-epsilon2*ddx(i)
         endif
         dxp=(xp1-xp2)
         if(abs(dxp).gt.epsilon) then
            dxc=(xp1-xc)/dxp
            zc=zp1-dxc*(zp1-zp2)
            yc=yp1-dxc*(yp1-yp2)
         else
            write(iptty,*)'Error in face_intersect_line. stop'
            write(iptty,*)'cant have xp1=xp2 for abs(ipc)=1 case'
            stop
         endif
      elseif(abs(ipc).eq.+2) then
         if(ipc.eq.+2) then
            yc=(1.+epsilon2)*ddy(i)
         else
            yc=-epsilon2*ddy(i)
         endif
         dyp=(yp1-yp2)
         if(abs(dyp).gt.epsilon) then
            dyc=(yp1-yc)/dyp
            zc=zp1-dyc*(zp1-zp2)
            xc=xp1-dyc*(xp1-xp2)
         else
            write(iptty,*)'Error in face_intersect_line. stop'
            write(iptty,*)'cant have yp1=yp2 for abs(ipc)=2 case'
            stop
         endif
      elseif(abs(ipc).eq.+3) then
         if(ipc.eq.+3) then
            zc=(1.+epsilon2)*ddz(i)
         else
            zc=-epsilon2*ddz(i)
         endif
         dzp=(zp1-zp2)
         if(abs(dzp).gt.epsilon) then
            dzc=(zp1-zc)/dzp
            xc=xp1-dzc*(xp1-xp2)
            yc=yp1-dzc*(yp1-yp2)
         else
            write(iptty,*)'Error in face_intersect_line. stop'
            write(iptty,*)'cant have zp1=zp2 for abs(ipc)=3 case'
            stop
         endif
      endif

      return

      end

c..................................................................

