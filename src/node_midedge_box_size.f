      subroutine node_midedge_box_size(
     1     x_range,y_range,z_range,
     2     xmin,   ymin,   zmin,
     3     xmax,   ymax,   zmax, i)
!***********************************************************************
! Copyright 2006 Los Alamos National Security, LLC  All rights reserved
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
!D1 Generate the bounding box of the midpoint of all edges connected to
!D1 node i.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 14-Sep-06, Programmer: C. Gable
!D2
!D2 $Log:   /pvcs.config/fehm90/src/node_midedge_box_size.f_a  $
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
C
C     INPUT: i - node number
C
C     OUTPUT: x_range,y_range,z_range
C             The range output is the bounding box
C             of the midpoint of all edges connected to
C             node i.
C
      use combi
      implicit none
      real*8 xi, yi, zi
      real*8 xj, yj, zj
      real*8 xmin, ymin, zmin
      real*8 xmax, ymax, zmax
      real*8 x_mid, y_mid, z_mid
      real*8 x_range, y_range, z_range
      integer i, j, n_connections, i1, i2, i_edge_j

      real*8 big_number
      data big_number / 1.e20 /

      integer if_debug
      data if_debug / 0 /
            
      xmin = big_number
      ymin = big_number
      zmin = big_number
      xmax =-big_number
      ymax =-big_number
      zmax =-big_number
      
      xi = cord(i,1)
      yi = cord(i,2)
      zi = cord(i,3)

      i1 = nelm(i)+1
      i2 = nelm(i+1)
C
C     Loop over all the connections to node i
C     Since it does not matter, the loop below will include
C     the i=j connection.
C
      n_connections = 0
      do j = i1,i2
         i_edge_j = nelm(j)

         if(i_edge_j .ne. i) then

         xj = cord(i_edge_j,1)
         yj = cord(i_edge_j,2)
         zj = cord(i_edge_j,3)
         
         x_mid = (xi + xj)/2.0d0
         y_mid = (yi + yj)/2.0d0
         z_mid = (zi + zj)/2.0d0
         
         xmin = min(xmin, x_mid)
         ymin = min(ymin, y_mid)
         zmin = min(zmin, z_mid)
         
         xmax = max(xmax, x_mid)
         ymax = max(ymax, y_mid)
         zmax = max(zmax, z_mid)

         n_connections = n_connections + 1

         if(if_debug .ge. 2)then
           write(6,*)'node i=', i, ' node j= ', i_edge_j
           write(6,*)xi, yi, zi, xj, yj, zj
           write(6,*)xmin,ymin,zmin,xmax,ymax,zmax
         endif

         endif

      enddo
      
      x_range = xmax - xmin
      y_range = ymax - ymin
      z_range = zmax - zmin
      
      if(if_debug .ge. 1)then
      write(6,*)'node =', i, ' # connections = ', n_connections
      write(6,*)' x min/max ', xmin, xmax, x_range
      write(6,*)' y min/max ', ymin, ymax, y_range
      write(6,*)' z min/max ', zmin, zmax, z_range
      endif
      
      return
      end
