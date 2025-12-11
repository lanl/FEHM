subroutine setup_stressperm
!***********************************************************************
! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
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

  
  use comai
  use combi
  use comdi
  use comdti
  use comsi
  implicit none
  integer jj1,i,iispmd,i1,i2,jj, kb
  integer kbxmin,kbxmax,kbymin,kbymax,kbzmin,kbzmax
  real*8 xdmin,xdmax,ydmin,ydmax,zdmin,zdmax,xi,yi,zi
  real*8 xdmint,xdmaxt,ydmint,ydmaxt,zdmint,zdmaxt
  real*8 xkb,ykb,zkb,xd1,yd1,zd1,dis,dist,disx,disy,disz
  real*8 dis_tol
  
  parameter (dis_tol = 1.d-8)
  
  jj1 =0 
  do i = 1,n0
     iispmd = ispm(i)
     ispmd = ispmt(iispmd)
     if(ispmd.eq.3.or.ispmd.eq.5.or.ispmd.eq.7 &
          .or.ispmd.eq.8) then
        xi = cord(i,1)
        yi = cord(i,2)
        zi = cord(i,3)
        xdmax = -1.d30
        xdmin = -1.d30
        ydmax = -1.d30
        ydmin = -1.d30
        zdmax = -1.d30
        zdmin = -1.d30
        xdmaxt = -1.d30
        xdmint = 1.d30
        ydmaxt = -1.d30
        ydmint = 1.d30
        zdmaxt = -1.d30
        zdmint = 1.d30                  
        kbxmax = i
        kbxmin = i
        kbymax = i
        kbymin = i
        kbzmax = i
        kbzmin = i
        i1 = nelm(i)+1
        i2 = nelm(i+1)
        do jj = i1,i2
           kb = nelm(jj)
           xkb = cord(kb,1)
           ykb = cord(kb,2)
           zkb = cord(kb,3)
           if(xkb-xi.gt.0) then
              dis = abs(xkb-xi)
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              if(dis.gt.xdmax) then
                 kbxmax = kb
                 xdmax = dis
                 xdmaxt = dist
              else if(abs(dis-xdmax).lt.dis_tol) then
                 if(dist.lt.xdmaxt) then
                    xdmaxt = dist
                    kbxmax = kb
                 endif
              endif
           else if(xkb-xi.lt.0) then
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              dis = abs(xkb-xi)
              if(dis.gt.xdmin) then
                 kbxmin = kb
                 xdmin = dis
                 xdmint = dist
              else if(abs(dis-xdmin).lt.dis_tol) then
                 if(dist.lt.xdmint) then
                    xdmint = dist
                    kbxmin = kb
                 endif
              endif
           endif
           if(ykb-yi.gt.0) then
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              dis = abs(ykb-yi)
              if(dis.gt.ydmax) then
                 kbymax = kb
                 ydmax = dis
                 ydmaxt = dist
              else if(abs(dis-ydmax).lt.dis_tol) then
                 if(dist.lt.ydmaxt) then
                    ydmaxt = dist
                    kbymax = kb
                 endif
              endif
           else if(ykb-yi.lt.0) then
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              dis = abs(ykb-yi)
              if(dis.gt.ydmin) then
                 kbymin = kb
                 ydmin = dis
                 ydmint = dist
              else if(abs(dis-ydmin).lt.dis_tol) then
                 if(dist.lt.ydmint) then
                    ydmint = dist
                    kbymin = kb
                 endif
              endif
           endif
           if(zkb-zi.gt.0) then
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              dis = abs(zkb-zi)
              if(dis.gt.zdmax) then
                 kbzmax = kb
                 zdmax = dis
                 zdmaxt = dist
              else if(abs(dis-zdmax).lt.dis_tol) then
                 if(dist.lt.zdmaxt) then
                    zdmaxt = dist
                    kbzmax = kb
                 endif
              endif
           else if(zkb-zi.lt.0) then
              dist = abs(xkb-xi) + abs(ykb-yi) + abs(zkb-zi)
              dis = abs(zkb-zi)
              if(dis.gt.zdmin) then
                 kbzmin = kb
                 zdmin = dis
                 zdmint = dist
              else if(abs(dis-zdmin).lt.dis_tol) then
                 if(dist.lt.zdmint) then
                    zdmint = dist
                    kbzmin = kb
                 endif
              endif
           endif
        enddo
        !     
        !     store max and min of each direction (may create (slightly)larger stencil)
        !     
        ipermx(i,1) = kbxmin
        ipermx(i,2) = kbxmax
        ipermy(i,1) = kbymin
        ipermy(i,2) = kbymax	     
        ipermz(i,1) = kbzmin
        ipermz(i,2) = kbzmax	  
        if(cord(kbxmax,1)-cord(kbxmin,1).le.0.0) then
           jj1 =1
           write(ierr,*) & 
                'dis(x) failed, node ',i,' model 3 or 5 sub stres_perm'
        endif
        if(cord(kbymax,2)-cord(kbymin,2).le.0.0) then
           jj1 =1
           write(ierr,*) & 
                'dis(y) failed, node ',i,' model 3 or 5 sub stres_perm' 
        endif
        if(cord(kbzmax,3)-cord(kbzmin,3).le.0.0) then
           jj1 =1
           write(ierr,*) & 
                'dis(z) failed, node ',i,' model 3 or 5 sub stres_perm'
        endif
        !     stop for zero distances               
        if(jj1.ne.0) stop     
     endif
  enddo
  go to 2001
  !     
  !     test code
  !     
  open(99,file='stress_perm',status = 'unknown')
  do i = 1,n0
     iispmd = ispm(i)
     ispmd = ispmt(iispmd)
     if(ispmd.eq.3.or.ispmd.eq.5) then
        write(99,*) 'node = ',i
        write(99,*)  & 
             ' x node pair', ipermx(i,2), ipermx(i,1),  &
             ' dis ', cord(ipermx(i,2),1)-cord(ipermx(i,1),1)
        if(cord(ipermx(i,2),1)-cord(ipermx(i,1),1).le.0.001)  &
             write (99,*) '>>> zero distance x <<<<'
        write(99,*)  & 
             ' y node pair', ipermy(i,2), ipermy(i,1),  &
             ' dis ', cord(ipermy(i,2),2)-cord(ipermy(i,1),2)
        if(cord(ipermy(i,2),2)-cord(ipermy(i,1),2).le.0.001)  &
             write (99,*) '>>> zero distance y <<<<'
        write(99,*)  & 
             ' z node pair', ipermz(i,2), ipermz(i,1),  &
             ' dis ', cord(ipermz(i,2),3)-cord(ipermz(i,1),3)
        if(cord(ipermz(i,2),3)-cord(ipermz(i,1),3).le.0.001)  &
             write (99,*) '>>> zero distance z <<<<'    
     endif
  enddo
  close(99)
  
2001 continue   
  
  return
  
end subroutine setup_stressperm
!......................................................................
