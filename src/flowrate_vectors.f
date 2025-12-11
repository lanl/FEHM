      subroutine flowrate_vectors(iflg)  
c
c calculate flowrate vectors
c gaz initiated 112022
c
      use combi
      use comci
      use comdi
      use comei
      use comflow
      use davidi
      use comgi
      use comfi
      use comdti
      use comxi
      use comai
      use comco2, only : icarb, c_axy, c_vxy
      implicit none
      
      integer i, i1, i2, if, iff1, ii, iii, iz, j, jj, jjj, kb
c gaz 020723
c gaz 060323  dimensions changed 20 to 30
      integer  open_file, iflow, izone, nzones, kb_num, inc_flow, itestf
      integer zone_list(10)
      integer iflg,neqp1,nmatavw
      real*8  cosx, cosy, cosz, flow_zonex, flow_zoney
      real*8 disx,disy,disz, dis,dis2
      real*8 xi,yi,zi,xkb,ykb,zkb
      real*8 axy,vxy,flowratex,flowratey,flowratez,tol_flow,axy_t
      real*8 flow1(30)
      real*8 coskb(30,2),anorm
      integer i_kb(30),inorm 
      real*8 facflow 
      real*8, allocatable ::  flow2(:,:)
      parameter(tol_flow=1.d-14, facflow = 1.0, inorm = 1)
      save iflow  
      save flow2 
      save inc_flow
      nzones = 1
      itestf =0
c gaz 021324 turned off
      if(itestf.eq.0) return
      zone_list(1) = 5
      zone_list(2) = 6
      zone_list(3) = 7
      neqp1 = neq+1
      nmatavw = 0                        
      if(iflg.eq.-1) then
c gaz 020723 create  file for flux vectors output
       iflow = open_file('flow_vectors', 'unknown')
       if(.not.allocated(flow2)) allocate(flow2(neq,3))
       flow2 = 0.0
      else if(iflg.eq.-2) then
        close(iflow)
      else if(iflg.eq.1) then
       write(iflow,100)   days                                                 
100    format( 'flow vectors, time = ',1p,g14.5)
       do jj =  1, nzones
c       write(iflow,*)'zone ', zone_list(jj)
        flow_zonex = 0.0
        flow_zoney = 0.0
       do i = 1, neq
       izone = izonef(i)
c gaz 021423 force all nodes
       zone_list(jj)  = izone
       if(izone.eq.zone_list(jj))  then
        xi = cord(i,1)
        yi = cord(i,2)
        if(icnl.eq.0) then
         zi = cord(i,3)
        else
         zi = 0.0
        endif
        i1 = nelm(i) + 1
        i2 = nelm(i+1)
        flowratex = 0.0
        flowratey = 0.0  
        flowratez = 0.0  
        kb_num = 0   
        axy_t = 0.0  
        do j = i1,i2
         kb = nelm(j)
        if(i.ne.kb) then
         kb_num = kb_num+1
         xkb = cord(kb,1)
	   ykb = cord(kb,2)             
         disx = (xkb-xi)
         disy = (ykb-yi)    
        
          if(icnl.eq.0) then
           zkb = cord(kb,3)
           disz = (zkb-zi)
           dis2 = disx**2 + disy**2 + disz**2
          else
           zkb = 0.0
           dis2 = disx**2 + disy**2          
          endif
         dis =sqrt(dis2+1.d-12)
         cosx = (disx/dis)
         cosy = (disy/dis)
c gaz 060323 if heat conduction only, zero internodal flow
         if(idoff.lt.0) then
          axy = 0.0
         else
          axy = a_axy(j-neqp1+nmatavw)
         endif
c         axy = a_axy(j-neqp1+nmatavw)
         axy_t = axy_t +axy
         flow1(kb_num) = axy
         coskb(kb_num,1) = cosx
         coskb(kb_num,2) = cosy
         i_kb(kb_num) = kb                                
c         if(icnl.eq.0.and.axy.gt.tol_flow) then
         if(icnl.eq.0) then    
          cosz = (disz/dis)
          flowratex = flowratex + cosx*axy*facflow 
          flowratey = flowratey + cosy*axy*facflow   
          flowratez = flowratez + cosz*axy*facflow   
          flow2(i,1) =  cosx*axy*facflow 
          flow2(i,2) =  cosy*axy*facflow
          flow2(i,3) =  cosz*axy*facflow 
c         else if(axy.gt.tol_flow) then
         else
          flowratex = flowratex + cosx*axy*facflow  
          flowratey = flowratey + cosy*axy*facflow   
          flowratez = 0.0    
         endif
       endif
       enddo
        flow2(i,1) =  flowratex 
        flow2(i,2) =  flowratey 
c gaz 021523
        if(inorm.eq.0) then
          pnx(n+i) = flowratex 
          pny(n+i) = flowratey  
        else
          anorm =sqrt(flowratex**2+flowratey**2)+tol_flow                    
          pnx(n+i) = flowratex/anorm 
          pny(n+i) = flowratey/anorm  
        endif     
       flow_zonex = flow_zonex + flowratex
       flow_zoney = flow_zoney + flowratey
       if(itestf.ne.0) then
       write(iflow,101) i, izone, flowratex, flowratey, flowratez, 
     &                 xi, yi, zi        
101    format ('node ',i6,' zone ',i3,' flow vector x,y,z,',1p, 
     &         3(1x,g14.6),' coor x y z ',3(1x,g14.6))                           
       write(iflow,102) i, (i_kb(i1), i1 = 1,kb_num)
       write(iflow,104)'flow ',axy_t, (flow1(i1), i1 = 1,kb_num)
       write(iflow,103)'cosx*flow ',(flow1(i1)*coskb(i1,1),
     & i1 = 1,kb_num)
       write(iflow,103)'cosy*flow ',(flow1(i1)*coskb(i1,2),
     & i1 = 1,kb_num)
102    format ('node and neigh ',i6,t22, 10(1x,i14))
103    format (a10,t30,1p,10(1x,g14.6))
104    format (a5,1p,t8,g14.6,t30,10(1x,g14.6))
       endif
          write (iflow,106) i, flow2(i,1),flow2(i,2)
106    format(i9,1x,1p, g14.6,1x,g14.6)
       endif

       enddo
       if(itestf.ne.0) then
       write (iflow,105) zone_list(jj), flow_zonex, flow_zoney
       write (iflow,*) '*********************************************'
105    format('zone ',i5,' flow_zonex ',1p,g14.6,' flow_zoney ',g14.6)
       endif
       enddo
      else
      endif
      return
      end