      subroutine varchk_simple_awh(iflg,i1,i2,ndummy,dumawh1,
     &                             dumawh2,strd_simple)
c 
c simple phase change algorithm (two pass: thermal and material)
c
c gaz 091022 initial implementation
c
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use davidi      
      use com_exphase
      use commass_AWH
      use com_prop_data, only : pl_last, tl_last, pcl_last
      implicit none
     
      integer iflg,id1,id,i,j,i1,i2,ij,jmia,neqp1
      integer nr1,nr2,nr3 
      integer ieosd,ieosdc,ndummy
      integer i3, i4, kb
      real*8  dumawh1,dumawh2
      real*8 pl,tl,pcl,sl,strd_simple, step_nr
      real*8 pvapor,dpsatt,dpsats,dtsatp,psatl,tboil
      real*8 pv_diff,s_low,s_high,s_12,s_32,s_21,s_23
c gaz 030323 moved  phi_max  to comfi
      real*8 phi_fac, aiped_awh
      real*8 phi_temp,del1_min
      parameter (del1_min = 1.d-10)
      parameter (s_low=-1.d-2,s_high = 1.0)	   
      parameter (s_12=0.9,s_32 = 0.1,step_nr=0.98, aiped_awh = 1.d+3)	
      parameter (s_21 = 1.0, s_23 = 0.0, phi_fac = 10. )
   
      if(iflg.eq.1) then
       do id1 = i1,i2
        id = id1 + ndummy
        pl = phi(id)
        tl = t(id)
        pcl = pci(id)
        sl = s(id)
        ieosd = ieos(id)
        tboil=psatl(pl-pcl,0.0,0.0,dtsatp,
     &                             dpsats,1,0.0)
        pvapor = psatl(tl,0.,0.,dpsatt,dpsats,0,0.0)  
        pv_diff = (pl-pcl) - pvapor 
        ieosd = ieos(id)
        ieosdc = ieosd
c thermal phase test      
         if(ieosd.eq.1) then
c liquid phase 
          if(pv_diff.lt.0.0) then                                
           ieosdc = 2
           s(id) = s_12
          endif 
         elseif(ieosd.eq.2) then
c two phase 
c         s_high = 1.0
          if(sl.ge.s_high) then    
            ieosdc = 1
            if(pl-pcl.lt.pvapor) then
             pci(id) = pl - pvapor
             s(id) = 0.9999d0
             ieosdc = 2
            else                               
              pcl = pl-pvapor
              pci(id)=pcl 
              ieosdc=1
              s(id) = 1.0d0          
            endif           
           endif
           if(sl.le.s_low) then
            ieosdc = 3
            s(id) = s_23
          endif
         elseif(ieosd.eq.3) then
c gas phase 
          if(pv_diff.gt.0.0) then
           ieosdc = 2 
           s(id) = s_32
          endif
         elseif(ieosd.eq.4) then
c supercritical phase  
         endif
        pci(id)=max(0.0d00,pci(id))
        if(phi(id).lt.0.0d0) phi(id) = 0.0d0  
        if(ieosdc.ne.ieosd) then
         strd_simple=step_nr
         ieos(id) = ieosdc
        endif
       enddo
       continue
      else if(iflg.eq.2) then
c               
               do i=i1,i2
c
                  ieosd=ieos(i)
c gaz 032722                  
                  pl_last = phi(i)
                  tl_last = t(i)
                  pcl_last = pci(i) 
                  sl = s(i)
                  if(ps(i).eq.0.0) then
c     gaz 10-18-2001
                     t(i)=t(i)-del(2)*strd
c gaz 121718 added ieosd 4                     
                  elseif(ieosd.eq.1.or.ieosd.eq.4) then
                     phi(i)=phi(i)-del(1)*strd
                     t(i)=t(i)-del(2)*strd
                     pci(i)=pci(i)-del(3)*strd
                  elseif(ieosd.eq.2) then
                     phi(i)=phi(i)-del(1)*strd
                     s(i)=s(i)-del(2)*strd
c  GAZ 5/1/98          
                     t(i)=t(i)-del(3)*strd
c gaz 122418                      
                         pvapor=psatl(t(i),pcp(i),dpcef(i),dpsatt,
     &                           dpsats,0,an(i))  
c gaz debug 093019 stays in!                       
                     if(pvapor.ge.phi(i)) then
                          t(i)=psatl(phi(i),pcp(i),dpcef(i),
     &                          dtsatp,dpsats,1,an(i)) 
                        pvapor = phi(i)
                     endif     
                     pci(i) = phi(i)-pvapor
                  elseif(ieosd.eq.3) then
                     phi(i)=phi(i)-del(1)*strd
                     t(i)=t(i)-del(2)*strd
                     pci(i)=pci(i)-del(3)*strd
                  endif

                  pci(i)=max(0.0d00,pci(i))
                  s(i)=min(1.d00,max(0.0d00,s(i)))
c gaz 032722  load variable changes into bp
                  if(ieosd.ne.2) then
                   del(1) = phi(i)-pl_last
                   del(2) = t(i)-tl_last
                   del(3) = pci(i) - pcl_last                  
                  else
                   del(1) = phi(i)-pl_last
                   del(3) = t(i)-tl_last
                   del(2) = s(i) - sl                     
                  endif                  
               enddo 
c gaz 082622 added nr_completed = 3              
              nr_completed = 3

      else if(iflg.eq.3) then
c form src term with phi_max
        if(idelp.ne.0) then         
         dumawh1 = aiped_awh*(phi(i1)-phi_max)
         dumawh2 = aiped_awh
         idelp = 2
        endif
      else if(iflg.eq.4) then
c perform sucessive relaxation on neighbor nodes
c new residual for node id (with just updated variables)
       id = i1 + ndummy
       call awh_accumulation_calc(1,id,id,0)
       do id1 = i1,i2
        id = id1 + ndummy
c generate neighbor list and residuals
         i3 = nelmdg(id)+1
         i4 = nelm(id+1)
        do j = i3, i4
         kb = nelm(j)
         call awh_accumulation_calc(2,kb,kb,0)
        enddo
       enddo
      endif     
      return
      end

