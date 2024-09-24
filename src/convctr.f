      subroutine convctr(iflg) 
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
!D1  PURPOSE
!D1
!D1  This subroutine manages the conversion from one type of physics to 
!D1  another
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 21-APR-03, Programmer: George Zyvoloski
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/convctr.f_a  $
!D2
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4  Equation variable
!D4
C!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
**********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comrxni
      use comii
      use davidi
      implicit none


c     
      integer iflg,icode, izone, inode, idir
      integer mi,neqp1,i,i1,i2,j,jj,ja,ii
      integer maxiconv, input_unit
      parameter (maxiconv=100)
      real*8 dist,dumconv,dumconv1,pdum,tdum,rolconv 
      real*8 head0_dum
c================================================================
      if(iconv.eq.0) return
c================================================================
      if(iflg.eq.0) then
c     
c     read input       
c     iconvf=1 : initial conditions         
c     iconvf=2 : boundary conditions         
c     conv1 is the reference pressure 
c     conv2 is the reference temperature
c     cordc is the reference coordinate
c     idirc is the gradient direction   
c     var1 is the linear temperature gradient
c     
         if(nconv.eq.0) then
            allocate(isconv(maxiconv))
            allocate(ifconv(maxiconv))
            read(inpt,*) nconv
            iconv = 1
            nall = n0*nconv
            isconv(iconv)=1
            ifconv(iconv)=nconv
            if(.not.allocated(izone_conv)) then
               allocate(izone_conv(max(1,nconv)))
               allocate(iconvf(max(1,nconv)))
               allocate(cordc(max(1,nconv)))
               allocate(varc(max(1,nconv)))
               allocate(conv1(max(1,nconv)))
               allocate(conv2(max(1,nconv)))
               allocate(headconv_val(max(1,nconv)))
               allocate(idirc(max(1,nconv)))
               allocate(izone_conv_nodes(nall))   
               izone_conv_nodes = 0
            end if
            headconv_val(iconv) = headconv
            read(inpt,*) 
     &           (izone_conv(i),iconvf(i),conv1(i),conv2(i),
     &           cordc(i),idirc(i),varc(i),i=1,nconv)

c     Loop over each zone for determining izone_conv array

            do izone = 1, nconv
               do inode = 1, n0
                  if(izonef(inode).eq.izone_conv(izone)) then
                     izone_conv_nodes(inode+n0*(izone-1)) =
     &               izone_conv(izone)
                  end if
               end do
            end do
           continue
         else 
            nconv0 = nconv
            read(inpt,*) nconv
            iconv = iconv+1
            nconv = nconv + nconv0
c gav 062223
            nall = n0*nconv
            isconv(iconv)=nconv0 + 1
            ifconv(iconv)=nconv

            allocate(izone_conv_temp(max(1,nconv)))
            allocate(iconvf_temp(max(1,nconv)))
            allocate(cordc_temp(max(1,nconv)))
            allocate(varc_temp(max(1,nconv)))
            allocate(conv1_temp(max(1,nconv)))
            allocate(conv2_temp(max(1,nconv)))
            allocate(headconv_val_temp(max(1,nconv)))
            allocate(idirc_temp(max(1,nconv)))
            allocate(izone_conv_nodes_temp(nall)) 
            

            izone_conv_temp(1:nconv0) = izone_conv(1:nconv0)
            iconvf_temp(1:nconv0) = iconvf(1:nconv0)
            cordc_temp(1:nconv0) = cordc(1:nconv0)
            varc_temp(1:nconv0) = varc(1:nconv0)
            conv1_temp(1:nconv0) = conv1(1:nconv0)
            conv2_temp(1:nconv0) = conv2(1:nconv0)
            headconv_val_temp(1:nconv0) = headconv_val(1:nconv0)
            idirc_temp(1:nconv0) = idirc(1:nconv0)
c gaz 062223
c            izone_conv_nodes_temp(1:nall-n0) = 
c     &           izone_conv_nodes(1:nall-n0)
            izone_conv_nodes_temp(1:nconv0*n0) =
     &      izone_conv_nodes(1:nconv0*n0)  
            deallocate(izone_conv,iconvf,cordc,varc,conv1)
            deallocate(conv2,idirc,izone_conv_nodes,headconv_val)    

            allocate(izone_conv(max(1,nconv)))
            allocate(iconvf(max(1,nconv)))
            allocate(cordc(max(1,nconv)))
            allocate(varc(max(1,nconv)))
            allocate(conv1(max(1,nconv)))
            allocate(conv2(max(1,nconv)))
            allocate(headconv_val(max(1,nconv)))
            allocate(idirc(max(1,nconv)))
            allocate(izone_conv_nodes(nall))   
            izone_conv_nodes(1:nconv0*n0) =
     &      izone_conv_nodes_temp(1:nconv0*n0)

            izone_conv = izone_conv_temp
            iconvf = iconvf_temp
            cordc = cordc_temp
            varc = varc_temp
            conv1 = conv1_temp
            conv2 = conv2_temp
            headconv_val = headconv_val_temp
            idirc = idirc_temp
c            izone_conv_nodes = izone_conv_nodes_temp

            deallocate(izone_conv_temp,iconvf_temp)
            deallocate(cordc_temp,varc_temp,conv1_temp)
            deallocate(conv2_temp,idirc_temp)
            deallocate(izone_conv_nodes_temp)
            deallocate(headconv_val_temp)  
            
            headconv_val(iconv) = headconv
            read(inpt,*) 
     &           (izone_conv(i),iconvf(i),conv1(i),conv2(i),
     &           cordc(i),idirc(i),varc(i),i=nconv0+1,nconv)

c     Loop over each zone for determining izone_conv array

            do izone = nconv0+1, nconv
               do inode = 1, n0
                  if(izonef(inode).eq.izone_conv(izone)) then
c gaz 062223 
c                     izone_conv_nodes(inode+n0*(iconv-1)) =                                               
                          izone_conv_nodes(inode+n0*(izone-1)) =
     &                    izone_conv(izone)
                  end if
               end do
            end do
         endif
         continue
      else if(iflg.eq.1.and.ico2.lt.0) then
c
c   modify initial values and BC's   
c   isothermal fluid physics 
c     
         if(iread.le.0) then
c no initial values from restart file
            do ii = 1,iconv
c           loop on conv calls
               do izone=isconv(ii),ifconv(ii)
c              loop on zones in each conv call
c gaz 110123 added 'conv'
                 call zone_saved(2,'conv',izone_conv(izone),izone,
     &           input_unit)
                  do inode=1,n0
c                 loop over nodes
                     if(izone_conv_nodes(inode+n0*(izone-1)).eq.
     &                    izone_conv(izone)) then
c                    identify node for conversion from head to pressure
                        idir = max(1,idirc(izone))
                        dist = cord(inode,idir)-cordc(izone)
                        headconv = headconv_val(izone)
                        if(iconvf(izone).eq.1) then
c                       check for initial value type
                           if(pho(inode).le.cord(inode,igrav)) then
c                          check for head le  elevation
                              pho(inode)= conv1(izone)
                              s(inode) = 0.
                           else
c                          convert head to pressure
                              ihead=1
                              dumconv = crl(1,1)
                              dumconv1 = crl(4,1)
                              pdum =conv1(izone)+rol0*headconv*(-grav)
                              tdum = conv2(izone)
                              call water_density(tdum,pdum,rolconv)
                              crl(1,1)=rolconv
                                 if(pres0.le.0.0) then
                                  crl(4,1)=pdum
                                 else
                                  crl(4,1)=pres0
                                 endif
                              head0_dum = head0
                              head0 = headconv
                              rho1grav = rolconv*9.81d-6
                              call headctr(5,inode,pho(inode),
     &                             pho(inode))
                              head0 = head0_dum
                              crl(1,1)= dumconv
                              crl(4,1)= dumconv1
                              ihead=0
                              if(varc(izone).ne.0) then
c                              can have variable temperature if it does not change
                                 to(inode) = tdum+ varc(izone)*dist
                              endif
                              s(inode) = 1.d0
c                          end IC check
                           endif
                        else if(iconvf(izone).eq.2) then
c                       check for head boundary value type
                           if(ka(inode).lt.0) then
                              if(pflow(inode).le.cord(inode,igrav)) then
c                             boundary value le elevation
                                 pflow(inode)= conv1(izone)
                                 esk(inode) = -1.
                                 ka(inode) = 0
                                 pho(inode) = conv1(izone)
                              else
                                 ihead=1
                                 dumconv = crl(1,1)
                                 dumconv1 = crl(4,1)
c gaz 052423
c                                 headconv = pflow(inode)-
c     &                           cord(inode,igrav)
                                  headconv = 0.0
                                 pdum =conv1(izone)+rol0*headconv*
     &                                (-grav)
                                 tdum = conv2(izone)
                                 call water_density(tdum,pdum,rolconv)
                                 crl(1,1)=rolconv
                                 if(pres0.le.0.0) then
                                  crl(4,1)=pdum
                                 else
                                  crl(4,1)=pres0
                                 endif
                                 
                                 head0_dum = head0
                                 head0 = headconv
                                 rho1grav = rolconv*9.81d-6
                                 call headctr(5,inode,pflow(inode),
     &                                pflow(inode))
                                 head0 = head0_dum
                                 crl(1,1)= dumconv
                                 crl(4,1)= dumconv1
                                 ihead=0
                              endif
c                          end BC check
                           endif
c                       end conv check
                        endif
c                    end identify node for conversion from head to pressure
                     endif
c                 end loop on nodes
                  enddo
c              end loop on zones for head conversion
               enddo
c           end loop on calls for conv
            enddo
c        end  conv section for no restart file
         else
c        restart file is used
            do ii = 1,iconv
c           loop on conv calls
               do izone=isconv(ii),ifconv(ii)
c              loop on zones
c gaz 110123 added 'conv'
                  call zone_saved(2,'conv',izone_conv(izone),izone,
     &            input_unit)
                  do inode=1,n0
c                 loop on nodes
                     if(izone_conv_nodes(inode+n0*(izone-1)).eq.
     &                    izone_conv(izone)) then
c                    check for zone match
                        idir = max(1,idirc(izone))
                        dist = cord(inode,idir)-cordc(izone)
                        headconv = headconv_val(izone)
c gaz 080523 intial values already in Mpa so skip iconvf(izone).eq.1 
                        if(iconvf(izone).eq.1.and.iread.eq.0) then
c                       check for ic or bc type
                           if(pho(inode).le.cord(inode,igrav)) then
                              pho(inode)= conv1(izone)
                              s(inode) = 0.
                           else
                              ihead=1
                              dumconv = crl(1,1)
                              dumconv1 = crl(4,1)
                              pdum =conv1(izone)+rol0*headconv*(-grav)
                              tdum = conv2(izone)
                              call water_density(tdum,pdum,rolconv)
                              crl(1,1)=rolconv
                                 if(pres0.le.0.0) then
                                  crl(4,1)=pdum
                                 else
                                  crl(4,1)=pres0
                                 endif
                              head0_dum = head0
                              head0 = headconv
                              rho1grav = rolconv*9.81d-6
                              call headctr(5,inode,pho(inode),
     &                             pho(inode))
                              head0 = head0_dum
                              crl(1,1)= dumconv
                              crl(4,1)= dumconv1
                              ihead=0
                              if(varc(izone).ne.0) then
                                 to(inode) = tdum+ varc(izone)*dist
                              endif
                              s(inode) = 1.d0
                           endif
                        else if(iconvf(izone).eq.2) then
                           if(ka(inode).lt.0) then
c                          check for bc nodes                     
                              if(pflow(inode).le.cord(inode,igrav)) then
                                 pflow(inode)= conv1(izone)
                                 esk(inode) = -1.
                                 ka(inode) = 0
                                 pho(inode) = conv1(izone)
                              else
                                 ihead=1
                                 dumconv = crl(1,1)
                                 dumconv1 = crl(4,1)
c gaz 052423
c                                 headconv = pflow(inode)-
c     &                           cord(inode,igrav)
                                  headconv = 0.0
                                 pdum =conv1(izone)+rol0*headconv*
     &                                (-grav)
                                 tdum = conv2(izone)
                                 call water_density(tdum,pdum,rolconv)
                                 crl(1,1)=rolconv
                                 if(pres0.le.0.0) then
                                  crl(4,1)=pdum
                                 else
                                  crl(4,1)=pres0
                                 endif
                                 
                                 head0_dum = head0
                                 head0 = headconv
                                 rho1grav = rolconv*9.81d-6
                                 call headctr(5,inode,pflow(inode),
     &                                pflow(inode))
                                 head0 = head0_dum
                                 crl(1,1)= dumconv
                                 crl(4,1)= dumconv1
                                 ihead=0
                              endif
c                          end bc  node check
                           endif
c                       end check for ic or bc
                        endif
c                    end check zone match
                     endif
c                 end loop on nodes
                  enddo
c              end loop on zones
               enddo
c           end loop on conv calls
            enddo
c        end restart check
         endif

         if(ico2.lt.0) then
            rho1grav = crl(1,1)*(9.81d-6)
         else
            rho1grav = rol0*9.81d-6
         endif

         if(allocated(izone_conv)) then
            deallocate(isconv)
            deallocate(ifconv)
            deallocate(izone_conv)
            deallocate(iconvf)
            deallocate(varc)
            deallocate(conv1)
            deallocate(conv2)
            deallocate(cordc)
            deallocate(idirc)
            deallocate(izone_conv_nodes)   
         end if

      else if(iflg.eq.1.and.ico2.eq.0) then
c     
c     modify initial values and BC's
c     
         if(iread.le.0) then
            do ii = 1,iconv
               do izone=isconv(ii),ifconv(ii)
c gaz 060523
c gaz 110123 added 'conv'
                 call zone_saved(2,'conv',izone_conv(izone),ii,jj)
                  do inode=1,n0
                     if(izone_conv_nodes(inode+n0*(ii-1)).eq.
     &                    izone_conv(izone)) then
                        idir = max(1,idirc(izone))
                        dist = cord(inode,idir)-cordc(izone)
                        headconv = headconv_val(izone)
                        if(iconvf(izone).eq.1) then
                           ihead=1
                           dumconv = crl(1,1)
                           dumconv1 = crl(4,1)
                           pdum =conv1(izone)+rol0*headconv*(-grav)
                           tdum = conv2(izone)
                           call water_density(tdum,pdum,rolconv)
                           crl(1,1)=rolconv
                           crl(4,1)=pres0
                           head0_dum = head0
                           head0 = headconv
                           rho1grav = rolconv*9.81d-6
                           call headctr(5,inode,pho(inode),pho(inode))
                           head0 = head0_dum
                           crl(1,1)= dumconv
                           crl(4,1)= dumconv1
                           ihead=0
                           if(varc(izone).ne.0) then
                              to(inode) = tdum + varc(izone)*dist
                           else
                              to(inode) = tdum       
                           endif
                        else if(iconvf(izone).eq.2) then
                           if(ka(inode).lt.0) then
                              ihead=1
                              dumconv = crl(1,1)
                              dumconv1 = crl(4,1)
                              pdum =conv1(izone)+rol0*headconv*(-grav)
                              tdum = conv2(izone)
                              call water_density(tdum,pdum,rolconv)
                              crl(1,1)=rolconv
                              crl(4,1)=pres0
                              head0_dum = head0
                              head0 = headconv
                              rho1grav = rolconv*9.81d-6
                              call headctr(5,inode,pflow(inode),
     &                             pflow(inode))
                              head0 = head0_dum
                              crl(1,1)= dumconv
                              crl(4,1)= dumconv1
                              ihead=0
                              if(varc(izone).ne.0) then
                                 esk(inode) = -(tdum+ varc(izone)*dist)
                              else
                                 esk(inode) = -tdum
                              endif
                           else if(ka(inode).gt.0) then
                              if(varc(izone).ne.0) then
                                 esk(inode) = -(tdum+ varc(izone)*dist)
                              else
                                 esk(inode) = -tdum
                              endif
                           endif
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            do ii = 1,iconv
               do izone=isconv(ii),ifconv(ii)
c gaz 060523
c gaz 110123 added 'conv'
                 call zone_saved(2,'conv',izone_conv(izone),0,0)
                  do inode=1,n0
                     if(izone_conv_nodes(inode+n0*(ii-1)).eq.
     &                    izone_conv(izone)) then
                        idir = max(1,idirc(izone))
                        dist = cord(inode,idir)-cordc(izone)
                        headconv = headconv_val(izone)
                        if(iconvf(izone).eq.2) then
                           if(ka(inode).lt.0) then
                              ihead=1
                              dumconv = crl(1,1)
                              dumconv1 = crl(4,1)
                              pdum =conv1(izone)+rol0*headconv*(-grav)
                              tdum = conv2(izone)
                              call water_density(tdum,pdum,rolconv)
                              crl(1,1)=rolconv
                              crl(4,1)=pres0
                              head0_dum = head0
                              head0 = headconv
                              rho1grav = rolconv*9.81d-6
                              call headctr(5,inode,pflow(inode),
     &                             pflow(inode))
                              head0 = head0_dum
                              crl(1,1)= dumconv
                              crl(4,1)= dumconv1
                              ihead=0
                              if(varc(izone).ne.0) then
                                 esk(inode) = -(tdum+ varc(izone)*dist)
                              else
                                 esk(inode) = -tdum
                              endif
                           else if(ka(inode).gt.0) then
                              if(varc(izone).ne.0) then
                                 esk(inode) = -(tdum+ varc(izone)*dist)
                              else
                                 esk(inode) = -tdum
                              endif
                           endif
                        endif
                     endif
                  enddo
               enddo
            enddo
         endif
    
         if(ico2.lt.0) then
            rho1grav = crl(1,1)*(9.81d-6)
         else
            rho1grav = rol0*9.81d-6
         endif
         
         if(allocated(izone_conv)) then
            deallocate(isconv)
            deallocate(ifconv)
            deallocate(izone_conv)
            deallocate(iconvf)
            deallocate(varc)
            deallocate(conv1)
            deallocate(conv2)
            deallocate(cordc)
            deallocate(idirc)
            deallocate(izone_conv_nodes)   
         end if

      else if(iflg.eq.2) then
c     
c     convert surfer output heads for temperature dependence
c     


      else if(iflg.eq.3) then
c     

      endif
c     
      return
      end                
      subroutine water_density(t,p,rol)
c     
c     calculate water density(assume (*,1) coefficients
c gaz 112819 hardwired (*,1) coefs to avoid possible memory issues      
c     
      use comii
      use comai, only : itsat
      implicit none
      real*8 t,p,pd,rol,prop,dpropt,dpropp
      real*8 rnwn1,rnwn2,rnwn3,rnwn
      real*8 rnwd1,rnwd2,rnwd3,rnwd
c gaz 112819 hardwired crl(,1)  coefs  
      real*8 crl01
      real*8 crl02
      real*8 crl05
      real*8 crl03
      real*8 crl06
      real*8 crl08
      real*8 crl04
      real*8 crl07
      real*8 crl10
      real*8 crl09
      real*8 crl11
      real*8 crl12
      real*8 crl15
      real*8 crl13
      real*8 crl16
      real*8 crl18
      real*8 crl14
      real*8 crl17
      real*8 crl20
      real*8 crl19

      real*8 pmin1
      real*8 pmax1
      real*8 tmin1
      real*8 tmax1
 

c
c**** density of liquid ****
c
c**** date created 10/31/88, max error(percent) 0.055 ****
c
      pmin1 =    0.001
      pmax1 =  110.000
      tmin1 =   15.000
      tmax1 =  360.000

      crl01 =  0.10000000d+01
      crl02 =  0.17472599d-01
      crl05 =  0.49564109d-02
      crl03 = -0.20443098d-04
      crl06 = -0.40757664d-04
      crl08 =  0.50330978d-04
      crl04 = -0.17442012d-06
      crl07 =  0.50676664d-07
      crl10 = -0.18383009d-06
      crl09 =  0.33914814d-06
      crl11 =  0.10009476d-02
      crl12 =  0.16812589d-04
      crl15 =  0.48841156d-05
      crl13 = -0.24582622d-07
      crl16 = -0.32967985d-07
      crl18 =  0.53249055d-07
      crl14 = -0.17014984d-09
      crl17 =  0.28619380d-10
      crl20 = -0.12221899d-09
      crl19 =  0.30456698d-09        
      if(itsat.le.10) then
c
c**** density of liquid ****
c
c**** date created 10/31/88, max error(percent) 0.055 ****
c
      pmin1 =    0.001
      if(p.le.pmin1) p = pmin1
      pmax1 =  110.000
      if(p.ge.pmax1) p = pmax1
      tmin1 =   15.000
      if(t.le.tmin1) t = tmin1
      tmax1 =  360.000
      if(t.ge.tmax1) t = tmax1
c      
       rnwn1=crl01+crl02*p+
     &     crl03*p*p+
     &     crl04*p*p*p
       rnwn2=crl05*t+crl06*t*t+
     &     crl07*t*t*t
       rnwn3=crl08*t*p+
     &     crl10*t*t*p+
     &     crl09*t*p*p
       rnwn=rnwn1+rnwn2+rnwn3
       rnwd1=crl11+crl12*p+
     &     crl13*p*p+
     &     crl14*p*p*p
       rnwd2=crl15*t+crl16*t*t+
     &     crl17*t*t*t
       rnwd3=crl18*t*p+
     &     crl20*t*t*p+
     &     crl19*t*p*p
       rnwd=rnwd1+rnwd2+rnwd3
       rol=rnwn/rnwd
      else
c density and derivatives          
          call eos_aux(itsat,t,p,0,1,prop,dpropt,dpropp)
          rol = prop
      endif
      return
      end
