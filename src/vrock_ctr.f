      subroutine vrock_ctr(iflg,ndummy)
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
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
CD1
CD1 PURPOSE
CD1
CD1 To calculate variable rock density and heat capacity.
CD1     
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Started 080417 George Zyvoloski (Consultant) 
CD2 Based on subroutine vcon for variable thermal conductivity     
CD2 gaz 080117
CDA
CDA REFERENCES
CDA
CDA
C**********************************************************************
CPS
C**********************************************************************

      use combi
      use comci
      use comdi
      use comei
      use comfi
      use comgi
      use comii
      use comdti
      use comai
      use comki
      implicit none
      
      integer iflg,ndummy,i,ivrc,mid,mi,it,itp, i1, i2, itbl,iparam
      integer ntblines_roc, lu, itrocd
      real*8 vr1,vr2,vr3,vr4, vr5, vr6, vr7, vr8, vr12, sqrsat, tmpPor
      real*8 t_dum, cpr_t1, cpr_t2, denr_t1, denr_t2
      real*8 ddenr1t,ddenr2t,dcpr1t,dcpr2t
      real*8 cprmi,dcprmit,denrmi,ddenrmit
      real*8 term_t1, term_t2, diff_term
      real*8 tl, tmelt0, tmelt1, tmelt2, tmeltdt, heat_latent
      character*200 chdum, file_flux
      integer max_vrock_model, max_vrock_tabl 
      parameter (max_vrock_model = 2000, max_vrock_tabl = 100)
      real*8  strd_vroc
      parameter (strd_vroc = 0.9)
      integer open_file

      logical null1

c read in data
      if(ivrock.ne.0) then
         if(iflg.eq.0) then
          allocate(ivroc(max_vrock_model), table_vroc(max_vrock_tabl))
          allocate(vroc1f(max_vrock_model), vroc2f(max_vrock_model))
          allocate(vroc3f(max_vrock_model), vroc4f(max_vrock_model))
          allocate(vroc5f(max_vrock_model), vroc6f(max_vrock_model))
          allocate(vroc7f(max_vrock_model), vroc8f(max_vrock_model))
          allocate(vroc9f(max_vrock_model),itroc(n0))
          allocate(ntable_vroc(max_vrock_model))
          allocate(tblindx_roc(max_vrock_tabl,2))
          allocate(ddenrt(n0),dcprt(n0))
          allocate(ivrn(n0))
          if(.not.allocated(urock)) allocate(urock(n0),durockt(n0))
            ntable_roc = 0 
            i=0
 10         continue
            read(inpt,'(a80)') wdd1
            if(.not. null1(wdd1)) then
               backspace inpt
               read(inpt,*) ivrc
               backspace inpt
               if(isalt.ne.0) then
                if(ivrc.ne.6) then
                 write(ierr,*) 
     &            'warning non salt vrock model entered '
                 if(iout.ne.0) write(iout,*) 
     &            'warning non salt vrock model entered '
                 if(iptty.ne.0) write(iptty,*) 
     &            'warning non salt vrock model entered ' 
                endif
               endif
               i=i+1
               if(ivrc .eq.  1) then
c constant model for rock density and heat capacity                   
                  read(inpt,*) ivroc(i), 
     &              vroc1f(i),vroc2f(i)              
               else if(ivrc .eq.  2) then
c linear expansion model for rock density and heat capacity                   
                  read(inpt,*) ivroc(i), 
     &              vroc1f(i),vroc2f(i),vroc3f(i),vroc4f(i),vroc5f(i), 
     &              vroc6f(i)    
               else if(ivrc .eq.  3) then
c linear model for rock density and quadratic model for heat capacity                     
                  read(inpt,*) ivroc(i), 
     &              vroc1f(i),vroc2f(i),vroc3f(i),vroc4f(i),vroc5f(i), 
     &              vroc6f(i),vroc7f(i)
               elseif(ivrc .eq. 4) then
c tablular data for density and heat capacity
                  ntable_roc = ntable_roc + 1
                  read(inpt,'(a)') chdum  
                   chdum=trim(chdum)
                   read(chdum,*) ivroc(i), file_flux
                   file_flux = trim(file_flux)
                   table_vroc(ntable_roc) = file_flux
                   ntable_vroc(i) = ntable_roc
                   lu = open_file(file_flux,'old')
                   call manage_rock_tables(0,lu,
     &               ntable_roc,0,0.d0,0.d0,0.d0,0.d0,0.d0)
                   
               else if(ivrc .eq.  5) then
c linear model for rock density and heat capacity, with melting  
c reference temperature is 20 C for solid phase
c  vroc1f(i) density  denr
c  vroc2f(i) derivative denr wrt temperature
c  vroc3f(i) specific heat capacity for solid phase   (Cps)
c  vroc4f(i) derivative Cps wrt temperature
c reference temperature is melt temperature for melt phase
c  vroc5f(i) specific heat capacity for liquid (melt) phase   (Cpl, set = Cps at tmelt)
c  vroc6f(i) derivative Cpl wrt temperature                   
c  vroc7f(i) melt temperature (tmelt0)
c  vroc8f(i) latent_heat
c  vroc9f(i) spead temperature                  

                   read(inpt,*) ivroc(i), 
     &              vroc1f(i),vroc2f(i),vroc3f(i),vroc4f(i),
     &              vroc6f(i),vroc7f(i),vroc8f(i),vroc9f(i)     
c the internal energy derivative wrt temperature must be positive  
c therfore set liquid cpr(tmelt) tp solid(tmelt)                   
                    vr3=vroc3f(i)                  
                    vr4=vroc4f(i)
                    tmelt0 = vroc7f(i)
                    cpr_t1 =(vr4*(tmelt0-20.)+vr3)  
                    vroc5f(i) = cpr_t1                                                                    
               else if(ivrc .eq.  6) then   
c tables and latent heat    
c tablular data for density and heat capacity
                  ntable_roc = ntable_roc + 1
                  read(inpt,'(a)') chdum  
                   chdum=trim(chdum)
                   read(chdum,*) ivroc(i),file_flux,vroc7f(i),vroc8f(i),
     &                  vroc9f(i)
                   file_flux = trim(file_flux)
                   table_vroc(ntable_roc) = file_flux
                   ntable_vroc(i) = ntable_roc
                   lu = open_file(file_flux,'old')
                   call manage_rock_tables(0,lu,ntable_roc,0,
     &                  0.d0,0.d0,0.d0,0.d0,0.d0)
               endif
            else
               if(i.eq.0) ivrc=0
               go to 20
            endif
            go to 10
c     read in nodal capillary type
 20         continue

            narrays = 1
            itype(1) = 4
            default(1) = 0
            macro = "vroc"
            igroup = 2
            
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(9), macro, igroup, ireturn,
     &           i4_1=ivrn(1:n0))
            
            macroread(9) = .TRUE.

        if(ntable_roc.ne.0) then
          ntblines_roc = tblindx_roc(ntable_roc,2) 
          allocate(roc_table(ntblines_roc,3))
           do i = 1, ntblines_roc
             roc_table(i,1:3) = temp_table(i,1:3)
           enddo

c xhua error#6631 
c A non-optional actual argument must be present when invoking a procedure with an explicit interface
c should have 4 integers followed by 5 reals 
c         call manage_rock_tables(3,0,0,0.d0,0.d0,0.d0,0.d0,0.d0) original
          call manage_rock_tables(3,0,0,0,0.d0,0.d0,0.d0,0.d0,0.d0)
        endif

            
         else if(iflg.eq.1) then

c load heat capacity
            do mi=1,n0
               it=ivrn(mi)
               if(it.ne.0) then
               	itp=ivroc(it)
               	if(itp.eq.1) then
c     constant rock density and heat capacity
                   tl = t(mi)
                   denr(mi) = vroc1f(it)
                   cpr(mi) = vroc2f(it)*energy_conv
                   ddenrt(mi) = 0.0
                   dcprt(mi) = 0.0
                   urock(mi) = denr(mi)*cpr(mi)*(tl)
                   durockt(mi) = denr(mi)*cpr(mi)                
               	elseif(itp.eq.2) then
c     linear variation with temperature
c     vroc1f=reference temperature,vroc2f=reference rock density
c     vroc3f=d(density)/d(temp) at reference conditions
                  tl = t(mi)
                  vr1=vroc1f(it)
                  vr2=vroc2f(it)
                  vr3=vroc3f(it)
                  denr(mi)=(vr3*(t(mi)-vr1) + vr2)
                  ddenrt(mi) = vr3
                  vr4=vroc4f(it)
                  vr5=vroc5f(it)
                  vr6=vroc6f(it)
                  cpr(mi)=(vr6*(t(mi)-vr4) + vr5)*energy_conv
                  dcprt(mi) = vr6*energy_conv
                  urock(mi) = denr(mi)*cpr(mi)*(tl)
                  durockt(mi) = denr(mi)*(dcprt(mi)*(tl)+
     &                 cpr(mi)) + ddenrt(mi)*cpr(mi)*(tl)
               	else if(itp.eq.3) then 
c     linear variation with temperature
c     vroc1f=reference temperature,vroc2f=reference rock density
c     vroc3f=d(density)/d(temp) at reference conditions
                  tl = t(mi)
                  vr1=vroc1f(it)
                  vr2=vroc2f(it)
                  vr3=vroc3f(it)
                  denr(mi)=(vr3*(t(mi)-vr1) +vr2)
                  vr4=vroc4f(it)
                  vr5=vroc5f(it)
                  vr6=vroc6f(it)
                  vr7=vroc7f(it)
                  cpr(mi)=(vr6*(t(mi)-vr4) +vr7*(t(mi)-vr4)**2 +vr5)*
     &               energy_conv
                  dcprt(mi) = (vr6 + 2.*vr7*(t(mi)-vr4))*energy_conv
                  urock(mi) = denr(mi)*cpr(mi)*(tl)
                  durockt(mi) = denr(mi)*(dcprt(mi)*(tl)+
     &                 cpr(mi)) + ddenrt(mi)*cpr(mi)*(tl)                  
               	else if(itp.eq.4) then
                    tl = t(mi)
                    ntable_roc = ntable_vroc(it) 
                   call manage_rock_tables(1,0,ntable_roc,10,
     &                tl,cpr(mi),dcprt(mi),denr(mi),ddenrt(mi) )
                  urock(mi) = denr(mi)*cpr(mi)*(tl)
                  durockt(mi) = denr(mi)*(dcprt(mi)*(tl)+
     &                 cpr(mi)) + ddenrt(mi)*cpr(mi)*(tl)                   
                  else if(itp.eq.5) then 
c internal energy formulation
c linear model for rock density and heat capacity, with melting  
c reference temperature is 20 C for solid phase
c  vroc1f(i) density  denr
c  vroc2f(i) derivative denr wrt temperature
c  vroc3f(i) specific heat capacity for solid phase   (Cps)
c  vroc4f(i) derivative Cps wrt temperature
c reference temperature is melt temperature for melt phase
c  vroc5f(i) specific heat capacity for liquid (melt) phase   (Cpl) reference value
c  vroc6f(i) derivative Cpl wrt temperature                   
c  vroc7f(i) melt temperature
c  vroc8f(i) latent_heat (aready input in Mj/kg, typical value 1 to 2)
c  vroc9f(i) spead temperature          
                  tl = t(mi)
                  vr1=vroc1f(it)
                  vr2=vroc2f(it)
                  denr(mi)=vr2*(tl-20.) + vr1    
                  ddenrt(mi) = vr2
                  tmelt0 = vroc7f(it)  
c      center the temperature ramp                  
                  tmeltdt = vroc9f(it)/2. 
                  tmelt1 = tmelt0 - tmeltdt
                  tmelt2 = tmelt0 + tmeltdt                  
                  heat_latent = vroc8f(it)
                   if(tl.lt.tmelt1) then
                    vr3=vroc3f(it)                  
                    vr4=vroc4f(it)
                    cpr(mi)=(vr4*(tl-20.)+vr3)*energy_conv
                    dcprt(mi) = vr4*energy_conv                     
                    urock(mi) = 0.0 + denr(mi)*cpr(mi)*(tl)
                    durockt(mi) = denr(mi)*(dcprt(mi)*(tl)+
     &                 cpr(mi)) + ddenrt(mi)*cpr(mi)*(tl)
                     continue
                   else if(tl.gt.tmelt2) then
                    vr5=vroc5f(it)                  
                    vr6=vroc6f(it)
                    cpr(mi)=(vr6*(tl-tmelt0)+vr5)*energy_conv
                    dcprt(mi) = vr6*energy_conv                        
                    urock(mi)=denr(mi)*(heat_latent+cpr(mi)*
     &                         (tl))
                    durockt(mi)=denr(mi)*(dcprt(mi)*(tl)+cpr(mi))
     &                 + ddenrt(mi)*(cpr(mi)*(tl)+heat_latent)
                     continue
                   else
c average both density and internal energy 
                    vr3=vroc3f(it)                  
                    vr4=vroc4f(it)
                    cpr_t1 =(vr4*(tmelt1-20.)+vr3)*energy_conv   
                    denr_t1  = vr2*(tmelt1 -20.) + vr1
                    vr5=vroc5f(it)                  
                    vr6=vroc6f(it)                    
                    cpr_t2 =(vr6*(tmelt2-tmelt0)
     &                        +vr5)*energy_conv   
                    denr_t2  = vr2*(tmelt2 -20.) + vr1    
                    term_t1 = denr_t1*cpr_t1*(tmelt1)
                    term_t2 = denr_t2*(cpr_t2*(tmelt2)+
     &                      heat_latent)
                    diff_term = term_t2-term_t1
                    urock(mi)= term_t1 + diff_term*(tl-tmelt1)/
     &                 (2.*tmeltdt)   
                    durockt(mi) = diff_term/(2.*tmeltdt)
                    continue
                  endif                                         
                   
                else if(itp.eq.6) then 
c internal energy formulation with tables and latent heat
c  denrmi density  denr
c  ddenrmit derivative denr wrt temperature
c  cprmi specific heat capacity for solid phase   (Cps)
c  dcprmit derivative Cps wrt temperature                
c  vroc7f(i) melt temperature
c  vroc8f(i) latent_heat (aready input in Mj/kg, typical value 1 to 2
c  vroc9f(i) spread temperature difference   
                   tl = t(mi)
                   ntable_roc = ntable_vroc(it) 
                   call manage_rock_tables(1,0,
     &               ntable_roc,10,tl,denrmi,ddenrmit,cprmi,dcprmit)
c cpr and denr available from table 
                    denr(mi) = denrmi
                    ddenrt(mi) = ddenrmit
                    cpr(mi) = cprmi
                    dcprt(mi) = dcprmit
                  tmelt0 = vroc7f(it)  
c      center the temperature difference                  
                  tmeltdt = vroc9f(it)/2. 
                  tmelt1 = tmelt0 - tmeltdt
                  tmelt2 = tmelt0 + tmeltdt                  
                  heat_latent = vroc8f(it)
                   if(tl.lt.tmelt1) then
                    urock(mi) = 0.0 + denr(mi)*cpr(mi)*(tl)
                    durockt(mi) = denr(mi)*(dcprt(mi)*(tl)+
     &                 cpr(mi)) + ddenrt(mi)*cpr(mi)*(tl)
                     continue
                   else if(tl.gt.tmelt2) then                    
                     continue
                    urock(mi)=denr(mi)*(heat_latent+cpr(mi)*tl)
                    durockt(mi)=denr(mi)*(dcprt(mi)*(tl)+cpr(mi))
     &                 + ddenrt(mi)*(cpr(mi)*(tl)+heat_latent)           
                   else
c tmelt1 - tmelt - dt
c denr_t1 = rock density (tmelt1)
c ddenr1t = deriv wrt T denr(tmelt1)
c cpr_t1 = heat capacity (tmelt1) 
c dcpr1t = deriv wrt T heat capacity(tmelt1)
                    ntable_roc = ntable_vroc(it) 
                   call manage_rock_tables(1,0,
     &               ntable_roc,10,tmelt1,denr_t1,ddenr1t,cpr_t1,dcpr1t)
                 
                   call manage_rock_tables(1,0,
     &               ntable_roc,10,tmelt2,denr_t2,ddenr2t,cpr_t2,dcpr2t)
  
                    term_t1 = denr_t1*cpr_t1*(tmelt1)
                    term_t2 = denr_t2*(cpr_t2*(tmelt2)+
     &                      heat_latent)                    
                    diff_term = term_t2-term_t1
                    urock(mi)= term_t1 + diff_term*(tl-tmelt1)/
     &                 (2.*tmeltdt)   
                    durockt(mi) = diff_term/(2.*tmeltdt)
                    continue
                  endif
                 end if                   
                endif
              enddo
            else if(iflg.eq.2) then
c calculate rock phase state
              do mi = 1, n0 
               it=ivrn(mi)
                if(it.ne.0) then
                 itp=ivroc(it)
                 if(itp.ge.5) then
                  tmelt0 = vroc7f(it)  
                  tmeltdt = vroc9f(it)/2. 
                  tmelt1 = tmelt0 - tmeltdt
                  tmelt2 = tmelt0 + tmeltdt                        
                  tl = t(mi)
                  itrocd = itroc(mi)
                  if(tl.lt.tmelt1.and.itrocd.ne.1) then
                      itroc(mi) = 1
                      strd = min(strd,strd_vroc)
                  else if(tl.gt.tmelt2.and.itrocd.ne.3) then
                      itroc(mi) = 3
                      strd = min(strd,strd_vroc)                        
                  else if(tl.ge.tmelt1.and.tl.le.tmelt2.
     &                     and.itrocd.ne.2) then
                      itroc(mi) = 2
                      strd = min(strd,strd_vroc)                        
                 endif
                 endif
                endif                
              enddo
              else if(iflg.eq.3) then
c calculate rock phase state
              do mi = 1, n0 
               it=ivrn(mi)
                if(it.ne.0) then
                 itp=ivroc(it)
                 if(itp.ge.5) then
                  tmelt0 = vroc7f(it)                   
                  tmeltdt = vroc9f(it)/2. 
                  tmelt1 = tmelt0 - tmeltdt
                  tmelt2 = tmelt0 + tmeltdt                        
                  tl = t(mi)
                  if(tl.lt.tmelt1) then
                      itroc(mi) = 1
                  else if(tl.gt.tmelt2) then
                      itroc(mi) = 3                     
                  else if(tl.ge.tmelt1.and.tl.le.tmelt2) then
                      itroc(mi) = 2                      
                 endif
                 endif
                endif
              enddo
            end if
            end if
     
c      endif
      return
      end

      subroutine manage_rock_tables(iflg,table_unit,i_table,iparam,
     &           t_dum,var1_dum, dvar1_dumt, var2_dum, dvar2_dumt)
c 
c manage rock tables   
c     
c Data is found on the following lines
      use comai 
      use comdi
      implicit none

c     declare the 9 parameters
      integer iflg,table_unit,i_table,iparam 
      real*8 t_dum, dvar1_dumt,var1_dum,var2_dum, dvar2_dumt

c     variables
      real*8 var_dum, dvar_dumt 
      integer lu, lasttbl,i1,i2,mi,i
      integer ndx,cn,nparams,maxlines

      parameter(nparams = 3,maxlines = 100000)
      character*300 chdum
      
      if(iflg.eq.0) then
c read table input  
c table are appended  into one large table        
        if(.not.allocated(temp_table)) allocate(temp_table(maxlines,3))
c read 1 title line  
          read (table_unit,'(a)') chdum
           ndx = 0    
            if (i_table .eq. 1) then
c               first line is a text header (hence table data starts at line 2)
               tblindx_roc(i_table,1) = 2
               ndx = 0
            else
               tblindx_roc(i_table, 1) = tblindx_roc(i_table - 1, 2) + 1
               ndx = tblindx_roc(i_table - 1, 2)
            end if           
            do
               read (table_unit, '(a)', end = 5) chdum
c Input is terminated with a blank line or 'end' or end-of-file)
               if (len_trim(chdum) .eq.0 .or. 
     &               chdum(1:3) .eq. 'end') exit
               ndx = ndx + 1
               read (chdum, *) (temp_table(ndx,cn), cn = 1, nparams)
            end do
 5          lasttbl = i_table
            tblindx_roc(i_table, 2) = ndx
            if (table_unit .ne. inpt) close (table_unit)  
      else  if(iflg.eq.1) then    
c   Extract Tabular data
c iparam = 2: rock density  iparam = 3: rock specific heat 
c 
      i1 = tblindx_roc(i_table , 1)
      i2 = tblindx_roc(i_table , 2)
       if(iparam.ne.10) then
        do i = i1, i2 - 1   
         if (t_dum .le. roc_table(i, 1) .and. i .eq. i1) then
            var_dum = roc_table(i1, iparam)
            dvar_dumt = 0.
            go to 10
         else if (t_dum.ge.roc_table(i + 1, 1).and.i + 1 .eq. i2)then
            var_dum = roc_table(i2, iparam)
            dvar_dumt = 0.
            go to 10
         else if (t_dum .ge. roc_table(i, 1) .and. 
     &           t_dum .lt. roc_table(i + 1, 1)) then
            dvar_dumt =(roc_table(i + 1,iparam)-roc_table(i, iparam))/
     &           (roc_table(i + 1,1) - roc_table(i, 1))
            var_dum  = roc_table(i,iparam) + 
     &                 dvar_dumt* (t_dum - roc_table(i, 1))
            go to 10
         end if
        end do 
       else
        do i = i1, i2 - 1   
         if (t_dum .le. roc_table(i, 1) .and. i .eq. i1) then
            var1_dum = roc_table(i1, 2)
            dvar1_dumt = 0.
            var2_dum = roc_table(i1, 3)
            dvar2_dumt = 0.            
            go to 10
         else if (t_dum.ge.roc_table(i + 1, 1).and.i + 1 .eq. i2)then
            var1_dum = roc_table(i2, 2)
            dvar1_dumt = 0.
            var2_dum = roc_table(i2, 3)
            dvar2_dumt = 0.              
            go to 10
         else if (t_dum .ge. roc_table(i, 1) .and. 
     &           t_dum .lt. roc_table(i + 1, 1)) then
            dvar1_dumt =(roc_table(i + 1,2)-roc_table(i, 2))/
     &           (roc_table(i + 1,1) - roc_table(i, 1))
            var1_dum  = roc_table(i,2) + 
     &                 dvar1_dumt* (t_dum - roc_table(i, 1))
            dvar2_dumt =(roc_table(i + 1,3)-roc_table(i, 3))/
     &           (roc_table(i + 1,1) - roc_table(i, 1))
            var2_dum  = roc_table(i,3) + 
     &                 dvar2_dumt* (t_dum - roc_table(i, 1))            
            go to 10
         end if
        end do           
       endif
          return
10        if(iparam.eq.2) then
           var1_dum  = var_dum
           dvar1_dumt = dvar_dumt          
          else if(iparam.eq.3) then
           var2_dum  = var_dum*energy_conv
           dvar2_dumt = dvar_dumt*energy_conv 
          else if(iparam.eq.10) then            
           var2_dum  = var2_dum*energy_conv
           dvar2_dumt = dvar2_dumt*energy_conv           
          endif    

      else if (iflg.eq.3) then
c release memory for temp_table
       if(allocated(temp_table)) deallocate(temp_table)
      endif         
         return  
         end

