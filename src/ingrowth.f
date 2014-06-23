      subroutine ingrowth(begin_time,end_time)
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
!D1 PURPOSE
!D1 
!D1 To calculate the ingrowth of radio nuclides.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 Initial implementation: 19-JUL-1997, Programmer: cli
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/ingrowth.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:08:36   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:08   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:40   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3  2.3.5 Cell-based particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!D4 Subroutine used to calculate the ingrowth of radio nuclides
!D4                                                            
!D4 nreac: # of decay-ingrowth reactions.                      
!D4 ori:   parent species                                      
!D4 obj:   daughter species				      
!D4 nsegs: # of segments in i th parent species                
!D4 lsport: starting position of decaying parent species in a  
!D4         segment                                            
!D4 nsport: starting and ending index of decaying species in   
!D4         a segment                                          
!D4 time1/2: starting and ending time for decay calculation    
!D4 ivdt: the inverse of particle injection time interval      
!D4 ndaughter: total number of daughter particles generated at 
!D4            the end of end_time                             
!D4 nprevd: total number of daughter particles generated at the
!D4         end of previous time step                          
!D4 ingrow: number of new daughter particles generated within  
!D4         this time step                                     
!D4 tmsport: recording the strting and ending time for a       
!D4          segment                                           
!D4 msort: a subroutine used to sort the particles of the newly
!D4        generated species in terms of ascending time        
!D4                                                            
!D4 This subroutine uses an integration to approximate the      
!D4 summation for ingrowth calculation by assuming that        
!D4 particles are injected at equal time intervals. Thus, the   
!D4 do loop                                                    
!D4        sum=sum+(1-exp(-kfact*1(t-ti)))                     
!D4 end do                                                     
!D4 where t is time and ti is the injection time of the i th   
!D4 particle, is replaced by                                   
!D4        sum=ivdt*{[t2-t1]+[exp(-kfact*t2)-exp(-kfact*t1)]}  
!D4 Therefore, the number of operations are greatly reduced    
!D4 from (# of particles) * (1+1+1+1+1) to approximately 8.     
!D4                                                            
!D4 For particle tracking purposes, the origin injection time  
!D4 of the parent species is passed to the daughter species. If
!D4 the species continues to decay, then particles of the new     
!D4 species will be sorted in ascending time for the decay     
!D4 calculation purpose, the related parameters will be        
!D4 calculated and stored.                                     
!D4                               7/19/1997, cli               
!D4 Modified the code to handle variable confactors for decay  
!D4 ingrowth calculations.
!D4 
!*******************************************************************************

      use comai
      use comdti
      use comdi
      use compart
      use comsk

      implicit none


      integer i,ip,id
      real*8 begin_time,end_time
      
      do 100 i=1,nreac
         ip=ori(i)
         id=obj(i)
         if(id.eq.-1)then
            call decay
         else
            sumdecayed(ip)=num_particles(id)
            call decayingrowth
            sumdecayed(ip)=num_particles(id)-sumdecayed(ip)
         endif
 100  continue

      return
      contains

c     decay and decayingrowth are internal subroutines
C     subroutine for decay calculation only
C     Begin subroutine decay

      subroutine decay
	
      implicit none
      integer ingrow,lstart,j,k,ndaughter,lns,lns1,lns2
      integer i2,j2
      real*8  time1,time2
      
      do 150 j=1,nsegs(ip)
	   lns1=aidex(ip)+j
	   lns=lns1-1
         ndaughter=0
         if((lsport(lns)).gt.nsport(lns1))goto 150
	   lns2=lns+lns
         time1=(end_time-tmsport(lns2))/86400.
         time2=(end_time-tmsport(lns2-1))/86400.
         if(time2.lt.0.)time2=0.
         if(time1.lt.0.)time1=0.
         ndaughter=ivdt(lns)*(time2-time1+
     &        (exp(-kfact(ip)*time2)-exp(-kfact(ip)*time1))/
     &        kfact(ip))+0.5
         ingrow=ndaughter-nprevd(lns)
         
         if(ingrow.gt.0)then
            lstart=lsport(lns)
            do 200 i2=1,ingrow
               if(box(lstart,ip).gt.0)box(lstart,ip)=0
c bhl_6/4/08
c zvd 7/1/08 Only change the filtered box if mass will be output
               if (abs(prnt_rst) .ge. 40 .and. ripfehm .eq. 1) then
                  if(box(lstart,ip).lt.0 .and. abs(box(lstart,ip)) .lt. 
     &                 ibox_offset) box(lstart,ip)=box(lstart,ip) 
     &                 - ibox_offset
               end if
c bhl_6/4/08
               lstart=lstart+1
 200        continue	  
            lsport(lns)=lstart
            nprevd(lns)=ndaughter
         endif
 150  continue

      return
      end	subroutine decay

C     subroutine for decay-ingrowth calculation
C     Begin subroutine decayingrowth

      subroutine decayingrowth

      implicit none 
      
c bhl_12/8/05
c     integer i2,ingrow,sumingrow,lstart,itmp,lnsd,lnsd1,lnsd2
      integer i2,ingrow,sumingrow,lstart,itmp,lnsd,lnsd1,lnsd2,itmp2
c bhl_12/8/05
      integer lidex,lns,lns1,lns2,j2,ndaughter,j,k
      integer np_parent_j,np_dghtr_j
      real*8 time1,time2
      
      lidex=1
      sumingrow=0
c bhl_12/8/05
      itmp2=0
c bhl_12/8/05
      do 150 j=1,nsegs(ip)
	   lns1=aidex(ip)+j
	   lns=lns1-1
         itmp=num_particles(id)
         ndaughter=0
         if(lsport(lns).gt.nsport(lns1))goto 150
         lns2=lns+lns
         time1=(end_time-tmsport(lns2))/86400.
         time2=(end_time-tmsport(lns2-1))/86400.
         if(time2.lt.0.)time2=0.
         if(time1.lt.0.)time1=0.
         ndaughter=ivdt(lns)*(time2-time1+
     &        (exp(-kfact(ip)*time2)-exp(-kfact(ip)*time1))/
     &        kfact(ip))+0.5
         ingrow=ndaughter-nprevd(lns)  
         if(ingrow.gt.0)then
            lstart=lsport(lns)
            do 200 i2=1,ingrow
c bhl_6/4/08
c zvd 7/1/08 Only change the filtered box if mass will be output
               if (abs(prnt_rst) .ge. 40 .and. ripfehm .eq. 1) then
                  if(box(lstart,ip).lt.0 .and. abs(box(lstart,ip)) .lt. 
     &                 ibox_offset) box(lstart,ip)=box(lstart,ip) 
     &                 - ibox_offset
               end if
c bhl_6/4/08
               if(box(lstart,ip).gt.0)then
                  itmp=itmp+1
c bhl_12/8/05
                  itmp2=itmp2+1
c bhl_12/8/05
                  theta(itmp,id)=theta(lstart,ip)
                  start_time(itmp,id)=start_time(lstart,ip)
                  frac_done(itmp,id)=frac_done(lstart,ip)
                  box(itmp,id)=box(lstart,ip)
c bhl_5/5/08
c zvd 6/9/08 changed prnt_rst flag to >= 40
                  if (abs(prnt_rst) .ge. 40 .and. ripfehm .eq. 1) then
                     bconf_sav(itmp,id)=bconfactor(lns)
                  endif
c bhl_5/5/08
                  timeleft(itmp,id)=timeleft(lstart,ip)
                  if(ncolsizes.gt.0) then
                     if(sizes(id).ne.0) then
                        partsize(itmp,sizes(id)) =
     2                       partsize(lstart,sizes(ip))
                     end if
                  end if

c s kelkar March 22 07
                  if(flag_col_irrev(ip)) then
                  if(flag_col_irrev(id))then
c                     np_parent_j=lstart-lsport(lns)+1
c                     np_dghtr_j=itmp-nprevd(lns)+1
                     rcoll_div(itmp,irrevs(divs(id))) =
     2                    rcoll_div(lstart,irrevs(divs(ip)))
                     ret_weight(itmp,irrevs(divs(id))) =
     2                    ret_weight(lstart,irrevs(divs(ip)))
                  else
                     ret_weight_daughter(itmp,divs_d(id)) =
     2                    ret_weight(lstart,irrevs(divs(ip)))

                  endif
                  endif
                  if(flag_col_daughter(ip).eq.1) then
                  if(flag_col_daughter(id).eq.1)then
                     ret_weight_daughter(itmp,divs_d(id)) =
     2                    ret_weight_daughter(lstart,divs_d(ip))

                  endif
                  endif

                  box(lstart,ip)=0
                  sumingrow=sumingrow+1
               endif
               lstart=lstart+1
 200        continue	  
            lsport(lns)=lstart
            nprevd(lns)=ndaughter
         endif
         if(sumingrow.gt.0)then
            num_particles(id)=itmp
c bhl_12/8/05
            num_part_mem(id)=num_part_mem(id)+itmp2
            itmp2=0
c bhl_12/8/05
            nsegs(id)=nsegs(id)+1
            lnsd1=nsegs(id)+aidex(id)            
            lnsd=lnsd1-1
	      lsport(lnsd)=nsport(lnsd)+1
            nsport(lnsd1)=num_particles(id)
            confactor(lnsd)=confactor(lns)
            bconfactor(lnsd)=bconfactor(lns)
            do k=1,nreac
               if(ori(k).eq.id)then
                  lnsd2=lnsd+lnsd
                  nprevd(lnsd)=0
                  tmsport(lnsd2-1)=begin_time
                  tmsport(lnsd2)=end_time
                  ivdt(lnsd)=sumingrow/(tmsport(lnsd2)-
     &                 tmsport(lnsd2-1))
                  ivdt(lnsd)=86400.*ivdt(lnsd)
               endif
            enddo
            sumingrow=0
         endif
 150  continue
      
      return
      end     subroutine decayingrowth
      end  subroutine  ingrowth
