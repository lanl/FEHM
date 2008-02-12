      subroutine inptrk(simnum)
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
CD1
CD1  PURPOSE
CD1
CD1  This subroutine is actived with particle tracking.  It is called
CD1  from subroutine input to read the data in the ptrk macro.
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/inptrk.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:10   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:26   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:58   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:02   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:52 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.11   Wed Jun 12 15:21:08 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.10   Mon Jun 03 11:18:04 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.9   Fri May 31 11:00:18 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.8   Fri Feb 16 09:46:00 1996   zvd
CD2 Modified requirement.
CD2 
CD2    Rev 1.7   Thu Jan 18 12:39:44 1996   hend
CD2 Updated Comments
CD2 
CD2    Rev 1.6   Wed Jan 10 13:06:38 1996   hend
CD2 Added Prolog
CD2 
CD2    Rev 1.5   09/06/95 10:00:58   robinson
CD2 Included checks to eliminate dispersivity values of 0
CD2 
CD2    Rev 1.4   04/25/95 09:44:22   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.3   03/15/95 17:05:06   robinson
CD2 Added diffusion and dispersion to particle tracking model
CD2 
CD2    Rev 1.2   02/03/95 10:09:52   llt
CD2 added logical null
CD2 
CD2    Rev 1.1   02/02/95 15:21:16   llt
CD2 defined null as a logical - required on IBM
CD2
CD2    Rev 1.0   01/28/95 14:00:22   llt
CD2 new particle tracking module
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.3.5 Cell-based particle-tracking module
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      use combi
      use comdi
      use compart
      use compfrac, only : ffmax
      use comdti
      use comai
      use comki
      use comsk
      use comsptr, only : tprpflag

c jj is used as a counter, dummy because a line of zero's will be taken
c as no input, whereas here it is valid input

      implicit none

      integer jj, jjj, transflag, irfile, i, j, numparams
      integer max1d,max1d2,lnwds,nwds,ii
      real dummy,half_life, r, ran_sp
      logical null1, afm
      character*3 string3
      logical used, ready
      character*80 dummy_string
      character*100 tfilename, pfilename
      real*8 ptrparam(500)
      real*8 simnum, confread, kf_avg
      integer isimnum, idpr, ninputs, set_value, open_file
      logical done, dfree_flag, nulldum
      integer icount, idummy, ipart
      real*8 rect, sum_rcoll, kf_max, kf_min, h2o_diff, tort_diff
      real*8 dret, dretlog, ret, retlog, retmin, retminlog, probret
      real*8 r_avg, retp, retn, drlog, rlog_avg, rlogmin, rlogmax
      real*8 , allocatable :: ret_weight_temp(:)
      real*8 , allocatable :: rcoll_div_temp(:)

      integer, parameter:: input_max=20,i2=2,i4=4
c bhl_11/3/06 
      integer, parameter::i5=5
c bhl_11/3/06         
      real, parameter::   add0_1=0.1

      integer imsg(input_max),msg(input_max)
      real*8  xmsg(input_max)
      character*32 cmsg(input_max)
      character*100 input_line

      integer li,lns,lns1,itemp,np_temp
      integer :: np_temp_max = 0

c  dummy is number of particles, already read in
 !cli to be compatible with V2.2, the code will use nstep when
 !max1d is not present in the input file

      msg=0
      imsg=0
      xmsg=0.
      lnwds=0
      max1d=0
      min_part=1

      read(inpt,'(a80)')input_line
      call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)

      do ii=1,nwds
         if(msg(ii) == 3)then
	    lnwds=lnwds+1
	    exit
         end if
      end do

      nwds=nwds-lnwds
      if(nwds <2 )then
         write(ierr,*)'Number of inputs is:',nwds
         write(ierr,*)'For ptrk macro input line'
         write(ierr,*)input_line
         write(ierr,*)'Program stopped due to input error'
         stop
      elseif(nwds == i2)then
         write(ierr,*)'Using FEHM V2.22 ptrk input format'
      elseif(nwds == i4)then
         write(ierr,*)'Using FEHM V2.23 MPTR input format'
c bhl_11/3/06 
      elseif(nwds == i5)then
         write(ierr,*)'Using FEHM V2.25 MPTR input format'
c bhl_11/3/06 
      else
         write(ierr,*)'Number of inputs: ',nwds
         write(ierr,*)'Expecting either 2,4 or 5 input parameters'
         write(ierr,*)'Please check the number of inputs for ',
     &        'ptrk macro'
         write(ierr,*)'Make sure the input values are correct' 
         write(ierr,*)'Code continues with the above warnings'
      end if

      dummy=imsg(1)+xmsg(1)
      rseed=imsg(2)+xmsg(2)
      if (nwds .ge. i4) then
         rseed_release=imsg(3)+xmsg(3)
         max1d = imsg(4)+xmsg(4)
      else
         rseed_release=rseed
         max1d=nstep
      end if
c bhl_11/3/06 
      if (nwds .eq. i5) then
         min_part = imsg(5)+xmsg(5)
         if (min_part .lt. 1) min_part=1
      end if
c bhl_11/3/06 

c*** water table rise modification
      read(inpt,'(a4)') macro
      if (macro.eq.'wtri') then
         read(inpt,*) water_table
         wtrise_flag = .true.
      else
         water_table = -1.d+10
         wtrise_flag = .false.
         backspace (inpt)
      end if
c*** water table rise modification

! initialize arrays

      nspeci=1
      nsp = 1
      if(.not.allocated(confactor)) then
         allocate(confactor(max1d),bconfactor(max1d), gmol(1))
         allocate(aidex(nspeci + 1))
         bconfactor = 1
         gmol = 1
         aidex(1) = 1
      end if

c     Check to see if rip option is being used, read in
c     input line accordingly

      read(inpt,'(a3)') string3
      if(string3(1:3) .eq. 'rip') then
         read (inpt, *) daycs,daycf,dayhf,dayhs, 
     2        ripfehm, confread,gmol(1)
         confactor = confread
      else
         backspace (inpt)
         read (inpt, *) daycs,daycf,dayhf,dayhs
         ripfehm = 0
         confactor = 1.
      end if

c bhl_11/3/06 (moved below option to read ripfehm)
! If this a GoldSIM simulation
      if (ripfehm .eq. 1) then
! Set random seeds to zero since seeds will be set by GoldSIM
         rseed = 0
         rseed_release = 0
      end if
c bhl_11/3/06 

c If the time when then heat and mass solution is enabled is
c less than the time when it was to be disabled set them =
      if (dayhs.lt.dayhf) then
         dayhs=dayhf
      endif

c Trak_type: 1 for liquid only, 2 for vapor only 
c Radioactive decay half life in seconds (kfact)
c Output type: 1 for #/total volume, 2 for #/kg fluid, 3 for #/kg vapor
c              anything else for # of particles
c Prnt_rst=0 for no restart info, *1 or *2 for yes
c prnt_rst also controls output to the .ptrk and .ptrk_fin files
c zvd 11/09/06 -- option to specify which parameters to output in .ptrk
c 1=Number Having Entered System,2=Number Currently In System,
c 3=Number Having Left System,4=Number Having Decayed,
c 5=Number Having Been Filtered,6=Number That Left This Time
      read (inpt, '(a80)') input_line
      read (input_line, *) trak_type(1), half_life, pout, prnt_rst
      if (abs(prnt_rst) .ge. 20) then
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
         if (nwds .gt. 4) then
            prnt_var = .false.
            do i = 5, nwds
               if (imsg(i) .gt. 0 .and. imsg(i) .lt. 7)
     &              prnt_var(imsg(i)) = .true.
            end do   
         else
            prnt_var = .true.
         end if
      else
         write (ierr, *) 'Particle statistics are no longer written to',
     &        ' the .out file.'
         write (ierr, *) 'Use new prnt_rst options if you need this ',
     &        'data.'
      end if

c     ZVD 10-31-2002 Recognize keyword for type curve data input and 
c        read in data if used
      read(inpt,'(a80)') dummy_string
      if(dummy_string(1:6).ne.'tcurve') then
         backspace inpt
      else
         read(inpt,'(a80)') input_line
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
         numparams = imsg(1)
         if (numparams .eq. 3 .and. nwds .gt. 1) then
            ffmax = xmsg(2)
         else
!Use default value for maximum fracture flow in 
         end if
         read(inpt,'(a100)') tfilename
         call read_tcurve (tfilename, numparams, ierr)
      end if

       read(inpt,'(a4)')macro
       if(macro.eq.'zptr')then
         read(inpt,*)ipzone
         ipzone1=ipzone+1
         if(idpdp.ne.0)then
 	   ipzone2=ipzone1+ipzone1
         else
	   ipzone2=ipzone1
         endif
         if(.not.allocated(pconcmax)) then
            allocate(pconcmax(n0))
         end if
         pconcmax=0
c     One species, allocate pcount as 1*ipzone1 = ipzone1

         if(.not.allocated(idzone)) then
            allocate(idzone(ipzone1),pcount(ipzone2,1))
            allocate(idcmax(ipzone2,1))
            allocate(idcavg(ipzone2,1))
            allocate(idflux(ipzone2))
         end if
         pcount = 0
         read(inpt,*)(idzone(idpr),idpr=1,ipzone)
       else
         ipzone = 0
         if(.not.allocated(pconcmax)) then
            allocate(pconcmax(1))
         end if
         pconcmax=0
         ipzone1 = 1
	 if(idpdp.ne.0)then
   	   ipzone2=ipzone1+ipzone1
	 else
	   ipzone2=ipzone1
	 endif
         if(.not.allocated(idzone)) then
            allocate(idzone(ipzone1),pcount(ipzone2,1))
         end if
         pcount = 0
         backspace(inpt)
       endif
c compute decay constant from half-life
      if (half_life .ne.0.) then
         kfact(1)=.6931472/half_life
      else
         kfact(1) = 0.
      end if

c     For pout option 5, need array for particle concentration

      if(pout.eq.5) then
         if(.not.allocated(partconc)) then
            allocate(partconc(num_particles(1)),nodepconc(n0))
         end if
      end if
      read(inpt,'(a80)') wdd1
      if(wdd1(1:4).eq.'file') then
         read(inpt,'(a100)') pfilename

         irfile = open_file(pfilename,'old')
         read(irfile,*) ninputs
c     Comment out to make consistent with inmptr - BAR 1-16-00
c            read(irfile,'(a4)') wdd1

         isimnum = int(simnum+0.01)

c     Read through the file until we get to the one we want for this realization
         do j = 1, isimnum
            read(irfile,*)(ptrparam(i),i=1,ninputs)
         end do
         close(irfile)
      else
         backspace inpt
      end if

      afm = .false.
      read(inpt,'(a80)') wdd1
      if(wdd1(1:3).eq.'afm') then
         afm = .true.
      else
         backspace inpt
      end if

c zvd 02/28/07 Add for free water diffusion coefficient and tortuosity
      dfree_flag = .false.
      read(inpt,'(a80)') wdd1
      if (wdd1(1:4).eq.'dfre') then
         dfree_flag = .true.
      else
         backspace inpt
      end if

c     Determine if particle size distribution is being used
c     If so, read in the distribution

      read(inpt,'(a80)') wdd1
      if(wdd1(1:4).eq.'size') then
         sizes(1) = 1
         jj = 0
 1991    continue
         read(inpt,'(a80)') wdd1
         if (null1(wdd1)) then
            goto 1992
         else
            jj = jj + 1
            read(wdd1,*) porsize(jj,1),probsize(jj,1)
            goto 1991
         end if
 1992    continue

         if(probsize(jj,1).ne.1.) then
            write(ierr,*)'Fatal error, the colloid particle size'
            write(ierr,*)'distribution must end at 1'
            if(iptty.ne.0) then
               write(ierr,*)'Fatal error, the colloid particle size'
               write(ierr,*)'distribution must end at 1'
            end if
            stop
         end if

c     Assign particle size for each particle

         do i = 1, max_particles

            done = .false.
            j=1
            r=ran_sp(rseed)
            do while(probsize(j,1).lt.r)
               j=j+1
            end do
c	Interpolate to get exact size
            partsize(i,1)=porsize(j-1,1)+
     2           (porsize(j,1)-porsize(j-1,1))*
     3           (r-probsize(j-1,1))/(probsize(j,1)-probsize(j-1,1))

         end do

      else
         backspace inpt
      end if

C     Hari 01-Nov-06 read in parameters for the colloid diversity model

      read(inpt,'(a80)') wdd1
      if(wdd1(1:4).eq.'dive') then
         flag_col_daughter(1) = 0
         divs_d(1) = 0
         flag_diversity(1) = .true.
         divs(1) = 1
         jj = 1

         read(inpt,'(a100)') input_line
         if (input_line(1:4) .eq. 'file') then
c     Hari read blank for filename since scanin read colloid 
c     distribution filename and opened file
            read (inpt,*)
         else
            backspace inpt
         end if
c     Use simnum that is passed in for the multiple realization
         isimnum = int(simnum+0.01)				
         if(isimnum<=0) isimnum = 1
         read(inpt,'(a100)') input_line
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
         tprpflag(1) = imsg(1)
         if (ripfehm .eq. 0 .and. nwds .ge. 2) isimnum = imsg(2)
c     Hari read in cdf of colloid retardation factors

c     Hari choose whether to read in cdf or equation
         if (iread_rcoll .ne. 0) rewind(iread_rcoll)
         if(tprpflag(1).eq.11.or.tprpflag(1).eq.12)then 
            read(iread_rcoll,'(a100)',end=19931) input_line
            read(iread_rcoll,*,end=19931) idummy
            do while (idummy.ne.isimnum)
               nulldum=.false.
               do while(.not.nulldum)
                  read(iread_rcoll,'(a100)',end=19931) input_line
                  nulldum=null1(input_line)
               enddo
               read(iread_rcoll,*,end=19931) idummy
            enddo
            jjj=0
            read(iread_rcoll,'(a100)',end=19931) input_line
            nulldum=null1(input_line)
            do while(.not.nulldum)
               jjj=jjj+1
               backspace(iread_rcoll)
               read(iread_rcoll,*,end=19931)rcdiv(jjj,divs(1)),
     2              probdiv(jjj,divs(1)) 
               read(iread_rcoll,'(a100)',end=19931) input_line
               nulldum=null1(input_line)
            enddo
            nprobdivs(divs(1)) = jjj
            k_rev(divs(1))=1.
            r_min(divs(1))=rcdiv(1,divs(1))
            r_max(divs(1))=rcdiv(nprobdivs(divs(1)),divs(1))
c     Hari importance sampling will be done at the time of breakthrough
c     for each particle
c     call impsample(tprpflag(divs(1)))
            goto 19932

19931       write(iptty,*)'Error in inmptr for tprpflag=',
     1           tprpflag(jj)
            write(iptty,*)'realization_num ',isimnum,
     1           ' not found in iread_coll file. STOP.'
            write(ierr,*)'Error in inmptr for tprpflag=',
     1           tprpflag(jj)
            write(ierr,*)'realization_num ',isimnum,
     1           ' not found in iread_coll file. STOP.'
            stop
            
19932       continue
            
         elseif(tprpflag(1).eq.13.or.tprpflag(1).eq.14) then
c     Reading an equation
            if (iread_rcoll .eq. 0) then
c     Read equations parameters from macro file
               read(inpt,*,end=19941) k_rev(divs(1)), r_min(divs(1)),
     1              r_max(divs(1)), slope_kf(divs(1)), cint_kf(divs(1)) 

            else
c     Read equation parameters from external file
               read(iread_rcoll,'(a100)',end=19941) input_line
               read(iread_rcoll,*,end=19941)idummy
               do while (idummy.ne.isimnum)
                  read(iread_rcoll,*,end=19941)idummy
               enddo
               backspace(iread_rcoll)
               read(iread_rcoll,*,end=19941)idummy,k_rev(divs(1)),
     1              r_min(divs(1)),r_max(divs(1)),
     2              slope_kf(divs(1)),cint_kf(divs(1)) 
            end if
c     Hari ptrparam option
c     If k_rev is negative we are getting all parameters from ptrparams
            if(k_rev(divs(1))<0.)then
               set_value=int(add0_1-k_rev(divs(1)))
               k_rev(divs(1))=ptrparam(set_value)
               set_value=int(add0_1-r_min(divs(1)))
               r_min(divs(1))=ptrparam(set_value)
               set_value=int(add0_1-r_max(divs(1)))
               r_max(divs(1))=ptrparam(set_value)
               set_value=int(add0_1-slope_kf(divs(1)))
               slope_kf(divs(1))=ptrparam(set_value)
               set_value=int(add0_1-cint_kf(divs(1)))
               cint_kf(divs(1))=ptrparam(set_value)
            end if  
            
c     s kelkar importance sampling done in part_track
c     for each time step , each species
c     call impsample(tprpflag(jj))
            goto 19942

19941       write(iptty,*)'Error in inptrk for tprpflag=',
     1           tprpflag(1)
            write(iptty,*)'realization_num ',isimnum,
     1           ' not found in iread_coll file. STOP.'
            write(ierr,*)'Error in inptrk for tprpflag=',
     1           tprpflag(1)
            write(ierr,*)'realization_num ',isimnum,
     1           ' not found in iread_coll file. STOP.'
            stop

19942       continue
         endif
         
         read(inpt,'(a80)') wdd1
         call parse_string(wdd1,imsg,msg,xmsg,cmsg,nwds)
         if(wdd1(1:4).eq.'reve')then
c     Hari calculate mean retardation factor for reversible colloid species
c     s kelkar flag_log=0 gives linear average, =1 is log average 
            if (nwds .ge. 2 .and. msg(2) .eq. 1) flag_log(1) = imsg(2)
            reves(divs(1)) = 1
            
            if(tprpflag(1).eq.11.or.tprpflag(1).eq.12)then
               sum_rcoll = 0.0
               do jj = 1, nprobdivs(divs(1))-1
                  if(flag_log(1).eq.0) then
                     rect = (probdiv(jj+1,divs(1))-
     2                    probdiv(jj,divs(1)))*(
     3                    rcdiv(jj+1,divs(1))+
     4                    rcdiv(jj,divs(1)))/2.
                  else
c     average of the log(R)
                     rect = (probdiv(jj+1,divs(1))-
     2                    probdiv(jj,divs(1)))*(
     3                    dlog10(rcdiv(jj+1,divs(1)))+
     4                    dlog10(rcdiv(jj,divs(1))))/2.
                  endif
                  sum_rcoll= sum_rcoll+rect
               enddo
               if(flag_log(1).eq.0) then
                  mean_rcoll_reve(reves(divs(1))) = sum_rcoll
               else
                  mean_rcoll_reve(reves(divs(1))) = 10.**sum_rcoll
               endif
            elseif(tprpflag(1).eq.13.or.tprpflag(1).eq.14)then
c     s kelkar do a quick numerical integration of Kf for the equation
c     log10(CDF)=b+m(log10(Kf)), ie 
c     prob for interval d(K_f) =[m*(K_f**(m-1))*(10.**b)]*d(K_f)
               kf_avg = 0.         
               if(flag_log(1).eq.0) then
                  dret=(r_max(divs(1))-r_min(divs(1)))*0.1
                  retmin=r_min(1)
               else
                  if(flag_log(1).eq.1) then
c     average of log(R) using the equation
                     rlog_avg=0.
                     rlogmin=dlog10(r_min(divs(1)))
                     rlogmax=dlog10(r_max(divs(1)))
                     drlog=(rlogmax-rlogmin)/100.
                     retp=r_min(divs(1))
                     do ii=1,100
c     using the value of R at the mid-point
                        retlog=rlogmin+drlog*(ii-0.5)
                        retn = 10**(rlogmin+drlog*ii)
                        ret=10.**retlog
                        dret = retn - retp
                        probret = dret * slope_kf(divs(1)) *
     1                       (10.**(cint_kf(divs(1)))) *
     2                       (ret-1)**(slope_kf(divs(1))-1.) *
     3                       (k_rev(divs(1))**slope_kf(divs(1)))
                        rlog_avg= rlog_avg+probret*retlog
                        retp = retn
                     enddo
                     r_avg = 10.**rlog_avg
                     mean_rcoll_reve(reves(divs(1))) = r_avg
                  else
c     linear average using the equation
                     r_avg = 0.
                     dret=(r_max(divs(1))-r_min(divs(1)))/100.
                     do ii=1,100
c     using the value of R at the mid-point
                        ret=r_min(divs(1))+dret*(ii-0.5)
                        probret= dret * slope_kf(divs(1)) *
     1                       (10.**(cint_kf(divs(1)))) *
     2                       (ret-1)**(slope_kf(divs(1))-1.) *
     3                       (k_rev(divs(1))**slope_kf(divs(1)))
                        r_avg= r_avg+probret*ret
                     enddo
                     mean_rcoll_reve(reves(divs(1))) = r_avg
                  endif
               end if
c     calculate mean of equation
            endif
         else
c....................................................................
c     4/2/07 s kelkar arrays for retardation -irreversible colloid 
c     4/11/07 sampling with equal divisions in the CDF space
c     diversity model, ret_weight and rcoll_div 
c     assign sampling method (equal_weight, random, importance)
            if (nwds .ge. 2 .and. msg(2) .eq. 3) then
               select case (cmsg(2)(1:2))
               case ('eq','EQ')
                  flag_method(1) = 1
               case ('ra','RA')
                  flag_method(1) = 0
               case ('im','IM')
                  flag_method(1) = 2
               end select
            end if
            irrevs(divs(1)) = 1
            flag_col_irrev(1) = .true.
            allocate(ret_weight(max_particles,1))
            allocate(rcoll_div(max_particles,1))
            if (ripfehm.eq.0) then
               if(flag_diversity(1).and.flag_col_irrev(1))then
                  np_temp=num_particles(1)
                  if (np_temp .gt. 0) then
                     if( np_temp .gt. np_temp_max) np_temp_max=
     &                    np_temp
                     allocate(rcoll_div_temp(np_temp))
                     allocate(ret_weight_temp(np_temp))
                     rcoll_div_temp=0.
                     ret_weight_temp=0.
                     call impsample_ptrk(tprpflag(1),1,
     1                    np_temp,rcoll_div_temp,
     2                    ret_weight_temp,np_temp_max)
                     itemp=0
                     do i=1, np_temp
                        itemp=itemp+1
                        ret_weight(i,irrevs(divs(1)))=
     1                       ret_weight_temp(itemp)
                        rcoll_div(i,irrevs(divs(1)))=
     1                       rcoll_div_temp(itemp)
                     enddo
                     deallocate (ret_weight_temp)
                     deallocate (rcoll_div_temp)
                  end if
               endif
            end if
         endif

      else
         backspace(inpt)
      endif
c     Hari end colloid diversity code
      jj = 0
 991  continue
      read(inpt,'(a80)') wdd1
      if (null1(wdd1)) then
         goto 99
      else
         jj = jj + 1
         backspace inpt
         read(inpt,*) transflag
         backspace inpt
         if(transflag.gt.0) then
c     No colloids
            if(afm .and. .not. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1),
     1              tclx(1,jj), tcly(1,jj), tclz(1,jj),
     2              diffmfl(1,jj), rd_frac(jj,1),
     3              matrix_por(jj),aperture(jj),
     4              sresidual(jj),gamma_afm(jj)
            else if(afm .and. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1),
     1              tclx(1,jj), tcly(1,jj), tclz(1,jj),
     2              h2o_diff, tort_diff, rd_frac(jj,1),
     3              matrix_por(jj),aperture(jj),
     4              sresidual(jj),gamma_afm(jj)
            else if(.not. afm .and. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1),
     1              tclx(1,jj), tcly(1,jj), tclz(1,jj),
     2              h2o_diff, tort_diff, rd_frac(jj,1),
     3              matrix_por(jj),aperture(jj)
            else
               read(inpt,*) transflag, kd(jj,1),
     1              tclx(1,jj), tcly(1,jj), tclz(1,jj),
     2              diffmfl(1,jj), rd_frac(jj,1),
     3              matrix_por(jj),aperture(jj)
            end if
            kcoll(jj,1) = 0. 
            rcoll(jj,1) = 1.
            fcoll(jj,1) = 1.
         else
c     Colloids
            if(afm .and. .not. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1), 
     2              tclx(1,jj), tcly(1,jj), tclz(1,jj), 
     2              diffmfl(1,jj), rd_frac(jj,1),
     3              matrix_por(jj), aperture(jj), kcoll(jj,1),
     3              rcoll(jj,1), fcoll(jj,1), sresidual(jj),
     4              gamma_afm(jj)
            else if (afm .and. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1), 
     2              tclx(1,jj), tcly(1,jj), tclz(1,jj), 
     2              h2o_diff, tort_diff, rd_frac(jj,1),
     3              matrix_por(jj), aperture(jj), kcoll(jj,1),
     3              rcoll(jj,1), fcoll(jj,1), sresidual(jj),
     4              gamma_afm(jj)
            else if (.not. afm .and. dfree_flag) then
               read(inpt,*) transflag, kd(jj,1), 
     2              tclx(1,jj), tcly(1,jj), tclz(1,jj), 
     2              h2o_diff, tort_diff, rd_frac(jj,1),
     3              matrix_por(jj), aperture(jj), kcoll(jj,1),
     3              rcoll(jj,1), fcoll(jj,1)
            else
               read(inpt,*) transflag, kd(jj,1), 
     2              tclx(1,jj), tcly(1,jj),
     2              tclz(1,jj), diffmfl(1,jj), rd_frac(jj,1),
     3              matrix_por(jj),aperture(jj),kcoll(jj,1),
     3              rcoll(jj,1), fcoll(jj,1)
            end if
            transflag=-transflag
         end if
         if(tclx(1,jj).lt.0.) then
            tclx(1,jj) = ptrparam(int(0.1+abs(tclx(1,jj))))
         end if  
         if(tcly(1,jj).lt.0.) then
            tcly(1,jj) = ptrparam(int(0.1+abs(tcly(1,jj))))
         end if  
         if(tclz(1,jj).lt.0.) then
            tclz(1,jj) = ptrparam(int(0.1+abs(tclz(1,jj))))
         end if  
         if (tclx(1,jj).lt.1e-30) tclx(1,jj)=1e-30
         if (tcly(1,jj).lt.1e-30) tcly(1,jj)=1e-30
         if (tclz(1,jj).lt.1e-30) tclz(1,jj)=1e-30
         if(aperture(jj).lt.0.) then
            aperture(jj) = ptrparam(int(0.1+abs(aperture(jj))))
         end if
         if (.not. dfree_flag) then
            if(diffmfl(1,jj).lt.0.) then
               diffmfl(1,jj) = ptrparam(int(0.1+abs(diffmfl(1,jj))))
            end if
         else
            if (h2o_diff .lt. 0. .and. tort_diff .lt. 0.) then
               diffmfl(1,jj) = ptrparam(int(0.1+abs(h2o_diff))) *
     &              ptrparam(int(0.1+abs(tort_diff)))
            else if (h2o_diff .lt. 0.) then
               diffmfl(1,jj) = ptrparam(int(0.1+abs(h2o_diff))) *
     &              tort_diff
            else if (tort_diff .lt. 0.) then
               diffmfl(1,jj) = h2o_diff * 
     &              ptrparam(int(0.1+abs(tort_diff)))
            else
               diffmfl(1,jj) = h2o_diff * tort_diff
            end if 
         end if
         if(kcoll(jj,1).lt.0.) then
            kcoll(jj,1) = ptrparam(int(0.1+abs(kcoll(jj,1))))
         end if  
         if(fcoll(jj,1).lt.0.) then
            fcoll(jj,1) = ptrparam(int(0.1+abs(fcoll(jj,1))))
         end if  
         if(rcoll(jj,1).lt.0.) then
            rcoll(jj,1) = ptrparam(int(0.1+abs(rcoll(jj,1))))
         end if  
         if(rd_frac(jj,1).lt.0.) then
            rd_frac(jj,1) = ptrparam(int(0.1+abs(rd_frac(jj,1))))
         end if  
         if(kd(jj,1).lt.0.) then
            kd(jj,1) = ptrparam(int(0.1+abs(kd(jj,1))))
         end if  
         if(matrix_por(jj).lt.0.) then
            matrix_por(jj) = ptrparam(int(0.1+abs(matrix_por(jj))))
         end if  
         if(afm) then
            if(sresidual(jj).lt.0.) then
               sresidual(jj) = ptrparam(int(0.1+abs(sresidual(jj))))
            end if  
            if(gamma_afm(jj).lt.0.) then
               gamma_afm(jj) = ptrparam(int(0.1+abs(gamma_afm(jj))))
            end if  
         end if

c     set transport simulation type

         set_transport_type:select case(transflag)
         case(1)
            dispflag(jj,1)=0
            diffflag(jj,1)=0
         case(2)
            dispflag(jj,1)=1
            diffflag(jj,1)=0
         case(3)
            dispflag(jj,1)=0
            diffflag(jj,1)=1
         case(4)
            dispflag(jj,1)=1
            diffflag(jj,1)=1
         case(5)
            dispflag(jj,1)=0
            if(afm) then
               diffflag(jj,1)=-2
            else
               diffflag(jj,1)=-1
            end if
         case(6)
            dispflag(jj,1)=1
            if(afm) then
               diffflag(jj,1)=-2
            else
               diffflag(jj,1)=-1
            end if
         case(7)
            if(afm) then
               diffflag(jj,1)=5
            else
               diffflag(jj,1)=3
            end if
            dispflag(jj,1)=1
         case(8)
            if(afm) then
               diffflag(jj,1)=-5
            else
               diffflag(jj,1)=-7
            end if
            dispflag(jj,1)=1
         end select set_transport_type
      endif
      
      goto 991

 99   continue

      narrays=1
      itype(1)=4
      default(1)=1
      macro="ptrk"
      igroup=5
      call initdata2(inpt,ischk,n0,narrays,itype,default,
     +     macroread(21),macro,igroup,ireturn,
     2     i4_1=itrc(1:n0))

      narrays=3
      itype(1)=8
      itype(2)=8
      itype(3)=8
      default(1)=0.
      default(2)=0.
      default(3)=0.
      macro="ptrk"
      igroup=6
      call initdata2(inpt,ischk,n0,narrays,itype,default,
     +     macroread(21),macro,igroup,ireturn,
     2     r8_1=pcnsk(1:n0),r8_2=t1sk(1:n0),r8_3=t2sk(1:n0))

c     Option for reading in particle concentrations
      if(pout.eq.5) then
         narrays=1
         itype(1)=8
         default(1)=0.
         macro="ptrk"
         igroup=7
         call initdata2(inpt,ischk,n0,narrays,itype,default,
     +        macroread(21),macro,igroup,ireturn,
     2        r8_1=nodepconc(1:n0))
      end if

      macroread(21) = .TRUE.

      ptrak=.true.
      restarting = .FALSE.
c     initialize variables
      box(1,1)=0
      npt(2)=n0
      npn=0
      max1d2=max1d+max1d
      if(.not.allocated(lsport)) then
         allocate(lsport(max1d+1),nsport(max1d+1))
         allocate(nsegs(nspeci+1),nprevd(max1d+1),astep(nspeci+1))
         allocate(ivdt(max1d+1),tmsport(max1d2+2))
      end if
      nsegs(1) = 1
      astep(1)=nstep
      lsport(1) = 1
      nsport(2)=num_particles(1)
      
      return
      end

