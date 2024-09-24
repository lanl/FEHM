      subroutine insptr(ptitle)
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
!D1 This subroutine is activated with streamline particle tracking. It
!D1 is called from subroutine input to read the data in the sptr macro.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2 Initial implementation: 30-APR-1998, Programmer: G. Roemer
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/insptr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:22   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Replaced random_seed / random_number with ran_sp
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:02   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:27:10   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:04:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:43:00 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.6 Streamline particle-tracking module
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
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

      use combi
      use comdi
      use comdti
      use comai
      use comki
      use comsptr
      use comxi
      use compart
      use compfrac, only : ffmax
      use comsk

      implicit none
      integer i, iroot, idum, ifile, open_file, lf
      integer jj, flag_count, icliff1, icliff2, icliff3
      character*80 dummy_string, ptitle
      logical null1, done
      character*150 input_msg
      integer imsg(16)
      character*32 cmsg(16)
      integer msg(16)
      integer nwds
      real*8 xmsg(16)
      integer itabs
      integer nrandseeds, nzbtctmp
      integer, allocatable :: randseeds(:)
      integer, allocatable :: izonetmp(:,:), zbtctmp(:), tottmp(:)
      integer numparams
      integer numwell,ja,jb,jc,inode, nnwell
      real*8 rwell
      character*100 tfilename, root
      character*100 clfg1name, clfg2name, clfilename
      character*200 formstring
!      integer ncoef, max_con

c............... s kelkar Jan 10 07, for colloid diversity model
      integer tprpdum,jjj, realization_num
      logical nulldum, techead
c..........................................................


      macro = 'sptr'
      btc_flag = .false.
      exclude_particle = .false.
      output_end = .false.
c Assign file unit numbers and determine  output filenames
      isptr1 = nufilb(17)
      isptr2 = nufilb(18)
      isptr3 = nufilb(19)

c note that nufilb(20) amd (21) are already defined 
c and so are suffix(20) and (21)
c isptr7, is used for "sptr_save"
      isptr4 = nufilb(20)
      isptr5 = nufilb(21)
      isptr6 = nufilb(22)
      isptr7 = nufilb(23)
      isptr8 = nufilb(25)
      isptr9 = nufilb(26)
c.................................................

      root = ''
      if (null1(root_name)) then
! Find root of tracer file name or if not present, output file name
         if (nmfil(9) .ne. nmfily(3) .and. nmfil(9) .ne. ' ') then
c     Prefix from tracer file name
            call file_prefix(nmfil(9), iroot)
            if (iroot .gt. 94) iroot = 94
            root =  nmfil(9)(1:iroot)

         else if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
c     Prefix from output file name
            call file_prefix(nmfil(5), iroot)
            
            if (iroot .gt. 94) iroot = 94
            root = nmfil(5)(1:iroot)

         else 
c     Use default filenames "fehmn.sptr1", "fehmn.sptr2", "fehmn.sptr3", etc
            root = "fehmn"
            iroot = len_trim (root)
         endif
      else
         iroot = len_trim (root_name)
         if (iroot .gt. 94) iroot = 94
         root = root_name(1:iroot)
      end if

      if (.not. null1(root)) then
         nmfil(17) = root(1:iroot)
         nmfil(18) = root(1:iroot)
         nmfil(19) = root(1:iroot)
         nmfil(17)(iroot+1:iroot+6) = suffix(17)
         nmfil(18)(iroot+1:iroot+6) = suffix(18)
         nmfil(19)(iroot+1:iroot+6) = suffix(19)

c     s kelkar march 1 02.............................
         nmfil(20) = root(1:iroot)
         nmfil(21) = root(1:iroot)
         nmfil(22) = root(1:iroot)
!         nmfil(23) = root(1:iroot)
         nmfil(25) = root(1:iroot)
         nmfil(26) = root(1:iroot)
         nmfil(20)(iroot+1:iroot+6) = suffix(20)
         nmfil(21)(iroot+1:iroot+6) = suffix(21)
         nmfil(22)(iroot+1:iroot+6) = suffix(22)
!         nmfil(23)(iroot+1:iroot+6) = suffix(23)
         nmfil(25)(iroot+1:iroot+6) = suffix(25)
         nmfil(26)(iroot+1:iroot+6) = suffix(26)
      end if

c     Read in line, distinguish between 3 inputs (old way) and
c     more than 3. If there are more than 3, the fourth input
c     is the random number seed
c     BAR 12-15-2000
c      read(inpt,*) dtmx, iprt, iprto

! Group 1 data
      read (inpt,'(a)') input_msg
      call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)

c     assign dtmx (ensure that either a real or integer
c     input is read properly)
      if (msg(1).eq.1) then
         dtmx=imsg(1)
      else
         dtmx=xmsg(1)
      end if
! ZVD 05/05 - added for increasing btc time steps
!     set default value for dtmn
      dtmn = abs(dtmx)

c     assign iprt (only integers allowed)
      iprt = imsg(2)
c     assign iprto (only integers allowed)
      iprto = imsg(3)
c     Look for rseed, if it's there, assign it and
c     seed the random number generator
! ZVD eliminate random_seeds, implement use of ran_sp
      if(nwds.gt.3) then
         rseed = imsg(4)
! Can't have 0 for seed value in ran_sp
         if (rseed .eq. 0) rseed = 466201
c     get number of seeds the computer uses
c         call random_seed(size=nrandseeds)
c     allocate array of integer seed values
c         if(.not.allocated(randseeds)) 
c     2        allocate(randseeds(nrandseeds))
c     Assign values of array to input value
c         randseeds = rseed
c         call random_seed(put=randseeds)
      else
         rseed = 466201
      end if

! ZVD 06/05 added control for path information output
! 5th parameter is minimum time between output and 6th parameter
! is minimum distance (iprto > 0)
      if(nwds.gt.4) then
         if (msg(5).eq.1) then
            tplim=imsg(5)
         else
            tplim=xmsg(5)
         end if
      else
! Default no time limit for particle path output
         tplim = 0.
         ttpo = 0.
      end if
      if(nwds.gt.5) then
         if (msg(6).eq.1) then
            distlim=imsg(6)
         else
            distlim=xmsg(6)
         end if
      else
! Default no distance limit for particle path output
         distlim = 0.
      end if
      if(nwds.gt.6) then
         if (msg(7).eq.1) then
            linelim=imsg(7)
         else
            linelim=xmsg(7)
         end if
         if(linelim.eq.0) linelim = 1000000000
      else
! Default set line# limit for particle path output very large
         linelim = 1000000000
      end if
! weight for the Del(D) terms in random walk dispersion
      if(nwds.gt.7) then
         if (msg(8).eq.1) then
            divd_weight=1.*imsg(8)
         else
            divd_weight=xmsg(8)
         end if
      else
! Default set weight =1.
         divd_weight=1.
      end if

c Section to read keywords and associated data (independent of input
c order) ZVD 16-Oct-2006
! The courant factor line (Group 2) should be moved before keyword 
! input, but is parsed below to maintain compatibility with older 
! input decks

      done = .false.
      write_prop = 0
      flag_count = 0
      sptrx_flag = .false.
      xyz_flag = .false.
      alt_btc = .false.
      xyz_flag2 = .false.
      ip_flag = .false.
      trans_flag = .false.
      itensor = -999
      nzbtc = 0
      techead = .false.
      part_mult = 2.
      part_steps = 10
      part_frac = .1 * num_part
      delta_part = .0001 * num_part
      time_btc = 0.d0
      nplum = 0
      cliff_flag=.false.
      corner_flag = .false.
      epsilon_corner=0.0

! Default values for wtsi simulation
      smin=0.01
      deltawt=0.001

      do while(.not.done)
         done = .true.
         read(inpt,'(a80)') dummy_string
         read_flags:select case(dummy_string(1:2))
         case default
! The first non-character line read should contain courant factor, etc.
c     Read in line, distinguish between 2 inputs (old way) and
c     more than 2. If there are more than 2, the third input
c     is the flag for the dispersion tensor option
c     BAR 12-12-2000
            if (itensor .eq. -999) then
! We are reading Group 2 data
               done = .false.
               call parse_string(dummy_string, imsg, msg, xmsg, cmsg, 
     &              nwds)
c     assign courant_factor (ensure that either a real or integer
c     input is read properly)
               if (msg(1).eq.1) then
                  courant_factor=imsg(1)
               else
                  courant_factor=xmsg(1)
               endif
               if(courant_factor.eq.0) courant_factor = 0.25
c     assign iprtr (only integers allowed)
               iprtr = imsg(2)
c     assign itensor - set to input value (integer only) otherwise
c     if it's not in the input file set it to 0, it will be set to 5 
c     under keyword tprp if it has a zero value
               if(nwds.ge.3) then
                  itensor = imsg(3)
               else
                  itensor = 0
               end if

c......s kelkar march 27 2001  assign irevers- if it is not in the
c input file, set to 0, otherwise set to input value(integer only)
               if(nwds.ge.4) then
                  irevers = imsg(4)
               else
                  irevers = 0
               end if
c................................

c......s kelkarNov 13 2002  assign freez_time- if it is not in the
c input file, set to .0, otherwise set to input value(positive only)
c if freez_time is gt.0, then ptrac3 is called only at the end of the 
c flow calculations, and for a velocity frozen at that time,  
c particle tracks are calculated for freez_time days
               if(nwds.ge.5) then
                  if (msg(5).eq.1) then
                     freez_time  = imsg(5)
                  else
                     freez_time  = xmsg(5)
                  end if
               else
                  freez_time  = 0
               end if
            else if (itensor .eq. 0 ) then
! We are done with keywords and haven't ready any transport properties
! This is really old style input, so we needed to read it
               backspace (inpt)
c     Old style input
               read(inpt,*) ial,dispersivity1(1),
     2              dispersivity2(1),dispersivity3(1),
     2              dispersivity4(1),vratio(1)
               tprpflag(1) = 2
               random_flag=.true.
               kd(1,1) = 0.
               diffmfl(1,1) = 1.e-30
               rd_frac(1,1) = 1.
               matrix_por(1) = 1.e-30
               aperture(1) = -0.001
               itrc(1:n0) = 1
! Now make itensor 5
               itensor = 5
            else
               backspace (inpt)
            end if
c................................
         case('tc', 'TC')
c     ZVD 10-31-2002 Recognize keyword for type curve data input and 
c        read in data if used
            done = .false.
            read (inpt,'(a)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            numparams = imsg(1)
            if (numparams .eq. 3 .and. nwds .gt. 1) then
               ffmax = xmsg(2)
            else
!Use default value for maximum fracture flow in 
            end if
            read(inpt,'(a100)') tfilename
            call read_tcurve (tfilename, numparams, ierr)
         case ('om', 'OM')
c     ZVD 03-17-2004 Recognize keyword for grid with OMR
! Number of omr_nodes were read in scanin 
! (a value of zero represents a non-omr grid)
            done = .false.
            read(inpt, *) idum
! If a sptr restart is wanted (name, etc were read in scanin)
            if (save_omr) then
               read(inpt,'(a80)') dummy_string
               read(inpt,'(a80)') dummy_string
               read(inpt,'(a80)') dummy_string
            end if
         case ('tp', 'TP')
            done = .false.
c     ZVD 06-21-2005 Recognize keyword for transport porosities
            if (dummy_string(3:3) .eq. 'o' .or. 
     &           dummy_string(3:3) .eq. 'O') then
! We are reading transport porosities 'tpor'
               read(inpt,'(a80)') dummy_string
               if (dummy_string(1:4).eq.'file') then
                  read(inpt,'(a100)') tfilename
                  ifile = open_file(tfilename,'old')
               else
                  backspace inpt
                  ifile = inpt
               end if
            
!  zvd 06/02/03  Taken from rdcon
! Set default porosity to bogus value, assign all unassigned
! nodes the value set for rock macro after all data has been read 
               narrays=1
               itype(1) = 8
               macro = "sptr"
               igroup = 99
               default(1) = 10.
               call initdata2( ifile, ischk, n0, narrays, itype,
     2              default, macroread(21), macro, igroup, ireturn,
     3              r8_1=ps_trac(1:n0) )
               if (ifile .ne. inpt) close (ifile)
            else
! We are reading the transport property 'tprp' models
c     Section to find transport properties 
c     keyword and read in values if they are present
c     Note: the keyword tprp stands for "transport properties'
               if(itensor.eq.0) then
                  itensor = 5
                  if(iptty.gt.0) then
                     write(iptty,*)'***********************************'
                     write(iptty,*)' WARNING: in sptr macro input file,'
                     write(iptty,*)' no value was entered for itensor'
                     write(iptty,*)' defaulted to  itensor = 5'
                     write(iptty,*)'***********************************'
                  endif
                  write(ierr,*)'************************************'
                  write(ierr,*)' WARNING: in sptr macro input file,'
                  write(ierr,*)'no value was entered for itensor'
                  write(ierr,*)'defaulted to  itensor = 5'
                  write(ierr,*)'************************************'
               end if
               itabs = abs(itensor)

               jj = 0
               do
                  read(inpt,'(a150)') input_msg
                  if (null1(input_msg)) exit
                  jj = jj + 1
                  ial = 1
                  read(input_msg,*) tprpflag(jj)

                  if(tprpflag(jj).eq.1) then
                     read(input_msg,*) tprpflag(jj), kd(jj,1)
                     dispersivity1(jj) = 0.
                     dispersivity2(jj) = 0.
                     dispersivity3(jj) = 0.
                     dispersivity4(jj) = 0.
                     dispersivity6(jj) = 0.
                     dispersivity7(jj) = 0.
                     dispersivity8(jj) = 0.
                     vratio(jj) = 0.

                  elseif(tprpflag(jj).eq.2) then
                     random_flag=.true.
                     call parse_string(input_msg, imsg, msg, xmsg, 
     &                    cmsg, nwds)
                     read_disp1:select case(itabs)

c ...s kelkar aug 27 01  including case(0), otherwise, if itensor=0
c then itabs=0 and code gets stuck in an infinite do loop.It needs 
c to read this line so that it will go on to the next line
                     case(0)
! using input_msg so not backspacing
!                       read(inpt,*)
c....................................

                     case(1)
c     Modified axisymmetric form of Lichtner et al.
c     dispersivity1 = molecular diffusion
c     dispersivity2,3,4,6  =  alpha1 thru 4
c     dispersivity7 and 8 = direction cosines in x and y dir for 
c     the axis of symmetry 
                        call read2_itensor1

                     case(2)
c     Burnett and Frind
                        call read2_itensor2

                     case(3)
c     modified Burnett and Frind derived by Lichtner et al.
                        call read2_itensor3

                     case(4)
c     Isotropic tensor - Tompson et al.
                        call read2_itensor4

                     case(5)
c     Old version of fehm
                        call read2_itensor5

                     end select read_disp1
                     
                     diffmfl(1,jj) = 1.e-30
                     rd_frac(jj,1) = 1.
                     matrix_por(jj) = 1.e-30
                     aperture(jj) = -0.001

                  elseif(tprpflag(jj).eq.3) then
                     call parse_string(input_msg, imsg, msg, xmsg, 
     &                 cmsg, nwds)
c zvd 05/21/07 Added parameter for optional input of secondary spacing
                     if (nwds .ge. 7 .and. msg(7) .ne. 3) then
                        read(input_msg,*) tprpflag(jj), kd(jj,1), 
     2                       diffmfl(1,jj), rd_frac(jj,1),
     3                       matrix_por(jj),aperture(jj),secondary(jj)
                     else 
                        read(input_msg,*) tprpflag(jj), kd(jj,1), 
     2                       diffmfl(1,jj), rd_frac(jj,1),
     3                       matrix_por(jj),aperture(jj)
                     end if                        
                     dispersivity1(jj) = 0.
                     dispersivity2(jj) = 0.
                     dispersivity3(jj) = 0.
                     dispersivity4(jj) = 0.
                     dispersivity6(jj) = 0.
                     dispersivity7(jj) = 0.
                     dispersivity8(jj) = 0.
                     vratio(jj) = 0.

                  elseif(tprpflag(jj).eq.4) then
                     random_flag=.true.
                     call parse_string(input_msg, imsg, msg, xmsg, 
     &                 cmsg, nwds)
                     read_disp2:select case(itabs)

c ...s kelkar aug 27 01  including case(0), otherwise, if itensor=0
c then itabs=0 and code gets stuck in an infinite do loop.It needs 
c to read this line so that it will go on to the next line
                     case(0)
! using input_msg so not backspacing
!                        read(inpt,*)
c....................................

                     case(1)
c     Modified axisymmetric form of Lichtner et al.
c     dispersivity1 = molecular diffusion
c     dispersivity2,3,4,6  =  alpha1 thru 4
c     dispersivity7 and 8 = direction cosins in x and y dir for the axis
c     of symmetry 
                        call read4_itensor1
                     case(2)
c     Burnett and Frind
                        call read4_itensor2

                     case(3)
c     modified Burnett and Frind derived by Lichtner et al.
                        call read4_itensor3

                     case(4)
c     Isotrpoic tensor - Tompson et al.
                        call read4_itensor4

                     case(5)
c     Old version of fehm
                        call read4_itensor5

c     In the old version, vratio<0 was the flag to omit
c     the gradd term. Now set the itensor term to reflect that
                        if(vratio(jj).lt.0.) itensor = -5
                     end select read_disp2
c ...s kelkar jan 10 07
cHari 01-Nov-06 implement colloid diversity model

               elseif(tprpflag(jj).eq.11.or.tprpflag(jj).eq.12) then
                  div_flag=.true.
                  read(input_msg,*) tprpdum,realization_num
                  read(inpt,'(a80)') dummy_string
                  if (dummy_string(1:4) .eq. 'file') then
c Skip past file name
                     read(inpt,*)
                  else
                     stop
                  endif
                  read(iread_rcoll,'(a80)',end=19931) dummy_string
                  read(iread_rcoll,*,end=19931)idum
                  do while (idum.ne.realization_num)
                     nulldum=.false.
                     do while(.not.nulldum)
                        read(iread_rcoll,'(a80)',end=19931) dummy_string
                        nulldum=null1(dummy_string)
                     enddo
                     read(iread_rcoll,*,end=19931)idum
                  enddo
                  jjj=0
                  read(iread_rcoll,'(a80)',end=19931) dummy_string
                  nulldum=null1(dummy_string)
                  do while(.not.nulldum)
                     jjj=jjj+1
                     backspace(iread_rcoll)
                     read(iread_rcoll,*,end=19931)rcdiv(jjj,1),
     2                    probdiv(jjj,1)
                     read(iread_rcoll,'(a80)',end=19931) dummy_string
                     nulldum=null1(dummy_string)
                  enddo
                  nprobdivs(1) = jjj
                  k_rev=1.
                  r_min=rcdiv(1,1)
                  r_max=rcdiv(nprobdivs(1),1)
                  if(.not.allocated(rcoll_div))
     &               allocate(rcoll_div(num_part,1))
                  if(.not.allocated(ret_weight))
     &               allocate(ret_weight(num_part,1))
                   call impsample(tprpflag(jj))
c                  goto 19932
                   exit
19931             write(iptty,*)'Error in insptr for tprpflag=',
     1                 tprpflag(jj)
                  write(iptty,*)'realization_num ',realization_num,
     1                 ' not found in iread_coll file. STOP.'
                  write(ierr,*)'Error in insptr for tprpflag=',
     1                 tprpflag(jj)
                  write(ierr,*)'realization_num ',realization_num,
     1                 ' not found in iread_coll file. STOP.'
                  stop

c19932             continue

               elseif(tprpflag(jj).eq.13.or.tprpflag(jj).eq.14) then
                  div_flag=.true.
c s kelkar specifying log(CDF) vs log(Kf) as a st line equation
c with multiple realizations specified in the file iread_rcoll
                  read(inpt,'(a80)') dummy_string
                  if (dummy_string (1:4) .eq. 'file') then
c Skip past file name
                     read(inpt,*)
                     read(input_msg,*) tprpdum,realization_num
                     read(iread_rcoll,'(a80)',end=19941) dummy_string
                     read(iread_rcoll,*,end=19941)idum
                     do while (idum.ne.realization_num)
                        read(iread_rcoll,*,end=19941)idum
                     enddo
                     backspace(iread_rcoll)
                     read(iread_rcoll,*,end=19941)idum,k_rev,
     1                    r_min,r_max,slope_kf,cint_kf,
     2                    dispersivity1(jj),
     3                    dispersivity2(jj),dispersivity3(jj),
     4                    dispersivity4(jj),vratio(jj)
                     if(.not.allocated(rcoll_div))
     &                    allocate(rcoll_div(num_part,1))
                     if(.not.allocated(ret_weight))
     &                    allocate(ret_weight(num_part,1))
                     call impsample(tprpflag(jj))
                     goto 19942
19941                write(iptty,*)'Error in insptr for tprpflag=',
     1                    tprpflag(jj)
                     write(iptty,*)'realization_num ',realization_num,
     1                    ' not found in iread_coll file. STOP.'
                     write(ierr,*)'Error in insptr for tprpflag=',
     1                    tprpflag(jj)
                     write(ierr,*)'realization_num ',realization_num,
     1                    ' not found in iread_coll file. STOP.'
                     stop

19942                continue

                  else
c s kelkar specifying log(CDF) vs log(Kf) as a st line equation
                     backspace inpt
                     backspace inpt
                     read(inpt,*) tprpdum,k_rev,
     1                    r_min,r_max,slope_kf,cint_kf,
     2                    dispersivity1(jj),
     3                    dispersivity2(jj),dispersivity3(jj),
     4                    dispersivity4(jj),vratio(jj)
                     if(.not.allocated(rcoll_div))
     &                    allocate(rcoll_div(num_part,1))
                     if(.not.allocated(ret_weight))
     &                    allocate(ret_weight(num_part,1))
                     call impsample(tprpflag(jj))
                  end if
                  exit
               end if



               end do
               
               narrays=1
               itype(1)=4
               default(1)=1
               macro="sptr"
               igroup=4
               call initdata2(inpt,ischk,n0,narrays,itype,default,
     +              macroread(21),macro,igroup,ireturn,
     2              i4_1=itrc(1:n0))

            end if
         case('un','UN')
c     BAR 12-08-2000 Recognize keyword for unstructured particle
c     tracking run
! Currently does nothing
            done = .false.
         case('wt', 'WT')
            done = .false.
c keyword for specifying if initial locations are to be moved deltawt
c meters below the water table.
c default values set above in case wtdt keyword is not entered
!            smin=0.01
!            deltawt=0.001
            call parse_string(dummy_string, imsg, msg, xmsg, cmsg, nwds)
            if(nwds.gt.1) then
               if (msg(2).eq.1) then
                  deltawt=imsg(2)
               else
                  deltawt=xmsg(2)
               end if
               if (deltawt .eq. 0.) then
                  deltawt = .001
                  write(ierr,*) 'deltawt must be greater than 1.e-12, ',
     &                 'set to default value of 0.001'
               end if
               if(nwds.gt.2) then
                  if (msg(3).eq.1) then
                     smin=imsg(3)
                  else
                     smin=xmsg(3)
                  end if    
               endif
            else
               write(ierr,*)"check sptr input file. Can't use keyword"
               write(ierr,*)'wtdt without following it with a real #'
               write(ierr,*)' that equals deltawt.'
            endif
         case('vo','VO')
            done = .false.
c     keyword for recognizing if control volume values are output 
c     to .sptrx file
            sptrx_flag = .true.
            call parse_string(dummy_string, imsg, msg, xmsg, cmsg, nwds)
! cform default is formatted
            if (nwds .ge. 2) then
               if (msg(2) .eq. 3) then
                  if (cmsg(2)(1:4) .eq. 'unfo') cform(26)='unformatted'
               end if
            endif
c     Keywords to specify which parameters to output along the particle path
         case('xy', 'XY')
! Flag to indicate coordinates should be included for abbreviated output
            xyz_flag = .true.
            done = .false.
         case('ip', 'IP')
! Flag to indicate initial position should be included for abbreviated output
            ip_flag = .true.
            done = .false.
         case('tr', 'TR')
! Flag to indicate particle start time should be included for abbreviated output -- the initial node will be set to 0 if this option is used
            trans_flag = .true.
            done = .false.
! Flags to indicate whether particles that are outside the model domain should be included in the simulation
         case ('in', 'IN')
! Include all particles that are input (this is the default)
            exclude_particle = .false.
            done = .false.
         case ('ex', 'EX')
! Exclude particles with locations outside the model domain
            exclude_particle = .true.
            done = .false.
         case ('ou', 'OU')
! Exclude particles with locations outside the model domain
            output_end = .true.
            done = .false.
         case('po', 'PO')
            done = .false.
! If POD basis funtion derivatives are needed, zvd 16-Aug-04
c     ZVD 09-16-2004 Recognize keyword for output of derivatives for
c     model reduction POD basis functions
            if (dummy_string(3:3).eq.'d' .or. 
     &           dummy_string(3:3).eq.'D') then
! Do nothing, flag was set in scanin
            else
! Else po represents porosity
               flag_count = flag_count + 1  
               write_prop(1) = flag_count
            end if
         case('sa', 'SA')
! Saturation
            flag_count = flag_count + 1  
            write_prop(2) = flag_count
            done = .false.
         case('pe', 'PE')
! Permeability
            flag_count = flag_count + 1  
            write_prop(3) = flag_count
            done = .false.
         case('de', 'DE')
! Density
            flag_count = flag_count + 1  
            write_prop(4) = flag_count
            done = .false.
         case('pr', 'PR')
! Pressure
            flag_count = flag_count + 1  
            write_prop(5) = flag_count
            done = .false.
         case('te', 'TE')
! Temperature
            flag_count = flag_count + 1  
            write_prop(6) = flag_count
            done = .false.
         case('zo', 'ZO')
! Zone
            flag_count = flag_count + 1  
            write_prop(7) = flag_count
            done = .false.
         case('id', 'ID')
! Particle ID
            flag_count = flag_count + 1  
            write_prop(8) = flag_count
            done = .false.
         case('zb','ZB')
! ZVD 05/05 - added flagging alternate output, where
! a line will be written to the sptr3 file every time a particle
! enters a breakthrough zone
            done = .false.
            call parse_string(dummy_string, imsg, msg, xmsg, cmsg, nwds)
            if (nwds .gt. 1) then
               if (msg(2) .eq. 3 .and. cmsg(2) .eq. 'alt') then
                  alt_btc = .true.
                  if (nwds .gt. 2) then
                     if (msg(3) .eq. 3 .and. cmsg(3) .eq. 'xyz') then
                        xyz_flag2 = .true.
                     end if
                  end if
               else if (msg(2) .eq. 3 .and. (cmsg(2) .eq. 'tec' .or.
     &            cmsg(2) .eq. 'tecplot')) then
c     Output btc data with tecplot style header
                  techead = .true.
               end if
            end if
            
! ZVD 05/05 - added input for increasing btc time steps
            read(inpt,'(a80)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
c     assign nzbtc (ensure that either a real or integer
c     input is read properly)
            if (msg(1).eq.1) then
               nzbtc=imsg(1)
            else
               nzbtc=xmsg(1)
            endif
! ZVD 05/05 - if present read in minimum time step size, dtmn
            if(nwds.ge.2) then
               if (msg(2).eq.1) then
                  dtmn = imsg(2)
               else
                  dtmn = xmsg(2)
               endif
            end if
! ZVD 05/05 - if present read in time step multiplication factor
            if(nwds.ge.3) then
               if (msg(3).eq.1) then
                  part_mult = imsg(3)
               else
                  part_mult = xmsg(3)
               endif
            end if
! ZVD 05/05 - if present read in fraction of particles that should 
! break through before increasing the time step
            if(nwds.ge.4) then
               if (msg(4).eq.1) then
                  part_frac = imsg(4) * num_part
               else
                  part_frac = xmsg(4) * num_part
               endif
            end if
! ZVD 05/05 - if present read in fraction of particles that should 
! break through during a time step so that the time step is not 
! increased
            if(nwds.ge.5) then
               if (msg(5).eq.1) then
                  delta_part = imsg(5) * num_part
               else
                  delta_part = xmsg(5) * num_part
               endif            
            end if
! ZVD 05/05 - if present read in number of time steps that should be 
! checked for delta_part before increasing the time step
! increased
            if(nwds.ge.6) then
               if (msg(6).eq.1) then
                  part_steps = imsg(6)
               else
                  part_steps = xmsg(6)
               endif
            end if
! ZVD 05/05 - if present read in time (in days) small breakthrough 
! time steps should be started (for late initial breakthrough)
            if(nwds.ge.7) then
               if (msg(7).eq.1) then
                  time_btc = imsg(7)
               else
                  time_btc = xmsg(7)
               endif
            end if
            if(nzbtc.gt.0) then
               read(inpt,*)(zbtc(jj),jj=1,nzbtc)
            end if
         case('pl','PL')
c.... s kelkar 7/23/01.................
            done = .false.
            read(inpt,*) nplum
            if(nplum.gt.0) then
               read(inpt,*)(plum(jj),jj=1,nplum)
            end if
c......................................
         case('cl','CL')
c     S Kelkar 1-5-2006 Recognize keyword for cliff nodes option 
            done = .false.
            cliff_flag=.true.
            if (dummy_string(1:6) .eq. 'cliffg') then
! Read file names for generating & saving cliff node information
               read(inpt,'(a100)') clfg1name
               read(inpt,'(a100)') clfg2name
               read(inpt,'(a100)') clfilename
               icliff1 = open_file(clfg1name,'old')
               icliff2 = open_file(clfg2name,'old')
               icliff3 = open_file(clfilename,'unknown')
               call setup_cliffnodes(icliff1,icliff2,icliff3)
            else if (dummy_string(1:6) .eq. 'cliffr') then
! Read name of file containing cliff node information
               read(inpt,'(a100)') clfilename
               icliff1 = 0
               icliff2 = 0
               icliff3 = open_file(clfilename,'unknown')
               call setup_cliffnodes(icliff1,icliff2,icliff3)
            end if
         case('co','CO')
c.........
c.........s kelkar jan 31 06
c keyword for adjusting the velocities in corner nodes
            done = .false.
            epsilon_corner=0.01
            corner_flag=.true.
            call parse_string(dummy_string, imsg, msg, xmsg, cmsg, nwds)
            if(nwds.ge.2)  epsilon_corner= xmsg(2)         
c.........
         case('ca','CA')
c     S Kelkar 1-10-2005 Recognize keyword for pumping nodes for 
c     particle capture
c        read in data if used
            done = .false.
 9191       read(inpt,'(a80)') wdd1
            if(null1(wdd1)) then
               goto 9192
            elseif(wdd1(1:4).eq.'file') then
               read (inpt,'(a100)') tfilename
               ifile = open_file(tfilename,'old')
c     open(ifile,file=wdd1,form = 'formatted',
c     *           status = 'old')
               read(ifile,*) ja
c     assign array for storing well radii
               allocate(well_radius(ja))
c using irray(*,1) as scratch storage
               do jb=1,ja 
                  read(ifile,'(a80)') input_msg
                  call parse_string(input_msg, imsg, msg, xmsg, cmsg, 
     &                 nwds)
                  if (msg(1).eq.1) then
                     nnwell=imsg(1)
                  else
                     nnwell=xmsg(1)
                  endif
                  if(nwds.ge.2) then
                     if (msg(2).eq.1) then
                        rwell = imsg(2)
                     else
                        rwell = xmsg(2)
                     endif
                  else
c     set flag for default well radius=0.1*min(ddx,ddy)
                     rwell=1.1e+20
                  end if
                  irray(nnwell,0)=-(10000000+jb)
                  well_radius(jb)=rwell
               enddo
               close(ifile)
               goto 9192
            else
               backspace inpt
               read(inpt,*)ja,jb,jc
               if(ja.eq.0) goto 9192
c     Input by zones is first, or else input is by node
               if( ja .lt. 0 ) then
                  do inode = 1, n0
                     if( izonef(inode) .eq. abs(ja) ) then
                        irray(inode,0)=-(inode)
                     endif
                  enddo
               elseif(ja.le.n0.and.jb.le.n0) then
                  do inode=ja,jb,jc
                     irray(inode,0)=-(inode)
                  enddo
               else
                  write(ierr,*)'Error in keyword captur in sptr'
                  write(ierr,*)' bad integer values in sptr file'
                  write(ierr,*)'STOP'
                  stop
               endif
            endif
            goto 9191
 9192       continue
c.......................
         case('sp','SP')
c     S Kelkar 1-10-2005 Recognize keyword for spring nodes  
c     
c        read in data if used
            done = .false.
 9291       read(inpt,'(a80)') wdd1
            if(null1(wdd1)) then
               goto 9292
            elseif(wdd1(1:4).eq.'file') then
               read (inpt,'(a100)') tfilename
               ifile = open_file(tfilename,'old')
               read(ifile,*)ja
c using irray(*,1) as scratch storage
               read(ifile,*)(irray(jb,1),jb=1,ja)
               do jb=1,ja
                  inode=irray(jb,1)
                  irray(inode,0)=-(inode+2000000)
                  irray(jb,1)=0
               enddo
               close(ifile)
               goto 9292
            else
               backspace inpt
               read(inpt,*)ja,jb,jc
               if(ja.eq.0) goto 9292
c     Input by zones is first, or else input is by node
               if( ja .lt. 0 ) then
                  do inode = 1, n0
                     if( izonef(inode) .eq. abs(ja) ) then
                        irray(inode,0)=-(inode+2000000)
                     endif
                  enddo
               elseif(ja.le.n0.and.jb.le.n0) then
                  do inode=ja,jb,jc
                     irray(inode,0)=-(inode+2000000)
                  enddo
               else
                  write(ierr,*)'Error in keyword spring in sptr'
                  write(ierr,*)' bad integer values in sptr file'
                  write(ierr,*)'STOP'
                  stop
               endif
            endif
            goto 9291
 9292       continue

         end select read_flags
      end do

c ZVD 12-08-2005 if there are no flags set -- don't output anything 
c and issue warning so people will know the default has changed
c Old default was to output porosity and saturation

      if(flag_count.eq.0) then
         if (iptty .ne. 0 ) write (iptty, *) 
     &        'No parameter data will be output to sptr2 file'
      end if

      if(icnl.ne.0) then
         dispersivity3 = 0.
         dispersivity4 = 0.
      end if

c.......................

c      read(inpt,*) itm, ist
      read(inpt,'(a80)') dummy_string
      call parse_string(dummy_string, imsg, msg, xmsg, cmsg, nwds)
      count_steps_max = 1000000
      if (msg(1).eq.1) then
         itm=imsg(1)
      else
         itm=xmsg(1)
      end if
      if (msg(2).eq.1) then
         ist=imsg(2)
      else
         ist=xmsg(2)
      end if
      if (nwds .gt. 2) then
         if (msg(3) .eq. 1) then
            count_steps_max=imsg(3)
         else if (msg(3) .eq. 2) then
            count_steps_max=xmsg(3)
         else
            count_steps_max = 1000000
         end if
      end if   
! sptr_flag is set in scanin
      read(inpt,*) nx, ny, nz
      read(inpt,*) x10, y10, z10
      read(inpt,*) xdim, ydim, zdim
      part_id = 0
      ttp1 = 0.
      sptr_snode = .false.

      read(inpt,'(a80)') dummy_string
      if(.not. null1(dummy_string)) then
         if (dummy_string(1:4) .eq. 'file') then
            read (inpt,'(a100)') tfilename
            ifile = open_file(tfilename,'old')
            read(ifile,'(a150)') input_msg
            if (input_msg(1:4) .eq. 'FEHM') then
! This file was generated by FEHM and last seed and 
! ending time were written to file
               read(ifile,'(a80)') dummy_string
               read(ifile, *) rseed
! In case file was modified or manually generated
               if (rseed .eq. 0) rseed = 466201 
               sptr_snode = .true.
            else
               backspace ifile
            end if
         else
            backspace inpt
            ifile = inpt
         end if
         do i = 1,num_part
            read(ifile,'(a150)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            if (msg(1).eq.1) then
               ijkv(i)=imsg(1)
            else
               ijkv(i)=xmsg(1)
            end if
            if (msg(2).eq.1) then
               x1(i)=imsg(2)
            else
               x1(i)=xmsg(2)
            end if
            if (msg(3).eq.1) then
               y1(i)=imsg(3)
            else
               y1(i)=xmsg(3)
            end if
            if (msg(4).eq.1) then
               z1(i)=imsg(4)
            else
               z1(i)=xmsg(4)
            end if
            if (nwds .gt. 4 .and. msg(5) .ne. 3) then
               if (sptr_snode) then
! If a code generated file is being used for particle location input
! field 5 will be the particle number, field 6 will be the node, and
! field 7 will be particle travel time
                  if (msg(5).eq.1) then
                     part_id(i,1)=imsg(5)
                  else
                     part_id(i,1)=xmsg(5)
                  end if
                  part_id(i,2) = ijkv(i)
                  if (nwds .gt. 5 .and. abs(ist) .eq. 1) then
                     if (msg(6).eq.1) then
                        ijkv(i)=imsg(6)
                     else
                        ijkv(i)=xmsg(6)
                     end if
                     if (ist .gt. 0) then
                        if (msg(7).eq.1) then
                           ttp1(i)=imsg(7)
                        else
                           ttp1(i)=xmsg(7)
                        end if
                     else
                        ttp1(i) = 0.d0
                     end if
                  end if
               else
! If a manually generated file is being used for particle location input
! field 5 will be particle travel time
                  if (ist .gt. 0) then
                     if (msg(5).eq.1) then
                        ttp1(i)=imsg(5)
                     else
                        ttp1(i)=xmsg(5)
                     end if
                  else
                     ttp1(i) = 0.d0
                  end if
                  part_id(i,1) = i
                  part_id(i,2) = ijkv(i)
               end if
            else
               part_id(i,1) = i
               part_id(i,2) = ijkv(i)
            end if
         enddo
! If zone breakthrough was being done
         read(ifile,'(a80)') dummy_string
         if (dummy_string(1:4) .eq. 'zbtc') then
            call parse_string(dummy_string, imsg, msg, xmsg,
     &           cmsg, nwds)
            nzbtctmp = imsg(2)
            allocate (zbtctmp(nzbtctmp), tottmp(nzbtctmp))
            allocate (izonetmp(nzbtctmp, num_part))
            read(ifile, *) (zbtctmp(i), i = 1, nzbtctmp)
            read(ifile, *) (tottmp(i), i = 1, nzbtctmp)
            do i = 1, num_part
               read (ifile, *) (izonetmp(jj, i), jj = 1,nzbtctmp)
            end do 
            do i = 1, nzbtc
               do jj = 1, nzbtctmp
                  if (zbtc(i) .eq. zbtctmp(jj)) then
                     totalpart(i) = tottmp(jj)
                     izonebtc(i,1:num_part) = izonetmp(jj,1:num_part)
                     exit
                  end if
               end do
            end do
            deallocate (zbtctmp, tottmp, izonetmp) 
         else
            backspace ifile
         end if

         x3 = x1
         y3 = y1
         z3 = z1
! Read past blank line at end of sptr macro input
         read(ifile,'(a80)') dummy_string
         if(.not. null1(dummy_string)) write (ierr, 110)
         if (ifile .ne. inpt) close(ifile)
       endif
 110   format ('Error at end of sptr macro, input not terminated',
     &      ' with blank line, ist = ', i2)

! Open "isptr*" files if used, and write version identifier to file
      if (iprt .ne. 0) then
         open(isptr1, file = nmfil(17), status = cstats(17))
         write(isptr1 , 1000)  verno, jdate, jtime
         write(isptr1 , 1001)  ptitle
      end if
! Modifications by BAR for different output file formats for sptr2
      if (abs(iprto) .eq. 1) then
         open(isptr2, file = nmfil(18), status = cstats(18), 
     &        form = cform(18))
         write(isptr2 , 1000)  verno, jdate, jtime
         write(isptr2 , 1001)  ptitle
         if (iprto .lt. 0) write(isptr2, *) num_part
      else if (abs(iprto) .eq. 2) then
         open(isptr2, file = nmfil(18), status = cstats(18),
     2        form='unformatted')
         write(isptr2)  verno, jdate, jtime
         write(isptr2)  ptitle
         if (iprto .lt. 0) write(isptr2) num_part
      else if (abs(iprto) .eq. 3) then

! comment this binary code for linux/mac compiles with gfortran
         open(isptr2, file = nmfil(18), status = cstats(18),
     2        form='binary')
         write(isptr2)  verno, jdate, jtime
         write(isptr2)  ptitle         
         if (iprto .lt. 0) write(isptr2) num_part
! uncomment this Warning for linux/mac compile with gfortran
!        write(ierr, 120)
!        if (iout .ne. 0)  write(iout, 120)
!        if (iptty .ne. 0)  write(iptty, 120)
!        stop
 120     format ('Binary output not supported in this version ',
     &        '-- STOPPING (try unformatted)')

      end if
      if (nzbtc .gt. 0) then
         open(isptr3, file = nmfil(19), status = cstats(19))
         if (techead) then
            write(isptr3 , 2001)  trim(ptitle)
            formstring = 'VARIABLES = "Time (days)"'
            lf = 26
            do jj = 1, nzbtc
               write (formstring(lf:lf+12), 2002) zbtc(jj)
               lf = lf + 12
            end do
            write(isptr3 , '(a)') trim(formstring)
            write(isptr3 , 2000)  verno, jdate, jtime
         else
            write(isptr3 , 1000)  verno, jdate, jtime
            write(isptr3 , 1001)  ptitle
            if (alt_btc) then
               if (xyz_flag2) then
                  if (write_prop(3) .ne. 0) then
                     write(isptr3 , 1006)
                  else
                     write(isptr3 , 1005)
                  end if
               else
                  write(isptr3 , 1004)
               end if
            else
               write(isptr3 , 1002)
            end if
         end if
      end if
      if (sptr_flag) then
! File for saving particle locations at end of a simulation
! so that a restart can be made
         open(isptr8, file = nmfil(25), status = cstats(25))
         write(isptr8 , 1000)  verno, jdate, jtime
         write(isptr8 , 1001)  ptitle
      end if
      if(nplum.eq.-1) then
         open(isptr4, file = nmfil(20), status = cstats(20))
         write(isptr4 , 1000)  verno, jdate, jtime
         write(isptr4 , 1001)  ptitle
         write(isptr4, 1003) 
 
      endif
      if(nplum.eq.-2) then
         open(isptr5, file = nmfil(20), status = cstats(20))
         write(isptr5 , 1000)  verno, jdate, jtime
         write(isptr5 , 1001)  ptitle
      endif
      if (nplum .gt. 0) then
         open(isptr6, file = nmfil(20), status = cstats(20))
         write(isptr6 , 1000)  verno, jdate, jtime
         write(isptr6 , 1001)  ptitle
         write(isptr6 , 1002)
      end if
      if (sptrx_flag) then
! File for saving control volumes for plume_calc 
! s kelkar sep 12 05
!         ncoef=0
!         max_con=0
! Open file and write headers
         open(isptr9, file = nmfil(26), status = cstats(26),
     1           form=cform(26))
         if(cform(26).eq.'unformatted') then
            write(isptr9) verno, jdate, jtime
            write(isptr9) ptitle
! Moved to sptr_volume_out, nelm(neq_primary+1) not defined yet
!            write(isptr9) iw, neq_primary, nelm(neq_primary+1), 
!     &           ncoef, max_con
         elseif(cform(26).eq.'formatted') then
            write(isptr9, 1000) verno, jdate, jtime
            write(isptr9, 1001) ptitle
!            write(isptr9, '(5i10)' ) iw, neq_primary, 
!     &           nelm(neq_primary+1), ncoef, max_con
         endif
      end if

 1000 format(a30, 3x, a11, 3x, a8)
 2000 format('TEXT = "', a30, 3x, a11, 3x, a8, '"')
 1001 format(a80)
 2001 format('TITLE = "', a, '"')
 1002 format(' Time (days)      Zone1 Particles   . . .')
 2002 format(' "Zone ', i4, '"')
 1003 format(' Time (days)',4x, 'x(m)',4x,'y(m)',4x,'z(m)',4x,
     &     'node#',4x,'particle#',4x,'timestep')
 1004 format(2x, 'Time (days)', 7x, 'Particle#', 6x, 'ID', 
     &     4x, 'Zone', 4x, 'Node')
 1005 format(3x, 'x(m)', 13x, 'y(m)', 13x, 'z(m)', 13x, 'Time (days)', 
     &     7x, 'Particle#', 6x, 'ID', 7x, 'Zone', 4x, 'Node')
 1006 format(3x, 'x(m)', 13x, 'y(m)', 13x, 'z(m)', 13x, 'Time (days)', 
     &     7x, 'Particle#', 6x, 'ID', 4x, 'Zone', 4x, 'Node', 3x,
     &     'Log xperm')

      call plotc1(0,0)
      return

      contains
c.......................................................................
      subroutine read2_itensor1
c     Modified axisymmetric form of Lichtner et al.
c     dispersivity1 = molecular diffusion
c     dispersivity2,3,4,6  =  alpha1 thru 4
c     dispersivity7 and 8 = direction cosines in x and y dir for 
c     the axis of symmetry 

      use comai, only : inpt
      use compart
      use comsptr

      implicit none

      if(nwds.ge.10) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            dispersivity1(jj) =imsg(3)
         else
            dispersivity1(jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            dispersivity2(jj) =imsg(4)
         else
            dispersivity2(jj) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            dispersivity3(jj) =imsg(5)
         else
            dispersivity3(jj) =xmsg(5)
         end if
         if (msg(6).eq.1) then
            dispersivity4(jj) =imsg(6)
         else
            dispersivity4(jj) =xmsg(6)
         end if
         if (msg(7).eq.1) then
            dispersivity6(jj) =imsg(7)
         else
            dispersivity6(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity7(jj) =imsg(8)
         else
            dispersivity7(jj)  =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity8(jj) =imsg(9)
         else
            dispersivity8(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            vratio(jj) =imsg(10)
         else
            vratio(jj) =xmsg(10)
         end if                
      else
         call tprp_error ('read2_itensor1')
      endif

      if(nwds.ge.11) then
         call tprp_warn ('read2_itensor1', 11, 2)
         if (msg(11).eq.1) then
            divdwt(jj) =imsg(11)
         else
            divdwt(jj) =xmsg(11)
         end if
      endif

c     read(inpt,*) tprpflag(jj), kd(jj,1),
c     2                    dispersivity1(jj),dispersivity2(jj),
c     2                    dispersivity3(jj),dispersivity4(jj),
c     2                    dispersivity6(jj),dispersivity7(jj),
c     2                    dispersivity8(jj),vratio(jj)
      
      return
      
      end subroutine read2_itensor1

c...................................................................  
      subroutine read2_itensor2
c     Burnett and Frind tensor

      use comai, only : inpt
      use compart
      use comsptr

      implicit none

      if(nwds.ge.7) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            dispersivity1(jj) =imsg(3)
         else
            dispersivity1(jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            dispersivity2(jj) =imsg(4)
         else
            dispersivity2(jj) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            dispersivity3(jj) =imsg(5)
         else
            dispersivity3(jj) =xmsg(5)
         end if
         if (msg(6).eq.1) then
            dispersivity4(jj) =imsg(6)
         else
            dispersivity4(jj) =xmsg(6)
         end if
         if (msg(7).eq.1) then
            vratio(jj) =imsg(7)
         else
            vratio(jj) =xmsg(7)
         end if                
      else
         call tprp_error ('read2_itensor2')
      endif

      if(nwds.ge.8) then
         call tprp_warn ('read2_itensor2', 8, 2)
         if (msg(8).eq.1) then
            divdwt(jj) =imsg(8)
         else
            divdwt(jj) =xmsg(8)
         end if
      endif
      
c     read(inpt,*) tprpflag(jj), kd(jj,1),
c     2                       dispersivity1(jj),
c     2                       dispersivity2(jj),dispersivity3(jj),
c     2                       dispersivity4(jj),vratio(jj)
      
      return
      
      end subroutine read2_itensor2

c...................................................................
      subroutine read2_itensor3
c     modified Burnett and Frind derived by Lichtner et al.

      use comai, only : inpt
      use compart
      use comsptr

      implicit none

      if(nwds.ge.8) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            dispersivity1(jj) =imsg(3)
         else
            dispersivity1(jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            dispersivity2(jj) =imsg(4)
         else
            dispersivity2(jj) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            dispersivity3(jj) =imsg(5)
         else
            dispersivity3(jj) =xmsg(5)
         end if
         if (msg(6).eq.1) then
            dispersivity4(jj) =imsg(6)
         else
            dispersivity4(jj) =xmsg(6)
         end if
         if (msg(7).eq.1) then
            dispersivity6(jj) =imsg(7)
         else
            dispersivity6(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            vratio(jj) =imsg(8)
         else
            vratio(jj) =xmsg(8)
         end if                
      else
         call tprp_error ('read2_itensor3')
      endif

      if(nwds.ge.9) then
         call tprp_warn ('read2_itensor3', 9, 2)
         if (msg(9).eq.1) then
            divdwt(jj) =imsg(9)
         else
            divdwt(jj) =xmsg(9)
         end if
      endif
      
c     read(inpt,*) tprpflag(jj), kd(jj,1),
c     2                       dispersivity1(jj),
c     2                       dispersivity2(jj),dispersivity3(jj),
c     2                       dispersivity4(jj),dispersivity6(jj),
c     3                       vratio(jj)
      
      return
      
      end subroutine read2_itensor3

c...................................................................
      subroutine read2_itensor4
c     Isotrpoic tensor - Tompson et al.

      use comai, only : inpt
      use compart
      use comsptr

      implicit none

      if(nwds.ge.6) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            dispersivity1(jj) =imsg(3)
         else
            dispersivity1(jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            dispersivity2(jj) =imsg(4)
         else
            dispersivity2(jj) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            dispersivity3(jj) =imsg(5)
         else
            dispersivity3(jj) =xmsg(5)
         end if
         if (msg(6).eq.1) then
            vratio(jj) =imsg(6)
         else
            vratio(jj) =xmsg(6)
         end if                
      else
         call tprp_error ('read2_itensor4')
      endif

      if(nwds.ge.7) then
         call tprp_warn ('read2_itensor4', 7, 2)
         if (msg(7).eq.1) then
            divdwt(jj) =imsg(7)
         else
            divdwt(jj) =xmsg(7)
         end if
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     dispersivity1(jj),
c     2     dispersivity2(jj),dispersivity3(jj),
c     3     vratio(jj)
      
      return
      
      end subroutine read2_itensor4

c...................................................................
      subroutine read2_itensor5
c     Old version of fehm

      use comai, only : inpt
      use compart
      use comsptr

      implicit none

      if(nwds.ge.7) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            dispersivity1(jj) =imsg(3)
         else
            dispersivity1(jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            dispersivity2(jj) =imsg(4)
         else
            dispersivity2(jj) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            dispersivity3(jj) =imsg(5)
         else
            dispersivity3(jj) =xmsg(5)
         end if
         if (msg(6).eq.1) then
            dispersivity4(jj) =imsg(6)
         else
            dispersivity4(jj) =xmsg(6)
         end if
         if (msg(7).eq.1) then
            vratio(jj) =imsg(7)
         else
            vratio(jj) =xmsg(7)
         end if                
      else
         call tprp_error ('read2_itensor5')
      endif

      if(nwds.ge.8) then
         call tprp_warn ('read2_itensor5', 8, 2)
         if (msg(8).eq.1) then
            divdwt(jj) =imsg(8)
         else
            divdwt(jj) =xmsg(8)
         end if
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     dispersivity1(jj),
c     2     dispersivity2(jj),dispersivity3(jj),
c     2     dispersivity4(jj),vratio(jj)
      
      return
      
      end subroutine read2_itensor5

c...................................................................
      subroutine read4_itensor1
c     Modified axisymmetric form of Lichtner et al.
c     dispersivity1 = molecular diffusion
c     dispersivity2,3,4,6  =  alpha1 thru 4
c     dispersivity7 and 8 = direction cosines in x and y dir for 
c     the axis of symmetry 

      use comai, only : inpt
      use comdi
      use compart
      use comsptr

      implicit none

      if(nwds.ge.14) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            diffmfl(1,jj) =imsg(3)
         else
            diffmfl(1,jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            rd_frac(jj,1) =imsg(4)
         else
            rd_frac(jj,1) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            matrix_por(jj)= imsg(5)
         else
            matrix_por(jj)= xmsg(5)
         end if
         if (msg(6).eq.1) then
            aperture(jj)= imsg(6)
         else
            aperture(jj)= xmsg(6)
         end if

         if (msg(7).eq.1) then
            dispersivity1(jj) =imsg(7)
         else
            dispersivity1(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity2(jj) =imsg(8)
         else
            dispersivity2(jj) =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity3(jj) =imsg(9)
         else
            dispersivity3(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            dispersivity4(jj) =imsg(10)
         else
            dispersivity4(jj) =xmsg(10)
         end if
         if (msg(11).eq.1) then
            dispersivity6(jj) =imsg(11)
         else
            dispersivity6(jj) =xmsg(11)
         end if
         if (msg(12).eq.1) then
            dispersivity7(jj) =imsg(12)
         else
            dispersivity7(jj)  =xmsg(12)
         end if
         if (msg(13).eq.1) then
            dispersivity8(jj) =imsg(13)
         else
            dispersivity8(jj) =xmsg(13)
         end if
         if (msg(14).eq.1) then
            vratio(jj) =imsg(14)
         else
            vratio(jj) =xmsg(14)
         end if
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         if (nwds .ge. 15 .and. msg(15) .ne. 3) then
            call tprp_warn ('read4_itensor1', 15, 1)
            if (msg(15).eq.1) then
               secondary(jj) =imsg(15)
            else
               secondary(jj) =xmsg(15)
            end if
         end if
         if(nwds .ge. 16 .and. msg(16) .ne. 3) then
            call tprp_warn ('read4_itensor1', 16, 2)
            if (msg(16).eq.1) then
               divdwt(jj) =imsg(16)
            else
               divdwt(jj) =xmsg(16)
            end if
         endif
      else
         call tprp_error ('read4_itensor1')
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     diffmfl(1,jj), rd_frac(jj,1),
c     3     matrix_por(jj),aperture(jj),
c     4     dispersivity1(jj),dispersivity2(jj),
c     5     dispersivity3(jj),dispersivity4(jj),
c     6     dispersivity6(jj),dispersivity7(jj),
c     7     dispersivity8(jj),vratio(jj)
      
      return
      
      end subroutine read4_itensor1

c...................................................................
      subroutine read4_itensor2
c     Burnett and Frind tensor

      use comai, only : inpt
      use comdi
      use compart
      use comsptr

      implicit none

      if(nwds.ge.11) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            diffmfl(1,jj) =imsg(3)
         else
            diffmfl(1,jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            rd_frac(jj,1) =imsg(4)
         else
            rd_frac(jj,1) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            matrix_por(jj)= imsg(5)
         else
            matrix_por(jj)= xmsg(5)
         end if
         if (msg(6).eq.1) then
            aperture(jj)= imsg(6)
         else
            aperture(jj)= xmsg(6)
         end if
         if (msg(7).eq.1) then
            dispersivity1(jj) =imsg(7)
         else
            dispersivity1(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity2(jj) =imsg(8)
         else
            dispersivity2(jj) =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity3(jj) =imsg(9)
         else
            dispersivity3(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            dispersivity4(jj) =imsg(10)
         else
            dispersivity4(jj) =xmsg(10)
         end if
         if (msg(11).eq.1) then
            vratio(jj) =imsg(11)
         else
            vratio(jj) =xmsg(11)
         end if                
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         if (nwds .ge. 12 .and. msg(12) .ne. 3) then
            call tprp_warn ('read4_itensor2', 12, 1)
            if (msg(12).eq.1) then
               secondary(jj) =imsg(12)
            else
               secondary(jj) =xmsg(12)
            end if
         end if
         if(nwds .ge. 13. and. msg(13) .ne. 3) then
            call tprp_warn ('read4_itensor2', 13, 2)
            if (msg(13).eq.1) then
               divdwt(jj) =imsg(13)
            else
               divdwt(jj) =xmsg(13)
            end if
         endif
      else
         call tprp_error ('read4_itensor2')
      endif

c     read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     diffmfl(1,jj), rd_frac(jj,1),
c     3     matrix_por(jj),aperture(jj),
c     4     dispersivity1(jj),dispersivity2(jj),
c     5     dispersivity3(jj),dispersivity4(jj),
c     6     vratio(jj)
      
      return
      
      end subroutine read4_itensor2

c...................................................................
      subroutine read4_itensor3
c     modified Burnett and Frind derived by Lichtner et al.

      use comai, only : inpt
      use comdi
      use compart
      use comsptr

      implicit none

      if(nwds.ge.12) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            diffmfl(1,jj) =imsg(3)
         else
            diffmfl(1,jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            rd_frac(jj,1) =imsg(4)
         else
            rd_frac(jj,1) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            matrix_por(jj)= imsg(5)
         else
            matrix_por(jj)= xmsg(5)
         end if
         if (msg(6).eq.1) then
            aperture(jj)= imsg(6)
         else
            aperture(jj)= xmsg(6)
         end if

         if (msg(7).eq.1) then
            dispersivity1(jj) =imsg(7)
         else
            dispersivity1(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity2(jj) =imsg(8)
         else
            dispersivity2(jj) =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity3(jj) =imsg(9)
         else
            dispersivity3(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            dispersivity4(jj) =imsg(10)
         else
            dispersivity4(jj) =xmsg(10)
         end if
         if (msg(11).eq.1) then
            dispersivity6(jj) =imsg(11)
         else
            dispersivity6(jj) =xmsg(11)
         end if
         if (msg(12).eq.1) then
            vratio(jj) =imsg(12)
         else
            vratio(jj) =xmsg(12)
         end if                
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         if (nwds .ge. 13 .and. msg(13) .ne. 3) then
            call tprp_warn ('read4_itensor3', 13, 1)
            if (msg(13).eq.1) then
               secondary(jj) =imsg(13)
            else
               secondary(jj) =xmsg(13)
            end if
         end if
         if(nwds .ge. 14 .and. msg(14) .ne. 3) then
            call tprp_warn ('read4_itensor3', 14, 2)
            if (msg(14).eq.1) then
               divdwt(jj) =imsg(14)
            else
               divdwt(jj) =xmsg(14)
            end if
         endif      
      else
         call tprp_error ('read4_itensor3')
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     diffmfl(1,jj), rd_frac(jj,1),
c     3     matrix_por(jj),aperture(jj),
c     4     dispersivity1(jj),dispersivity2(jj),
c     5     dispersivity3(jj),dispersivity4(jj),
c     6     dispersivity6(jj),vratio(jj)

      return
      
      end subroutine read4_itensor3

c...................................................................
      subroutine read4_itensor4
c     Isotrpoic tensor - Tompson et al.

      use comai, only : inpt
      use comdi
      use compart
      use comsptr

      implicit none

      if(nwds.ge.10) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            diffmfl(1,jj) =imsg(3)
         else
            diffmfl(1,jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            rd_frac(jj,1) =imsg(4)
         else
            rd_frac(jj,1) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            matrix_por(jj)= imsg(5)
         else
            matrix_por(jj)= xmsg(5)
         end if
         if (msg(6).eq.1) then
            aperture(jj)= imsg(6)
         else
            aperture(jj)= xmsg(6)
         end if
         if (msg(7).eq.1) then
            dispersivity1(jj) =imsg(7)
         else
            dispersivity1(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity2(jj) =imsg(8)
         else
            dispersivity2(jj) =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity3(jj) =imsg(9)
         else
            dispersivity3(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            vratio(jj) =imsg(10)
         else
            vratio(jj) =xmsg(10)
         end if                
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         if (nwds .ge. 11 .and. msg(11) .ne. 3) then
            call tprp_warn ('read4_itensor4', 11, 1)
            if (msg(11).eq.1) then
               secondary(jj) =imsg(11)
            else
               secondary(jj) =xmsg(11)
            end if
         end if
         if(nwds.ge.12.and. msg(12) .ne. 3) then
            call tprp_warn ('read4_itensor4', 12, 2)
            if (msg(12).eq.1) then
               divdwt(jj) =imsg(12)
            else
               divdwt(jj) =xmsg(12)
            end if
         endif
      else
         call tprp_error ('read4_itensor4')
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     diffmfl(1,jj), rd_frac(jj,1),
c     3     matrix_por(jj),aperture(jj),
c     4     dispersivity1(jj),dispersivity2(jj),
c     5     dispersivity3(jj),vratio(jj)

      return
      
      end subroutine read4_itensor4

c...................................................................
      subroutine read4_itensor5
c     Old version of fehm

      use comai, only : inpt
      use comdi
      use compart
      use comsptr

      implicit none

      if(nwds.ge.11) then
         tprpflag(jj)=imsg(1)
         if (msg(2).eq.1) then
            kd(jj,1)=imsg(2)
         else
            kd(jj,1)=xmsg(2)
         end if
         if (msg(3).eq.1) then
            diffmfl(1,jj) =imsg(3)
         else
            diffmfl(1,jj) =xmsg(3)
         end if
         if (msg(4).eq.1) then
            rd_frac(jj,1) =imsg(4)
         else
            rd_frac(jj,1) =xmsg(4)
         end if
         if (msg(5).eq.1) then
            matrix_por(jj)= imsg(5)
         else
            matrix_por(jj)= xmsg(5)
         end if
         if (msg(6).eq.1) then
            aperture(jj)= imsg(6)
         else
            aperture(jj)= xmsg(6)
         end if

         if (msg(7).eq.1) then
            dispersivity1(jj) =imsg(7)
         else
            dispersivity1(jj) =xmsg(7)
         end if
         if (msg(8).eq.1) then
            dispersivity2(jj) =imsg(8)
         else
            dispersivity2(jj) =xmsg(8)
         end if
         if (msg(9).eq.1) then
            dispersivity3(jj) =imsg(9)
         else
            dispersivity3(jj) =xmsg(9)
         end if
         if (msg(10).eq.1) then
            dispersivity4(jj) =imsg(10)
         else
            dispersivity4(jj) =xmsg(10)
         end if
         if (msg(11).eq.1) then
            vratio(jj) =imsg(11)
         else
            vratio(jj) =xmsg(11)
         end if                
c zvd 05/21/07  Added parameter for optional input of secondary spacing
         if (nwds .ge. 12 .and. msg(12) .ne. 3) then
            call tprp_warn ('read4_itensor5', 12, 1)
            if (msg(12).eq.1) then
               secondary(jj) =imsg(12)
            else
               secondary(jj) =xmsg(12)
            end if
         end if
         if(nwds .ge. 13 .and. msg(13) .ne. 3) then
            call tprp_warn ('read4_itensor5', 13, 2)
            if (msg(13).eq.1) then
               divdwt(jj) =imsg(13)
            else
               divdwt(jj) =xmsg(13)
            end if
         endif
      else
         call tprp_error ('read4_itensor5')
      endif

c      read(inpt,*) tprpflag(jj), kd(jj,1),
c     2     diffmfl(1,jj), rd_frac(jj,1),
c     3     matrix_por(jj),aperture(jj),
c     4     dispersivity1(jj), dispersivity2(jj),
c     5    dispersivity3(jj),dispersivity4(jj),
c     6     vratio(jj)
      
      return
      
      end subroutine read4_itensor5

c...................................................................
      subroutine tprp_error (subname)

      use comai, only : ierr, iptty 
      use comsptr, only : tprpflag

      implicit none
      character*14 subname

      write(ierr, 100) subname
      write(ierr, 200)
      write(ierr, 300) jj, tprpflag(jj)
      if (iptty .ne. 0) then
         write(iptty, 100) subname
         write(iptty, 200)
         write(iptty, 300) jj, tprpflag(jj)
      end if

 100  format (a14, ': error in sptr macro tprp input')
 200  format ('# of entries do not match tprpflag: STOPPING')
 300  format ('model = ', i4, ' tprpflag = ', i1)

      end subroutine  tprp_error  

c...................................................................
      subroutine tprp_warn (subname, ith, iopt)

      use comai, only : ierr, iptty 
      use comsptr, only : tprpflag

      implicit none
      integer ith, iopt
      character*14 subname

      write(ierr, 100) subname, ith
      write(ierr, 200)
      write(ierr, 300) jj, tprpflag(jj)
      if (iptty .ne. 0) then
         write(iptty, 100) subname, ith
         if (iopt .eq. 1) then
            write(iptty, 200)
         else if (iopt .eq. 2) then
            write(iptty, 250)
         end if
         write(iptty, 300) jj, tprpflag(jj)
      end if

 100  format (a14, ' CAUTION :', i2,'th entry on input line for tprp')
 200  format ('interpreted as secondary spacing')
 250  format ('interpreted as multiplying weight for Div.(D)')
 300  format ('model = ', i4, ' tprpflag = ', i1)

      end subroutine tprp_warn 
 
      end subroutine insptr
