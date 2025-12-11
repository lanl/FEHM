      subroutine inmptr(simnum)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine is actived with particle tracking.  It is called
CD1  from subroutine part_trak to read the data in the ptrk macro and
CD1  then call set_ptrk to initialize the position of the particles.
CD1
CD1  To reduce the use of memory and number of calculations under
CD1  multiple species simulation conditions. The ptrk macro was
CD1  read in at here instead of in input subroutine. Thus, we need
CD1  to reopen the input data file and find the ptrk macro.
CD1  added the following for multiple species,pout,and prnt_rst.
CD1  first read in number of species, then read in
CD1  rock properties: alpha_x,_y,_z,f_spacing, and M_porosity.
C***********************************************************************
CD2
CD2  REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-JAN-96    S. Henderson   22      Add prolog.
CD2              S. Henderson           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/inmptr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:09:02   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:10:18   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:24:50   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:03:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:42:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2 modified for rip-fehm runs using variable confactors.   2/6/98, cli 
CD2
CD2 modified for ingrowth model.     July 19, 1997,cli
CD2
CD2 Modified for multiple species    Apr 18, 1997, cli
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
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
c
c Date: Summer '94
c Author: Stephen Henderson, EES-5
c         Contact: Bruce Robinson, EES-4
c                  George Zyvoloski, EES-5

      use combi
      use comdi
      use compart
      use compfrac, only : ffmax
      use comsptr, only : tprpflag
      use comdti
      use comai
      use comki
      use comxi
      use comsk
      
      implicit none

      integer open_file, size_counter, jj, jjj, ipart, ith1, itemp
      integer:: transflag,ithp,inpt1, idpr, ninputs, nsizep
      integer:: cnum,ith0,ja,jb,jc,nloc,ldt,idaughter,itmp
      integer:: i, j, k,ith, layers, layeri, istp, numtot
      integer:: dump10,idummy,max1d2,isimnum,irfile,set_value
      integer:: failed_nodes,itmp_shrink,numparams,temp_aidex
      integer:: temp_aidex_1,total_1d_size,max1d, n_col_irv, ithn
      integer:: np_temp, lns, lns1, li, rseed_old
      integer::  np_temp_max = 0
      integer:: lnwds, nwds, ii   !cli
      integer:: icount, colloid_counter, rev_count, irrev_count

      real:: half_life,tmpcnsk,dt,statime,endtime,period
      real:: tmpli, omega, xmag, tmptime,dt1,dt9,dt10
      real:: pinmass,sumass,temp_confactor

      real*8 simnum, rln2,ptrparam(500), kf_avg	
      real*8 rect, sum_rcoll, kf_max, kf_min, h2o_diff, tort_diff
      real*8 dret, ret, r_avg, probret, retp, retn
      real*8 rlogmin, rlogmax, drlog, retlog, rlog_avg
      character*132 input_line
      character*10 dummy_string
      
      real r, ran_sp

      logical:: null1, used, ready, afm, done, dfree_flag, nulldum
      
      character(len=100):: pfilename
      character(len=100):: tfilename

!define parameters 

      integer, parameter::istp10=20, input_max=20,i4=4,i5=5,i9=9,i10=10
c bhl_12/8/05 
      integer, parameter::i11=11
c bhl_12/8/05         
c bhl_11/3/06 
      integer, parameter::i6=6
c bhl_11/3/06         
       
      real, parameter::   add0_1=0.1, xhalf=0.5, xquarter=0.25
      real, parameter::   pai2=6.2831852,max_p_limit=0.98
      real*8, parameter:: mintcl=1.E-30, p_init=-5.0000
      real*8 , allocatable :: ret_weight_temp(:)
      real*8 , allocatable :: rcoll_div_temp(:)

      integer imsg(input_max),msg(input_max)
      real*8  xmsg(input_max)
      character*32 cmsg(input_max)
      logical it_is_open

! initialize variables
      size_counter = 0
      cnum=0
      max1d=0
c bhl_11/3/06         
      min_part=1
c bhl_11/3/06         
      ipzone=0		
      rev_count = 0
      irrev_count = 0
      colloid_counter = 0
						
      ipzone1=1
      rln2=log(2.)
! Set nsp to 1 for mptr simulations
      nsp = 1

      if(idpdp.ne.0)then
         ipzone2=ipzone1+ipzone1
      else
         ipzone2=ipzone1
      endif
	
!find mptr macro in FEHM data file(inpt) for input ptrk data 
      inquire(unit=inpt,opened=it_is_open)
      if(.not.it_is_open) then
c	If file got closed, open it again
         open(inpt,file=nmfil(2))
      else
         rewind(inpt)
      end if
      do
         read(inpt,'(a4)')macro
         locate_zone_mptr: select case (macro)
c 02-Jan-12 zvd - add zonn (any type of zone definition needs to be re-read)
         case ('zone', 'zonn')
	    call start_macro(inpt,inpt1,macro)
	    cnum=cnum+1
	    call zone(cnum,inpt1)
	    if(inpt1 /= inpt)then	
               call done_macro(inpt1)
	    endif	
         case ('mptr')
	    call start_macro(inpt,inpt1,macro)
            exit
         end select locate_zone_mptr
      end do
!set up parameters for monitoring zones designated by zptr
	   
!cli, keep code compatible with V2.22, in case max1d is not read      
      msg=0
      imsg=0
      xmsg=0.
      lnwds=0
      read(inpt1,'(a132)') input_line
      call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)

      do ii=1,nwds
         if(msg(ii) == 3)then
            lnwds=lnwds+1
            exit
         end if
      end do

      nwds=nwds-lnwds
      if(nwds <=3 )then
         write(ierr,*)'Number of inputs is:',nwds
         write(ierr,*)'For MPTR macro input line'
         write(ierr,*)input_line
         write(ierr,*)'Program stopped due to input error'
         stop
      elseif(nwds == i4)then
         write(ierr,*)'Using FEHM V2.22, you may need to adjust nstep'
         write(ierr,*)'to run the problem successfully'
         write(ierr,*)''
      elseif(nwds == i5)then
         write(ierr,*)'Using FEHM V2.23 MPTR input format'
         write(ierr,*)''
         max1d=imsg(5)+xmsg(5)
c bhl_11/3/06 
      elseif(nwds == i6)then
         write(ierr,*)'Using FEHM V2.25 MPTR input format'
         write(ierr,*)''
         max1d=imsg(5)+xmsg(5)
         min_part=imsg(6)+xmsg(6)
c bhl_11/3/06 
      else
         write(ierr,*)'Number of inputs: ',nwds
         write(ierr,*)'Expecting either 4, 5 or 6 input parameters'
         write(ierr,*)'Please check the number of inputs for MPTR macro'
         write(ierr,*)'Make sure the input values are correct' 
         write(ierr,*)'Code continues with the above warnings'
         write(ierr,*)''
      end if

      nspeci=imsg(1)+xmsg(1)
      maxlayers=imsg(2)+xmsg(2)
      max_particles=imsg(3)+xmsg(3)
      ripfehm=imsg(4)+xmsg(4)
c bhl_11/3/06 
      if (min_part<1) min_part=1
c bhl_11/3/06 

! cli continue read the rest of the inputs
c prnt_rst also controls output to the .ptrk and .ptrk_fin files
c zvd 11/09/06 -- option to specify which parameters to output in .ptrk
c 1=Number Having Entered System,2=Number Currently In System,
c 3=Number Having Left System,4=Number Having Decayed,
c 5=Number Having Been Filtered,6=Number That Left This Time
      read (inpt1, '(a132)') input_line
      read (input_line, *) pout, prnt_rst
      if (abs(prnt_rst) .ge. 20) then
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
         if (nwds .gt. 2) then
            prnt_var = .false.
            do i = 3, nwds
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

! cli, check the values of max1d, to be compatible with V2.22
!max1d will be set to nstep if not given.

      if(max1d==0)then
         max1d=nspeci*nstep
         write(ierr,*)'------------------- Note ---------------------'
         write(ierr,*)'Using V2.22 input format. Assign 1-D array size'
         write(ierr,*)'max1d to nstep. Please adjust nstep to reflect'
         write(ierr,*)'maximum number of segments expected'
         write(ierr,*)'----------------------------------------------'
      else if(max1d <=1)then
         max1d = nspeci*nstep
         write(ierr,*)'------------------- Note ---------------------'
         write(ierr,*)'Using V2.23 input format'
         write(ierr,*)'Input value for max1d is:', max1d,
     &        ' which is small'
         write(ierr,*)'Adjust max1d to:',max1d
         write(ierr,*)'You may want to adjust max1d in the future to'
         write(ierr,*)'a suitable value'
         write(ierr,*)'----------------------------------------------'
      end if

c     ZVD 10-31-2002 Recognize keyword for type curve data input and 
c        read in data if used
      read(inpt1,'(a4)') macro
      if(macro(1:4).ne.'tcur') then
         backspace inpt1
      else
         read(inpt1,'(a132)') input_line
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
         numparams = imsg(1)
         if (numparams .eq. 3 .and. nwds .gt. 1) then
            ffmax = xmsg(2)
         else
!Use default value for maximum fracture flow in 
         end if
         read(inpt1,'(a100)') tfilename
         call read_tcurve (tfilename, numparams, ierr)
      end if

      read(inpt1,'(a4)')macro
      if(macro.eq.'zptr')then						
         read(inpt1,*)ipzone						
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
         if(.not.allocated(idzone)) then
            allocate (idzone(ipzone1),pcount(ipzone2,nspeci))
            allocate(idcmax(ipzone2,nspeci))
            allocate(idcavg(ipzone2,nspeci))
            allocate(idflux(ipzone2))
            pcount = 0
         end if
         read(inpt1,*)(idzone(idpr),idpr=1,ipzone)
      else
         ipzone=0
         ipzone1=1
         if(idpdp.ne.0)then
	    ipzone2=ipzone1+ipzone1
         else
	    ipzone2=ipzone1
         endif
         if(.not.allocated(pconcmax)) then
            allocate(pconcmax(1))
         end if
         pconcmax=0
         if(.not.allocated(idzone)) then
            allocate (idzone(ipzone1),pcount(ipzone2,nspeci))
            pcount = 0
         end if
         backspace(inpt1)
      endif

!read in particle tracking control variables and allocate memory

! zvd, keep code compatible with V2.22, in case rseed_release is not read      
      msg=0
      imsg=0
      xmsg=0.
      lnwds=0
      read(inpt1,'(a132)') input_line
      call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
c     zvd 02-May-08
c     Save rseed if GoldSim has passed in a value since  we need to use
c     it for particle size distribution
      if (ripfehm == 1) rseed_old = rseed
      rseed = imsg(1)+xmsg(1)
      if (nwds .eq. 2) then
         rseed_release=imsg(2)+xmsg(2)
      else
         rseed_release = rseed
      end if
! If this a GoldSIM simulation
      if (ripfehm == 1) then
! Set random seeds to zero since seeds will be set by GoldSIM
         rseed = rseed_old
         rseed_release = 0
      end if
      
c*** water table rise modification
      read(inpt1,'(a4)') macro
      if (macro.eq.'wtri') then
         read(inpt1,*) water_table
      else
         water_table = -1.d+10
         backspace (inpt1)
      end if
c*** water table rise modification

      read (inpt1, *) daycs,daycf,dayhf,dayhs
      if(.not.allocated(ptindex)) then
         if (nstep .lt. istp10) nstep = istp10
         max1d2=max1d+max1d
         allocate(ptindex(n0),insnode(neq),dum_p(min(10000,nstep)))
         allocate(lsport(max1d+nspeci),nsport(max1d+nspeci+1))
         allocate(nsegs(nspeci),nprevd(max1d))
         allocate(ivdt(max1d),tmsport(max1d2))
         allocate(confactor(max1d),gmol(nspeci))
         allocate(bconfactor(max1d))
         allocate(ioconfactor(nspeci),sumdecayed(nspeci))
         allocate(newconfactor(nspeci),aidex(nspeci))
         allocate(conftime(nspeci),astep(nspeci))
!cli added for fraction release in decay-ingrowth
         allocate(p_fraction(nspeci)) 
         sumdecayed=0
!set default p_fraction to 0.25
         p_fraction=xquarter
      else
c    Not the first iteration, so we have to temporarily put the array
c    sizes of a few arrays back to their original sizes
         deallocate(dum_p)
         deallocate(insnode)
         deallocate(ptindex)
         deallocate(pcnsk)
         allocate(dum_p(min(10000,nstep)))
         allocate(insnode(neq))
         allocate(ptindex(n0))
         allocate(pcnsk(n7))
      end if
!find an unused unit for read in randomly generated parameter file

      !initialize lsport,nprevd
      lsport=0
      nsport=0
      nprevd=0

      read(inpt1,'(a80)') wdd1
      if(wdd1(1:4).eq.'file') then
         read(inpt1,'(a100)') pfilename

         irfile = open_file(pfilename,'old')
         read(irfile,*) ninputs

	  
!Read through the file until we get to the one 
!we want for this realization 
        
         isimnum = int(simnum+0.01)				
         if(isimnum<=0)isimnum = 1				 
         do j = 1, isimnum						
            read(irfile,*)(ptrparam(i),i=1,ninputs)
         end do
         close(irfile)
      else
         backspace inpt1
      end if

      afm = .false.
      read(inpt1,'(a80)') wdd1
      if(wdd1(1:3).eq.'afm') then
         afm = .true.
      else
         backspace inpt1
      end if

c zvd 02/28/07 Add for free water diffusion coefficient and tortuosity
      dfree_flag = .false.
      read(inpt1,'(a80)') wdd1
      if (wdd1(1:4).eq.'dfre') then
         dfree_flag = .true.
      else
         backspace inpt1
      end if

!read in x,y,z dispersivities, f_aperture, and m_porosities
!if a value<0,use the randomly generated data in ptrparam
	    
      do j=1,maxlayers
         read(inpt1,'(a80)') wdd1
         if(null1(wdd1)) then
            exit
         else
	    backspace(inpt1)
            if(afm) then
               read(inpt1,*)tclx(1,j),tcly(1,j),tclz(1,j),
     &              aperture(j),matrix_por(j),
     3              sresidual(j),gamma_afm(j)
            else
               read(inpt1,*)tclx(1,j),tcly(1,j),tclz(1,j),
     &              aperture(j),matrix_por(j)
            end if

            if(tclx(1,j)<0.)then
               tclx(1,j)=ptrparam(int(add0_1+abs(tclx(1,j))))
            end if  
            if(tcly(1,j)<0.)then
               tcly(1,j)=ptrparam(int(add0_1+abs(tcly(1,j))))
            end if  
            if(tclz(1,j)<0.)then
               tclz(1,j)=ptrparam(int(add0_1+abs(tclz(1,j))))
            end if  
            if(aperture(j)<0.)then
               aperture(j)=ptrparam(int(add0_1+abs(aperture(j))))
            end if  
            if(matrix_por(j)<0.)then
               matrix_por(j)=ptrparam(int(add0_1+abs(matrix_por(j))))
            end if  
            if(afm) then
               if(sresidual(j)<0.)then
                  sresidual(j)=ptrparam(int(add0_1+abs(sresidual(j))))
               end if  
               if(gamma_afm(j)<0.)then
                  gamma_afm(j)=ptrparam(int(add0_1+abs(gamma_afm(j))))
               end if  
            end if
	    if(tclx(1,j)>=0..and.tclx(1,j)<mintcl)tclx(1,j)=mintcl
	    if(tcly(1,j)>=0..and.tcly(1,j)<mintcl)tcly(1,j)=mintcl
	    if(tclz(1,j)>=0..and.tclz(1,j)<mintcl)tclz(1,j)=mintcl
         endif
      end do
	  
      narrays=1
      itype(1)=4
      default(1)=1
      macro="mptr"
      igroup=5
      call initdata2(inpt1,ischk,n0,narrays,itype,default,
     &     macroread(21),macro,igroup,ireturn,i4_1=itrc(1:n0))


!read in ith species transport properties:trak_type,half_life
!daughter species, conversion factor(# particles/mole),
!new conversion factor, time to change conversion factor, and
!molecule weight.
  
      ori=0
      obj=0
      nloc=0
      itmp=1
      nreac=0 
      astep=0 
      total_1d_size=0
      insnode(itmp)=0

      do i=1, nspeci
         ithp=0
         sumass=0.
         msg=0
         imsg=0
         xmsg=0.
         lnwds=0
         read(inpt1,'(a132)')input_line
         call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)

         do ii=1,nwds  
            if(msg(ii) == 3)then
               lnwds=lnwds+1
               exit
	    end if
         end do

         nwds=nwds-lnwds
         if(nwds <i9 )then
            write(ierr,*)'Number of inputs is:',nwds
            write(ierr,*)'For MPTR macro for species: ',imsg(1)
            write(ierr,*)input_line
            write(ierr,*)'Program stopped for wrong number of inputs'
            stop
         elseif(nwds == i9)then
            ith=imsg(1)+xmsg(1)
            write(ierr,*)'Using FEHM V2.22 format for species:', ith
            write(ierr,*)'you may need to adjust nstep to run the ' 
            write(ierr,*)'problem successfully'
         elseif(nwds == i10)then
            ith=imsg(1)+xmsg(1)
            write(ierr,*)'Using FEHM V2.23 MPTR format for speci:',ith
            astep(ith)=imsg(i10)+xmsg(i10)
c     bhl_12/8/05
         elseif(nwds == i11)then
            ith=imsg(1)+xmsg(1)
            write(ierr,*)'Using FEHM V2.24 MPTR format for speci:',ith
            astep(ith)=imsg(i10)+xmsg(i10)
            cfraction(ith)=imsg(i11)+xmsg(i11)
c     bhl_12/8/05
         else
            ith=imsg(1)+xmsg(1)
            astep(ith)=imsg(i10)+xmsg(i10)
         end if
         
         trak_type(ith)=imsg(2)+xmsg(2)
         half_life=imsg(3)+xmsg(3)
         idaughter=imsg(4)+xmsg(4)
         temp_confactor=imsg(5)+xmsg(5)
         newconfactor(ith)=imsg(6)+xmsg(6)
         conftime(ith)=imsg(7)+xmsg(7)
         gmol(ith)=imsg(8)+xmsg(8)
         p_fraction(ith)=imsg(9)+xmsg(9)   

! cli to be compatible with V2.22 input structure, astep is assigned
! default values of nstep. If p_fraction>=1, p_fraction is set to 0.25
         
         if(astep(ith)<=1)then
            write(ierr,*)'Particle injection segments for species:',ith
            write(ierr,*)'input is:',astep(ith)
            astep(ith)=nstep
            write(ierr,*)'is reset to:',astep(ith)
            write(ierr,*)''
         end if

         if(p_fraction(ith)>1.0)then
	    write(ierr,*)'For species: ',ith
	    write(ierr,*)'Input p_fraction value is:',p_fraction(ith)
            p_fraction(ith)=xquarter
	    write(ierr,*)'is reset to: ',xquarter
	    write(ierr,*)''
         end if

         total_1d_size=total_1d_size+astep(ith)

         if(total_1d_size > max1d)then
	    write(ierr,*) 
	    write(ierr, 30) max1d, total_1d_size
	    write(ierr, 40)
            if (iptty .ne. 0) then
               write(iptty, 30) max1d, total_1d_size
               write(iptty, 40)
            end if
	    stop
         end if
 30      format ('Error, max array size max1d:', i8, 
     &        ' smaller than the sum of indi. species', i8) 
 40      format ('Program stopped')

         if(ith == 1)then
            aidex(ith)=1
         end if

         ith1=ith+1
         if(ith1 <= nspeci)then
            aidex(ith1)=aidex(ith)+astep(ith)
         end if
         
         confactor(aidex(ith))=temp_confactor
         
         kfact(ith)=0.
         nsegs(ith)=0
         if(half_life /= 0.)then
            kfact(ith)=rln2/half_life
            if(idaughter.ne.0)then
               nreac=nreac+1
               ori(nreac)=ith
               obj(nreac)=idaughter
            endif
         endif
         

c     Determine if particle size distribution is being used
c     If so, read in the distribution


         read(inpt1,'(a80)') wdd1
         if(wdd1(1:4).eq.'size') then
            size_counter = size_counter + 1
            sizes(ith) = size_counter
            jj = 0
 1991       continue
            read(inpt1,'(a80)') wdd1
            if (null1(wdd1)) then
               goto 1992
            else
               jj = jj + 1
               read(wdd1,*) porsize(jj,sizes(ith)),
     2              probsize(jj,sizes(ith))
               goto 1991
            end if
 1992       continue
            
            if(probsize(jj,sizes(ith)).ne.1.) then
               write(ierr, 100)
               write(ierr, 110)
               if (iout .ne. 0) then
                  write(iout, 100)
                  write(iout, 110)
               end if
               if(iptty.ne.0) then
                  write(iptty, 100)
                  write(iptty, 110)
               end if
               stop
            end if
 100        format ('Fatal error, the colloid particle size')
 110        format ('distribution must end at 1')

c     Assign particle size for each particle
            
            do ipart = 1, max_particles
               
               done = .false.
               j=1
               r=ran_sp(rseed)
               do while(probsize(j,sizes(ith)).lt.r)
                  j=j+1
               end do
c     Interpolate to get exact size
               partsize(ipart,sizes(ith))=porsize(j-1,sizes(ith))+
     2              (porsize(j,sizes(ith))-porsize(j-1,sizes(ith)))*
     3              (r-probsize(j-1,sizes(ith)))/
     4              (probsize(j,sizes(ith))-probsize(j-1,sizes(ith)))

            end do
         else
            backspace inpt1
         end if

C     Hari read in parameters for the colloid diversity model

         read(inpt1,'(a80)') wdd1
         if(wdd1(1:4).eq.'dive') then
            read(inpt1,*) flag_col_daughter(ith)
            if (flag_col_daughter(ith).eq.0) then
c     Daughter product is a colloid
               flag_diversity(ith) = .true.
               colloid_counter = colloid_counter + 1
               divs(ith) = colloid_counter
               jj = 1

               read(inpt1,'(a132)') input_line
               if (input_line(1:4) .eq. 'file') then
c     Hari read blank for filename since scanin read colloid 
c     distribution filename and opened file
                  read (inpt1,*)
               else
                  backspace inpt1
               end if
c     Use simnum that is passed in for the multiple realization
               isimnum = int(simnum+0.01)				
               if(isimnum<=0) isimnum = 1
               read(inpt1,'(a132)') input_line
               call parse_string(input_line,imsg,msg,xmsg,cmsg,nwds)
               tprpflag(ith) = imsg(1)
               if (ripfehm .eq. 0 .and. nwds .ge. 2) isimnum = imsg(2)

c     Hari choose whether to read in cdf or equation
               if (iread_rcoll .ne. 0) rewind(iread_rcoll)
               if(tprpflag(ith).eq.11.or.tprpflag(ith).eq.12)then 
c Reading a table from external file
                  read(iread_rcoll,'(a132)',end=19931) input_line
                  read(iread_rcoll,*,end=19931) idummy
                  do while (idummy.ne.isimnum)
                     nulldum=.false.
                     do while(.not.nulldum)
                        read(iread_rcoll,'(a132)',end=19931) input_line
                        nulldum=null1(input_line)
                     enddo
                     read(iread_rcoll,*,end=19931) idummy
                  enddo
                  jjj=0
                  read(iread_rcoll,'(a132)',end=19931) input_line
                  nulldum=null1(input_line)
                  do while(.not.nulldum)
                     jjj=jjj+1
                     backspace(iread_rcoll)
                     read(iread_rcoll,*,end=19931)rcdiv(jjj,divs(ith)),
     2                    probdiv(jjj,divs(ith)) 
                     read(iread_rcoll,'(a132)',end=19931) input_line
                     nulldum=null1(input_line)
                  enddo
                  nprobdivs(divs(ith)) = jjj
                  k_rev(divs(ith))=1.
                  r_min(divs(ith))=rcdiv(1,divs(ith))
                  r_max(divs(ith))=rcdiv(nprobdivs(divs(ith)),divs(ith))
c     Hari importance sampling will be done at the time of breakthrough
c     for each particle
c     call impsample(tprpflag(divs(ith)))
                  goto 19932

19931             write(iptty,*)'Error in inmptr for tprpflag=',
     1                 tprpflag(jj)
                  write(iptty,*)'realization_num ',isimnum,
     1                 ' not found in iread_coll file. STOP.'
                  write(ierr,*)'Error in inmptr for tprpflag=',
     1                 tprpflag(jj)
                  write(ierr,*)'realization_num ',isimnum,
     1                 ' not found in iread_coll file. STOP.'
                  stop
                  
19932             continue
                  
               elseif(tprpflag(ith).eq.13.or.tprpflag(ith).eq.14) then
c Reading an equation
                  if (iread_rcoll .eq. 0) then
c Read equations parameters from macro file
                     read(inpt1,*,end=19941)k_rev(divs(ith)),
     1                    r_min(divs(ith)),r_max(divs(ith)),
     2                    slope_kf(divs(ith)),cint_kf(divs(ith)) 
                     
                  else
c Read equation parameters from external file
                     read(iread_rcoll,'(a132)',end=19941) input_line
                     read(iread_rcoll,*,end=19941)idummy
                     do while (idummy.ne.isimnum)
                        read(iread_rcoll,*,end=19941)idummy
                     enddo
                     backspace(iread_rcoll)
                     read(iread_rcoll,*,end=19941)idummy,
     1                    k_rev(divs(ith)),r_min(divs(ith)),
     2                    r_max(divs(ith)),slope_kf(divs(ith)),
     3                    cint_kf(divs(ith)) 
                  end if
c     Hari ptrparam option
c     If k_rev is negative we are getting all parameters from ptrparams
                  if(k_rev(divs(ith))<0.)then
                     set_value=int(add0_1-k_rev(divs(ith)))
                     k_rev(divs(ith))=ptrparam(set_value)
                     set_value=int(add0_1-r_min(divs(ith)))
                     r_min(divs(ith))=ptrparam(set_value)
                     set_value=int(add0_1-r_max(divs(ith)))
                     r_max(divs(ith))=ptrparam(set_value)
                     set_value=int(add0_1-slope_kf(divs(ith)))
                     slope_kf(divs(ith))=ptrparam(set_value)
                     set_value=int(add0_1-cint_kf(divs(ith)))
                     cint_kf(divs(ith))=ptrparam(set_value)
                  end if  
                  
c     s kelkar importance sampling done in part_track
c     for each time step , each species
c     call impsample(tprpflag(jj))
                  goto 19942

19941             write(iptty,*)'Error in inmptr for tprpflag=',
     1                 tprpflag(jj)
                  write(iptty,*)'realization_num ',isimnum,
     1                 ' not found in iread_coll file. STOP.'
                  write(ierr,*)'Error in inmptr for tprpflag=',
     1                 tprpflag(jj)
                  write(ierr,*)'realization_num ',isimnum,
     1                 ' not found in iread_coll file. STOP.'
                  stop

19942             continue
               endif
            else
               ncoll_daugh=ncoll_daugh+1
               divs_d(ith) = ncoll_daugh
               goto 19943
            end if
            
            read(inpt1,'(a80)') wdd1
            call parse_string(wdd1,imsg,msg,xmsg,cmsg,nwds)
            if(wdd1(1:4).eq.'reve')then
c Hari calculate mean retardation factor for reversible colloid species
c s kelkar flag_log=0 gives linear average, =1 is log average 
               if (nwds .ge. 2 .and. msg(2) .eq. 1) 
     &              flag_log(ith) = imsg(2)
               rev_count = rev_count+1
               reves(divs(ith)) = rev_count
               
               if(tprpflag(ith).eq.11.or.tprpflag(ith).eq.12)then
                  sum_rcoll = 0.0
                  do jj = 1, nprobdivs(divs(ith))-1
                     if(flag_log(ith).eq.0) then
                        rect = (probdiv(jj+1,divs(ith))-
     2                       probdiv(jj,divs(ith)))*(
     3                       rcdiv(jj+1,divs(ith))+
     4                       rcdiv(jj,divs(ith)))/2.
                     else
c average of the log(R)
                        rect = (probdiv(jj+1,divs(ith))-
     2                       probdiv(jj,divs(ith)))*(
     3                       dlog10(rcdiv(jj+1,divs(ith)))+
     4                       dlog10(rcdiv(jj,divs(ith))))/2.
                     endif
                     sum_rcoll= sum_rcoll+rect
                  enddo
                  if(flag_log(ith).eq.0) then
                     mean_rcoll_reve(reves(divs(ith))) = sum_rcoll
                  else
                     mean_rcoll_reve(reves(divs(ith))) = 10.**sum_rcoll
                  endif
               elseif(tprpflag(ith).eq.13.or.tprpflag(ith).eq.14)then
c s kelkar do a quick numerical integration of Kf for the equation
c log10(CDF)=b+m(log10(Kf)), ie prob that R is in [R,R+dR] 
c = [m*(10.**b)*(K_rev**m)*({R-1}**{m-1})]*d(R)
                  kf_avg = 0.
                  if ((r_max(divs(ith)) - r_min(divs(ith))) .lt.
     &                 1.d-22) then
                     mean_rcoll_reve(reves(divs(ith))) = 
     &                    r_min(divs(ith))
                  else
                     if(flag_log(ith).eq.1) then
c     average of log(R) using the equation
                        rlog_avg=0.
                        rlogmin=dlog10(r_min(divs(ith)))
                        rlogmax=dlog10(r_max(divs(ith)))
                        drlog=(rlogmax-rlogmin)/100.
                        retp=r_min(divs(ith))
                        do ii=1,100
c     using the value of R at the mid-point
                           retlog=rlogmin+drlog*(ii-0.5)
                           retn = 10**(rlogmin+drlog*ii)
                           ret=10.**retlog
                           dret = retn - retp
                           probret = dret * slope_kf(divs(ith)) *
     1                          (10.**(cint_kf(divs(ith)))) *
     2                          (ret-1)**(slope_kf(divs(ith))-1.) *
     3                          (k_rev(divs(ith))**slope_kf(divs(ith)))
                           rlog_avg= rlog_avg+probret*retlog
                           retp = retn
                        enddo
                        r_avg = 10.**rlog_avg
                        mean_rcoll_reve(reves(divs(ith))) = r_avg
                     else
c     linear average using the equation
                        r_avg = 0.
                        dret=(r_max(divs(ith))-r_min(divs(ith)))/100.
                        do ii=1,100
c     using the value of R at the mid-point
                           ret=r_min(divs(ith))+dret*(ii-0.5)
                           probret= dret * slope_kf(divs(ith)) *
     1                          (10.**(cint_kf(divs(ith)))) *
     2                          (ret-1)**(slope_kf(divs(ith))-1.) *
     3                          (k_rev(divs(ith))**slope_kf(divs(ith)))
                           r_avg= r_avg+probret*ret
                        enddo
                        mean_rcoll_reve(reves(divs(ith))) = r_avg
                     endif
                  end if
c     calculate mean of equation
               endif
            else
c     assign sampling method (equal_weight, random, importance)
               if (nwds .ge. 2 .and. msg(2) .eq. 3) then
                  select case (cmsg(2)(1:2))
                  case ('eq','EQ')
                     flag_method(ith) = 1
                  case ('ra','RA')
                     flag_method(ith) = 0
                  case ('im','IM')
                     flag_method(ith) = 2
                  end select
               end if

               irrev_count = irrev_count+1
               irrevs(divs(ith)) = irrev_count
               
               flag_col_irrev(ith) = .true.
            endif

         else
            backspace(inpt1)
         endif
c     Hari end colloid diversity code
19943    continue
         if(i.eq.1)then
                                !initialize variables
            do j=1,maxlayers
               kd(j,ith)=0.
               rd_frac(j,ith)=1.
               diffmfl(ith,j)=mintcl
               dispflag(j,ith)=0
               diffflag(j,ith)=0
            enddo
         else
            do j=1,maxlayers
               kd(j,ith)=kd(j,ith0)
               rd_frac(j,ith)=rd_frac(j,ith0)
               diffmfl(ith,j)=diffmfl(ith0,j)
               dispflag(j,ith)=dispflag(j,ith0)
               diffflag(j,ith)=diffflag(j,ith0)
            enddo
         endif
         
c read in the # of layers we need to alternate,followed by 
c layer transport properties (transflag,kd,rd,mdif).
         
         read(inpt1,*)layers
         do j=1,layers
            read(inpt1,*) layeri, transflag
            backspace inpt1
            if(transflag.gt.0) then
               if (.not. dfree_flag) then
                  read(inpt1,*)layeri,transflag,kd(layeri,ith),
     &                 rd_frac(layeri,ith),diffmfl(ith,layeri)
               else
                  read(inpt1,*)layeri,transflag,kd(layeri,ith),
     &                 rd_frac(layeri,ith),h2o_diff,tort_diff
               end if
               kcoll(layeri,ith) = 0. 
               rcoll(layeri,ith) = 1.
               fcoll(layeri,ith) = 1. 
            else
               if (.not. dfree_flag) then
                  read(inpt1,*)layeri,transflag,kd(layeri,ith),
     &                 rd_frac(layeri,ith),diffmfl(ith,layeri),
     2                 kcoll(layeri,ith),rcoll(layeri,ith),
     3                 fcoll(layeri,ith)
               else
                  read(inpt1,*)layeri,transflag,kd(layeri,ith),
     &                 rd_frac(layeri,ith),h2o_diff,tort_diff,
     2                 kcoll(layeri,ith),rcoll(layeri,ith),
     3                 fcoll(layeri,ith)
               end if
               transflag = -transflag
            end if

c     if value<0, use randomly generated parameters in ptrparam

            if(kd(layeri,ith)<0.)then
               set_value=int(add0_1-kd(layeri,ith))
               kd(layeri,ith)=ptrparam(set_value)
            end if  
            if(kcoll(layeri,ith)<0.)then
               set_value=int(add0_1-kcoll(layeri,ith))
               kcoll(layeri,ith)=ptrparam(set_value)
            end if  
            if(fcoll(layeri,ith)<0.)then
               set_value=int(add0_1-fcoll(layeri,ith))
               fcoll(layeri,ith)=ptrparam(set_value)
            end if  
            if(rcoll(layeri,ith)<0.)then
               set_value=int(add0_1-rcoll(layeri,ith))
               rcoll(layeri,ith)=ptrparam(set_value)
            end if  
            if(rd_frac(layeri,ith)<0.)then
               set_value=int(add0_1-rd_frac(layeri,ith))
               rd_frac(layeri,ith)=ptrparam(set_value)
            end if  
            if (.not. dfree_flag) then
               if(diffmfl(ith,layeri)<0.)then
                  set_value=int(add0_1-diffmfl(ith,layeri))
                  diffmfl(ith,layeri)=ptrparam(set_value)
               end if  
            else
               if (h2o_diff .lt. 0. .and. tort_diff .lt. 0.) then
                  set_value=int(add0_1-h2o_diff)
                  diffmfl(ith,layeri) = ptrparam(set_value)
                  set_value=int(add0_1-tort_diff)
                  diffmfl(ith,layeri) = diffmfl(ith,layeri) * 
     &                 ptrparam(set_value)
               else if (h2o_diff .lt. 0.) then
                  set_value=int(add0_1-h2o_diff)
                  diffmfl(ith,layeri) = ptrparam(set_value) * tort_diff
               else if (tort_diff .lt. 0.) then
                  set_value=int(add0_1-tort_diff)
                  diffmfl(ith,layeri) = h2o_diff * ptrparam(set_value)
               else
                  diffmfl(ith,layeri) = h2o_diff * tort_diff
               end if
	    end if

c     set transport simulation type

            set_transport_type:select case(transflag)
            case(1)
               dispflag(layeri,ith)=0
               diffflag(layeri,ith)=0
            case(2)
               dispflag(layeri,ith)=1
               diffflag(layeri,ith)=0
            case(3)
               dispflag(layeri,ith)=0
               diffflag(layeri,ith)=1
            case(4)
               dispflag(layeri,ith)=1
               diffflag(layeri,ith)=1
            case(5)
               dispflag(layeri,ith)=0
               if(afm) then
                  diffflag(layeri,ith)=-2
               else
                  diffflag(layeri,ith)=-1
               end if
            case(6)
               dispflag(layeri,ith)=1
               if(afm) then
                  diffflag(layeri,ith)=-2
               else
                  diffflag(layeri,ith)=-1
               end if
            case(7)
               if(afm) then
                  diffflag(layeri,ith)=5
               else
                  diffflag(layeri,ith)=3
               end if
               dispflag(layeri,ith)=1
            case(8)
               if(afm) then
                  diffflag(layeri,ith)=-5
               else
                  diffflag(layeri,ith)=-7
               end if
               dispflag(layeri,ith)=1
            end select set_transport_type

         end do

c     read in source terms

         read(inpt1,*)ns
         do j=1,ns
            istp=0
            itmp=itmp+1
            read(inpt1,*)ja,jb,jc,tmpcnsk
            call loadpcnsk(ja,jb,jc,tmpcnsk,nloc)
            insnode(itmp)=nloc
            
            read(inpt1,'(a80)')wdd1
            do while(.not.null1(wdd1))
               backspace(inpt1)
               set_mass_release_pattern:select case(wdd1(1:3))
               case('sin')      !sin function release
                  read(inpt1,*)wdd1,statime,endtime,period,numtot,dt
                  if(dt.gt.0.)then
                     ldt=dt
                     dt=(endtime-statime)/ldt
                  else
                     ldt=(endtime-statime)/dt
                  endif
                  omega=pai2/period
                  xmag=numtot*omega/(1-cos(omega*(endtime-statime)))
                  t1sk(1)=statime
                  t2sk(1)=t1sk(1)+dt
                  tmptime=xhalf*dt+t1sk(1)
                  dum_p(1)=xmag*sin(omega*tmptime)
                  istp=istp+1
                  do k=2,ldt
                     t1sk(k)=t2sk(k-1)
                     t2sk(k)=t1sk(k)+dt
                     tmptime=tmptime+dt
                     dum_p(k)=xmag*sin(omega*tmptime)
                     if(dum_p(k) > 0)istp=istp+1
                  enddo
               case('exp')   !exponential release
                  read(inpt1,*)statime,endtime,period,numtot,dt
                  if(dt < 0.)then
                     dt=-dt
                     ldt=(endtime-statime)/dt
                  else
                     ldt=dt
                     dt=(endtime-statime)/ldt
                  endif
                  xmag=numtot/(1-exp(-period*(endtime-statime))/period)
                  t1sk(1)=statime
                  t2sk(1)=t1sk(1)+dt
                  tmptime=xhalf*dt+t1sk(1)
                  dum_p(1)=xmag*exp(-period*tmptime)
                  istp=istp+1
                  do k=2,ldt
                     tmptime=tmptime+dt
                     t1sk(k)=t2sk(k-1)
                     t2sk(k)=t1sk(k)+dt
                     dum_p(k)=xmag*exp(-period*tmptime)
                     istp=istp+1
                  enddo
               case default     !discrete release
                  istp=istp+1
                  read(inpt1,*,err=999)pinmass,t1sk(istp),t2sk(istp)
                  if(confactor(aidex(ith)) /= 0.)then
                     dum_p(istp)=pinmass*confactor(aidex(ith))+xhalf
                     sumass=pinmass+sumass
c     cli added to make sure we do not have cases with mass but no particle
                     if(pinmass /= 0..and.dum_p(istp) ==0)then
                        write(ierr,*)'In inmptr, mass input line:',istp
                        write(ierr,*)'confactor too small, 0 particle ',
     &                       'injected'
                     end if
                  else
                     dum_p(istp)=pinmass
                  endif
               end select set_mass_release_pattern
               read(inpt1,'(a80)')wdd1
            enddo
                  
	!set up parameters for daughter species
 999        continue
            if(idaughter /= 0)then
               if(istp == 1.and.dum_p(istp) > istp10)then
                  dt10=(t2sk(istp)-t1sk(istp))/real(istp10)
                  dump10=dum_p(istp)/istp10
                  dt1=dt10/dump10
                  dt9=dt10-dt1
                  do k=1,istp10
                     t1sk(k)=t1sk(1) + (k-1)*dt10
                     t2sk(k)=t1sk(k)+dt9
                     dum_p(k)=dump10
                  enddo
                  t2sk(istp10)=t2sk(istp10)+dt1
                  istp=istp10
               endif
            endif
            failed_nodes=insnode(itmp)-insnode(itmp-1)
           call set_mptr(itmp,failed_nodes,ith,istp,ithp,idaughter)
         end do
	
         if(i == 1)then
            itmp_shrink=itmp
         endif

         !set up conversion factors

         !converting confactor and bconfactor arrays from 2-D to 1D    	

         ith0=ith
         if(confactor(aidex(ith)) == 0.)then
            if(ripfehm == 1)then
               ioconfactor(ith)=1
            else
               confactor(aidex(ith))=1.
               bconfactor(aidex(ith))=1.
               itemp=aidex(ith)
               do k=2,astep(ith)
                  itemp=itemp+1
                  confactor(itemp)=confactor(aidex(ith))
                  bconfactor(itemp)=bconfactor(aidex(ith))
               enddo
            endif
         else    
            ioconfactor(ith)=0
            temp_aidex_1=aidex(ith)+astep(ith)-1
          !cli modified this part to include switch option for confactor
            confactor(temp_aidex_1)=confactor(aidex(ith))
            bconfactor(temp_aidex_1)=1.D0/abs(confactor(aidex(ith)))
            if(sumass /= 0.)confactor(aidex(ith))=ithp/sumass
            bconfactor(aidex(ith))=1.D0/confactor(aidex(ith))
            itemp=aidex(ith)
            do k=2,astep(ith)-1
               itemp=itemp+1
               confactor(itemp)=abs(confactor(temp_aidex_1))
               bconfactor(itemp)=bconfactor(temp_aidex_1)
            enddo
         endif
         npt(ith+1)=ith*n0
      end do
      npt(ith+1)=ith*n0
         
      npt(1)=0
c     zvd 02-May-08
c     For GoldSim reset the random seed to zero for start of simulation
c     It will be set in part_track from the value passed through in(4)
      if (ripfehm == 1) rseed = 0
         
      !check the size of max1d vs. the sum of astep to make sure
      !enough space is allocated for the 1-D arrays, lsport,nprevd,
      !confactor, bconfactor, ivdt, and tmsport.

      if(total_1d_size>max1d)then
         write(ierr,*)''
         write(ierr,*)'Error, maximum array size max1d:',max1d,
     &        'smaller' 
         write(ierr,*)'than the sum of individual species',
     &        total_1d_size
         write(iptty,*)'Error, maximum array size max1d:',max1d,
     &        'smaller' 
         write(iptty,*)'than the sum of individual species',
     &        total_1d_size
         write(ierr,*)'Program stopped'
         stop
      endif

      macroread(21) = .TRUE.
      restarting = .FALSE.

      !set the safe thresh hold for max # of particles
      if(ripfehm==1)max_particles=max_particles*max_p_limit
	
      !initialize variables, p(), and flow_ot()
      npn=0
      if(idpdp /= 0)then
         nsizep=(nelm(neq+1)-neq-1)+(nelm(neq+1)-neq-1)
      else
         nsizep=nelm(neq+1)-neq-1
      endif
      if(.not. allocated(p)) then
         allocate(p(nsizep))
      end if
      p = p_init		
      flow_ot=0.		

      if(pout<0) then
         write(istrc,*) 'Node Number, Time A Particle ',
     &        'Left This Node (sec)'
      else
         call plotc1(0,1)
      endif

      !finishing reading ptrk data and close the unit inpt1

      if(inpt1.ne.inpt)then						
         call done_macro(inpt1)					
      endif										
	
      !cli_08/26/1999:shrink insnode, ptindex, pcnsk, and dum_p arrays

      call shrink_array_ptindex_pcnsk
      
c....................................................................
c 3/15/07 s kelkar irreversible colloid diversity model
c 4/11/07 sampling with equal divisions in the CDF space
c we need to allocate ret_weight and rcoll_div as arrays
c first see if irreversible colloid species are present

               n_col_irv=0
               do ith=1,nspeci
                  if(flag_diversity(ith).and.flag_col_irrev(ith))then
                     n_col_irv=n_col_irv+1
                  endif
               enddo
               if(n_col_irv.gt.0.) then
               if (.not.allocated(ret_weight)) then
                 allocate(ret_weight(max_particles,n_col_irv))
                 allocate(rcoll_div(max_particles,n_col_irv))
               endif
               endif
               if(ncoll_daugh.gt.0.) then
               if (.not.allocated(ret_weight_daughter)) then
                 allocate(ret_weight_daughter(max_particles,
     1                                        ncoll_daugh))
                 ret_weight_daughter=1.0
               endif
               endif

               if (ripfehm.eq.0) then
                  do ith=1,nspeci
                     if(flag_diversity(ith).and.flag_col_irrev(ith))then
                        do li=1,nsegs(ith)
                           lns1=aidex(ith)+li
                           lns=lns1-1
                           np_temp=nsport(lns1)-lsport(lns)+1
                           if (np_temp .gt. 0) then
                              if( np_temp .gt. np_temp_max) np_temp_max=
     &                             np_temp
                              allocate(rcoll_div_temp(np_temp))
                              allocate(ret_weight_temp(np_temp))
                              rcoll_div_temp=0.
                              ret_weight_temp=0.
                              call impsample_ptrk(tprpflag(ith),ith,
     1                             np_temp,rcoll_div_temp,
     2                             ret_weight_temp,np_temp_max)
                              itemp=0
                              do i=lsport(lns),nsport(lns1)
                                 itemp=itemp+1
                                 ret_weight(i,irrevs(divs(ith)))=
     1                                ret_weight_temp(itemp)
                                 rcoll_div(i,irrevs(divs(ith)))=
     1                                rcoll_div_temp(itemp)
                              enddo
                              deallocate (ret_weight_temp)
                              deallocate (rcoll_div_temp)
                           end if
                        enddo   
                     endif
                  enddo
               endif
c end importance sampling changes...................................


      return

	
      contains
	 
!cli_added:The following subroutine is used to shrink the  
!array size of insnode, pcnsk(), ptindx(), and dum_p()
	
        subroutine shrink_array_ptindex_pcnsk

	implicit none

	integer:: i, error, nloc1
c	integer,allocatable, dimension(:):: tmp_ptindex,tmp_insnode
        integer, allocatable :: tmp_ptindex(:)
        integer, allocatable :: tmp_insnode(:)

	real(8), allocatable:: tmp_pcnsk(:)
	
!allocate temporary storage for ptindex and pcnsk

	nloc=insnode(itmp_shrink)
	nloc1=nloc+1
	if(.not.allocated(tmp_ptindex))then
	  allocate(tmp_ptindex(nloc1),tmp_pcnsk(nloc1))
	  allocate(tmp_insnode(itmp_shrink))
	endif
	
!store insnode in temporary storage

	do i=1,itmp_shrink
	  tmp_insnode(i)=insnode(i)
	enddo

!store ptindex and pcnsk values in temporory storage

	do i=1,nloc
	  tmp_ptindex(i)=ptindex(i)
	  tmp_pcnsk(i)=pcnsk(i)
	enddo

!deallocate  old insnode, ptindex, pcnsk, and dum_p 
	
	error=0
	deallocate(dum_p,insnode,ptindex, pcnsk, stat=error)
	if (error == 0) then
           if(iout .ne. 0) write(iout, 120)
           if(iptty.ne.0) write(iptty, 120)
	else
           if(iout .ne. 0) write(iout, 130)
           write(ierr, 130)
           if(iptty.ne.0) write(iptty, 130)
	endif
 120    format ('Arrays:dum_p,insnode,ptindex,& pcnsk deallocated')
 130    format ('Failed to deallocate:dum_p,insnode,ptindex,pcnsk')

!allocate new dum_p,insnode, ptindex, pcnsk.

	error=0
        allocate(ptindex(nloc1))
        allocate(pcnsk(nloc1))
        allocate(insnode(itmp_shrink))
        allocate(dum_p(1))
c	allocate(dum_p(1),insnode(itmp_shrink),ptindex(nloc1),
c     &         pcnsk(nloc1), stat=error)
	if(error == 0)then
           if (iout .ne. 0) then
              write(iout,*)''
              write(iout,*)'dum_p()   re_allocated with size:',
     .             size(dum_p)
              write(iout,*)'insnode() re_allocated with size:',
     .             size(insnode)
              write(iout,*)'ptindex() re_allocated with size:',
     .             size(ptindex)
              write(iout,*)'pcnsk()   re_allocated with size:',
     .             size(pcnsk)
           end if
           if(iptty.ne.0) then
              write(iptty,*)'dum_p()   re_allocated with size:',
     .             size(dum_p)
              write(iptty,*)'insnode() re_allocated with size:',
     .            size(insnode)
              write(iptty,*)'ptindex() re_allocated with size:',
     .            size(ptindex)
              write(iptty,*)'pcnsk()   re_allocated with size:',
     .            size(pcnsk)
	  endif
       else
          if (iout .ne. 0) then
             write(iout,*)''
             Write(iout,*)'Failed to re_allocate new array size for:'
             write(iout,*)'dum_p, insnode, ptindex, and pcnsk'
          end if
	  write(ierr,*)''
	  write(ierr,*)'Failed to re_allocate new array size for:'
	  write(ierr,*)'dum_p, insnode, ptindex, and pcnsk'
	  write(ierr,*)''
	  if(iptty.ne.0) then
             Write(iptty,*)'Failed to re_allocate new array size for:'
             write(iptty,*)'dum_p, insnode, ptindex, and pcnsk'
	  endif
	endif

!transfer values from the temporory storage 
!back to the new arrays

	do i=1,itmp_shrink
	  insnode(i)=tmp_insnode(i)
	enddo
	
	do i=1,nloc
	  ptindex(i)=tmp_ptindex(i)
	  pcnsk(i)=tmp_pcnsk(i)
	enddo

	
!deallocate temporary storage

	deallocate(tmp_ptindex,tmp_pcnsk,tmp_insnode)
	
	

	return
	end subroutine shrink_array_ptindex_pcnsk

      end subroutine inmptr
