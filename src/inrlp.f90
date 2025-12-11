      subroutine inrlp
!*************************************************************************
! Copyright  2015.   Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
! National  Laboratory  (LANL),  which is operated by  Los Alamos National
! Security, LLC  for the U. S. Department of Energy.  The U. S. Government
! has rights to use, reproduce, and distribute this software.  Neither the
! U. S. Government nor Los Alamos National Security, LLC or persons acting
! on their behalf,  make any warranty,  express or implied, or assumes any
! liability for the accuracy, completeness, or usefulness of the software,
! any information pertaining to the software,  or  represents that its use
! would not infringe privately owned rights.

! The  software  being licensed  may  be Export Controlled.  It may not be
! distributed  or  used by individuals  or entities prohibited from having
! access to the software package, pursuant to United States export control
! laws and regulations. An export control review and determination must be
! completed before LANS will provide access to the identified Software.
!*************************************************************************
!
      use comai
      use comco2, only : icarb
      use comcomp, only : ioil
      use comdi, only : irlp,xfperm,yfperm,zfperm
      use comdti, only : n0
      use comki
      use comrlp

      implicit none

      integer :: lastmodel = 0
      integer :: lastrun = 0
      integer :: lasttbl = 0
      integer i, it, j, j2, jt, k, k2, maxphase, maxcpl, ndx, nparams,mm
      integer table_unit, open_file, cn, tblnum, ip, ip2, couple,mi,ir,kf
      integer nwds, imsg(20), msg(20),ij,kl,water,ic,ictype(3),nrlp_phases
      integer flag
      integer, allocatable :: irlptmp(:)
      real*8, allocatable :: xfptmp(:), yfptmp(:), zfptmp(:)
      real*8 xmsg(20),su_cut,smcutf,slcut,hmin,alamdam,alpha,alpham,alamda
      real*8 scutm, smcut,smax,alambda,amladam,beta,fac,facf,smcutm,d1
      parameter(scutm = 1.d-03, maxphase = 30, maxcpl = 30)
      parameter(hmin = 1.d-8)
      parameter(su_cut = 0.99d00)
            
      character*32 cmsg(20), dum_string
      character*180 chdum
      character*200 table_file
      logical null1, null_new, cap

      save lastmodel, lastrun, lasttbl
      phase_comp(20,1)=3;phase_comp(20,2)=6  !water/co2_liquid
      phase_comp(21,1)=3;phase_comp(21,2)=4   !water/co2_gas
      phase_comp(22,1)=6;phase_comp(22,2)=4   !co2_liquid/co2_gas 
      phase_comp(23,1)=1;phase_comp(23,2)=7   !oil/water
      phase_comp(24,1)=1;phase_comp(24,2)=8    !gas/water
      phase_comp(25,1)=2;phase_comp(25,2)=1     !air/water
      phase_comp(26,1)=1;phase_comp(26,2)=6;phase_comp(26,3)=4    !water/co2_liquid/co2_gas
      phase_comp(28,1)=2;phase_comp(28,2)=5  !water/vapor
  	   nrlp_phases=0   ;it=0 ;k2=0 ! initialize table and rlp indices
      maxrp = 30
      maxcp = 30
      msg = 0
      imsg = 0
      xmsg = 0.
      cmsg = ''
      macro = 'rlpm'
      if (lastrun .ne. irun) then
         lastrun = irun
         lastmodel = 0
      end if

      allocate (irlptmp(n0))
      if(.not.allocated (rlp_phase)) then
         allocate (rlp_phase(nrlp, maxrp))
      endif
      if (.not. allocated (rlp_group)) then
         allocate (rlp_group(nrlp))
      endif
      if(.not.allocated(rlp_type2)) then
         allocate (rlp_type(nrlp))
         allocate (rlp_type2(nrlp, maxrp), rlp_pos(nrlp, maxphase))
         allocate (rlp_param(nrlp, maxrp * max_rp))
         allocate (rlp_fparam(nrlp, maxrp * max_rpf))
         rlp_pos = 0
         rlp_phase = 0
         rlp_type = 0
         rlp_type2 = 0
         rlp_param = 0.
         rlp_fparam = 0.
         allocate (cap_type(nrlp))
         allocate (cap_coupling(nrlp, maxcp), cap_type2(nrlp, maxcp))
         allocate (cap_param(nrlp, max_cp*maxcp)) 
         allocate (cap_fparam(nrlp, max_cp*maxcp))
         allocate (cap_pos(nrlp, maxcpl))
         cap_pos = 0
         cap_coupling = 0
         cap_type = 0
         cap_type2 = 0
         cap_param = 0.
         if (ntable .ne. 0) then
            allocate (tblindx(ntable, 2), rlp_table(ntblines, 5))
            tblindx = 0
            rlp_table = 0.
         end if
      end if

         if (.not. allocated(xfperm)) then
            allocate(xfperm(nrlp),yfperm(nrlp),zfperm(nrlp))
            xfperm = 1.d0
            yfperm = 1.d0
            zfperm = 1.d0
         else if (fperm_flag) then
!     Assign values from temporary fperm arrays
            allocate (xfptmp(nrlp), yfptmp(nrlp), zfptmp(nrlp))
            xfptmp = 1.d0
            yfptmp = 1.d0
            zfptmp = 1.d0
            do mi = 1, n0
               ir = irlp(mi)
               if (ir  .ne. 0) then
                  xfptmp(ir) = xfperm(mi)
                  yfptmp(ir) = yfperm(mi)
                  zfptmp(ir) = zfperm(mi)
               end if
            end do
            deallocate (xfperm, yfperm, zfperm)
            allocate(xfperm(nrlp),yfperm(nrlp),zfperm(nrlp))
            xfperm = xfptmp 
            yfperm = yfptmp
            zfperm =  zfptmp 
            deallocate (xfptmp, yfptmp, zfptmp)
         end if

     
      i = lastmodel
      j = 0
      it = lasttbl
      do
         read (inpt, '(a80)') chdum
         if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or.   &
              chdum(1:3) .eq. 'END') exit
 25         call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
!         write(*,*) 'line 93 ',chdum
         select case (cmsg(1))
         case ('group', 'GROUP')
         
!     i is the current working group number
!    first, do some checking on completeness of specifications in the previous group
            if(lastmodel>0) then
               if(j>1) then
!     call ch(water,i,j-1)
               else
                  if (j .eq. 1) then
                     if (rlp_phase(i,j) .lt. 9) then
                        write (ierr, 50) rlp_group(lastmodel)
!                        stop
                     end if
                  else
                     write (ierr, 50) rlp_group(lastmodel)
!                     stop
                  end if
               endif
            endif
!     increment group number (i)      		  
            i = i + 1            
            lastmodel = i
!     store the group number            
            rlp_group(i) = imsg(2)
            cycle
!     now we are inside the group, and the choices are rlp, table, or cap   
!*****************TABLE ***********************************************      
         case ('table', 'TABLE')  ! expecting 'table  table_index phase_couple'
            j=j+1   ! always increment rlp number
            k2=k2+1  ! always increment cap number
            if (nwds .lt. 2) write(ierr, 60) 1, cmsg(2), rlp_group(i)  ! expect at least phase_couple
			         nparams=4            
            tblnum = imsg(2) ! first param is table index  ; not sure this is used    
            if(nwds==6) then   ! older style input; ignore all parameters except last (phase couple)
             couple=ic(cmsg(6), ictype)
! gaz debug 120317            
!            couple = 25
			         else
             couple=ic(cmsg(2), ictype)  ! couple will be >20; ictype will include each phase            
            endif
!     increment table index            
            it = it + 1
            if (it .eq. 1) then
               tblindx(it,1) = 1
               ndx = 0
            else
!     keep track of how many lines are in the table            
               tblindx(it, 1) = tblindx(it - 1, 2) + 1
               ndx = tblindx(it - 1, 2)   ! last line of previous table
            end if
 
            rlp_phase(i,it)=couple 
            rlp_pos(i,couple)=it 
            rlp_pos(i, ictype(1)) = it  
            rlp_pos(i, ictype(2)) = it  
            
            rlp_type(i) = 11
            rlp_type2(i,couple) = 11
            rlp_type2(i, ictype(1)) = 11
            rlp_type2(i, ictype(2)) = 11
            k = (j - 1) * max_rp
       
            rlp_param(i, k + 1) = it
               cap_type2(i, couple) = 9  
               cap_pos(i, couple) = it
               k = (k2 - 1) * max_cp
!     don't think we need this one            		
               cap_param(i, k + 1) = it
            read (inpt, '(a80)') chdum
            if (chdum(1:4) .eq. 'file') then
!     Data is read from file
               read (inpt, '(a200)') table_file
               table_unit = open_file(table_file,'old')
!     Read past any header lines in the table (header lines should start with a character)
               do 
                  read (table_unit,'(a)') chdum(1:40)
                  call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
                  if (msg(1) .ne. 3) then
                     backspace (table_unit)
                     exit
                  end if
100             continue
!                write(*,'(a)') chdum(1:40)
               end do
            else
!     Data is found on the following lines
               table_unit = inpt
               backspace (inpt)
            end if
            do
               read (table_unit, '(a80)', end = 5) chdum
!     Input is terminated with a blank line or 'end' or end-of-file)
               if (null_new(chdum) .or. chdum(1:3) .eq. 'end') exit
               ndx = ndx + 1
               read (chdum, *) (rlp_table(ndx,cn), cn = 1, nparams)  ! read 4 parameters on each line
!     Capillary pressure is put into position 5
                rlp_table(ndx,5) = rlp_table(ndx,nparams)
            end do
 5          lasttbl = it
            tblindx(it, 2) = ndx
            if (table_unit .ne. inpt) close (table_unit)
            cmsg(2) = 'tabular'
!     End of case 'tabular'
         case ('rlp', 'RLP')
         if(cmsg(3).ne.'same') then
!  		 write(*,*) 'not same'
!     **************** REL PERMS ******************************
! increment rlp index           
            nrlp_phases = nrlp_phases+1
            j=nrlp_phases
            rlp_phase(i,j)=ic(cmsg(2), ictype)  
! for single phase: ic is an integer  < 10 corresponding to the phase, ictype=0).  
! for a phase couple: ic is an integer > 10 corresponding to the couple, ictype(1-3) has each phase
            rlp_pos(i, abs(rlp_phase(i, j))) = j
!     Read parameters associated with rlp phase

            k = (j - 1) * max_rp
            kf = (j -1)* max_rpf
            select case (cmsg(3))
            case ('constant','Cons')
               if (nwds .lt. 3) write(ierr, 60) 1, cmsg(2), rlp_group(i)
               rlp_type2(i, rlp_phase(i,j)) = 1
               rlp_param(i, k + 1) = xmsg(4) + imsg(4)
            case ('linear')
               if (nwds .lt. 4) write(ierr, 60) 2, cmsg(2), rlp_group(i)
               rlp_type2(i, rlp_phase(i,j)) = 2
               rlp_param(i, k + 1) = xmsg(4) + imsg(4)
               rlp_param(i, k + 2) = xmsg(5) + imsg(5)
            case ('exponential')
               if (nwds .lt. 6) write(ierr, 60) 3, cmsg(2), rlp_group(i)
               rlp_type2(i, rlp_phase(i,j)) = 3
               do mm=1,3
               rlp_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do
               mm=4
               if(nwds==7) then
                  rlp_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)      ! this is an optional fourth parameter         
               else
                  rlp_param(i, k + mm)= 1
               endif            
            case ('corey')
!     Restrict to two-phase models
               if (nwds .lt. 4) write(ierr, 60) 2, cmsg(2), rlp_group(i)
               smax=xmsg(5) + imsg(5)
               if(smax.le.0.1) then ! assume this is old-style input
               	write(ierr,'(a80,f5.2)') 'Assuming second parameter of rlp corey model should be 1.- ',smax
                smax=1.-smax
               endif              	
               rlp_type2(i, rlp_phase(i,j)) = 4
               rlp_param(i, k + 1) = xmsg(4) + imsg(4)
               rlp_param(i, k + 2) = smax
            case ('brooks-corey')
! ekeating edit 12/2012.  Now require 4 parameters:    srl srg lambda and 'l' (from Miller et al., pg 88)          
               if (nwds .lt. 5) write(ierr, 60) 3, cmsg(2), rlp_group(i)
               rlp_type2(i, rlp_phase(i,j)) = 5
               if (nwds.eq.5) then
                  write(ierr,*) nwds
                  write(ierr,'(a50)') 'warning:  using obsolete' 
                  write(ierr,'(a50)') ' brooks-corey input'
               endif
               do mm=1,3
               	rlp_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do
			   mm=4
			   if (nwds.eq.6) then 
			   	rlp_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
			   else
                  rlp_param(i, k + mm) =  2.0
               endif
            case('vg_cap','vg_cp','vg', 'vg_1','vg_1_cap', 'vg_ek')
            	rlp_fparam(i,kf+1)=-99;cap=.false.
!            	if(cmsg(3).ne.'vg') then
            		cap=.true.
!            	endif
            	if (nwds .lt. 6) write(ierr, 60) 3, cmsg(2), rlp_group(i)  ! require at least 3 parameters
            	if(cap) then
                		if(nwds.eq.7) then   ! alpha was specified, will be ignored
               	  			write(ierr,*) 'Alpha parameter in ', xmsg(6), i, k,' rlp macro will be ignored.'  &
               	  		  //' It should be specified in cap keyword'
							mm=6
							xmsg(mm)=xmsg(mm+1);imsg(mm)=imsg(mm+1)
			   			endif  
			   	endif 
               do mm=1,3
               	rlp_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do            

			   	if(cmsg(3).eq.'vg_cap'.or.cmsg(3).eq.'vg_cp') then
			   			rlp_type2(i, rlp_phase(i,j)) = 9
			   	elseif  (cmsg(3).eq.'vg') then 
			   			rlp_type2(i, rlp_phase(i,j)) = 7
			   	elseif  (cmsg(3).eq.'vg_1'.or.cmsg(3).eq.'vg_ek') then 	
			   			rlp_type2(i, rlp_phase(i,j)) = 6
			   	elseif  (cmsg(3).eq.'vg_1_cap') then 	
			   			rlp_type2(i, rlp_phase(i,j)) = 8
			   	elseif  (cmsg(3).eq.'vg_corey') then 	
			   			rlp_type2(i, rlp_phase(i,j)) = 10
			   	endif
			   						   	

! check to see if fracture model is specified            	
               read (inpt, '(a80)') chdum
               if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or.   &
                    chdum(1:3) .eq. 'END') exit
               call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
               if (cmsg(1) .eq. 'fracture') then
               if (nwds .lt. 8) write(ierr, 60) 7, cmsg(2), rlp_group(i)  ! make sure at least 7 parameters were specified
               if(nwds.eq.9) then   ! alpha was specified, will be ignored
				  if(cap) then
               	  write(ierr,*) 'Alpha parameter ', xmsg(4),imsg(4),' in fracture rlp macro will be ignored. '   &
               	  //' It should be specified in cap keyword'
						do mm=4,8
						xmsg(mm)=xmsg(mm+1);imsg(mm)=imsg(mm+1)
						end do
				  else
				      write(ierr,*) 'Too many parameters for fracture keyword after rlp keyword'
				      stop
				  endif
			   endif
			   do mm=1,7
			   	rlp_fparam(i,kf+mm)=xmsg(mm+1) + imsg(mm+1)			   
			   end do
               else
                  backspace inpt
               end if	
            case ('stone')
!     For 3-phase water/oil/gas, oil rel perm in terms of krow and krog
!     Table corresponding to oil_water, oil_gas phase and
!     connate water saturation
               rlp_type2(i, rlp_phase(i,j)) = 12
               rlp_param(i, k+1) = it - 1
               rlp_param(i, k+2) = it
               rlp_param(i, k+3) = xmsg(4) + imsg(4)
            case default
               write (ierr, 20) trim(cmsg(3)), rlp_group(i)
               stop
            end select
		 else
!		 	write(*,*) 'about to err'
		 	write(ierr,20) 'option "same" is obselete and is ignored'
		 endif
            
         case('cap','CAP')
!*************CAPILLARY PRESSURE *******************************************	
            k2 = k2 + 1				
            couple = ic(cmsg(2), ictype)
            cap_pos(i, couple) = k2
! Read parameters associated with capillary model
            k = (k2 - 1) * max_cp
            select case (cmsg(3))   ! number of parameters = nwds-3
            case ('linear')
               if (nwds .lt. 5) write(ierr, 60) 2, cmsg(3), rlp_group(i)
               cap_type2(i, couple) = 1
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
            case ('linear_for')
               if (nwds .lt. 5) write(ierr, 60) 2, cmsg(3), rlp_group(i)
               cap_type2(i, couple) = 2
               cap_param(i, k + 1) = xmsg(4) + imsg(4)
               cap_param(i, k + 2) = xmsg(5) + imsg(5)
               if (cap_param(i, k + 2) .gt. 0.) then
                  cap_param(i, k + 1) = cap_param(i, k + 1) /    &
                       cap_param(i, k + 2)
               else
                  cap_param(i, k + 1) = 0.
               end if
            case ('exponential')
               if (nwds .lt. 5) write(ierr, 60) 3, cmsg(3), rlp_group(i)
               cap_type2(i, couple) = 3
               do mm=1,3
               cap_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do
               mm=4
               if(nwds.eq.6) then
               cap_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               else
               cap_param(i, k + mm) = 1.
               endif
            case ('brooks-corey')
               if (nwds .lt. 9) write(ierr, 60) 6, cmsg(3), rlp_group(i)
!     ekeating 2012
!     cap_type(i, j) = 4
               cap_type(i) = 4
               cap_type2(i, couple) = 4
               do mm=1,6
               cap_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do 
            case ('vg', 'vg_cap', 'vg_cp', &
                   'vg_ek')
               cap_fparam(i,k+1)=-99    
               if (nwds .lt. 9) write(ierr, 60) 6, cmsg(3), rlp_group(i)
               if(cmsg(3).ne.'vg_ek') then
                  cap_type2(i, couple) = 5
               else
                  cap_type2(i, couple) = 6
                endif
               
               do mm=1,6
               cap_param(i, k + mm) = xmsg(mm+3) + imsg(mm+3)
               end do 
! Smin Smax, alpha (1/m),N,S1,S2
               alpha = cap_param(i, k + 3)
               beta = cap_param(i, k + 4)
               alamda = 1.0-1.0/beta
               smcut=(cap_param(i, k + 6)-cap_param(i, k + 1))/(cap_param(i, k + 2)-cap_param(i, k + 1))
               smcutm = max(smcut,scutm)
               slcut = smcutm*(cap_param(i, k + 2)-cap_param(i, k + 1)) + cap_param(i, k + 1)
               fac=cap_param(i, k + 5)
         		if (.not. allocated(vg1)) then
            		allocate (vg1(nrlp,22),vg2(nrlp,22),vg3(nrlp,22),vg4(nrlp,22))
            		allocate (cp1(nrlp,22),cp2(nrlp,22))
         		end if
         			flag=0
                  call vgcap_fit3(flag,cap_param(i, k + 1),cap_param(i, k + 2),slcut,smcutm,fac,  &
                      alpha,alamda,vg1(i,k2),vg2(i,k2),vg3(i,k2),vg4(i,k2)  &
                      ,hmin)
                     
! not sure if this needs to be set or not
               cap_param(i, k + 6) = smcutm
!     get fit at saturated end(star=su_cut)
               slcut = su_cut*(cap_param(i, k + 2)-cap_param(i, k + 1)) + cap_param(i, k + 1)
               flag=4;d1=0.
               call vgcap_fit3(flag,cap_param(i, k + 1),cap_param(i, k + 2),slcut,su_cut ,fac,alpha  &
                   ,alamda, d1,d1, cp1(i,k2),cp2(i,k2),hmin)
               read (inpt, '(a80)') chdum
               if (null1(chdum) .or. chdum(1:3) .eq. 'end' .or.    &
                chdum(1:3) .eq. 'END') exit
               call parse_string2(chdum,imsg,msg,xmsg,cmsg,nwds)
               
               if (cmsg(1) .eq. 'fracture') then
                  if (nwds .lt. 7) write(ierr, 60) 6, dum_string,  rlp_group(i)
                  do mm=1,6
                  cap_fparam(i, k + mm) = xmsg(mm+1) + imsg(mm+1)
                  end do
 
                  smcut = (cap_fparam(i, k + 6) - cap_fparam(i, k + 1))  &
                         / (cap_fparam(i, k + 2) -   cap_fparam(i, k + 1))
                   	cap_fparam(i, k + 6) = max (smcut, scutm)
				  facf=cap_fparam(i, k + 5)
                  alpham = alpha
                  alamdam = alamda
                  alpha = cap_fparam(i, k + 3)
                  beta = cap_fparam(i, k + 4)
                  alamda = 1.0-1.0/beta
                  smcut=(cap_fparam(i, k + 6)-cap_fparam(i, k + 1))/(cap_fparam(i, k + 2)-cap_fparam(i, k + 1))
                  smcutf = max(smcut,scutm)
                  slcut = smcutf*(cap_fparam(i, k + 2)-cap_fparam(i, k + 1)) + cap_fparam(i, k + 1)
                  if(fac.ge.0.0d00) then
! next call matches cutoff capillary pressure
                     facf = cap_fparam(i, k + 5)
                     call vgcap_match(2,smcutm,alpham,alamdam,fac,smcutf,alpha,alamda,facf)  ! this routine may change facf
                     fac=facf
                     cap_fparam(i, k + 5)=facf
                  else
                     fac = cap_fparam(i, k + 5)
                  endif

         		if (.not. allocated(vg1f)) then
            		allocate (vg1f(nrlp,22),vg2f(nrlp,22),vg3f(nrlp,22),vg4f(nrlp,22))
            		allocate (cp1f(nrlp,22),cp2f(nrlp,22))
         		end if
! need to check to make sure this call should be based on fac, not facf         		
                 call vgcap_fit3(0,cap_fparam(i, k + 1),cap_fparam(i, k + 2),slcut,smcutf,fac,   &
                       alpha,alamda,vg1f(i,k2),vg2f(i,k2),vg3f(i,k2),vg4f(i,k2)   &
                       ,hmin)
                  cap_fparam(i, k + 6) = smcutf
!     get fit at saturated end(star=su_cut)
                  slcut = su_cut*(cap_fparam(i, k + 2)-cap_fparam(i, k + 1)) + cap_fparam(i, k + 1)
!                 call vgcap_fit3(4,cap_param(i, k + 1),cap_param(i, k + 2),slcut,su_cut,fac   &
!                      ,alpha,alamda,0.0d0,0.0d0,cp1f(i,k2),cp2f(i,k2),hmin)
                flag=4;d1=0.
                  call vgcap_fit3(flag,cap_param(i, k + 1),cap_param(i, k + 2),slcut,su_cut,fac   &
                       ,alpha,alamda,d1,d1,cp1f(i,k2),cp2f(i,k2),hmin)
 
                 else
                  backspace inpt
                end if
        	    end select  !  finished looking at cap options
            case default
!  ************ if not group, table, rlp, or cap, assume it's 'rlp' for backwards compatability***************************
               write (ierr, '(3a15)') 'Assuming ',trim(cmsg(1)), 'is an rlp model'
               k=len(trim(chdum))
               chdum(6:k+5)=trim(chdum)
               chdum(1:5)='rlp '
               go to 25
!            end select
         end select
! take note of the fact if water has been defined            
!            if(rlp_phase(i,j)==2.or.rlp_phase(i,j)==3) water=j
	  
      end do

      narrays = 1
      itype(1) = 4
      default(1) = 0
      macro = "rlpm"
      igroup = 2
      call initdata2( inpt, ischk, n0, narrays,   &
           itype, default, macroread(7), macro, igroup, ireturn,   &
           i4_1=irlptmp(1:n0) )
!     write(49,*) 'inrlp :nrlp: ',nrlp
      do i = 1, n0
         do j = 1, nrlp
            if (rlp_group(j) .eq. irlptmp(i)) then
               irlp(i) = j
               exit
            end if
         end do
      end do
       
      macroread(7) = .TRUE.

      deallocate (irlptmp)	
 10   format ('Unrecognized phase ', a, ' in rlpm, check model ', i4) 
 20   format ('Unrecognized type ', a, ' in rlpm, check model ', i4)
 30   format ('Invalid model type ', a, ' in rlpm, check model ', i4)
 40   format ('Invalid coupling ', a, ' in rlpm, check model ', i4)
 50   format ('At least two phases must be specified in rlpm, ',  &
           'check model ', i4)
 60   format (i1, 'parameters are required for model type', a,  &
       ' in rlpm, check model ', i4)
 65   format ('At least ', i1, 'parameters are required for model type',  &
           a, ' in rlpm, check model ', i4)
      end subroutine inrlp
      subroutine ch(water,i,j)
      use comrlp, only : rlp_param,max_rp
      implicit none
      integer water,i,j,kl,max_rop,ij,k
! i is group number, j is number of phases
! now make sure smin and smax are consistent for all phases
!
	if(water==0) then
           write(*,*) 'no rlp model specified for water phase'
           stop
	endif
	kl=(water-1) * max_rp
	write(*,*) 'water = ',water,' group = ',i,   &
             ' number of phases = ',j
	do ij=1,j
           k = (ij - 1) * max_rp
           if(ij.ne.water) then
              rlp_param(i,k+1)=1.-rlp_param(i,kl+2)
              rlp_param(i,k+2)=1.-rlp_param(i,kl+1)   
           endif         
	end do    
	return
	end

                                                                                                                      