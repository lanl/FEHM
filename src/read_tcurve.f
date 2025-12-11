      subroutine read_tcurve (tfilename, numpars, ierr)
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
!D1 Read type curves for parallel fracture dispersion interpolation 
!D1 calculations.
!D1 
!***********************************************************************
!D2 
!D2 REVISION HISTORY
!D2 
!D2 FEHM Version 2.20 [10086-STN-2.20-00]
!D2 Initial implementation: 28-Oct-02, Programmer: Z. Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/read_tcurve.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:48   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:13:30   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.3.5 Cell-based particle-tracking module
!D3 2.3.6 Streamline particle-tracking module
!D3 
!***********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Sandia data curves
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!***********************************************************************
 
      use compfrac
      use comai, only : neq
      use comdi, only : nspeci
      implicit none
      character*3 output_flag
      integer i, ierr, j, k, l, m, mm, n, dummy, numpars, nump_max
      integer ncurves
      real*8 dum(3)
      real*8 xconvert
      character*100 tfilename
      character*4 dummy_string
      character*80 input_msg
      integer msg(10)
      integer nwds
      real*8 xmsg(10)
      integer imsg(10),ditnumlines,i2,index
      character*32 cmsg(10)


      numparams = numpars
      open (121, file = tfilename, status = 'old', action = 'read')

      sigma_low = 1.d10
      sigma_high = 0.d0
      omega_low = 1.d10
      omega_high = 0.d0
      par3_low = 1.d10
      par3_high = 0.d0

      read(121,'(a4)') dummy_string
      if(dummy_string(1:4).eq.'free') then
         read(121,*) (log_flag(i),i=1,numparams)
         read(121,*) (normal_param(i),i=1,numparams)
         read(121,*) curve_structure

c     Input used to be 1 parameter, now it's 2. This section of code
c     checks for the number of parameters input. If only 1,
c     set the second parameter weight_factor to the old hardwired
c     value of 1.e-3 to make code behave as before. Otherwise,
c     set the value

         if(curve_structure.gt.1) then
            backspace 121
            read(121,'(a)') input_msg
            call parse_string(input_msg, imsg, msg, xmsg, cmsg, nwds)
            curve_structure = imsg(1)
            if(nwds.gt.1) then

c     Second number could be input as an integer, so make sure the
c     code accounts for that

               if(msg(2).eq.1) then
                  weight_factor = imsg(2)
               else
                  weight_factor = xmsg(2)
               end if
            else
               weight_factor = 1.e-3
            end if
         end if

         allocate(itf_curve(nspeci,neq,curve_structure))
         allocate(wt_curve(nspeci,neq,curve_structure))
         itf_curve = 0
         wt_curve = 0.
         read(121,*) ncurves
         allocate(param1(ncurves))
         allocate(param2(ncurves))
         param1 = 0.
         param2 = 0.
         nump1 = ncurves
         nump2 = 1
         if (numparams .gt. 2) then
            allocate (param3 (ncurves))
            param3 = 0.
            d4 = 4
         else
            d4 = 1
         end if
         nump3 = 1
         
      else
         backspace 121
         curve_structure = 0
         read (121, *) nump1
         allocate (param1 (nump1))
         read (121,*) (param1(i), i= 1, nump1) 
         do i = 1, nump1
            sigma_low(1) = min (sigma_low(1), param1(i))
            sigma_high(1) = max (sigma_high(1), param1(i))
         end do
         read (121, *) nump2
         allocate (param2 (nump2))
         read (121, *) (param2(i), i= 1, nump2)
         do i = 1, nump2
            omega_low(1) = min(omega_low(1), param2(i))
            omega_high(1) = max(omega_high(1), param2(i))
         end do
         if (numparams .gt. 2) then
            read (121, *) nump3
            allocate (param3 (nump3))
            read (121, *) (param3(i), i= 1, nump3)
            do i = 1, nump3
               par3_low(1) = min(par3_low(1), param3(i))
               par3_high(1) = max(par3_high(1), param3(i))
            end do
            d4 = 4
         else
            nump3 = 1
            allocate (param3 (nump3))
            param3 = 0.
            par3_low = 0.
            par3_high = 0.
            d4 = 1
         end if
      end if
      read (121, *) nump_max
      allocate (nump(nump1,nump2,nump3,d4))
      allocate (dtime(nump1,nump2,nump3,d4,nump_max), 
     .     conc(nump1,nump2,nump3,d4,nump_max))
      dtime = 0.
      conc = 0.
!      backspace (121)
      do l = 1, d4
         do i = 1, nump1
            do j = 1, nump2
               do k = 1, nump3
                  read (121, *) nump(i,j,k,l), (dum(n), 
     .                 n = 1, numparams)
                  if(curve_structure.gt.0) then
                     param1(i) = dum(1)
                     sigma_low(1) = min (sigma_low(1), dum(1))
                     sigma_high(1) = max (sigma_high(1), dum(2))
                     param2(i) = dum(2)
                     omega_low(1) = min (omega_low(1), dum(1))
                     omega_high(1) = max (omega_high(1), dum(2))
                     if(numparams.gt.2) then
                        param3(i) = max(1.d-10,dum(3))
                        par3_low(1) = min (par3_low(1), dum(3))
                        par3_high(1) = max (par3_high(1), dum(3))
                     else
                        par3_low = 0.
                        par3_high = 0.
                     end if
                  end if
                  do m = 1, nump(i,j,k,l)
                     read (121, *) dtime(i,j,k,l,m), conc(i,j,k,l,m)
                  enddo
                  if(numparams.gt.2) then
c                     xconvert = (param1(i)+param2(i))/
                     xconvert = param2(i)*(1.-param3(i))/
     2                    (param1(i)+param2(i)*param3(i))
                     dtime(i,j,k,l,1:nump(i,j,k,l)) = 
     2                    (dtime(i,j,k,l,1:nump(i,j,k,l))-1.)/xconvert
                  end if
                  mm=1
                  do m=2,nump(i,j,k,l)	
                     if(conc(i,j,k,l,m)-conc(i,j,k,l,mm)
     2                    .lt.-1.e-20 )then
                        write(ierr,*)'Decreasing data found in type ',
     .                       'curve, stop'
                        write(ierr,*)'point:',m,' conc=',
     .                       conc(i,j,k,l,m)
                        stop
                     end if
                     if(conc(i,j,k,l,m).ne.conc(i,j,k,l,mm))then
                        mm=mm+1
                        if(mm.ne.m)then
                           conc(i,j,k,l,mm)=conc(i,j,k,l,m)
                           dtime(i,j,k,l,mm)=dtime(i,j,k,l,m)
                        end if
                     end if
                  end do
                  nump(i,j,k,l)=mm
               enddo
            enddo
         enddo
      enddo
      sigma_low(2) = sigma_high(1)
      sigma_high(2) = sigma_low(1)
      omega_low(2) = omega_high(1)
      omega_high(2) = omega_low(1)
      par3_low(2) = par3_high(1)
      par3_high(2) = par3_low(1)
      ipartout = 0
      read(121,'(a3)',end=2000) output_flag
      if(output_flag(1:3).eq.'out') then
         ipartout = 1
         allocate(param_density(nump1,nump2,nump3))
         param_density = 0
      end if
 2000 continue
      close (121)

      pfrac_read = .true.
      return
      end
