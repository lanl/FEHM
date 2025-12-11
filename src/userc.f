      subroutine userc(iz, i, rc_ss, drc_ss)
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
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To allow for user specified concentrations or fluxes as
CD1 a function of time.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 FEHM Version 2.0, SC-194
CD2
CD2 Initial implementation: 01-OCT-1993, Programmer: B. Robinson
CD2
CD2 $Log:   /pvcs.config/fehm90/src/userc.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:24   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:21:42   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:28:54   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:26   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:06   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:47:56 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2
C***********************************************************************
CD3 
CD3 REQUIREMENTS TRACEABILITY
CD3 
CD3 N/A
CD3 
C***********************************************************************
CD4 
CD4 SPECIAL COMMENTS AND REFERENCES
CD4 
CD4 This user subroutine has two options:
CD4    1 -- Use a given set of concentration v. time data to control
CD4         the injection concentration with time. To use this option
CD4         the time data must range from 0 seconds ago to x seconds 
CD4         ago, and the negative inlet concnetration must be the time
CD4         ago in years to start the simulation. 
CD4    2 -- Use a given set of tracer source sink (mol/s) data v. time
CD4         to set rc. Those nodes for which this option is to be used
CD4         must have an inlet concentration of -9876.
CD4
C***********************************************************************

      use combi
      use comdi
      use comdti
      use comai
      use comci
      use compart
      use comrxni, only : scl, dsccl
      use comuserc
      implicit  none

      integer jl,ju,jm,iz,i2
      real*8 slope,getconc,timein,getflux
      integer, optional :: i
      real(8), optional :: rc_ss, drc_ss
      integer ispecies, j, k, k1, curcolumn
      integer, allocatable :: nsindex(:)
      integer nsnodes, iread1, nszones, mi, mim
c     SPC
      integer flag_kd2, irip, jjj, sflag
      real*8  deltat, in3_old
      
      real*8 srmiml, rc_ss_erosion, drc_ss_erosion
      real*8 inflrate, molesin, molesout
c     SPC  add for transformation factor
c     real*8 transform_kd(8,10)
      real*8 transform_kd(100,100)
      integer flag_kd, i1,j1, ii, jj,  tf1
      
      logical used
      integer icount, open_file
      character*100 filename
c     SPC
      save flag_kd2,in3_old,tf1
      save transform_kd,flag_kd,  ii, jj
      
      
      if(ripfehm.EQ.0) then
         
         if ((iz.eq.0).and.(iuserc.ne.-9324)) then
C     Initialization
            iuserc=-9324
            if (iout .ne. 0) write (iout,*)
     &           'Solute transport user subroutine is invoked'
            if (iatty.gt.0) write(iatty,*)
     &           'Solute transport user subroutine is invoked'
            

            filename = ''
c     zvd 04/26/2007 Allow user to enter name of file for data
            read (inpt,'(a100)') filename
            if (filename(1:4) .eq. 'file') then
               read (inpt,'(a100)') filename
            else
               backspace inpt
               filename = ''
               filename(1:14) = 'userc_data.dat'
            end if
            iread1 = open_file(trim(filename),'old')

            read(iread1,*) usroption
            select case (usroption)
         case (1,2)
c     usroption = 2 Time-varying solute mass flux input at prescribed nodes
            read(iread1,*) nu_points
            if(.not.allocated(userconc)) then
               allocate(userconc(nspeci,nu_points), usertime(nu_points))
            end if
            read(iread1,*)
            do i2=1,nu_points
               read(iread1,*) (userconc(j,i2),j=1,nspeci)
            enddo
            read(iread1,*)
            do i2=1,nu_points
               read(iread1,*) usertime(i2)
            enddo
         case (3)
            read(iread1,*) nu_points, nsnodes
            if(.not.allocated(nodeindex)) then
               allocate(nodeindex(n0), nsindex(nsnodes))
            end if
            nodeindex = 0
            read(iread1,*)(nsindex(j),j=1,nsnodes)
            do j = 1, nsnodes
               nodeindex(nsindex(j)) = j
            end do
            if(.not.allocated(userconc3)) then
               allocate(userconc3(nspeci,nsnodes,nu_points),
     2              usertime(nu_points))
            end if
            do k = 1, nspeci
               do k1 = 1, nu_points
                  read(iread1,*) usertime(k1),
     2                 (userconc3(k,j,k1),j=1,nsnodes)
               end do
            end do
         case (4)
            read(iread1,*) nu_points, nszones
            if(.not.allocated(nodeindex)) then
               allocate(nodeindex(n0), nsindex(nszones))
            end if
            nodeindex = 0
            read(iread1,*)(nsindex(j),j=1,nszones)
            do j = 1, nszones
               do k = 1, n0
                  if(izonef(k).eq.nsindex(j)) then
                     nodeindex(k) = j
                  end if
               end do
            end do
            if(.not.allocated(userconc3)) then
               allocate(userconc3(nspeci,nszones,nu_points),
     2              usertime(nu_points))
            end if
            do k = 1, nspeci
               do k1 = 1, nu_points
                  read(iread1,*) usertime(k1),
     2                 (userconc3(k,j,k1),j=1,nszones)
               end do
            end do
         case(5)
c     usroption = 5 Solute recirculation model
c     read number of production/infiltration zones
            read (iread1, *) erosion_factor
            if (.not. allocated(recycle_factor)) then
               allocate (recycle_factor(nspeci))
            end if
            read (iread1, *) (recycle_factor(j), j = 1, nspeci)
            read (iread1, *) nu_zones
            if (.not. allocated(prodzones)) then
               allocate (prodzones(nu_zones) , inflzones(nu_zones))
               allocate (moles_recycle(nu_zones,nspeci))
            end if
            read (iread1, *) (prodzones(j), j = 1, nu_zones)
            read (iread1, *) (inflzones(j), j = 1, nu_zones)
c     convert erosion factor fromm 1/yr to 1/s
            erosion_factor = erosion_factor / (365.25 * 86400.)
         end select
         close(iread1)
      else if (iuserc.eq.-9324) then
C     transient changes
         if ((iz.eq.1).and.(usroption.eq.1)) then
            if (l.ne.0) then
               if (iuchk.ne.-9876) then 
                  do i2=1,n0
                     if (cnsk(i2).lt.0.) timeCl36start=-cnsk(i2)
                  enddo
                  iuchk=-9876
                  timeCl36start=timeCl36start*365.25*8.64e4
               endif
               timein=timeCl36start-days*8.64e4
               if (timein.le.0.) then 
                  getconc=userconc(1,1)
                  goto 1
               else if (timein.ge.usertime(nu_points)) then
                  i2=nu_points
                  goto 2
               endif
               
               jl=0
               ju=nu_points+1
 10            if (ju-jl.gt.1) then
                  jm=(ju+jl)/2
                  if(timein.gt.usertime(jm))then
                     jl=jm
                  else
                     ju=jm
                  endif
                  goto 10
               endif
               i2=jl+1
               
 2             if(i2.eq.1) then
                  getconc=userconc(1,1)
               else
                  slope=userconc(1,i2)-userconc(1,i2-1)
                  slope=slope/(usertime(i2)-usertime(i2-1))
                  getconc=userconc(1,i2-1)+slope*
     2                 (timein-usertime(i2-1))
               endif
               
 1             do i2=1,n0
                  if (cnsk(i2).lt.0.) cnsk(i2)=-getconc
               enddo
            endif
            
         else if ((iz.eq.2).and.(usroption.eq.2)) then
            if (l.ne.0) then
               if (iuchk.ne.-9876) then
                  ispecies = 1 + npn/n0
                  if(cnsk(i+npn) .eq. -9876) then
                     timein=days*8.64e4
                     if (timein.le.usertime(1)) then 
                        getflux=userconc(ispecies,1)
                        goto 13
                     else if (timein.ge.usertime(nu_points)) then
                        i2=nu_points
                        goto 23
                     endif
                     
                     jl=0
                     ju=nu_points+1
 103                 if (ju-jl.gt.1) then
                        jm=(ju+jl)/2
                        if(timein.gt.usertime(jm))then
                           jl=jm
                        else
                           ju=jm
                        endif
                        goto 103
                     endif
                     i2=jl+1
                     
 23                  if (i2.eq.1) then
                        getflux=userconc(ispecies,1)
                     else
                        slope=userconc(ispecies,i2)-
     2                       userconc(ispecies,i2-1)
                        slope=slope/(usertime(i2)-usertime(i2-1))
                        getflux=userconc(ispecies,i2-1)+slope*
     2                       (timein-usertime(i2-1))
                     endif
 13                  continue
                     rc_ss = -getflux
                     drc_ss = 0.
                  end if
               endif               
            endif

         else if ((iz.eq.2).and.(usroption.eq.3)) then
C---------------------------------------------------------
c-------phil added fix to bypass nodes not in the nodeindex
c     2/18/04
c---------------------------------------------------------
            if(nodeindex(i).NE.0) then
               if (l.ne.0) then
                  ispecies = 1 + npn/n0
                  curcolumn = nodeindex(i)
                  if(pcnsk(i+npn) .lt. 0.) then
                     timein=days*8.64e4
                     if (timein.le.usertime(1)) then 
                        getflux=userconc3(ispecies,curcolumn,1)
                        goto 43
                     else if (timein.ge.usertime(nu_points)) then
                        i2=nu_points
                        goto 33
                     endif
                     
                     jl=0
                     ju=nu_points+1
 303                 if (ju-jl.gt.1) then
                        jm=(ju+jl)/2
                        if(timein.gt.usertime(jm))then
                           jl=jm
                        else
                           ju=jm
                        endif
                        goto 303
                     endif
                     i2=jl+1
                     
 33                  if (i2.eq.1) then
                        getflux=userconc3(ispecies,curcolumn,1)
                     else
                        slope=userconc3(ispecies,curcolumn,i2)-
     2                       userconc3(ispecies,curcolumn,i2-1)
                        slope=slope/(usertime(i2)-usertime(i2-1))
                        getflux=userconc3(ispecies,curcolumn,i2-1)+
     2                       slope*(timein-usertime(i2-1))
                     endif
 43                  continue
                     cnsk(i+npn) = -getflux
                  end if
               end if               
            endif
c     (nodeindex.NE.0)  phs 2/18/04  
            
         else if ((iz.eq.2).and.(usroption.eq.4)) then
            if (l.ne.0) then
               ispecies = 1 + npn/n0
               curcolumn = nodeindex(i)
               if(pcnsk(i+npn) .lt. 0.) then
                  timein=days
                  if (timein.le.usertime(1)) then 
                     getflux=userconc3(ispecies,curcolumn,1)
                     goto 143
                  else if (timein.ge.usertime(nu_points)) then
                     i2=nu_points
                     goto 133
                  endif
                  
                  jl=0
                  ju=nu_points+1
 403              if (ju-jl.gt.1) then
                     jm=(ju+jl)/2
                     if(timein.gt.usertime(jm))then
                        jl=jm
                     else
                        ju=jm
                     endif
                     goto 403
                  endif
                  i2=jl+1
                  
 133              if (i2.eq.1) then
                     getflux=userconc3(ispecies,curcolumn,1)
                  else
                     slope=userconc3(ispecies,curcolumn,i2)-
     2                    userconc3(ispecies,curcolumn,i2-1)
                     slope=slope/(usertime(i2)-usertime(i2-1))
                     getflux=userconc3(ispecies,curcolumn,i2-1)+
     2                    slope*(timein-usertime(i2-1))
                  endif
 143              continue
c     
c     IF this is an isothermal air-water simualtion
c     
                  mi = i+npn
                  mim=mi-npn
                  if( ico2 .lt. 0 ) then
c     
c     Set liquid and vapor source/sink flow rates
c     
                     srmiml = sk(mim)
                     
c     ELSE
c     
                  elseif( ico2 .gt. 0 ) then
c     
c     Set liquid and vapor source/sink flow rates
c     
                     srmiml = sk(mim)
                     
                  else
                     srmiml = sk(mim) * s(mim)
                     
c     
c     ENDIF
c     
                  end if
c     
c     IF liquid is entering the system
c     
                  if( srmiml .lt. 0. ) then
c     
c     IF solute can enter with the liquid
c     
                     if( icns(nsp) .gt. 0 ) then
c     
c     IF the current time is during the solute injection period
c     
                        if( days .gt. abs(t1sk(mi)) .and.
     2                       days .le. abs(t2sk(mi)) ) then
c     
c     Compute contribution to solute source term
c     
                           
                           rcss(mi) = getflux * srmiml
                           rc_ss = rcss(mi)
                           drc_ss = 0.
c     
c     ENDIF
c     
                        end if
c     
c     ENDIF
c     
                     end if
                     
                  end if
               endif
               
               
               
            endif
         else if ((iz.eq.2).and.(usroption.eq.5)) then
c     Recirculation option
            if (l .ne. 0) then
               mi = i + npn
               ispecies = 1 + npn/n0
               do k = 1, nu_zones
                  if (izonef(i) .eq. inflzones(k)) then
                     if (abs(sk(i)) .gt. zero_t) then
                        rc_ss = rc_ss + moles_recycle(k, ispecies) 
     &                       * sk(i)
                     end if
c     Erosion model
                     if (erosion_factor .ne. 0 .and. scl(mi) .gt.
     &                    zero_t) then
C     Using rate/s 
                        rc_ss_erosion = 1.d-3 * erosion_factor 
     &                       * sx1(i) * rolf(i) * scl(mi)
                        drc_ss_erosion = 1.d-3 * erosion_factor 
     &                       * sx1(i) * rolf(i) * dsccl(mi)
                        rc_ss = rc_ss + rc_ss_erosion
                        drc_ss = drc_ss + drc_ss_erosion
                     end if
c     if (mi .eq. 7022) then
c     write(iptty,*) 'userc 7022:', rc_ss
c     end if
                     exit
                  end if
               end do
            end if

         else if ((iz.eq.3).and.(usroption.eq.5)) then
c     usroption = 5 Solute recirculation model
c     foreach production zone
            moles_recycle = 0.
            do k = 1, nu_zones
c     determine if the zone is currently producing
c     if producing
               inflrate = 0.
               do ispecies = 1, nspeci
                  molesout = 0.
                  npn = npt(ispecies)
                  do j = 1, n0
                     mi = j + npn
                     if (izonef(j) .eq. prodzones(k)) then
c     get the solute production rate (molesout -> moles/s)
                        if (sk(j) .gt. zero_t) then
                           molesout = molesout + sk(j) * anlo(mi)
                        end if
                     end if
                  end do
                  if (molesout .gt. zero_t .and. ispecies .eq. 1) then
c     determine infiltration for zone (not dependent on species)
                     do j = 1, n0
                        if (izonef(j) .eq. inflzones(k)) then
                           if (abs(sk(j)) .gt. zero_t) then
                              inflrate = inflrate + abs(sk(j))
                           end if
                        end if
                     end do
                  end if
c     molesin -> moles/kg
                  if (inflrate .gt. zero_t) then
                     moles_recycle(k, ispecies) = (molesout / inflrate)
     &                    * recycle_factor(ispecies)
                  end if
               end do
            end do
         endif
         
      end if
      
      else
c---------------------------------------------------------
c     Phil and Hari  6/2004 - 4/2011
c     This section is for Tracer transport coupled 
c     with Goldsim.  The input flux is calculated from
c     the mass/time 
C     i  also, feed in kd from in(*) array
c     AND check to see if deltat is negative and reset
c     meaning that a new realization is underway!
c---------------------------------------------------------
c---  Changes to SUPER INDEXING  9 28 04  IN(8+in(7)-1+nsp)

         if(in(4).EQ.666) then 	 

c     - - - - - - only open and read the transform matrix once.
            if (flag_kd.NE.666) then
               open(222,file='transformation_matrix.txt')
               read(222,*) 
c     ii=cluster  jj=flowfield 
               read(222,*) ii,jj
               write(iaunit,*) 'Transform Matrix  column=cluster row=FF'
               do j1 = 1,jj
                  read(222,*) (transform_kd(i1,j1) , i1=1,ii)
                  write(iaunit,771) (transform_Kd(i1,j1), i1=i,ii)
               end do
               close(222)
               flag_kd=666
            end if

            if(in3_old.NE.in(3)) flag_kd2 = 0

            if((i.eq.1).AND.(in(1).GT.0)) then

               if(flag_kd2.NE.1) then
                  write(iaunit,*) 'Transformation in USERC.f '
                  write(iaunit,*) 'Cluster Flowfield NSP Kd Factor'
                  tf1 = 666
               end if

               if((iadsfl(nsp,1).eq.66).and.(flag_kd2.ne.1))then
		  a1adfl(nsp,1) = in(7+int(in(7))+nsp)*
     &                 transform_kd(int(in(6)),int(in(2)))
                  write(iaunit,772) int(in(6)),int(in(2)),nsp,
     &                 a1adfl(nsp,1),transform_kd(int(in(6)),int(in(2)))
                  if(nsp.eq.nspeci) then
                     flag_kd2 = 1
                     in3_old = in(3)
                  end if
               endif
               
 771           format(10(F6.3)) 
 772           format(3I6,2x,2G9.3)

c--------4/27/2011  in(5) is now delta T direct from GS in yrs
               
               deltat = in(5)*3600.*24.*365.25

c     SPC comment: this is for areaG run with Kd specified
c     SPC comment: With no Kd option, shall use follow line
c     getflux = in(11+ nspeci+sflag)/deltat    

               sflag = int(in(7+nsp))
               getflux = in(11+2*nspeci+sflag)/deltat

               rc_ss   = -getflux
               drc_ss  = 0.

            endif               ! i EQ 1  
	 endif                  ! in(4) EQ 666
      endif
      
      return
      end
