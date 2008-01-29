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
!!D2 
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
      implicit  none

      integer n_points,jl,ju,jm,iz,iuserc,ichk,i2,usroption
      real*8 slope,getconc,timein,timeCl36start,getflux
      integer, optional :: i
      real(8), optional :: rc_ss, drc_ss
      integer ispecies, j, k, k1, curcolumn
c      pointer(ipuserconc,userconc)
c      real*8 userconc(99999)
      real(8), allocatable :: userconc(:,:)
      real(8), allocatable :: userconc3(:,:,:)
c      pointer(ipusertime,usertime)
c      real*8 usertime(99999)
      real(8), allocatable :: usertime(:)
      integer, allocatable :: nsindex(:)
      integer, allocatable :: nodeindex(:)
c      pointer(ipsrcnodes,srcnodes)
c      integer srcnodes(99999)
      integer nsnodes, iread1, nszones, mi, mim
      real*8 srmiml
      logical used
      integer icount, open_file, flag_user
      character*100 filename

      save timeCl36start,iuserc,ichk,usroption,n_points
      save usertime, userconc, userconc3
      save nodeindex, nsindex

      if ((iz.eq.0).and.(iuserc.ne.-9324)) then
         iuserc=-9324
         if (iout .ne. 0) 
     &      write (iout,*) 'Solute transport user subroutine is invoked'
         if (iatty.gt.0) 
     &      write(iatty,*)'Solute transport user subroutine is invoked'
         

         do icount = 1, 100
            filename(icount:icount) = ' '
         end do
c zvd 04/26/2007 Allow user to enter name of file for data
         read (inpt,'(a100)') filename
         if (filename(1:4) .eq. 'file') then
            read (inpt,'(a100)') filename
         else
            backspace inpt
            filename = ''
            filename(1:14) = 'userc_data.dat'
         end if
         iread1 = open_file(trim(filename),'old')

         
c-------------------------------------------------------------
c-------  temp change to allow list of nodes then time conc
c-------  that applies to entire list, for NTS tritium 
c-------   phs  2/19/2004 
c-------------------------------------------------------------
            read(iread1,*) usroption
            if(usroption.eq.3) then
               read(iread1,*) n_points, nsnodes
               if(.not.allocated(nodeindex)) then
                  allocate(nodeindex(n0), nsindex(nsnodes))
               end if
               nodeindex = 0
               read(iread1,*)(nsindex(j),j=1,nsnodes)
               do j = 1, nsnodes
                  nodeindex(nsindex(j)) = j
               end do
               if(.not.allocated(userconc3)) then
                  allocate(userconc3(nspeci,nsnodes,n_points),
     2                 usertime(n_points))
               end if
               do k = 1, nspeci
                  do k1 = 1, n_points
                     read(iread1,*) usertime(k1),
     2                    userconc3(k,1,k1)                  ! phs
                     do j = 2,nsnodes                        ! phs
                       userconc3(k,j,k1) = userconc3(k,1,k1) ! phs
                     enddo                                   ! phs
                  end do
               end do
            elseif(usroption.eq.4) then
               read(iread1,*) n_points, nszones
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
                  allocate(userconc3(nspeci,nszones,n_points),
     2                 usertime(n_points))
               end if
               do k = 1, nspeci
                  do k1 = 1, n_points
                     read(iread1,*) usertime(k1),
     2                    (userconc3(k,j,k1),j=1,nszones)
                  end do
               end do
            else
               read(iread1,*) n_points
               if(.not.allocated(userconc)) then
                allocate(userconc(nspeci,n_points), usertime(n_points))
               end if
               read(iread1,*)
               do i2=1,n_points
                  read(iread1,*) (userconc(j,i2),j=1,nspeci)
               enddo
               read(iread1,*)
               do i2=1,n_points
                  read(iread1,*) usertime(i2)
               enddo
            end if
            close(iread1)
         else if (iuserc.eq.-9324) then
            if ((iz.eq.1).and.(usroption.eq.1)) then
               if (l.ne.0) then
                  if (ichk.ne.-9876) then 
                     do i2=1,n0
                        if (cnsk(i2).lt.0.) timeCl36start=-cnsk(i2)
                     enddo
                     ichk=-9876
                     timeCl36start=timeCl36start*365.25*8.64e4
                  endif
                  timein=timeCl36start-days*8.64e4
                  if (timein.le.0.) then 
                     getconc=userconc(1,1)
                     goto 1
                  else if (timein.ge.usertime(n_points)) then
                     i2=n_points
                     goto 2
                  endif
                  
                  jl=0
                  ju=n_points+1
 10               if (ju-jl.gt.1) then
                     jm=(ju+jl)/2
                     if(timein.gt.usertime(jm))then
                        jl=jm
                     else
                        ju=jm
                     endif
                     goto 10
                  endif
                  i2=jl+1
                  
 2                if(i2.eq.1) then
                     getconc=userconc(1,1)
                  else
                     slope=userconc(1,i2)-userconc(1,i2-1)
                     slope=slope/(usertime(i2)-usertime(i2-1))
                     getconc=userconc(1,i2-1)+slope*
     2                    (timein-usertime(i2-1))
                  endif
                  
 1                do i2=1,n0
                     if (cnsk(i2).lt.0.) cnsk(i2)=-getconc
                  enddo
               endif
               
            else if ((iz.eq.2).and.(usroption.eq.2)) then
               if (l.ne.0) then
                  if (ichk.ne.-9876) then
                     ispecies = 1 + npn/n0
                     if(cnsk(i+npn) .eq. -9876) then
                        timein=days*8.64e4
                        if (timein.le.usertime(1)) then 
                           getflux=userconc(ispecies,1)
                           goto 13
                        else if (timein.ge.usertime(n_points)) then
                           i2=n_points
                           goto 23
                        endif
                        
                        jl=0
                        ju=n_points+1
 103                    if (ju-jl.gt.1) then
                           jm=(ju+jl)/2
                           if(timein.gt.usertime(jm))then
                              jl=jm
                           else
                              ju=jm
                           endif
                           goto 103
                        endif
                        i2=jl+1
                        
 23                     if (i2.eq.1) then
                           getflux=userconc(ispecies,1)
                        else
                           slope=userconc(ispecies,i2)-
     2                          userconc(ispecies,i2-1)
                           slope=slope/(usertime(i2)-usertime(i2-1))
                           getflux=userconc(ispecies,i2-1)+slope*
     2                          (timein-usertime(i2-1))
                        endif
 13                     continue
                        rc_ss = -getflux
                        drc_ss = 0.
                     end if
                  endif
                  
                  
                  
               endif
            else if ((iz.eq.2).and.(usroption.eq.3)) then
C---------------------------------------------------------
c------- phil added fix to bypass nodes not in the nodeindex
c                       2/18/04
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
                     else if (timein.ge.usertime(n_points)) then
                        i2=n_points
                        goto 33
                     endif
               
                     jl=0
                     ju=n_points+1
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
                  
               endif
              endif   ! (nodeindex.NE.0)  phs 2/18/04
            else if ((iz.eq.2).and.(usroption.eq.4)) then
               if (l.ne.0) then
                  ispecies = 1 + npn/n0
                  curcolumn = nodeindex(i)
                  if(pcnsk(i+npn) .lt. 0.) then
                     timein=days
                     if (timein.le.usertime(1)) then 
                        getflux=userconc3(ispecies,curcolumn,1)
                        goto 143
                     else if (timein.ge.usertime(n_points)) then
                        i2=n_points
                        goto 133
                     endif
                     
                     jl=0
                     ju=n_points+1
 403                 if (ju-jl.gt.1) then
                        jm=(ju+jl)/2
                        if(timein.gt.usertime(jm))then
                           jl=jm
                        else
                           ju=jm
                        endif
                        goto 403
                     endif
                     i2=jl+1
                     
 133                 if (i2.eq.1) then
                        getflux=userconc3(ispecies,curcolumn,1)
                     else
                        slope=userconc3(ispecies,curcolumn,i2)-
     2                       userconc3(ispecies,curcolumn,i2-1)
                        slope=slope/(usertime(i2)-usertime(i2-1))
                        getflux=userconc3(ispecies,curcolumn,i2-1)+
     2                       slope*(timein-usertime(i2-1))
                     endif
 143                 continue
c     
c   IF this is an isothermal air-water simualtion
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
c   ENDIF
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
     2                          days .le. abs(t2sk(mi)) ) then
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
            endif
            

         end if
             
         return
         end



