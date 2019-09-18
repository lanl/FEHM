      subroutine pest(iflg)
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
!D1 Parameter ESTimation IO routine
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM Version 2.0, SC-194
!D2
!D2 Initial implementation: 22-JAN-1997, Programmer: G. Zyvoloski
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/pest.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:38   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:11:52   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:11:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:25:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:05:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:44:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2
!***********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3 
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS
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
c GAZ 1/22/97
c Major modification GAZ 7/30/98
C***********************************************************************

      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comflow
      use comxi
      use davidi

      implicit none

      real*8 value,valuex,valuey,valuez
      real*8 tolp, head_ini 

      real*8, allocatable :: xpest(:)
      real*8, allocatable :: ypest(:)
      real*8, allocatable :: zpest(:)
      real*8, allocatable :: fperms(:,:)
      real*8 sumfout, sumsink, sumsource, sumboun, xc, yc, zc
c gaz 080419      
      real*8 avgsat, volsat, s_dum, vel_mag
      integer, allocatable :: nelm_pest(:)
      integer, allocatable :: npestq(:)
      integer, allocatable :: npestc(:)
      integer, allocatable :: iparam(:,:)
      integer i, iflg, inode, inneq, iroot, izone, k, nodew
      integer idummy, addnode, i1, i2, iconn, indexa_axy, nparam
      integer itp,inump,mpestq,mpestc
      logical matrix_node, null1
      character*80 line, pest_param_info
      character*80 pest_deriv_ini, pest_deriv_fin           

      save xpest,ypest,zpest,tolp,fperms        
      save nelm_pest,npestq,npestc              
      save iparam,itp,nparam,mpestq,mpestc
      if(ipest.eq.0) return
      if(iflg.eq.0) then

c     Assign file unit numbers and determine output filenames
         ispest = nufilb(15)
         ispst1 = nufilb(16)
         if (null1(root_name)) then
            if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
c     Prefix from output file name
               call file_prefix(nmfil(5), iroot)
               if (iroot .gt. 94) iroot = 94
               nmfil(15) = nmfil(5)(1:iroot) // suffix(15)
               nmfil(16) = nmfil(5)(1:iroot) // suffix(16)
            else 
c     Use default filenames "fehmn.pest", "fehmn.pest1"
               if (nmfil(15)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil15 file: ',
     .                 nmfil(15), ' unable to determine pest file name'
                  stop
               end if
               if (nmfil(16)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil16 file: ',
     .                 nmfil(16), ' unable to determine pest1 file name'
                  stop
               end if
            endif
         else
            iroot = len_trim (root_name)
            if (iroot .gt. 94) iroot = 94
            nmfil(15) = root_name(1:iroot) // suffix(15)
            nmfil(16) = root_name(1:iroot) // suffix(16)
         end if
         open(ispest, file = nmfil(15), status = cstats(15))
         write(ispest, 1000) verno, jdate, jtime
         open(ispst1, file = nmfil(16), status = cstats(16))
         write(ispst1, 1000) verno, jdate, jtime
 1000    format('PEST Output: ', a30, 3x, a11, 3x, a8)
         
c     input section
         read (inpt,'(a80)') wdd1
         mpestq = 0
         mpestc = 0
         read (wdd1,*,end=10) mpest,mpestq,mpestc    
 10      continue
         if(.not.allocated(npest)) then
            allocate(npest(mpest))
            if(mpestq.ne.0) allocate(npestq(mpestq))
            if(mpestc.ne.0) allocate(npestc(mpestc))
            allocate(xpest(mpest))
            allocate(ypest(mpest))
            allocate(zpest(mpest))
            if(icnl.eq.0) then
               allocate(nelm_pest(mpest*8))
            else
               allocate(nelm_pest(mpest*4))
            endif
            if(.not.allocated(izoneflxz)) then
               allocate(izoneflxz(n0))
            endif
         end if
         nelm_pest=0

         read (inpt, *) (npest(i), i =1, mpest)
c**** read in coordinates
         do  i = 1,mpest 
            if (npest(i) .eq.-999) then
               read (inpt, *) xc, yc, zc
               call near3 (xc, yc, zc, nodew, 0)
               npest(i) =  nodew
            else if (npest(i) .lt. 0) then
               read (inpt, *) xpest(i), ypest(i), zpest(i)
            end if
         end do
c**** read in zones for flux targets
         if( mpestq .ne. 0 ) then
            read (inpt, *) (npestq(i), i =1, mpestq)
         end if
c**** read in zones for concentration targets
         if( mpestc .ne. 0 ) then
            read (inpt, *) (npestc(i), i =1, mpestc)
c**** read in coordinates for concentration
            do  i = 1,mpestc 
               if (npestc(i) .eq.-999) then
                  read (inpt, *) xc, yc, zc
                  call near3 (xc, yc, zc, nodew, 0)
                  npestc(i) =  nodew
               else if(npestc(i).lt.0) then
                  write (ierr, 667) 
                  if (iptty .ne. 0) write(iptty, 667) 
                  stop
               end if
            end do
         end if
 667     format ('>>> only -999 option allowed for concentration',
     &        ' (stopping)')
c**** read in info on derivatives
         nparam = 0
         read(inpt,'(a80)',end=999) line
         if(line(1:4).eq.'deri') then
            read(inpt,*,end=999) nparam    
            allocate(iparam(nparam,2))
            allocate(fperms(nparam,3))
            backspace inpt
            read(inpt,*,end=999) nparam
            do i =1,nparam
               read(inpt,*,end=999) iparam(i,1),
     &              iparam(i,2),fperms(i,1),fperms(i,2),fperms(i,3) 
            enddo
            do i= 1,80
               pest_deriv_ini(i:i) = ' '
               pest_deriv_fin(i:i) = ' '
            enddo
            read(inpt,*) tolp,itp,inump        
            read(inpt,'(a80)') pest_param_info
            read(inpt,'(a80)',end=998) pest_deriv_ini 
            read(inpt,'(a80)',end=998) pest_deriv_fin 
            go to 999
 998        continue                 
         else
            backspace inpt 
         endif
 999     continue

      else if(iflg.eq.-1) then

c**** if npest(i) < 0 and find node number(or element) ****
         do  i = 1,mpest 
            if (npest(i) .lt. 0.and.npest(i).ne.-999) then
               call pest_elem(0, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,value,phi)
               if(k.eq.0) then
                  npest(i)=0
               endif
            end if
         end do

      else if(iflg.eq.1) then
         if(istea_pest.eq.0) then
c     output section
            write(ispest,1010) days
            write(ispst1,1010) days
 1010       format('Time = ',g24.12,' days', ' Steady_State')
            write(ispest,1011) mpest
            write(ispst1,1011) mpest
 1011       format('Nodes = ',i8)
         else
c     output section
            write(ispest,1012) days
            write(ispst1,1012) days
 1012       format('Time = ',g24.12,' days',' Transient ')
            write(ispest,1011) mpest
            write(ispst1,1011) mpest
	 endif
         if(ihead.eq.0) then
            if (istea_pest.eq.0) then
               write(ispest,*) 'pressure     zone'
            else
               write(ispest,*) 'pressure  pressure difference  zone'
	    endif
            do  i = 1,mpest 
               if(npest(i).gt.0) then
                  if (istea_pest.eq.0) then
                     write(ispest,'(i10,f40.20,1x,i7)') npest(i),
     &                phi(npest(i)),izonef(npest(i))
                  else
                   write(ispest,'(i10,2f40.20,1x,i7)') npest(i),
     &              phi(npest(i)),phi(npest(i))-phini(npest(i)),
     &              izonef(npest(i))
                  endif
                  if (rlp_flag .eq. 1) then
                     write(ispst1,'(i8,1x,i4,1x,3(g20.12,1x))') 
     &                    npest(i),irlp(npest(i)), phi(npest(i)), 
     &                    s(npest(i)), t(npest(i))
                  else
                     write(ispst1,'(i8,1x,i4,1x,3(g20.12,1x))') 
     &                    npest(i),0, phi(npest(i)), s(npest(i)), 
     &                    t(npest(i))
                  end if
               else if(npest(i).lt.0) then
                  call pest_elem(1, i, nelm_pest, xpest(i), 
     &                 ypest(i), zpest(i),k,value,phi)
                  write(ispest,'(i10,f40.20)') npest(i),value
               else
                  write(ispest,'(i10,f40.20,5x,
     &                 21hobservation not found)') 0,0.d00
               endif
            end do
         elseif(ihead.ne.0) then
            if (istea_pest.eq.0) then
               write(ispest,*) 'heads'
            else
               write(ispest,*) 'heads and head differences'
            endif	
            do  i = 1,mpest 
               if(npest(i).gt.0) then
                  if (istea_pest.eq.0) then
                     call headctr(4,npest(i),phi(npest(i)),
     &                    head(npest(i)))
                     write(ispest,'(i10,f40.20)') npest(i),
     &                    head(npest(i))
                  else
                     call headctr(4,npest(i),phi(npest(i)),
     &                    head(npest(i)))
                     call headctr(4,npest(i),phini(npest(i)),head_ini)
                     write(ispest,'(i10,2f30.15)') npest(i),
     &                    head(npest(i)),head(npest(i))-head_ini    
                  endif
                  if (rlp_flag .eq. 1) then
                     if (ifree .ne. 0) then
                        write(ispst1,'(i8,1x,i4,1x,3(g20.12,1x))') 
     &                       npest(i), irlp(npest(i)), head(npest(i)), 
     &                       s(npest(i)), t(npest(i))
                     else
                        write(ispst1,'(i8,1x,i4,1x,3(g20.12,1x))') 
     &                       npest(i), irlp(npest(i)), head(npest(i)), 
     &                       1.0d0, t(npest(i))
                     end if
                  else
                     write(ispst1,'(i8,1x,i4,1x,3(g20.12,1x))') 
     &                    npest(i),0, head(npest(i)), 1.0d0, t(npest(i))
                  end if
               else if(npest(i).lt.0) then
                  call pest_elem(1, i, nelm_pest, xpest(i), 
     &                 ypest(i), zpest(i),k,value,head)
                  write(ispest,'(i10,f40.20)') npest(i),value
               else
                  write(ispest,'(i10,f40.20,5x,
     &                 21hobservation not found)') 0,0.d00
               endif
            end do
         endif                       
         write(ispest,*) 'saturations'
         do  i = 1,mpest 
            if(npest(i).gt.0) then
               if (irdof .ne. 13 .or. ifree .ne. 0) then
                  write(ispest,'(i10,f40.20)') npest(i),s(npest(i))
               else
                  write(ispest,'(i10,f40.20)') npest(i),1.0d0
               end if
            else if(npest(i).lt.0) then
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,value,s)
               write(ispest,'(i10,f40.20)') npest(i),value
            else
               write(ispest,'(i10,f40.20,5x,21hobservation not found)')
     &              0,0.d00               
            endif
         end do
         write(ispest,*) 'temperature'
         do  i = 1,mpest 
            if(npest(i).gt.0) then
               write(ispest,'(i10,f40.20)') npest(i),t(npest(i))
            else if(npest(i).lt.0) then
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,value,t)
               write(ispest,'(i10,f40.20)') npest(i),value
            else
               write(ispest,'(i10,f40.20,5x,21hobservation not found)')
     &              0,0.d00               
            endif
         end do
         write(ispest,*) 'permeabilities'
         do  i = 1,mpest 
            if(npest(i).gt.0) then
               write(ispest,'(i10,1p,3g18.10)') npest(i),
     &              pnx(npest(i))*1.d-6,pny(npest(i))*1.d-6,
     &              pnz(npest(i))*1.d-6
            else if(npest(i).lt.0) then
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuex,pnx)
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuey,pny)
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuez,pnz)
               write(ispest,'(i10,1p,4g18.10)') npest(i),valuex*1.d-6,
     &              valuey*1.e-6,valuez*1.d-6               
            else
               write(ispest,'(i10,f40.20,5x,21hobservation not found)')
     &              0,0.d00               
            endif
         end do
         write(ispest,*) 'velocities'
         do  i = 1,mpest 
            if(npest(i).gt.0) then
               vel_mag =sqrt(pnx(npest(i)+n)**2 + pny(npest(i)+n)**2 + 
     &          pnz(npest(i)+n)**2)
               write(ispest,'(i10,1p,4g18.10)') npest(i),
     &              pnx(npest(i)+n),pny(npest(i)+n),
     &              pnz(npest(i)+n), vel_mag
            else if(npest(i).lt.0) then
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuex,pnx(n+1))
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuey,pny(n+1))
               call pest_elem(1, i, nelm_pest, xpest(i), 
     &              ypest(i), zpest(i),k,valuez,pnz(n+1))
               vel_mag =sqrt(valuex**2 + valuey**2 + 
     &          valuez**2)               
               write(ispest,'(i10,1p,4g18.10)') npest(i),valuex,
     &              valuey, valuez                
            else
               write(ispest,'(i10,f40.20,5x,21hobservation not found)')
     &              0,0.d00               
            endif
         end do

c     Compute flux passing thorugh a zone

         if(nflxz.ne.0) then

            write(ispest,*)
            write(ispest,*)
            write(ispest,1040) nflxz
            write(ispst1,1041) nflxz
 1040       format('Zones = ',i8, '      Total Flux (kg/s) ', 
     &           'Leaving Zone (flxz macro option)')
            write(ispest,*)
 1041       format('Zones = ',i8)
            write(ispest,1044)
 1044      format(t3,'Zone',t14,'Source',t29,'Sink',t44,'Net',
     &      t59,'Boundary',t74,'Avg Saturation')               
            write(ispest,*)
c     Compute fluxes out of zone (>0), leaving out fluxes
c     into other parts of the zone

            do izone = 1, nflxz
               sumfout = 0.
               sumsink = 0.
               sumsource = 0.
               sumboun = 0.0
               avgsat = 0.0
               volsat = 1.d-30
c     Loop over all nodes
               do inode = 1, n0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
                  if(inode.gt.neq) then
                     inneq = inode-neq
                     matrix_node = .true.
                     addnode = nelm(neq+1)-neq-1
                     idummy = neq
                  else
                     matrix_node = .false.
                     inneq = inode
                     addnode = 0
                     idummy = 0
                  end if
c gaz 051519  (repeated 080419)       
                     if(irdof.eq.13) then
                         s_dum = 1.
                     else
                         s_dum = s(inneq)
                     endif
c     Determine if node is part of the zone being summed
                  if(izoneflxz(inode).eq.izone) then
c     Add boundary condition sources
                     sumboun=sumboun + sk(inode)
c gaz 080419 avg sat calc
                     volsat = volsat + volume(inneq)*ps(inneq)
                     avgsat = avgsat + volume(inneq)*ps(inneq)*s_dum
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
                     i1 = nelm(inneq)+1
                     i2 = nelm(inneq+1)
c     loop over every connecting node
                     do iconn = i1, i2
                        indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
                        if(a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                           if(izoneflxz(idummy+nelm(iconn))
     2                          .ne.izone.or.nelm(iconn)
     3                          .eq.inneq) then
                              sumfout = sumfout + a_axy(indexa_axy)
                              if(nelm(iconn).eq.inneq) then
                                 sumsink = sumsink +
     2                                a_axy(indexa_axy)
                              end if
                           end if
                        elseif(a_axy(indexa_axy).lt.0.) then
c     add to source sum
                           if(nelm(iconn).eq.inneq) then
                              sumsource = sumsource +
     2                             a_axy(indexa_axy)
                           end if
                        end if
                     end do
                  end if
               end do
c     Write results
c gaz 080419
               avgsat = avgsat/volsat
               write(ispest,1045) iflxz(izone),
     2              sumsource, sumsink, sumfout, sumboun, avgsat
               write(ispst1,1045) iflxz(izone),
     2              sumsource, sumsink, sumfout, sumboun, avgsat
 1045          format(1x,i5,5x,5(2x,e15.7))
            end do

         end if
c     
c     print out concentrations
c     

         if(mpestc.ne.0.and.iccen.ne.0) then
            if(allocated(anl)) then
               write(ispest,*) 'concentrations'
               do  i = 1,mpestc 
                  if(npestc(i).gt.0) then
                     write(ispest,'(i10,3g18.10)') npestc(i),
     $                    max(anl(npestc(i)),1.d-90),
     &                    max(anv(npestc(i)),1.d-90)
                  else if(npestc(i).lt.0) then
                     call pest_elem(1, i, nelm_pest, xpest(i), 
     &                    ypest(i), zpest(i),k,valuex,anl)
                     call pest_elem(1, i, nelm_pest, xpest(i), 
     &                    ypest(i), zpest(i),k,valuey,anv)
                     write(ispest,'(i10,3g18.10)') npest(i),
     &                    max(valuex,1.d-90),max(valuey,1.d-90)  
                  else
                     write(ispest,'(i10,f40.20,5x,
     &                    21hobservation not found)')
     &                    0,0.d00               
                  endif
               end do
            endif
         endif

      else if(iflg.eq.2) then
c     
c     calculate derivatives if necessary
c     
         if(nparam.ne.0) then
            if(inump.lt.0) then
               if (iout .ne. 0) then
                  write(iout,*)
                  write(iout,*) '>>> calculating ',
     &                 nparam,' analytical sensitivities'        
                  write(iout,*)
               end if
               if (iatty .ne. 0) then
                  write(iatty,*)
                  write(iatty,*)  '>>> calculating ',
     &                 nparam,' analytical sensitivities'
                  write(iatty,*)
               end if
            else
               if (iout .ne. 0) then
                  write(iout,*)
                  write(iout,*) '>>> calculating ',
     &                 nparam,' central difference sensitivities'
                  write(iout,*)
               end if
               if (iatty .ne. 0) then
                  write(iatty,*)
                  write(iatty,*) '>>> calculating ',
     &                 nparam,' central difference sensitivities'
                  write(iatty,*)
               end if
            endif
            do i = 1,mpestq
               do k=1,nflxz
                  if(npestq(i).eq.iflxz(k)) then
                     npestq(i)=k
                     go to 5
                  endif
               enddo
 5             continue
            enddo
            if(inump.lt.0) then
               call parameter_sensitivity_analyt(nparam,
     &              pest_param_info,pest_deriv_ini,pest_deriv_fin,
     &              iparam,fperms,tolp,itp,mpestq,npestq)
            else
               call parameter_sensitivity_numerical(nparam, 
     &              pest_param_info,pest_deriv_ini,pest_deriv_fin,
     &              iparam,fperms,tolp,itp,mpestq,npestq)
            endif
            if (iout .ne. 0) then
               write(iout,*)
               write(iout,*) '>>> finished parameter sensitivities'
               write(iout,*)
            end if
            if (iatty .ne. 0) then
               write(iatty,*)
               write(iatty,*) '>>> finished parameter sensitivities'
               write(iatty,*)
            end if
         endif
         if(allocated(iparam)) deallocate(iparam)
         if(allocated(fperms)) deallocate(fperms)
         close(ispest)
         close(ispst1)
      endif
      return
      end

      subroutine pest_elem(iflg,ii,nelm_pest, xg, yg, zg, k, value,
     &     value_array)
c     this subroutine interpolates x,y,z from element data
C***********************************************************************
      use combi
      use comdti
      use comai
      implicit none

      logical eleb
      integer i, ij, j, in, izone, jz, maxitz, nin, nsl
      integer iflg,ii,k,nelm_pest(*)
      real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8 delx, dely, delz, detja
      real*8 etad, excid, rx, ry, rz, sid, sxd, sy, sz
      real*8 sa11, sa12, sa13, sa21, sa22, sa23, sa31, sa32, sa33
      real*8 tols, tolt, xl, xmn, xmx, yl, ymn, ymx, zl, zmn, zmx
      real*8 xg, yg, zg
      real*8 xz(8), yz(8), zz(8)
      real*8 value,value_array(*)

      maxitz = 10
      tols = 1.0e-07
      tolt = 1.0e-07
      nsl=ns
      k = 0

      if(iflg.eq.0) then
         if (icnl .eq. 0)  then
c**** 3-d calculation ****
c**** find coordinates contained in zone ****
c**** determine element type           ****
c**** find extreme coordinates in zone ****
c**** loop over elements
            do j =1, nei
               xmn = 1.0e+20
               xmx = -1.0e+20
               ymn = 1.0e+20
               ymx = -1.0e+20
               zmn = 1.0e+20
               zmx = -1.0e+20
               do in = 1, nsl
                  ij=nelm((j-1)*ns+in)
                  xmn = min(cord(ij,1), xmn)
                  xmx = max(cord(ij,1), xmx)
                  ymn = min(cord(ij,2), ymn)
                  ymx = max(cord(ij,2), ymx)
                  zmn = min(cord(ij,3), zmn)
                  zmx = max(cord(ij,3), zmx)
               end do
               
               nin=0

c**** find nodal coordinates belonging in this zone ****
               xl = xg            
               if (xl .ge. xmn-tols .and. xl .le. xmx+tols)  then
                  yl = yg         
                  if (yl .ge. ymn-tols .and. yl .le. ymx+tols)  then
                     zl = zg           
                     if (zl .ge. zmn-tols.and. zl .le.zmx+tols)  then
                        nin =  1
                     end if
                  end if
               end if

c**** do calculations on nodal coordinates in this zone ****
               if(nin.ne.0) then

c**** use newton raphson iteration to find local coordinates ****
c**** set local coordinates initially to zero ****
                  sid = 0.0
                  etad = 0.0
                  excid = 0.0
                  
                  do jz = 1, nsl
                     ij = nelm((j-1)*ns+jz)
                     xz(jz)=cord(ij,1)
                     yz(jz)=cord(ij,2)
                     zz(jz)=cord(ij,3)
                  enddo

                  do iad = 1, maxitz
c**** evaluate shape functions ****
                     call  sfn3r   (sid, etad, excid)

c**** form jacobian ****
c**** define geometric coefficients ****
                     sxd = 0.0
                     sy = 0.0
                     sz = 0.0
                     sa11 = 0.0
                     sa12 = 0.0
                     sa22 = 0.0
                     sa21 = 0.0
                     sa13 = 0.0
                     sa23 = 0.0
                     sa31 = 0.0
                     sa32 = 0.0
                     sa33 = 0.0

                     do jz = 1, nsl
                        sxd = sxd +w (jz, 1)*xz(jz)
                        sy = sy  +w (jz, 1)*yz(jz)
                        sz = sz  +w (jz, 1)*zz(jz)
                        sa11 = sa11+wx(jz, 1)*xz(jz)
                        sa12 = sa12+wy(jz, 1)*xz(jz)
                        sa13 = sa13+wz(jz, 1)*xz(jz)
                        sa21 = sa21+wx(jz, 1)*yz(jz)
                        sa22 = sa22+wy(jz, 1)*yz(jz)
                        sa23 = sa23+wz(jz, 1)*yz(jz)
                        sa31 = sa31+wx(jz, 1)*zz(jz)
                        sa32 = sa32+wy(jz, 1)*zz(jz)
                        sa33 = sa33+wz(jz, 1)*zz(jz)
                     end do

c**** save jacobian information ****
                     a11 = sa22*sa33-sa23*sa32
                     a21 = -sa21*sa33+sa23*sa31
                     a31 = sa21*sa32-sa22*sa31
                     a12 = -sa12*sa33+sa13*sa32
                     a22 = sa11*sa33-sa13*sa31
                     a32 = -sa11*sa32+sa12*sa31
                     a13 = sa12*sa23-sa13*sa22
                     a23 = -sa11*sa23+sa13*sa21
                     a33 = sa11*sa22-sa12*sa21
                     detja = sa11*sa22*sa33+sa12*sa23*sa31+
     .                    sa13*sa21*sa32-sa13*sa22*sa31-
     .                    sa23*sa32*sa11-sa33*sa21*sa12

c**** form residuals ****
                     rx = sxd-xl
                     ry = sy -yl
                     rz = sz -zl

c**** get correction ****
                     delx = (a11*rx+a12*ry+a13*rz)/detja
                     dely = (a21*rx+a22*ry+a23*rz)/detja
                     delz = (a31*rx+a32*ry+a33*rz)/detja
                     sid = sid - delx
                     etad = etad - dely
                     excid = excid - delz

c**** if correction small end iteration ****
                     if (sqrt(delx*delx+dely*dely+delz*delz) 
     .                    .le. tols) go to 10
                  end do

c**** iteration did not converge so stop ****
                  izone = ii
                  write(ierr, 6000) izone
                  if (iout .ne. 0) write(iout, 6000)  izone
                  if (iptty .gt. 0)  write(iptty, 6000)  izone
 6000             format(/, 1x, 'iteration in zone did not ', 
     *                 'converge, izone = ', i10, ' please check',
     *                 ' icnl in macro CTRL ')
                  
                  stop

 10               continue

c**** check if local coordinates within the element ****
                  eleb = .TRUE.
                  if ((abs(sid) .gt. 1.0 + tolt) .or. 
     .                 (abs(etad) .gt. 1.0 + tolt) .or.  
     .                 (abs(excid) .gt. 1.0 + tolt))  then
                     eleb = .FALSE.
                  end if

c**** the node i is in zone izone, save information ****`
                  if (eleb) then
                     k= j                
                     do in=1,nsl
                        ij=nelm((k-1)*ns+in)
                        nelm_pest((ii-1)*ns+in)=ij
                     enddo
                     xg=sid
                     yg=etad
                     zg=excid
                  else
                     k= 0           
                     do in=1,nsl
                        nelm_pest((ii-1)*ns+in)=0 
                     enddo
                  end if

               end if

            end do

         else if (icnl .ne. 0)  then
c**** 2-d calculation ****
c**** find coordinates contained in zone ****
c**** determine element type           ****
c**** find extreme coordinates in zone ****
c**** loop over elements
            do j =1, nei
               xmn = 1.0e+20
               xmx = -1.0e+20
               ymn = 1.0e+20
               ymx = -1.0e+20
               do in = 1, nsl
                  ij=nelm((j-1)*ns+in)
                  xmn = min(cord(ij,1), xmn)
                  xmx = max(cord(ij,1), xmx)
                  ymn = min(cord(ij,2), ymn)
                  ymx = max(cord(ij,2), ymx)
               end do

               nin=0

c**** find nodal coordinates belonging in this zone ****
               xl = xg            
               if (xl .ge. xmn .and. xl .le. xmx)  then
                  yl = yg         
                  if (yl .ge. ymn .and. yl .le. ymx)  then
                     nin =  1
                  end if
               end if

c**** do calculations on nodal coordinates in this zone ****
               if(nin.ne.0) then

c**** use newton raphson iteration to find local coordinates ****
c**** set local coordinates initially to zero ****
                  sid = 0.0
                  etad = 0.0

                  do jz = 1, nsl
                     ij = nelm((j-1)*ns+jz)
                     xz(jz)=cord(ij,1)
                     yz(jz)=cord(ij,2)
                  enddo
                  
                  do iad = 1, maxitz
c**** evaluate shape functions ****
                     call  sfn2r   (sid, etad)

c**** form jacobian ****
c**** define geometric coefficients ****
                     sxd = 0.0
                     sy = 0.0
                     sa11 = 0.0
                     sa12 = 0.0
                     sa22 = 0.0
                     sa21 = 0.0

                     do jz = 1, nsl
                        sxd = sxd +w (jz, 1)*xz(jz)
                        sy = sy  +w (jz, 1)*yz(jz)
                        sa11 = sa11+wx(jz, 1)*xz(jz)
                        sa12 = sa12+wy(jz, 1)*xz(jz)
                        sa21 = sa21+wx(jz, 1)*yz(jz)
                        sa22 = sa22+wy(jz, 1)*yz(jz)
                     end do

c**** save jacobian information ****
                     a11 = sa22
                     a22 = sa11
                     a12 = -sa12
                     a21 = -sa21
                     detja = sa11*sa22-sa12*sa21

c**** form residuals ****
                     rx = sxd-xl
                     ry = sy -yl

c**** get correction ****
                     delx = (a11*rx+a12*ry)/detja
                     dely = (a21*rx+a22*ry)/detja
                     sid = sid -delx
                     etad = etad-dely

c**** if correction small end iteration ****
                     if (sqrt(delx*delx+dely*dely) .le. tols) go to 20

                  end do

c**** iteration did not converge so stop ****
                  write(ierr, 6000)  izone
                  if (iout .ne. 0) write(iout, 6000)  izone
                  if (iptty .gt. 0)  write(iptty, 6000)  izone
                  
                  stop
                  
 20               continue

c**** check if local coordinates within the element ****
                  eleb = .TRUE.
                  if ((abs(sid) .gt. 1.0 + tolt) .or. 
     .                 (abs(etad) .gt. 1.0 + tolt))  then
                     eleb = .FALSE.
                  end if

c**** the node i is in zone izone, save information ****`
                  if (eleb) then
                     k= j                
                     do in=1,nsl
                        ij=nelm((k-1)*ns+in)
                        nelm_pest((ii-1)*ns+in)=ij
                        xg=sid
                        yg=etad
                     enddo
                  else
                     k= 0           
                     do in=1,nsl
                        nelm_pest((ii-1)*ns+in)=0 
                     enddo
                  end if

               end if

            end do
         endif

      else if(iflg.eq.1) then
c     evaluate interpolated value for pest output
         if(icnl.eq.0) then
            call sfn3r(xg,yg,zg)
         else
            call sfn2r(xg,yg)
         endif
         value = 0.0
         do i=1,nsl
            value = value +w(i,1)*value_array(nelm_pest((ii-1)*ns+i))
         enddo
      endif
      return
      end
      subroutine parameter_sensitivity_analyt
     &     (nparam,pest_param_info,pest_deriv_ini,pest_deriv_fin,
     &     iparam,fperms,tolp,itp,mpestq,npestq)
C**********************************************************************
C     D1
C     D1 PURPOSE
C     D1
C     D1 To generate sensitivity matrix for parameters for PEST
C     D1
C**********************************************************************
      use comai
      use combi
      use comci
      use comdi
      use comgi
      use comji
      use comxi
      use comei
      use comfi
      use comdti
      use davidi
      implicit none


      real*8, allocatable :: piv_a(:)
      real*8, allocatable :: deriv(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: dum2(:,:)
      real*8, allocatable :: bp1(:)
      real*8 pive,fdum2,facr,tollr
      real*8 apiv,aij,tolp
      real*8 bpid,adbl,bdbl
      integer neqp1,cnum,np,ijj,ij,ipar,kb,ityp,nzmap
      integer ispst2,i,j,ii,id,i1,i2,jj,idg,ikd,nparam,itp
      integer ispst3,ispst4,nparam1,mpest1,lispst2,mpestq
      integer iparam(nparam,*)
      real*8 fperms(nparam,*)
      integer, allocatable :: node_param(:)
      integer, allocatable :: izone_map(:)
      integer npestq(*)
      character*80 pest_deriv_file
      character*80 parameter_file
      character*80 pest_param_info
      character*80 pest_deriv_ini, pest_deriv_fin
      logical cdev,cdev1,cdev2
      integer max_zones
      parameter(max_zones=1000)
      cnum = 1

c     coding : this routine should probably be called from
c     subroutine pest just before exiting run
c     loop on number of parameters
c     read parameter file
c     need to get parameter file name form pest macro
      allocate (deriv(neq))
      allocate (node_param(neq))
      ispst2 =ispst1 + 20
      inquire(file=pest_param_info,opened=cdev2)
      if(cdev2) open(ispst2,file=pest_param_info,status='old')
      if(pest_deriv_ini(1:10).ne.'          ') then
         ispst3 =ispst2 + 1
         ispst4 =ispst3 + 1
         open (ispst4,file='deriv.temp',status='unknown',
     &        form='unformatted')
         write(ispst4) nparam,mpest
         open (ispst3,file=pest_deriv_ini,status='unknown',
     &        form='unformatted')
 882    format(/,'>>> reading re-start (initial) sensitivities, file: ',
     &        a40)
 782     format(/,'>>> no re-start sensitivities available from file: ',
     &        a40)
         read(ispst3,end=887) nparam1,mpest1 
         if (iout .ne. 0) write(iout,882) pest_deriv_ini(1:40)
         if (iptty .ne. 0) write(iptty,882) pest_deriv_ini(1:40)
         cdev = .true.
         backspace ispst3
      endif
      go to 888
 887  cdev = .false.
      if (iout .ne. 0) write(iout,782) pest_deriv_ini(1:40)
      if (iptty .ne. 0) write(iptty,782) pest_deriv_ini(1:40)
c     
 888  continue
c     
c     try reading parameter info from pest_param_info file first
c     at pesent, it is used to relate additional zone numbers for a single 
c     parameter
c     
      if(cdev2) then
         allocate(izone_map(max_zones))
         do i=1,nparam
            izone_map(iparam(i,1)) = iparam(i,1)
         enddo
c     first line is text
         read(ispst2,'(a80)',end=770)
c     second line is number of relationships
c     next lines zone j is repaced by zone k
         read(ispst2,*) nzmap
         do i=1,nzmap
            read(ispst2,*) j,ii    
            izone_map(j) = ii
         enddo
         do i = 1,neq
            node_param(i) = izone_map(izonef(i))
         enddo
         deallocate(izone_map)
         go to 771
      endif
c     
c     assume form of zone file
c     
c     izonef(i) is the parameter asociated with node i
c     
 770  do i = 1,neq
         node_param(i) = izonef(i)
      enddo
c     
 771  if(cdev2) close(ispst2)
c     
      do i=1,80
         if(nmfil(16)(i:i).eq.'.') then
            pest_deriv_file(1:i)=nmfil(16)(1:i)
            pest_deriv_file(i+1:i+5)='deriv'
            go to 5
         endif
      enddo
 5    continue
      lispst2 = i+5
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
c     
c     get and keep jacobian matrix for flow solution
c     
      call airctr(3,0)
c     
      allocate(dum(4*neq))
      allocate(piv_a(neq))
c     
      a=0.0
c     
      neqp1 = neq+1
      do 101 id=1,neq
c     
c     decide on equation type
c     
         call geneq2(id)
         a(nelmdg(id)-neqp1+nmat(1))=
     &        a(nelmdg(id)-neqp1+nmat(1)) + dq(id)
 101  continue

c     
c     coding for 1dof waterflow only
c     
c     
      do id=1,neq
         np=nelmdg(id)
         i1=nelm(id)+1
         i2=nelm(id+1)
         apiv=a(np-neqp1+nmat(1))
         if(apiv.le.0.0) then
            if (iptty .ne. 0) then
               write(iptty,*) "notice: apiv le 0, stopping"
               write(iptty,*) id,np,apiv
            end if
            write(ierr,*) "notice: apiv le 0, stopping"
            write(ierr,*) id,np,apiv
            stop
         endif
         do ijj=i1,i2
            ij=nelm(ijj)
            aij=a(ijj-neqp1+nmat(1))/apiv
            a(ijj-neqp1)=aij
         enddo
         piv_a(id) = 1./apiv
      enddo

      write(ispst2,*) nparam,mpest
      if(cdev) then
         read(ispst3) nparam1,mpest1
         if(nparam1.ne.nparam.or.mpest1.ne.mpest) then
            cdev1=.false.
         else
            cdev1=.true.
            allocate(bp1(neq))
         endif
      endif
c     
c     save ILU factorization and use for all parameters 
c     just need iback = 0 for parameter 1
c     iback = 1 for all other parameters
c     set iback = 1 so factorization is performed once

      iback = 0
      do ii =1,nparam
c     loop on number of grid blocks(ie max number of observations)
         ipar = iparam(ii,1)
         ityp = iparam(ii,2)
         if(ityp.gt.0) then
            if (iout .ne. 0) write(iout,*)
     &           'calculating sensitivities for permeability ',ipar
            if (iatty .ne. 0) write(iatty,*)
     &           'calculating sensitivities for permeability ',ipar
         else if(ityp.lt.0) then
            if (iatty .ne. 0) write(iatty,*)
     &           'calculating sensitivities for perm factor  ',ipar
            if (iout .ne. 0) write(iout,*)
     &           'calculating sensitivities for perm factor  ',ipar
         else 
            if (iatty .ne. 0) write(iatty,*)
     &           'not calculating sensitivities for parameter  ',ipar
            if (iout .ne. 0) write(iout,*)
     &           'not calculating sensitivities for parameter  ',ipar
            if(cdev) read(ispst3,end = 535) (bp1(i),i= 1,neq)
 535        bp = -1.11000000000000d33
            iter = 0
            go to 995
         endif
         bp = 0.0
         do jj =1,neq
c     calculate new source term (reflecting d(obs)/d(par))
            call node_sensitivity
     &           (nparam,jj,node_param,ipar,deriv,fperms,ii)
         enddo
c     
c     call linear equation solver for sensitivities wrt parameter ii
c     
c     need to get unnormalized array a 
c     
         fdum2 = 0.0
         do id=1,neq
            pive = piv_a(id)        
            bp(id)=-bp(id)*pive
            fdum2=fdum2+bp(id)*bp(id)
         enddo
         fdum=sqrt(fdum2)

c     
c     if previous derivative file exists 
c     read previous sensitivities
c     
         if(cdev1) then
            read(ispst3,end = 540) (bp1(i),i= 1,neq)
            go to 550
 540        bp1 = 0.0d00
            if (iout .ne. 0) write(iout,*) 
     &           'end-of-file encountered in derivative ini file',
     &           ' (deriv set =0)'
            if (iatty .ne. 0) write(iatty,*)
     &           'end-of-file encountered in derivative ini file',
     &           ' (deriv set =0)'
 550        continue
c     
c     matrix multply bp= bp-a*bp1
c     to form augmented residual 
c     
            fdum2 = 0.0
            neqp1 = neq +1
            do i = 1, neq
               i1=nelm(i)+1
               i2=nelm(i+1)
               bpid = bp(i)
               do j=i1,i2
                  kb = nelm(j)
                  adbl = a(j-neqp1)
                  bdbl = bp1(kb)
c     bpid = bpid - a(j-neqp1)*bp1(kb)
                  bpid = bpid - adbl*bdbl                   
               enddo
               bp(i) = bpid
               fdum2=fdum2+bp(i)*bp(i)
            enddo
            fdum=sqrt(fdum2)
         endif

         facr=1.0
         if(tolp.ge.0.0) then
            tollr = tolp*fdum    
         else
            tollr = abs(tolp)
         endif
         
         
         if(fdum.gt.tollr) then
            iter=itp           
            call solve_new(neq,a,b,bp,nmat,nb,nrhs,nelm,nop,north
     *           ,tollr,irb,iirb,npvt,gmres,dum,piv
     *           ,h,c,ss,g,y,iter,iback,1,iptty,maxor,accm)
c     set iback = 1 so same factorization is used
            iback = 1
         else
            iter = 0
            if (iout .ne. 0) then
               write(iout,*) 
     &              'sensitivity residual smaller than tolerance: ',tolp
               write(iout,*) 
     &              'residual = ',fdum,' calculation not performed'
            end if
            if (iatty .ne. 0) then
               write(iatty,*) 
     &              'sensitivity residual smaller than tolerance: ',tolp
               write(iatty,*) 
     &              'residual = ',fdum,' calculation not performed'
            end if
         endif

c     
c     printout analytic derivatives to file
c     the natural form is d(obs)/d(par(ii))
c     
         if(cdev1) then
            do i=1,neq
               bp(i)= bp1(i)+bp(i)
            enddo 
         endif
 995     continue
         if(pest_deriv_ini(1:10).ne.'          ') then
            write(ispst4) (bp(i),i=1,neq)
         endif
         if (iatty .ne. 0) write(iatty,*) 'solver iterations = ', iter
         if (iout .ne. 0) write(iout,*) 'solver iterations = ', iter
         write(ispst2,10) (bp(npest(j)),j=1,mpest)
 10      format(5(1x,g21.14))

      enddo
      close(ispst2)
c     might need to re-order the derivative file
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
      read(ispst2,*) nparam1,mpest1
      allocate(dum2(nparam1,mpest1))
      do i=1,nparam1
         read(ispst2,*) (dum2(i,j),j=1,mpest1)
      enddo
      close(ispst2)
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
      write(ispst2,*) nparam1,mpest1
      do j=1,mpest1 
         write(ispst2,772) (dum2(i,j),i=1,nparam1)
         write(ispst2,*) ' '
 772     format(5(1x,g21.14))
      enddo
      deallocate(dum2)
      close(ispst2)
 880  format
     &     (/,'>>> writing re-start sensitivities, file: ',a30)
      write(iatty,881) pest_deriv_file(1:lispst2)
      write(iout,881) pest_deriv_file(1:lispst2)
 881  format
     &     (/,'>>> writing pest usable sensitivities, file: ',a30)
      if(pest_deriv_ini(1:10).ne.'          ') then
         if(cdev1) deallocate(bp1)
         close(ispst4)
      endif
      if(pest_deriv_fin(1:10).ne.'          ') then
         if (iatty .ne. 0) write(iatty,880) pest_deriv_fin
         if (iout .ne. 0) write(iout,880) pest_deriv_fin
         open(ispst2,file=pest_deriv_fin,status='unknown',
     &        form='unformatted')
         open(ispst4,file='deriv.temp',status='old',
     &        form='unformatted')
         read(ispst4) nparam1,mpest1
         write(ispst2) nparam1,mpest1
         do j=1,nparam1
            read(ispst4) (bp(i),i=1,neq)
            write(ispst2) (bp(i),i=1,neq)
         enddo
         close(ispst4,status='delete')
         close(ispst2)
      else
         if (iatty .ne. 0) write(iatty,880) 'deriv.temp'
         if (iout .ne. 0) write(iout,880) 'deriv.temp'
      endif
      return
      end

      subroutine parameter_sensitivity_numerical 
     &     (nparam,pest_param_info,pest_deriv_ini,pest_deriv_fin,
     &     iparam,fperms,tolp,itp,mpestq,npestq)
C**********************************************************************
C     D1
C     D1 PURPOSE
C     D1
C     D1 To generate sensitivity matrix for parameters for PEST
C     D1
C**********************************************************************
      use comai
      use combi
      use comci
      use comdi
      use comgi
      use comji
      use comxi
      use comei
      use comfi
      use comdti
      use davidi
      implicit none


      real*8, allocatable :: piv_a(:)
      real*8, allocatable :: deriv(:)
      real*8, allocatable :: dum(:)
      real*8, allocatable :: dum2(:,:)
      real*8, allocatable :: headp(:,:)
      real*8, allocatable :: flxp(:,:)
      real*8, allocatable :: bp1(:)
      real*8 pive,fdum2,facr,tollr,perm_fac
      real*8 apiv,aij,tolp, dfactor, pnx1
      real*8 sumfout, sumsink, sumsource 
      real*8 bpid,adbl,bdbl
      integer neqp1,cnum,np,ijj,ij,ipar,kb,ityp,nzmap
      integer ispst2,i,j,ii,id,i1,i2,jj,idg,ikd,nparam,itp
      integer ispst3,ispst4,nparam1,mpest1,lispst2
      integer mpestq,izone
      integer, allocatable :: node_param(:)
      integer, allocatable :: izone_map(:)
      integer iparam(nparam,*)
      integer npestq(*)
      real*8 fperms(nparam,*)
      character*80 pest_deriv_file
      character*80 parameter_file
      character*80 pest_param_info
      character*80 pest_deriv_ini, pest_deriv_fin
      logical cdev,cdev2
      integer max_zones, iterp
      parameter(perm_fac=1.d06)
      parameter(max_zones=1000)
      cnum = 1

c     coding : this routine should probably be called from
c     subroutine pest just before exiting run
c     loop on number of parameters
c     read parameter file
c     need to get parameter file name form pest macro
      allocate (deriv(neq))
      allocate (node_param(neq))
      ispst2 =ispst1 + 20
      inquire(file=pest_param_info,opened=cdev2)
      if(cdev2) open(ispst2,file=pest_param_info,status='old')
c     
c     try reading parameter info from pest_param_info file first
c     at pesent, it is used to relate additional zone numbers for a single 
c     parameter
c     
      if(cdev2) then
         allocate(izone_map(max_zones))
         do i=1,nparam
            izone_map(iparam(i,1)) = iparam(i,1)
         enddo
c     first line is text
         read(ispst2,'(a80)',end=770)
c     second line is number of relationships
c     next lines zone j is repaced by zone k
         read(ispst2,*) nzmap
         do i=1,nzmap
            read(ispst2,*) j,ii    
            izone_map(j) = ii
         enddo
         do i = 1,neq
            node_param(i) = izone_map(izonef(i))
         enddo
         deallocate(izone_map)
         go to 771
      endif
c     
c     assume form of zone file
c     
c     izonef(i) is the parameter asociated with node i
c     
 770  do i = 1,neq
         node_param(i) = izonef(i)
      enddo
c     
 771  if(cdev2) close(ispst2)
c     
      do i=1,80
         if(nmfil(16)(i:i).eq.'.') then
            pest_deriv_file(1:i)=nmfil(16)(1:i)
            pest_deriv_file(i+1:i+5)='deriv'
            go to 5
         endif
      enddo
 5    continue
      lispst2 = i+5
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
c     
      write(ispst2,*) nparam,mpest+mpestq
c     
      allocate(headp(2,n))
      if(mpestq.ne.0) then
         allocate(flxp(2,mpestq))
      endif
      maxsolve = itp
      if(tolp.lt.0) then
         fdum1 = -1.0
         tmch = abs(tolp)
      else
         fdum1 = 1.0
         tmch = tolp
      endif

      iback = 0
      do ii =1,nparam
c     loop on number of grid blocks(ie max number of observations)
         ipar = iparam(ii,1)
         ityp = iparam(ii,2)
         if(mpestq.ne.0) flxp = 0.0d00
         if(ityp.gt.0) then
            if (iatty .ne. 0) write(iatty,*)
     &           '>>>     calculating sensitivities for permeability ',
     &           ipar
            if (iout .ne. 0) write(iout,*)
     &           '>>>     calculating sensitivities for permeability ',
     &           ipar
         else if(ityp.lt.0) then
            if (iatty .ne. 0) write(iatty,*)
     &           '>>>     calculating sensitivities for perm factor  ',
     &           ipar
            if (iout .ne. 0) write(iout,*)
     &           '>>>     calculating sensitivities for perm factor  ',
     &           ipar
         else 
            if (iatty .ne. 0) write(iatty,*)
     &           'not calculating sensitivities for parameter  ', ipar
            if (iout .ne. 0) write(iout,*)
     &           'not calculating sensitivities for parameter  ', ipar
 535        bp = -1.11000000000000d33
            flxp = -1.11000000000000d33
            iter = 0
            go to 995
         endif
c     
c     generate head values for permeabilities
c     do 2 times for central differences
c     set porosity = 0   
         ps = 1.d-12
c     
         iterp = 0
         do jj=1,2
            if(jj.eq.1) then
               dfactor = -1.0d00/(10.d00)**(ityp)
            else
               dfactor = 2.0d00/(10.d00)**(ityp)
            endif
            kb=0
            do i = 1,n
               if(node_param(i).eq.ipar) then
                  kb=kb+1
                  if(kb.eq.1) then
                     i1 = i
                     pnx1 = pnx(i)
                  else
                     if(pnx(i).ne.pnx1) then
                        if (iatty .ne. 0) then
                           write(iatty,*) 'perms not equal in same ',
     &                          'zone : stopping'
                           write(iatty,*) 'zone = ' ,ipar, 'nodes ',i1,i
                        end if
                        if (iout .ne. 0) then
                           write(iout,*) 'perms not equal in same ',
     &                          'zone : stopping'
                           write(iout,*) 'zone = ' ,ipar, 'nodes ',i1,i
                        end if
                     endif
                  endif
                  pnx(i)=pnx(i)*(1.d00+dfactor)*fperms(ii,1)
                  pny(i)=pnx(i)*fperms(ii,2)
                  pnz(i)=pnx(i)*fperms(ii,3)
               endif
            enddo
            call bnswer
            iterp =iterp +itert
            do i=1,n
               pho(i) = phi(i)
            enddo
            call headctr(2,0,0.0d00,0.0d00)
            do i=1,n
               headp(jj,i) = head(i)
            enddo
            if(mpestq.gt.0) then
               do i=1,mpestq
                  izone =npestq(i)
                  call pest_fluxz
     &                 (izone,sumfout, sumsink, sumsource, flxp(jj,i))
               enddo
            endif
         enddo
c     
c     form sensitivity
c     
         do i = 1,n
            bp(i) = 
     &           (headp(2,i)-headp(1,i))/(pnx(i)*dfactor)*perm_fac
         enddo
         do i = 1,mpestq
            flxp(1,i) = 
     &           (flxp(2,i)-flxp(1,i))/(pnx(i)*dfactor)*perm_fac
         enddo
c     
c     calculate flux targets if necessary
c     
 995     continue
         if (iatty .ne. 0) write(iatty,60) kb, iterp
         if (iout .ne. 0) write(iout,60) kb, iterp
         kb=0
 60      format ('nodes with this parameter ',
     &        i8,' solver iterations = ',i8)
         write(ispst2,10) (bp(npest(j)),j=1,mpest),
     &        (flxp(npestq(j),1),j=1,mpestq)
 10      format(5(1x,g21.14))
         iter=iter-iterp
         iterp = 0

      enddo
      deallocate (headp,flxp)
      close(ispst2)
c     might need to re-order the derivative file
c     so open file again
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
      read(ispst2,*) nparam1,mpest1
      allocate(dum2(nparam1,mpest1))
      do i=1,nparam1
         read(ispst2,*) (dum2(i,j),j=1,mpest1)
      enddo
      close(ispst2)
      open(ispst2,file=pest_deriv_file(1:lispst2),status='unknown')
      write(ispst2,*) nparam1,mpest1
      do j=1,mpest1 
         write(ispst2,772) (dum2(i,j),i=1,nparam1)
         write(ispst2,*) ' '
 772     format(5(1x,g21.14))
      enddo
      deallocate(dum2)
      close(ispst2)
      if (iatty .ne. 0) write(iatty,881) pest_deriv_file(1:lispst2)
      if (iout .ne. 0) write(iout,881) pest_deriv_file(1:lispst2)
 881  format
     &     (/,'>>> writing pest usable sensitivities, file: ',a30)
      close(ispst2)
      return
      end

      subroutine node_sensitivity
     &     (nparam,i,node_param,iparameter,deriv,fperms,ii)
C**********************************************************************
C     D1
C     D1 PURPOSE
C     D1
C     D1 To generate sensitivities for each node wrt a parameter   
C     D1 At present only coded for 1-phase saturated zone problems
C     D1 At present only coded for parameters that are permeabilities
C     D1
C**********************************************************************

      use comai
      use combi
      use comci
      use comdi
      use comgi
      use comji
      use davidi
      implicit none

      integer i
      integer ii
      integer icd
      integer ii1
      integer ii2
      integer idg
      integer iq
      integer jmi
      integer ij 
      integer ij1 
      integer ij2 
      integer jm
      integer neqp1
      integer iz
      integer kb
      integer neighc
      integer kz
      integer nmatavw
      real*8 sx1d
      real*8 axi
      real*8 ayi
      real*8 azi
      real*8 alxi
      real*8 alyi
      real*8 alzi
      real*8 avxi
      real*8 avyi
      real*8 avzi
      real*8 pvii
      real*8 phii
      real*8 swi
      real*8 dili
      real*8 dilkb
      real*8 divi
      real*8 divkb
      real*8 axkb
      real*8 aykb
      real*8 azkb
      real*8 alxkb
      real*8 alykb
      real*8 alzkb
      real*8 sx2c
      real*8 sx4d
c     real*8 sxzc
      real*8 pvikb
      real*8 phikb
      real*8 radi
      real*8 radkb
      real*8 fid
      real*8 fid1
      real*8 axyd
      real*8 axy
      real*8 axyf
      real*8 pvxy
      real*8 pxy
      real*8 pxyh
      real*8 pxyi
      real*8 sx4h
      real*8 vxyd
      real*8 vxy
      real*8 vxyf
      real*8 dpvti
      real*8 sx2t
      real*8 sx3t
      real*8 sxzt
      real*8 ti
      real*8 dis2
      real*8 delx2
      real*8 dely2
      real*8 delz2
      integer isl,nparam
      real*8 grav_air

      logical bit
      integer iz4m1
      integer imd,iwd
      integer node_param(*)
      integer iparameter
      real*8 derivii
      real*8 deriv(*)
      real*8 fperms(nparam,*)
      real*8 sx_tol
      real*8 dkxkbdkk,dkykbdkk,dkzkbdkk
      real*8 dkxidkk,dkyidkk,dkzidkk
      real*8 dkxdki,dkydki,dkzdki    
      real*8 dkxdk,dkydk,dkzdk    
      real*8 dkxdkkb,dkydkkb,dkzdkkb    
      real*8 perm_fac, dis_tol

      parameter (dis_tol=1.d-20,sx_tol=1.d-30,perm_fac=1.d06)
      

c     changed by avw -- entered here by seh
      neqp1=neq+1
      if(i.gt.neq) then
         nmatavw=ldna
      else
         nmatavw=0
      endif
      if(icons.le.maxit) then
         grav_air=0.0
      else
         grav_air=grav
      endif

c     
c     storage for upwind
c     

      sx1d=sx1(i)
      axi=pnx(i)
      ayi=pny(i)
      azi=pnz(i)
      alxi=axi
      avxi=axi
      alyi=ayi
      avyi=ayi
      alzi=azi
      avzi=azi
      pvii=phi(i)
      phii=pvii-pcp(i)
      dpvti=dpcef(i)
      dili=dil(i)
      divi=div(i)
      ti=t(i)
      swi=s(i)
c     
c     form constants for i>neq
c     
      if(i.gt.neq.and.idualp.eq.0) then
         icd=neq
      else
         icd=0
      endif
      iz=i-icd
c     
      iz4m1 = 4*(iz-1)+1
c     
      ii1=nelm(i-icd)+1
      ii2=nelm(i-icd+1)
      idg=nelmdg(i-icd)-ii1+1
      neqp1=neq+1
      iq=0
      jmi=nelmdg(i-icd)

      do 58 jm=ii1,ii2
         if(nelm(jm).ne.iz) then
            iq=iq+1
            it8(iq)=nelm(jm)+icd
            it9(iq)=iq
            it11(iq)=jm-neqp1
         endif
 58   continue
      iq=0
      do jm=ii1,jmi-1
         kb = nelm(jm)
         ij1=nelmdg(kb)+1
         ij2=nelm(kb+1)
         do  ij=ij1,ij2
            if(nelm(ij).eq.iz) then
               iq=iq+1
               it10(iq)=istrw(ij-neqp1)
               if(imdnode.ne.0) then
                  imd = mdnodes(i) + mdnodes(kb)
                  if(imd.lt.2) it10(iq) = -abs(it10(iq))
               endif
            endif
         enddo
      enddo
      do jm=jmi+1,ii2
         iq=iq+1
         it10(iq)=istrw(jm-neqp1)
      enddo
c     
c     determine if derivative exist fo node and neighbors
c     perm_fac accounts for perm scaling of 1.e6 (to correct viscosity)
c     
      if(node_param(i).eq.iparameter) then
         derivii  = 1.0*perm_fac
      else
         derivii  = 0.0
      endif
      do  jm=1,iq
         kb=it8(jm)
         if(node_param(kb).eq.iparameter) then
            deriv(jm) = 1.0*perm_fac
         else
            deriv(jm) = 0.0
         endif
      enddo
      dkxidkk = derivii*fperms(ii,1)
      dkyidkk = derivii*fperms(ii,2)
      dkzidkk = derivii*fperms(ii,3)
c     
c     3-d geometry
c     
      if(icnl.eq.0) then
         do 59 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            neighc=it9(jm)
            iwd=it10(jm)
            iw =abs(iwd)
            axkb=pnx(kb)
            aykb=pny(kb)
            azkb=pnz(kb)
            alxkb=axkb
            alykb=aykb
            alzkb=azkb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            perml(3)=2.*alzkb*alzi/(alzkb+alzi)
            sx2c=sx(iw,isox)+sx(iw,isoy)+sx(iw,isoz)
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            delz2=(cord(kz,3)-cord(iz,3))**2
            dis2=delx2+dely2+delz2
            if(dis2.gt.dis_tol.or.iwd.lt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2)+delz2/perml(3))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2),perml(3))
            endif
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            t1(neighc)=pxyi
            t2(neighc)=pxyh
            dkxkbdkk = deriv(jm)*fperms(ii,1)
            dkykbdkk = deriv(jm)*fperms(ii,2)
            dkzkbdkk = deriv(jm)*fperms(ii,3)
            if(abs(sx2c).ge.sx_tol) then
               dkxdkkb = 2.*((alxi/(alxkb+alxi)
     &              -alxkb*alxi/(alxkb+alxi)**2)*dkxkbdkk)
               dkxdki = 2.*((alxkb/(alxkb+alxi)
     &              -alxkb*alxi/(alxkb+alxi)**2)*dkxidkk)
               dkydkkb = 2.*((alyi/(alykb+alyi)
     &              -alykb*alyi/(alykb+alyi)**2)*dkykbdkk)
               dkydki = 2.*((alykb/(alykb+alyi)
     &              -alykb*alyi/(alykb+alyi)**2)*dkyidkk)
               dkzdkkb = 2.*((alzi/(alzkb+alzi)
     &              -alzkb*alzi/(alzkb+alzi)**2)*dkzkbdkk)
               dkzdki = 2.*((alzkb/(alzkb+alzi)
     &              -alzkb*alzi/(alzkb+alzi)**2)*dkzidkk)
               t3(neighc)=-sx2c*dis2/
     &              (delx2/perml(1)+dely2/perml(2)+delz2/perml(3))**2
     &              *(delx2/perml(1)**2*dkxdki+dely2/perml(2)**2*dkydki
     &              +delz2/perml(3)**2*dkzdki)          
               t4(neighc)=-sx2c*dis2/
     &              (delx2/perml(1)+dely2/perml(2)+delz2/perml(3))**2*
     &              (delx2/perml(1)**2*dkxdkkb+dely2/perml(2)**2*dkydkkb
     &              +delz2/perml(3)**2*dkzdkkb)          
            else
               t3(neighc)=0.0
               t4(neighc)=0.0
            endif
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 59      continue
         if(irdof.ne.11) then
c     
c     liquid phase calculations
c     
            do 60 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyi=t1(neighc)
               sx4d=t6(neighc)
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=axyd
 60         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
            isl=1
            do 61 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c     add coding to save upwind position
               if(iad.le.iad_up) then
                  fid=0.5
                  axyd=t8(neighc)
                  if(axyd.lt.0.0) fid=dnwgt
                  if(axyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c     
                  call setbit(nbits,neighc,upwind_l(iz4m1),fid)
c     
               else
                  if(bit(nbits,neighc,upwind_l(iz4m1))) then 
                     t9(neighc)=1.0
                  else
                     t9(neighc)=0.0
                  endif
               endif
 61         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 62 jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  axyd=t8(neighc)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  dilkb=dil(kb)
                  axyf=(fid*dilkb+fid1*dili)
                  axy=axyd*axyf
c     

c     bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
c     
c     new source term for analytic derivatives
c     
                  bp(iz+nrhs(1))=bp(iz+nrhs(1))+
     &                 (t3(jm)+t4(jm))*axyf*(phi(kb)-phii)
     &                 +0.5*(-grav)*(t3(jm)+t4(jm))*(rolf(i)+rolf(kb))
     &                 *(cord(kz,igrav)-cord(iz,igrav))

 62            continue
            endif
         endif
         if(irdof.ne.13) then   
c     
c     vapour phase calculations
c     
            do 63 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
 63         continue
c     
c     determine upwind nodes and if vapour phase exists
c     
            isl=1
            do 64 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c     add coding to save upwind position
               if(iad.le.iad_up) then
                  fid=0.5
                  vxyd=t8(neighc)
                  if(vxyd.lt.0.0) fid=dnwgt
                  if(vxyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c     
                  call setbit(nbits,neighc,upwind_v(iz4m1),fid)
c     
               else
                  if(bit(nbits,neighc,upwind_v(iz4m1))) then 
                     t9(neighc)=1.0
                  else
                     t9(neighc)=0.0
                  endif
               endif
 64         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 65 jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  vxyd=t8(neighc)
                  divkb=div(kb)
                  vxyf=(fid*divkb+fid1*divi)
                  vxy=vxyd*vxyf
c     

                  bp(iz+nrhs(2))=bp(iz+nrhs(2))+vxy
 65            continue
            endif
         endif
c     
c     
c     2-d geometry
c     
      elseif(icnl.ne.0) then
         radi=cord(iz,3)
         do 69 jm=1,iq
            kb=it8(jm)
            kz=kb-icd
            iw=it10(jm)
            neighc=it9(jm)
            axkb=pnx(kb)
            aykb=pny(kb)
            alxkb=axkb
            alykb=aykb
            perml(1)=2.*alxkb*alxi/(alxkb+alxi)
            perml(2)=2.*alykb*alyi/(alykb+alyi)
            radkb=0.5*(radi+cord(kz,3))
            sx2c=radkb*(sx(iw,isox)+sx(iw,isoy))
            pvikb=phi(kb)
            phikb=pvikb-pcp(kb)
            delx2=(cord(kz,1)-cord(iz,1))**2
            dely2=(cord(kz,2)-cord(iz,2))**2
            dis2=delx2+dely2
            if(dis2.gt.dis_tol.or.iwd.lt.0) then
               pxy=sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2))
            else
               pxy=sx2c*sx_mult*max(perml(1),perml(2))
            endif
            pxyi=pxy*(phikb-phii)
            pxyh=pxy*(pvikb-pvii)
            dkxkbdkk = deriv(jm)*fperms(ii,1)
            dkykbdkk = deriv(jm)*fperms(ii,2)
            if(abs(sx2c).ge.sx_tol) then
               dkxdkkb = 2.*((alxi/(alxkb+alxi)
     &              -alxkb*alxi/(alxkb+alxi)**2)*dkxkbdkk)
               dkxdki = 2.*((alxkb/(alxkb+alxi)
     &              -alxkb*alxi/(alxkb+alxi)**2)*dkxidkk)
               dkydkkb = 2.*((alyi/(alykb+alyi)
     &              -alykb*alyi/(alykb+alyi)**2)*dkykbdkk)
               dkydki = 2.*((alykb/(alykb+alyi)
     &              -alykb*alyi/(alykb+alyi)**2)*dkyidkk)
               dkxdk = 2.*((alxkb*dkxidkk+ alxi*dkxkbdkk)/(alxkb+alxi)
     &              -alxkb*alxi/((alxkb+alxi)**2)*(dkxidkk+dkxkbdkk))
               dkydk = 2.*((alykb*dkyidkk+ alyi*dkykbdkk)/(alykb+alyi)
     &              -alykb*alyi/((alykb+alyi)**2)*(dkyidkk+dkykbdkk))
               t3(neighc)=-sx2c*dis2/
     &              (delx2/perml(1)+dely2/perml(2))**2
     &              *(delx2/perml(1)**2*dkxdki+dely2/perml(2)**2*dkydki)
               t3(neighc)=-sx2c*dis2/
     &              (delx2/perml(1)+dely2/perml(2))**2
     &              *(delx2/perml(1)**2*dkxdk+dely2/perml(2)**2*dkydk)
               t4(neighc)=-sx2c*dis2/(delx2/perml(1)+
     &              dely2/perml(2))**2*(delx2/perml(1)**2
     &              *dkxdkkb+dely2/perml(2)**2*dkydkkb)
               t4(neighc)=0.0         
            else
               t3(neighc)=0.0
               t4(neighc)=0.0
            endif
            t6(neighc)=-grav*t3(neighc)
            t7(neighc)=-grav_air*t4(neighc)
 69      continue
         if(irdof.ne.11) then
c     
c     liquid phase calculations
c     
            do 70 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyi=t1(neighc)
               sx4d=t6(neighc)
               axyd=pxyi+0.5*sx4d*(rolf(i)+rolf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=axyd
 70         continue
c     
c     determine upwind nodes and if liquid phase exists
c     
            isl=1
            do 71 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c     add coding to save upwind position
               if(iad.le.iad_up) then
                  fid=0.5
                  axyd=t8(neighc)
                  if(axyd.lt.0.0) fid=dnwgt
                  if(axyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c     
                  call setbit(nbits,neighc,upwind_l(iz4m1),fid)
c     
               else
                  if(bit(nbits,neighc,upwind_l(iz4m1))) then 
                     t9(neighc)=1.0
                  else
                     t9(neighc)=0.0
                  endif
               endif
 71         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 72 jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  axyd=t8(neighc)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  dilkb=dil(kb)
                  axyf=(fid*dilkb+fid1*dili)
                  axy=axyd*axyf
c     
c     

c     bp(iz+nrhs(1))=bp(iz+nrhs(1))+axy
c     
c     new source term for analytic derivatives
c     
                  bp(iz+nrhs(1))=bp(iz+nrhs(1))+
     &                 (t3(jm)+t4(jm))*axyf*(phi(kb)-phii)
     &                 +0.5*(-grav)*(t3(jm)+t4(jm))*(rolf(i)+rolf(kb))
     &                 *(cord(kz,igrav)-cord(iz,igrav))

 72            continue
            endif
         endif
         if(irdof.ne.13) then
c     
c     vapour phase calculations
c     
            do 73 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
               pxyh=t2(neighc)
               sx4h=t7(neighc)
               vxyd=pxyh+0.5*sx4h*(rovf(i)+rovf(kb))
     *              *(cord(kz,igrav)-cord(iz,igrav))
               t8(neighc)=vxyd
 73         continue
c     
c     determine upwind nodes and if vapour phase exists
c     
            isl=1
            do 74 jm=1,iq
               kb=it8(jm)
               kz=kb-icd
               neighc=it9(jm)
c     add coding to save upwind position
               if(iad.le.iad_up) then
                  fid=0.5
                  vxyd=t8(neighc)
                  if(vxyd.lt.0.0) fid=dnwgt
                  if(vxyd.gt.0.0) fid=upwgt
                  t9(neighc)=fid
c     
                  call setbit(nbits,neighc,upwind_v(iz4m1),fid)
c     
               else
                  if(bit(nbits,neighc,upwind_v(iz4m1))) then 
                     t9(neighc)=1.0
                  else
                     t9(neighc)=0.0
                  endif
               endif
 74         continue
c     
c     form equations
c     
            if(isl.ne.0) then
               do 75 jm=1,iq
                  kb=it8(jm)
                  kz=kb-icd
                  neighc=it9(jm)
                  fid=t9(neighc)
                  fid1=1.0-fid
                  vxyd=t8(neighc)
                  divkb=div(kb)
                  vxyf=(fid*divkb+fid1*divi)
                  vxy=vxyd*vxyf
c     
                  bp(iz+nrhs(2))=bp(iz+nrhs(2))+vxy
 75            continue
            endif
         endif
c     
      endif
      bp(iz+nrhs(1))=bp(iz+nrhs(1))+sx1d*deni(i)
      if(irdof.ne.13) then
         bp(iz+nrhs(2))=bp(iz+nrhs(2))+sx1d*denei(i)
      endif
      
      return
      end

      subroutine pest_fluxz
     &     (izone,sumfout, sumsink, sumsource, sumboun)
c     
c     subroutine to calculate flux observations 
c     
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comei
      use comgi
      use comfi
      use comflow
      use comxi
      use davidi

      implicit none

      real*8 sumfout, sumsink, sumsource, sumboun
      integer i, iflg, inode, inneq, iroot, izone, k
      integer idummy, addnode, i1, i2, iconn, indexa_axy
      logical matrix_node


c     Compute fluxes out of zone (>0), leaving out flues
c     into other parts of the zone

c     Loop over all nodes
      do inode = 1, n0
c     Determine if node is fracture or matrix, set indexes
c     and flags accordingly
         if(inode.gt.neq) then
            inneq = inode-neq
            matrix_node = .true.
            addnode = nelm(neq+1)-neq-1
            idummy = neq
         else
            matrix_node = .false.
            inneq = inode
            addnode = 0
            idummy = 0
         end if
c     Determine if node is part of the zone being summed
         if(izoneflxz(inode).eq.izone) then
c     Add boundary condition sources
            sumboun=sumboun + sk(inode)
c     Set index for looping through a_axy depending on whether
c     the node is a fracture or matrix node
            i1 = nelm(inneq)+1
            i2 = nelm(inneq+1)
c     loop over every connecting node
            do iconn = i1, i2
               indexa_axy = iconn-neq-1+addnode
c     add to sum if it is flow out of the node
               if(a_axy(indexa_axy).gt.0.) then
c     add to sum only if the connecting node is not also
c     in the zone or else the connecting node is itself, i.e.
c     the value is a sink term
                  if(izoneflxz(idummy+nelm(iconn))
     2                 .ne.izone.or.nelm(iconn)
     3                 .eq.inneq) then
                     sumfout = sumfout + a_axy(indexa_axy)
                     if(nelm(iconn).eq.inneq) then
                        sumsink = sumsink +
     2                       a_axy(indexa_axy)
                     end if
                  end if
               elseif(a_axy(indexa_axy).lt.0.) then
c     add to source sum
                  if(nelm(iconn).eq.inneq) then
                     sumsource = sumsource +
     2                    a_axy(indexa_axy)
                  end if
               end if
            end do
         end if
      end do

      return
      end
