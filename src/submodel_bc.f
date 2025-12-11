      subroutine submodel_bc(iflg)
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
!!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To create new "flow" macro to represent boundary conditions on 
!D1 extracted submodel.
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 9-Jan-02, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/submodel_bc.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!D3 2.3.7 Sources and sinks
!D3 2.6   Provide Input/Output Data Files
!D3 3.0   INPUT AND OUTPUT REQUIREMENTS
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************

      use comflow
      use davidi
      use comji
      use comfi
      use comgi
      use comxi
      use comei
      use comdi
      use comii
      use comci
      use combi
      use comdti
      use comki  
      use comai

      implicit none

      integer, allocatable :: kq_dum(:), icount(:)
      integer i,j,ii,jj,kb,i1,i2,neqp1,max_subboun
      integer izone1,izone2,iflg,ibnd,iroot,idsubm 
      integer idsubmc,isubmd,isubmodel0,open_file   
      integer mi,ik,ityps,itemp,ic
      integer nmatavw,izik
c gaz 110819 removed tref, pref (now global)       
      real*8 subflux,aiped, flux_gh,  pres_gh
      real*8 head_value, parm1, parm2, parm3
      logical null1
      character*4 keyword
      character*9 temp_name
      character*4 dsubm
      character*100 submod_root
      character*120 submod_name
      save keyword,izone1,izone2
      parameter(aiped=1.d02,max_subboun = 50)
      parameter(temp_name = 'subm_temp')
c

      save submod_root, iroot, icount
      if(isubbc.eq.2) then
         if(iflg.eq.0) then
c
c create filenamne
c  
            if (.not. allocated(submodfile)) then
c 
c this should be the fist time called
c
               allocate (submodfile(max_subboun))
               allocate (isubmodelfile(max_subboun))
               allocate (isubmodnamlen(max_subboun))
               allocate (submod_filename(max_subboun))
               allocate (izonesub1(max_subboun))
               allocate (itypsd(max_subboun))         
               allocate (keyms1(max_subboun),keyms2(max_subboun))
               allocate (keyms3(max_subboun),keyms4(max_subboun))
               allocate (iflux_list(n0))
               allocate (icount(max_subboun))
               icount = 0
               isubmodel = 0
               submodfile = 0
               isubmodelfile = 0
               isubmodnamlen = 0
               submod_filename = ''
               if (null1(root_name)) then
! Use  file root name
                  if (nmfil(9) .ne. nmfily(3) .and. nmfil(9) .ne. ' ') 
     &                 then
                     call file_prefix(nmfil(9), iroot)
                     if (iroot .gt. 100) iroot = 100
                     submod_root(1:iroot) = nmfil(9)(1:iroot)
                  else
                     if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) 
     &                    .ne. ' ') then
                        call file_prefix(nmfil(5), iroot)
                        if (iroot .gt. 100) iroot = 100
                        submod_root(1:iroot) = nmfil(5)(1:iroot)
                     else
                        if (nmfil(2)(1:1) .eq. ' ' ) then
                           write (ierr, *) 'FILE ERROR: nmfil2 file: ', 
     &                          nmfil(2), ' unable to ',
     &                          'determine submod file prefix'
                           stop
                        else
                           call file_prefix(nmfil(2), iroot)
                           if (iroot .gt. 100) iroot = 100
                           submod_root(1:iroot) = nmfil(2)(1:iroot)
                        end if
                     end if
                  endif
               else
                  iroot = len_trim (root_name)
                  if (iroot .gt. 100) iroot = 100
                  submod_root(1:iroot) = root_name(1:iroot)
               end if        
            endif
c
c          
c        read input parameters
c
            isubmd = 0
            isubmodel0 = isubmodel
            itypsd = 0
            keyms1 = ''
            keyms2 = ''
            keyms3 = ''
            keyms4 = ''

            do
               read(inpt,'(a80)') wdd1
               if(null1(wdd1)) exit
c create a file for each submodel type         
               isubmd = isubmd + 1
               isubmodel = isubmodel + 1
               write(dsubm,'(i4.4)') isubmodel
               submod_name = ''
               submod_name(1:iroot) = submod_root(1:iroot)
               submod_name(iroot+1:iroot+1) ='.'
               submod_name(iroot+2:iroot+5) = dsubm(1:4)
               submod_name(iroot+6:iroot+11) = '.wflow'
               isubmodnamlen(isubmodel)= iroot+11
               submod_filename(isubmodel) = 
     &              submod_name(1:isubmodnamlen(isubmodel))
c         
c complete name here
c
               isubmodelfile(isubmodel) = 
     &              open_file(submod_filename(isubmodel), 'unknown')
               write(isubmodelfile(isubmodel),*) ' '
               read(wdd1,*) keyms1(isubmd),keyms2(isubmd),
     &              keyms3(isubmd),keyms4(isubmd)
               if(keyms4(isubmd).eq.'type') then
                  read(wdd1,*) keyms1(isubmd),keyms2(isubmd),
     &                 keyms3(isubmd),keyms4(isubmd),itypsd(isubmd)
               endif
            end do
            
	    allocate (kq_dum(n0))
            kq_dum = 0
            igroup = 1
            narrays = 1
            itype(1) = 4
            default(1) = 0
            igroup = 1
            call initdata2( inpt, ischk, n0, narrays,
     &           itype, default, macroread(8), macro, igroup, ireturn,
     &           i4_1 = kq_dum(1:n0)) 
c 
c  write out information to the opened file 
c
            do i = 1, n0
               isubmd = kq_dum(i)  
               if(isubmd.ne.0) then
                  j = isubmodel0 + isubmd
c icount > 0 will indicate that this model has been used
                  icount(j)  = icount(j) + 1
                  if(keyms4(isubmd).ne.'type') then
                     write(isubmodelfile(j),301) i,isubmd,keyms1(isubmd)
     &                   ,keyms2(isubmd),keyms3(isubmd),keyms4(isubmd)
                  else
                     write(isubmodelfile(j),301) i,isubmd,
     &                    keyms1(isubmd),keyms2(isubmd),keyms3(isubmd),
     &                    keyms4(isubmd),itypsd(isubmd)
                  endif
               endif
            enddo
            do j = isubmodel0+1,isubmodel 
               write(isubmodelfile(j),'(a4)') 'end '
               write(isubmodelfile(j),'(8(1x,i9))') (izonef(ik),ik=1,n0)
               close (isubmodelfile(j))
            enddo
 301        format(1x,i9,1x,i4,4(1x,a5),1x,2i5)       
c close file
            deallocate(kq_dum)

         else if(iflg.eq.2) then
c
c read data type and printout new flow macros
c 
            neqp1 = neq+1
            nmatavw=ldna
            isubmd = 1 
            do mi = 1, isubmodel
               isubmodelfile(mi) = open_file(submod_filename(mi), 'old')
               if (icount(mi) .eq. 0) then
c The model wasn't used, close and delete file
                  close (isubmodelfile(mi), status = 'DELETE')
                  write (ierr, 320) mi
                  if (iout .ne. 0) write (iout, 320) mi
                  if (iptty .ne. 0) write (iptty, 320) mi
 320              format (/, ' WARNING : wflow model ', i2, 
     &                 ' was not assigned to any nodes - file deleted') 
               else
c Can't use file identifiers from before
c               open(isubmodelfile(mi), 
c     &              file = subm_name(1:j),status = 'unknown')  
c
c read zonefile
c   
                  ityps = 0
                  j = 0
                  read(isubmodelfile(mi),'(a4)') keyword
                  read(isubmodelfile(mi),301)  i,isubmd,keyms1(1),
     &                 keyms2(1),keyms3(1),keyms4(1)
                  if(keyms4(1).eq.'type') then
                     read(isubmodelfile(mi),301) i,isubmd,keyms1(1),
     &                    keyms2(1),keyms3(1),keyms4(1),ityps
                  endif
                  rewind isubmodelfile(mi)
                  do
                     read(isubmodelfile(mi),'(a4)') keyword
                     if(keyword.eq.'end ') then
                        read (isubmodelfile(mi),'(8(1x,i9))') 
     &                       (izonef(ik),ik=1,n0)
                        exit
                     endif
                  enddo
                  rewind isubmodelfile(mi)     
                  itemp = open_file(temp_name, 'unknown') 
                  read(isubmodelfile(mi),'(a4)') keyword
                  if(ityps.eq.0) then
                     write(itemp,'(a4)')'flow'
                  else if(ityps.ne.0) then
                     write(itemp,'(a4)')'flo3'
                  endif
                  ic = 1
                  do       
                     read(isubmodelfile(mi),'(a4)') keyword
                     if(keyword.eq.'end ') then
                        exit
                     else
                        backspace isubmodelfile(mi)
                        if(ityps.eq.0) then
                           read(isubmodelfile(mi),301) ik,isubmd,
     &                          keyms1(isubmd),keyms2(isubmd),
     &                          keyms3(isubmd),keyms4(isubmd)
                        else
                           read(isubmodelfile(mi),301) ik,isubmd,
     &                          keyms1(isubmd),keyms2(isubmd),
     &                          keyms3(isubmd),keyms4(isubmd),ityps
                        endif
c     Key words may be identified by a unique subset of the string
                        if(keyms1(isubmd)(1:1).eq.'p' .or. 
     &                       keyms1(isubmd)(1:1).eq.'P') then
c     Pressure
                           parm1 = phi(ik)
                        else if(keyms1(isubmd)(1:1).eq.'h' .or.
     &                          keyms1(isubmd)(1:1).eq.'H') then
c     Head
                           call headctr(4,ik,phi(ik),head_value)
                           parm1 = head_value
                        else if(keyms1(isubmd)(1:1).eq.'f' .or.
     &                          keyms1(isubmd)(1:1).eq.'F') then
c     Flux
                           if(ka(ik).eq.0) then
                              izik = izonef(ik)
                              i1 = nelm(ik)+1
                              i2 = nelm(ik+1)
                              parm1 = 0.0
                              do jj = i1,i2
                                 kb = nelm(jj)
                                 if(izonef(kb).ne.izik) then
                                    parm1 = parm1 + 
     &                                   a_axy(jj-neqp1+nmatavw)
                                 endif
                              enddo
                           else
                              parm1 = sk(ik)
                           endif
                        endif
                        if(keyms2(isubmd)(1:1).eq.'s' .or.
     &                       keyms2(isubmd)(1:1).eq.'S') then
c     Saturation
                           if(ifree.eq.0.and.irdof.eq.13) then
                              parm2 = 1.0
                           else
                              parm2 = s(ik)
                           endif
                        else if(keyms2(isubmd)(1:1).eq.'t' .or.
     &                          keyms2(isubmd)(1:1).eq.'T') then
c     Temperature
                           parm2 = -t(ik)
                        else if(keyms2(isubmd)(1:1).eq.'e' .or.
     &                          keyms2(isubmd)(1:1).eq.'E') then
c     Enthalpy
                           parm2 = enlf(ik)
                        else if(keyms2(isubmd)(1:1).eq.'w' .or.
     &                          keyms2(isubmd)(1:1).eq.'W') then
c     Water only source
                              parm2 = 1.0
                        else if(keyms2(isubmd)(1:1).eq.'a' .or.
     &                          keyms2(isubmd)(1:1).eq.'A') then
c     Water only source
                              parm2 = -1.0                          
                        else
                           write(ierr, 350) 2, trim(keyms2(isubmd)), mi
                           if (iout .ne. 0) write(iout, 350) 2, 
     &                          trim(keyms2(isubmd)), mi
                           if (iptty .ne. 0) write(iptty, 350) 2, 
     &                          trim(keyms2(isubmd)), mi
c The model wasn't used, close and delete file
                           close (isubmodelfile(mi), status = 'DELETE')
                           go to 399
                        endif
c     Impedance key words need 4 characters, to identify type 
                        if (keyms1(isubmd).eq.'flux') then
                           parm3 = 0.e0
                        else if(keyms3(isubmd).eq.'imph') then
                           parm3 = 1.e00
                        else if(keyms3(isubmd).eq.'impl') then
                           parm3 = 1.e-4
                        else if(keyms3(isubmd).eq.'impn') then
                           parm3 = -1.e00
                        else
                           parm3 = 1.e02
                        endif    
                        ic = ic+1
                        write(itemp,401)ik,ik,1,parm1,parm2,parm3,ityps
                     endif
                  enddo
 350              format (/, ' WARNING : Invalid Keyword ', i1, ' "', 
     &                 a, '", output not written for wflow model ', i2)
                  write(itemp,*)
                  ic = ic+1
                  rewind isubmodelfile(mi)
                  rewind itemp
                  do i = 1,ic
                     read(itemp,'(a80)') wdd1(1:80)
                     write(isubmodelfile(mi),'(a80)') wdd1(1:80)
                  enddo
                  close(isubmodelfile(mi))
 399              close(itemp, status = 'DELETE')
               endif
            enddo 
         endif
 401     format(2(1x,i8),1x,i3,1p,1x,g15.7,0p,f10.4,1p,g12.4,i5)
         return
      endif
c tradional submodel code      
      if(iflg.eq.0) then
c     
c     set up output files and unit numbers
c     
c     Assign file unit numbers and determine output filenames
         isubm  = nufilb(24)
         if (null1(root_name)) then
            if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
c     Prefix from output file name
               call file_prefix(nmfil(5), iroot)
               if (iroot .gt. 94) iroot = 94
               nmfil(24) = nmfil(5)(1:iroot) // suffix(24)
            else 
c     Use default filenames "fehmn.subbc"
               if (nmfil(24)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil24 file: ',
     .                 nmfil(24), 
     .                 ' unable to determine submodel file name'
                  stop
               end if
            endif
         else
            iroot = len_trim (root_name)
            if (iroot .gt. 94) iroot = 94
            nmfil(24) = root_name(1:iroot) // suffix(24)
         end if
         open(isubm, file = nmfil(24), status = cstats(24))

c     read input parameters
         read(inpt,'(a80)') wdd1
         izone1 = 0
         izone2 = 0
         read(wdd1,*,end=10) keyword,izone1,izone2
         go to 20
 10      izone2=0
 20      if(keyword.ne.'flux'.and.keyword.ne.'head'
     &        .and.keyword.ne.'pres'.and.keyword.ne.'flow'.and.
     &        keyword.ne.'init'.and.keyword.ne.'flgh') then
            write (ierr, *) '>>>> error in keyword for macro subm <<<<'
            if (iout .ne. 0) then
               write(iout,*) 
               write(iout,*) '>>>> error in keyword for macro subm <<<<'
            end if
            if(iptty.ne.0) then
               write(iptty,*)'>>>> error in keyword for macro subm <<<<'
            endif
            stop
         endif
         if(keyword.eq.'init') then
            write(isubm, 1000) 'pres', verno, jdate, jtime, icnl
         else if(keyword.eq.'flgh') then
            write(isubm, 1000) 'flow', verno, jdate, jtime, icnl
         else 
            write(isubm, 1000) 'flow', verno, jdate, jtime, icnl
         endif
 1000    format (a4,'  Submodel BCs ', a30, 3x, a11, 3x, a8, i4)

      else if(iflg.eq.1) then
c     
c     save zonefile for submodel
c     
         allocate(izonesubm(neq))
         izonesubm = 0
         do i=1,neq
            if(izonef(i).eq.izone1) izonesubm(i)=izone1
            if(izone2.ne.0.and.izonef(i).eq.izone2) izonesubm(i)=izone2
         enddo
      else if(iflg.eq.2) then
c     
c     create flow macro for submodel
c     
         neqp1 = neq+1
	 if(keyword.eq.'flgh') then
c     print out all generalized head bc
c gaz 110819 tref now global       
c            tref = crl(6,1)
            do ii = 1,ngh 
               i = node_gh(ii)
               call headctr(5,i,pres_gh,pflow_gh(ii))   
	       flux_gh = wellim_gh(ii)*(phi(i)-pres_gh)
               write(isubm,1998) 
     &              i,i,flux_gh,tref,idir_gh(ii),cord(i,1),cord(i,2),
     &              cord(i,3),wellim_gh(ii)         
            enddo
 1998       format(2i10,' 1 ',1x,1pg18.9,1x,g9.3,' 1. ',i3,'  ',
     &           3(g12.6,1x),' k*den/vis*(A/d) ',g10.4)
            write(isubm,*)
            write(isubm,'(a4)') 'stop'

         else if(keyword.eq.'init') then
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
                  if(ihead.ne.0) then
                     write(isubm,1999) i,i,head(i),t(i),
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  else 
                     write(isubm,1999) i,i,phi(i),t(i),
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  end if
               end if
            enddo
 1999       format(2i10,' 1 ',1x,1pg18.9,1x,g9.2,' 1 # ',
     &           3(g12.6,1x),'z ',i5)
            write(isubm,*)
            write(isubm,'(a4)') 'stop'
         elseif(keyword.eq.'flow') then
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
                  if(ihead.ne.0) then
                     write(isubm,2000) i,i,head(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  else 
                     write(isubm,2000) i,i,phi(i),
     &                    aiped,cord(i,1),cord(i,2),cord(i,3),izone1
                  end if
               end if
            enddo
         else if(keyword.ne.'flux') then
c     set pressure boundary conditions
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
                  i1=nelm(i)+1
                  i2=nelm(i+1)
                  do ii=i1,i2
                     kb=nelm(ii)
                     if(izonesubm(kb).ne.izone1) then
                        go to 30
                     endif
                  enddo
                  go to 40
 30               continue
                  if(keyword.eq.'head') then
                     write(isubm,2000) i,i,head(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  else if(keyword.eq.'pres') then
                    if(ico2.lt.0) then
                     write(isubm,2000) i,i,phi(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                    else
                     if(itsat.le.10) then
                      write(isubm,2001) i,i,phi(i),-t(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1  
                     else if(itsat.gt.10) then
                       write(isubm,2002) i,i,phi(i),t(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1      
                     endif
                    endif
                  end if
 2000             format(2i10,' 1 ',1x,1pg18.9,' 1. ',g9.2,'  # ',
     &                 3(g12.6,x), 'z ',i5)
 2001             format(2i10,' 1 ',1x,1pg18.9,1x,g11.4,1x,g9.2,' # ',
     &                 3(g12.6,x), 'z ',i5)    
 2002             format(2i10,' 1 ',1x,1pg18.9,1x,g11.4,1x,g9.2,' # ',
     &                 3(g12.6,x), 'z ',i5,' ice')         
               endif
 40            continue
            enddo
            write(isubm,*)
            write(isubm,'(a4)') 'stop'
         else if (keyword.eq.'flux') then
c zvd 16-Mar-09 remove sum over a_axy, just want boundary source/sink
c     and don't need to consider if neighbor is in excluded zone
c         else if(izone2.ne.0) then
c     code when submodel zone and other domain is specified
c            do i=1,neq
c               if(izonesubm(i).eq.izone1) then
c     first add existing flux(recharge)
c     gaz    subflux = sk(i)
c                  subflux = 0.0d00
c                  ibnd = 0
c                  i1=nelm(i)+1
c                  i2=nelm(i+1)
c                  do ii=i1,i2
c                     kb=nelm(ii)
c                     if(izonesubm(kb).eq.izone2) then
c                        subflux = subflux+a_axy(ii-neqp1)
c                        ibnd = 1
c                     endif
c                  enddo
c                  if(ibnd.ne.0) then
c                     write(isubm,2000) i,i,subflux,0.0d00,
c     &                    cord(i,1),cord(i,2),cord(i,3), izone1
c                  endif
c               endif
c            enddo
c            write(isubm,*)
c            write(isubm,'(a4)') 'stop'
c         else if(izone2.eq.0) then
c     code when only submodel zone is specified
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
c     first add existing flux(recharge)
c                  ibnd = 0
                  ibnd = 1
                  subflux = sk(i)
c                  i1=nelm(i)+1
c                  i2=nelm(i+1)
c                  do ii=i1,i2
c                     kb=nelm(ii)
c                     if(izonesubm(kb).ne.izone1) then
c                        subflux = subflux+a_axy(ii-neqp1)
c                        ibnd = 1
c                     endif
c                  enddo
                  if(ibnd.ne.0) then
                     write(isubm,2000) i,i,subflux,0.0d00,
     &                    cord(i,1),cord(i,2),cord(i,3), izone1
                  endif
               endif
            enddo
            write(isubm,*)
            write(isubm,'(a4)') 'stop'
         endif
c     end of major if block
         deallocate(izonesubm)
      endif
      
      return
      end
