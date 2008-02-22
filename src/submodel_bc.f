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
!D2    Rev 2.5   06 Jan 2004 10:44:16   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:20:28   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
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
      use comai

      implicit none

      integer i,j,ii,jj,kb,i1,i2,neqp1
      integer izone1,izone2,iflg,ibnd,iroot     
      real*8 subflux,aiped, flux_gh, tref, pres_gh
      character*4 keyword
      logical null1
      save keyword,izone1,izone2
      parameter(aiped=1.d02)
c     
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
            tref = crl(6,1)
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
                     write(isubm,1999) i,i,head(i),20.,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  else 
                     write(isubm,1999) i,i,phi(i),20.,
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
                     write(isubm,2000) i,i,phi(i),aiped,
     &                    cord(i,1),cord(i,2),cord(i,3),izone1
                  end if
 2000             format(2i10,' 1 ',1x,1pg18.9,' 1. ',g9.2,'  # ',
     &                 3(g12.6,x), 'z ',i5)
               endif
 40            continue
            enddo
            write(isubm,*)
            write(isubm,'(a4)') 'stop'
         else if(izone2.ne.0) then
c     code when submodel zone and other domain is specified
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
c     first add existing flux(recharge)
c     gaz    subflux = sk(i)
                  subflux = 0.0d00
                  ibnd = 0
                  i1=nelm(i)+1
                  i2=nelm(i+1)
                  do ii=i1,i2
                     kb=nelm(ii)
                     if(izonesubm(kb).eq.izone2) then
                        subflux = subflux+a_axy(ii-neqp1)
                        ibnd = 1
                     endif
                  enddo
                  if(ibnd.ne.0) then
                     write(isubm,2000) i,i,subflux,0.0d00,
     &                    cord(i,1),cord(i,2),cord(i,3)
                  endif
               endif
            enddo
            write(isubm,*)
            write(isubm,'(a4)') 'stop'
         else if(izone2.eq.0) then
c     code when only submodel zone is specified
            do i=1,neq
               if(izonesubm(i).eq.izone1) then
c     first add existing flux(recharge)
                  ibnd = 0
                  subflux = sk(i)
                  i1=nelm(i)+1
                  i2=nelm(i+1)
                  do ii=i1,i2
                     kb=nelm(ii)
                     if(izonesubm(kb).ne.izone1) then
                        subflux = subflux+a_axy(ii-neqp1)
                        ibnd = 1
                     endif
                  enddo
                  if(ibnd.ne.0) then
                     write(isubm,2000) i,i,subflux,0.0d00,
     &                    cord(i,1),cord(i,2),cord(i,3)
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
