      subroutine write_avs_node_s(icall,
     2     neq,
     3     nscalar,
     4     lu,
     5     ifdual,
     6     iriver2)
!***********************************************************************
!  Copyright, 1993, 2004,  The  Regents of the University of California.
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
CD1 PURPOSE
CD1
CD1 Output AVS scalar node information for FEHM
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 10-SEP-93    Carl Gable     22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/write_avs_node_s.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:44   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:16   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:36:40   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.2   06 Jun 2001 08:26:26   pvcs
!D2 Update for extended dispersion tensor model
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:20 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.7   Wed Apr 09 09:38:50 1997   gaz
CD2 changes made to accomodate head output
CD2 
CD2    Rev 1.6   Thu Oct 24 15:11:52 1996   zvd
CD2 Increased output precision
CD2 
CD2    Rev 1.5   Mon Aug 05 11:23:06 1996   hend
CD2 Fix to Rev 1.4
CD2 
CD2    Rev 1.4   Mon Aug 05 09:56:02 1996   hend
CD2 Added check to print 1d-20 for smaller values
CD2 
CD2    Rev 1.3   Fri Feb 02 14:23:02 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.2   04/14/95 14:26:56   zvd
CD2 Corrected erroneous warning conditions for 1 scalar output
CD2 
CD2    Rev 1.1   01/20/95 13:31:38   tam
CD2 Changed format for strings from * to a56, kept length to 80 so left justified
CD2 
CD2    Rev 1.0   08/23/94 15:34:18   llt
CD2 Original version
CD2
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   None
CD3
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 None
CD4
CD4   
C*******************************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5   None
CD5
CD5 Local Types
CD5
CD5   None
CD5
CD5 Local variables
CD5
CD5   Identifier      Type     Description
CD5
CD5 Local Subprograms
CD5
CD5   None
CD5
C***********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6
CD6
CD6
C***********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7 None
CD7
C***********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.6 Provide Input/Output Data Files
CD9 3.0 INPUT AND OUTPUT REQUIREMENTS
CD9
C***********************************************************************
CDA
CDA REFERENCES
CDA
CDA None
CDA
C************************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN 
CPS Liquid Pressure   phi(i)-pcp(i)
CPS Temperature       t(i)
CPS Vapor Pressure    phi(i)
CPS Saturation        max(s(i),sattol)
CPS Hydraulic head    head(i) 
CPS Porosity          ps(i) 
CPS Source            sk(i) 
CPS Density (liquid)  rolf(i) 
CPS Density (vapor)   rovf(i)
CPS
CPS ERROR Checking
CPS error if nscalar.gt.maxscalar
CPS error if iopressure.ne.0  and  (ioliquid.eq.0).and.(iovapor.eq.0)
CPS error if (ioliquid.ne.0).or.(iovapor.ne.0) and iopressure.eq.0
CPS
CPS Output all 6 scalars
CPS   write phi(i)-pcp(i) t(i) phi(i) max(s(i),sattol) head(i) ps(i)
CPS Output 5 scalars
CPS   Lp, Vp, t, s, h 
CPS   Lp, Vp, t, s, ps 
CPS   Lp, Vp, t, h, ps 
CPS   Lp, t, s, h, ps 
CPS   Vp, t, s, h, ps 
CPS   Vp, Lp, s, h, ps 
CPS
CPS Output 3 scalars
CPS   if ioliquid = 0      - Vapor P, temperature, saturation
CPS   elseif iotemperature = 0 - Liquid P, Vapor P, saturation
CPS   elseif iovapor = 0       - Liquid P, temperature, saturation
CPS   elseif iosaturation = 0  - Liquid P, Vapor P, temperature
CPS   plus head combinations 
CPS   else error on 3 scalar output
CPS
CPS Output 2 scalars
CPS   if ioliquid, iotemperature = 0     - Vapor P, saturation
CPS   elseif iotemperature, iovapor  = 0     - Liquid P, saturation
CPS   elseif iosaturation, iovapor   = 0     - Liquid P, temperature
CPS   elseif ioliquid, iosaturation  = 0     - Vapor P, temperature
CPS   elseif iosaturation, iotemperature = 0 - Vapor P, Liquid P
CPS   elseif iovapor, ioliquid = 0           - temperature, saturation
CPS   plus head combinations 
CPS   else error on 2 scalar output
CPS
CPS Output 1 scalar
CPS   if ioliquid, iopressure   = 1          - Liquid P
CPS   elseif iotemperature          = 1      - temperature
CPS   elseif iovapor, iopressure    = 1      - Vapor P
CPS   elseif iosaturation    =1              - saturation
CPS   elseif iohead    =1                    - head
CPS   else error on 1 scalar output
CPS
CPS   
CPS END 
CPS 
C***********************************************************************
C---  Phil Adding Water Vapor Pressure to output for No Head Case
c---   pci coming into this subroutine from avs_io.f  
c----------------------------------------------------------------------

      use avsio
      use comai, only : altc, days, iadif, icnl, idof, nei_in,
     &     ns_in, phi_inc, istrs, neq_primary
      use combi, only : corz, izonef, nelm, nelmdg, sx1
      use comci, only : rolf, rovf
      use comdi
      use comfi, only : pci
      use comflow, only : a_axy, a_vxy
      use comwt, only : sattol, head_id, rlptol
      use comsi
      use davidi
c     RJP 1/12/07 added following
      use comco2
      use comriv, only : npoint_riv, nnelm_riv, nelm_riv
      implicit none

      integer maxscalar
      parameter (maxscalar = 32)
      integer neq,nscalar,lu,ifdual,icall,open_file,offset,iriver2
      integer i,j,iolp,iovp,nout,iz,iendz,il,idz, i1, i2, index, iaxy, k
      integer size_head, size_pcp, istart, iend, ic1, ic2, length, nadd
      integer icord1, icord2, icord3, ns_in0, irivp 
      integer nelm2(ns_in)
      real*8 hdum, sdum, px, py, pz, flxdum
      character*80 title(2*maxscalar+3)
      character*150 :: tecstring = ''
      character*150 :: tecstring_riv = ''
      character*500 tstring2
      character*500 string
      character*20 vstring
      character*43 tstring
      character*5 char_type
      character*3 dls

      save tecstring, tecstring_riv
C     BEGIN
      size_head = size(head)
      size_pcp = size(pcp)
C     nscalar=(iovapor+ioliquid)*iopressure+iosaturation+iotemperature
C     calculation done in avs_io()
      nout = nscalar + iocord
      ioLP = ioliquid*iopressure
      ioVP = iovapor*iopressure

C     ERROR checking:
      
      if(nscalar .gt. maxscalar)then
         write(lu,*)'--------------------------------------------'
         write(lu,*)'ERROR: WRITE_AVS_NODE_S'
         write(lu,*)'nscalar   = ',nscalar,' is greater than '
         write(lu,*)'maxscalar = ',maxscalar
         write(lu,*)'--------------------------------------------'
         return
      endif
      irivp = 0
      if(ifdual .ne. 0)then
         istart = neq + 1
         iend = neq * 2
         nadd = nelm(neq+1)-neq-1
         offset = neq
      else 
         if (iriver2 .ne. 0) then
c     Output for river/well nodes         
            istart = neq_primary + 1
            iend = neq_primary + npoint_riv
            nadd = 0
            offset = 0
            if(iriver2.eq.2) then
               irivp = 2
               ns_in0 = ns_in
               ns_in = 2
            endif
         else
            istart = 1
            iend = neq_primary
            nadd = 0
            offset = 0
            irivp = 0
         end if
      endif
      
      if (icall .eq. 1.and.irivp.eq.0) tecstring = ''
      if (icall .eq. 1.and.irivp.ne.0) tecstring_riv = ''

c     Surfer headers need to be written to each zone file
      if (altc(1:3) .ne. 'sur') then
         call write_avs_head_s(icall,nscalar,lu,ifdual,0,iriver2)
      end if

      if (iozone .ne. 0 ) then
         iendz = nsurf
      else 
         iendz = 1
         idz = iozone
      end if

      do iz = 1, iendz
c     Zone loop
         if (iozone .ne. 0) then
            idz = izone_surf(iz)
            if (altc(1:3) .eq. 'tec') then
               if (icall .gt. 1 .and. iocord .ne. 0) then
                  string = ''
                  if (icnl .eq. 0) then
                     if (iozid .eq. 0) then
                        write (string, 125) '1-3', iz
                     else
                        write (string, 125) '1-3,5', iz
                     end if
                  else
                     if (iozid .eq. 0) then
                        write (string, 125) '1-2', iz
                     else
                        write (string, 125) '1-2, 4', iz
                     end if
                  end if
                  if(irivp.eq.0) then
                     tecstring = trim(string)
                  else
                     tecstring_riv = trim(string)
                  endif
               end if
               write (lu, 118) trim(timec_string)
               if(irivp.eq.0) then
                  write (lu, 120) idz, trim(tecstring)
               else
                  write (lu, 120) idz, trim(tecstring_riv)
               endif
            end if
         else
            if (altc(1:3) .eq. 'tec') then
               string = ''
               if (icall .eq. 1 .and. iogeo .eq. 1) then
                  select case (ns_in)
               case (5,6,8)
                  write (string, 135) neq, nei_in, 'FEBRICK'
               case (4)
                  if (icnl .eq. 0) then
                     write (string, 135) neq, nei_in, 
     &                    'FETETRAHEDRON'
                  else
                     write (string, 135) neq, nei_in, 
     &                    'FEQUADRILATERAL'
                  end if
               case (3)
                  write (string, 135) neq, nei_in, 'FETRIANGLE'
               case (2)
                  if(irivp.eq.0) then
                     write (string, 135) neq, nei_in, 'FELINESEG'
                     ns_in=ns_in0
                  else if(irivp.eq.2)then                    
                     write (string, 135) npoint_riv, npoint_riv-1,
     &                    'FELINESEG'
                     ns_in=ns_in0
                  endif
               case (0)
c     fdm grid
                  write (string, '(a)') ''
               end select
               if(irivp.eq.0) then
                  tecstring = trim(string)
                  write (lu, 130) trim(timec_string), trim(tecstring)
               else
                  tecstring_riv = trim(string)
                  write (lu, 130) trim(timec_string),
     &                 trim(tecstring_riv)
               endif                  
               if (ns_in .eq. 0) then
                  if (iozid .eq. 0) then
                     write (string, 125) '1-3', iz
                  else
                     write (string, 125) '1-3, 5', iz
                  end if
               else if (icnl .eq. 0) then
                  if (iozid .eq. 0) then
                     write (string, 140) '1-3', iz
                  else
                     write (string, 140) '1-3, 5', iz
                  end if
               else
                  if (iozid .eq. 0) then
                     write (string, 140) '1-2', iz
                  else
                     write (string, 140) '1-2, 4', iz
                  end if
               end if
               if(irivp.eq.0) then
                  tecstring = trim(tecstring) // trim(string)
               else
                  tecstring_riv = trim(tecstring_riv) // trim(string)
               endif
            else if (icall .eq. 1 .and. iocord .ne. 0) then
               write (lu, 130) trim(timec_string)
               if (icnl .eq. 0) then
                  if (iozid .eq. 0) then
                     write (string, 125) '1-3', iz
                  else
                     write (string, 125) '1-3, 5', iz
                  end if
               else
                  if (iozid .eq. 0) then
                     write (string, 125) '1-2', iz
                  else
                     write (string, 125) '1-2, 4', iz
                  end if
               end if
               if(irivp.eq.0) then
                  tecstring = trim(string)
                  write (lu, 130) trim(timec_string), trim(tecstring) 
               else
                  tecstring_riv = trim(string)
                  write (lu, 130) trim(timec_string),
     &                 trim(tecstring_riv)
               endif                             
            else
               if(irivp.eq.0) then
                  write (lu, 130) trim(timec_string), trim(tecstring) 
               else
                  write (lu, 130) trim(timec_string),
     &                 trim(tecstring_riv)
               endif
            end if
         end if
      end if

      if (altc(1:3) .eq. 'sur') then
         call write_avs_head_s(icall,nscalar,lu,ifdual,idz,iriver2)
         dls = ', '
         k = 2
      else if (altc(1:4) .eq. 'avsx') then
         dls = ' : '
         k = 3
      else
         dls = ' '
         k = 1
      end if
      
      do i = istart, iend
                                ! Node loop          
         string = ''
         if (iozone .ne. 0) then
            if (izone_surf_nodes(i).ne.idz) goto 200
         end if
                                ! Node number will be written first for avs and sur files
         if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'sur') then
            if (ifdual .eq. 0) then
               write(string, 100) i
            else 
               write(string, 100) i - neq
            end if
            ic1 = 11
         else
            ic1 = 1
         end if
         if (iocord .ne. 0) then
            if (altc(1:4) .eq. 'avsx') then
c     For avsx all three coordinates are always output
               do j = 1,3
                  write(vstring,110) dls(1:k), corz(i - offset,j)
                  ic2 = ic1 + len_trim(vstring)
                  string(ic1:ic2) = vstring
                  ic1 = ic2 + 1
               end do
            else if (altc(1:4) .eq. 'avs' .or. altc(1:3) .eq. 'sur'
     &              .or. icall .eq. 1) then
c     Coordinates will only be output in the first file for tecplot
               select case (icnl)
            case (1, 4)
               icord1 = 1
               icord2 = 2
               icord3 = 1
            case (2, 5)
               icord1 = 1
               icord2 = 3
               icord3 = 2
            case(3, 6)
               icord1 = 1
               icord2 = 3
               icord3 = 1
            case default
               icord1 = 1
               icord2 = 3
               icord3 = 1
            end select
            do j = icord1, icord2, icord3
               write(vstring,110) dls(1:k), corz(i - offset,j)
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
            end do
         end if
      end if
c     Node numbers are written after coordinates for tec files
      if (altc(1:3) .eq. 'tec') then
         if (ifdual .eq. 0) then
            write(vstring, 105) dls(1:k), i
         else 
            write(vstring, 105) dls(1:k), i - neq
         end if
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iozid .eq. 1) then
         if (altc(1:4) .eq. 'avs' .or. altc(1:3) .eq. 'sur'
     &        .or. icall .eq. 1) then
            write(vstring, 115) dls(1:k), izonef(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
      end if
      if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
         if (size_pcp .ne. 1) then
            write(vstring,110) dls(1:k), phi(i)-pcp(i)
         else
            write(vstring,110) dls(1:k), phi(i)-phi_inc
         end if
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iopressure .eq. 1 .and. iovapor .eq. 1) then
         write(vstring,110) dls(1:k), phi(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         if (iadif .eq. 1) then
            write(vstring,110) dls(1:k), phi(i)-pci(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
      end if
      if (iocapillary .eq. 1) then
         write(vstring,110) dls(1:k), pcp(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iotemperature .eq. 1) then
         write(vstring,110) dls(1:k), t(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iosaturation .eq. 1) then
         if (ps(i) .le. 0.) then
            sdum = 0.d0
         else if (irdof .ne. 13 .or. ifree .ne. 0) then
c     sdum = max(s(i), sattol)
c     if (sdum .le. sattol) sdum = 0.d0
c     saturations are never zeroed out, report what is in array
            sdum = s(i)
         else
            sdum = 1.0d0
         end if
         write(vstring,110) dls(1:k), sdum
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (ioco2 .eq. 1) then
         if (ps(i) .le. 0.) then
            sdum = 0.d0
            write(vstring,112) dls(1:k), sdum
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
            write(vstring,112) dls(1:k), sdum
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
            write(vstring,112) dls(1:k), sdum
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
            write(vstring,112) dls(1:k), sdum
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
            write(vstring,115) dls(1:k), int(sdum)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         else
                                ! Water volume fraction
            write(vstring,112) dls(1:k), fw(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
                                ! Liquid co2 fraction
            write(vstring,112) dls(1:k), fl(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
                                ! Gaseous co2 fraction
            write(vstring,112) dls(1:k), fg(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
                                ! Dissolved co2 mass fraction
            write(vstring,112) dls(1:k), yc(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1                 
                                ! Phase state of co2
            write(vstring,115) dls(1:k), ices(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
      end if
      if (iohead .eq. 1 .and. size_head .ne. 1) then
         if (ps(i) .le. 0.) then
            hdum = 0.d0
         else
            hdum = head(i)
c     might need help in the 
            if (irdof .ne. 13 .or. ifree .ne. 0) then
               if (s(i).lt.sattol+rlptol) hdum = head_id
            endif 
         end if
         write(vstring,110) dls(1:k), hdum
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (ioporosity .eq. 1) then
         write(vstring,110) dls(1:k), ps(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iosource .eq. 1) then
         write(vstring,110) dls(1:k), sk(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
         write(vstring,110) dls(1:k), rolf(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iodensity .eq. 1 .and. iovapor .eq. 1) then
         write(vstring,110) dls(1:k), rovf(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      end if
      if (iopermeability .eq. 1) then
         if (idof .ne. 0) then
            px=log10(pnx(i)*1.d-06)
            py=log10(pny(i)*1.d-06)
            pz=log10(pnz(i)*1.d-06)
         else
            px = 0.
            py = 0.
            pz = 0.
         end if
         write(vstring,110) dls(1:k), px
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         write(vstring,110) dls(1:k), py
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         write(vstring,110) dls(1:k), pz
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      endif
      if (ioflx .eq. 1 .and. ioliquid .eq. 1) then
         if (.not. net_flux) then
            iaxy = nelmdg (i) - (neq + 1) + nadd
            if (vol_flux) then
               if (sx1(i) .gt. 0.) then
                  flxdum = a_axy(iaxy) / sx1(i)
               else
                  flxdum = 0.d0
               end if
            else
               flxdum = a_axy(iaxy)
            end if
            write(vstring,110) dls(1:k), flxdum
         else
            i1 = nelm(i) + 1
            i2 = nelm(i+1)
            flxdum = 0.
            do index = i1, i2
               iaxy = index - neq - 1 + nadd
               if(a_axy(iaxy).gt.0.) then
                  flxdum = flxdum + a_axy(iaxy)
               end if
            end do
            if (vol_flux) then
               if (sx1(i) .gt. 0.) then
                  flxdum = flxdum / sx1(i)
               else
                  flxdum = 0.d0
               end if
            end if
            write(vstring,110) dls(1:k), flxdum
         end if
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1                  
      end if
      if (ioflx .eq. 1 .and. iovapor .eq. 1) then
         iaxy = nelmdg (i) - (neq + 1) + nadd
         write(vstring,110) dls(1:k), a_vxy(iaxy)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1                  
      end if
      if (iodisp .eq. 1.and.idisp_rel.ne.0) then
         write(vstring,110) dls(1:k), du(i)-du_ini(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         write(vstring,110) dls(1:k), dv(i)-dv_ini(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         if (icnl .eq. 0) then
            write(vstring,110) dls(1:k), dw(i)-dw_ini(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
      else if (iodisp .eq. 1) then
         write(vstring,110) dls(1:k), du(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         write(vstring,110) dls(1:k), dv(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         if (icnl .eq. 0) then
            write(vstring,110) dls(1:k), dw(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
      endif 
      if (iostress .ne. 0) then
         write(vstring,110) dls(1:k), str_x(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         write(vstring,110) dls(1:k), str_y(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
         if(icnl.eq.0) then
            write(vstring,110) dls(1:k), str_z(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         end if
         write(vstring,110) dls(1:k), str_xy(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1   
         if(icnl.eq.0) then
            write(vstring,110) dls(1:k), str_xz(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
            write(vstring,110) dls(1:k), str_yz(i)
            ic2 = ic1 + len_trim(vstring)
            string(ic1:ic2) = vstring
            ic1 = ic2 + 1
         endif            
      endif 
      if (iostrain .eq. 1) then
         write(vstring,110) dls(1:k), vol_strain(i)
         ic2 = ic1 + len_trim(vstring)
         string(ic1:ic2) = vstring
         ic1 = ic2 + 1
      endif  
      length = len_trim(string)
      write(lu,'(a)') string(1:length)           
 200  enddo
      call flush(lu)
      if (altc(1:3) .eq. 'sur') close (lu)
      enddo

      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. iogeo .eq. 1)
     &     then
c     Read the element connectivity and write to tec file
         if(irivp.eq.0) then
            il = open_file(geoname,'old')
! avsx geometry file has an initial line that starts with neq_primary
            read(il,*) i
            if (i .ne. neq_primary) backspace il
            do i = 1, neq
               read(il,*)
            end do
            do i = 1, nei_in
               read (il,*) i1,i2,char_type,(nelm2(j), j=1,ns_in)
               write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
            end do
            close (il)
         else
c     river segments (2 node elements)
            do i = 1,nnelm_riv
               write(lu,'(2(i8))') nelm_riv(i,1)-neq,
     &              nelm_riv(i,2)-neq
            enddo
         endif
      end if
      if (altc(1:3) .ne. 'sur') close (lu)

 100  format(i10.10)
 105  format(a, i10.10)
 110  format(a, g16.9)
 112  format(a, f10.4)
 115  format(a, i4)
c     120  format('ZONE T = "',i4.4,' Simulation time ',1p,g16.9,' days"', a)
c     118  format('TEXT T = "Simulation time ',1p,g16.9,' days"')
 118  format('TEXT T = ', a)
 120  format('ZONE T = "',i4.4, '"', a)
 125  format(', VARSHARELIST = ([', a,'] = ', i4, ')')
c     130  format('ZONE T = "Simulation time ',1p,g16.9,' days"', a)
 130  format('ZONE T =', a, a) 
 135  format(', N = ', i8, ', E = ', i8, ', DATAPACKING = POINT',
     &     ', ZONETYPE = ', a)
 140  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')

      return
      end

