      subroutine write_avs_node_mat_s(icall,
     2     neq,
     3     nscalar,
     4     lu,
     5     ifdual,
     6     iriver2, mout)
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
!D2 
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
      use comai, only : altc, days, grav, iadif, icnl, ico2, idof, 
     &     ichead, ihead, nei_in, ns_in, phi_inc, istrs, ivf,
     &     neq_primary, rho1grav, ifdm_elem, igrav, ns, gdkm_flag,
     &     iccen, igrav, wdd, jdate, jtime, verno, ice, rlp_flag, idoff 
      use combi, only: corz, izonef, nelm, nelmdg, sx1, ncord, 
     &     ncord_inv, elem_temp, elem_geo, nact_elem_geo
      use comci, only: rolf, rovf
      use comdi
      use comfi, only : pci
      use comflow, only : a_axy, a_vxy
      use comfem
      use comii, only : crl 
      use comwt, only : sattol, head_id, rlptol
      use comsi
      use davidi
c     RJP 1/12/07 added following
      use comco2
      use comriv, only : npoint_riv, nnelm_riv, nelm_riv, iriver
      use comrxni, only : rxn_flag
      use comchem, only: cpntnam
      implicit none

      integer maxscalar
      parameter (maxscalar = 37)
      integer neq,nscalar,lu,ifdual,icall,open_file,offset,iriver2
      integer i,j,iolp,iovp,nout,iz,iendz,il,idz, i1, i2, index, iaxy, k
      integer size_head, size_pcp, istart, iend, ic1, ic2, length, nadd
      integer icord1, icord2, icord3, ns_in0, irivp, iocord_tmp 
      integer, allocatable :: nelm2(:)
      integer izunit,nin,ii,n_elem,ns_elem,ie, e_mem(8)
      integer neq_sv, nei_in_sv, icall_sv, neq_p, iblanking_value
      integer i_pri,i_sec
      logical zone_saved
c gaz 060822
      logical opnd
      real*8 hdum, sdum, px, py, pz, flxdum
c gaz 050321,060822
      integer mout_dim, num_unit
      real*8 vx_dum,vy_dum,vz_dum
      real*8 pdum, tdum, rolconv, dumconv, dumconv1
      real*8 sxx_dum, syy_dum, szz_dum, sxy_dum, sxz_dum, syz_dum
      real*8 pi,sat_out_tol
c      character*80 title(2*maxscalar+3)
      character*150 :: tecstring = ''
      character*150 :: tecstring_riv = ''
      character*500 string
c      character*20 vstring
      character*43 tstring
      character*5 char_type
      character*3 dls
      character*30 zonesavename, char_temp
      character*6 zonestring
c gaz 050221
      integer maxtitle, mout, neq_write
      parameter(maxtitle = 22)
      logical :: xon = .false., yon = .false., zon = .false.
c gaz 061322
      logical exists
      character*5 dual_char
      character*42 title(maxtitle), units(maxtitle), pstring
      character*42, allocatable :: title_kd(:)
      character*600 print_title, vstring
      character*300 temp_string
      character*15 nform
      parameter(iblanking_value = -9999)
      parameter (pi=3.1415926535)
      parameter (sat_out_tol=1.d-98)
      save tecstring, tecstring_riv
C     BEGIN
      size_head = size(head)
      size_pcp = size(pcp)
C     nscalar=(iovapor+ioliquid)*iopressure+iosaturation+iotemperature
C     calculation done in avs_io()
      ioLP = ioliquid*iopressure
      ioVP = iovapor*iopressure

c--------------------------------------------
c  Shaoping add  10-19-2017
      zone_saved =.false.
c-------------------------------------------
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
c gaz 050921 added below after avs output crashed      
      temp_string = ''
      print_title = ''
      pstring = ''
      dual_char = '' 
      iocord_tmp = iocord
c gaz 050621      
      if (iodual .eq. 1.and.ifdual.ne.0) dual_char = 'Dual '
      if (iodual .ne. 1.and.ifdual.ne.0) dual_char = 'GDKM '    
      if (iogdkm .ne. 0 .and. ifdual .ne. 0) then
c     Output for gdkm nodes
        if(iogdkmblank.eq.0) then
         istart = neq_primary + 1
         iend = neq
         offset = 0
         nadd = 0
        else
c gaz 040917 gdkm blanking            
         istart = 1
         iend = neq_primary
         offset = 0
         nadd = 0            
        endif
        irivp = 0 
c         if (icnl .eq. 0) then
c            iocord = 3
c         else
c            iocord = 2
c         end if
c gaz 051421 added ckeck for coordinate output        
        if(iocord_tmp.ne.0) then 
         if (icnl .eq. 0) then
            iocord = 3
         else
            iocord = 2
         end if
        endif     
      else if (ifdual .ne. 0)then
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
c gaz 050221 material file heading
      title(1) = trim(dual_char) // 'Permeability (m**2) in X'
      title(2) = trim(dual_char) // 'Permeability (m**2) in Y'
      title(3) = trim(dual_char) // 'Permeability (m**2) in Z'
      title(4) = trim(dual_char) // 'Thermal Conductivity (W/m*K) in X'
      title(5) = trim(dual_char) // 'Thermal Conductivity (W/m*K) in Y'
      title(6) = trim(dual_char) // 'Thermal Conductivity (W/m*K) in Z'
      title(7) = trim(dual_char) // 'Porosity'
      title(15)= trim(dual_char) // 'Rock bulk density (kg/m**3)'
      title(8) = trim(dual_char) // 'Rock specific heat (MJ/kg*K)'
      title(9) = trim(dual_char) // 'Capillary pressure (MPa)'
      title(10)= trim(dual_char) // 'Relative permeability model'
      title(11)= trim(dual_char) // 'Capillary pressure model'
      title(12) = 'X coordinate (m)'
      title(13) = 'Y coordinate (m)'
      title(14) = 'Z coordinate (m)'         
      
      units(1) = '(m**2)'
      units(2) = '(m**2)'
      units(3) = '(m**2)'
      units(4) = '(W/m*K)'
      units(5) = '(W/m*K)'
      units(6) = '(W/m*K)'
      units(7) = '(non dim)'
      units(15) = '(kg/m**3)'
      units(8) = '(MJ/kg*K)'
      units(9) = '(MPa)'
      units(10)= '(flag)'
      units(11)= '(flag)'
      ic1 = 1
      
      select case (icnl)
      case (1, 4)
         xon = .true.
         yon = .true.
      case (2, 5)
         xon = .true.
         zon = .true.
      case(3, 6)
         yon = .true.
         zon = .true.
      case default
         xon = .true.
         yon = .true.
         zon = .true.
      end select
      if(altc(1:3).eq.'avs' .and. altc(4:4) .ne. 'x') then
        if(icnl.eq.0) then
          mout_dim = 3
        else
          mout_dim = 2
        endif
         write (temp_string, '(i2)') mout + mout_dim
         ic2 = len_trim(temp_string)
         pstring(ic1:ic2) = temp_string
         write (temp_string, '(a2)') ' 1'
         do i = 1, mout+ mout_dim
            ic1 = ic2 + 1
            ic2 = ic2 + len_trim(temp_string)
            pstring(ic1:ic2) = temp_string          
         end do
         length = len_trim(pstring)
         write (lu, '(42a)') pstring(1:length)
c gaz  051521 fixed unit problem    
! Coordinates will be written  
         if(iocord_tmp. ne. 0) then
            if (xon) write (lu, 101) trim(title(12))
            if (yon) write (lu, 101) trim(title(13)) 
            if (zon) write (lu, 101) trim(title(14))
         endif   
         if (idoff .ne. -1) then
! Permeability will be written
            if (xon) write (lu, 101) trim(title(1))
            if (yon) write (lu, 101) trim(title(2))
            if (zon) write (lu, 101) trim(title(3))
         end if     
         if (ico2 .gt. 0 .or. ice .ne. 0) then
! Conductivity will be written
            if (xon) write (lu, 101) trim(title(4))
            if (yon) write (lu, 101) trim(title(5))
            if (zon) write (lu, 101) trim(title(6))
         end if
! Porosity, bulk density and specific heat will be written
         write (lu, 101) trim(title(7)), trim(units(7))
         write (lu, 101) trim(title(15))
         write (lu, 101) trim(title(8))
         if (irdof .ne. 13) then
! Capillary pressure will be written
            write(lu, 101) trim(title(9))
         end if
         if (rlp_flag .eq. 1) then
! rlp and cap model flags will be written if rlp_flag .eq. 1
            write (lu, 101) trim(title(10)), trim(units(10))
            write (lu, 101) trim(title(11)), trim(units(11))
         end if
         if (iccen .eq. 1 .and. iokd .eq. 1 .and. rxn_flag .eq. 0) then
! Kd will be output for each transport specie and model
            allocate (title_kd(nspeci))
            do i = 1, nspeci
               title_kd(i) = trim(cpntnam(i)) // ' (Kd l/kg) (l/kg)'
            end do
         end if
      else
         ic1 = 1
         ic2 = 0
         if (altc(1:4) .eq. 'avsx') then
            dls = ' : '
            write (nform, 200) dls
            temp_string = 'node'
         else if (altc(1:3) .eq. 'tec') then
!            nform = '(' // "'" // ' "' // "', a, '" // '"' // "')"
            write (nform, 300)
            if (iocord .ne. 0) then
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
                  icord1 = 2
                  icord2 = 3
                  icord3 = 1
               case default
                  icord1 = 1
                  icord2 = 3
                  icord3 = 1
               end select
! Write X coordinate
               if (icnl .ne. 3 .and. icnl .ne. 6) then
                  write(temp_string,nform) trim(title(12))
                  ic2 = ic2 + len_trim(temp_string)
                  print_title(ic1:ic2) = temp_string
                  ic1 = ic2 + 1
               end if
! Write Y coordinate
               if (icnl .ne. 2 .and. icnl .ne. 5) then
                  write(temp_string,nform) trim(title(13))
                  ic2 = ic2 + len_trim(temp_string)
                  print_title(ic1:ic2) = temp_string
                  ic1 = ic2 + 1
               end if
! Write Z coordinate
               if (icnl .ne. 1 .and. icnl .ne. 4) then
                  write(temp_string,nform) trim(title(14))
                  ic2 = ic2 + len_trim(temp_string)
                  print_title(ic1:ic2) = temp_string
                  ic1 = ic2 + 1
               end if
               write (temp_string, fmt=nform) 'node'
            else
               write (temp_string, fmt=nform) 'node'
            end if
         else if (altc(1:3) .eq. 'sur') then
            dls = ', '
            write (nform, 200) dls
            temp_string = 'node'
         end if
         ic2 = ic2 + len_trim(temp_string)
         print_title(ic1:ic2) = temp_string
         ic1 = ic2 + 1
         if (idoff .ne. -1) then
! Permeability will be written
            if (xon) then
               write (temp_string, fmt=nform) trim(title(1))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
            if (yon) then
               write (temp_string, fmt=nform) trim(title(2))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
            if (zon) then
               write (temp_string, fmt=nform) trim(title(3))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
         end if     
         if (ico2 .gt. 0 .or. ice .ne. 0) then
! Conductivity will be written
            if (xon) then
               write (temp_string, fmt=nform) trim(title(4))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
            if (yon) then
               write (temp_string, fmt=nform) trim(title(5))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
            if (zon) then
               write (temp_string, fmt=nform) trim(title(6))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
         end if
! Porosity, bulk density, and specific heat will be written
         write (temp_string, fmt=nform) trim(title(7))
         ic2 = ic2 + len_trim(temp_string)
         print_title(ic1:ic2) = temp_string
         ic1 = ic2 + 1
         write (temp_string, fmt=nform) trim(title(15))
         ic2 = ic2 + len_trim(temp_string)
         print_title(ic1:ic2) = temp_string
         ic1 = ic2 + 1
         write (temp_string, fmt=nform) trim(title(8))
         ic2 = ic2 + len_trim(temp_string)
         print_title(ic1:ic2) = temp_string
         ic1 = ic2 + 1
         if (irdof .ne. 13) then
! Capillary pressure will be written
            write(temp_string, fmt=nform) trim(title(9))
            ic2 = ic2 + len_trim(temp_string)
            print_title(ic1:ic2) = temp_string
            ic1 = ic2 + 1
         end if
         if (rlp_flag .eq. 1) then
! rlp and cap model flags will be written if rlp_flag .eq. 1
            write (temp_string, fmt=nform) trim(title(10))
            ic2 = ic2 + len_trim(temp_string)
            print_title(ic1:ic2) = temp_string
            ic1 = ic2 + 1
            write (temp_string, fmt=nform) trim(title(11))
            ic2 = ic2 + len_trim(temp_string)
            print_title(ic1:ic2) = temp_string
            ic1 = ic2 + 1
         end if
! Kd will be written
         if (iccen .eq. 1 .and. iokd .eq. 1 .and. rxn_flag .eq.0) then
            allocate (title_kd(nspeci))
            do i = 1, nspeci
               title_kd(i) = trim(cpntnam(i)) // ' (Kd l/kg)'
               write (temp_string, fmt=nform) trim(title_kd(i))
               ic2 = ic2 + len_trim(temp_string)
               print_title(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end do
         end if
         length = len_trim(print_title)
         if (altc(1:3) .ne. 'tec') then
            write (lu, '(a)') print_title(1:length)
         else
            write (lu, 124) verno, jdate, jtime, trim(wdd)
            if (iogrid .eq. 1 .and. iocord .eq. 0) then
               write (lu, 139)
            end if
            write (lu, '("VARIABLES = ", a)') print_title(1:length)
            if (iogeo .eq. 1) then
               if(iriver.ne.2.and.gdkm_flag.ne.1) then
                 neq_write = neq
               else
                 neq_write = neq_primary
               endif                   
               select case (ns_in)
               case (5,6,8)
                  write (temp_string, 134) neq_write, nei_in, 'FEBRICK'
               case (4)
                  if (icnl .eq. 0) then
                     write (temp_string, 134) neq_write, nei_in, 
     &                    'FETETRAHEDRON'
                  else
                     write (temp_string, 134) neq_write, nei_in, 
     &                    'FEQUADRILATERAL'
                  end if
               case (3)
                  write (temp_string, 134) neq_write, nei_in, 
     &                    'FETRIANGLE'
               case (2)
                  write (temp_string, 134) neq_write, nei_in, 
     &                    'FELINESEG'
               case (0)
! fdm grid
                  write (temp_string, '(a)') ''
               end select
               write (lu, 129) trim(temp_string)
            end if
         endif
      end if
 124  format('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
 129  format('ZONE T = "Material properties"', a)
 134  format(', N = ', i8, ', E = ', i8, ', DATAPACKING = POINT',
     &     ', ZONETYPE = ', a)
 139  format('FILETYPE = "SOLUTION"')
 200  format("( ' ", a, "', a)")
 300  format("('",' "', "', a, '",'"',"')")

      temp_string = ''
      vstring = ''

c gaz 050221 end material file heading         
      nout = nscalar + iocord
 
      if (icall .eq. 1.and.irivp.eq.0) tecstring = ''
      if (icall .eq. 1.and.irivp.ne.0) tecstring_riv = ''

c     Surfer headers need to be written to each zone file
      if (altc(1:3) .ne. 'sur') then
c gaz  043021 see what is there for a header        
c        call write_avs_head_s(icall,nscalar,lu,ifdual,0,iriver2)
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
c open and read saved zone file if they exist
            call zone_saved_manage(1,izunit,idz,nin,n_elem,zone_saved)
c            
            if(zone_saved) then
             zonestring = ' '   
             write(zonestring(1:5),'(i5)') idz
             neq_sv = neq
             nei_in_sv = nei_in
             icall_sv = icall
             neq = nin
             nei_in = n_elem
             icall = 1
             iogeo = 1
c copy FETYPE stuff in here   
            if (altc(1:3) .eq. 'tec') then
               string = ''
               if (icall .eq. 1 .and. iogeo .eq. 1) then
                  select case (ns_in)
                  case (5,6,8)
                     write (string, 135) neq, nei_in, 'FEBRICK'
                  case (4)
                     if (icnl .eq. 0) then
                        write (string, 135) neq, nei_in, 
     &                       'FETETRAHEDRON'
                     else
                        write (string, 135) neq, nei_in, 
     &                       'FEQUADRILATERAL'
                     end if
                  case (3)
                     write (string, 135) neq, nei_in, 'FETRIANGLE'
                  case (2)
                     if(irivp.eq.0) then
                        write (string, 135) neq, nei_in, 'FELINESEG'
                        ns_in=ns_in0
                     else if(irivp.eq.2)then                    
                        write (string, 135) npoint_riv, npoint_riv-1,
     &                       'FELINESEG'
                        ns_in=ns_in0
                     endif
                  case (0)
c     fdm grid
                     write (string, '(a)') ''
                  end select   
                  if(irivp.eq.0) then
c gaz 111516 mods to write zone number                      
                     if (iogeo .eq. 1) then
                        tecstring = trim(string)
                        write (lu, 131) trim(zonestring),  
     &                       trim(timec_string), trim(tecstring)
                     else if (iogrid .eq. 1) then
                        tecstring = trim(string)
                        write (lu, 130) trim(timec_string), 
     &                       trim(gridstring), trim(times_string)
                     else
                        write (lu, 130) trim(timec_string), 
     &                       trim(tecstring)
                     end if
                  else
                     tecstring_riv = trim(string)
                     write (lu, 130) trim(timec_string),
     &                    trim(tecstring_riv), trim(gridstring), 
     &                    trim(times_string)
                  endif   
               endif
c end copy FETYPE             
            endif
c         endif
c  an else here????         
         else if (altc(1:3) .eq. 'tec') then
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
c gaz 040517 gdkm needs small modification
               if(gdkm_flag.ne.0) then
                neq_p = neq_primary
               else
                neq_p = neq
               endif 
               string = ''
               if (icall .eq. 1 .and. iogeo .eq. 1) then

               else if (icall .eq. 1 .and. iocord .ne. 0) then
c gaz 050221 matertial file does not need time(_string)                  
c         
               else if (icall .eq. 1 .and. iozid .ne. 0) then
c   gaz 050221    write (lu, 130) trim(timec_string)
                  write (string, 125) '2', iz
                  if(irivp.eq.0) then
                     tecstring = trim(string)
                  else
                     tecstring_riv = trim(string)
                  endif
               else
                  if(irivp.eq.0) then
                     if (iogeo .eq. 1) then
c gaz 050221       write (lu, 130) trim(timec_string), 
c     &                       trim(tecstring)
                     else if (iogrid .eq. 1) then
c                        write (lu, 130) trim(timec_string), 
c     &                       trim(tecstring),
c     &                       trim(gridstring), trim(times_string)
                     else
c                        write (lu, 130) trim(timec_string), 
c   &                       trim(tecstring)                        
                     end if   
                  else
c                     write (lu, 130) trim(timec_string),
c     &                    trim(tecstring_riv), trim(gridstring), 
c     &                    trim(times_string)
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
      
         if (iohead .eq. 1 .and. ichead .eq. 1) then
            ihead=1
            dumconv = crl(1,1)
            dumconv1 = crl(4,1)
            pdum = pres0+rol0*head0*(-grav)
            tdum = temp0        
            call water_density(tdum,pdum,rolconv)
            crl(1,1)=rolconv
            crl(4,1)=pres0
            rho1grav = rolconv*9.81d-6
         end if
c if saved zone exists use 1, nin form
         if(zone_saved) then
          istart = 1
          iend = nin
         endif
         do ii = istart, iend
c     Node loop          
            string = ''
            if (iozone .ne. 0) then
               if(zone_saved) then
                i = ncord(ii)
               else
                i = ii
                if (izone_surf_nodes(i).ne.idz) goto 201
               endif
            else
              i = ii
            end if
c gaz 052221 need this for coordinate output            
          i_pri = i  
c     Node number will be written first for avs and sur files
c gaz 040517 this is where gdkm number is changed            
            if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'sur') then  
c gaz 052221 need i_pri = i for gdkm coordinates
             if(ifdual.ne.0.and. iogdkm .eq. 1.
     &                        and. iogdkmblank .eq. 1) then
c gaz  040817 if blanking use primary grid node number (and blank variable) 
c gaz identify secondary node variable (i_sec)                   
                   write(string, 100) i
                   i_pri = i
                   if(nelm(nelm(i_pri+1)).gt.neq_primary) then
                    i_sec = nelm(nelm(i_pri+1))
                    i = i_sec
                   else  
c gaz use a blanking value
                    i_sec = iblanking_value
                    i = i_sec
                   endif                   
               else if (ifdual .eq. 0 .or. iogdkm .eq. 1) then
                    write(string, 100) i
               else 
                  write(string, 100) i - neq_primary    
               end if
c ic1 positions the column for printout               
              ic1 = 11
c
            else
               ic1 = 1
            end if
            if (iocord .ne. 0) then
c Only output coordinates that are used
               if (altc(1:3) .eq. 'tec' .and. icall .ne. 1) then
c     Coordinates will only be output in the first file for tecplot
c     (do nothing)
c gaz 051121 added 'tec' condition            
               elseif (altc(1:3) .eq. 'tec' .and. 
     &            icall .eq. 1) then
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
               if(ifdual.ne.0.and. iogdkm .eq. 1.
     &                        and. iogdkmblank .eq. 1) then
c gaz  040817 if blanking use primary grid node number (and blank variable) 
c gaz identify secondary node variable (i_sec)                   
                   write(vstring, 105) dls(1:k), i
                   i_pri = i
                   if(nelm(nelm(i_pri+1)).gt.neq_primary) then
                    i_sec = nelm(nelm(i_pri+1))
                    i = i_sec
                   else  
c gaz use a blanking value
                    i_sec = iblanking_value
                    i = i_sec
                   endif                   
               else if (ifdual .eq. 0 .or. iogdkm .eq. 1) then
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
     &              .or. icall .eq. 1 .or. iogrid .eq. 1) then
                  write(vstring, 115) dls(1:k), izonef(i)
                  ic2 = ic1 + len_trim(vstring)
                  string(ic1:ic2) = vstring
                  ic1 = ic2 + 1
               end if
            end if


c gaz 050321 perm            
            if (iomaterial .eq. 1.and.idoff .ne. -1) then
c add coordinates here  
c gaz 052221 always use i_pri for coordinates               
            if(iocord_tmp.ne.0.and.altc(1:3).ne.'tec') then
             if(xon) then
               write(vstring,110) dls(1:k), corz(i_pri,1)
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
             endif
             if(yon) then
               write(vstring,110) dls(1:k), corz(i_pri,2)
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
              endif 
             if(zon) then
               write(vstring,110) dls(1:k), corz(i_pri,3)
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
             endif               
            endif
c permeability                  
             if(xon) then    
              if(i.ne.iblanking_value)then
               write(vstring,110) dls(1:k), pnx(i)*1.e-6
              else
               write(vstring,110) dls(1:k), iblanking_value
              endif
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
             endif
             if(yon) then  
              if(i.ne.iblanking_value)then 
               write(vstring,110) dls(1:k), pny(i)*1.e-6
              else
               write(vstring,110) dls(1:k), iblanking_value
              endif 
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
             endif
             if(zon) then  
              if(i.ne.iblanking_value)then 
               write(vstring,110) dls(1:k), pny(i)*1.e-6
              else
               write(vstring,110) dls(1:k), iblanking_value
              endif 
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
             endif
            endif 
c thermal conductivity  
c gaz 051421            
            if (ico2 .gt. 0 .or. ice .ne. 0) then
               if(i.ne.iblanking_value)then
                 vx_dum = thx(i)*1.e06
                 vy_dum = thy(i)*1.e06
                 vz_dum = thz(i)*1.e06
                  write(vstring,110) dls(1:k), vx_dum
               else
               write(vstring,110) dls(1:k), iblanking_value
               endif
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
              if(i.ne.iblanking_value)then 
               write(vstring,110) dls(1:k), vy_dum
              else
               write(vstring,110) dls(1:k), iblanking_value
              endif 
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
              if(i.ne.iblanking_value)then 
               write(vstring,110) dls(1:k), vz_dum
              else
               write(vstring,110) dls(1:k), iblanking_value
              endif 
               ic2 = ic1 + len_trim(vstring)
               string(ic1:ic2) = vstring
               ic1 = ic2 + 1
            endif
       
            if (iomaterial .eq. 1) then      
              if(i.ne.iblanking_value)then
! Porosity, bulk density, and specific heat will be written
               write (temp_string, '(3(a,1p,g14.6))') dls, ps(i), 
     &           dls, denr(i), dls, cpr(i)
              else
               write (temp_string, '(3(a,i14))') dls, iblanking_value,
     &           dls, iblanking_value, dls, iblanking_value
              endif
            ic2 = ic2 + len_trim(temp_string)
            string(ic1:ic2) = temp_string
            ic1 = ic2 + 1
            endif
            if (iomaterial .eq. 1) then  
! Capillary pressure will be written             
             if (irdof .ne. 13) then               
              if(i.ne.iblanking_value)then
               write (temp_string, '(a,1p,g14.6)') dls, pcp(i)
              else
               write (temp_string, '(a,1p,i14)') dls, iblanking_value
              endif 
               ic2 = ic2 + len_trim(temp_string)
               string(ic1:ic2) = temp_string
               ic1 = ic2 + 1
             end if
            endif
            if (iomaterial .eq.1.and.rlp_flag .eq.1) then  
! rlp and cap model flags will be written if rlp_flag .eq. 1
              if(i.ne.iblanking_value) then  
               write (temp_string, '(2(a,i14))') dls, irlp(i), 
     &              dls, icap(i) 
              else
               write (temp_string, '(2(a,i14))') dls, iblanking_value, 
     &              dls, iblanking_value   
              endif
               ic2 = ic2 + len_trim(temp_string)
               string(ic1:ic2) = temp_string
               ic1 = ic2 + 1
            end if
! Kd will be written
            if(iomaterial .eq.1.and.iccen .eq. 1 .and. iokd .eq. 1 .and.
     &           rxn_flag .eq.0) then
               do i = 1, nspeci                  
                if(i.ne.iblanking_value) then  
                 write (temp_string, '(a,g14.6)') dls, a1adfl(i,itrc(j))
                else
                 write (temp_string, '(a,i14)') dls, iblanking_value
               endif
               ic2 = ic2 + len_trim(temp_string)
               string(ic1:ic2) = temp_string
               ic1 = ic2 + 1                                    
              end do
            end if                             
            length = len_trim(string)
            write(lu,'(a)') string(1:length)                
 201        enddo
c            
c add element information here for saved zones and then exit
c
         if(zone_saved) then
           do i = 1, n_elem
            write(lu,'(9(1x,i7))') 
     &        (ncord_inv(elem_temp(i,j)),j = 1,ns_in)
           enddo
           neq = neq_sv
           nei_in = nei_in_sv
           icall = icall_sv
           deallocate(elem_temp)
         endif
         if (iohead .eq. 1 .and. ichead .eq. 1) then
            crl(1,1)= dumconv
            crl(4,1)= dumconv1
            ihead=0
            if(ico2.lt.0) then
               rho1grav = crl(1,1)*(9.81d-6)
            else
               rho1grav = rol0*9.81d-6
            endif
         end if

         call flush(lu)
         if (altc(1:3) .eq. 'sur') close (lu)
      enddo
      if(zone_saved) return
      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. iogeo .eq. 1)
     &     then
c     Read the element connectivity and write to tec file
         if (ifdual .eq. 1 .and. iogdkm .eq. 1) then
c     Do nothing unless blanking used
          if(iogdkmblank.ne.0) then
c gaz 040817 attach geometry to gdkm file
            il = open_file(geoname,'old')
c     avsx geometry file has an initial line that starts with neq_primary
            read(il,*) i
            if (i .ne. neq_primary) backspace il
            do i = 1, neq_primary
               read(il,*)
            end do
            allocate (nelm2(ns_in))
            do i = 1, nei_in
               read (il,*) i1,i2,char_type,(nelm2(j), j=1,ns_in)
               write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
            end do
            deallocate(nelm2)
            close (il)              
          endif
          else if(irivp.eq.0) then
            if(nact_elem_geo.eq.0) then
             inquire(file = geoname, number = num_unit)
             close(num_unit)
             il = open_file(geoname,'old')
c     tec geometry file has 4 initial lines t
             read(il,*) 
             read(il,*)
             read(il,*)
             read(il,'(a20,i9)') wdd(1:20),i
             if (i .ne. neq_primary) backspace il
             do i = 1, neq
               read(il,*)
             end do
             allocate (nelm2(ns_in))
             do i = 1, nei_in
               read (il,*) (nelm2(j), j=1,ns_in)
               write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
             end do
             deallocate(nelm2)
             close (il)
            else
             do i = 1, nei_in
               write(lu,'(8(i8))') (elem_geo((i-1)*ns_in+j),j=1,ns_in)
            end do
            endif
         else
c     river segments (2 node elements)
            do i = 1,nnelm_riv
               write(lu,'(2(i8))') nelm_riv(i,1)-neq,
     &              nelm_riv(i,2)-neq
            enddo
         endif
      end if
c gaz added element output  (hex only) for fdm generated grid   
      if (icall .eq. 1 .and. altc(1:3) .eq. 'tec' .and. ivf .eq. -1
     &     .and. ifdm_elem. eq. 1) then
c first generate elements   
c gaz 061422 simplified  element        
          allocate (nelm2(ns_in))
          do i = 1, nei_in
            write(lu, '(8(i8))') (elem_geo((i-1)*ns_in+j), j=1,ns_in)
          end do
          deallocate(nelm2)          

      end if        
      if (altc(1:3) .ne. 'sur') close (lu)
      iocord = iocord_tmp

 100  format(i10.10)
c gaz 051021 (format 101 is the old 100 in  write_avs_node_mat.f  
c used in avs output   
C gaz 051521      
 101  format(a, ' ', a)
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
 130  format('ZONE T =', a, a, a, a) 
 131  format('ZONE ',a,',',' T =', a, a, a, a) 
 135  format(', N = ', i8, ', E = ', i8, ', DATAPACKING = POINT',
     &     ', ZONETYPE = ', a)
 140  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')

      return
      end

