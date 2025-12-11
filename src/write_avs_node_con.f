      subroutine write_avs_node_con(icall,npt,neq,nspeci,
     .     lu, ifdual)
c gaz debug 110914
c      subroutine write_avs_node_con(icall,an,anv,npt,neq,nspeci,
c     .     lu, ifdual)
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
CD1 Output concentration fields from FEHM
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
CD2 $Log:   /pvcs.config/fehm90/src/write_avs_node_con.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:24:40   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:14   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:42   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:30   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.9   Fri Feb 02 14:19:38 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.8   09/11/95 17:30:46   awolf
CD2 Set so all concentration values less than 1e-20 are set to 1e-20
CD2 for use with avs
CD2 
CD2    Rev 1.7   08/07/95 11:49:36   awolf
CD2 Fixed for avs output of dpdp
CD2 
CD2    Rev 1.6   04/10/95 15:16:48   llt
CD2 changed max function for IBM compatibility
CD2 
CD2    Rev 1.5   04/10/95 14:44:04   robinson
CD2 Minimum output concentration is now 0 rather than 1.e-20
CD2 
CD2    Rev 1.4   01/28/95 16:14:06   llt
CD2 needed endif
CD2 
CD2    Rev 1.3   01/28/95 14:21:10   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.2   01/20/95 13:29:12   tam
CD2 Changed format for strings from * to a56, kept length to 80 so left justified
CD2 
CD2    Rev 1.1   11/29/94 13:42:32   llt
CD2 Removed protying, so could run under cc compile, instead of c compiler. 
CD2 Used #ifdef to determine format for routine name, depending on machine.
CD2 (Changes made by tam.)
CD2 
CD2    Rev 1.0   08/23/94 15:33:54   llt
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
C***********************************************************************
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
C***********************************************************************
CPS
CPS PSEUDOCODE
CPS
CPS BEGIN 
CPS   
CPS END 
CPS 
C***********************************************************************

      use avsio
      use comai, only : altc, days, grav, iadif, icnl, ico2, idof, 
     &     ichead, ihead, nei_in, ns_in, phi_inc, istrs, ivf,
     &     neq_primary, rho1grav, ifdm_elem, igrav, ns, gdkm_flag,
     &     verno,jdate,jtime,wdd, icconc 
      use combi, only : corz, izonef, nelm, ncord, ncord_inv, elem_temp,
     &     contour_conc_files
      use comchem
      use comdi, only : nsurf, izone_surf, izone_surf_nodes, icns,
     &     an, anv

      use compart, only : ptrak, pout
      use comriv, only : npoint_riv, nnelm_riv, nelm_riv, iriver
      use comrxni
      use comdti
      implicit none

      integer add_dual, maxcon, iz, idz, iendz, il, open_file
      integer neq,nspeci,lu,ifdual,icall,length,i1,i2
      integer icord1, icord2, icord3, iaq, ivap, isolid
      integer npt(*), ip1, ip2, ic1, ic2
      integer, allocatable ::  nelm2(:)
      parameter (maxcon = 100)
      real*8, allocatable :: an_dum(:,:)
      real*8, allocatable :: anv_dum(:,:)
      real*8, allocatable :: antmp(:,:)
c gaz 062717      
c      character*60, allocatable :: title(:)
      character*14 tailstring
      character*8 dual_char
c gaz 062717      
c      character*3 dls, snum
       character*3 snum
c      character*5 char_type
      character*60 fstring
      character*35 tmpname
      character*30 cordname(3)
c      character*150 :: string = '', tecstring = '', sharestring = ''
      real*8 write_array(maxcon)
      integer i, ic, im, in, iv, ix, istep, j, k, n, iblanking_value
      integer t1(maxcon),itotal2,write_total, iocord_temp,ns_in0

      logical zone_saved
      integer izunit,nin,ii,n_elem,ns_elem,ie, e_mem(8)
      integer neq_sv, nei_in_sv, icall_sv, neq_p
      integer i_pri,i_sec,nscalar

c      character*80 title(2*maxscalar+3)
      character*150 :: tecstring = ''
      character*150 :: tecstring_riv = ''
      character*500 string
      character*20 vstring
      character*43 tstring
      character*5 char_type
      character*3 dls
      character*30 zonesavename, char_temp
      character*6 zonestring
      character*500 sharestring
      character*6 share_string
      character*150 :: tec_string = ''
c      character*45 title(3*maxvector), title2(2)
      integer maxvector
      parameter (maxvector = 3)
      character*60, allocatable :: title(:)
      parameter(iblanking_value = -9999)
      real*8 pi
      parameter (pi=3.1415926535)
      integer nadd, istart, iend, irivp, iocord_tmp , offset, iriver2
      integer irxn_title
      integer i_wrt

      real*8 complex_conc
      real*8 minc, maxc
c gaz 072917      
      character*200 file_flux

      parameter (minc = 1.0d-90, maxc = 1.0d+20)

      save tecstring, sharestring, tecstring_riv

c--------------------------------------------
c  Shaoping add  10-19-2017
      zone_saved =.false.
c-------------------------------------------

      data cordname(1) / 'X (m)' /, 
     &     cordname(2) / 'Y (m)' /,
     &     cordname(3) / 'Z (m)' /  

      allocate(an_dum(n0,nspeci))
      allocate(anv_dum(n0,nspeci))
      allocate(antmp(n0,nspeci))
      do i = 1, nspeci
       ip1 = (i-1)*n0 + 1
       ip2 = ip1 + n0 -1
       an_dum(1:n0,i) = an(ip1:ip2)
       anv_dum(1:n0,i) = anv(ip1:ip2)
      enddo
      iocord_temp = iocord
   
      irivp = 0
      iriver2 = iriver
      
      iocord_tmp = iocord
      
      if (iogdkm .ne. 0 .and. ifdual .ne. 0) then
c     Output for gdkm nodes
        if(iogdkmblank.eq.0) then
         istart = neq_primary + 1
c gaz 070118         
c         iend = n0
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
         if (icnl .eq. 0) then
            iocord = 3
         else
            iocord = 2
         end if
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
 
      if (icall .eq. 1.and.irivp.eq.0) tecstring = ''
      if (icall .eq. 1.and.irivp.ne.0) tecstring_riv = ''
 
   
      if(ifdual .eq. 0)then
         istep = 0
         add_dual=0
         dual_char = ''
         tailstring = '_con_node'
      else
         istep = maxcon
c gaz       add_dual=neq   
         if (iodual .eq. 1) then
            add_dual=neq
            dual_char = 'Dual '
            tailstring = '_con_dual_node'
         else if (iogdkm .eq. 1) then
            add_dual=0
            dual_char = 'GDKM '
            tailstring = '_con_gdkm_node'
            if (icnl .eq. 0) then
               iocord = 3
            else
               iocord = 2
            end if
         end if
      endif

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
            icord1 = 1
            icord2 = 3
            icord3 = 1
         case default
            icord1 = 1
            icord2 = 3
            icord3 = 1
         end select
      end if

      if (iozone .ne. 0 ) then
         iendz = nsurf
      else 
         iendz = 1
         idz = iozone
      end if
c different than scalars 
       if (altc(1:4) .eq. 'avsx') then
         dls = ' : '
      else if (altc(1:3) .eq. 'sur') then
         dls = ', '
      else
         dls = ' '
      end if     
c       
c       
      if (altc(1:3) .ne. 'sur') then
         call namefile2(icall,lu,ioformat,tailstring,0)
! file will be opened in zone loop for surfer
      end if
         icall_sv = icall
c gaz 062717
         do iz = 1, iendz
c     Zone loop 
c if block 1             
            if (iozone .ne. 0) then
               idz = izone_surf(iz) 
c open and read saved zone file if they exist
            call zone_saved_manage(1,izunit,idz,nin,n_elem,zone_saved)
c   gaz 062717         
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
c             irivp = 0
c gaz 112716 FE geometry goes here   
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
c               endif       
        

               endif
c end copy FETYPE             
            endif
c         endif
c gaz end  FE geometry  
c  if block 1
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
c end id block 1             
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
                  select case (ns_in)
                  case (5,6,8)
                     write (string, 135) neq_p, nei_in, 'FEBRICK'
                  case (4)
                     if (icnl .eq. 0) then
                        write (string, 135) neq_p, nei_in, 
     &                       'FETETRAHEDRON'
                     else
                        write (string, 135) neq_p, nei_in, 
     &                       'FEQUADRILATERAL'
                     end if
                  case (3)
                     write (string, 135) neq_p, nei_in, 'FETRIANGLE'
                  case (2)
                     if(irivp.eq.0) then
                        write (string, 135) neq_p, nei_in, 'FELINESEG'
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
c (gaz put back in 072419)not needed below 
c use this info to write header later.
                if(icall.eq.1) then
                 if(irivp.eq.0) then
                     if (iogeo .eq. 1) then
                        tecstring = trim(string)
                        tecstring = trim(timec_string)//trim(tecstring)
                     else if (iogrid .eq. 1) then

                        tecstring = trim(timec_string)// 
     &                       trim(gridstring)//trim(times_string)
                     else
                        tecstring = trim(timec_string)// 
     &                       trim(tecstring)
                     end if
                  else
                     tecstring_riv = trim(string)
c                     write (lu, 130) trim(timec_string),
c     &                    trim(tecstring_riv), trim(gridstring), 
c     &                    trim(times_string)
                  endif 
                endif  
                if(icall .eq. 1 .and. iocord .ne. 0) then
                  sharestring = ''
                  if (icnl .eq. 0) then
                     if (iozid .eq. 0) then
                        write (sharestring, 140) '1-3', iz
                     else
                        write (sharestring, 140) '1-3,5', iz
                     end if
                  else
                     if (iozid .eq. 0) then
                        write (sharestring, 140) '1-2', iz
                     else
                        write (sharestring, 140) '1-2, 4', iz
                     end if
                  end if  
                endif

c not needed above  (not need for icall = 1), save for later                
                  
               else if (icall .eq. 1 .and. iocord .ne. 0) then
                  

               else if (icall .eq. 1 .and. iozid .ne. 0) then
                  write (lu, 130) trim(timec_string)
                  write (string, 125) '2', iz
                  if(irivp.eq.0) then
                     tecstring = trim(string)
                  else
                     tecstring_riv = trim(string)
                  endif
               else
                  if(irivp.eq.0) then
                     if (iogeo .eq. 1) then
c                        write (lu, 130) trim(timec_string), 
c     &                       trim(tecstring)
                     else if (iogrid .eq. 1) then
c                        write (lu, 130) trim(timec_string), 
c     &                       trim(tecstring),
c     &                       trim(gridstring), trim(times_string)
                     else
c gaz 063018 should not write
c                         if(icall.ne.1) write (lu, 130)  
c     &                       trim(timec_string), trim(tecstring)      
                     end if   
                  else
c                     write (lu, 130) trim(timec_string),
c     &                    trim(tecstring_riv), trim(gridstring), 
c     &                    trim(times_string)
                  end if
                continue
              end if
            continue
           end if
          continue         
         end if


        if(rxn_flag.eq.0)then
         if(nspeci .gt. maxcon)then
            write(lu,*)'--------------------------------------------'
            write(lu,*)'ERROR: WRITE_AVS_NODE_CON'
            write(lu,*)'nspeci = ',nspeci,
     2           ' is greater than maxcon = ',maxcon
            write(lu,*)'Subroutine only able to handle up to', maxcon,
     &           'tracers'
            write(lu,*)'--------------------------------------------'
            deallocate(anv_dum,an_dum,antmp)
            return
         endif


c     Formatted write to accomodate the different way fortran 90
c     does unformatted writes (f77 way no longer worked because
c     it put commas between the numbers, avs didn't like that)

c----------------------------------------
c     PHS 4/27/00   Added new format for avsxpress!!
c---------------------------------------
         fstring = ''
         if (altc(1:3) .eq. 'tec') then
            itotal2 = nspeci+iocord+iozid+1
         else
            itotal2 = nspeci+iocord+iozid
         end if
         allocate (title(itotal2))
         do i = 1, iocord
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(i) = trim(cordname(i)) // ', (m)'
            else
               title(i) = trim(cordname(i))
            end if
         end do
         if (altc(1:3) .eq. 'tec') then
            title(i) = 'Node'
            i = i + 1
         end if
         if (iozid .eq. 1) then
            title(i) = 'Zone'
            i = i + 1
         end if
         j = i
         iaq = 0
         ivap = 0
         isolid = 0
         antmp = an_dum
         do i = 1, nspeci
            select case (icns(i))
            case (1, 2, -2)
               iaq = iaq + 1
               tmpname = cpntnam(iaq)
               if (icns(i) .eq. -2) then
                  if (tmpname(1:7) .eq. 'Aqueous') then
                     tmpname(1:5) = 'Vapor'
                     tmpname(6:17) = tmpname(8:19)
                     tmpname(18:20) = ''
                     cpntnam(iaq) = tmpname
                  end if
                  do k = 1, n0
                     antmp(k+add_dual,i) = anv_dum(k+add_dual,i)
                  end do
               end if
               title(j) = trim(dual_char)//trim(cpntnam(iaq))
            case (0)
               if (ptrak) then
                  write (snum, '(i3)') i
                  select case (abs(pout))
                  case (0)
                     tmpname = "Particles/Fluid mass Species "
                  case (1)
                     tmpname = "Particles/Total volume Species "
                  case (2)
                     tmpname = "Particles/Fluid volume Species "
                  case (3, 7)
                     tmpname = "Number of Particles Species "
                  case (4)
                     tmpname = "Cl-36 Species "
                  case (5)
                     tmpname = "Mixed Mean Concentration Species "
                  case (6)
                     tmpname = "C-14 Species "
                  end select
                  title(j) = trim(dual_char) // trim(tmpname) // snum
               else
                  isolid = isolid + 1
                  title(j) = trim(dual_char)//trim(immnam(isolid))
               end if
            case (-1)
               ivap = ivap + 1
               title(j) = trim(dual_char)//trim(vapnam(ivap))
            end select
            j = j+ 1
         end do
      
         if(altc(1:4).eq.'avsx') then
            write(fstring, 664) itotal2
            write(lu,fstring) days,(trim(title(i)),i=1,itotal2)
         else if (altc(1:3) .eq. 'avs') then
            write(lu,'(i3, 100(1x, i1))') itotal2,(1,i=1,nspeci+iocord)
            write(lu, '(a)') (trim(title(i)),i=1,itotal2)
c first write             
         else if (altc(1:3) .eq. 'tec') then
c iz is the increment of the zone loop             
            if (icall_sv .eq. 1 .or. iogrid .eq. 1) then
              if(iz.eq.1) then  
               write(lu, 98) verno, jdate, jtime, trim(wdd)
c gaz 071022 removed solution line - need to check               
c               if (iogrid .eq. 1) write(lu, 100) 
               write (fstring, 99) itotal2
               write(lu, fstring) 'VARIABLES = ', (trim(title(i)), 
     &              i=1,itotal2)
              endif
            end if
              if(iozone.ne.0) then
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
      
c if saved zone exists use 1, nin form
c start ifblock 1
         if(zone_saved) then
          istart = 1
          iend = nin
         endif
c start do loop 1   
c gaz 072519  mod icall = 1 to tecstring
            if (altc(1:3) .eq. 'tec') then
               if (iozone .eq. 0 .or. iogrid .eq. 1) then
                  if(icall.eq.1) then
c                   write (lu, 94) trim(timec_string)
                   write (lu, 94) trim(tecstring)                      
c gaz 070118 this could be the problem                 
                  else
                    write (lu, 94) trim(timec_string),trim(gridstring),
     &              trim(sharestring)
c                   write (lu, 94) trim(timec_string), trim(tecstring),
c     &                 trim(gridstring), trim(times_string)                      
                  endif
               else
                  if (icall .gt. 1 .and. iozone .ne. 0) then
                     write (tecstring, 125) trim(sharestring), iz
c gaz debug 072717
                     write (lu, 95) idz, trim(tecstring)
                  end if
               end if
               if (icall .eq. 1 .and. iz .eq. iendz) then
                  if (iogeo .eq. 1) then
                     tecstring = trim(tecstring) // trim(string)
                  else if (iocord .ne. 0 .or. iozid .ne. 0) then
                     write (tecstring, 125) trim(string), iz
                  end if   
               end if
            else if (altc(1:3) .eq. 'sur') then
               call namefile2(icall,lu,ioformat,tailstring,idz)
               fstring = ''
               write (fstring, 96) itotal2
               write(lu, fstring) (trim(title(i)), i=1,itotal2)
            end if
            fstring = ''
            if (iozid .eq. 0 .and. iocord .eq. 0) then
               write (fstring, 333) itotal2, dls
            else if (iozid .eq. 0 .and. iocord .ne. 0) then
               if (altc(1:3) .ne. 'tec') then
                  write (fstring, 333) itotal2, dls
               else
                  if (icall .eq. 1) then
                     write (fstring, 104) iocord, nspeci
                  else
                     write (fstring, 333) itotal2 - iocord, dls
                  end if
               end if
            else if (iozid .eq. 1 .and. iocord .eq. 0) then
               if (altc(1:3) .ne. 'tec' .or. icall .eq. 1 .or. 
     &              iogrid .eq. 1) then
                  write (fstring, 334) dls, nspeci, dls
               else
                  write (fstring, 333) itotal2, dls
               end if
            else if (iozid .eq. 1 .and. iocord .ne. 0) then
               if (altc(1:3) .ne. 'tec') then
                  write (fstring, 335) iocord, dls, nspeci, dls,
     &                 nspeci, dls
               else
                  if (icall .eq. 1) then
                     write (fstring, 106) iocord, nspeci
                  else
                     write (fstring, 333) itotal2, dls
                  end if
               end if
            end if             
         
         
         
         
         do ii = istart, iend
c     Node loop          
            string = ''
            if (iozone .ne. 0) then
               if(zone_saved) then
                i = ncord(ii)
               else
                i = ii
                if (izone_surf_nodes(i).ne.idz) goto 199
               endif
            else
              i = ii
            end if
c     Node number will be written first for avs and sur files
c gaz 040517 this is where gdkm number is changed            
            if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'sur') then
               if (ifdual .eq. 0) then
c------------------------------------
c  Shaoping change  10-20-2017
c                 write(string, 100) i
                  write(string, 109) i 
c------------------------------------
               else if (iogdkm .eq. 1.and.iogdkmblank.ne.0) then
c------------------------------------
c  Shaoping change  10-20-2017
c                 write(string, 100) i
                  write(string, 109) i 
c------------------------------------
               else
c------------------------------------
c  Shaoping change  10-20-2017
c                 write(string, 100) i - neq                 
                  write(string, 109) i - neq             
c------------------------------------
               end if
c ic1 positions the column for printout               
               ic1 = 11
            else
               ic1 = 1
            end if
            if (iocord .ne. 0) then
c Only output coordinates that are used
               if (altc(1:3) .eq. 'tec' .and. icall .ne. 1) then
c     Coordinates will only be output in the first file for tecplot
c     (do nothing)
               else 
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
c start endif *****
            if (altc(1:3) .eq. 'tec') then
               if(ifdual.ne.0.and. iogdkm .eq. 1.
     &                        and. iogdkmblank .eq. 1) then
c gaz  040817 if blanking use primary grid node number (and blank variable) 
c gaz identify secondary node variable (i_sec)                   
                   write(vstring, 105) dls(1:k), i
                   i_wrt = i
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
                  i_wrt = i
               else 
                  write(vstring, 105) dls(1:k), i - neq
                  i_wrt = i -neq
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
c end endif *****
           
            

         if (altc(1:3) .eq. 'tec' .and. iocord .ne. 0) then
                  if (icall .eq. 1 .and. iozid .eq. 0) then
                   if(i.ne.iblanking_value) then
                    write(lu, fstring) (corz(i_wrt,j), j =icord1,icord2,
     +                    icord3), i_wrt, (min(maxc, max(minc,
     +                    antmp(i+add_dual,n))), n=1,nspeci)
                   else
                    write(lu, fstring) (corz(i_wrt,j), j =icord1,icord2,
     +                    icord3), i_wrt, (iblanking_value, n=1,nspeci)
                   endif
                  else if (icall .eq. 1 .and. iozid .eq. 1) then
                   if(i.ne.iblanking_value) then
                    write(lu, fstring) (corz(i_wrt,j), j =icord1,icord2,
     +                   icord3), i_wrt, izonef(i), (min(maxc, max(minc,
     +                   antmp(i+add_dual,n))), n=1,nspeci)
                   else
                    write(lu, fstring) (corz(i_wrt,j), j =icord1,icord2,
     +                    icord3), i_wrt, izonef(i),
     +                     (iblanking_value, n=1,nspeci) 
                   endif
                  else
c gaz debug 072819  don't print zone for icall gt 1                    
                     if(i.ne.iblanking_value) then
                      write(lu, fstring) i_wrt, 
     +                 (min(maxc, max(minc,antmp(i+add_dual,n))), 
     +                 n=1,nspeci)
                     else
                      write(lu, fstring) i_wrt, 
     +                    (iblanking_value, n=1,nspeci)       
                     endif
                  end if


               else if (altc(1:3) .eq. 'tec' .and. iozid .ne. 0) then
                  if (icall .eq. 1 .or. iogrid .eq. 1) then
                    if(i.ne.iblanking_value) then
                     write(lu, fstring) i_wrt, izonef(i), (min(maxc,
     +                    max(minc, antmp(i+add_dual,n))), n=1,nspeci)
                    else
                     write(lu, fstring) i_wrt, izonef(i), 
     +                    (iblanking_value, n=1,nspeci) 
                    endif
                  else
                    if(i.ne.iblanking_value) then  
                     write(lu, fstring) i_wrt, (min(maxc, max(minc,
     +                    antmp(i+add_dual,n))), n=1,nspeci)
                    else
                     write(lu, fstring) i_wrt, 
     +                   (iblanking_value, n=1,nspeci)        
                    endif
                  end if
                  
               else if (iocord .ne. 0) then
                  if (iozid .eq. 0) then
                   if(i.ne.iblanking_value) then
                    write(lu, fstring) i_wrt,(corz(i_wrt,j), j = icord1,
     +                    icord2, icord3), (min(maxc, max(minc,
     +                    antmp(i+add_dual,n))), n=1,nspeci)
                   else 
                    write(lu, fstring) i_wrt,(corz(i_wrt,j), j = icord1,
     +                    icord2, icord3), (iblanking_value, n=1,nspeci)
                   endif
                  else
                   if(i.ne.iblanking_value) then                      
                    write(lu, fstring) i_wrt, (corz(i_wrt,j), j =icord1,
     +                    icord2, icord3), izonef(i), (min(maxc, 
     +                    max(minc, antmp(i+add_dual,n))), n=1,nspeci)
                   else
                    write(lu, fstring) i_wrt, (corz(i_wrt,j), j =icord1,
     +                    icord2, icord3), izonef(i), 
     +                    (iblanking_value, n=1,nspeci)
                   endif
                  end if

               else
c-------------------------------
c     Shaoping temporary  fix                  
                 i_wrt =i
c----------------------------------------
                  if (iozid .eq. 0) then
                    if(i.ne.iblanking_value) then  
                     write(lu, fstring) i_wrt,(min(maxc, max(minc,
     +                    antmp(i+add_dual,n))), n=1,nspeci)
                    else
                     write(lu, fstring) i_wrt,
     +                    (iblanking_value, n=1,nspeci)           
                    endif
                  else
                    if(i.ne.iblanking_value) then  
                     write(lu, fstring) i_wrt, izonef(i),
     +                 (min(maxc, max(minc, antmp(i+add_dual,n))),
     +                   n=1,nspeci)
                    else
                     write(lu, fstring) i_wrt, izonef(i),
     +                    (iblanking_value, n=1,nspeci)           
                    endif                     
                  end if
               end if
 199   enddo
            if (altc(1:3) .eq. 'sur') close (lu)

         
         deallocate (title)

c---------------------------------------------------
c     PHS 4/28/2000    end of changes for xpress format
c     PHS 8/10/2000    more changes for xpress below. . .
c     for when rxn_flag NE zero
c---------------------------------------------------
c  middle ifblock 2
      else                      ! IF rxn_flag  NE zero

         fstring = ''
         if(ncplx.ne.0)then
            itotal2 = ncpntprt*2+nimmprt+ncplxprt+nvapprt+iocord+iozid
         else
            itotal2 = ncpntprt+nimmprt+nvapprt+iocord+iozid
         endif

         allocate (title(itotal2))
         do i = 1,itotal2
            t1(i)=1
         enddo
c=======================================================
c     Titles for the rxn stuff  (like Np[aq] ,(Moles/kg H20))
c=======================================================
         do i = 1, iocord
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(i) = trim(cordname(i)) // ', (m)'
            else
               title(i) = trim(cordname(i))
            end if
         end do

         if (iozid .eq. 1) then
            title(i) = 'Zone'
         end if


         irxn_title = iocord + iozid         
         if(altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x')
     2        write(lu,'(i3, 100(1x, i1))')itotal2,(t1(i),i=1,itotal2)


c------CHU 12/02/2015  change so when rxn on, anv output correctly in con files
         do i = 1,ncpntprt
            ic = cpntprt(i)
            irxn_title = irxn_title + 1
            if (icns(i) .eq. -2) then
              do j=1, nspeci
                   do k = 1, n0
                      an_dum(k+add_dual,j) = anv_dum(k+add_dual,j)
                   end do
              end do
               title(irxn_title)= trim(dual_char)
     &              //trim(cpntnam(ic))//' (Moles/kg Vapor)'

            else

            	if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               		title(irxn_title) = trim(dual_char)//trim(cpntnam(ic))
               		write(lu,190) trim(title(irxn_title))
            	else
               		title(irxn_title)= trim(dual_char)
     &              //trim(cpntnam(ic))//' (Moles/kg H20)'
            	end if
            end if
         enddo
c------CHU 12/02/2015  change so when rxn on, anv output correctly in con files

         do i = 1,nimmprt
            im = immprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(immnam(im))
               write(lu,200) trim(title(irxn_title))
            else
               title(irxn_title)= trim(dual_char)
     &              //trim(immnam(im))//' (Moles/kg Rock)'
            end if
         enddo

         do i = 1,nvapprt
            iv = vapprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(vapnam(iv))
               write(lu,210) trim(title(irxn_title))
            else
               title(irxn_title)= trim(dual_char)
     &              //trim(vapnam(iv))//' (Moles/kg Vapor)'
            end if
         enddo

         if(ncplx.ne.0)then
            do i = 1,ncpntprt
               irxn_title = irxn_title + 1
               ic = cpntprt(i)
               if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
                  title(irxn_title) = trim(dual_char)//'Free Ion '
     &                 //trim(cpntnam(ic))
                  write(lu,220) trim(title(irxn_title))
               else
                  title(irxn_title)=trim(dual_char) //'Free Ion '
     &                 //trim(cpntnam(ic))//' (Moles/kg H20)'
               end if
            enddo
         endif

         do i = 101, ncplxprt+100
            ix = cplxprt(i)
            irxn_title = irxn_title + 1
            if (altc(1:3) .eq. 'avs' .and. altc(4:4) .ne. 'x') then
               title(irxn_title) = trim(dual_char)//trim(cplxnam(ix))
               write(lu,190) trim(title(irxn_title))
            else
               title(irxn_title)=trim(dual_char)
     &              //trim(cplxnam(ix))//' (Moles/kg H20)'
            end if
         enddo
         
         fstring = ''
         if(altc(1:4).eq.'avsx') then
            write(fstring, 664) itotal2
            write(lu, fstring) days,(trim(title(i)),i=1,itotal2)

         else if (altc(1:3) .eq. 'tec') then
            if (icall .eq. 1 .or. iogrid .eq. 1) then
               write (fstring, 99) itotal2+1
               write(lu, 98) verno, jdate, jtime, trim(wdd)
c gaz 071022 removed 'SOLUTION'               
c               write(lu, 100)
               write(lu, fstring) 'VARIABLES = ', (trim(title(i)), 
     &              i=1,iocord), 'node', (trim(title(i)), 
     &              i=iocord+1,itotal2)
            end if
            if (iozone .ne. 0) write (lu, 97) trim(timec_string)
         end if

c=================================================
c     Write out node info to the  _con file
c=================================================

c start do loop 1  
c         istart = 1
c         iend = 1
         do ii = istart, iend
c     Node loop          
            string = ''
            if (iozone .ne. 0) then
               if(zone_saved) then
                i = ncord(ii)
               else
                i = ii
                if (izone_surf_nodes(i).ne.idz) goto 299
               endif
            else
              i = ii
            end if
            if (altc(1:3) .eq. 'tec') then
               if (iozone .eq. 0 .or. iogrid .eq. 1) then
c gaz writes   every ii should be ii = start  
                if(ii.eq.istart) then
c gaz 071122 tecstring duplicates other strings                    
                 if(icall.eq.1) then
                  write (lu, 94) trim(timec_string), 
     &                 trim(gridstring), trim(times_string) 
                 else
                  write (lu, 94) trim(timec_string), 
     &                 trim(gridstring), trim(sharestring)                  
                 endif
                 endif
               else
                  if (icall .gt. 1 .and. iozone .ne. 0) then
                     write (tecstring, 125) trim(sharestring), iz
                  end if
                  write (lu, 95) idz, trim(tecstring)
               end if
               if (icall .eq. 1 .and. iz .eq. iendz) then
                  if (iogeo .eq. 1) then
                     tecstring = trim(tecstring) // trim(string)
                  else if (iocord .ne. 0 .or. iozid .ne. 0) then
                     write (tecstring, 125) trim(string), iz
                  end if   
               end if
            else if (altc(1:3) .eq. 'sur') then
               call namefile2(icall,lu,ioformat,tailstring,idz)
               write (fstring, 96) itotal2
               write(lu, fstring) (trim(title(i)), i=1,itotal2)
            end if

c--------------------------------------------
c  Shaoping add  10-20-2017
c           do in = 1,neq
c gaz 063018
            do in = i,i
c-------------------------------------------
               j=0
               if (iocord .ne. 0) then
                  do i = icord1, icord2, icord3
                     j = j+1
                     write_array(j)=corz(in,i)
                  end do
               end if
               do i = 1,ncpntprt
                  ic = cpntprt(i)
                  j=j+1
                  write_array(j)=an_dum(in+add_dual,pcpnt(ic))
               enddo
               do i = 1,nimmprt
                  im = immprt(i)
                  j=j+1
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,an_dum(in+add_dual,pimm(im))))
               enddo
               do i = 1,nvapprt
                  iv = vapprt(i)
                  j=j+1            
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,an_dum(in+add_dual,pvap(iv))))
               enddo
               if(ncplx.ne.0)then
                  do i = 1,ncpntprt
                     ic = cpntprt(i)
                     j=j+1
                     write_array(j)=min(1.0d+40,
     2                    max(1.0d-90,cpntsv(ic,in+add_dual)))
                  enddo
               endif
               do i = 101,ncplxprt+100
                  ix = cplxprt(i)
                  j=j+1
                  call cplxcalc(in+add_dual,ix,complex_conc)
                  write_array(j)=min(1.0d+40,
     2                 max(1.0d-90,complex_conc))
               enddo
               
               write_total = j
               fstring = ''
               
               if (iozid .eq. 0 .and. iocord .eq. 0) then
                  write (fstring, 333) write_total, dls
               else if (iozid .eq. 0 .and. iocord .ne. 0) then
                  if (altc(1:3) .ne. 'tec') then
                     write (fstring, 333) write_total, dls
                  else
                     if (icall .eq. 1) then
                        write (fstring, 104) iocord, write_total-iocord
                     else
                        write (fstring, 333) write_total - iocord, dls
                     end if
                  end if
               else if (iozid .eq. 1 .and. iocord .eq. 0) then
                  if (altc(1:3) .ne. 'tec' .or. icall .eq. 1 .or.
     &                 iogrid .eq. 1) then
                     write (fstring, 334) dls, write_total, dls
                  else
                     write (fstring, 333) write_total, dls
                  end if
               else if (iozid .eq. 1 .and. iocord .ne. 0) then
                  if (altc(1:3) .ne. 'tec') then
                     write (fstring, 335) iocord, dls, dls, 
     &                    write_total - iocord, dls
                  else
                     if (icall .eq. 1) then
                        write (fstring, 106) iocord, write_total-iocord
                     else
                        write (fstring, 333) write_total - iocord, dls
                     end if
                  end if
               end if
       
               if (altc(1:3) .eq. 'tec' .and. iocord .ne. 0) then
                  if (icall .eq. 1 .and. iozid .eq. 0) then
                     write(lu,fstring) (write_array(i),i=1,iocord), 
     &                    in, (write_array(i),i=iocord+1,write_total)
                  else if (icall .eq. 1 .and. iozid .eq. 1) then
                     write(lu,fstring) (write_array(i),i=1,iocord), in,
     &                    izonef(in), (write_array(i),
     &                    i=iocord+1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i),
     &                    i=iocord+1,write_total)
                  end if

               else if (altc(1:3) .eq. 'tec' .and. iozid .ne. 0) then
                  if (icall .eq. 1 .or. iogrid .eq. 1) then
                     write(lu,fstring) in, izonef(in), 
     &                    (write_array(i),i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i),
     &                    i=iocord+1,write_total)
                  end if

               else if (iocord .ne. 0) then
                  if (iozid .eq. 0) then
                     write(lu,fstring) in, (write_array(i),
     &                    i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i), 
     &                    i=1,iocord), izonef(in), 
     &                    (write_array(i),i=iocord+1,write_total)
                  end if
               else
                  if (iozid .ne. 0) then
                     write(lu,fstring) in, izonef(in), 
     &                    (write_array(i),i=1,write_total)
                  else
                     write(lu,fstring) in, (write_array(i), 
     &                    i=1,write_total)
                  end if
               end if
            enddo
 299   enddo
c  end ifblock 2
          endif
c end do loop 1
c          endif
c end if block 1           
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
         if (altc(1:3) .eq. 'sur') close (lu)
      enddo
      if(zone_saved) then 
       if(.not.allocated(contour_conc_files)) 
     &     allocate(contour_conc_files(100))
        icconc = icconc + 1
        inquire(unit=lu,name = file_flux) 
        contour_conc_files(icconc) = file_flux
      endif
c gaz 073017 don't think flush(lu) is needed      
      call flush(lu)
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
            il = open_file(geoname,'old')
c     tec geometry file has an initial line that starts TITLE
c gaz 050522 adding tec and other if blocks            
           if(altc(1:3).eq.'tec') then
            read(il,*) 
            read(il,*) 
            read(il,*) 
            read(il,*) 
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
            read(il,*) i
            if (i .ne. neq_primary) backspace il
            do i = 1, neq
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
         call structured(4)
         il = open_file('fdm_elem.macro','old')
         read(il,*) 
         read(il,*)  ns_in , nei_in
         allocate (nelm2(ns_in))
         do i = 1, nei_in
            read (il,*) i1, (nelm2(j), j=1,ns_in)
            write(lu, '(8(i8))') (nelm2(j), j=1,ns_in)
         end do
         deallocate(nelm2)
         close (il)
         end if   

      if (altc(1:3) .ne. 'sur') close (lu)
      iocord = iocord_tmp

c 94   format('ZONE T = "Simulation time ',1p,g16.9,' days"', a)
 94   format('ZONE T = ', a, a, a, a)
 95   format('ZONE T = "',i4.4,'"', a)
 96   format("('node', ", i3, "(', ', a))")
c 97   format('TEXT T = "Simulation time ',1p,g16.9,' days"')
 97   format('TEXT T = ', a)
 98   format ('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
 99   format ('(a11, ', i3, "('",' "',"', a, '",'"',"'))")
 100  format ('FILETYPE = "SOLUTION"')
 109  format ('FILETYPE = "SOLUTION"', i6)
 104  format ('(1x, ', i3, '(g16.9, 1x), i10.10, ', i3, '(1x, g16.9))')

 125  format(', VARSHARELIST = ([', a,'] = ', i4, ')')
 135  format(', N = ', i8, ', E = ', i8, ', DATAPACKING = POINT',
     &     ', ZONETYPE = ', a)
 140  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')
 664  format("('nodes at ',e16.9,' days ',",i3,"(' : ', a))")
 190  format(a,', (Moles/kg H20)')
 200  format(a,', (Moles/kg Rock)')
 210  format(a,', (Moles/kg Vapor)')
 220  format('Free Ion ',a,', (Moles/kg H20)')
C     No zid
 333  format("(1x, i10.10, ", i3, "('", a, "', g16.9))")
C     Zid no coordinates
 334  format("(1x, i10.10, '", a, "', i4, ", i3, "('", a, "', g16.9))")
C     Coordinates and zid
 335  format("(1x, i10.10, ", i3, "('", a, "', g16.9), '", a, "', i4,",
     &     i3, "('", a, "', g16.9))")
 130  format('ZONE T =', a, a, a, a) 
 131  format('ZONE ',a,',',' T =', a, a, a, a)     
 303  format(', VARSHARELIST = ([', a,'] = ', i4, ')')
 304  format('FILETYPE = "SOLUTION"')
 305  format(', VARSHARELIST = ([', a,'] = ', i4, '), ',
     &     'CONNECTIVITYSHAREZONE = 1')   
 118  format('TEXT T = ', a)
 120  format('ZONE T = "',i4.4, '"', a)     
 110  format(a, g16.9)
 115  format(a, i4)
 106  format ('(1x, ', i3, '(g16.9, 1x), i10.10, 1x, i4, ', i3, 
     &     '(1x, g16.9))')
 105  format(a, i10.10)     
c      deallocate(anv_dum,an_dum,antmp)
      if(allocated(anv_dum)) deallocate(anv_dum)
      if(allocated(an_dum)) deallocate(an_dum)
      if(allocated(antmp)) deallocate(antmp)
      return
      end
