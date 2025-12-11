      subroutine avs_io(inj)
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
CD1 Produce FEHM output in AVS UCD format
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
CD2 $Log:   /pvcs.config/fehm90/src/avs_io.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:18   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:54:18   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:48   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:54 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.10   Wed Feb 14 10:47:16 1996   zvd
CD2 Modified requirements
CD2 
CD2    Rev 1.9   Thu Feb 08 13:35:28 1996   llt
CD2 changed ifdual = 4 to ifdual = 1
CD2 
CD2    Rev 1.8   Mon Jan 29 13:10:54 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.7   Mon Jan 29 10:18:12 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.6   09/08/95 11:55:42   zvd
CD2 Corrected unit designator from ierr to iptty for terminal output.
CD2 Changed use of avs_prefix to file_prefix so a single file prefix routine is used.
CD2 
CD2    Rev 1.5   09/05/95 13:24:40   llt
CD2 changed "" string comparisons to ' ', for Cray
CD2 
CD2    Rev 1.4   01/28/95 14:19:52   llt
CD2 modified for the revised reactive transport module
CD2 
CD2    Rev 1.3   01/26/95 14:21:04   tam
CD2 added icnl to call avs_write_cord, to allow 0's in z for 2-dim.
CD2 
CD2    Rev 1.2   12/12/94 16:20:54   tam
CD2 changed include from ' ' to " "
CD2 so includes can be in current or a defined directory
CD2 
CD2    Rev 1.1   11/29/94 13:42:12   llt
CD2 Removed protying, so could run under cc compile, instead of c compiler. 
CD2 Used #ifdef to determine format for routine name, depending on machine.
CD2 (Changes made by tam.)
CD2 
CD2    Rev 1.0   08/23/94 15:32:40   llt
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
c
c       inj = 0 on first call, this will produce header files
c               and grid coordinate and connectivity information.
c       inj > 0 output node scalar and vector data
c
c       inj < 0 no action
c
c***references (none)
c***routine called  namefile1, namefile2, write_avs_node_con, 
C      write_avs_node_s,
c      write_avs_node_v, avs_write_cord, avs_write_struc,
c      write_avs_node_mat, write_avs_ucd_header
      use avsio
      use comai
      use combi
      use comchem
      use comdi
      use comdti, only : n0
      use comriv, only : iriver
      use comrxni
      use comxi
      use commeth
      use comwt, only : isconwt, col
      use davidi

      implicit none

      character*110 logname
      character*14 head_mat
      character*14 head_mat_dual
      character*14 head_sca
      character*14 head_sca_dual
      character*14 head_vec
      character*14 head_vec_dual
      character*14 head_con
      character*14 head_con_dual
      character*14 head_hf
      character*14 geom_tail
      character*14 material_tail
      character*14 material_dual_tail
      character*14 scalar_tail
      character*14 duals_tail
      character*14 vector_tail
      character*14 dualv_tail
      character*14 concen_tail
      character*14 dualcon_tail
      character*14 heatflux_tail
      character*14 log_tail
      character*14 wt_tail
      character*19 tmp_tail

      integer itotal2, neq_tmp, neq_gdkm2
      integer inj, lu, lu_log, icall, num_cdata, num_mdata, io_err
      integer groot, len_scalar, len_vec, len_material, length
      integer ifdual, n1, ktmp, ktmp_head, ltmp, mout, iocord_temp
      integer nscalar, nscalar_dual, nvector, nvector_dual
      integer nconcen, nconcen_dual, nmaterial, nmaterial_dual
      integer nheatflux
c gaz 060722 
      integer file_size_bytes, lu_test
C
      real*8, allocatable :: dum(:)
      logical null1, exists, opnd
      character*10 response
      character*20 dum_geo
C
      save lu_log

      iocord_temp = 0
      data head_mat /'mat_head'/
      data head_sca /'sca_head'/
      data head_vec /'vec_head'/
      data head_con /'con_head'/
      data head_hf /'hf_head'/
      data geom_tail /'_geo'/
      data material_tail /'mat_node'/
      data scalar_tail /'_sca_node'/
      data vector_tail /'_vec_node'/
      data concen_tail /'_con_node'/
      data heatflux_tail /'_hf_node'/
      data log_tail /'avs_log'/
      data wt_tail /'_wt'/
      data icall /1/
      data num_cdata / 0 /, num_mdata / 0 /
      data io_err /0/

      if (iogdkm .eq. 1) then
         head_mat_dual = 'mat_gdkm_head'
         head_sca_dual = 'sca_gdkm_head'
         head_vec_dual = 'vec_gdkm_head'
         head_con_dual = 'con_gdkm_head'
         material_dual_tail = 'mat_gdkm_node'
         duals_tail  = '_sca_gdkm_node'
         dualv_tail  = '_vec_gdkm_node'
         dualcon_tail  = '_con_gdkm_node'
         neq_gdkm2 = neq - neq_primary
         neq_tmp = neq
         if (icnl .eq. 0) then
            iocord_temp = 3
         else
            iocord_temp = 2
         end if
      else
         head_mat_dual = 'mat_dual_head'
         head_sca_dual = 'sca_dual_head'
         head_vec_dual = 'vec_dual_head'
         head_con_dual = 'con_dual_head'
         material_dual_tail = 'mat_dual_node'
         duals_tail = '_sca_dual_node'
         dualv_tail = '_vec_dual_node'
         dualcon_tail = '_con_dual_node'
         neq_tmp = neq_primary
      end if
C
      if(ico2.lt.0) then
         allocate(dum(neq))
         dum=0.0d00
      endif

      if (inj .eq. 0) then
C
C     Substitute file_prefix to find root of contour file name or if not
C     present, output or input file name
         if (null1(root_name)) then
            if (nmfil(10) .ne. nmfily(3) .and. nmfil(10) .ne. ' ') then
               call file_prefix(nmfil(10), iaroot)
               if (iaroot .gt. 94) iaroot = 94
               avs_root(1:iaroot) = nmfil(10)(1:iaroot)
            else if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ')
     &              then
               call file_prefix(nmfil(5), iaroot)
               if (iaroot .gt. 94) iaroot = 94
               avs_root(1:iaroot) = nmfil(5)(1:iaroot)
            else 
               if (nmfil(2)(1:1) .eq. ' ' ) then
                  write (ierr, *) 'FILE ERROR: nmfil2 file: ', nmfil(2),
     .                 ' unable to determine contour file prefix'
                  stop
               else
                  call file_prefix(nmfil(2), iaroot)
                  if (iaroot .gt. 94) iaroot = 94
                  avs_root(1:iaroot) = nmfil(2)(1:iaroot)
               end if
            endif
         else
            iaroot = len_trim (root_name)
            if (iaroot .gt. 94) iaroot = 94
            avs_root(1:iaroot) = root_name(1:iaroot)
         end if
         iaroot = iaroot + 1
         avs_root(iaroot:iaroot) = '.'
      end if

C     figure out how many  components in each file, and
C     the length of each component.
      
      len_scalar = 1
      len_vec = 3
      len_material = 1

      if(rxn_flag.eq.0)then
         itotal2 = nspeci
      else
         if(ncplx.ne.0)then
            itotal2=ncpntprt*2+nimmprt+ncplxprt+nvapprt
         else
            itotal2=ncpntprt+nimmprt+nvapprt
         endif
      endif      
      nscalar = (ioliquid+iovapor+iovapor*iadif)*iopressure
     >     +iotemperature+iosaturation+iohead+ioporosity
     >     +iosource+(ioliquid+iovapor)*iodensity
     >     +iopermeability*3+iozid+(ioliquid+iovapor)*ioflx
     >     +iocapillary+ioco2*7+iodisp*3+iostress+iostrain
      if (nscalar .ne. 0) then
         nscalar_dual   = nscalar * iodual + nscalar * iogdkm + 
     &        iocord_temp
c gaz 031022 added ioscalar to indicate if scalar variables are available
        ioscalar = 1 
      else
          ioscalar = 0
          nscalar_dual   = 0
      end if
      nvector        = (ioliquid+iovapor)*iovelocity
      if (nvector .ne. 0) then
         nvector_dual   = nvector * iodual
c     No velocities for gdpm or gdkm
         if (iogdkm .eq. 1) then
            write (ierr, *) 'Velocities are not computed for gdkm nodes'
            write (ierr, *) ' and will not be output'
         end if
      else
         nvector_dual = 0
      end if
      if (idoff .eq. -1) then
c     Heat conduction only
         nheatflux = ioheatflux
      else
c     Advection and conduction
         nheatflux = 2 * ioheatflux
      end if
      if (iomaterial .ne. 0) then
         mout = 0
         if (idoff .ne. -1) then
! Permeability will be written
            mout = mout + 3
         end if     
         if (ico2 .gt. 0 .or. ice .ne. 0) then
! Conductivity will be written
            mout = mout + 3
         end if
! Porosity and specific heat will be written
         mout = mout + 2
         if (irdof .ne. 13) then
! Capillary pressure will be written
            mout = mout + 1
         end if
         if (rlp_flag .eq. 1) then
! rlp and cap model flags will be written if rlp_flag .eq. 1
            mout = mout + 2
         end if
         nmaterial    = iomaterial * mout
         nmaterial_dual  = nmaterial * iodual + nmaterial * iogdkm
      else
         nmaterial    = iomaterial
         nmaterial_dual  = iodual * iomaterial
      end if
      nconcen        = itotal2
      if (itotal2 .ne. 0) then
         nconcen_dual   = itotal2 * iodual + itotal2 * iogdkm + 
     &        iocord_temp
      else
         nconcen_dual   = 0
      end if
c special case for hydrate (gaz 10-29-03)
c
      if(iohyd.ne.0) nscalar = 4
      
      if (inj .eq. 0) then
         icall = 1
         if (iogdkm .eq. 1) neq_tmp = neq_gdkm2
C     first time through open a file for logging of avs output
C     Open with ascii format, not binary
         
         call namefile1(lu_log,2,avs_root,log_tail,
     &        iaroot,ierr)
         write(lu_log,300) verno, jdate
 300     format('# ',a30,3x,a11)
         write(lu_log,305)
 305     format('#   LOG AVS OUTPUT')
         write(lu_log,310)wdd
 310     format('# ',a80)
 
         gridstring = ''
         if (iogeo .eq. 1 .or. iogrid .eq. 1) then
! Find geometry file name and write if it does not exist
            geoname = ' '
            lu = 9999
! Use grid file name root as default (if different from input file)
            if (nmfil(3) .ne. nmfil(2) .and. nmfil(3) .ne. ' ') then
               call file_prefix(nmfil(3), groot)
               if (groot .gt. 100) groot = 100
               geoname(1:groot) = nmfil(3)(1:groot)
               if (iogeo .eq. 1) then
                  geoname(groot+1:groot+4) = ".geo"
               else
                  geoname(groot+1:groot+9) = "_grid.dat"
               end if
               if (geoname .eq. nmfil(3)) then
! In case coordinate file has ".geo" or "_grid.dat" as its suffix
                  length = len_trim(geoname)
                  geoname(length+1:length+4) = "_" // "altc(1:3)"
               end if
               inquire(file = geoname, exist = exists)
! If the file exists, we do not need to do anything else 
               if (.not. exists) then
! Open the file to be written
                  response = ''
                  open (lu, file = geoname, iostat = io_err, err=400)
                  inquire (file = geoname, write = response)
! In case we can't write to directory with grid file
 400              if (response .eq. 'NO' .or. io_err .ne. 0) then
                     if (response .eq. 'NO') close (lu)
                     if (iout .ne. 0) then 
                        write (iout, *) "Can't write to file:", geoname
                     end if
                     if (iptty .ne. 0) then 
                        write (iptty, *) "Can't write to file:", geoname
                     end if
                     geoname = ' '
                  end if
               else
! Check to make sure we can open the existing file
                  open (lu, file = geoname, iostat = io_err, err=410)
               end if
            end if
            
 410        inquire (file = geoname, OPENED=opnd)
c            inquire(file = geoname, exist = exists)
c            inquire(file = geoname, recl = file_size_bytes)
c            inquire(file = geoname, number = lu_test)
            if (.not. opnd .or. geoname .eq. ' ') then
! Use root name determined above
               geoname = ''
               geoname(1:iaroot) = avs_root(1:iaroot)
               if (iogeo .eq. 1) then
                  geoname(iaroot+1:iaroot+3) = "geo"
               else if (iogrid .eq. 1) then
                  geoname(iaroot:iaroot+9) = "_grid.dat"
               end if
               if (geoname .eq. nmfil(3)) then
! In case coordinate file has ".geo" as its suffix
                  length = len_trim(geoname)
                  geoname(length + 1:length + 4) = "_" // altc(1:3)
               end if
               inquire(file = geoname, exist = exists)
! If the file exists, we do not need to do anything else 
               if (.not. exists) then
! Open the file to be written
                  inquire(file = geoname, OPENED=opnd)
                  if(.not.opnd) open (lu, file = geoname)
                   inquire(file = geoname, OPENED=opnd)
                   inquire(file = geoname, exist = exists)
                   inquire(file = geoname, size = file_size_bytes)
                   inquire(file = geoname, number = lu_test)
               else if(file_size_bytes.gt.0) then
c gaz 050922 check if geo file matches contour option   
                  if(.not.opnd) open (lu, file = geoname)
                  read(lu,'(a7)') dum_geo(1:7)
                  if(dum_geo(1:7).eq.'TITLE ='.and.
     &               altc(1:3).ne.'tec') then
                   close(lu,dispose = 'delete')
                   open (lu, file = geoname)
                   exists = .false.
                  else if(dum_geo(1:7).ne.'TITLE ='.and.
     &               altc(1:3).eq.'tec') then
                   close(lu,dispose = 'delete')
                   open (lu, file = geoname)
                   exists = .false.
                  endif
               end if
               end if
c gaz 050922 check for mismatch  with tec and avs geo files
            if (exists.and.file_size_bytes.gt.0) then 
                   read(lu,'(a7)') dum_geo(1:7)
                  if(dum_geo(1:7).eq.'TITLE ='.and.
     &               altc(1:3).ne.'tec') then
                   close(lu,dispose = 'delete')
                   open (lu, file = geoname)
                   exists = .false.
                  else if(dum_geo(1:7).ne.'TITLE ='.and.
     &               altc(1:3).eq.'tec') then
                   close(lu,dispose = 'delete')
                   open (lu, file = geoname)
                   exists = .false.
                  else
                   backspace lu
                  endif             
            endif    
            if (exists.and.file_size_bytes.gt.0) then
               if (iout .ne. 0) write (iout, *) 
     &              "Using existing geometry file: ", geoname
               if (iptty .ne. 0) write (iptty, *) 
     &              "Using existing geometry file: ", geoname
c gaz 041822  050922 commented out             
c               call tec_write_grid(lu)
! We will reopen the file when it needs to be read
               if (opnd)  close (lu)
c gaz 060822 if opened fill file               
c            else if (exists) then
          else
! Write geometry file
               if(.not.opnd) open (lu, file = geoname)
               
               if (iout .ne. 0) write (iout, *) 
     &              "Geometry written to file: ", geoname
               if (iptty .ne. 0) write (iptty, *) 
     &              "Geometry written to file: ", geoname
              
               if (ioformat .eq. 1) then
C No binary option
               else
                  if (iogeo .eq. 1.and.altc(1:3).eq.'avs') then
                     call avs_write_cord(corz(1,1),
     1                    corz(1,2),
     2                    corz(1,3),
     3                    neq_primary,lu,ioformat,icnl,ierr)

                     call avs_write_struc(nelm,ns,icnl,nei,
     1                    corz(1,1),
     2                    corz(1,2),
     3                    corz(1,3),
     4                    neq_primary,lu,ioformat,ierr)
                     close(lu)
c 050922 gaz commented out                     call tec_write_grid(0)
                  elseif (iogeo .eq. 1.and.altc(1:3).eq.'tec') then
                     call tec_write_grid(lu)
                     close(lu)                     
                  end if
               endif
            end if
         end if
         
C     output material properties
c====================================================================
c   added altc to pass to write_avs_node_mat so avsx can get zone info
c    PHS 8/11/2000
c=====================================================================
         if (nmaterial .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_mat,
     &           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' )
     &              call write_avs_ucd_header(lu,verno,jdate,wdd,
     1              neq_primary,nei,nmaterial,num_cdata,num_mdata)
               
               if (altc(1:3) .eq. 'tec') then
                  tmp_tail = trim(material_tail) // '.dat'
               else if (altc(1:3) .eq. 'sur') then
                  tmp_tail = trim(material_tail) // '.csv'
               else if (altc(1:4) .eq. 'avsx') then
                  tmp_tail = trim(material_tail) // '.avsx'
               else 
                  tmp_tail = trim(material_tail) // '.avs'
               end if
               call namefile1(lu,ioformat,avs_root,tmp_tail,
     &              iaroot,ierr)
               ifdual = 0
c gaz 042921 first step to make mat file like scalar file               
c               call write_avs_node_mat (lu,ifdual,nmaterial)
              call write_avs_node_mat_s(icall, neq_primary, nscalar, lu,
     &                ifdual, 0, nmaterial)
               close(lu)
            endif
         endif
         
         if (nmaterial_dual .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_mat_dual,
     .           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' )
     &              call write_avs_ucd_header
     1              (lu,verno,jdate,wdd,neq_tmp,nei,nmaterial_dual,
     2              num_cdata,num_mdata)
               n1 = neq_primary + 1
               if (altc(1:3) .eq. 'tec') then
                  tmp_tail = trim(material_dual_tail) // '.dat'
               else if (altc(1:3) .eq. 'sur') then
                  tmp_tail = trim(material_dual_tail) // '.csv'
               else if (altc(1:4) .eq. 'avsx') then
                  tmp_tail = trim(material_dual_tail) // '.avsx'
               else
                  tmp_tail = trim(material_dual_tail) // '.avs'
               end if
               call namefile1(lu,ioformat,avs_root,tmp_tail,
     &              iaroot,ierr)
               ifdual = 1
c gaz 043021                
c               call write_avs_node_mat (lu,ifdual,nmaterial_dual)
c               write(lu,*) 'gdkm material test '
              call write_avs_node_mat_s(icall, neq, nscalar, lu,
     &                ifdual, 0, nmaterial_dual)
               close(lu)
            endif
         endif
         
         if (nscalar .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_sca,iaroot,
     &           ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' )
     &              call write_avs_ucd_header(lu,verno,jdate,wdd,
     &              neq_primary,nei,nscalar,num_cdata,num_mdata)
            endif
         endif
         
         if (nscalar_dual .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_sca_dual,
     .           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' ) 
     &              call write_avs_ucd_header
     1              (lu,verno,jdate,wdd,neq_tmp,nei,nscalar_dual,
     2              num_cdata,num_mdata)
            endif
         endif
         
C     nvector needs to be multiplied by len_vec for the header 
         if (nvector .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_vec,
     &           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
               continue
            else
               if (altc(1:3) .eq. 'avs' ) 
     &              call write_avs_ucd_header(lu, verno,jdate,wdd,
     1              neq_primary,nei,(nvector*len_vec),num_cdata,
     2              num_mdata)
            endif
         endif
         
         if (nvector_dual .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_vec_dual,
     &           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' ) 
     &              call write_avs_ucd_header(lu,verno,jdate,wdd,
     1              neq_tmp,nei,(nvector_dual*len_vec),num_cdata,
     2              num_mdata)
            endif
         endif
         
         if (nconcen .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_con,
     &           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' ) 
     &              call write_avs_ucd_header(lu,verno,jdate,wdd,
     1              neq_primary,nei,nconcen,num_cdata,num_mdata)
            endif
         endif
         
         if (nconcen_dual .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_con_dual,
     &           iaroot,ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' ) 
     &              call write_avs_ucd_header
     1              (lu,verno,jdate,wdd,neq_tmp,
     &              nei,nconcen_dual,num_cdata,num_mdata)
            endif
         endif

         if (nheatflux .ne. 0) then
            if (altc(1:3) .eq. 'avs' ) 
     &           call namefile1(lu,ioformat,avs_root,head_hf,iaroot,
     &           ierr)
            if (ioformat .eq. 1) then
C No binary option
            else
               if (altc(1:3) .eq. 'avs' )
     &              call write_avs_ucd_header(lu,verno,jdate,wdd,
     &              neq_primary,nei,nheatflux*len_vec,num_cdata,
     &              num_mdata)
            endif
         end if

         if (iogdkm .eq. 1) neq_tmp = neq
         
      else if (inj .lt. 0) then
         
C     ignore this call, write everything on the inj .gt. 0 call
         
      else if (inj .gt. 0) then
         
C     write out field information
         
c     use user defined root appended with _v for vector fields
c     use user defined root appended with _s for scalar fields
c     use user defined root appended with _dp for dual perm. fields

         timec_string = ''
         times_string = ''
         if (timec_flag .eq. 1) then
            contour_time = abs(days / 365.25d00)
            time_units = 'years'
         else if (timec_flag .eq. 2) then
            contour_time = abs(days)
            time_units = 'days'
         else if (timec_flag .eq. 3) then
            contour_time = abs(days * 86400.d00)
            time_units = 'seconds'
         else if (timec_flag .eq. 4) then
            contour_time = abs(days * 24.d00)
            time_units = 'hours'
         end if
         write (timec_string, 500) contour_time, trim(time_units)
         
         if (iogrid .ne. 0) write (times_string, 510) contour_time
c If the STRANDID is not zero, the SOLUTIONTIMES must be different
c for each zone in the STRAND
 500     format ('"Simulation time ', 1p, g16.9, a, '"')
 510     format (', STRANDID = 0, SOLUTIONTIME = ', 1p, g16.9)

         if (nscalar .ne. 0) then
            nscalar = nscalar
            if (.not. allocated(head)) then
               allocate (head(1))
               head = 0.
            end if
            if (ioformat .eq. 1) then
C No binary option
            else
               ifdual = 0

c  PHS 5/03/2000   added altc and days to the pass to write_avs_node_s
C zvd 10/23/2007 remove altc, days from call, they are in comai

              if(iohyd.eq.0) then
                 call write_avs_node_s(icall, neq_primary, nscalar, lu,
     &                ifdual, 0)
c     RJP 04/30/2007 added for wellbore nodes below
                 if (iriver .ne. 0) 
     &                call write_avs_node_s(icall, neq_primary, nscalar,
     &                lu, ifdual, iriver)
              else
                 call namefile2(icall,lu,ioformat,scalar_tail,iaroot)
                 call write_avs_node_h(
     1              phi(1),t(1),fracw(1),frachyd(1),
     2              neq_primary,nscalar,lu,
     3              ioliquid,
     4              iovapor,
     5              iopressure,
     6              iotemperature,
     6              iofw, 
     7              iofh,
     8              ifdual)
                 call write_rlp_hyd(
     &                neq_primary,lu,ifdual) 
              endif
           endif
        endif
         
C     ktmp is the index for vapor values
C     ltmp is the index for liquid values
         
        if (nvector .ne. 0) then
           ltmp = n + 1
           if (irdof .eq. 13) then
              ktmp = n + 1
           else
              ktmp = n + n + 1
           end if
           ifdual = 0
           if (ioformat .eq. 1) then
C No binary option
           else
              call write_avs_node_v(icall,pnx(ktmp),pny(ktmp),pnz(ktmp),
     1             pnx(ltmp),pny(ltmp),pnz(ltmp),
     2             neq_primary,nvector,lu,
     3             ifdual)
              close(lu)
           endif
        endif
         
        if (nscalar_dual .ne. 0) then
           ktmp = neq_primary + 1
           ifdual = 1
           if (ioformat .eq. 1) then
C     No binary option
           else

c     PHS 5/03/2000   added altc and days to the pass to write_avs_node_s
C     zvd 10/23/2007 remove altc, days from call, they are in comai
              if (.not. allocated(head)) then
                 allocate (head(1))
                 head = 0.
                 ktmp_head = 1
              else
                 ktmp_head = size (head)
                 if (ktmp .lt. ktmp_head) ktmp_head = ktmp
              end if

              if (iohyd.eq.0) then 
                 call write_avs_node_s(icall, neq_tmp, nscalar,
     &                lu, ifdual, 0)
c     RJP 04/30/2007 added for wellbore nodes below
                 if (iriver .ne. 0) 
     &                call write_avs_node_s(icall, neq_tmp, 
     &                nscalar, lu, ifdual, iriver)

              else
                 call namefile2(icall,lu,ioformat,
     &                duals_tail,iaroot)
                 call write_avs_node_h(
     1                phi(ktmp),t(ktmp),fracw(ktmp),frachyd(ktmp),
     2                neq_tmp,nscalar,lu,
     3                ioliquid,
     4                iovapor,
     5                iopressure,
     6                iotemperature,
     6                iofw, 
     7                iofh,
     8                ifdual)
                 call write_rlp_hyd(
     &                neq_tmp,lu,ifdual) 
                 close(lu)
              endif
           endif
        endif
        
C     ktmp is the index for vapor values
C     ltmp is the index for liquid values
        
C     zvd 10/23/2007 remove altc, days from call, they are in comai
        if (nvector_dual .ne. 0) then
           ltmp = neq_primary + n + 1
c           ltmp = neq_tmp + n + 1
           if (irdof .eq. 13) then
              ktmp = neq_primary + n + 1
c              ktmp = neq_tmp + n + 1
           else
              ktmp = neq_primary + n + n + 1
c              ktmp = neq_tmp + n + n + 1 
           end if
           ifdual = 1
           if (ioformat .eq. 1) then
C     No binary option
           else
              call write_avs_node_v(icall,
     1             pnx(ktmp),pny(ktmp),pnz(ktmp),
     2             pnx(ltmp),pny(ltmp),pnz(ltmp),
     3             neq_tmp,nvector,lu,
     4             ifdual)
              close(lu)
           endif
        endif
        
        if (iowt .ne. 0) then
           call namefile2(icall,lu,ioformat,wt_tail,0)
           isconwt = lu
           if (ifree .ne. 0) then
              if (.not. allocated(col)) then
                 if (.not. allocated(izone_free_nodes)) then
                    allocate(izone_free_nodes(n0))
                    izone_free_nodes=0
                 end if
                 call wtsi_column
              end if
              call airctr(14,0)
           else
              call wtsictr(14)
           end if
           close(lu)
        end if
         
        if ((nconcen .ne. 0) .and. (ioconcentration .ne. 0)) then
            ifdual = 0
            if (ioformat .eq. 1) then
C No binary option
            else

c  PHS 4/27/2000   added altc and days to the pass to write_avs_node_con
C zvd 12/20/2002 pass in ntp array [sometimes only dimensioned (1)]
C zvd 10/23/2007 remove altc, days from call, they are in comai
c               call write_avs_node_con(an,anv,iovapor,npt(2),neq
c gaz 111414
               dum_geo(1:10) = wdd(1:10)
               call write_avs_node_con(icall,npt,neq_primary,
     &              nspeci,lu,ifdual)

            endif
            
            if (iodual .ne. 0 .or. iogdkm .ne. 0) then
               ifdual = 1
               if (ioformat .eq. 1) then
C No binary option
               else

c  PHS 4/27/2000   added altc and days to the pass to write_avs_node_con
c gaz 111414 and 070118
c               call write_avs_node_con(icall,npt,neq_primary,
c     &              nspeci,lu,ifdual)
               call write_avs_node_con(icall,npt,neq,
     &              nspeci,lu,ifdual)               
               endif
               
            endif
         endif

         if (nheatflux .ne. 0) then
! No dual options for now
            ifdual = 0
            if (ioformat .eq. 1) then
C     No binary option
            else
               call write_avs_node_hf(icall,neq,nheatflux,lu,ifdual)
            end if
         end if
         
C     Write output to logfile
         
         if (iptty .ne. 0 ) then
            write (iptty, *) 'Finished writing files for: ',
     .           avs_root(1:iaroot)
            write (iptty, *) 'Number of times called: ',icall
         end if
         if (iout .ne. 0 ) then
            write (iout, *) 'Finished writing files for: ',
     .           avs_root(1:iaroot)
            write (iout, *) 'Number of times called: ',icall
         end if
         if (icall .eq. 1) then
            write(lu_log, 315) trim(time_units)
         end if

         write(lu_log, 320) avs_root(1:iaroot), icall, contour_time
         icall_max = max(icall,icall_max)
         call flush(lu_log)
          
 315     format('# Root filename', 15x, 'Output Time (', a, ')')
        
 320     format(1x, a, i5.5, 5x, g16.9)
         
C     Increment ncall when node output is written
         
         icall = icall + 1     
      endif

      if(allocated(dum)) deallocate(dum)
      
      end
