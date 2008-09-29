      subroutine write_avs_head_s(icall,
     2     nscalar,
     3     lu,
     4     ifdual, 
     5     idz,
     6     iriver2)
!***********************************************************************
!  Copyright, 2006,  The  Regents of the University of California.
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
CD2 $Log:   /pvcs.config/fehm90/src/write_avs_head_s.f_a  $
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

      use avsio
      use comai, only : altc, contim, days, iadif, icnl, jdate, jtime, 
     &     verno, wdd
      use comdi, only : head
      use davidi
      implicit none

      integer maxscalar
      parameter (maxscalar = 22)
      integer neq,nscalar,lu,ifdual
      integer i,j,iolp,iovp,nout,icall,idz,iriver2
      integer size_head, size_pcp, istart, iend, ic1, ic2, length
      character*5 dual_char
      character*14 tailstring
      character*20 tailstring2
      character*20 formstring
      character*40 title(maxscalar+3)
      character*20 units(maxscalar+3)
      character*1200 tstring2
      character*280 string
      character*50 tstring

      data units( 1) /'(MPa)'/,
     &     units( 2) /'(MPa)'/,
     &     units( 3) /'(deg C)'/,
     &     units( 4) /'(no dim)'/,
     &     units( 5) /'(m)'/,
     &     units( 6) /'(MPa)'/,
     &     units( 7) /'(no dim)'/,
     &     units( 8) /'(kg/s)'/,
     &     units( 9) /'(kg/m**3)'/,
     &     units(10) /'(kg/m**3)'/,
     &     units(11) /'(log m**2)'/,
     &     units(12) /'(log m**2)'/,
     &     units(13) /'(log m**2)'/,
     &     units(17) /'(no dim)'/,
     &     units(18) /'(kg/s)'/,
     &     units(19) /'(kg/s)'/,
     &     units(20) /'(MPa)'/,
     &     units(21) /'(no dim)'/,
     &     units(22) /'(no dim)'/,
     &     units(23) /'(no dim)'/,
     &     units(24) /'(no dim)'/,
     &     units(25) /'(no dim)'/,
     &     units(14) /'(m)'/,
     &     units(15) /'(m)'/,
     &     units(16) /'(m)'/

C     BEGIN
      size_head = size(head)
C     nscalar=(iovapor+ioliquid)*iopressure+iosaturation+iotemperature
C     calculation done in avs_io()
      nout = nscalar + iocord
      if(ifdual .ne. 0)then
         dual_char = 'Dual '
         if (iriver2 .eq. 0) then
            tailstring = '_sca_dual_node'
         else
            tailstring = '_wel_dual_node'
         end if
      else
         dual_char = ''
         if (iriver2 .eq. 0) then
            tailstring = '_sca_node'
         else
            tailstring = '_wel_node'
         end if
      endif

      call namefile2(icall,lu,ioformat,tailstring,idz)

c     Header is only written to the first tecplot file
      if (icall .gt. 1 .and. altc(1:3) .eq. 'tec') return
      
      title( 1) = trim(dual_char) // 'Liquid Pressure (MPa)'
      title( 2) = trim(dual_char) // 'Vapor Pressure (MPa)'
      title( 3) = trim(dual_char) // 'Temperature (deg C)'
      title( 4) = trim(dual_char) // 'Saturation'
      title( 5) = trim(dual_char) // 'Hydraulic Head (m)'
      title( 6) = trim(dual_char) // 'Water Vapor Pressure (MPa)'
      title( 7) = trim(dual_char) // 'Porosity'
      title( 8) = trim(dual_char) // 'Source (kg/s)'
      title( 9) = trim(dual_char) // 'Liquid Density (kg/m**3)'
      title(10) = trim(dual_char) // 'Vapor Density (kg/m**3)'
      title(11) = trim(dual_char) // 'X Permeability (log m**2)'
      title(12) = trim(dual_char) // 'Y Permeability (log m**2)'
      title(13) = trim(dual_char) // 'Z Permeability (log m**2)'
      title(17) = trim(dual_char) // 'Zone'
      if (vol_flux .and. net_flux) then
         title(18) = trim(dual_char) // 'Net Liquid Flux (kg/s/m**3)'
         units(18) = '(kg/s/m**3)'
      else if (vol_flux) then
         title(18) = trim(dual_char) // 'Liquid Flux (kg/s/m**3)'
         units(18) = '(kg/s/m**3)'
      else if (net_flux) then
         title(18) = trim(dual_char) // 'Net Liquid Flux (kg/s)'
      else
         title(18) = trim(dual_char) // 'Liquid Flux (kg/s)'
      end if
      title(19) = trim(dual_char) // 'Vapor Flux (kg/s)'
      title(20) = trim(dual_char) // 'Capillary Pressure (MPa)'
      title(21) = trim(dual_char) // 'Water Fraction'
      title(22) = trim(dual_char) //'Super-Critical/Liquid CO2 Fraction'
      title(23) = trim(dual_char) // 'Gaseous CO2 Fraction'
      title(24) = trim(dual_char) // 'Dissolved CO2 Mass Fraction'
      title(25) = trim(dual_char) // 'CO2 Phase State'
      title(14) = 'X coordinate (m)'
      title(15) = 'Y coordinate (m)'
      title(16) = 'Z coordinate (m)'
      
      if(altc(1:3).eq.'avs' .and. altc(4:4) .ne. 'x') then
C---  Max number of scalars is 9 + 3 coordinates
         string = ''
         write (string, '(i2.2)') nout
         ic1 = 3
         ic2 = ic1 + 3
         do i = 1, nout
            write (string(ic1:ic2), '(2x,i1)') 1
            ic1 = ic2 + 1
            ic2 = ic1 + 3
         end do
         length = len_trim(string)
         write(lu, '(a)') string(1:length)

         string = ''
         if (iocord .ne. 0) then
c     Write X coordinate
            if (icnl .ne. 3 .and. icnl .ne. 6) 
     &           write(lu,200) trim(title(14)), trim(units(14))
c     Write Y coordinate
            if (icnl .ne. 2 .and. icnl .ne. 5) 
     &           write(lu,200) trim(title(15)), trim(units(15))
c     Write Z coordinate
            if (icnl .ne. 1 .and. icnl .ne. 4) 
     &           write(lu,200) trim(title(16)), trim(units(16))
         end if
         if (iozid .eq. 1) then
            write(lu,200) trim(title(17)), trim(units(17))
         end if
         if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,200) trim(title(1)), trim(units(1))
         end if
         if (iopressure .eq. 1 .and. iovapor .eq. 1) then
            write(lu,200) trim(title(2)), trim(units(2))
            if (iadif .eq. 1) then
               write(lu,200) trim(title(6)), trim(units(6))
            end if
         end if
         if (iocapillary .eq. 1) then
            write(lu,200) trim(title(20)), trim(units(20))
         end if         
         if (iotemperature .eq. 1) then
            write(lu,200) trim(title(3)), trim(units(3))
         end if
         if (iosaturation .eq. 1) then
            write(lu,200) trim(title(4)), trim(units(4))
         end if
         if (ioco2 .eq. 1) then
            write(lu,200) trim(title(21)), trim(units(21)) 
            write(lu,200) trim(title(22)), trim(units(22)) 
            write(lu,200) trim(title(23)), trim(units(23)) 
            write(lu,200) trim(title(24)), trim(units(24)) 
            write(lu,200) trim(title(25)), trim(units(25)) 
         end if
         if (iohead .eq. 1 .and. size_head .ne. 1) then
            write(lu,200) trim(title(5)), trim(units(5))
         end if
         if (ioporosity .eq. 1) then
            write(lu,200) trim(title(7)), trim(units(7))
         end if
         if (iosource .eq. 1 ) then
            write(lu,200) trim(title(8)), trim(units(8))
         end if
         if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,200) trim(title(9)), trim(units(9))
         end if
         if (iodensity .eq. 1 .and. iovapor .eq. 1) then
            write(lu,200) title(10)
         end if
         if (iopermeability .eq. 1) then
            write(lu,200) trim(title(11)), trim(units(11))
            write(lu,200) trim(title(12)), trim(units(12))
            write(lu,200) trim(title(13)), trim(units(13))
         end if
         if (ioflx .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,200) trim(title(18)), trim(units(18))
         end if
         if (ioflx .eq. 1 .and. iovapor .eq. 1) then
            write(lu,200) trim(title(19)), trim(units(19))
         end if         
         
      else
C--   altc is 'avsx', 'tec', 'sur'
         tstring = ''
         if(altc(1:4) .eq. 'avsx') then
            write (formstring, 100) ' : '
            write(tstring2, '("nodes at ", g16.9, " days ")') days
         else if (altc(1:3) .eq. 'tec') then
            write (formstring, 300) 
            tstring2 = 'VARIABLES = '
         else if (altc(1:3) .eq. 'sur') then
            write (formstring, 100) ', '
            tstring2 = 'node'
         end if

         ic1 = 1
         ic2 = len_trim(tstring2)
         if (iocord .ne. 0) then
c     Write X coordinate
            if (icnl .ne. 3 .and. icnl .ne. 6) then
               write(tstring,formstring) trim(title(14))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
c     Write Y coordinate
            if (icnl .ne. 2 .and. icnl .ne. 5) then
               write(tstring,formstring) trim(title(15))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
c     Write Z coordinate
            if (icnl .ne. 1 .and. icnl .ne. 4) then
               write(tstring,formstring) trim(title(16))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
         end if
         if (altc(1:3) .eq. 'tec') then
            write(tstring,formstring) 'node'
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if  
         if (iozid .eq. 1) then
            write(tstring,formstring) trim(title(17))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
            write(tstring,formstring) trim(title(1))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iopressure .eq. 1 .and. iovapor .eq. 1) then
            write(tstring,formstring) trim(title(2))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            if (iadif .eq. 1) then
               write(tstring,formstring) trim(title(6))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
         end if
         if (iocapillary.eq. 1) then
            write(tstring,formstring) trim(title(20))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iotemperature .eq. 1) then
            write(tstring,formstring) trim(title(3))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iosaturation .eq. 1) then
            write(tstring,formstring) trim(title(4))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (ioco2 .eq. 1) then
            write(tstring,formstring) trim(title(21))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(22))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(23))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(24))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(25))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iohead .eq. 1 .and. size_head .ne. 1) then
            write(tstring,formstring) trim(title(5))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (ioporosity .eq. 1) then
            write(tstring,formstring) trim(title(7))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iosource .eq. 1 ) then
            write(tstring,formstring) trim(title(8))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
            write(tstring,formstring) trim(title(9))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iodensity .eq. 1 .and. iovapor .eq. 1) then
            write(tstring,formstring) trim(title(10))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (iopermeability .eq. 1) then
            write(tstring,formstring) trim(title(11))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(12))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(13))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (ioflx .eq. 1 .and. ioliquid .eq. 1) then
            write(tstring,formstring) trim(title(18))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (ioflx .eq. 1 .and. iovapor .eq. 1) then
            write(tstring,formstring) trim(title(19))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
         end if
         if (altc(1:3) .eq. 'tec') then
            write(lu, 400) verno, jdate, jtime, trim(wdd)
         end if
         write(lu, '(a)') tstring2(ic1:ic2)         
      end if

 100  format("('", a, "', a)")
 200  format(a, ', ', a)
 300  format("('",' "', "', a, '",'"',"')")
 400  format('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
 500  format('VARIABLES = "node" ', a)
      return
      end

