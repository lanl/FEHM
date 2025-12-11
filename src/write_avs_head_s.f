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
     &     verno, wdd, i_vtk
      use comdi, only : head
      use comsi, only : iPlastic,flag_excess_shear, flag_principal
      use davidi
      implicit none

      integer maxscalar
      parameter (maxscalar = 37)
      integer neq,nscalar,lu,ifdual
      integer i,j,iolp,iovp,nout,icall,idz,iriver2,iocord_temp
      integer size_head, size_pcp, istart, iend, ic1, ic2, length
      character*5 dual_char
      character*14 tailstring
      character*20 tailstring2
      character*20 formstring
      character*40 title(maxscalar+4)
      character*20 units(maxscalar+4)
      character*1200 tstring2
      character*400 string
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
     &     units(14) /'(m)'/,
     &     units(15) /'(m)'/,
     &     units(16) /'(m)'/,     
     &     units(17) /'(no dim)'/,
     &     units(18) /'(kg/s)'/,
     &     units(19) /'(kg/s)'/,
     &     units(20) /'(MPa)'/,
     &     units(21) /'(no dim)'/,
     &     units(22) /'(no dim)'/,
     &     units(23) /'(no dim)'/,
     &     units(24) /'(no dim)'/,
     &     units(25) /'(no dim)'/,
     &     units(26) /'(m)'/,
     &     units(27) /'(m)'/,    
     &     units(28) /'(m)'/, 
     &     units(29) /'(MPa)'/,          
     &     units(30) /'(MPa)'/,
     &     units(31) /'(MPa)'/,
     &     units(32) /'(no dim)'/,     
     &     units(33) /'(MPa)'/,           
     &     units(34) /'(MPa)'/,
     &     units(35) /'(MPa)'/,
     &     units(36) /'(no dim)'/,     
     &     units(37) /'(MPa)'/,
     &     units(38) /'(MPa)'/,
     &     units(39) /'(deg)'/,
     &     units(40) /'(kg/m**3)'/,
     &     units(41) /'(kg/m**3)'/

C     BEGIN
      size_head = size(head)
C     nscalar=(iovapor+ioliquid)*iopressure+iosaturation+iotemperature
C     calculation done in avs_io()
      iocord_temp = iocord

      if(ifdual .ne. 0)then
          if (iriver2 .eq. 0) then
            if (iodual .eq. 1) then
               dual_char = 'Dual '
               tailstring = '_sca_dual_node'
            else if (iogdkm .eq. 1) then
               dual_char = 'GDKM '
               tailstring = '_sca_gdkm_node'
c gaz 051821 debug               
c               if (icnl .eq. 0) then
c                  iocord = 3
c               else
c                  iocord = 2
c               end if
            end if               
         else
            dual_char = 'Dual '
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
      nout = nscalar + iocord

      call namefile2(icall,lu,ioformat,tailstring,idz)

c     Header is only written to the first tecplot file
      if (icall .gt. 1 .and. altc(1:3) .eq. 'tec' .and. (iogrid .eq. 0 
     &   .or. (ifdual .eq. 1 .and. iogdkm .eq. 1))) return
      
      title( 1) = trim(dual_char) // 'Liquid Pressure (MPa)'
      title( 2) = trim(dual_char) // 'Vapor Pressure (MPa)'
      if (altc(1:3) .eq. 'tec'.and. i_vtk.ne.0) then
c gaz 122722 added heading change for vtk
         title( 3) = trim(dual_char) // 'Temperature (deg C)'
      else if (altc(1:3) .eq. 'tec') then
         title( 3) = trim(dual_char) // 'Temperature (<sup>o</sup>C)'
      else
         title( 3) = trim(dual_char) // 'Temperature (deg C)'
      end if
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
      title(14) = 'X coordinate (m)'
      title(15) = 'Y coordinate (m)'
      title(16) = 'Z coordinate (m)'      
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
      title(21) = trim(dual_char) // 'Water Saturation'
      title(22) = trim(dual_char) //
     &     'Super-Critical/Liquid CO2 Saturation'
      title(23) = trim(dual_char) // 'Gaseous CO2 Saturation'
      title(24) = trim(dual_char) // 'Dissolved CO2 Mass Fraction'
      title(40) = trim(dual_char) // 'CO2 Liquid Density (kg/m**3)'
      title(41) = trim(dual_char) // 'CO2 Gas Density (kg/m**3)'
      title(25) = trim(dual_char) // 'CO2 Phase State'
      title(26) = 'X displacement (m)'
      title(27) = 'Y displacement (m)'
      title(28) = 'Z displacement (m)'  
c gaz 052317      
      if(flag_principal.eq.0) then
       title(29) = 'X stress (MPa)'
       title(30) = 'Y stress (MPa)'
       title(31) = 'Z stress (MPa)'      
       title(33) = 'XY stress (MPa)'
       title(34) = 'XZ stress (MPa)'
       title(35) = 'YZ stress (MPa)'  
      else if(flag_principal.eq.1) then
       title(29) = 'Sigma Max (MPa)'
       title(30) = 'Sigma 2   (MPa)'
       title(31) = 'Sigma Min (MPa)'      
       title(33) = 'Ang Sigma Max-Z'
       title(34) = 'Ang Sigma Max-X'
       title(35) = 'Ang Sigma Max-Y'           
      endif
      title(32) = 'Volume Strain'  
      title(36) = 'Plastic strain (no dim)'  
      title(37) = 'Youngs Mod (MPa)'  
      title(38) = 'Excess Shear (MPa)'  
      title(39) = 'Shear Angle (deg)'  

      
      if(altc(1:3).eq.'avs' .and. altc(4:4) .ne. 'x') then
C---  Max number of scalars is number of parameters + 3 coordinates
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
c gaz 051621 modfying use of units() to avoid duplication for avs
c will leave variables with no units as the unit array is relevant             
c     Write X coordinate
            if (icnl .ne. 3 .and. icnl .ne. 6) 
     &           write(lu,201) trim(title(14))
c     Write Y coordinate
            if (icnl .ne. 2 .and. icnl .ne. 5) 
     &           write(lu,201) trim(title(15))
c     Write Z coordinate
            if (icnl .ne. 1 .and. icnl .ne. 4) 
     &           write(lu,201) trim(title(16))
         end if
         if (iozid .eq. 1) then
            write(lu,200) trim(title(17)), trim(units(17))
         end if
         if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,201) trim(title(1))
         end if
         if (iopressure .eq. 1 .and. iovapor .eq. 1) then
            write(lu,201) trim(title(2))
            if (iadif .eq. 1) then
               write(lu,201) trim(title(6))
            end if
         end if
         if (iocapillary .eq. 1) then
            write(lu,201) trim(title(20))
         end if         
         if (iotemperature .eq. 1) then
            write(lu,201) trim(title(3))
         end if
         if (iosaturation .eq. 1) then
            write(lu,200) trim(title(4)), trim(units(4))
         end if
         if (ioco2 .eq. 1) then
            write(lu,200) trim(title(21)), trim(units(21)) 
            write(lu,200) trim(title(22)), trim(units(22)) 
            write(lu,200) trim(title(23)), trim(units(23)) 
            write(lu,200) trim(title(24)), trim(units(24)) 
            write(lu,200) trim(title(40)), trim(units(40)) 
            write(lu,200) trim(title(41)), trim(units(41)) 
            write(lu,200) trim(title(25)), trim(units(25)) 
         end if
         if (iohead .eq. 1 .and. size_head .ne. 1) then
            write(lu,201) trim(title(5))
         end if
         if (ioporosity .eq. 1) then
            write(lu,200) trim(title(7)), trim(units(7))
         end if
         if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,201) trim(title(9))
         end if
         if (iodensity .eq. 1 .and. iovapor .eq. 1) then
            write(lu,201) title(10)
         end if
         if (iopermeability .eq. 1) then
            write(lu,201) trim(title(11))
            write(lu,201) trim(title(12))
            write(lu,201) trim(title(13))
         end if
         if (iosource .eq. 1 ) then
            write(lu,201) trim(title(8))
         end if
         if (ioflx .eq. 1 .and. ioliquid .eq. 1) then
            write(lu,201) trim(title(18))
         end if
         if (ioflx .eq. 1 .and. iovapor .eq. 1) then
            write(lu,201) trim(title(19))
            
         end if         
         if (iodisp .eq. 1) then
            write(lu,201) trim(title(26))
            write(lu,201) trim(title(27))
            if (icnl .eq. 0) 
     &           write(lu,201) trim(title(28))
         end if 
         if (iostress .ne. 0) then
            write(lu,201) trim(title(29))
            write(lu,201) trim(title(30))
            if (icnl .eq. 0) then
             write(lu,201) trim(title(31))
            end if
            write(lu,201) trim(title(33))
            if (icnl .eq. 0) then
               write(lu,201) trim(title(34))
               write(lu,201) trim(title(35))
            end if
            if(iPlastic.eq.1) then
               write(lu,201) trim(title(36))
            endif
            if (flag_excess_shear.eq.1) then
               write(lu,201) trim(title(36))
               write(lu,201) trim(title(37))
               write(lu,201) trim(title(38))
            end if
         end if  	    
         if (iostrain .eq. 1) then
            write(lu,200) trim(title(32)), trim(units(32))
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
            write(tstring,formstring) trim(title(40))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(41))
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
         if (iosource .eq. 1 ) then
            write(tstring,formstring) trim(title(8))
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
         if (iodisp .eq. 1) then
            write(tstring,formstring) trim(title(26))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(27))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            if (icnl .eq. 0) then
               write(tstring,formstring) trim(title(28))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
         end if 
         if (iostress .ne. 0) then
            write(tstring,formstring) trim(title(29))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            write(tstring,formstring) trim(title(30))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            if (icnl .eq. 0) then
               write(tstring,formstring) trim(title(31))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
            write(tstring,formstring) trim(title(33))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
            if (icnl .eq. 0) then
               write(tstring,formstring) trim(title(34))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
               write(tstring,formstring) trim(title(35))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
            if(iPlastic.eq.1) then
               write(tstring,formstring) trim(title(36))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            endif
            if (flag_excess_shear.eq.1) then
               write(tstring,formstring) trim(title(37))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
               write(tstring,formstring) trim(title(38))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
               write(tstring,formstring) trim(title(39))
               tstring2 = tstring2(ic1:ic2) // tstring
               length = len_trim(tstring)
               ic2 = ic2 + length
            end if
         end if 
         if (iostrain .eq. 1) then
            write(tstring,formstring) trim(title(32))
            tstring2 = tstring2(ic1:ic2) // tstring
            length = len_trim(tstring)
            ic2 = ic2 + length
	 endif         
         if (altc(1:3) .eq. 'tec') then
            write(lu, 400) verno, jdate, jtime, trim(wdd)
            if (iogrid .eq. 1 .and. iocord .eq. 0) then
               write (lu, 600)
            end if
         end if
         write(lu, '(a)') tstring2(ic1:ic2)         
      end if
      iocord = iocord_temp

 100  format("('", a, "', a)")
 200  format(a, ', ', a)
 201  format(a)    
 300  format("('",' "', "', a, '",'"',"')")
 400  format('TITLE = "', a30, 1x, a11, 1x, a8, 1x, a, '"')
 500  format('VARIABLES = "node" ', a)
 600  format('FILETYPE = "SOLUTION"')

      return
      end

