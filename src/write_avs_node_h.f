      subroutine write_avs_node_h(
     2                 phi,t,fracw,frachyd,neq,nscalar,lu,
     3                 ioliquid,
     4                 iovapor,
     5                 iopressure,
     6                 iotemperature,
     7                 iofw,
     7                 iofh,
     8                 ifdual)
!***********************************************************************
!  Copyright, 2004,  The  Regents of the University of California.
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
!D1 Output AVS scalar node information for FEHM (hydrate)
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 FEHM V2.30
!D2
!D2 $Log:   /pvcs.config/fehm90/src/write_avs_node_h.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************

      use comai, only : altc, days
      implicit none

      integer maxscalar
      parameter (maxscalar = 8)
      integer neq,nscalar,lu,ioliquid,iovapor,iopressure,iotemperature
      integer iofw,iofh,ifdual
      integer i,iolp,iovp,istep,nout
      real*8 phi(neq), t(neq), fracw(neq)
      real*8 frachyd(neq)
      character*80 title(maxscalar)

      data title( 1) /'Liquid Pressure (MPa), (MPa)'/,
     1     title( 2) /'Temperature (deg C), (deg C)'/,
     2     title( 3) /'Fraction Water (no dim)'/,
     3     title( 4) /'Fraction Hydrate (no dim)'/,
     5     title( 5) /'Dual Liquid Pressure (MPa), (MPa)'/,
     6     title( 6) /'Dual Temperature (deg C), (deg C)'/,
     7     title( 7) /'Dual Fraction Water (no dim)'/,        
     8     title( 8) /'Dual Fraction Hydrate (no dim)'/

      if(ifdual .ne. 0)then
         istep = 4
      else
         istep = 0
      endif
    

C---Max number of scalars is 4, output all 4 variables
C---Output all of these
      nout=4
      if(nout .eq. 4)then

        if(altc(1:4).NE.'avsx') then
          write(lu,90)nout,1,1,1,1,1
          write(lu,'(a56)')title(1+istep)
          write(lu,'(a56)')title(2+istep)
          write(lu,'(a56)')title(3+istep)
          write(lu,'(a56)')title(4+istep)
          do i = 1,neq
            write(lu,100)i,phi(i),
     &           t(i),
     &           max(fracw(i),1d-20),
     &           max(frachyd(i),1d-20)
          enddo
        else
          write(lu,667) days,
     x     (title(i+istep),i=1,4)
          do i = 1,neq
            write(lu,666)i,phi(i),
     &           t(i),
     &           max(fracw(i),1d-20),
     &           max(frachyd(i),1d-20)
          enddo
        end if
       end if
         

 90   format(i1,2x,5(i5,2x))
 100  format(i10.10,2x,4(e16.9,2x))
 666  format(i10.10,4(' : ',e16.9))
 667  format('nodes at ',e10.4,' days ',5(' : ',a40))

      return
      end

