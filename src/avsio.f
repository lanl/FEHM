      module avsio
!    avsio
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
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Include file containing definitions for the avsio library functions.
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2              T. Cherry      00022    Initial Implementation
!D2
!D2 $Log:   /pvcs.config/fehm90/src/avsio.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:20   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:55:08   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:04:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:21:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:55:54   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:38:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   12/12/94 16:31:54   tam
!D2 added iodisplacement to common iocontrol data
!D2 
!D2    Rev 1.1   06/20/94 11:01:52   zvd
!D2 Added prolog and typed variables.
!D2 
!D2    Rev 1.0   01/20/94 10:22:22   pvcs
!D2 original version in process of being certified
!D2
!***********************************************************************
!D3
!D3 INTERFACES
!D3
!D3 None
!D3
!***********************************************************************
!D4
!D4 GLOBAL OBJECTS
!D4
!D4 Global Constants
!D4
!D4   None
!D4
!D4 Global Types
!D4
!D4   None
!D4
!D4 Global Variables
!D4
!D4                            COMMON
!D4   Identifier      Type     Block     Description
!D4
!D4   ***** COMMON Block iocontrol variables *****
!D4   ioconcentration INT      iocontrol AVS output usage flag
!D4                                        (0 - not set, 1 - set)
!D4   iodual          INT      iocontrol AVS output usage flag
!D4                                        (0 - not set, 1 - set)
!D4   ioformat        INT      iocontrol AVS output usage flag 
!D4                                        (1 - unformatted, 2 - ascii format)
!D4   ioliquid        INT      iocontrol AVS output usage flag 
!D4                                        (0 - not set, 1 - set)
!D4   iomaterial      INT      iocontrol AVS output usage flag
!D4                                        (0 - not set, 1 - set)
!D4   iopressure      INT      iocontrol AVS output usage flag
!D4                                        (0 - not set, 1 - set)
!D4   iohead          INT      iocontrol AVS output usage flag
!D4                                        (0 - not set, 1 - set)
!D4   iosaturation    INT      iocontrol AVS output usage flag 
!D4                                        (0 - not set, 1 - set)
!D4   iotemperature   INT      iocontrol AVS output usage flag 
!D4                                        (0 - not set, 1 - set)
!D4   iovapor         INT      iocontrol AVS output usage flag 
!D4                                        (0 - not set, 1 - set)
!D4   iovelocity      INT      iocontrol AVS output usage flag 
!D4                                        (0 - not set, 1 - set)
!D4   
!***********************************************************************
!D5
!D5 LOCAL IDENTIFIERS
!D5
!D5 None
!D5
!***********************************************************************
!D6
!D6 FUNCTIONAL DESCRIPTION
!D6
!D6 None
!D6
!***********************************************************************
!D7
!D7 ASSUMPTIONS AND LIMITATIONS
!D7
!D7 None
!D7
!***********************************************************************
!D8
!D8 SPECIAL COMMENTS
!D8
!D8 None
!D8
!***********************************************************************
!D9
!D9 REQUIREMENTS TRACEABILITY
!D9
!D9 None
!D9
!***********************************************************************
!DA
!DA REFERENCES
!DA
!DA None
!DA
!***********************************************************************
!PS
!PS PSEUDOCODE
!PS 
!PS None
!PS
!***********************************************************************

      integer iomaterial,ioliquid,iovapor,iodual,iohyd,ioflx
      integer iovelocity,iopressure,iotemperature,iosaturation
      integer ioconcentration,iodisplacement,iohead,ioformat
      integer iofw,iofh,ioporosity,iosource,iodensity,iocord
      integer iopermeability,iogeo,iozone,iowt,iokd,iozid,iogrid
      integer iocapillary,ioco2,iogdkm,iogdkmblank,ioheatflux
      integer iodisp, iostrain, iostress      

      integer iaroot, timec_flag
c gaz 020522 need max calls for vtk conversion   
      integer icall_max
c gaz 0310522 need ioscalar (created this date) for vtk conversion 
      integer ioscalar
      real*8  contour_time
      character*8 time_units
      character*50 timec_string, times_string
      character*120 avs_root, geoname, gridstring
      logical net_flux, vol_flux, dit_flag

      end module avsio
