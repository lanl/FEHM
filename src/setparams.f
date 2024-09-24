        subroutine setparams
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
CD1 Initialize/set parameter values.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY 
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 06-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/setparams.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:52   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:15:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:32   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:34   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:36   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:22 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.12   Fri Feb 09 08:54:46 1996   robinson
CD2 Initialized integer flags ichng and izeolites
CD2 
CD2    Rev 1.11   Wed Feb 07 12:16:32 1996   gaz
CD2 initialized constants before scanin search
CD2 
CD2    Rev 1.10   Thu Feb 01 16:06:56 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.9   Tue Jan 09 14:11:24 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.8   12/15/95 11:07:36   gaz
CD2 fixed ichng value to 1 when rlp macro is not called
CD2 
CD2    Rev 1.6   12/13/95 08:48:10   gaz
CD2 changes to accodate variable rlp number of data
CD2 
CD2    Rev 1.5   09/15/95 09:01:48   gaz
CD2 increased nn to 200 to accomdate trease/gable grids
CD2 
CD2    Rev 1.4   01/28/95 13:55:38   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.3   06/20/94 11:09:04   zvd
CD2  
CD2 
CD2    Rev 1.2   03/23/94 14:42:02   robinson
CD2  
CD2 
CD2    Rev 1.1   03/18/94 15:55:42   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:27:30   pvcs
CD2 original version in process of being certified
CD2 
c 16-mar-94
c created n7a
c 12/9/94 gaz changed definition of nnop
c 12/16/94 gaz set nnop back
c 12/23/94 gaz nr,ldn,and nbd don't need to be set 
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   None
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   incoor                   I    Input coordinate data file
CD3   inpt                     I    Input data file
CD3   
C***********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4   None
CD4
CD4 Global Types
CD4
CD4   None
CD4
CD4 Global Variables
CD4
CD4                            COMMON
CD4   Identifier      Type     Block  Description
CD4
CD4   incoor          INT      faai   Unit number for coordinate input file
CD4   inpt            INT      faai   Unit number for input file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   ldn             INT      param  Maximum array space needed for jacobian
CD4                                     array matrix
CD4   ldn1            INT      param  ldn - ldn2
CD4   ldn2            INT      param  ldn / 2 + 1
CD4   lenintg         INT      param  Converts bits to words for allocating
CD4                                     memory
CD4   lenreal         INT      param  Converts bits to words for allocating
CD4                                     memory
CD4   maxor           INT      param  Length of gmres working arrays
CD4   n0              INT      param  Maximum number of nodes allowed
CD4   n044            INT      param  n0 * 44
CD4   n2              INT      param  2 * n0, storage parameter
CD4   n3              INT      param  3 * n0, storage parameter
CD4   n4              INT      param  Array storage parameter for
CD4                                     noncondensible gas solution
CD4   n5              INT      param  Array storage parameter for dual
CD4                                     porosity
CD4                                     solution
CD4   n6              INT      param  Array storage for ice solution
CD4   n7              INT      param  Array storage for tracer solution
CD4   n7a             INT      param  Array storage for drc variable     
CD4   n8              INT      param  Array storage for variable  porosity
CD4                                     solution
CD4   nbd             INT      param  180 * n0 maximum array space for
CD4                                     incomplete lu decomposition matrix
CD4   ne1             INT      param  Array size of common block /febb/,
CD4                                     nelmd + nnop + 4n0
CD4   ne2             INT      param  Array size of common /feb/,
CD4                                     n0 + 9nr + 6nq
CD4   ne3             INT      param  Array size of common block /fbs/,
CD4                                     3n0 + 336
CD4   ne5             INT      param  Array size of common block /fcc/,
CD4                                     39n0 + 5n7
CD4   ne6             INT      param  Array size of common block /fdd/,
CD4                                     36n0
CD4   ne7             INT      param  Array size of common block /fdd1/,
CD4                                     14n7 + 870
CD4   ne8             INT      param  Array size of common block /ffd2/,
CD4                                     7n8 + 4
CD4   ne9             INT      param  Array size of common block /fddi/,
CD4                                     5n0 + 1
CD4   ne10            INT      param  Array size of common block /fhh/,
CD4                                     14n0
CD4   ne11            INT      param  Array size of common block /co2/,
CD4                                     32n4 + 1
CD4   ne12            INT      param  Array size of common block /fgg/,
CD4                                     9nn + 2n3
CD4   ne13            INT      param  Array size of common block /dualp/,
CD4                                     23n5
CD4   ne14            INT      param  Array size of common block /fice/
CD4                                     2*n6 +1
CD4   ne15            INT      param  Array size of common block /ice/
CD4                                     n6
CD4   nelmd           INT      param  Maximum array space for element
CD4                                     connectivity array and (later) the nodal
CD4                                     connectivity array
CD4   nelucm          INT      param  nbd / 64
CD4   neq             INT      faai   Number of nodes, not including dual
CD4   nn              INT      param  Maximum number of connected elements
CD4   nnop            INT      param  Maximum array space for lu decomposition
CD4                                     connectivity array
CD4   nq              INT      param  Maximum array space for each finite
CD4                                     element coefficient array associated
CD4                                     with the stress solution
CD4   nr              INT      param  Maximum space allowed for each finite
CD4                                     element coefficient array
CD4   nspeci          INT      fdd1i  Number of species for tracer solution
CD4   zero_t          REAL*8   param  Value for real parameter zero
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4   
CD4   scanin                    Scan input file for parameters needed prior to
CD4                               data input        
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
CD5   macro           CHAR     Input variable for reading 4 character keywords
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
CD9 2.8 Provide Multiple Realization Option
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
CPS BEGIN setparams
CPS 
CPS   call scanin to read necessary parameters from input file
CPS   
CPS   IF number of nodes hasn't been read
CPS      WHILE input macro is not 'coor' 
CPS        search coordinate input file for 'coor' macro
CPS        IF macro 'coor' is found
CPS           read value for n0
CPS           rewind coordinate input file
CPS        EXIT IF macro 'coor' was found
CPS        END IF
CPS      END WHILE
CPS   END IF 
CPS   
CPS   if macro 'coor' was not found issue error message and terminate 
CPS    program
CPS   
CPS   IF idualp is set
CPS      set n0 to 3*neq
CPS   ELSE IF idpdp is set
CPS      set n0 to 2*neq
CPS   ELSE
CPS      set n0 to neq
CPS   END IF
CPS   
CPS   IF idualp is set
CPS      set n5 to neq
CPS   ELSE
CPS      set n5 to 1
CPS   ENDIF
CPS   
CPS   IF terminal output is being used
CPS      write n0 value to tty
CPS   END IF
CPS           
CPS   IF nspeci is not 0 
CPS      compute parameter n7
CPS   ELSE
CPS      set parameter to 1
CPS   END IF 
CPS   compute parameter n7a
CPS
CPS   compute/set parameters [n2,n3, nr, ldn, nbd, nelmd, nnop, nn, nq, 
CPS    n4, n5, n6, n8, ne1, ne2, ne3, ne5, ne6, ne7, ne8, ne9, ne10, 
CPS    ne11, ne12, ne13, ne14, ne15, maxor, zero_t, ldn2, ldn1, 
CPS    nelucm, n044, lenreal, lenintg]
CPS 
CPS END setparams
CPS
C***********************************************************************
c gaz050223 setparameters_2b is thesame as 2c
      use comchem
      use comdi
      use comzeoli
      use comdti
      use comai
      use combi
      use comwellphys
      use comxi, only : cform
c RJP12/13/06 ADDED following
      use comriv
      implicit none

      character*4 macro
      integer idumm
      integer nx,ny,nz
      
c set some parameters
      ico2 = 0
      icons = 0
      idpdp = 0
      idualp = 0
      izeolites = 0
      ichng = 0
c gaz 061513
      ich_pebi = 0
      if(irun.eq.1) then
         isox=1
         isoy=2
         isoz=3
      end if
c	gdpm nodes
      ngdpmnodes = 0
      gdpm_flag = 0
      gdkm_flag = 0
      iwellp_flag = 0
      nwellphy = 0
      i_rlp = 0

      call scanin
c gaz 101023 added  rewind incoor
      rewind incoor
      if (neq .eq. 0) then
         macro = '    '   
         do while (macro .ne. 'coor' .and. macro .ne. 'fdm ')
            if (cform(3)(1:3) .eq. 'unf' .or. 
     *           cform(3)(1:3) .eq. 'UNF') then
               read (incoor, end = 4103, err = 4103) macro
               if (macro .eq. 'coor')  then
                  read (incoor) neq
                  rewind (incoor)
                  go to 4105
               end if
            else
               read (incoor, 4102, end = 4103, err = 4103) macro
               if (macro .eq. 'coor' .or. macro .eq. 'fdm ') then
                  if(macro .eq. 'coor') then
                     read (incoor, *) neq
                  else if(macro .eq. 'fdm ') then
                     read (incoor,'(a80)') wdd1
                     if(wdd1(1:4).ne.'modf') then
                        read (incoor, *) nx,ny,nz
                        neq = nx*ny*nz
                     else 
                        read (incoor,'(a80)') wdd1
                        open(99,file=wdd1,form = 'formatted',
     *                       status = 'old')
                        read(99,*) nz,ny,nx
                        neq = nx*ny*nz
                        close(99)
                     endif
                  endif
                  rewind (incoor)
                  go to 4105
               end if
            end if
         end do
 4102    format (a4)
 4103    write(ierr, 4104)
         if (iptty .ne. 0 ) write(iptty, 4104)
 4104    format
     *   ( ' ****    COOR or FDM Required Input    **** ', /,  
     *        ' ****---------------------------****', /,
     *        ' ****       JOB STOPPED         ****', /,
     *        ' ****---------------------------****', /)
         stop
      end if

 4105 continue
c gaz 0226 first discovery of neq      
      if(ngdpmnodes.eq.-999) ngdpmnodes = neq
c      
      if (nei .eq. 0) then
         macro = '    '   
         do while (macro .ne. 'elem' .and. macro .ne. 'fdm ')
            if (cform(3)(1:3) .eq. 'unf' .or. 
     *           cform(3)(1:3) .eq. 'UNF') then
               read (incoor, end = 4203, err = 4203) macro
               if(macro .eq. 'elem') then
                  read (incoor) idumm, nei
                  rewind (incoor)
                  go to 4205
               end if
            else
               read (incoor, 4202, end = 4203, err = 4203) macro
               if (macro .eq. 'elem' .or. macro .eq. 'fdm ') then
                  if(macro .eq. 'elem') then
                     read (incoor, *) idumm, nei
                  else if(macro .eq. 'fdm ') then
                     nei = 0
                  endif
                  rewind (incoor)
                  go to 4205
               end if
            end if
         end do
 4202    format (a4)
 4203    write(ierr, 4204)
         if (iptty .ne. 0 ) write(iptty, 4204)
 4204    format
     *   ( ' ****    ELEM or STRU Required Input    **** ', /,  
     *        ' ****---------------------------****', /,
     *        ' ****       JOB STOPPED         ****', /,
     *        ' ****---------------------------****', /)
         stop
      end if
 4205 continue

c after the grid blocks are counted, set ichng=neq
c if appropriate(ichng is the size of the rlperm parameters
c ie if ichng = -1, set ichng = neq
      if(ichng.eq.-1) ichng=neq 

c at least set nrlp to 1 (might not be set if rlp macro not called)
      if(nrlp.eq.0 ) nrlp=1   
      
      if (idualp .eq. 1) then
         n0 = 3 * neq
         neq_primary = neq
      else if (idpdp .ne. 0) then
         n0 = 2 * neq
         neq_primary = neq
      else
c     gdpm nodes must be added to neq
         if (irun .eq. 1) then
            neq_primary = neq
            n0 = neq_primary + ngdpmnodes
            neq=n0
         end if
      end if

      if (idualp .eq. 1) then
         n5 = neq
      else
         n5 = 1
      end if

c RJP 12/13/06 added following for wellbore
c check for river or well nodes
c 
      if(nriver.ne.0) then

         neq_primary = neq
         neq = neq_primary + npoint_riv

         n0 = neq 
        
      endif
      if(.not.gdkm_new) then
       if (iptty .gt. 0) write(iptty, *) 'n0 = ', n0
       if(ngdpmnodes.ne.0.and.iptty.gt.0)
     &  write(iptty,*) neq_primary,' primary nodes'
       if(ngdpmnodes.ne.0.and.iptty.gt.0)
     &  write(iptty,*) ngdpmnodes,' gdpm nodes'
         
       if (iout .gt. 0) write(iout, *) 'n0 = ', n0
       if(ngdpmnodes.ne.0.and.iout.gt.0)
     &  write(iout,*) neq_primary,' primary nodes'
       if(ngdpmnodes.ne.0.and.iout.gt.0)
     &  write(iout,*) ngdpmnodes,' gdpm nodes'
      else
c gaz 092822 gdkm new calculates gdkm nodes later          
c       if (iptty .gt. 0) write(iptty, *) 'n0 = ', n0
        if (iptty .gt. 0) write(iptty, *)
     &  'total nodes(neq_primary + gdkm nodes calculated later)'
       if(iptty.gt.0)
     &  write(iptty,'(t1,i10,1x,a20)') neq_primary,'primary nodes'
       if(iptty.gt.0)
     &  write(iptty,'(t12,a)') 
     &  'gdkm node count is calculated after',
     &  'reading gdkm macro see "input title : gdkm"'       
c gaz 092822 gdkm new calculates gdkm nodes later          
c       if (iout .gt. 0) write(iptty, *) 'n0 = ', n0
        if (iout .gt. 0) write(iptty, *)
     &  'total nodes(neq_primary + gdkm nodes calculated later)'
       if(iout.gt.0)
     &  write(iout,'(t1,i10,1x,a)') neq_primary,'primary nodes'
       if(iout.gt.0)
     &  write(iout,'(t12,a)') 
     &  'gdkm node count is calculated after',
     &  'reading gdkm macro see "input title : gdkm"'    
      endif
      if (nspeci .ne. 0) then
         n7 = nspeci * n0
      else if(ico2.ge.0) then
         n7 = n0
      else
         n7 = 1
      end if
      
      n7a=max(n7,n0,n0*dimdrc)

c**** calculate parameters ****
       n2 = 2 * n0
       n3 = 3 * n0
       nr = 80000
       ldn = 70 * n0
       nbd = 2 * ldn
       nelmd = 40 * n0
       nnop = nelmd * 3 / 2
       nn = 324
       nq = 1
       n4 = n0
       n6 = n0
       n8 = n0
       ne1 = 4 * n0 + nelmd + nnop
       ne2 = n0 + 6 * nr + 6 * nq
       ne3 = 3 * n0 + 336
       ne5 = 32 * n0 + 5 * n7 + 8 * n0
       ne6 = 36 * n0
       ne7 = 11 * n7 + 3890
       ne8 = 7 * n8 + 4
       ne9 = 5 * n0 + 1
       ne10 = 14 * n0
       ne11 = 32 * n4 + 1
       ne12 = 9 * nn + 2 * n3
       ne13 = 23 * n5
       ne14 = 2 * n6 + 1
       ne15 = n6
c zvd - 28-May-09, set in scanin based on north
c       maxor = 201
       zero_t = 1.e-30
       ldn2 = n0*28+1       
       ldn1 = 0
       nelucm = nei           
       n044 = n0 * 44
 
c*** lenreal = converts bits to words for allocating memory
c = 4 if single precision real
c = 8 if double precision real
      lenreal = 8
 
c*** lenintg = converts bits to words for allocating memory
c = 4 if integer*4
      lenintg = 4
c
c allocate memory for boun macro if enabled
c consider doing the same for rlperms       
c
      call flow_boundary_conditions(0)

      return
      end

