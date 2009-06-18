      subroutine startup(tajj, tasii)
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
CD1 Perform miscellaneous startup calculations.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 07-DEC-93    Z. Dash        22      Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/startup.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:00   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:17:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:14:44   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:27:00   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:10:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:58 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.33   Fri May 31 15:28:24 1996   gaz
CD2 decrease allocation for b array in 6dof dpdp
CD2 
CD2    Rev 1.32   Fri May 17 13:00:24 1996   hend
CD2 Added optional parameter to time macro for initial time
CD2 
CD2    Rev 1.31   Fri Apr 26 16:14:36 1996   gaz
CD2 changes for mdnodes
CD2 
CD2    Rev 1.30   Wed Apr 03 15:15:52 1996   hend
CD2 Removed New Allocation Unless Tracer Problem
CD2 
CD2    Rev 1.29   Thu Mar 21 13:19:24 1996   hend
CD2 Fixed allocation of sehvariables for trac
CD2 
CD2    Rev 1.28   Mon Mar 04 16:15:46 1996   hend
CD2 Removed uneccesary calculations from coneq1 and added trac input 
CD2 option
CD2 
CD2    Rev 1.27   Fri Feb 16 11:36:12 1996   zvd
CD2 Modified requirements.
CD2 
CD2    Rev 1.26   Fri Feb 02 12:01:30 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.25   Fri Jan 12 17:58:50 1996   llt
CD2 changed mmgetblk arguments
CD2 
CD2    Rev 1.24   Thu Jan 11 12:51:58 1996   gaz
CD2 fixed requirements for reduced degree of freedom
CD2 
CD2    Rev 1.23   Tue Jan 09 14:14:06 1996   llt
CD2 gaz changes
CD2 
CD2    Rev 1.22   12/13/95 10:31:18   robinson
CD2 Changed allocation of a matrix to handle multi-species problems
CD2 
CD2    Rev 1.21   12/11/95 15:28:50   gaz
CD2 more rdof changes
CD2 
CD2    Rev 1.20   12/11/95 14:17:48   gaz
CD2 changed mdof ,idof definition to accomdate reduced degree of 
CD2 freedom dpdp
CD2 
CD2    Rev 1.19   11/15/95 15:33:18   gaz
CD2 changes to handle irdof with 4dof
CD2 
CD2    Rev 1.18   09/18/95 08:49:22   gaz
CD2 corrected (again) estimate for nnop
CD2 
CD2    Rev 1.17   09/15/95 09:23:56   gaz
CD2 corrected error in nnop estimate
CD2 
CD2    Rev 1.16   09/15/95 09:08:52   gaz
CD2 new estimates for nnop and nbd
CD2 
CD2    Rev 1.15   08/16/95 11:00:08   zvd
CD2 Corrected write to iatty when unassigned.
CD2 
CD2    Rev 1.14   08/08/95 08:57:26   awolf
CD2 rarng taken out by SEH and moved to infiles.f
CD2 
CD2    Rev 1.13   08/03/95 17:08:10   gaz
CD2 added calls to md_nodes for multiply defined nodes
CD2 
CD2    Rev 1.12   08/03/95 08:23:02   robinson
CD2 Fixed problem with the writing of the AVS geo file
CD2 
CD2    Rev 1.11   08/02/95 16:47:50   llt
CD2 changed allocation of ib array from nnop to nbd
CD2 
CD2    Rev 1.10   08/02/95 13:00:02   awolf
CD2 Allocate a_axy here instead of in allocmem now
CD2 
CD2    Rev 1.9   06/23/95 14:44:06   gaz
CD2 made nbd and nemx larger
CD2 
CD2    Rev 1.8   06/01/95 16:45:50   gaz
CD2 made change to allow idof=6 for h-m-a
CD2 
CD2    Rev 1.7   03/29/95 12:35:44   llt
CD2 changed allocation of bp array
CD2 
CD2    Rev 1.6   03/20/95 13:34:26   gaz
CD2 moved call to rarng before call to anonp
CD2 
CD2    Rev 1.2   03/28/94 16:41:44   robinson
CD2 Removed unneeded array.
CD2 
CD2    Rev 1.1   03/18/94 15:55:50   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.0   01/20/94 10:28:14   pvcs
CD2 original version in process of being certified
CD2 
c 17-mar-93
c because we are only allocating space for b at the 
c call to the solvers,we must allocate space for array
c ib here
c also the realloc routines are commented out
c we need to address this
c 12/16/94 gaz alloc space for nop just before slvesu
c 12/22/94 gaz read directly in nar for fill-in level
c 12/23/94 gaz got rid of ldn1,ldn2
c          cal resize nop after slvesu
c 12/23/94 gaz printout resized bmatrix size
c 1/3/95   gaz called thikness(1) for 2-d problems l 410 
c 1/3/95   gaz commented out x3>x2 l 310
C***********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3
CD3   Identifier      Type     Use  Description
CD3
CD3   tajj            REAL*8   I    Elapsed cpu time (for reading input and
CD3                                   coefficient generation)
CD3   tasii           REAL*8   O    Cpu time at start of solution computations 
CD3
CD3 Interface Tables
CD3
CD3   None
CD3
CD3 Files
CD3
CD3   Name                     Use  Description
CD3
CD3   iout                     O    File used for general code output
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
CD4   cord            REAL*8   fbs    Contains the coordinates of each node
CD4   cpr             REAL*8   fdd    Rock specific heat at each node
CD4   denei           REAL*8   fcc    Energy accumulation term
CD4   deneh           REAL*8   fdd    Last time step energy accumulation term
CD4                                     at each node
CD4   denej           REAL*8   fdd    Last time step energy accumulation time
CD4                                     derivative at each node
CD4   denh            REAL*8   fdd    Last time step mass accumulation term at
CD4                                     each node
CD4   deni            REAL*8   fcc    Mass accumulation term
CD4   denj            REAL*8   fdd    Last time step mass accumulation time
CD4                                     derivative at each node
CD4   denr            REAL*8   fdd    Rock density at each node
CD4   dit             REAL*8   fdd1   Array containing time step changes
CD4   dstm            REAL*8   fcc    Steam mass
CD4   icnl            INT      faai   Problem dimension
CD4   ig              INT      fhh    Variable order lu decomposition
CD4                                     information
CD4   iieos           INT      fddi1  Thermodynamics set at each node
CD4   iirb            INT      fhh    Inverse of irb
CD4   iout            INT      faai   Unit number for output file
CD4   ipa             POINTER  fee    Pointer to variable array a
CD4   ipistrw         POINTER  fbb    Pointer to variable array istrw
CD4   ipnelm          POINTER  fbb    Pointer to variable array nelm
CD4   ipnop           POINTER  fbb    Pointer to variable array nop
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   irb             INT      fhh    Array containing the reordered node
CD4                                     numbers
CD4   irdof           INT      david1 Reduced degree of freedom model used
CD4   istrw           INT      fbb    Starting positions in sx(nr,9) array of
CD4                                     finite element coefficients for each
CD4   mdof            INT      comrxni The maximum number of degrees
CD4                                    of freedom necessary to solve 
CD4                                    the heat & mass and the 
CD4                                    tracer solutions    
CD4 
CD4   nar             INT      fbb    Array containing gauss elimination order
CD4                                     for each node
CD4   neigh           INT      faai   Maximum number of neighbors occur in
CD4                                     tracer solution
CD4   nelm            INT      fbb    Initially information about nodes in each
CD4                                     element, later nodal connectivity
CD4   nemx            INT      faai   Number of unique (geometrically) elements
CD4   neq             INT      faai   Number of nodes, not including dual
CD4                                     porosity nodes
CD4   nop             INT      fbb    Matrix sparsity structure for lu
CD4                                     decomposition
CD4   nopt            INT      fhh    Array indicating active variables
CD4   npvt            INT      fhh    Pivot information for the lu
CD4                                     decomposition matrix
CD4   pflow           REAL*8   fdd    Flowing pressure at each source node
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   phini           REAL*8   fdd2   Initial pressure at each node
CD4   pho             REAL*8   fdd    Last time step pressure at each node
CD4   pnx             REAL*8   fdd    Permeability in the x-direction, liquid
CD4                                     velocity in the x-direction, vapor
CD4                                     velocity in the x-direction
CD4   pny             REAL*8   fdd    Permeability in the y-direction liquid
CD4                                     velocity in the y direction, vapor
CD4                                     velocity in the y-direction
CD4   pnz             REAL*8   fdd    Permeability in the z-direction, liquid
CD4                                     velocity in the z-direction, vapor
CD4                                     velocity in the z-direction
CD4   ps              REAL*8   fdd    Porosity at each node
CD4   psini           REAL*8   fdd2   Initial porosity at each node
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   so              REAL*8   fdd    Last time step saturation at each node
CD4   sx1             REAL*8   fbc    Contains volume associated with each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4   thx             REAL*8   fdd    Thermal conductivity x-direction
CD4   thy             REAL*8   fdd    Thermal conductivity y-direction
CD4   thz             REAL*8   fdd    Thermal conductivity z-direction
CD4   tini            REAL*8   fdd2   Initial temperature at each node
CD4   to              REAL*8   fdd    Last time step temperature at each node
CD4   volume          REAL*8   fdd    Volume associated at each node
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   airctr                   Manage the isothermal air-water calculations
CD4   anonp                    ?
CD4   bcon                     Adjust/manage boundary conditions

CD4   co2ctr                   Control isothermal air-water simulation
CD4   coeffc                   Change the coefficients of the equation of
CD4                              state fits
CD4   contr                    Write contour plot data
CD4   disk                     Read/write initial/final state data
CD4   dpdp                     Control dual porosity/dual permeability solution
CD4   dual                     Find dual porosity contributions to nodes
CD4   peint                    Initializes pressures and temperatures
CD4   plot                     Write history plot data
CD4   porosi                   Calculate pressure dependant porosity
CD4   radius                   Figure out radius in radial problem
CD4   rarng                    Rearrange coefficients in 2-d problem
CD4   setord                   Set up the order of solution for the equations
CD4                              at each node
CD4   sice                     Control ice simulation
CD4   slvesu                   Perform symbolic factorization and calculate
CD4                              computer storage necessary for the solvers
CD4   split                    Split rectangles(bricks) into 4 triangles(12
CD4                              tetrahedrals)
CD4   steady                   Calculate steady state solution
CD4   storsx                   Store and retrieve element coefficients
CD4   thickness                Adjust thicness in 2-d problems
CD4   tyming          REAL*8   Calculate elapsed cpu time
CD4   varchk                   Determine variable set and make n-r corrections
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
CD5   caz             REAL*8   Dummy argument to function tyming
CD5   i               INT      Loop index
CD5   ib              INT      Temporary array for slvesu
CD5   idum1           INT      Temporary array for slvesu
CD5   idum2           INT      Temporary array for slvesu
CD5   icode            INT      variable used in memory management
CD5   imm             INT      Counter
CD5   nbytes          INT      Space to be allocated for variable arrays
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
CD8 This routine initializes various startup conditions which are 
CD8 necessary for the code.
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.7   Provide Restart Capability
CD9 2.7.3 Resume the calculation
CD9 2.8   Provide Multiple Realization Option
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
CPS BEGIN startup
CPS   ???
CPS END startup
CPS
C***********************************************************************

      use comflow
      use comcouple

      use davidi
      use comgi
      use comei
      use comfi
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comii
      use comwt
      use compart, only : ptrak
      use comsptr, only : sptrak
      use comzone
      use commeth
      use comco2
      use comsi, only : idof_stress
      implicit none

      integer icall, dummyint, j, iconv_tmp
      integer mdof,ldnmax,ldnmaxtmp,ncont,neqp1,igmax
      integer nstorederiv,nstoresolve,i,icode,irdum
      integer ncon_size,ncon_size_orig
      integer open_file
      integer idb_file, if_debug
      real*8 tajj, tasii, tyming, dummyreal
      real*4 caz(2)
      real*8 rnwn1, rnwn2, rnwn3, rnwn, rnwd1, rnwd2, rnwd3, rnwd
      real*8 xmin_bbox, ymin_bbox, zmin_bbox
      real*8 xmax_bbox, ymax_bbox, zmax_bbox, scale_factor
      integer, allocatable :: ib(:)
      integer, allocatable :: idum1(:)
      integer, allocatable :: idum2(:)
      integer, allocatable :: ncon_temp(:)
      integer, allocatable :: nop_temp(:)      
      integer i1,i2,ipiv,icnt
      integer count, num_zone
      integer nsbb,idofdum
      real*8 very_small, fac_nop, fac_mult
      parameter (very_small= 1.d-30)

c     Calculate rho1grav
c      rho1grav = crl(1,1)*(9.81d-6)
c      rho1grav = 997.*9.81d-6
c**** Deallocate printout information array for variable porosity model
c**** if allocated and not a variable porosity problem
      if (iporos .ne. -4) then
         if (allocated(nskw3)) deallocate(nskw3)
      endif

c**** calculate number of neighbors and connected elements ****
      allocate(idum1(n0),idum2(n0))
      if (icnl .eq. 0)  then
c**** three dimensional problem ****
         nemx = 200
         neigh = 200
      else
c**** two   dimensional problem ****
         nemx = 200
         neigh = 200
      end if

      neigh = max0(neigh, nemx)
      
c**** set z coordinate = y coordinate for 2-d problems ****
c**** also insure gravity direction is correct ****
      if (icnl .ne. 0)  then
         if (igrav .eq. 3)  igrav = 2
      end if

c**** set idof = 3 if stress solution enabled and icnl = 0 (3-d) ****
      if (icnl .eq. 0 .and. idof_stress .lt. 0)  idof = 3

      if (n0 .lt. neq)  then
         write(ierr, 6000)  n0, neq
         if (iout .ne. 0) write(iout, 6000)  n0, neq
         if (iptty .gt. 0)  write(iptty, 6000)  n0, neq
 6000    format(/, 1x, '**** n0(', i6, ') .lt. neq(', i7, ') **** ', 
     *        'check parameter statements ****')
         stop
      end if
c**** set idof = 3 if air-water-heat
      if (ico2 .gt. 0) idof = max(idof,3)
      if (idof_co2 .gt. 0) idof = max(idof,idof_co2)
      if (icnl.eq.0) i1 = 3
      if (icnl.ne.0) i1 = 2
      if (idof_stress.gt.3) idof = max(idof + i1,5)
      if (istrs.ne.0) idof = max(idof,i1)
c**** 

      if(idpdp.ne.0.and.idof_stress.ne.0) then
        if(iptty.ne.0) then
         write(iptty,*) '>>>>> dpdp not implemented with stress <<<<'
         write(iptty,*) 'stopping in startup'
        endif
        if(iout.ne.0) then
         write(iout,*) '>>>>> dpdp not implemented with stress <<<<'
         write(iout,*) 'stopping in startup'
        endif
        stop
      endif
c********* Make dpdp, two-phase work properly
c check for dpdp(ngas or 2dof)
      if(idpdp .ne. 0) then    
c reorder id irdof gt 0
         if(irdof.gt.0.and.islord.eq.0) islord=1
         if(ico2 .gt. 0) then
            idof=6
         else
            idof=4
         endif
      endif

c calculate new idof based on irdof
      mdof=idof
      if(idof.eq.2) then
         if(irdof.gt.0) mdof=1
      else if(idof.eq.3) then
         if(irdof.eq.1) mdof=1
         if(irdof.eq.2) mdof=2
      else if(idof.eq.4) then
         if(irdof.gt.0) mdof=2
      else if(idof.eq.6) then
         if(irdof.eq.1) mdof=2
         if(irdof.eq.2) mdof=4
      endif

c Find the maximum degrees of freedom necessary by taking the 
c max of the degrees of freedom needed for the heat and mass solution
c and the max degrees of freedom necessary for the tracer solution
      if (compute_flow) then
         mdof = max(mdof,mdof_sol)
      else
         if (iccen .eq. 1)  then
            mdof = mdof_sol
         else
            mdof = 1
         end if
      end if

c**** set eos cofficient set eq to 1 since fit was good ****
c****          the entire range only done if iieos(i) = 0 ****
      if (iwelb .ne. 1)  then
         do i = 1, n
            if (iieos(i) .eq. 0)  iieos(i) = 1
         end do
      end if

c**** compatability scaling (viscosity is in pa-sec, not in mpa-sec) ***
      do i = 1, n
         if (idoff .gt. 0) then
            pnx(i) = pnx(i) * 1.0d+06
            pny(i) = pny(i) * 1.0d+06
            pnz(i) = pnz(i) * 1.0d+06
         end if
         if(ico2.ge.0.or.ice.ne.0) then
           thx(i) = thx(i) * 1.0d-06
           thy(i) = thy(i) * 1.0d-06
           thz(i) = thz(i) * 1.0d-06
         endif
         if (cpr(i) .gt. 1.0)  cpr(i) = cpr(i) * 1.0d-06
      end do

c**** set permeabilities = 0 for porosity = 0.0 ****
      do i = 1, n
         if (abs(ps(i)) .le. zero_t .and. idoff .gt. 0)  then
            pnx(i) = zero_t
            pny(i) = zero_t
            pnz(i) = zero_t
         end if
      end do

      days = 0.0
      mink = neq

c**** complete element information ****
!      if(contim.ge.0) then
         call contr (0)
!      else
!         call contr_days (0)
!      endif

c**** complete pest information(if interpolation is used) ****
      call pest(-1)
      if (mlz .ne. 0. and. irun .eq. 1 .and. lda .le. 0) call split(0)
      if (lda .le. 0.and.irun.eq.1) then
         if(ivf.ne.-1) then
            if(ianpe.ne.2) call anonp
            if(ianpe.ne.0) call anonp_ani
         else
            call structured(1)
            call structured(3)
         endif
      endif
c mlz is now the size of sx(mlz, j)
c      if (mlz .ne. 0 .and. irun .eq. 1 .and. lda .le. 0) call split(0)
      if (mlz .ne. 0 .and. irun .eq. 1 .and. lda .le. 0) call split(1)
      if (lda .ne. 0.and.irun.eq.1) call storsx
c gaz 11-09-2001 
      if(lda.le.0.and.imdnode.ne.0.and.irun.eq.1) then
         ncon_size=nelm(neq_primary+1)
         ncon_size_orig = ncon_size
         nr=nr+1
         allocate(ncon_temp(ncon_size))
         do i=1,ncon_size
            ncon_temp(i)=nelm(i)
         enddo
         deallocate(nelm)
         call md_nodes(4,0,ncon_size)
c gaz 1-24-2003
         nelmd = ncon_size
         allocate(nelm(ncon_size))
         do i=1,ncon_size_orig
            nelm(i)=ncon_temp(i)
         enddo
         deallocate(ncon_temp)
         call md_nodes(5,0,ncon_size)
      endif
c gaz 11-09-2001 allocation of istrw_itfc and istrw_cold
c done here after call to anonp,storsx, or structured
c still could ne changed in add_gdpm
      ncon_size=nelm(neq_primary+1)
      if(idpdp.eq.0) then
         if (.not. allocated (istrw_itfc)) 
     &        allocate(istrw_itfc(ncon_size))
         if (.not. allocated (istrw_cold))
     &        allocate(istrw_cold(ncon_size))
      else
         if (.not. allocated (istrw_itfc)) 
     &        allocate(istrw_itfc(2*ncon_size))
         if (.not. allocated (istrw_cold))
     &        allocate(istrw_cold(2*ncon_size))
      end if
      istrw_itfc = 0
      istrw_cold = 0
c gaz 11-11-2001 moved this before call to add_gdpm
      if (icnl .ne. 0)  then
         call radius
         call thickness(1)
      end if
c
c     Call routine to adjust the connectivity array if needed
c     to add gdpm nodes
c     RJP 12/13/06 modified to add river_ctr for wellbore
c trying to maintain compatibility of add_gdpm and river_ctr
      if (irun.eq.1) then
c add gdpm connections
         if (gdpm_flag .ne. 0) call add_gdpm

c add river or well connections
         if (nriver .ne. 0) then
c river_ctr(1) call is made in incoord
c        call river_ctr(1)
            call river_ctr(2)
            call river_ctr(3)	
c printout to file volumes,areas etc.
            call river_ctr(-3)	
            call river_ctr(4)
c printout to file connectivities
            call river_ctr(-4)
         endif
      end if

      if(interface_flag.ne.0) call setconnarray
c call sx_combine to break connections to fixed type BCs
      if(ianpe.eq.0) then
         if (irun.eq.1.and.inobr.eq.0) call sx_combine(1)
      else
         if (irun.eq.1.and.inobr.eq.0) call sx_combine_ani(1)
      endif
c call fluxo now to calculate neighbors if necessary
      call flxo(1)

c allocate space for velocity and alpha connection terms (seh)
      if (iccen .ne. 0) then
         sehsize=0
         do i=1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            ipiv=nelmdg(i)
            do icnt=ipiv+1,i2
               sehsize=sehsize+1
            enddo
         enddo
         if (idpdp.ne.0) sehsize=2*sehsize
         if(irun.eq.1) then
            allocate(alphaconl(sehsize),alphaconv(sehsize))
            allocate(sehvell(sehsize),sehvelv(sehsize))
         end if
         alphaconl=0
         alphaconv=0
         sehvell=0
         sehvelv=0
         sehdonevel=0
      endif

      iad = 0
c     Moved the head calls to after the initialization flag w/ iflg=-1
c**** print coordinates and elements to graphics files ****
c     Moved this call until after rarng is called - BAR 8-3-95
c
c     Moved earlier so that a_axy can be read in from disk
c
c       PHS  3/24/2004  adding a_wvxy to store water vapor mass flux
c
      neqp1=neq+1
      ncont=nelm(neqp1)
      ldna=ncont-neqp1
      if(irun.eq.1) then

         if(idof_co2.gt.0) then
            allocate(c_axy(ldna))
            allocate(c_vxy(ldna))
            c_axy = 0.d0
            c_vxy = 0.d0
         endif

         if (idoff .gt. 0) then
            if (idpdp .ne. 0) then
               allocate(a_axy(2*ldna+neq))
            else
               allocate(a_axy(ldna))
            end if
            if (irdof .ne. 13) then
               if (idpdp .ne. 0) then
                  allocate(a_vxy(2*ldna+neq))
                  if (iadif .eq. 1) then
                     allocate(a_wvxy(2*ldna+neq))
                  else
                     allocate(a_wvxy(1))
                  end if
               else
                  allocate(a_vxy(ldna))
                  if (iadif .eq. 1) then
                     allocate(a_wvxy(ldna))
                  else
                     allocate(a_wvxy(1))
                  end if
               end if
            else
               allocate(a_vxy(1), a_wvxy(1))
            end if
         else
            allocate(a_axy(1),a_vxy(1),a_wvxy(1))
         end if
      end if
      a_axy=0
      a_vxy=0
      a_wvxy=0

c
      if (iread .gt. 0) call diskread
      close (1010)

c added new feature -- set initial time if requested in time macro
      if (irsttime.ne.0) days=rsttime
      if (iread .le. 0 .and. (tin0 .gt. 0.0 .or. tin1.gt.0.0))
     &     call peint
c      if (iread.le.0 .and. igrad. ne.0) call gradctr(1)
c     gaz 9-26-04
      if (igrad. ne.0) call gradctr(1)
      if (iconv. ne.0) call convctr(1)
c
c check that source terms are ok on restart for wellbore model
c      call welbor(-1)

c      set clathrate properties equal to water properties for a test
c
      call icectr(-2,0)

      call icectrco2(-2,0)

c gaz 11-11-2001 moved before call to add_gdpm
c     if (icnl .ne. 0)  then
c        call radius
c        call thickness(1)
c     end if

      dtot = day * 86400.0

      if (sx1(1) .gt. 0.0)  then
         upwgt = 1.0 - upwgt
         upwgta = 1.0 - upwgta
      end if
      dnwgt = 1.0 - upwgt
      dnwgta = 1.0 - upwgta

c Compute space for gmres

c check if bcgs solution requested
      if(accm.eq.'bcgs') then
         kgmres = (3 + 2 * ( abs(north)+1 ) )*mdof*neq
         if (iout .ne. 0) write(iout,*) 'bcgstab method chosen '
         if (iatty.ne.0) write(iatty,*) 'bcgstab method chosen '
      elseif(accm.eq.'gmre') then
         kgmres = (north + 1) * mdof * neq
         if (iout .ne. 0) write(iout,*) 'gmres method chosen '
         if (iatty.ne.0) write(iatty,*) 'gmres method chosen '
      endif
c adjust space for a  matrix when it is called
      ldn=ldna*idof**2

c**** set up equation solvers ****
c  estimate parameters for LU arrays
      igmax=0
      do i=1,neq
         igmax=max(igmax,nar(i))
      enddo
c      nnop = neq*(ldna/neq+1)*igmax**2.1 +neqp1
c new estimate GAZ 080707
       if(igmax.le.1) then
	  nnop = ncont
	 else
	  fac_nop = 1.5
	  fac_mult = 2.2
        nnop = fac_nop*ldna*fac_mult**(igmax-1) + neqp1
	 endif
      if(ncont.gt.nnop) nnop=ncont
      nbd  = nnop*mdof**2
      if(irdof.gt.0) then
         irdum=mdof
      else
         irdum=0
      endif
      ldnmax = max(ldn,ldna*mdof_sol**2)
      ldnmaxtmp = ldnmax

c     ib needs only nnop storage (dof independent)
      if(irun.eq.1) then
         allocate(ib(nnop),nop(nnop))
c
c check for input inconsistencies
c
         call diagnostics(0)
c
c Compute approximate bounding box of each node
c before modification/simplification of connectivity.
c
         allocate(dx_bbox(neq))
         allocate(dy_bbox(neq))
         allocate(dz_bbox(neq))
C         if((icnl .eq. 4).or.(icnl .eq. 5).or.(icnl .eq. 6))then
         allocate(da_bbox(neq))
C         endif
         do i = 1, neq
            call node_midedge_box_size
     &           (dx_bbox(i),dy_bbox(i),dz_bbox(i),
     &           xmin_bbox, ymin_bbox, zmin_bbox, 
     &           xmax_bbox, ymax_bbox, zmax_bbox, i)
c
c        Check if the prblem is radial r,z instead of xyz
c
            if((icnl .eq. 4) .or. (icnl .eq. 5))then
               da_bbox(i) = (xmax_bbox**2 - xmin_bbox**2)*acos(-1.0)
            elseif(icnl .eq. 6)then
               da_bbox(i) = (ymax_bbox**2 - ymin_bbox**2)*acos(-1.0)
            else
               da_bbox(i) = dx_bbox(i)*dy_bbox(i)
            endif
         enddo
C
         if_debug = 0
         if(if_debug .ne. 0)then
            idb_file = open_file('tmp_debug_dx_dy_dz.tbl','unknown')
            write(idb_file,*)'dx_bbox, real'
            write(idb_file,*)'dy_bbox, real'
            write(idb_file,*)'dz_bbox, real'
            write(idb_file,*)'da_bbox, real'
            write(idb_file,*)'dv_bbox, real'
            write(idb_file,*)'sx1    , real'
            write(idb_file,*)'vrat   , real'
            do i = 1, neq
               write(idb_file,'(i8,9e15.6)')
     1              i,dx_bbox(i),dy_bbox(i),dz_bbox(i),
     1              da_bbox(i), dx_bbox(i)*dy_bbox(i)*dz_bbox(i),
     2              sx1(i),(dx_bbox(i)*dy_bbox(i)*dz_bbox(i))/sx1(i)
            enddo
            close(idb_file)
         endif
         if_debug = 0
c
c s kelkar may 20 09 moved call to ptrac1 here to do the geometry
c before simplify_ncon for -ve porosities. But moved the call to 
c load_omr_flux_array from ptrac1 to fehmn.f where ptrac1 used to be called
         if(sptrak) then
            call ptrac1
         endif

c simplify connectivity based on porosity=0.0
c only for isotropic problems (ianpe=0)
c
         if(ianpe.eq.0) then
            call   simplify_ncon
     &           (0,nelm,nelmdg,nop,istrw,ib,neq,idof,ka,ps,i)
            if (connect_out) call connections_list
            if (igauss .le. 1) then
               deallocate (nop)
            else
               nop = 0
            end if
            if(i.gt.0) then
               if (iout .ne. 0) write(iout, 100) i
               if (ischk .ne. 0) write(ischk, 100) i
               if (iptty .ne. 0) write(iptty, 100) i
 100           format ('>>>> ', i8, ' nodes eliminated (porosity <0)')
               do i=1,neq
                  if(nelm(i)+1.eq.nelm(i+1)) then
                     ka(i) = 0.0
                     sk(i) = 0.0
                     qh(i) = 0.0
                     qc(i) = 0.0
                     if (idoff .gt. 0) then
                        pnx(i)=very_small
                        pny(i)=very_small
                        pnz(i)=very_small
                     end if
                     pho(i)=0.1
                     to(i)=crl(6,1)
                     if (irdof .ne. 13 .or. ifree .ne. 0) s(i)=1.0
                  endif
               enddo
            endif
         end if
      end if
c
c call reordering algorithm here
c
      call renum(1)
c
      if(irun.eq.1) then
         if (igauss .gt. 1) then
            if (gdpm_flag.ne.0) then
               call slvesu (neq_primary, ib, nelm_primary, nop,
     *              nelmdg_primary, npvt, nar, irb, iirb, 
     *              idum1, idum2, ldnmaxtmp, nbd, kgmres, nelmd, nnop, 
     *              mdof, abs(north), ireord, iout, iatty, 
     *              irdum, nbnd,accm)
            else
               call slvesu (neq, ib, nelm, nop, nelmdg, npvt, nar, irb,
     *              iirb, idum1, idum2, ldnmaxtmp, nbd, kgmres,  
     *              nelmd, nnop, mdof, abs(north), ireord, 
     *              iout, iatty, irdum, nbnd,accm)
            endif
            allocate (nop_temp(nnop))
            nop_temp(1:nnop) = nop(1:nnop)
            deallocate (nop)
            allocate (nop(nnop))
            nop = nop_temp
            deallocate (nop_temp)
         else
! Skip call to slvesu
            do j=1,neq
               iirb(j) =  j
               irb (j) =  j
            end do

            nbnd   =  0
            do j=2,neq
               nbnd =  max0( nelm(j+1)-nelm(j-1),nbnd )
            end do
            if(gdpm_flag.ne.0) then
               do j=1,neq_primary
                  npvt(j) = nelmdg_primary(j)
               enddo
            else
               do j=1,neq
	          npvt(j) = nelmdg(j)
               enddo
            endif          
         endif
         deallocate(ib)
      end if
c     reset some array sizes
c GAZ 10/9/98 adjust space for 1dof problem
c GAZ 02/6/01 adjust space for gdpm problem
c
      if(irdof.eq.13) then
         if (iout .ne. 0) write(iout, 6009) 
         if (iptty .ne. 0) write(iptty, 6009) 
 6009    format('SZ only simulation, true size of a matrix :')
      endif
      if(gdpm_flag.ne.0) then
         if (iout .ne. 0) write(iout, 6012) 
         if (iptty .ne. 0) write(iptty, 6012) 
 6012    format('>>>> gdpm invoked, true size of matrices :')
      endif
      if(gdpm_flag.ne.0) then
         ldnmax = max(ldn,ldna*mdof_sol**2)
         if (iout .ne. 0) write(iout, 6011) ldnmax
         if (iptty .ne. 0) write(iptty, 6011) ldnmax
      endif
      if(irdof.eq.13) then
         ldn=ldna
         ldnmax = max(ldn,ldna*mdof_sol**2)
         if (iout .ne. 0) write(iout, 6011) ldnmax
         if (iptty .ne. 0) write(iptty, 6011) ldnmax
      endif
 6011 format(1x,'storage available for a matrix resized to ',
     &     i10, '<<<<<<') 
c  gaz 032208 change to richards eq
      if(irun.eq.1) then
         if(jswitch.ne.0) then
            call airctr(-2,0)
         else if (compute_flow .or. iccen .eq. 1) then
            allocate(a(ldnmax))
         endif
      end if
c modify storage of LU factorization matrix for 6-2
c 3-1 full GMRES schemes
      if(compute_flow .or. iccen .eq. 1) then
         if (igauss .gt. 1) then
            nsbb = nop(neq_primary+1)-(neq_primary+1)
         else
            nsbb = nelm(neq_primary+1)-(neq_primary+1)
         end if
         if (jswitch.eq.1 .and. idpdp .eq. 0) then
            nbd = nsbb
         else if (jswitch.eq.1 .and. idpdp .ne. 0) then
            nbd = 4*nsbb
         else if(irdof.ge.0) then
            nbd = nsbb*mdof**2
         else if(irdof.eq.-3) then
            nbd = nsbb*1**2
         else if(irdof.eq.-6) then
            nbd = nsbb*2**2
         else
            nbd = (nelm(neqp1)-neqp1)*mdof**2
         end if
      end if
      if(.not.compute_flow .and. allocated(nop)) then
         if (iccen .eq. 1 .and. igauss .gt. 1) then
! Do nothing, we need the nop array
         else
            deallocate(nop)
         end if
      end if
      if(irun.eq.1) then
         if (compute_flow .or. iccen .eq. 1) then
            if (ice .eq.0) then
               allocate(bp(max((mdof*n0),(idof*n0),2*n0)))
               bp = 0.0
            else
               allocate(bp(max((mdof*n0),((idof+1)*n0),2*n0)))
               bp = 0.0
            end if
c 
c  make sure size for stress solution is at least 6*neq 
c  for CO2 app (hopefully less for sequentially coupled)
c
c             if (istrs.ne.0) then
c             if(icnl.ne.0.and.idof_co2.ge.2) idofdum = 5
c             if(icnl.eq.0.and.idof_co2.ge.2) idofdum = 6
c             if(icnl.eq.0.and.ico2.eq.0) idofdum = 5
c             if(icnl.eq.0.and.ico2.ge.1) idofdum = 6            
c             if(idof_stress.le.3)then
c              idofdum = max(idof_stress,idof_co2)
c             endif
c             if (max((mdof*n0),2*n0).lt.idofdum*n0) then
c	        deallocate(bp)
c	        allocate(bp(n0*idofdum))
c             endif
c             idof = idofdum
c	      endif
c
c            bp = 0
c
         
         end if
      end if
c     New allocation of bigblock
      if(compute_flow .or. iccen .eq. 1) then
         if(irdof.ne.13) then
            nstorederiv = 40*n0 + n7a
         else
c there are 15 arrays we can eliminate for saturated only
            nstorederiv = 25*n0 + 15*1 + n7a
         endif
         nstorepiv = n0*mdof**2
         nstoresolve = nbd + kgmres + nstorepiv
         icall = -1
      else
c    try only minimal memory
c         nstorederiv = 39*n0 + 3*n7 + n7a
         nstorederiv = 39*1 + 3*1 + 1
         nstorepiv = 1
         nstoresolve = 1
         icall = -2
      end if
      nbigblock = max(nstorederiv,nstoresolve)
      call storage_derivatives(icall,irun)
c     write out new size for b matrix 
      if(compute_flow .or. iccen .eq. 1) then
         if (iout .ne. 0) write(iout, 6010) nbd
         if (iatty .ne. 0) write(iatty, 6010) nbd
      else
         if (iout .ne. 0) write(iout, 6010) 1
         if (iatty .ne. 0) write(iatty, 6010) 
      end if
 6010 format(1x,'storage available for b matrix resized to ',
     &     i10, '<<<<<<') 
      
c     Initialize eos values as necessary (for air)
c     first set phase state for wtsi
c      call wtsictr(-1)           
      if (ico2.lt.0.and.ice.eq.0) then
         call airctr(-1,0)
!     else if (ico2.lt.0.and.ice.ne.0) then
!        call icectr(-1,0)
      endif
c     Calculate rho1grav
      rho1grav = crl(1,1)*(9.81d-6)
c      rho1grav = 997.*9.81d-6

c     Moved calls to airctr for head option to below the iflg = -1 call
c if head input has been used,convert to pressures
c
      if(ihead.ne.0.or.ichead.ne.0) then
         if(head0.gt.0) then
            phi_inc = head0*crl(1,1)*(-grav)
         end if

c convert from head to pressure if no initial value file
         if(iread.le.0 .or. .not. pres_read) then
c
c this is new(set IC heads equal to lowest fixed head
c
c            call wtsictr(12)
c
            call airctr(7,0)
         else
c now add pressure inc. equal to the head inc.
c in reading pressures from intial value field
            do i= 1,n
               pho(i) = pho(i) + phi_inc
            enddo
         endif
c convert from head to pressure for boundary nodes           
         call airctr(9,0)
      endif

      if (ichead.ne.0 .and. head0 .gt. 0.) then
         phi_inc = head0*rol0*(-grav)
      end if
      
c**** set array ordering if necessary ****
      if(compute_flow .or. iccen .eq. 1) then
         call setord
      end if
      
c**** initialize coeffients adjust volumes in dual porosity calcs ****
      call dual (2)

c**** initialize coefficients adjust volumes in dpdp calcs ****
      call dpdp (1)

      do i = 1, n
         volume(i) = sx1(i)
         if (to (i) .le. zero_t)   to (i) = tin0
         if(irdof.ne.13) then
            if (pho(i) .le. zero_t)   pho (i) = pein
         end if
         phi (i) = pho(i)
         t (i) = to(i)
         vtot = vtot + volume(i)
      end do

      if (mass_read) then
         iconv_tmp = iconv
         iconv = 1
         call convctr(3)
         deallocate (mass_var)
         iconv = iconv_tmp
      end if
c  initialize pressures and temperatures
      if(ipini.eq.0.or.iread.eq.0) then      
         phini = pho
         psini = ps
         tini = to           
      else
c phini and tini have be read in from disk, just set pressure
         psini = ps
      endif
c initialize porosity model -5 if necessary
      call porosi(5)
      call stressctr(2,0)
c this call handled in call to porosi above
c
c  find permeability nodes and startup operations(allocate memory)
c
      call stressctr(16,0)
c
c   combine multiple generalzed head BCs
c
      call inflogh(1)
c      
c
c if boundary node fix initial value to fixed value
c
      if(ico2.lt.0.and.iread.le.0.and.ice.eq.0) then
         do i = 1, n
            if(ka(i).eq.-1.or.ka(i).eq.-2) then
               pho(i) = pflow(i)
               phi (i) = pho(i)
            endif
         end do
      endif
c  Set uz correction for water table
c
	call uz_wtctr(100)
c	
c     These arrays are present only if the flow field is
c     being computed
      if(compute_flow) then
         do i = 1, n
            if (ps (i) .le. zero_t)   deni (i) = 0.0
            if(irdof.ne.13) then
               if (ps (i) .le. zero_t)   dstm (i) = 0.0
            endif
         end do
      end if

c modify volumes and heat capacities to affect boundary conditions
      call bcon(1)

c**** initialize ngas varibles ****
      call co2ctr (6)

c bookeeping for air-water or methane
      if(ico2.lt.0.and.ice.eq.0) then
c get cell lengths for wtsi if necessary
         call airctr(11, 0)
         call airctr(6, 0)
      else if(ico2.lt.0.and.ice.ne.0) then
         call icectr(6, 0)
      endif
      if(icarb.eq.1) then
         call icectrco2(6,0)
      endif

c change porosity and permeability if necessary if Gangi model is used
      call porosi(3)

c**** determine initial coefficients for thermo fits ****
      call coeffc

c**** determine initial variable state ****
c gaz 10-18-2001     call sice (1)
      if(compute_flow .or. iccen .eq. 1) then
         if(ice.eq.0) then
            if(icarb.eq.1) then
               call icectrco2(-1,0)
               call icectrco2(14,0)
               call icectrco2(-34,0)
               call icectrco2(3,0)
               call icectrco2(-3,0)
               call icectrco2(-33,0)
               call icectrco2(-35,0)
            else
               call varchk (0, 0)
            end if
         else
            call icectr (-1,0)
            call icectr (1,0)
c added check for hydrate line
            call icectr (-6,0)
c allocate space check properties, dellocate space
            if(idof_meth.ne.7) then
               call icectr (-34,0)
               call icectr (3,0)
               call icectr (-3,0)
               call icectr (-33,0)
               call icectr (-35,0)
	    else
	      call icectr (-34,0)
c     id mobile methane and water
              call hydrate_equil(2,0)
              call icectr (3,0)
              call icectr (-3,0)
              call icectr (-33,0)
              call  hydrate_equil(3,0)
              call icectr (-35,0)
	      
           endif
        endif
      else
         if(ico2.lt.0.and.ice.eq.0) then
c     air water problem
            do i = 1, n0
               rolf(i)=crl(1,1)*
     2              (1.0+crl(3,1)*(phi(i)-crl(4,1)))
            end do
         else if(ico2.lt.0.and.ice.ne.0) then
c gaz 10-18-2001 following used ase dummy parameters
c rnwn1, rnwn2, rnwn3, rnwn, rnwd1, rnwd2, rnwd3, rnwd
            call methane_properties(2,2,phi(i),t(i),rnwn1, rnwn2,
     &           rolf(i),rnwn3, rnwn, rnwd1, rnwd2)
         else
c     nonisothermal problem
            do i = 1, n0
               rnwn1=crl(1,1)+crl(2,1)*phi(i)+
     2              crl(3,1)*phi(i)*phi(i)+
     3              crl(4,1)*phi(i)*phi(i)*phi(i)
               rnwn2=crl(5,1)*t(i)+crl(6,1)*t(i)*t(i)+
     2              crl(7,1)*t(i)*t(i)*t(i)
               rnwn3=crl(8,1)*t(i)*phi(i)+
     2              crl(10,1)*t(i)*t(i)*phi(i)+
     3              crl(9,1)*t(i)*phi(i)*phi(i)
               rnwn=rnwn1+rnwn2+rnwn3
               rnwd1=crl(11,1)+crl(12,1)*phi(i)+
     2              crl(13,1)*phi(i)*phi(i)+
     3              crl(14,1)*phi(i)*phi(i)*phi(i)
               rnwd2=crl(15,1)*t(i)+crl(16,1)*t(i)*t(i)+
     2              crl(17,1)*t(i)*t(i)*t(i)
               rnwd3=crl(18,1)*t(i)*phi(i)+
     2              crl(20,1)*t(i)*t(i)*phi(i)+
     3              crl(19,1)*t(i)*phi(i)*phi(i)
               rnwd=rnwd1+rnwd2+rnwd3
               rolf(i)=rnwn/rnwd
            end do
         end if
      end if
      call dual (1)
      if(compute_flow) then
         call dpdp ( 2 )
      end if

      if(compute_flow ) then
c**** initialize some variables **** 
         am0 = 0.0
         ame = 0.0
         astmo = 0.0
c
c change volumes so mass and energy comes out correct (GAZ 1-29-09)
c
      call vboun(1,0)
c      
c**** calculate initial mass and energy ****
         do i = 1, n
            if(irdof.ne.13) then
               deneh(i) = denei(i) * dtot
               ame = ame + deneh(i) * volume(i)
               astmo = astmo + dstm (i)
               denei(i) = 0.0
            endif
	    if (ifree .ne. 0) then
               so(i) = rlxyf(i)
            else if (irdof.ne.13) then
               so (i) = s(i)
            end if
            denh (i) = deni (i) * dtot
            am0 = am0 + denh (i) * volume(i)
            deni (i) = 0.0
            nopt (i) = 0
         end do
         if(ice.ne.0) then
c**** calculate initial mass and energy for component system ****
            call icectr(7,0)
         elseif(icarb.eq.1) then
            call icectrco2(7,0)
         endif
         amass = am0
         aener = ame
         asteam = astmo
      end if

c**** initialize noncondensible and ice arrays if applicable ****
      call co2ctr (1)
      dtotdm = dtot
      if (icgts .gt. 0)  then
         do nicg = 1, icgts
            if (dit(nicg) .gt. days)  then
               ditnd = dit(nicg)
               go to 81
            end if
         end do
         nicg = 1
         icgts = 0
         ditnd = 1.0d+30
      else
         ditnd = 1.0d+30
      end if
 81   continue

!      if(contim.ge.0) then
         call contr (-1)
         call contr ( 1)
!      else
!         call contr_days (-1)
!         call contr_days ( 1)
!      endif

c**** initalize concentration ****
      if (iccen .ne. 0)  call concen (-1,0,dummyreal)

      call pest(1)

! Moved history file setup below zone volume calculations
!      if (hist_flag) then
!         call plot_new (0, 0.0d00, 0.0d00, 0.0d00, 0.0d00)
!      else
!         call plot (0, 0.0d00, 0.0d00)
!      end if

      tajj = tyming(caz) - tajj

      if (iout .ne. 0) write(iout, 6020)  tajj
      if (iptty .gt. 0)  write(iptty, 6020)  tajj
 6020 format(//, 1x, 'time for reading input, forming coefficients ', 
     *     g10.3, //)

      tasii = tyming(caz)

c**** initialize pflow if needed ****
      do i = 1, n
         if (abs(pflow(i)) .lt. zero_t.or.
     &        pflow(i).eq.-999.)  pflow(i) = pho(i)
         if (icarb.ne.0) then
            if(abs(pflowco2(i)) .lt. zero_t.or.
     &           pflowco2(i).eq.-999.)  pflowco2(i) = phoco2(i)
         endif
         if(ico2.ge.0.or.ice.ne.0) then
            if (qflux(i) .eq. -999 .and. qflxm(i).gt.0.0)
     &           qflux(i) = to(i)
         endif
      end do
c
c call calculation of areas,weights,etc for use in boundary
c conditions
c
      call area_vol_weightctr(1)
c
c apply areas, weights etc.
c
      call area_vol_weightctr(2)
c
c correct air pressures that are less than pref
c correct appropriate boundary pressures and change to seepage face
c conditions if necessary
c
      if(ihead.ne.0) call airctr(10,0)
c
      daynew = day
c release memory for original coordinates
c     BAR - leave around since multiple simulations is a possibility
c      call rarng(1)

      if(sptrak .and. .not. compute_flow) ihf = 1

! zvd - Compute zone volumes for average output
      if (out_zones) then
         allocate (zone_volume(node_azones))
         zone_volume = 0.
         node_ptr_num => node_head_num
         node_ptr => node_head
         num_zone = 0
         outer: do
            if (.not. associated(node_ptr_num)) exit
            num_zone = num_zone + 1
! Number of nodes in current zone
            j = node_ptr_num%node_number
            count = 0
            if (j .ne. 0) then
               sum_values: do 
                  if (.not. associated(node_ptr)) exit
! Number of current node in the zone
                  i = node_ptr%node_number
                  zone_volume(num_zone) = sx1(i)*pnx(i) + 
     &                 zone_volume(num_zone)
                  node_ptr =>node_ptr%nnp
                  count = count + 1
                  if (count .eq. j) exit
               end do sum_values
            else
               zone_volume(num_zone) = 0
            end if
            node_ptr_num =>node_ptr_num%nnp
         end do outer
      end if

c call gdpm_corr to add connections if necessary
c        
       call gdpm_corr(-1)

! zvd - moved history file setup after zone volume calculations
      if (hist_flag) then
! Initial call for file setup
         call plot_new (0, 0.0d00, 0.0d00, 0.0d00, 0.0d00)
      else
         call plot (0, 0.0d00, 0.0d00)
      end if

! zvd - moved initial call here to accomodate new output options
! for trac macro
      if (.not. ptrak .and. .not. sptrak) call plotc1(0,0)

! zvd - Oct 7, 2004 Need to set sk values after we are sure all 
! modifications have been made to nelmdg
      if (.not. compute_flow) then
         if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &        'liq') then
            do i = 1, neq
               sk(i) = a_axy(nelmdg(i) - neq - 1)               
            end do
         end if
         if (flux_flag(1:3) .eq. 'all' .or. flux_flag(1:3) .eq. 
     &        'vap') then
            do i = 1, neq
               if (ico2 .gt. 0 ) then
                  qc(i) = a_vxy(nelmdg(i) - neq - 1) 
               else if (ico2 .lt. 0) then
                  qh(i) = a_vxy(nelmdg(i) - neq - 1) 
               end if
            end do
         end if
      end if

      deallocate(idum1,idum2)
      
      end
