      subroutine fehmn(method, state, ing, out)
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
CD1 Finite Element Heat and Mass Transfer in porous media.
CD1 
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 06-DEC-93    Z. Dash        22      Add prolog.
CD2 1980         G. Zyvoloski           Initial implementation.
CD2
CD2 $Log:   /pvcs.config/fehm90/src/fehmn.f_a  $
CD2
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:58   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:03:24   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:03:50   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 Update the GoldSim / FEHM interface
!D2 
!D2    Rev 2.2   06 Jun 2001 13:23:44   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:01:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:41:16 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.16   Fri Jun 21 15:33:06 1996   hend
CD2 Fixed Ultimately Lame Mistake in the Last Revision
CD2 That Is Not Even Worth Describing
CD2 
CD2    Rev 1.15   Fri Jun 21 15:28:46 1996   hend
CD2 Fixed possible division by 0
CD2 
CD2    Rev 1.14   Tue May 14 14:32:18 1996   hend
CD2 Updated output
CD2 
CD2    Rev 1.13   Wed May 08 14:09:20 1996   hend
CD2 Rearranged and added output
CD2 
CD2    Rev 1.12   Thu Feb 15 10:20:50 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.11   Tue Jan 30 09:24:28 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.10   08/16/95 16:27:54   robinson
CD2 Changed name of variable to set print out interval
CD2 
CD2    Rev 1.9   06/01/95 17:04:18   gaz
CD2 added another call to diagnostics
CD2 
CD2    Rev 1.8   03/20/95 08:25:22   gaz
CD2 un commented call to velocity . this will allow the velocties to be used
CD2 in the transport solution.
CD2 
CD2    Rev 1.7   03/10/95 10:34:46   llt
CD2 commented out call to velocity - gaz
CD2 
CD2    Rev 1.6   02/22/95 10:13:24   llt
CD2 in cont macro, CONTIM was not recognized when avs was used.
CD2
CD2    Rev 1.5   11/30/94 12:20:18   llt
CD2 Added verno definition to dated_*.f routines, so would also tell which
CD2 platform running on.
CD2 
CD2    Rev 1.4   11/29/94 18:22:04   llt
CD2 Changed length of jdate to 11 characters for ibm
CD2 
CD2    Rev 1.3   08/23/94 15:24:12   llt
CD2 Made a subroutine, so that all array allocation can be done before used, 
CD2 required for the ibm.
CD2 
CD2    Rev 1.2   03/18/94 15:53:26   gaz
CD2 Added solve_new and cleaned up memory management.
CD2 
CD2    Rev 1.1   02/28/94 11:54:04   zvd
CD2 Corrected  problem so contour and restart files are always updated for last
CD2 time step.
CD2 
CD2    Rev 1.0   01/20/94 10:23:46   pvcs
CD2 original version in process of being certified
CD2 
c 3/20/95 gaz un-commented out call to velocity
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
CD4   b               REAL*8          array containing the incomplete lu
CD4                                     decomposition of the jacobian matrix
CD4   day             REAL*8   faar   Current time step size in days
CD4   deneh           REAL*8   fdd    Last time step energy accumulation term
CD4                                     at each node
CD4   denei           REAL*8   fcc    Energy accumulation term
CD4   denej           REAL*8   fdd    Last time step energy accumulation time
CD4                                     derivative at each node
CD4   denh            REAL*8   fdd    Last time step mass accumulation term at
CD4                                     each node
CD4   deni            REAL*8   fcc    Mass accumulation term
CD4   denj            REAL*8   fdd    Last time step mass accumulation time
CD4                                     derivative at each node
CD4   dstm            REAL*8   fcc    Steam mass
CD4   dtot            REAL*8   faar   Current time step size in seconds
CD4   eflow           REAL*8   fdd    Energy flow at each source node
CD4   enlf            REAL*8   fcc    Liquid enthalpy
CD4   envf            REAL*8   fcc    Vapor enthalpy
CD4   esk             REAL*8   fdd    Inlet enthalpy associated with a source
CD4   ex              LOGICAL  faay   Logical for existence check (T/F)
CD4   iout            INT      faai   Unit number for output file
CD4   iptty           INT      faai   Unit number for selected tty output
CD4   jdate           CHAR     faac1  Contains the date (mm/dd/yr)
CD4   jtime           CHAR     faac1  Contains the time (hr:mn:sc)
CD4   l               INT      faai   Current time step number
CD4   n               INT      faai   Total number of nodes
CD4   nmat            INT      david2 Array used in the reduced degree of
CD4                                      freedom method
CD4   nstep           INT      faai   Maximum number of time steps
CD4   phi             REAL*8   fdd    Pressure at each node
CD4   pho             REAL*8   fdd    Last time step pressure at each node
CD4   ps              REAL*8   fdd    Porosity at each node
CD4   qh              REAL*8   fdd    Energy source term at each node
CD4   rolf            REAL*8   fcc    Liquid density
CD4   rovf            REAL*8   fcc    Vapor density
CD4   s               REAL*8   fdd    Liquid saturation at each node
CD4   sk              REAL*8   fdd    Source strength of each node
CD4   so              REAL*8   fdd    Last time step saturation at each node
CD4   t               REAL*8   fdd    Temperature at each node
CD4   to              REAL*8   fdd    Last time step temperature at each node
CD4   verno           CHAR     faac   Contains version number of FEHMN code
CD4   volume          REAL*8   fdd    Volume associated at each node
CD4
CD4 Global Subprograms
CD4
CD4   Identifier      Type     Description
CD4
CD4   allocmem                 Allocate memory to dynamic variable arrays
CD4   bnswer                   Invoke thermodynamics and solution routines
CD4   bcon                     Adjust/manage boundary conditions
CD4   co2ctr                   Control isothermal air-water simulation
CD4   concen                   Control trace simulation
CD4   contr                    Write contour plot data
CD4   data                     Zero all arrays and load thermo coefficients
CD4   datchk                   Analyze input data and some generated quantities
CD4   dated                    Find current date and time
CD4   disk                     Read/write initial/final state data
CD4   dual                     Find dual porosity contributions to nodes
CD4   enthp           REAL*8   Calculate enthalpy at a node as a function of
CD4                              temperature and pressure
CD4   fimpf                    Calculate fraction of variables over a given
CD4                              tolerance
CD4   infiles                  Control reading of input data files
CD4   iofile                   Manage the opening of input and output files
CD4   plot                     Write history plot data
CD4   resetv                   Reset variables to old time step values
CD4   setparams                Initialize/set parameter values
CD4   sice                     Control ice simulation
CD4   startup                  Perform miscellaneous startup calculations
CD4   timcrl                   Procedure to adjust timestep
CD4   flow_boundary_conditions Procedure to adjust timestep and BCs 
CD4   tyming          REAL*8   Calculate the elapsed cpu time
CD4   user                     User programmed special calculations
CD4   varchk                   Determine variable set and make n-r corrections
CD4   veloc                    Calculate fluid velocities
CD4   wellbor                  Wellbore input and simulation
CD4   wrtout                   Control output to files and tty
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
CD5   deneht          REAL*8 
CD5   denht           REAL*8 
CD5   eskd            REAL*8 
CD5   enthp           REAL*8 
CD5   flemax          REAL*8 
CD5   flmax           REAL*8 
CD5   ichk            INT
CD5   im              INT
CD5   ja              INT
CD5   mi              INT
CD5   prav            REAL*8  
CD5   pravg           REAL*8 
CD5   tajj            REAL*8   Elapsed cpu time (for reading input and
CD5                              coefficient generation)
CD5   tas             REAL*8   Total elapsed cpu time
CD5   tassem          REAL*8   Elapsed cpu time for time step
CD5   tasii           REAL*8   Cpu time at start of solution computations 
CD5   teinfl          REAL*8 
CD5   tinfl           REAL*8 
CD5   tmav            REAL*8 
CD5   tmavg           REAL*8 
CD5   tyming          REAL*8
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
CD8 This routine handles the entire functioning of the code by calling
CD8 the necessary routines and computations.
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C***********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.5.1 Implement time-step mechanism
CD9 2.7   Provide Restart Capability
CD9 2.7.3 Resume the calculation
CD9 2.8   Provide Multiple Realization Option
CD9 2.9   Interface with GoldSim
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
CPS BEGIN fehmn
CPS 
CPS   call tyming and dated to setup starting time and date
CPS   call iofile to manage the opening of input and output files
CPS   call setparams to calculate/set parameter values
CPS   call allocmem to allocate memory to dynamic variable arrays
CPS   call data to zero all arrays and load thermodynamic coefficients
CPS   call startup to perform miscellaneous startup calculations
CPS   call datchk to analyze input data and some generated quantities
CPS   
CPS   FOR each time step
CPS   
CPS       call timcrl to adjust timestep
CPS       IF user subroutine should be called
CPS          call user to invoke user subroutine
CPS       END IF
CPS       call welbor to compute wellbore solution if enabled
CPS       
CPS       FOR each node
CPS           calculate inflowing enthalpy
CPS           IF source input is a temperature
CPS              convert to enthalpy
CPS           END IF
CPS       END FOR
CPS       
CPS       IF heat and mass transfer solution is disabled 
CPS       
CPS          set last time step size and compute elapsed CPU time
CPS          
CPS       ELSE if heat and mass transfer solution is enabled 
CPS       
CPS          compute CPU time
CPS          call bnswer to invoke thermodynamic and solution routines
CPS          
CPS          IF any variables are out of bounds
CPS             IF maximum number of iterations was exceeded
CPS                write timestep, iterations, and timestep size to tty 
CPS                 if being used
CPS             END IF
CPS             reset days and call resetv to set variables to old 
CPS              timestep values
CPS             IF dual porosity or double porosity/double permability is 
CPS              enabled
CPS                call resetv to set dp/dpdp variables to old timestep 
CPS                 values
CPS             END IF
CPS             call varchk to determine variable set and make n-r corrections
CPS             call dual to invoke varchk for dual porosity nodes
CPS             set new day parameter and increment global iteration count
CPS             return to beginning of loop to adjust timestep (timcrl)
CPS             
CPS          END IF
CPS          
CPS          set last time step size and compute elapsed CPU time for 
CPS           this timestep and total elapsed time
CPS          IF total elapsed time exceeds maximum allowed runtime
CPS             write message to tty if being used about terminating 
CPS              program and exit time step loop
CPS          END IF
CPS          
CPS          initalize time step parameters and calculate outflows and 
CPS           update coefficients for next time step 
CPS          call varchk to update variable set and make n-r corrections
CPS          call dual to invoke varchk for dual porosity nodes
CPS          call fimpf to calculate fraction of variables over a given 
CPS           tolerance
CPS          call bcon to adjust boundary conditions
CPS          
CPS          FOR each node
CPS              IF porosity is present
CPS                 compute timestep fluid outflow and energy outflow
CPS              ELSE
CPS                 compute timestep energy outflow
CPS              END IF
CPS              IF this is a source/sink
CPS                 compute fluid outflow, energy outflow, average 
CPS                  pressure and temperature, and cumulative fluid 
CPS                  and energy outflow
CPS              END IF
CPS              compute mass and energy accumulation 
CPS              set variable last timestep values
CPS          END FOR
CPS          
CPS          call co2ctr to update co2 arrays
CPS          call sice to update sice arrays
CPS          
CPS       END IF
CPS       
CPS       call veloc to calculate velocities
CPS       call concen to obtain concentration solution
CPS       
CPS       compute intermediate flow total
CPS       IF flow total is greater than 0
CPS          compute average temperature, pressure per flow volume
CPS          compute intermediate energy flow and enrgy flow for timestep
CPS       END IF
CPS       
CPS       set printout flag to 0 (no printout for this time step)
CPS       IF tty output is enabled
CPS          set printout flag to 1 (tty printout only)
CPS       END IF
CPS       IF this isn't the finish time and we haven't exceeded the 
CPS        maximum number of timesteps
CPS          increment print-out interval counter
CPS          IF the counter is greater than or equal to the number of 
CPS           timesteps per printout interval
CPS             set printout flag to 2 and reinitialize interval counter
CPS          END IF
CPS       ELSE
CPS          set printout flag to 2
CPS       END IF
CPS       
CPS       IF printout is enabled
CPS          FOR each node
CPS              set last timestep pressure
CPS          END FOR
CPS          calculate mass and energy balance errors
CPS          call co2ctr to calculate gas mass balance error
CPS          call wrtout to initiate output file writes
CPS       END IF
CPS       
CPS       call plot to write history plot data
CPS       IF contour plots enabled or using intervals or time for plot
CPS          call contr to write contour plot data
CPS          call disk to write a restart file at this time
CPS          IF this is the final timestep
CPS             set flag to indicate contour and restart files already 
CPS              written
CPS          END IF
CPS       END IF
CPS       
CPS       IF the final time has been reached
CPS          exit the time step loop
CPS       END IF
CPS       
CPS   END FOR
CPS   
CPS   write final solution (call disk, contr, plot) and messages to 
CPS    output files and tty if using
CPS   
CPS END fehmn
CPS
C***********************************************************************

      use comai
      use combi
      use comci
      use comco2
      use comcomp
      use comdi
      use comdti
      use comei
      use comevap, only : evaporation_flag
      use comfi, only : qtc, qtotc, pci, pcio
      use comflow, only : a_axy
      use compart
      use comriv
      use comrtd, only : maxmix_flag
      use comrxni
      use comsi
      use comsk, only : save_omr
      use comsplitts
      use comsptr
      use comuserc, only : in
      use comwt
      use comxi
      use davidi
      use comfem, only : edgeNum1, NodeElems, ifem, flag_element_perm
      use comfem, only : fem_strain, conv_strain, conv_pstrain  
      use property_interpolate

      use petsc_initialize_package
      use petsc_finalize_package

c     added combi and comflow to get izonef and a_axy arrays
c     in subroutine computefluxvalues
      implicit none

c     These are PC attributes used as compiler directives. They
c     should be set as follows for the PC-RIP version of fehm.

c     If the stand-alone pc version of the code is being used, then
c     the fehmn attribute line should be omitted and the method
c     variable should be passed by reference.

c     For UNIX versions, these lines are ignored as comments.
C!DEC$ ATTRIBUTES dllexport, c :: fehmn
C!DEC$ ATTRIBUTES value :: method
C!DEC$ ATTRIBUTES reference :: method
C!DEC$ ATTRIBUTES reference :: state
C!DEC$ ATTRIBUTES reference :: ing
C!DEC$ ATTRIBUTES reference :: out

      integer(4) method, state
      real(8) ing(*), out(*)

c     irun is a counter for each realization in a multiple simulation
c     run of fehm. It is initialized to 0 in comai

      character*80 filename
      integer open_file, ifail

      logical used, die
      real*8 tims_save, day_saverip, in3save
      real*8 deneht, denht, eskd, enthp, flemax, flmax, prav, 
     *     pravg, tajj, tas, teinfl, tinfl, tmav, tmavg, tyming
      real(8) :: inflow_thstime = 0., inen_thstime = 0.
      real(8) :: contr_riptot = 1.0d+30
      real(8) :: tasii = 0., tassem = 0.
      real*4 caz(2)
c*** water table rise modification
      real*8 water_table_old
      real*8 prop,dpropt,dpropp,p_energy
c*** water table rise modification
      logical it_is_open, intfile_ex
      integer im, ja, mi
      integer :: ichk = 0, tscounter = 0, flowflag = 0
      integer number_of_outbuffers, jpr
      integer :: n_input_arguments = 0
      integer index_N_large, index_in_species
      integer index_in_flag, index_out_buffer
      integer is_ch
      integer :: is_ch_t = 0
      integer :: out_flag = 0
      integer ntty_save
c
c gaz debug 121415 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
      real*8 rel_hum,qin,qin_ng,qin_h2o,qin_enth,enth_avg,pl_dum

      save flowflag, ichk, tassem, tasii, tscounter,
     &     contr_riptot, tims_save, day_saverip, in3save,
     &     water_table_old, ntty_save

c zvd 09-Sep-2011 change size of in array to be consistent with iofile
c modification for GoldSim 
      if (.not. allocated(in)) allocate(in(4))
      in = ing(1:4)

      inquire(unit=6,opened=it_is_open)
      if(method.eq.2) then
c	When run from Goldsim, the normal fehm screen output gets
c	written to fehmn.log instead by opening unit 6 as a file
c	with this name
         if(.not.it_is_open) then
            open(6,file='fehmn.log')
            iptty = 6
         end if
c     Version number (updated in subroutine dated)
         out(1) = vernum
         return
      elseif(method.eq.99) then
c     Cleanup - close all files before starting a new realization
	
         do i = 1, 99
            inquire(unit=i,opened=it_is_open)
            if(it_is_open) then
               close(i)
            end if
         end do
c     Release all dynamic array memory at the end of a GoldSim simulation

         call releasemem

      elseif(method.eq.3) then
c        write(6,*) 'You are running fehm as a dll from rip.'
c        write(6,*) 'This requires that you input the number'
c        write(6,*) 'of input arguments that you are using'
c        write(6,*) 'in this rip simulation.'
c        write(6,*)
c        write(6,*) 'The input arguments are time (mandatory)'
c        write(6,*) 'plus the number of index parameters you'
c        write(6,*) 'in the rip input.'
c        write(6,*) 'Please input the number of input arguments:'
c        read(*,*) n_input_arguments
c        out(1) = n_input_arguments
c        out(2) = 0
         ripfehm = 1
         do i = 1, 80
            filename(i:i) = ' '
         end do
         filename(1:10)='fehmn.gold'
         iread = open_file(filename,'old')

         read(iread,*) n_input_arguments, out_flag
         out(1) = n_input_arguments
c         out(2) = number_of_outbuffers
         out(2) = 0
         call ingold
         close(iread)
c	System call to run batch file fehmn_real.bat
         if(iptty.ne.0) then
            write(iptty,*) 'Running fehmn_real.bat'
         end if
         call system('fehmn_real')
c         out(1) = 3
c         n_input_arguments = 3
c         out(2) = 12
c         number_of_outbuffers = out(2)
      elseif(method.eq.0) then
c     Initialize
         tscounter = 0
         flowflag = 0
         tassem = 0.0
c*** water table rise modification
         water_table_old = -1.d+10
c*** water table rise modification
c bhl_5/15/08
         ipmbal = 0
c bhl_5/15/08

c     Increase counter for simulation number

         if (ripfehm .ne. 0) then
            irun = irun + 1
         else
            irun = in(1)
         end if   

c**** set version number ****
ccc      verno = 'FEHMN XX-XX-XX      '
         tajj = tyming(caz)
         call dated (jdate, jtime)
c--Add copyright write out
         call write_copyright (6)
         call iofile (ichk)

c**** initialize/set parameter values
         call setparams

c**** allocate memory ****
         if(irun.eq.1) call allocmem

c**** call data initialization routine ****
         call data

c**** call co2_properties_interpolation_lookup_table RJP 04/09/07
         if (icarb .ne. 0) then
            inquire(file=nmfil(29), exist=intfile_ex)
            if(.not.intfile_ex) then
               write(ierr, 6010) trim(nmfil(29))
               write(ierr, 6012)
               if (iout .ne. 0) then
                  write(iout, 6010) trim(nmfil(29))
                  write(iout, 6012)
               end if         
               if (iptty .ne. 0) then
                  write(iptty, 6010) trim(nmfil(29))
                  write(iptty, 6012)
               end if         
               stop
            endif
	
            call read_interpolation_data(ifail,nmfil(29))

         end if

 6010    format('CO2 Properties Interpolation Table File not found: ', 
     &        /, a, /, 'Stopping')
 6012    format('Input correct name in control file using, co2in : ',
     &        'filename')

c**** read and write data ****
         in3save = in(3)
         if(in(3).eq.0) then
            in(3) = irun + .0001
         end if
         
         call infiles(in(3))
         if (nriver .ne. 0) call river_ctr(33)

         in(3) = in3save
c**** modify gravity to reflect vector value ****

c transferred to fehmn.f(GAZ 2/19/97)
         grav = -abs(grav) * 1.0d-06
c 
c**** call time varing boundary conditions ****
c 
         call flow_boundary_conditions(2)             
c 
c  10/22/99
c  moving ,the volumes are not defined yet(need for distributed source)
c       call flow_boundary_conditions(3)             
c 
c *** intialize chemistry if needed **
         if(rxn_flag.eq.1) call initchem
c**** call startup calculations ****
         call startup (tajj, tasii)
c moved flow_boundary_conditions(3) from above(could be dangerous!)
         call flow_boundary_conditions(3)             
c**** call to set up area coefficients for md nodes
         call md_nodes(6,0,0)
c**** call data checking routine ****
         call datchk
c**** initial active base variables if necessarhy
         call active_nodes_ctr(-1) 
c gaz 050809 moved to startup
c calculate initial stress field and displacements
c 
c         call stress_uncoupled(1)
c 
c reset boundary conditions for principal stresses (fraction of lithostatic)
c
c         call stressctr(3,0) 
c               
	 if(ico2.lt.0) then
            if (iout .ne. 0) write(iout,834) ifree1
            if (iptty .ne. 0) write(iptty,834) ifree1
	 endif
 834     format('Initial number of partially filled cells: ', i8,/)

c     Store final simulation time so that for a rip simulation
c     we know when the simulation is over

         tims_save = tims

c     rip avs output flag  - initialize

         contr_riptot = contim_rip

         if(.not.compute_flow .and. iccen .ne. 1 .and. 
     &        .not. sptrak) then
            if(allocated(sx)) deallocate(sx)
            if(allocated(istrw)) deallocate(istrw)
         end if

         if(sptrak) then
! Moved opening of "isptr*" files to insptr
!        open(isptr1, file = nmfil(17), status = cstats(17))
!        open(isptr2, file = nmfil(18), status = cstats(18))
!        open(isptr3, file = nmfil(19), status = cstats(19))
c s kelkar may 20 09 moved call to load_omr_flux_array from ptrac1 here
c s kelkar may 28 09 moved call to init_sptr_params from ptrac1 here
c where ptrac1 used to be called
c zvd - 19Nov2010
c     Moved call to sptr_save here, needs to be called after call
c     to load_omr_flux_array
            call init_sptr_params
            if (.not. compute_flow) then
               if (.not. sptr_exists) then
                  call load_omr_flux_array
                  if (save_omr) call sptr_save (1)
               endif
               if(.not.random_flag) then
c                  if(allocated(sx)) deallocate(sx)
c                  if(allocated(istrw)) deallocate(istrw)
               end if
            end if
c            call ptrac1
         endif

c     If block only entered if the code is being called to 
c     perform the calculation (other options are initialization
c     only, return the version number, etc.)
c  change to 4 in new version of rip
      elseif(method.eq.1) then
         ex = .FALSE.
         if(maxmix_flag.ne.0) then
            call generate_flow
         end if

c  stop simulation after stress calc for certain stress input
c set up time-spaced coupling         
         if(istrs_coupl.eq.-4) then
           timestress0 = days
           timestress = timestress0 + daystress
         endif
        if(istrs.ne.0.and.istrs_coupl.eq.0) go to 170
c Before the time step loop create the partitions for zones

         call paractr(1)
c Check for submodel BCs
         call paractr(5)

c
c restart after steady state has been achieved
c 
 999     continue
         if(ifinsh.eq.2) then
            l = 0
            call flow_boundary_conditions(3) 
            daynew = day
            days =0.0
            qt = 0.0
            qtot = 0.0
            qte = 0.0
            qtote = 0.0
            qtc = 0.0
            qtotc = 0.0
            amass = 0.0
            asteam = 0.0
            aener = 0.0
            toutfl = 0.0
            teoutf = 0.0
            dtot_next = day*86400.
            if(iporos.ne.0) call porosi(4)
         endif      

c     Call evaporation routine if this is an evaporation problem
         if (evaporation_flag) call evaporation(1)
c
c set counter for restarted timesteps to zer0
c
         nrestart_ts = 0


! -----------------------------------------------------------
!                  Initialize / Allocate PETSc memory
! -----------------------------------------------------------

            neq = neq + 0

            call petsc_initialize
            

            
            
! -----------------------------------------------------------
!                  Begin Decomposition
! -----------------------------------------------------------
			
            call UGridCreateUGDM(unstructured_grid,ugdm,ndof,option)
            
            call UGridPartition(unstructured_grid,option,Dual_mat,is_new, &
                      num_cells_local_new)
                      
                      
            call UGridDMCreateVector(unstructured_grid,ugdm,global_vec,GLOBAL,option)
            
			call UGridMapIndices(grid%unstructured_grid,discretization%dm_1dof%ugdm, &
			          grid%nG2L,grid%nL2G,grid%nG2A,option)
			          
			          

                      
            call UGridDecompose(discretization%grid%unstructured_grid, &
                              option)                      
                              
            call UGridMapSideSet(grid%unstructured_grid, &
                            region%sideset%face_vertices, &
                            region%sideset%nfaces,region%name, &
                            option,region%cell_ids,region%faces)
                            region%num_cells = size(region%cell_ids)                              

                            
            call DivideVariableArray(data_array, unstructured_grid,face_vertices,n_ss_faces, &
                           region_name,option,cell_ids,face_ids)                            
! -----------------------------------------------------------


c ************** major time step loop ***************************
         do l = 1, nstep
c
c     counter that keeps accumulating when rip is the time step driver
c     Or in a conventional simulation with heat and mass solution
c           gaz debug 082714
            tscounter = tscounter + 1
c       Don't use input value of initial time step anymore
c         if(abs(tscounter).eq.1.and.in(3).ne.0.) then
c            day_saverip = overf
c         end if

c**** time step control via iteration count ****

 100        continue
            call riptime
c
c     Set current index for flow field catalog number (rip option)
c
            flowflag = int(in(2))
cHari 3/1/07
c*** water table rise modification
            if (ripfehm .ne. 0) water_table_old = in(7)
c*** adjust timestep size

            call timcrl

c
c  manage the stress calls when ihms = istrs_coupl = -4
c
           istresscall = 0
c           
            if(ihms.eq.-4) then
               if(days.ge.timestress) then
                  istrs_coupl = -3
                  timestress0 = timestress
                  timestress = timestress0 + daystress
               else
                  istrs_coupl = ihms
               endif
            endif
           
c     Call evaporation routine if this is an evaporation problem
            if (evaporation_flag) call evaporation(2)

            if (ichk .ne. 0)  call user (ichk)
c
c check for possible source movement
c
            if(move_wt.eq.1) call wtsictr(7)

            dtot = day * 86400.0
c     No longer do this because GoldSim enters with in(1) = 0
c     the first time, so we don't need to have the user input the
c     initial delta time  
c         if(abs(tscounter).eq.1.and.in(1).ne.0.) then
c            dtot = 86400.*day_saverip
c         end if

c**** call welbore simulator ****
c         call welbor (1)

c**** calculate inflowing enthalpy ****
            if(ico2.ge.0.and.ice.eq.0) then
               if(icarb.eq.1) then
                  call icectrco2(-4,0)
               else

                  do mi = 1, n
                     eskd = esk(mi)
c**** if source input is a temp convert to enthalpy ****
                     if (itsat.le.10 .and. eskd .lt. 0.0)  then
c potential energy added to inflow energy in function enthp
                        eskd = enthp(mi, -eskd)
c below lines commented out                        
c                        if(ps(mi).eq.0.0.or.idof.le.1) then
c                          eskd=cpr(mi)*(-eskd)
c                        endif
                     else if(itsat.gt.10.and.eskd.lt.0.0) then
c enthalpy and derivatives  
c itsats gt 10 always convert the temperature                     
                        if(igrav.ne.0) then
                           p_energy = -grav*cord(mi,igrav)
                        else
                           p_energy = 0.0d0
                        endif    
                        eskd = abs(eskd)
                        if(eskd.ne.0.0) then
                           call eos_aux(itsat,eskd,phi(mi),1,1,prop,
     &                          dpropt,dpropp)
                           eskd = prop + p_energy    
                           if(ps(mi).eq.0.0.or.idof.le.1) then
                              eskd=cpr(mi)*eskd
                           endif 
                        endif                
                     end if
                     eflow(mi) = eskd
                  end do
               end if
            else if (ice.ne.0) then
               call icectr(-4,0)
            endif

c**** if heat and mass transfer solution is disabled ****
            if ((ihf .ne. ihs) .or. .not. compute_flow ) then

               dtotdm = dtot
               tassem = tyming(caz) - tassem

c**** if heat and mass transfer solution is enabled ****
            else
               dtotdm = dtot
         
               tassem = tyming(caz)

c**** form equations, calculate corrections. ****
               irestart_ts = 0
               if(ipara.eq.0) then
                  call bnswer
               else
                  call bnswer_part
               endif
            
c**** decrease time step if necessary             ****
c**** mlz.ne.0 means thermo variable is out of bounds ****
         
               if (mlz .le. -1)  then
                  irestart_ts = 1
                  if (mlz .eq. -1)  then
                     if (iout .ne. 0) write(iout, 6020)  l, iad, day
                     if (iptty .gt. 0)  write(iptty, 6020)  l, iad, day
 6020                format(/, 1x, 'timestep = ', i6, ' iterations = ', 
     *                    i4, ' timestep size = ', g15.4)
                  else if(mlz.eq.-2) then
                     if (iout .ne. 0) write(iout, *) 
     &                  'restart -  normalization failed'
                     if (iptty .ne. 0) write(iptty, *) 
     &                  'restart -  normalization failed'
                  end if
                  days = days - day
c	          Check for submodel BCs
                  call resetv (0)
                  if(l.eq.1) call paractr(5)
c
c gaz added code so counter is correct for restarted time steps
c
C zvd commented out, now reset if necessary in timcrl
c                  if(icontr.ne.0) then
c                     nicg=nicg-1
c                     ditnd=dit(nicg)
c                  endif
                  if (idualp .ne. 0)  then
                     call resetv (neq)
                     call resetv (neq + neq)
c	          Check for submodel BCs
c gaz debug 102314 see "if(l.eq.1) call paractr(5)" above
                     if(l.eq.1) call paractr(5)
                  end if
                  if (idpdp .ne. 0)  then
                     call resetv (neq)
c	          Check for submodel BCs
                     if(l.eq.1) call paractr(5)
                  end if
                  call varchk (0, 0)
                  call dual (1)
                  daynew = day
                  iad = abs(maxit) + 1
                  ntty_save = ntty
                  ntty = 2
                  call diagnostics(-1)
                  call diagnostics(1)  
                  ntty =   ntty_save 
c
c count restarted timestep
c
                  nrestart_ts = nrestart_ts + 1
                   go  to  100
               end if
          
c
               if (mlz .ge. 1)  then
                  if (iout .ne. 0) write(iout, 6021) l
                  if (iptty .gt. 0)  write(iptty, 6021) l
 6021             format(/, 1x, 'timestep = ', i6, ' timestep ',
     &                 'restarted because of balance errors',
     &                 ' or variable out of bounds')
                  days = days - day
                  call resetv (0)
                  if(icontr.ne.0) then
                     nicg=nicg-1
                     ditnd=dit(nicg)
                  endif
                  if (idualp .ne. 0)  then
                     call resetv (neq)
                     call resetv (neq + neq)
                  end if
                  if (idpdp .ne. 0)  then
                     call resetv (neq)
                  end if
                  if(ice.eq.0) then
                     call varchk (0, 0)
                  endif
                  call dual (1)
                  daynew = day
                  iad = abs(maxit) + 1
                  ntty_save = ntty
                  ntty = 2
                  call diagnostics(-1)
                  call diagnostics(1)
                  ntty = ntty_save
c
c count restarted timestep
c
                  nrestart_ts = nrestart_ts + 1
                  go  to  100
               end if
         
c**** stress routine calls ****

               dtotdm = dtot

               tassem = tyming(caz) - tassem

c**** solve for heat and mass transfer solution ****
               tas = tyming(caz) - tasii

               if (tas .gt. rnmax)  then

                  if (iout .ne. 0) write(iout, 6030) trim(nmfil(7))
                  if (iptty .gt. 0)  write(iptty, 6030) trim(nmfil(7))
 6030             format(/, 1x, '**** allotted time gone, terminating ',
     *                 'run : restart = ', a, ' ****')
c Make sure fin file is written if we are stopping
                  if (isave .ne. 0) call diskwrite
                  if(isalt.ne.0) call saltctr(21,0,0.0d00,0.0d00)
                  go  to  170
               
               end if

c zero variables used for current time step data
               amass = 0.0
               asteam = 0.0
               aener = 0.0
               tmav = 0.0
               prav = 0.0
               tmavg = 0.0
               pravg = 0.0
               pow = 0.0
               inflow_thstime=0.0
               inen_thstime=0.0

c**** calculate outflows and update coefficients for next time step ****
c qtoti and qtotei will have this value subtracted to get values for
c this time step only
               qtoti = qtot
               qtotei = qtote
               aw = awt
               ay = 1.0 - aw

c**** update variables and parameters ****

c gaz 12-30-99
c           call varchk (0, 0)
               call dual (1)

c**** calculate AIM factor ****

               call fimpf

c delete infinite reservoir nodes from mass balance calcs

               call bcon(2)

c check for steady state solution
               if(isteady.ne.0) then
                  call steady(1,0.,0.)
               endif

c
c correct mass correction
c
               call bcon(4)
c
c  calculate phase change information
c
               is_ch = 0
               nphase_liq = 0
               nphase_2 = 0
               nphase_gas = 0
               do i=1,n
                  if(ieos(i).eq.1) nphase_liq = nphase_liq + 1
                  if(ieos(i).eq.2) nphase_2 = nphase_2 + 1
                  if(ieos(i).eq.3) nphase_gas = nphase_gas + 1
                  if (irdof .ne. 13 .or. ifree .ne. 0) then
                     if(s(i).lt.1.0.and.so(i).ge.1.0) then
                        is_ch=is_ch +1
                     else if(s(i).gt.0.0.and.so(i).le.0.0) then
                        is_ch=is_ch +1
                     else if(s(i).ge.1.0.and.so(i).lt.1.0) then
                        is_ch=is_ch +1
                     else if(s(i).le.0.0.and.so(i).gt.0.0) then
                        is_ch=is_ch +1
                     endif
                  endif
               enddo
               is_ch_t = is_ch_t + is_ch
               if(l.ne.1) then
                dnphase_liq = nphase_liq - nphase_liq_0
                dnphase_2 = nphase_2 - nphase_2_0
                dnphase_gas = nphase_gas - nphase_gas_0
               else
                dnphase_liq = 0
                dnphase_2 = 0
                dnphase_gas = 0
               endif
               nphase_liq_0 = nphase_liq
               nphase_2_0 = nphase_2
               nphase_gas_0 = nphase_gas
c
c call thermo because the solver is overwriting the deni and denei arrays
c    
               if(istrs_coupl.gt.0.and.ico2.eq.0) then
                  call thermw(0)
               endif

c dtotdm is the current time step size in seconds
               do ja = 1, n

                  if (abs(ps(ja)) .gt. zero_t)  then
c qt = kg out - kg in
c qte = MJ out - MJ in
                     qt = qt + sk(ja) * dtotdm
                     qte = qte + qh(ja) * dtotdm
                  else
c check with gaz -- it seems like this should be qh(ja) or the code
c below is wrong
                     qte = qte + sk(ja) * dtotdm
                  end if
                  if (sk(ja) .gt. 0.0 .and. ps(ja).gt.0.0)  then
c there is outflow at this node
c toutfl = summation of kg out of system
c qtot = summation of kg out of system
                     toutfl = toutfl + sk(ja) * dtotdm
                     prav = prav + sk(ja) * phi(ja) * dtotdm
                     tmav = tmav + sk(ja) * t (ja) * dtotdm
                     qtot = qtot + sk(ja) * dtotdm
                  else if (sk(ja).lt.0.0 .and. ps(ja).gt.0.0) then
c there is inflow at this node
                     inflow_thstime=inflow_thstime+sk(ja)*dtotdm
                  endif
                  if (qh(ja).gt.0.0 .and. ps(ja).gt.0.0) then
c there is outflow of energy at this node
c teoutf = summation of MJ out of system
c qtote = summation of MJ out of system
                     teoutf = teoutf + qh(ja) * dtotdm
                     qtote = qtote + qh(ja) * dtotdm
                  else if (qh(ja).lt.0.0 .and. ps(ja).gt.0.0) then
c there is inflow of energy at this node
                     inen_thstime=inen_thstime+qh(ja)*dtotdm
                  end if
                  denht = denh(ja)
c                  to (ja) = t(ja)
                  denh (ja) = denh (ja) + deni (ja) * dtot
c asteam = kg of steam in system
c amass = kg of mass in system
c aener = energy in system
    
                  amass = amass + denh(ja) * volume(ja)
                  if(irdof.ne.13) then
                     asteam = asteam + dstm(ja)
                     deneht = deneh(ja)
                     deneh(ja) = deneh(ja) + denei(ja) * dtot
                     aener = aener + deneh(ja) * volume(ja)
                     denei(ja) = 0.0
                  endif
c                  pho (ja) = phi(ja)
                  if(abs(irdof).eq.14) then
                     denej (ja) = denj(ja)
                  else
                     denj (ja) = 0.0
                  endif        
                  deni (ja) = 0.0
                  if (irdof .ne. 13) then
                     if(ifree.ne.0) then
                        so (ja) = rlxyf(ja)
                     else
                        so (ja) = max(0.0d00,min(1.0d00,s(ja)))
                     end if
                  end if
    
               end do
               if(istrs_coupl.gt.0.or.istrs_coupl.eq.-3) then
c save flow residuals
                  call stressctr(17,0) 
c**** update stress arrays ****
c solve for displacements
c....... s kelkar 22 Aug 2012
                  if(istrs_coupl.ge.5) then
                     if(ifem.eq.1) conv_strain = fem_strain
                     if(iPlastic.eq.1) conv_pstrain = plastic_strain
                  endif
c........
                  if(istrs_coupl.eq.-3) then
                     istresscall = 1
                     call stress_uncoupled(3)
c update volume strains
                     call stressctr(6,0)
c update porosity
c s kelkar 22 Aug 2012. if pore_factor>0 this is done in 
c porosity_wrt_displacements which is called from from 
c bnswer and gensl_stress_coupled_3D  
                     if(pore_factor.eq.0) call stressctr(-7,0)
                  endif
c add displacements to total displacements	        
                  call stressctr(10,0)
c update displacements
                  call stressctr(12,0)
c update volume strains
                  call stressctr(-6,0)
c calculate stresses
                  call stressctr(13,0)	
               endif   
c calculate subsidence
               call subsidence(1)         
c update peaceman term for wellbore pressure
               if(isubwd.ne.0) call wellimped_ctr(1)
               do ja = 1,n
                  to (ja) = t(ja)
                  pho (ja) = phi(ja)
               enddo
c**** update co2 arrays ****
               call co2ctr (3)

c**** update component arrays or mixtures ****
c added call for balance error
               call icectr (9,0)
               call icectr (8,0)
               if(icarb.eq.1) call icectrco2 (9,0)
               if(icarb.eq.1) call icectrco2 (8,0)

c**** update ice arrays ****
c gaz 10-18-2001 call sice (3)

c**** call diagnostics  ****
c            call diagnostics(1)
c            call diagnostics(2)

c end if block for heat and mass transfer solution
            endif
c
c find max residuals for flow (H +M) solution
c
            call diagnostics(-1)
c            
c**** calculate velocities ****
            if(compute_flow .or. iccen .eq. 1) call veloc

c**** obtain concentration solution ****
            in3save = in(3)
            if(in(3).eq.0) then
               in(3) = irun + .0001
            end if
c save permeability and porosity if this can change with chemical transport     
c now just for salt
c gaz debug 091414 average 
c
c  save old ps and pnx
c
            if(isalt.ne.0) call saltctr(3,0,0.0d00,0.0d00) 
c            
            call concen (1,tscounter)
c
c average and updates new porosities and perms if necessary
c
            if(isalt.ne.0) call saltctr(4,0,0.0d00,0.0d00) 
c save new ps and perms for restart file
            if(isalt.ne.0) call saltctr(5,0,0.0d00,0.0d00) 
c
            in(3) = in3save 
        
c compute kg out of system this time step
            qtoti = qtot - qtoti
c compute MJ out of system this time step
            qtotei=qtote-qtotei
c compute MJ/s out of system this time step
            pow=qtotei/dtotdm

            if (qtoti .gt. 0.0)  then
               tmavg = tmav / qtoti
               pravg = prav / qtoti
c     qtotei = qtote - qtotei
c     pow = qtotei / dtotdm
            endif

c**** printout at specified intervals ****

            ntty = 0
            if (iptty .eq. 6)  ntty = 1
            if (ifinsh .eq. 0 .and. l .lt. nstep)  then
               iac = iac + 1
               if (iac .ge. iprtout)  then
                  ntty = 2
                  iac = 0
               end if
            else
               ntty = 2
            end if

            do im = 1, n
               pho(im) = phi(im)
            end do

            if(ice.eq.0) then

c**** calculate mass balance error ****
c tinfl = -total kg in
c teinfl = -total MJ in
               tinfl = qt - toutfl
               flmax = max(abs(tinfl), abs(toutfl))
               teinfl = qte - teoutf
               flemax = max(abs(teinfl), abs(teoutf))

               if (idof .le. 1)  flmax = 0.0
               dife = 0.0
               if (flemax .ne. 0.0)  dife = (aener-ame + qte)/flemax
c [(kg in system)-(kg originally in system)+(total kg left system -
c    total kg entered system)]/ max(total kg in, total kg out)
               if (flmax  .ne. 0.0)  difm = (amass-am0 + qt)/flmax

c mass balance for air
               call co2ctr (4)
c material balance for component mixtures
            endif

c**** call output routine ****
            if (ntty .gt. 0 .and. iout .ne. 0)  then
c retrieve flow residuals
               if(istrs_coupl.gt.0.or.istrs_coupl.eq.-3) then
                  call stressctr(18,0) 
               endif


!               print *, "The rank is ", rank

               if (rank == 0) then    ! only print out once 
            
                  call wrtout(tassem,tas,dabs(tinfl),dabs(teinfl),
     &                 dabs(inflow_thstime),dabs(inen_thstime),
     &                 is_ch,is_ch_t)

               end if 

               if(istrs_coupl.gt.0.or.istrs_coupl.eq.-3) then
c output  displacements and stresses 
                  call stressctr(11,0) 
               endif 

            end if

c**** call history plot ****

            if (hist_flag) then
               call plot_new (1, dabs(tinfl), dabs(teinfl),
     &              dabs(inflow_thstime), dabs(qtoti))
c     &           dabs(inflow_thstime), dabs(inen_thstime))
                if(istrs.ne.0) call  stressctr(14,0)
            else
               call plot (1, tmavg, pravg)
            end if
c**** call wellbore pressures
           if(isubwd. ne.0) call wellimped_ctr(4)
c check for steady state solution
            if(isteady.ne.0) then
               call steady(2,dabs(inflow_thstime),dabs(inen_thstime))
            endif

c     Move call to before the contour output routines
c     BAR 7-20-99

            if(sptrak) then
c........... s kelkar nov 13, 02....................
c if freez_time is gt.0, then ptrac3 is called only at the end of the 
c flow calculations, and for a velocity frozen at that time,  
c particle tracks are calculated for freez_time days
               if(freez_time.eq.0.) then
                  call ptrac3
               endif
            endif
c........................................

c**** call contour plot ****
            dayscn = dayscn + day

            if (icontr .ne. 0 .or. mod(l, ncntr) .eq. 0 .or.
     *           dayscn .ge. abs(contim))  then

!               if(contim.ge.0) then
               call contr (1)
               call contr (-1)
               call river_ctr(6)
               call active_nodes_ctr(4)
!               else
!                  call contr_days (1)
!                  call contr_days (-1)
!               endif 
               if (allocated(itc) .and. nicg .gt. 1) then
                  if (itc(nicg-1).gt.0) then
!                 call disk (1)
                     if (isave .ne. 0) call diskwrite
                     if(isalt.ne.0) call saltctr(21,0,0.0d00,0.0d00)
                     if (iflxn .eq. 1) call flxn
                  endif
               else
!              call disk (1)
                  if (isave .ne. 0) call diskwrite
                  if(isalt.ne.0) call saltctr(21,0,0.0d00,0.0d00)
                  if (iflxn .eq. 1) call flxn
               end if
               call pest(1)
                if(isalt.ne.0) call saltctr(6,0,0.0d00,0.0d00)              
               if (ifinsh .ne. 0)  ex = .TRUE.
               if (dayscn .ge. abs(contim)) dayscn = 0.

            end if


            if (ifinsh .eq. 1) then
c standard simulation finish
               istea_pest = 0
	         
               go  to  170
            else if (ifinsh .eq. 2) then

! Make sure last steady state time step is output
               if (hist_flag) then
                  call plot_new (1, dabs(tinfl), dabs(teinfl),
     &                 dabs(inflow_thstime), dabs(qtoti))
               end if
               call contr (1)
               call contr (-1)
               call river_ctr(6)
               call active_nodes_ctr(4) 
c finished the steady state simulation, now doing transient
               days = 0.0
               istea_pest = 0
               call pest(1)
               istea_pest = 1
               call flow_boundary_conditions(4)
               go to 999
            endif


c EHK check for kill file
        inquire(file='kill.file',exist=die)
        if(die) goto 170

        end do      ! End major time step loop 




c**** write solution to plot tapes ****

         l = l - 1

 170     continue


c ******************* end major time step loop ****************


! -----------------------------------------------------------
!               Finalize / Deallocate PETSc memory
! -----------------------------------------------------------
      call petsc_finalize()

! ----------------------------------------------------------- 
 
c     check for steady state solution
         if(isteady.ne.0) then
            ntty = 2
            call steady(3,dabs(inflow_thstime),dabs(inen_thstime))
         endif

c......     s kelkar nov 13, 02........
c if freez_time is gt.0, then ptrac3 is called only at the end of the 
c flow calculations, and for a velocity frozen at that time,  
c particle tracks are calculated for freez_time days
         if(sptrak) then
            if(freez_time.gt.0.) then
               call ptrac3
            endif
         endif
c................................................

         if(sptrak) then
            close(isptr1)
            close(isptr2)
            close(isptr3)
            if (pod_flag) then
! Call to write out derivatives for model reduction basis functions
               call pod_derivatives
            end if
         endif
c call subsidence for last time, iflg = 2
         call subsidence(2)

c     
c     gaz 1-6-2002
c     printout submodel boundary conditions if necessary
c     
         if(isubbc.ne.0) call submodel_bc(2)
c     
c     For a rip simulation, write avs output info every
c     contim_rip days
         if(tims.ge.contr_riptot .and. ripfehm .ne. 0) then
            contr_riptot = contr_riptot + contim_rip
!         if(contim.ge.0) then
            call contr (1)
            call contr (-1)
            call river_ctr(6)
            call active_nodes_ctr(4)
!         else
!            call contr_days (1)
!            call contr_days (-1)
!         endif
         end if

      if(die) then
         if (iout .ne. 0) then
         write(iout, '(a40)') 'Kill file present; simulation terminated'
         endif
         if (iptty .gt. 0) then
         write(iptty, '(a40)')'Kill file present; simulation terminated'
         endif
      endif


      ! only print out once for MPI run
      if (rank == 0) then

         if (iout .ne. 0) write(iout, 6040)  days, l
         if (iptty .gt. 0)  write(iptty, 6040)  days, l
 6040    format(//, 1x, 'simulation ended: days ', 1pg30.23, 
     *        ' timesteps ', i5)
     
      end if 
 
c
c calculate final stress field and displacements
c output contour information
c 
         if(istresscall.eq.0.and.ihms.eq.-4)then
           istrs_coupl = -2
         endif 

         call stress_uncoupled(2)

c zvd 30-Jun-10
c Move call to disk_write and contr after call to stress_uncoupled
         nsave = 1

         if (.not. ex) then
            if(in(1).eq.0) then
               if( tscounter .eq. 1 .or. in(1) .eq. 0
     2              .or. tims .eq. tims_save) then
                  if (allocated(itc) .and. nicg .gt. 1) then
                     if (itc(nicg-1).gt.0) then
!                    call disk (nsave)
                        if (isave .ne. 0) call diskwrite
                        if(isalt.ne.0) call saltctr(21,0,0.0d00,0.0d00)
                        if (iflxn .eq. 1) call flxn
                     endif
                  else
!                 call disk (nsave)
                     if (isave .ne. 0) call diskwrite
                     if(isalt.ne.0) call saltctr(21,0,0.0d00,0.0d00)
                     if (iflxn .eq. 1) call flxn
                  end if
!               if(contim.ge.0) then
                  if(istrs_coupl.ne.-2.and.istrs_coupl.ne.-1) then
! contr will be called below in stress_uncoupled for a stress solution
c contr was called in stress_uncoupled above
                     call contr (1)
                     call contr (-1)
                  end if
                  call river_ctr(6)    
                  call active_nodes_ctr(4)           
!               else
!                  call contr_days (1)
!                  call contr_days (-1)
!               endif
                  istea_pest = 0
                  call pest(1)
c     **** call pest to calculate sensitivities if necessary
                  call pest(2)
               end if
            end if
         end if
c     New convention is to make days the - of its value to
c     write to the output history file, then change it back
c     after calling plot. days needs to be correct if the
c     code is being called in "rip" mode

         days = -days
c      
         if (hist_flag) then
            call plot_new (1, dabs(tinfl), dabs(teinfl),
     &           dabs(inflow_thstime), dabs(inen_thstime))
! Add call to plot_new so all dummy variables will be deallocated
! after final history data is output
            call plot_new (2, dabs(tinfl), dabs(teinfl),
     &           dabs(inflow_thstime), dabs(inen_thstime))
            if(istrs.ne.0) call  stressctr(14,0)
         else
            call plot (1, tmavg, pravg)
         end if

c     Change it back

         days = -days

         ! only print out once for MPI run
         if (rank == 0) then  

              if (iout .ne. 0) write(iout, 6041) itotal,itotals
              if (iptty .gt. 0) write(iptty, 6041)  itotal,itotals
 6041    format(//, 1x, 'total N-R iterations = ', i10
     &        ,/,1x, 'total solver iterations = ', i10)

              if (iout .ne. 0) write(iout, 6042) tyming(caz) - tasii 
              if (iptty .gt. 0) write(iptty, 6042) tyming(caz) - tasii
 6042    format(//, 1x, 'total code time(timesteps) = ', f13.6)

              call dated (jdate, jtime)

              if (iout .ne. 0) write(iout, 6052)  verno, jdate, jtime
              if (iptty .gt. 0) write(iptty, 6052)  verno, jdate, jtime

         end if 


         if (ripfehm .eq. 0) then
            close (inpt)
            if (iout .ne. 0) close (iout)
         end if
 6052    format(//, 1x, '****----------------------------------------', 
     *        '-----------------****', 
     *        /, 1x, '**** This program for ', 35x, '    ****', 
     *        /, 1x, '****   Finite Element Heat and Mass Transfer ', 
     *        'in porous media ****', 
     *        /, 1x, '****-------------------------------------------', 
     *        '--------------****', 
     *        /, 1x, '**** ', 12x, '  Version  : ', a30, ' ****', 
     *        /, 1x, '**** ', 12x, '  End Date : ', a11, 18x, '  ****', 
     *        /, 1x, '**** ', 12x, '      Time : ', a8, 20x, '   ****', 
     *        /, 1x, '****-------------------------------------------',
     *        '--------------****')

c     Add call to routine to transfer particle information to out(i)
c     array for rip simulations

         if(n_input_arguments .ne. 0) then
            if(int(in(n_input_arguments+1)).ne.0
     2           .and.int(in(n_input_arguments+4)).ne.0)
     3           call loadoutarray
         end if

c RJP 04/10/07 this is for co2 properties table look-up
         call interpolation_arrays_deallocate()
c
c gaz 031314 split hexes into tets
c
      if (sv_hex_tet.and.ns.eq.8) then
       rewind incoor
       call incoord
       ivf_sv = ivf
       ivf = -1
       call split(0)
       ivf = ivf_sv
      endif
c
c     no more method = 4
c     elseif(method.eq.4) then

c     routine computes fluid flux values exiting each output buffer
c     Routine no longer needed. These are computed in part_track
c     and stored in the array idflux

c     call computefluxvalues

c     End if around the part of code for performing the simulation
      else
         continue
      end if

c     Add statement to set return flag

c     For the RIP version, set so that no error
c     is passed back to rip
      state = 0

      return

      contains

c     Subroutine riptime - scope is local to fehmn

      subroutine riptime
      implicit none
      logical used
      real(8) daystmp
      integer ncall, ichloc
      character*23 string_call
      character*10 ffname
      character*5 ch5

      wtrise_flag = .false.
      if(in(1).ne.0.) then
         tims = abs(in(1))*365.25

c     Read in new flow field if the input flag has changed from the 
c     Previous time step

c         if(int(in(2)).ne.flowflag) then
c*** water table rise modification
         if( (int(in(2)).ne.flowflag) .or.
     &        (abs(in(7)-water_table_old).gt.1.d-6) ) then
c*** water table rise modification

c     Flag adjusted to tell particle tracker that new
c     flow field is being read in

            tscounter = -tscounter

c*** water table rise modification
c zvd 21-Jul-08 Always make water table adjustment when a new flow 
c field is read
            water_table = in(7)
            wtrise_flag = .true.
c*** water table rise modification

c     Define file name except for the number

            ffname = 'ff    .ini'

c     Create number such that 1 is 10001, 2 is 10002, etc.

            ncall = int(in(2)) + 10000
c     Write then number to ch5 character string

            write(ch5,'(i5)') ncall
c     Place the number in the empty space of ffname, one
c     character at a time. Start w/ digit 2, so that the 
c     file name for #1 is ff00001.ini, etc.

            do ichloc = 2, 5
               ffname(1+ichloc:1+ichloc)=ch5(ichloc:ichloc)
            end do

            if (iptty .gt. 0) write(iptty,*) 
     2           'Reading a new restart file: ', ffname
            write(iptty,*) 
     2           'Reading a new restart file: ', ffname


c     Check to see if file is already open, if so get file
c     number and rewind file

            used = .false.
            inquire(file=ffname,opened=used)
            if(used) then
               inquire(file=ffname, number = iread)
               rewind (iread)
            else


               do i = 1,80
                  filename(i:i) = ' '
               end do
               do i = 1, 10
                  filename(i:i) = ffname(i:i)
               end do
               iread = open_file(filename,'old')


            end if
            nmfil(6) = ''
            nmfil(6) = ffname

            daystmp = days
!            call disk(0)
            if (iread .gt. 0) call diskread
            days = daystmp
! zvd 22-Mar-02
! File is now closed after data is read in call disk(0)
!            close(iread)


         end if
      elseif(in(3).ne.0) then
c bhl 2005
c     Read in new flow field if the input flag has changed from the 
c     Previous time step

c         if(int(in(2)).ne.flowflag) then
c*** water table rise modification
         if( (int(in(2)).ne.flowflag) .or.
     &        (abs(in(7)-water_table_old).gt.1.d-6) ) then
c*** water table rise modification
            write(ierr,*)'in(2):',in(2)
            write(ierr,*)'flowflag:',flowflag


c     Flag adjusted to tell particle tracker that new
c     flow field is being read in

            tscounter = -tscounter

c*** water table rise modification
            if(abs(in(7)-water_table_old).gt.1.d-6) then
               water_table = in(7)
               wtrise_flag = .true.
            else
               wtrise_flag = .false.
            end if
c*** water table rise modification

c     Define file name except for the number

            ffname = 'ff    .ini'

c     Create number such that 1 is 10001, 2 is 10002, etc.

            ncall = int(in(2)) + 10000
c     Write then number to ch5 character string

            write(ch5,'(i5)') ncall
c     Place the number in the empty space of ffname, one
c     character at a time. Start w/ digit 2, so that the 
c     file name for #1 is ff00001.ini, etc.

            do ichloc = 2, 5
               ffname(1+ichloc:1+ichloc)=ch5(ichloc:ichloc)
            end do

            if (iptty .gt. 0) write(iptty,*) 
     2           'Reading a new restart file: ', ffname
            write(iptty,*) 
     2           'Reading a new restart file: ', ffname


c     Check to see if file is already open, if so get file
c     number and rewind file

            used = .false.
            inquire(file=ffname,opened=used)
            if(used) then
               inquire(file=ffname, number = iread)
               rewind (iread)
            else


               do i = 1,80
                  filename(i:i) = ' '
               end do
               do i = 1, 10
                  filename(i:i) = ffname(i:i)
               end do
               iread = open_file(filename,'old')


            end if
            nmfil(6) = ''
            nmfil(6) = ffname


            daystmp = days
!            call disk(0)
            if (iread .gt. 0) call diskread
            days = daystmp
! zvd 22-Mar-02
! File is now closed after data is read in call disk(0)
!            close(iread)


         end if

c bhl 2005
c     tims = overf
c     We are here if it is a GoldSim run and in(1) = 0
c     Here we want the code to take a very small time step, essentially
c     zero. This is done by the user setting a low value of daymin
c     in the ctrl macro
         tims = daymin
c     Run batch file on first timestep
         if(abs(tscounter).eq.1) then
            string_call(1:14) = 'fehmn_ts0.bat '

            ncall = int(in(2)) + 10000

c     Write then number to ch5 character string

            write(ch5,'(i5)') ncall
c     Place the number in the empty space of string, one
c     character at a time. Start w/ digit 15

            do ichloc = 2, 5
               string_call(13+ichloc:13+ichloc)=ch5(ichloc:ichloc)
            end do
            string_call(19:19) = ' '

c     Do the same with in(3)

            ncall = int(in(3)) + 10000

c     Write then number to ch5 character string

            write(ch5,'(i5)') ncall
c     Place the number in the empty space of string, one
c     character at a time. Start w/ digit 20

            do ichloc = 2, 5
               string_call(18+ichloc:18+ichloc)=ch5(ichloc:ichloc)
            end do


            if(iptty.ne.0) then
               write(iptty,*) 'Running fehmn_ts0.bat'
               write(iptty,*)
     2              'Calling arguments are ',
     3              string_call(15:18),' ',string_call(20:23)
            end if
            call system(string_call(1:23))
         end if
      end if
      return
      end subroutine riptime
c     Subroutine loadoutarray - scope is local to fehmn

      subroutine loadoutarray
      implicit none
      real*8, allocatable :: out_save(:)
      real*8, allocatable :: time_dump(:)
      integer ispecies
      integer number_of_species
      integer ns2,izones
      integer nflow_frac, number_of_zones,indexout,indexmzone
      integer add_spots, add_spots2
      integer index_temp
      real*8 cur_time, prev_time, del_time
      real*8 :: cur_time_save = 0.
      save out_save, time_dump, cur_time_save

cHari 3/1/07    
c     Before V 2.23, in(4) was the correct index, but now it is 
c     in(8) because two random number seeds, a flag, and the 
c     water table elevation were added to the interface before M_fine
c     Now, M_fine is in(8). The number added to get to index_N_large 
c     is now 9 instead of 5
c     BAR 2-9-2005

      index_N_large=int(in(8))*2+9

c     As of V 2.23, we now use a flag to decide whether there are
c     nflow_frac inputs of fracture fractional flow data to
c     handle, or if the array skips that input and proceeds
c     directly to number_of_species. This flag is used in the 
c     if block below. BAR 2-9-2005

      if(int(in(6)).eq.1) then
c     flow fraction data exists
         index_temp = index_N_large+int(in(index_N_large))+1
         nflow_frac = int(in(index_temp))
         index_in_species=index_N_large+int(in(index_N_large))+
     2        nflow_frac+2
         number_of_species = int(in(index_in_species))
      else
c     no flow fraction data
         index_in_species=index_N_large+int(in(index_N_large))+1
         number_of_species = int(in(index_in_species))
      end if
c     index_in_species=index_N_large+int(in(index_N_large))+1
c     number_of_species = int(in(index_in_species))
CHari if number of species is > 45 then we assume that flow
CHari fractions are in use and we are really being passed nspecies*nlarge

c     if(number_of_species.gt.45)then
c     nflow_frac = number_of_species
c     index_in_species=index_N_large+int(in(index_N_large))+
c     2      nflow_frac+2
c     number_of_species = int(in(index_in_species))
c     endif
      index_in_flag=index_in_species+1
      index_out_buffer=index_in_flag+2
      number_of_outbuffers = int(in(index_out_buffer))


c     Number of output buffers is the total number
c     The water table is split into zones, and for
c     dual perm. problems there is a fracture and matrix
c     exit buffer for each zone. Therefore, rip will pass
c     the total number of buffers, and we divide by 2 to
c     get the number of output zones

      if(idpdp.ne.0) then
         number_of_zones = number_of_outbuffers/2
      else
         number_of_zones = number_of_outbuffers
      end if

c     Allocate local array that saves the previous out array
c     time_dump is the time at which the previous particles from
c     this region have been dumped to the out array

      if(.not.allocated(out_save)) then
         allocate(out_save(number_of_outbuffers*number_of_species))
         allocate(time_dump(number_of_outbuffers*number_of_species))
         out_save = 0
         time_dump = 0
      end if

c     cur_time = current time in years
c     prev_time = time at previous time step (years)


      cur_time = days/365.25
      del_time = dtot/(365.25*86400.)
      prev_time = cur_time - del_time

c     Reinitialize the out_save and time_dump arrays
c     when a new realization is initiated (i.e. the cur_time
c     "clock" is reset to a low value

      if(cur_time_save.gt.cur_time) then
         out_save = 0
         time_dump = 0
      end if
      cur_time_save = cur_time


      if(out_flag.eq.0) then
         add_spots = 0
         indexout = 0
      else
         add_spots = 2*number_of_outbuffers*number_of_species
         add_spots2 = number_of_outbuffers*number_of_species
         indexout = add_spots
      end if

      do izones = 1, number_of_zones

c     Buffer associated with the fractures is done in the loop below

         do ispecies = 1, number_of_species
            indexout = indexout+1
            out(indexout)=0.
            if(ispecies <= nspeci)then

c     if block makes sure time_dump is set properly if this is the first
c     time step in a realization other than the first realization

               if(tscounter.eq.1) then
                  time_dump(indexout-add_spots) = prev_time
               end if
               
c     This if block checks that at least 2 particles have exited since
c     the last time particles have exited the system. If they have,
c     out array is updated and the dump time is set to the current time
c     The normal conversion to grams is gmol*pcount. This new
c     method accounts for the gradual trickling out of particles by 
c     multiplying this value by del_time/(cur_time-time_dump(indexout))
c     to account for the exiting of two particles over more than a single
c     time step. The out value remains at the previous value until at 
c     least 2 particles have left, then computes the average mass exiting
c     at the time step.

c     In part_track, the pcount value only gets reset when 2 particles
c     have accumulated in the exit bin.

               
                                !per the request of Dave Sevougian, cli removed the smoothing scheme implement
                                !below.

                                !cli      if(pcount(izones,ispecies).gt.0) then
               out(indexout) = del_time*gmol(ispecies)*
     2              pcount(izones,ispecies)/
     3              (cur_time-time_dump(indexout-add_spots))
               time_dump(indexout-add_spots) = cur_time
                                !cli      else

c     Set output mass to the previous value until at least 2 particles
c     have accumulated at the output bin

                                !cli         out(indexout) = out_save(indexout-add_spots)
                                !cli added this statement so that we do not smear the first non-zero mass output
                                !cli         if(out(indexout).eq.0.)then
                                !cli          time_dump(indexout-add_spots) = cur_time
                                !cli         end if
                                !cli      end if

c     Store new value of out in the out_save array

               out_save(indexout-add_spots) = out(indexout)
            endif

c     If output for max and avg concentration are to be passed back,
c     do that here

            if(out_flag.ne.0) then
c     average concentration
               out(indexout-add_spots) = idcavg(izones,ispecies)

c     maximum concentration
               out(indexout-add_spots2) = idcmax(izones,ispecies)

            end if

         end do
         
         

         if(idpdp.ne.0) then
c     Buffer associated with the matrix is done in the loop below
c     First, compute the index of the pcount array corresponding to
c     the matrix exiting particles. For this, we need to realize that
c     the first (number_of_zones+1) are the fracture data. The 1
c     is because the remaining particles not leaving any of the zones
c     are in the pcount array also. However, these are not passed
c     through in the out array in this version of the code.

            indexmzone = number_of_zones+1+izones

            do ispecies = 1, number_of_species
               indexout = indexout+1
               out(indexout)=0.
               if(ispecies <= nspeci)then

c     if block makes sure time_dump is set properly if this is the first
c     time step in a realization other than the first realization

                  if(tscounter.eq.1) then
                     time_dump(indexout-add_spots) = prev_time
                  end if

c     See comments above for fracture nodes for an explanation of this
cli   removed confactor from the out() calculations because, pcount is 
cli   already in mols.         

                                !per the request of Dave Sevougian, cli removed the smoothing scheme
                                !implemnted. 09-05-03

                                !cli         if(pcount(indexmzone,ispecies).gt.0) then
                  out(indexout) = del_time*gmol(ispecies)*
     2                 pcount(indexmzone,ispecies)/
     4                 (cur_time-time_dump(indexout-add_spots))
                  time_dump(indexout-add_spots) = cur_time
                                !cli         else

c     Set output mass to the previous value until at least 2 particles
c     have accumulated at the output bin

                                !cli            out(indexout) = out_save(indexout-add_spots)
                                !cli added this statement so that we do not smear the first non-zro mass output
                                !cli            if(out(indexout).eq.0.)then
                                !cli               time_dump(indexout-add_spots) = cur_time
                                !cli            end if
                                !cli         end if

c     Store new value of out in the out_save array

                  out_save(indexout-add_spots) = out(indexout)
               endif

c     If output for max and avg concentration are to be passed back,
c     do that here

               if(out_flag.ne.0) then
c     average concentration
                  out(indexout-add_spots) = idcavg(indexmzone,ispecies)

c     maximum concentration
                  out(indexout-add_spots2) = idcmax(indexmzone,ispecies)

               end if


            end do
         end if

      end do


      return
      end subroutine loadoutarray

c     computefluxvalues - passes fluid mass fluxes to rip
c     scope is local to fehmn

      subroutine computefluxvalues

      implicit none
      real*8 fluxfrac, fluxmat
      integer iznum, inode, nmedia, indexarray
      integer number_of_zones

c     number_of_outbuffers = int(in(n_input_arguments+4))

c     Number of output buffers is the total number
c     The water table is split into zones, and for
c     dual perm. problems there is a fracture and matrix
c     exit buffer for each zone. Therefore, rip will pass
c     the total number of buffers, and we divide by 2 to
c     get the number of output zones

      if(idpdp.ne.0) then
         number_of_zones = number_of_outbuffers/2
         nmedia = 2
      else
         number_of_zones = number_of_outbuffers
         nmedia = 1
      end if

c     zero out all flux values

      out(1:number_of_outbuffers) = 0.

c     We can only do this calculation if the particle tracking
c     has begun, otherwise idzone is not yet allocated and
c     is undefined

      if(allocated(idzone)) then

c     Loop over each output zone - both fracture and matrix
c     fluxes get computed and stored at the same time

         do iznum = 1, number_of_zones

c     indexarray is the position in the out array for the fracture flux
c     for this zone. indexarray+1 is the corresponding matrix
c     matrix

            indexarray = nmedia*(iznum-1)+1

c     Need to loop through each node in the modelto see if it is
c     in this output zone

            do inode = 1, neq

               if(izonef(inode).eq.idzone(iznum)) then

c     This is one of the fracture output zones - add to flux

                  fluxfrac = a_axy(nelmdg(inode)-neq-1)

c     Filter out the extremely large out flows, which denote
c     nodes connected to no other nodes in the grid

                  if(fluxfrac.lt.1.e8) then
                     out(indexarray) = out(indexarray) + fluxfrac
     2                    * 31557.600
                  end if

c     If dual perm., add matrix output flux to the next position in 
c     the out array

                  if(idpdp.ne.0) then
                     fluxmat = a_axy(nelm(neq+1)-neq-1+
     2                    nelmdg(inode)-neq-1)
                     if(fluxmat.lt.1.e8) then
                        out(indexarray+1) = out(indexarray+1) + fluxmat
     2                       * 31557.600
                     end if
                  end if

               end if

            end do


         end do

      end if
      return

      end subroutine computefluxvalues

c     Routine called fehmn for GoldSim runs

      subroutine ingold
      implicit none
      logical more
      character(80) single_line
      integer print_flag, iread1, iread2

      used = .true.
      iread = 1
      do while(used)
         inquire(unit = iread, opened = used)
         if(.not.used) then
            
c     Open file to read
            
            used = .false.
            open(iread,file = 'fehmn_real.bat')
            iread1=iread
         else
            used = .true.
            iread = iread + 1
         end if
      end do

      used = .true.
      iread = 1
      do while(used)
         inquire(unit = iread, opened = used)
         if(.not.used) then
            
c     Open file to read
            
            used = .false.
            open(iread,file = 'fehmn_ts0.bat')
            iread2=iread
         else
            used = .true.
            iread = iread + 1
         end if
      end do

      more = .true.
      print_flag = iread1
      do while(more)
         read(1,'(a80)',end=1000) single_line
         if(single_line(1:4).ne.'ts0') then
            write(print_flag,*) single_line
         else
            print_flag = iread2
         end if
      end do
 1000 more = .false.
      close(iread1)
      close(iread2)

      return
      end subroutine ingold

      end subroutine fehmn

