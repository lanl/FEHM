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
      use comdi
      use comdti
      use comei
      use comevap, only : evaporation_flag
      use comfi, only : qtc, qtotc
      use comflow, only : a_axy
      use compart
      use comriv
      use comrtd, only : maxmix_flag
      use comrxni
      use comsi
      use comsk, only : save_omr
      use comsplitts
      use comsptr
      use comuserc, only : in, iaunit
      use comwt
      use comxi
      use davidi
      use property_interpolate

      implicit none

c     These are PC attributes used as compiler directives. They
c     should be set as follows for the PC-RIP version of fehm.

c     If the stand-alone pc version of the code is being used, then
c     the fehmn attribute line should be omitted and the method
c     variable should be passed by reference.

c     For UNIX versions, these lines are ignored as comments.
!DEC$ ATTRIBUTES dllexport, c :: fehmn
!DEC$ ATTRIBUTES value :: method
!DEC$ ATTRIBUTES reference :: state
!DEC$ ATTRIBUTES reference :: ing
!DEC$ ATTRIBUTES reference :: out

      integer(4) method, state
      real(8) ing(*), out(*)
      real(8) qcout_old(100)
      
c     irun is a counter for each realization in a multiple simulation
c     run of fehm. It is initialized to 0 in comai

      character*80 filename
      character*40 logfile
      character*40 errfile
      character*40 ainfile
	character*2  clus(11)
      
      integer open_file, ifail

      logical used
      real*8 tims_save, day_saverip, in3save
      real*8 deneht, denht, eskd, enthp, flemax, flmax, prav, 
     *     pravg, tajj, tas, teinfl, tinfl, tmav, tmavg, tyming
      real(8) :: inflow_thstime = 0., inen_thstime = 0.
      real(8) :: contr_riptot = 1.0d+30
      real(8) :: tasii = 0., tassem = 0.
      real*4 caz(2)
      real*8 prop,dpropt,dpropp,p_energy
cSPC
      real*8 jjj, goldtime

      logical it_is_open, intfile_ex
      integer im, ja, mi, i, ii, j, k
      integer :: ichk = 0, tscounter = 0, flowflag = 0
      integer number_of_outbuffers, jpr
      integer :: n_input_arguments = 0
      integer index_N_large, index_in_species
      integer index_in_flag, index_out_buffer
      integer is_ch, skipflag
      integer :: is_ch_t = 0
      integer :: out_flag = 0
      integer :: size_old = 0
	integer jumpflag, loopflag, loopstart
cSPC
      integer cluster
      integer sflag, size_of_in
cSPC  add Aug 27, 2009 Goldsim species number for initialization
      save inflow_thstime, inen_thstime
      save qcout_old, size_of_in
      save jumpflag, loopflag, loopstart, logfile, ainfile, errfile
      
      save flowflag, ichk, tassem, tasii, tscounter,
     &     contr_riptot, tims_save, day_saverip, in3save,
     &     n_input_arguments, cluster, clus, skipflag
 
c-------------------------------------------------------------------------------- 
c - - - - - - - - - - - - - - - - - - - - PHS 1/26/2012    SKIP      Method = 1
c           Skip for no GW pathway                         SKIP
 
      if((method.eq.1).AND.(ing(4).EQ.667)) then
	
	  irun = 1
	  skipflag = 1
	  jjj = int(ing(int(7+2*ing(7)+1)))
	  do i = 1,jjj
	    out(i) = 0.0
	  end do 

	  goto 9999

      end if
	  
c------------------------------------------------------ Write Method State   
    
      if((method.NE.1).AND.(irun.NE.0).AND.(skipflag.NE.1)) then
	  write(iaunit,*) '--------------------------'
	  write(iaunit,*) 'method, state' , method, state
	  write(iaunit,*) '--------------------------'
      end if
			           
c--------------------------------------------------------------------------------
c---------------------------------------------                        Method = 2
      if(method.eq.2) then

         out(1) = 3.0

	end if

c--------------------------------------------------------------------------------
c - - - - - - - - - - - - - - - - - - - - - - - -                      Method = 3 
      if(method.eq.3) then

         ripfehm = 1
c - - - - - fehmn.gold has num of incoming and out (minus species mass) 	    
	   inquire(7777,opened=it_is_open)
         if(.not.it_is_open) then
		  open(7777,file='fehmn.gold', status='unknown')
         end if
	   rewind(7777)
	          
          read(7777,*) n_input_arguments, out_flag
	   close(7777)
                  
         out(1) = n_input_arguments
         out(2) = 0
      
	end if
c--------------------------------------------------------------------------------
c - - - - - - - - - - - - - - - - - - - - - - - -                     Method = 99     
      if(irun.NE.0.AND.method.EQ.99) then
c     Cleanup - close  files at the end of the realizations
         if (isave .ne. 0) call diskwrite     
c      Close log err and a_in_array.txt files
         
	  inquire(unit = ierr,opened=it_is_open)
        if(it_is_open) close (ierr)    

	  inquire(unit = iptty,opened=it_is_open)
        if(it_is_open) close (iptty)

	  inquire(unit = iaunit, opened=it_is_open)
        if(it_is_open) close (iaunit)	 

c     Release all dynamic array memory at the beginning of a realizatoin
	   call releasemem
	   irun = 0

c       Phil removed much YMP stuff in the following section
      end if

c-------------------------------------------------------------------------------- 
c - - - - - - - - - - - - - - - - - - - - PHS 4/20/2011                Method = 0        
      if(method.eq.0) then
   
        continue
                                        
      end if       ! END IF  Method = 2,3,99,0 
          
c-------------------------------------------------------------------------------- 
c - - - - - - - - - - - - - - - - - - - - PHS 4/20/2011     time=0   Method = 1  
c                                                   Beginning of new Realization      
      if((method.eq.1).AND.(ing(1).EQ.0)) then    

c     -------------------------------------------------------
c      GoldSim Initialization Stuff	

c       Allocate memory for the in array and 
c         load values of IN for Time=0

        size_of_in = n_input_arguments + 4 + ing(n_input_arguments + 1)
        
        if (not(allocated(in))) then
          allocate(in(size_of_in))
	  end if

	  do j = 1,size_of_in
          in(j)=ing(j)
	  end do

c  - - - - - PHS moving open err,log,a_inarray files to here.4/20/11

        cluster = int(in(6))

	  clus(1) = '01'
	  clus(2) = '02'
	  clus(3) = '03'
	  clus(4) = '04'
	  clus(5) = '05'
	  clus(6) = '06'
	  clus(7) = '07'
	  clus(8) = '08'
	  clus(9) = '09'
	  clus(10) = '10'
	  clus(11) = '11'

      
        logfile = 'clusterXX/fehmXX.log'
        errfile = 'clusterXX/fehmXX.err'
        ainfile = 'clusterXX/a_inXX_array.txt'
        logfile(8:9) = clus(cluster)
	  errfile(8:9) = clus(cluster)
	  ainfile(8:9) = clus(cluster)
        logfile(15:16) = clus(cluster)
	  errfile(15:16) = clus(cluster)
	  ainfile(15:16) = clus(cluster)

        
        
        ierr = 6555 + cluster
        inquire(unit = ierr,opened=it_is_open)
        if(.not.it_is_open) then
          open (ierr, file = errfile, status='unknown')
        end if

        iptty = 6666 + cluster
	  inquire(unit = iptty,opened=it_is_open)
        if(.not.it_is_open) then
	    open(iptty,file = logfile, status='unknown')
         end if

        iaunit = 6777 + cluster
	  inquire(unit = iaunit, opened=it_is_open)
        if(.not. it_is_open) then
	    open(iaunit,file = ainfile,status='unknown')
	  end if
	   
	  write(ierr,*) 'cluster ',cluster , 'method = ',method,
     x               'realiz ',int(in(3))
        write(iaunit,*)  'cluster ',cluster , 'method = ',method,
     x               'realiz ',int(in(3))
	  write(iptty,*) 'cluster ',cluster , 'method = ',method,
     x               'realiz ',int(in(3))

        write(iaunit,*) '-----------------------------------------'
	  write(iaunit,*) 'IN array values  -  size of in' , size_of_in
	  do j = 1,7+int(in(7))
          write(iaunit,720) j,int(in(j))
	  end do 
	  do j = 8+int(in(7)) , 7+2*int(in(7))
	    write(iaunit,721) j, in(j)
	  end do
	  do j = 8+2*int(in(7)),size_of_in-int(in(n_input_arguments + 1))
	    write(iaunit,720) j, int(in(j))
	  end do
	  do j = size_of_in-int(in(n_input_arguments+1))+ 1, size_of_in
	    write(iaunit,721) j, in(j)
	  end do 
        write(iaunit,*) '-----------------------------------------'

 720    format(I5,I5)
 721    format(I5,G12.3)


c      SPC add Aug 27, 2009, avoid odd out array values
   
	  do j = 1, in(8+in(7)*2)
           out(j) = 0.0
	  end do

c     PHS 4/29/2011 subroutine to get the name of the fehmn.files file
         
        call namefiles

c        End GoldSim Init Specific Stuff 
c     -------------------------------------- 
             
c     --------------------------------------------		     
c     Initialize  FEHM  stuff  
         tscounter = 0
         flowflag = 0
         tassem = 0.0

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
         call iofile (ichk)

	   write(iptty,*) 'After iofile call  cluster =',cluster

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
c gaz 050809 moved to startup
c
c calculate initial stress field and displacements
c 
c         call stress_uncoupled(1)
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
               end if
               if(.not.random_flag) then
c                  if(allocated(sx)) deallocate(sx)
c                  if(allocated(istrw)) deallocate(istrw)
               end if
            end if
c     call ptrac1
         endif

c        -------------------------------
c        Set GoldSim specific values for 
c        begining of new realiazation         
	   
	   qcout_old = 0.0
	   jumpflag = 0 
	   loopstart = 1
	   loopflag = 0
	   nstep = 100000	

c        Phil adding TS control for MDA G fixes
c        reset time parameters at new realization

	   days = 0.0
	   day  = 1.0
	   daymin = 1.e-4
c        ----------------------------------

      end if

c-----------------------------------------------------------------------------------
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - METHOD=1 MD=1
       if((method.EQ.1).AND.(ing(1).GT.0)) then

c      - - - - - - - - - - - - - - - - - -
c        load values of IN for Time GT 0, 
c         set goldtime and 
c         daymax = goldtime OR 10 yrs whichever is smaller

	   do j = 1,size_of_in
           in(j)=ing(j)
	   end do

         goldtime = in(1)*365.25
	   if(in(5).GT.10) then
	     daymax = 365.25*10
	   else  
	     daymax = in(5)*365.25
	   end if

c      - - - - - - - - - - -  - - - - - 

c      WRITE OUT to a_in_array the values of tracer coming in from GoldSim

         write(iaunit,*) '**********************************'
         write(iaunit,*) '  Species coming in from GoldSim  '
         write(iaunit,*) ' FEHM#   GS#   in_location    Mass    '

         do j = 1, int(in(7))
	     k =  size_of_in - int(in(8+in(7)*2)) + int(in(7+j))
           write(iaunit,737)  j, int(in(7+j)) , k,  in(k)
	   end do


 737     format(3I7,G12.3)
c
         ex = .FALSE.
         if(maxmix_flag.ne.0) then
            call generate_flow
         end if

c     stop simulation after stress calc for certain stress input
         if(istrs.ne.0.and.istrs_coupl.eq.0) go to 170
c     set up time-spaced coupling         
         if(istrs_coupl.eq.-4) then
           timestress0 = days
           timestress = timestress0 + daystress
         endif
         if(istrs.ne.0.and.istrs_coupl.eq.0) go to 170
c     Before the time step loop create the partitions for zones

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

c       -------------------------------
c       Maintain loop count when coming in and out of GoldSim  PHS 5/11/2011  

	   if(jumpflag.EQ.666) then
	      loopstart = loopflag
         else 
	      loopstart = 1
	   end if 
c       -------------------------------

c ************** major time step loop ***************************
         do l = loopstart, nstep
c
c     counter that keeps accumulating when rip is the time step driver
c     Or in a conventional simulation with heat and mass solution
c
            tscounter = tscounter + 1
	      

c**** time step control via iteration count ****

 100        continue
            
c
c     Set current index for flow field catalog number (rip option)
c
            flowflag = int(in(2))

c*** adjust timestep size
            call timcrl

c          -----------------------------------
c          PHS 5/11/2011  adding to force final time to = GoldSim timestep
		  if(days.GT.goldtime) then
		    day = day - (days-goldtime) 
	        days = goldtime
	        if(day.LT.1.) day = 1.0
		    daymax = day
	        write(iaunit,*) '---------------------------'
	        write(iaunit,*) 'day daymax years goldtime-yrs '
	        write(iaunit,778) day,daymax,days/365.25,goldtime/365.25
	        write(iaunit,*) '---------------------------'
	      end if
 778        format(4(F9.3,' -- '))
c          -----------------------------------

	
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
         
               if (mlz .ne. 0)  then
                  irestart_ts = 1
                  if (mlz .eq. -1)  then
                     if (iout .ne. 0) write(iout, 6020)  l, iad, day
                     if (iptty .gt. 0)  write(iptty, 6020)  l, iad, day
 6020                format(/, 1x, 'timestep = ', i6, ' iterations = ', 
     *                    i4, ' timestep size = ', g15.4)
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
                  call diagnostics(1)     
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
               do i=1,n
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
                  if(istrs_coupl.eq.-3) then
                     istresscall = 1
                     call stress_uncoupled(3)
c update volume strains
                     call stressctr(6,0)
c update porosity
                     call stressctr(-7,0)
                  endif
c add displacements to total displacements	        
                  call stressctr(10,0)
c update displacements
                  call stressctr(12,0)
c update volume strains
                  call stressctr(-6,0)
c calculate stresses
                  call stressctr(13,0)	
c allocate memory for permeability update if necessay
                  call stress_perm(-1,0)
c update permeabilities (explicit)
                  call stress_perm(1,0)			            
c deallocate memory for permeability update if necessay
                  call stress_perm(-2,0)               	        
               endif 
c calculate subsidence
               call subsidence(1)         
c update peaceman term for wellbore pressure
               if(isubwd.ne.0)call wellimped_ctr(1)
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

c**** calculate velocities ****
            if(compute_flow .or. iccen .eq. 1) call veloc

c**** obtain concentration solution ****
            in3save = in(3)
            if(in(3).eq.0) then
               in(3) = irun + .0001
            end if
         
            call concen (1,tscounter)
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
               call wrtout(tassem,tas,dabs(tinfl),dabs(teinfl),
     &              dabs(inflow_thstime),dabs(inen_thstime),
     &              is_ch,is_ch_t)
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
!               else
!                  call contr_days (1)
!                  call contr_days (-1)
!               endif
               if (allocated(itc) .and. nicg .gt. 1) then
                  if (itc(nicg-1).gt.0) then
!                 call disk (1)
                     if (isave .ne. 0) call diskwrite
                     if (iflxn .eq. 1) call flxn
                  endif
               else
!              call disk (1)
                  if (isave .ne. 0) call diskwrite
                  if (iflxn .eq. 1) call flxn
               end if
               call pest(1)
               if (ifinsh .ne. 0)  ex = .TRUE.
               if (dayscn .ge. abs(contim)) dayscn = 0.

            end if


            if (ifinsh .eq. 1) then
c standard simulation finish
               istea_pest = 0
	         
               go  to  170
            else if (ifinsh .eq. 2) then

c Make sure last steady state time step is output
               if (hist_flag) then
                  call plot_new (1, dabs(tinfl), dabs(teinfl),
     &                 dabs(inflow_thstime), dabs(qtoti))
               end if
               call contr (1)
               call contr (-1)
               call river_ctr(6)
c finished the steady state simulation, now doing transient
               days = 0.0
               istea_pest = 0
               call pest(1)
               istea_pest = 1
               call flow_boundary_conditions(4)
               go to 999
            endif

            write(iaunit,773) l,day,days,goldtime     		                       
 773        format('step tstep time Goldtime  ',i8,3G10.3)

c- - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - Jump Back to GS
	      if(days.GE.goldtime) then
	        call loadoutarray_trac
	        loopflag = l+1
	        jumpflag = 666
              goto 9999
	      end if   
c- - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - Jump Back to GS

         end do    ! l loop  

c ******************* end major time step loop ****************

c**** write solution to plot tapes ****

         l = l - 1

 170     continue
cSPC
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
        
         if (iout .ne. 0) write(iout, 6040)  days, l
         if (iptty .gt. 0)  write(iptty, 6040)  days, l
 6040    format(//, 1x, 'simulation ended: days ', 1pg30.23, 
     *        ' timesteps ', i5)
      
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
                        if (iflxn .eq. 1) call flxn
                     endif
                  else
!                 call disk (nsave)
                     if (isave .ne. 0) call diskwrite
                     if (iflxn .eq. 1) call flxn
                  end if
!               if(contim.ge.0) then
                  if(istrs_coupl.ne.-2.and.istrs_coupl.ne.-1) then
! contr will be called below in stress_uncoupled for a stress solution
                     call contr (1)
                     call contr (-1)
                  end if
                  call river_ctr(6)               
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
 	 
c RJP 04/10/07 this is for co2 properties table look-up
         call interpolation_arrays_deallocate()

       end if 
	 
	   
c     End if around the part of code for performing the simulation
 

c     Add statement to set return flag

c     For the RIP version, set so that no error
c     is passed back to rip

 9999   continue    
        state = 0

        if((method.EQ.1).AND.(iptty.GT.1000)) then 
          write(iptty,*) 'Return fehmn.f  Cluster Method Realiz ',
     x                   int(in(6)),method, int(in(3))
        end if


        do i = 1, 500
          inquire(unit=i,opened=it_is_open)
          if(it_is_open) then
             close(i)
          end if
        end do     
         
      return

c--------------------------------------------------------------------------------
c * * * * * * * * * * * * * * * * * * * *                   Final Return to Main
c * * * * * * * * * * * * * * * * * * * *             Below here lie subroutines
c--------------------------------------------------------------------------------
      contains

      subroutine loadoutarray_trac
cHari compute conc values to pass back to goldsim
      implicit none
      real*8, allocatable :: out_save(:)
      real*8, allocatable :: time_dump(:)
      integer ispecies
      integer number_of_species
      integer ns2,izones, sflag, node1, node393
      integer nflow_frac, number_of_zones,indexout,indexmzone
      integer add_spots, add_spots2
      real*8 cur_time, prev_time, del_time
      real*8 :: cur_time_save = 0.
      save out_save, time_dump, cur_time_save

c      index_in_species= 4 + 6 + 
c      number_of_species = int(in(index_in_species))
c      Changing SUPER INDEX CHANGE  9/28/04 in(8+in(7)*2))
cSPC note: for not using index 66 for Kd, here should be in(8+in(7))

      write(iaunit,*) 'loadoutarray - species loaded',in(8+in(7)*2)
      do jjj = 1, in(8+in(7)*2)
         out(jjj) = 0.0
      end do

      write(iaunit,*) '-------------------------------'
      write(iaunit,*) 'sflag goldsp mass_out'
      do sflag = 1, in(7)
         out(in(7+sflag)) = qcout(sflag) - qcout_old(sflag)
         qcout_old(sflag) = qcout(sflag)
         write(iaunit,*) sflag,int(in(7+sflag)),out(in(7+sflag))
  775    format(2I6,E12.4)
  776    format(I6,2E12.4)
      end do
	write(iaunit,*) '-------------------------------'
	write(iaunit,*) '----   YEARS   ', days/365.25
	write(iaunit,*) ' node1 conc    node393 conc  '
	do sflag = 1, int(in(7))
	  node1 = (sflag - 1 )* neq + 1
	  node393 = (sflag-1)*neq + 393
	  write(iaunit,776) sflag, anl(node1) , anl(node393)
	end do
      write(iaunit,*) '-------------------------------'

      return
      end subroutine loadoutarray_trac

c - - - - - Subroutine namefiles  Names the .files file
c- - - - - - SPC for multi cluster run

      subroutine namefiles
      implicit none
      
      integer ncase, open_file
	character*1  infil(10)
	character*2 species(19), cluster(11) 
      character*74 zzout
	

        infil(1) = 'A'
        infil(2) = 'B'
        infil(3) = 'C'
        infil(4) = 'D'
        infil(5) = 'E'
        infil(6) = 'F'
        infil(7) = 'G'
        infil(8) = 'H'
        infil(9) = 'I'
        infil(10) = 'J'

        cluster(1) = '01'
        cluster(2) = '02'
        cluster(3) = '03'
        cluster(4) = '04'
        cluster(5) = '05'
        cluster(6) = '06'
        cluster(7) = '07'
        cluster(8) = '08'
        cluster(9) = '09'
        cluster(10) = '10'
        cluster(11) = '11'

	  species(1)  = '01'
	  species(19) = '19'

        zzout(1:29) = 'clusterXX/A/fehmn_XX_YY.files '
        zzout(8:9) = cluster(int(in(6)))
	  zzout(11:11) = infil(int(in(2)))
	  zzout(19:20) = cluster(int(in(6)))
	  zzout(22:23) = species(int(in(7)))

        nmfil(1) = zzout

      return
      end subroutine namefiles

      end subroutine fehmn

