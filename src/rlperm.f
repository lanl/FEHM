      subroutine rlperm(ndummy,iz)
!***********************************************************************
!  Copyright, 2004,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
C**********************************************************************
CD1
CD1 PURPOSE
CD1
CD1 To compute the relative permeabilities for vapor and liquid.
CD1
C**********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Revision                    ECD
CD2 Date         Programmer     Number  Comments
CD2
CD2 3-4-94       G. Zyvoloski   N/A     Initial implementation
CD2                                     
CD2 $Log:   /pvcs.config/fehm90/src/rlperm.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:43:06   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:14:16   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:12:20   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:26:24   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:07:16   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:45:08 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.22   Fri Nov 21 13:26:00 1997   gaz
CD2 fine tuned some phase change parameters
CD2 
CD2    Rev 1.20   Mon Mar 31 08:42:10 1997   gaz
CD2 added new models
CD2 
CD2    Rev 1.19   Thu Jun 27 13:34:34 1996   zvd
CD2 Added error write to ierr
CD2 
CD2    Rev 1.18   Wed Jun 12 15:21:20 1996   zvd
CD2 Modified optional input file routines and reduced number of calls
CD2 
CD2    Rev 1.17   Mon Jun 03 11:18:30 1996   hend
CD2 Added macro name & comment capabi. to new input
CD2 
CD2    Rev 1.16   Fri May 31 11:04:50 1996   hend
CD2 Added optional input from specified file
CD2 
CD2    Rev 1.15   Fri Apr 26 16:00:02 1996   gaz
CD2 correction when upstream perm is used, another space for rlperm model
CD2 
CD2    Rev 1.14   Fri Feb 16 11:15:08 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.13   Wed Feb 07 11:28:56 1996   gaz
CD2 changed calling sequence for linear and corey functions
CD2 
CD2    Rev 1.12   Thu Feb 01 15:56:08 1996   hend
CD2 Updated Requirements Traceability
CD2 
CD2    Rev 1.11   12/13/95 08:45:24   gaz
CD2 now capillary pressure input for corey and linear rlp models
CD2 
CD2    Rev 1.10   11/15/95 16:25:08   gaz
CD2 fixes so upwind perms except fracture matrix interaction
CD2 
CD2    Rev 1.9   08/18/95 10:33:56   llt
CD2 irpd already defined, removed for cray
CD2 
CD2    Rev 1.8   04/27/95 15:25:34   robinson
CD2 Made sure 0.99 values were 0.99d0 in subroutine calls
CD2 
CD2    Rev 1.7   04/25/95 09:06:42   llt
CD2 retrieved lost log history information
CD2 
CD2    Rev 1.6   03/23/95 17:14:42   gaz
CD2 gaz : added changes for different smoothing techniques
CD2 both dry and wet end. Also put parameter iupk in common(It is now a macro)
CD2 
CD2    Rev 1.5   01/28/95 13:55:28   llt
CD2 water balance equation was modified
CD2 
CD2    Rev 1.4   08/23/94 08:23:42   llt
CD2 gaz changes
CD2
CD2    Rev 1.3   03/28/94 16:36:52   robinson
CD2 Fixed memory problem.
CD2
CD2    Rev 1.2   03/23/94 14:45:12   robinson
CD2 Revamped van Genutchen rel. perm. implementation
CD2
CD2    Rev 1.1   03/21/94 08:33:56   gaz
CD2 Added solve_new and cleaned up memory management.
CD2
CD2    Rev 1.0   01/20/94 10:27:16   pvcs
CD2 original version in process of being certified
CD2
c
c gaz 3/6/95 added calls to vgcap for linear models
C**********************************************************************
CD3
CD3 INTERFACES
CD3
CD3 Formal Calling Parameters
CD3 
CD3 Identifier   Type    Use     Description
CD3 
CD3 iz           int     I       Control parameter
CD3 ndummy       int     I       Parameter to set the correct node
CD3                                  number for dual porosity
CD3 
CD3 Interface Tables
CD3
CD3 None
CD3
CD3 Files
CD3
CD3 Name                  Use   Description
CD3 
CD3 inpt                   I    File contains input data
CD3
C**********************************************************************
CD4
CD4 GLOBAL OBJECTS
CD4
CD4 Global Constants
CD4
CD4 Identifier  Type     Description
CD4 
CD4
CD4 Global Types
CD4
CD4 NONE
CD4
CD4 Global Variables
CD4 
CD4 ipsx, ipdeef, ipdepf, ippdmef, ipdq, ipdqh, wdd1, inpt, irlpt, rp1f,
CD4 rp2f, rp3f, rp4f, rp5f, rp6f, rp7f, cp1f, rp18f, rp19f, rp20f, cp3f,
CD4 irlp, rp8f, rp9f, rp10f, rp11f, rp12f, rp13f, rp15f, rp16f, rp17f,
CD4 rp14f, cp2f, rp21f, rp22f, rp23f, cp4f, neq, pcp, dpcef, idpdp,
CD4 idualp, pnx, pny, pnz, drlef, drlpf, drvpf, s, narrays, pointer,
CD4 itype, default, igroup, ireturn,
CD4 macroread
CD4 
CD4                      COMMON
CD4 Identifier   Type    Block   Description
CD4
CD4  
CD4 Global Subprograms
CD4
CD4 Identifier    Type     Description
CD4 
CD4 null1         LOGICAL  Determines if there is additional data to
CD4                            read
CD4 welbor        N/A      Computes relative permeabilities for
CD4                            wellbore option
CD4 initdata      N/A      Reads and sets input values defined at each
CD4                              node
CD4 
CD4
C**********************************************************************
CD5
CD5 LOCAL IDENTIFIERS
CD5
CD5 Local Constants
CD5
CD5 Identifier   Type        Description
CD5 
CD5 scutm        real*8      Minimum saturation cutoff
CD5 iupk         int         Flag indicating upstream permeability
CD5                             weighting
CD5 hmin         real*8      Term used in van Genutchen calculation
CD5 darcyf       real*8      Term used in van Genutchen calculation
CD5
CD5 Local Types
CD5
CD5 None
CD5
CD5 Local variables
CD5
CD5 Identifier   Type        Description
CD5 
CD5 i            int         Counter for number of values read in
CD5 irlpd        int         Input flag denoting which relative
CD5                              permeability model is used
CD5 alpha        real*8      Term in van Genutchen calculation
CD5 beta         real*8      Term in van Genutchen calculation
CD5 alamda       real*8      Exponent in van Genutchen calculation
CD5 alpi         real*8      Term in van Genutchen calculation
CD5 smcut        real*8      Cutoff saturation in van Genutchen
CD5                              calculation(normalized)
CD5 slcut        real*8      Cutoff saturation in van Genutchen
CD5                              calculation
CD5 fac          real*8      Term in van Genutchen calculation
CD5 ds           real*8      Term in van Genutchen calculation
CD5 dhp          real*8      Term in van Genutchen calculation
CD5 izonef       int         Array of zone numbers for each node
CD5 izone        int         Current zone number for this node
CD5 ja           int         Do loop index for setting parameter values
CD5 jb           int         Do loop index for setting parameter values
CD5 jc           int         Do loop index for setting parameter values
CD5 ij           int         Do loop index over all nodes
CD5 mi           int         Current node number including dual
CD5                              porosity nodes
CD5 ieosd        int         Current equation of state
CD5 it           int         Current relative permeability model
CD5 irpd         int         Current region number
CD5 rl           real*8      Liquid relative permeability
CD5 rv           real*8      Vapor relative permeability
CD5 drls         real*8      Derivative of liquid relative
CD5                              permeability with saturation
CD5 drvs         real*8      Derivative of vapor relative permeability
CD5                              with saturation
CD5 rp1          real*8      Term used in van Genutchen calculation
CD5 rp2          real*8      Term used in van Genutchen calculation
CD5 rp3          real*8      Term used in van Genutchen calculation
CD5 rp4          real*8      Term used in van Genutchen calculation
CD5 denom        real*8      Denominator used in van Genutchen
CD5                              calculation
CD5 star         real*8      Dimensionless saturation used in van
CD5                              genutchen calculation
CD5 sl           real*8      Liquid saturation
CD5 hp           real*8      Term used in van Genutchen calculation
CD5 rl1          real*8      Liquid relative permeability
CD5 rv1          real*8      Liquid relative permeability
CD5 drls1        real*8      Derivative of liquid relative
CD5                              permeability with saturation
CD5 drvs1        real*8      Derivative of vapor relative
CD5                              permeability with saturation
CD5 akf          real*8      Term used in van Ganutchen calculation
CD5 akm          real*8      Term used in van Ganutchen calculation
CD5 porf         real*8      Term used in van Ganutchen calculation
CD5 permb        real*8      Permeability
CD5 
CD5 Local Subprograms
CD5
CD5 None
CD5
C**********************************************************************
CD6
CD6 FUNCTIONAL DESCRIPTION
CD6 
CD6 
C**********************************************************************
CD7
CD7 ASSUMPTIONS AND LIMITATIONS
CD7
CD7
C**********************************************************************
CD8
CD8 SPECIAL COMMENTS
CD8
CD8  Requirements from SDN: 10086-RD-2.20-00
CD8    SOFTWARE REQUIREMENTS DOCUMENT (RD) for the
CD8    FEHM Application Version 2.20
CD8
C**********************************************************************
CD9
CD9 REQUIREMENTS TRACEABILITY
CD9
CD9 2.4.4 Relative-permeability and capillary-pressure functions
CD9 2.6   Provide Input/Output Data Files
CD9 3.0   INPUT AND OUTPUT REQUIREMENTS
CD9
C**********************************************************************
CDA
CDA REFERENCES
CDA
CDA See FEHMN SRS, MMS, and SDD, Robinson's memo EES-4-92-354 for
CDA documentation.
CDA
C**********************************************************************
CPS
CPS PSEUDOCODE
CPS 
CPS BEGIN rlperm
CPS 
CPS IF data are read in on this call
CPS 
CPS   LOOP to read first group of data
CPS     read input
CPS   EXIT IF no more data (null1)
CPS     IF there is more data to read (not null1)
CPS       Read relative permeability model type
CPS       
CPS       IF the model is linear
CPS          Reread model type and 4 parameters
CPS       ELSE IF the model is Corey
CPS          Reread model type and 2 parameters
CPS       ELSE IF the model is van Genutchen
CPS          Reread model type and 6 parameters
CPS       ELSE IF the model is fracture/matrix van Genutchen
CPS          Reread model type and 15 parameters
CPS       END IF
CPS       
CPS       IF the model is either van Genutchen model
CPS         Calculate cutoff saturation and head difference
CPS         Calculate cefficients for spline fits
CPS         Calculate upper cutoff
CPS         Calculate additional cefficients for spline fits
CPS       ENDIF
CPS       
CPS       IF the model is fracture/matrix van Genutchen
CPS         Calculate cutoff saturation and head difference
CPS         Calculate coefficients for spline fits
CPS         Calculate upper cutoff
CPS         Calculate additional cefficients for spline fits
CPS       ENDIF
CPS       
CPS     ENDIF
CPS   
CPS   ENDLOOP
CPS   
CPS   Set parameters for reading region numbers
CPS     
CPS   initdata - read region numbers and set at each node
CPS   
CPS   Set flag to denote that the rlp macro has been called
CPS       
CPS ELSE we are calculating relative permeabilities on this call
CPS 
CPS   FOR each node
CPS   
CPS     IF this is a liquid only system
CPS       Set relative permeabilities and derivatives accordingly
CPS     ELSEIF this is a superheated vapor
CPS       Set relative permeabilities and derivatives accordingly
CPS     ELSEIF this is a two-phase system
CPS       
CPS       IF a linear model is used
CPS         Set parameters
CPS         IF liquid saturation is greater than cutoff
CPS           Set relative permeability to 1
CPS           IF liquid saturation is greater than other cutoff
CPS             Compute liquid relative permeability and derivative
CPS           ENDIF
CPS         ENDIF
CPS         IF vapor saturation is greater than cutoff
CPS           Set relative permeability to 1
CPS           IF vapor saturation is greater than other cutoff
CPS             Compute vapor relative permeability and derivative
CPS           ENDIF
CPS         ENDIF
CPS         
CPS       ELSEIF a Corey model is used
CPS       
CPS         Compute parameters
CPS         IF model results in only vapor flow
CPS           Set vapor relative permeability to 1
CPS         ELSEIF model results in only liquid flow
CPS           Set liquid relative permeability to 1
CPS         ELSE both vapor and liquid are moving
CPS           Compute liquid and vapor relative permeabilities and...
CPS           ... derivatives
CPS         ENDIF
CPS         
CPS       ELSEIF a van Genutchen single porosity relation is used
CPS       
CPS         vgcap - compute the capillary pressure
CPS         vgrlp - compute the relative permeability
CPS         Convert values for capillary pressure
CPS
CPS       ELSEIF the van Genutchen fracture model parameter are invoked
CPS         
CPS         Set permeability parameter values
CPS
CPS         IF this is a dpdp or dual porosity solution
CPS           IF this is a fracture node
CPS             vgcap - compute fracture capillary pressure
CPS             vgrlp - compute fracture relative permeability
CPS             Compute fracture permeability
CPS           ELSE this is a matrix node
CPS             vgcap - compute matrix capillary pressure
CPS             vgrlp - compute matrix relative permeability
CPS             Compute matrix permeability
CPS           ENDIF
CPS         ELSE it is a continuum solution
CPS           vgcap - compute matrix capillary pressure
CPS           vgrlp - compute matrix relative permeability
CPS           Compute normalized saturation
CPS           vgrlp - compute fracture relative permeability
CPS           Set fracture porosity weighted relative permeabilities
CPS         ENDIF
CPS           
CPS         IF upstream weighting of permeability is used
CPS           Correct relative permeability values
CPS         ENDIF
CPS           
CPS         Set permeabilities
CPS         Convert capillary pressure
CPS           
CPS       ENDIF
CPS       
CPS     ENDIF ends the various fluid types
CPS     
CPS     IF relative permeability models are employed
CPS       Set values for relative permeabilities and derivatives in...
CPS       ... arrays
CPS     ENDIF
CPS   
CPS   ENDFOR
CPS 
CPS   welbor - compute permeabilities for wellbore model option
CPS 
CPS ENDIF data are read in on this call
CPS 
CPS END rlperm
CPS
C**********************************************************************
c
c calculates relative permiabilities and derivatives
c

      use comhi
      use comfi
      use comei
      use comdi
      use comci
      use combi
      use comdti
      use comai
      use comki
      implicit none

      real*8 scutm,hmin,darcyf,tol_l,tol_u,su_cut 
      parameter(scutm = 1.d-03)
      parameter(hmin = 1.d-8)
      parameter(darcyf = 1.d12)
      parameter(tol_l  = 1.d-5)
      parameter(tol_u  = 1.d-5)
      parameter(su_cut = 0.7d00)


      integer iz,ndummy,i,irlpd,mi,ieosd,it,ir,j,num_models
      integer :: ireg = 1
      real*8 alpha,beta,alamda,alpi,smcut,slcut,fac,ds,dhp
      real*8 rp1,rp2,rp3,rp4,denom,star,hp,rl1,rv1,hp1,dhp1
      real(8) :: rl = 1., rv = 1., drls = 0., drvs = 0.
      real*8 drls1,drvs1,akf,akm,porf,permb,sl
      real*8 smcutm,smcutf,alpham,alamdam,facf
      real*8 rpa1, rpa2, rpa3, rpa4, rpa5
      real*8, allocatable :: xfptmp(:), yfptmp(:), zfptmp(:)
      logical null1,ex
      integer ireg0
      save ireg,ireg0
      if(l.eq.0) ireg0 = 0
      

      if(iz.eq.0) then
c     read in data

c     check for read from other file

         icapp = 0
         i = 0
         j = 0
         ex = .false.
         do
            read(inpt,'(a80)') wdd1
            if(null1(wdd1)) exit
            backspace inpt
            read(inpt,*) irlpd
            backspace inpt
            i = i+1
            
            if (irlpd .eq. 1) then
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp3f(i)
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               icapp = 1
               icapt(i)=1
            else if (irlpd .eq. 12 .or. irlpd .eq. 13) then
c linear with residual saturations as function of hydrate fraction
c for use with methane hydrate 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp5f is linear term for irreducible water saturation
c rp6f is linear term for irreducible gas saturation
c model 12 assumes the rv does not have a frac_hyd term
c model 13 assumes the rv has a frac_hyd term
c RJP 6/1/04 introduced two new terms to give an option 
c to specify sakamoto_a1 and sakamoto_b1 terms through 
c specifying rp7f(i) and rp8f(i) respectively.
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp3f(i),rp5f(i),rp6f(i),rp7f(i),
     &              rp8f(i)
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0. 
               endif
               icapp = 1
               icapt(i)=1
            else if (irlpd .eq. 14) then
c linear with residual saturations as function of hydrate fraction
c for use with methane hydrate 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp3f is maximum liquid rl perm
c rp4f is maximum vapor rv perm
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp3f(i)
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               icapp = 1
               icapt(i)=1
            else if (irlpd .eq. 15) then
c linear with residual saturations as function of hydrate fraction
c for use with methane hydrate 
c rp1f is constant term for irreducible water saturation
c rp2f is constant term for irreducible gas saturation
c rp3f is maximum liquid rl perm
c rp4f is maximum vapor rv perm
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp3f(i)
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               icapp = 1
               icapt(i)=1
c RJP 04/19/07 Introduced new rel-perm model for CO2.  
            elseif (irlpd .eq. 16 .or. irlpd .eq. 17) then
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),rp5f(i),rp6f(i),rp7f(i), rp8f(i),
     &              rp9f(i),rp10f(i),rp11f(i),rp12f(i), rp13f(i),
     &              rp14f(i)
            elseif (irlpd .eq. 18) then
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),rp5f(i),rp6f(i),rp7f(i), rp8f(i),
     &              rp9f(i),rp10f(i),rp11f(i),rp12f(i), rp13f(i),
     &              rp14f(i),rp15f(i),rp16f(i)
            else if (irlpd .eq. 21) then
c Thomeer capillary pressure and Corey relative perms
c rp1f is irreducible water saturation
c rp2f is irreducible gas saturation
c rp3f is pd1 in Boitnott's Thomeer formulation
c rp4f is fg1 in Boitnott's Thomeer formulation
c cp1f is minimum saturation cutoff for Thomeer formulation
c cp2f is maximum saturation cutoff for Thomeer formulation
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp2f(i)
               
            else if (irlpd .eq. 10) then
c linear with mim and max l and v saturations
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),rp5f(i),rp6f(i),cp1f(i),cp3f(i)
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               icapp = 1
               icapt(i)=1
            else if (irlpd .eq. 11) then
c Brooks-Corey capillary and rlperms
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),rp3f(i),
     &              rp4f(i),cp1f(i),cp3f(i)
            else if (irlpd .eq. -1) then
c constant rlperms
               read(inpt,*) irlpt(i),rp1f(i),rp2f(i),
     &              cp1f(i),cp3f(i)
c     If the permeability was entered as a negative value, use log permeability
               if (rp1f(i) .lt. 0.d0) 
     &              rp1f(i) =  1.0d00/10.0d00**(abs(rp1f(i)))
               if (rp2f(i) .lt. 0.d0) 
     &              rp2f(i) =  1.0d00/10.0d00**(abs(rp2f(i)))
               cp2f(i)=0.0
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               icapp = 1
               icapt(i)=1
            else if (irlpd .eq. 2) then
c Corey(geothermal) rlperms and linear capillary pressure
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i),
     *              cp1f(i), cp3f(i)
               cp2f(i)=0.0
               icapp = 1
               icapt(i)=1
               if((cp3f(i)-cp2f(i)) .gt. 0. ) then
                  cp1f(i) = cp1f(i) / (cp3f(i) - cp2f(i))
               else
                  cp1f(i) = 0.
               endif
               if(cp1f(i).gt.0.0) then
                  icapp = 1
                  icapt(i)=1
               endif
            else if (irlpd .eq. 3) then
c     van Genutchen capillary pressure and relative perms
c with rlp(h) 
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i), 
     *              rp4f(i), rp5f(i), rp6f(i)
            else if (irlpd .eq. -4) then
c     van Genutchen capillary pressure and relative perms fracture/matrix
c with rlp(h) non-equilibrium cap pressures
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i),
     *              rp4f(i), rp5f(i), rp6f(i), rp11f(i), rp12f(i), 
     *              rp13f(i), rp14f(i), rp15f(i), rp16f(i), rp17f(i),  
     *              rp18f(i), rp19f(i)
            else if (irlpd .eq. 4) then
c     van Genutchen capillary pressure and relative perms fracture/matrix
c with rlp(h) 
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i),
     *              rp4f(i), rp5f(i), rp6f(i), rp11f(i), rp12f(i), 
     *              rp13f(i), rp14f(i), rp15f(i), rp16f(i), rp17f(i),  
     *              rp18f(i), rp19f(i)      
            else if (irlpd .eq. 5) then
c     van Genutchen capillary pressure and relative perms
c with rlp(S) 
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i), 
     *              rp4f(i), rp5f(i), rp6f(i)
            else if (irlpd .eq. 6) then
c     van Genutchen capillary pressure and relative perms fracture/matrix
c with rlp(S) 
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i),
     *              rp4f(i), rp5f(i), rp6f(i), rp11f(i), rp12f(i), 
     *              rp13f(i), rp14f(i), rp15f(i), rp16f(i), rp17f(i),  
     *              rp18f(i), rp19f(i)
            else if (irlpd .eq. 7) then
c     van Genutchen capillary pressure and relative perms fracture/matrix
c with rlp(S)  and special fracture interaction term
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i),
     *              rp4f(i), rp5f(i), rp6f(i), rp11f(i), rp12f(i), 
     *              rp13f(i), rp14f(i), rp15f(i), rp16f(i), rp17f(i),  
     *              rp18f(i), rp19f(i), rp24f(i)
c allocate memory for f-m interaction and check for errors
               call rlp_frac(0,0,0.0d00,0.0d00,0.0d00,
     &              0.0d00,0.0d00,0.0d00,0.0d00,0.0d00,
     &              rp24f(i))
            else if (irlpd .eq. 8) then
c     van Genutchen capillary pressure and relative perms
c with rlp(h) 
c with relative perm of gas phase governed by VG parameters
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i), 
     &              rp4f(i), rp5f(i), rp6f(i)
            else if (irlpd .eq. 9) then
c     van Genutchen capillary pressure and relative perms
c with rlp(h) 
c with relative perm of gas phase governed by VG parameters
c with terms from Touma and Vauclin (as reported by Celia and Binney)
c rp13- Bw 
c rp14- Kas (multiplier of Ksat)
c rp15- Aa
c rp16- Ba
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i), 
     &              rp4f(i), rp5f(i), rp6f(i), 
     &              rp13f(i), rp14f(i), rp15f(i), rp16f(i) 
c keating 9/2008
            else if (irlpd .eq. 19) then
               read(inpt,*) irlpt(i), rp1f(i), rp2f(i), rp3f(i), 
     &              rp4f(i), rp5f(i), rp6f(i), rp7f(i)
            endif
c     If the permeability was entered as a negative value, use log permeability
            if (rp17f(i) .lt. 0.d0) 
     &           rp17f(i) =  1.0d00/10.0d00**(abs(rp17f(i)))
            if (rp18f(i) .lt. 0.d0) 
     &           rp18f(i) =  1.0d00/10.0d00**(abs(rp18f(i)))
                        
         end do
c     ****** end of input loop
         
 20      num_models=max(i,j-1)
         do i=1,num_models
            
            if(irlpt(i).ge.3 .and. irlpt(i)
     &        .le. 9. or. irlpt(i) .eq. -4 ) then
               
c     calculate cutoff saturation and head difference
c     GAZ 8-16-94 added slcut to work with liquid saturation
               alpha = rp3f(i)
               beta = rp4f(i)
               alamda = 1.0-1.0/beta
               smcut=(rp6f(i)-rp1f(i))/(rp2f(i)-rp1f(i))
               smcutm = max(smcut,scutm)
               slcut = smcutm*(rp2f(i)-rp1f(i)) + rp1f(i)
               fac=rp5f(i)
               if(fac.eq.0.0) then
                  call vgcap_fit(1,rp1f(i),rp2f(i),slcut,smcutm,fac,
     &                 alpha,alamda,rp7f(i),rp8f(i),rp9f(i),rp10f(i),
     &                 hmin)
               else if(fac.gt.0.0) then
                  call vgcap_fit(2,rp1f(i),rp2f(i),slcut,smcutm,fac,
     &                 alpha,alamda,rp7f(i),rp8f(i),rp9f(i),rp10f(i),
     &                 hmin)
               else if(fac.lt.0.0) then
                  call vgcap_fit(3,rp1f(i),rp2f(i),slcut,smcutm,fac,
     &                 alpha,alamda,rp7f(i),rp8f(i),rp9f(i),rp10f(i),
     &                 hmin)
               endif
               rp6f(i) = smcutm
c     get fit at saturated end(star=su_cut)
               slcut = su_cut*(rp2f(i)-rp1f(i)) + rp1f(i)
               call vgcap_fit(4,rp1f(i),rp2f(i),slcut,su_cut ,fac,alpha
     &              ,alamda, 0.0d0, 0.0d0, cp1f(i),cp2f(i),hmin)
            
               if (irlpt(i) .ne. 3 .and. irlpt(i) .ne. 5 
     &              .and. irlpt(i) .lt. 8 .or. irlpt(i) .eq. -4 ) then
               
C     terms for the fractures
c     calculate cutoff saturation and head difference
c     GAZ 8-16-94 added slcut to work with liquid saturation
                  alpham = alpha
                  alamdam = alamda
                  alpha = rp13f(i)
                  beta = rp14f(i)
                  alamda = 1.0-1.0/beta
                  smcut=(rp16f(i)-rp11f(i))/(rp12f(i)-rp11f(i))
                  smcutf = max(smcut,scutm)
                  slcut = smcutf*(rp12f(i)-rp11f(i)) + rp11f(i)
                  if(fac.ge.0.0d00) then
c next call matches cutoff capillary pressure
                     facf = rp15f(i)
                     call vgcap_match(2,smcutm,alpham,alamdam,fac
     &                    ,smcutf,alpha,alamda,facf)
                     fac=facf
                     rp15f(i)=facf
                  else
                     fac = rp15f(i)
                  endif
                  if(fac.gt.0.0) then
                     call vgcap_fit(2,rp11f(i),rp12f(i),slcut,smcutf,
     &                    fac,alpha,alamda,rp20f(i),rp21f(i),rp22f(i),
     &                    rp23f(i),hmin)
                  else if(fac.eq.0.0) then
                     call vgcap_fit(1,rp11f(i),rp12f(i),slcut,smcutf,
     &                 fac,alpha,alamda,rp20f(i),rp21f(i),rp22f(i),
     &                 rp23f(i),hmin)
                  else if(fac.lt.0.0) then
                     call vgcap_fit(3,rp11f(i),rp12f(i),slcut,smcutf,
     &                 fac,alpha,alamda,rp20f(i),rp21f(i),rp22f(i),
     &                 rp23f(i),hmin)
                  endif
                  rp16f(i) = smcutf
c     get fit at saturated end(star=su_cut)
                  slcut = su_cut*(rp12f(i)-rp11f(i)) + rp11f(i)
                  call vgcap_fit(4,rp1f(i),rp2f(i),slcut,su_cut,fac
     &                 ,alpha,alamda,0.0d0,0.0d0,cp3f(i),cp4f(i),hmin)
               endif
               if (irlpt(i) .ge. 5 .and. irlpt(i) .le. 7) then
                  rp1=rp1f(i)
                  rp2=rp2f(i)
                  rp3=rp3f(i)
                  rp4=rp4f(i)
                  call  vgrlps(0, sl, star, rp3, rp4, rp1, rp2, 
     2                 tol_l, tol_u, rl, drls, rv, drvs )
               endif
               if (irlpt(i) .eq. 9) then
c     calculate normalisation factor for saturation  
                  rp17f(i) = 1.d00/rp2f(i)**rp13f(i)
               endif
            end if
            if (irlpt(i) .eq.21) then
c    initialization of arrays for Thomeer 
               call thomeercap(0,0,0,0.,0.,0.,0.,0.,0.,0.,0.)
c    calculate cutoff values for cap pressure
               call thomeercap(1,0,i,0.d0,0.d0,
     &              0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
            endif 
            
         enddo                
c     
c     return if read data from a file
c     
         if(ex) return
         
c     read in nodal capillary type
         
c     read in nodal capillary type
      
         narrays = 1
         itype(1) = 4
         default(1) = 1
         macro = "rlp "
         igroup = 2
         call initdata2( inpt, ischk, n0, narrays,
     2        itype, default, macroread(7), macro, igroup, ireturn,
     3        i4_1=irlp(1:n0) )

         macroread(7) = .TRUE.
         
         do i=1,n0
            if(irlpt(irlp(i)).le.2) then
               icap(i)=irlp(i)
            else if(irlpt(irlp(i)).eq.21) then
               icap(i) = 1
            else
               icap(i)=0
            endif
         enddo
      
         if(idpdp.eq.0) then
            do mi=1,neq
               it = irlp(mi)
               ir = irlpt(it)
               if((ir.eq.4.and.rp19f(it).lt.0.0).or.
     &              (ir.eq.7.and.rp19f(it).lt.0.0)) then
                  write(ierr, 120)
                  write(ierr, 100) 
                  write(ierr, 110) 
                  write(ierr,*) "stopping"
                  write(ierr, 120)
                  if (iout .ne. 0) then
                     write(iout, 120)
                     write(iout, 100) 
                     write(iout, 110) 
                     write(iout,*) "stopping"
                     write(iout, 120)
                  end if
                  if (iptty .ne. 0) then
                     write(iptty, 120)
                     write(iptty, 100) 
                     write(iptty, 110)
                     write(iptty,*) "stopping"
                     write(iptty, 120)
                  end if
                  stop
               endif
            enddo
         endif

         if (.not. allocated(xfperm)) then
            allocate(xfperm(nrlp),yfperm(nrlp),zfperm(nrlp))
            xfperm = 1.d0
            yfperm = 1.d0
            zfperm = 1.d0
         else if (fperm_flag) then
c     Assign values from temporary fperm arrays
            allocate (xfptmp(nrlp), yfptmp(nrlp), zfptmp(nrlp))
            xfptmp = 1.d0
            yfptmp = 1.d0
            zfptmp = 1.d0
            do mi = 1, n0
               ir = irlp(mi)
               if (ir  .ne. 0) then
                  xfptmp(ir) = xfperm(mi)
                  yfptmp(ir) = yfperm(mi)
                  zfptmp(ir) = zfperm(mi)
               end if
            end do
            deallocate (xfperm, yfperm, zfperm)
            allocate(xfperm(nrlp),yfperm(nrlp),zfperm(nrlp))
            xfperm = xfptmp 
            yfperm = yfptmp
            zfperm =  zfptmp 
            deallocate (xfptmp, yfptmp, zfptmp)
         end if
            
 100     format ("cannot have anisotropic perms for rlp model 4")
 110     format ("or rlp model 7 with equivalent continuum")
 120     format ("*********************************************")
      else
c     
c     calculate relative perm and derivatives
c     
         do i = 1,neq
            ireg0 = ireg
            mi = i+ndummy
            ieosd = ieos(mi)
            if (rlp_flag .eq. 0) then
               it = 0
            else
               it = irlp(mi)
            end if
c     18-mar-94 gaz
c     if it=0 set irpd=0
            if(it.eq.0) then
               irpd=0
            else
               irpd = irlpt(it)
            endif
            if(irpd.ge.3.and.irpd.le.8) then 
               ieosd = 2
            endif
            if(ieosd.eq.1) then
c     liquid only
               rl = 1.0
               rv = 0.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.3) then
c     superheated vapor
               rl = 0.0
               rv = 1.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.4) then
c     no permeability, no porosity
               rl = 0.0
               rv = 0.0
               drls = 0.0
               drvs = 0.0
            elseif(ieosd.eq.2) then
               if (irpd .ne. 0) then
                  rp1 = rp1f(it)
                  rp2 = rp2f(it)
                  rp3 = rp3f(it)
                  rp4 = rp4f(it)
               else
! Use linear model as default
                  rp1 = 0.
                  rp2 = 0.
                  rp3 = 1.
                  rp4 = 1.
                  irpd = 1
               end if
               sl = s(mi)
               if(irpd.eq.1) then
c     two-phase mixture
c     linear function of saturations
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  sv = 1.0-sl
                  if(sl.ge.rp1) then
                     rl = 1.0
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom
                        drls = 1.0/denom
                     endif
                  endif
                  if(sv.ge.rp2) then
                     rv = 1.0
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom
                        drvs=-1.0/denom
                     endif
                  endif
               else if(irpd.eq.10) then
c     two-phase mixture
c     linear function of saturations
c     minimum rlperms
                  rl = rp5f(it)
                  rv = rp6f(it)
                  drls = 0.0
                  drvs = 0.0
                  sv = 1.0-sl
                  if(sl.ge.rp1) then
                     if(sl.le.rp3) then
                        denom = rp3-rp1
                        rl=(sl-rp1)/denom+rp5f(it)
                        drls = 1.0/denom
                     endif
                     if(sl.gt.rp3.or.rl.gt.1.0) then
                        rl=1.0
                        drls=0.0
                     endif
                  endif
                  if(sv.ge.rp2) then
                     if(sv.le.rp4) then
                        denom = rp4-rp2
                        rv=(sv-rp2)/denom+rp6f(it)
                        drvs=-1.0/denom
                     endif
                     if(sv.gt.rp4.or.rv.gt.1.0) then
                        rv=1.0
                        drvs=0.0
                     endif
                  endif
               elseif(irpd.eq.-1) then
c constant rl perms
                  rl = rp1
                  rv = rp2
                  drls = 0.0
                  drvs = 0.0
               elseif(irpd.eq.2) then
c     corey relationships(geothermal)
                  rl = 0.0
                  rv = 0.0
                  drls = 0.0
                  drvs = 0.0
                  denom = 1.0-rp1-rp2
                  star=(sl-rp1)/denom
                  ds = 1.0/denom
                  if(star.le.0.0) then
                     rv = 1.0
                  else if(star.ge.1.0) then
                     rl = 1.0
                  else
                     rl = star**4
                     drls = 4.0*star**3*ds
                     rv=(1.0-star)**2*(1.0-star**2)
                     drvs=(-2.0*(1.0-star)*(1.0-star**2)+
     2                    (1.0-star)**2*(-2.0*star))*ds
                  endif
               elseif(irpd.eq.3) then
c rlp(h)
c                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rl, drls, rv, drvs,0)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.8) then
c rlp(h) and vapor rlp
c                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rl, drls, rv, drvs,1)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.9) then
c rlp(h) and vapor rlp
c                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it),su_cut ,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  rpa1 = rp13f(it)
                  rpa2 = rp14f(it)
                  rpa3 = rp15f(it)
                  rpa4 = rp16f(it)
                  rpa5 = rp17f(it)
                  call vgrlpa(sl, star, rp3, rp4, hmin, hp,
     2                 dhp, rpa1, rpa2, rpa3, rpa4, rpa5,
     3                 rp1f(it), rp2f(it), rl, drls, rv, drvs)
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.5) then
c rlp(S)
c                  call vg_regions(1,ireg,mi,su_cut)
                  call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                 rp7f(it), rp8f(it), rp9f(it),
     3                 rp10f(it), rp6f(it), su_cut,
     4                 cp1f(it),cp2f(it),hp, dhp ,ireg    )
                  star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                  call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                 tol_l, tol_u, rl, drls, rv, drvs )
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.4) then
c     
c     akf-fracture saturated permeability(rp15)
c     akm-matrix saturated permeability(rp16)
c     porf-fracture fraction(rp17)
c     
                  akf = rp17f(it)
                  akm = rp18f(it)
                  porf = rp19f(it)
                  
                  if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                     
                     if( mi .le. neq ) then
c                        call vg_regions(2,ireg,mi,su_cut)
                        call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp, dhp, ireg    )
                        star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                       hp, dhp, rl, drls, rv, drvs,0)
                        permb = akf * porf
                        
                     else
c                        call vg_regions(1,ireg,mi,su_cut)
                        call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                       rp7f(it), rp8f(it), rp9f(it),
     3                       rp10f(it), rp6f(it), su_cut,      
     4                       cp1f(it),cp2f(it),hp, dhp, ireg    )
                        star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                        call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                       dhp, rl, drls, rv, drvs,0)
                        permb = akm * (1. - porf)
                     end if
                     
                  else
c gaz changed so that frac cap pressure is used with frac rel perm                  
c                     call vg_regions(1,ireg,mi,su_cut)
                     call vgcap( sl, rp1, rp2, rp3, rp4,  
     2                    rp7f(it), rp8f(it), rp9f(it),
     3                    rp10f(it), rp6f(it), su_cut,      
     4                    cp1f(it),cp2f(it),hp, dhp, ireg    )
                     star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                     call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                    dhp, rl, drls, rv, drvs,0)

                     star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                     call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                    hp, dhp, rl1, drls1, rv1, drvs1,0)
                     permb = akf * porf + akm * (1. - porf)
                     rl=(akf*rl1*porf+akm*rl*(1.0-porf))/permb
                     drls=(akf*drls1*porf+akm*drls*(1.0-porf))/permb
                     rv = 1. - rl
                     drvs = -drls
                  endif
                  if(iupk.gt.0) then
c     upstream weighting of intrinsic permeability
                     permb = permb*darcyf
                     rl = rl*permb
                     rv = rv*permb
                     drls = drls*permb
                     drvs = drvs*permb
                     permb = 1.0/darcyf
c     set saturated permeabilities(assume isotropic)
                     pnx(mi) = permb*1.e+6
                     pnz(mi) = pnx(mi) * zfperm(it)
                     pny(mi) = pnx(mi) * yfperm(it)
                     pnx(mi) = pnx(mi) * xfperm(it)
                  else
c     set saturated permeabilities(assume isotropic)
                     if(porf.ge.0.0) then
                        pnx(mi) = permb*1.e+6
                        pny(mi) = pnx(mi) * yfperm(it)
                        pnz(mi) = pnx(mi) * zfperm(it)
                        pnx(mi) = pnx(mi) * xfperm(it)
                     endif
                  endif
c     set capillary pressures
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd.eq.-4) then
c     
c     akf-fracture saturated permeability(rp15)
c     akm-matrix saturated permeability(rp16)
c     porf-fracture fraction(rp17)
c     
                  akf = rp17f(it)
                  akm = rp18f(it)
                  porf = rp19f(it)
                  
                  if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                     
                     if( mi .le. neq ) then
                        call vg_regions(2,ireg,mi,su_cut)
                        call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp, dhp, ireg    )
                        star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                       hp, dhp, rl, drls, rv, drvs,0)
                        permb = akf * porf
                        
                     else
                        call vg_regions(1,ireg,mi,su_cut)
                        call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                       rp7f(it), rp8f(it), rp9f(it), 
     3                       rp10f(it), rp6f(it), su_cut,      
     4                       cp1f(it),cp2f(it),hp, dhp, ireg    )
                        star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                        call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                       dhp, rl, drls, rv, drvs,0)
                        permb = akm * (1. - porf)
                     end if
                     
                  else
c gaz changed so that frac cap pressure is used with frac rel perm                  
                     call vgcap( sl, rp1, rp2, rp3, rp4,  
     2                    rp7f(it), rp8f(it), rp9f(it),
     3                    rp10f(it), rp6f(it), su_cut,      
     4                    cp1f(it),cp2f(it),hp, dhp, ireg    )
                     star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                     call vgrlp(sl, star, rp3, rp4, hmin, hp,
     2                    dhp, rl, drls, rv, drvs,0)
                     call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp1, dhp1, ireg    )   
                     star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
    
                     call vgrlp(sl, star, rp13f(it), rp14f(it), hmin,
     2                    hp1, dhp1, rl1, drls1, rv1, drvs1,0)
c                     hp = hp1 * porf + hp * (1. - porf) 
c                     dhp = dhp1 * porf + dhp * (1. - porf) 
                     permb = akf * porf + akm * (1. - porf)
                     rl=(akf*rl1*porf+akm*rl*(1.0-porf))/permb
                     drls=(akf*drls1*porf+akm*drls*(1.0-porf))/permb
                     rv = 1. - rl
                     drvs = -drls
                  endif
                  if(iupk.gt.0) then
c     upstream weighting of intrinsic permeability
                     permb = permb*darcyf
                     rl = rl*permb
                     rv = rv*permb
                     drls = drls*permb
                     drvs = drvs*permb
                     permb = 1.0/darcyf
c     set saturated permeabilities(assume isotropic)
                     pnx(mi) = permb*1.e+6
                     pnz(mi) = pnx(mi) * zfperm(it)
                     pny(mi) = pnx(mi) * yfperm(it)
                     pnx(mi) = pnx(mi) * xfperm(it)
                  else
c     set saturated permeabilities(assume isotropic)
                     if(porf.ge.0.0) then
                        pnx(mi) = permb*1.e+6
                        pny(mi) = pnx(mi) * yfperm(it)
                        pnz(mi) = pnx(mi) * zfperm(it)
                        pnx(mi) = pnx(mi) * xfperm(it)
                     endif
                  endif
c     set capillary pressures
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp
               elseif(irpd .eq. 6 .or. irpd .eq. 7) then
c     
c     akf-fracture saturated permeability(rp15)
c     akm-matrix saturated permeability(rp16)
c     porf-fracture fraction(rp17)
c     rlp(S)     
c     
                  akf = rp17f(it)
                  akm = rp18f(it)
                  porf = rp19f(it)
                  
                  if( idpdp .ne. 0 .or. idualp .ne. 0 ) then
                     
                     if( mi .le. neq )  then
c                        call vg_regions(2,ireg,mi,su_cut)
                        call vgcap( sl, rp11f(it), rp12f(it), rp13f(it),
     2                       rp14f(it), rp20f(it), rp21f(it), rp22f(it),
     3                       rp23f(it), rp16f(it), su_cut,    
     4                       cp3f(it),cp4f(it),hp, dhp, ireg    )
                        star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                        call  vgrlps(2, sl, star, rp13f(it),
     2                       rp14f(it), rp11f(it), rp12f(it),
     3                       tol_l, tol_u, rl, drls, rv, drvs)
                        permb = akf * porf
                        if(irpd.eq.7) then 
c calculate fracture term(from Sandia) if necessary
                           call rlp_frac(1,mi,sl,star,rp11f(it),
     &                          rp12f(it),rl,drls,1.-rl,-drls,rp24f(it))
                        endif
                        
                     else
c                        call vg_regions(1,ireg,mi,su_cut)
                        call vgcap( sl, rp1, rp2, rp3, rp4, 
     2                       rp7f(it), rp8f(it), rp9f(it),
     3                       rp10f(it), rp6f(it), su_cut,      
     4                       cp1f(it),cp2f(it),hp, dhp, ireg    )
                        star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                        call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                       tol_l, tol_u, rl, drls, rv, drvs )
                        permb = akm * (1. - porf)
                        if(irpd.eq.7) then
c calculate fracture term(from Sandia) if necessary
                           call rlp_frac(1,mi,sl,star,rp1,rp2,
     &                          rl,drls,1.-rl,-drls,rp24f(it))
                        endif
                     end if
                     
                  else
c                     call vg_regions(1,ireg,mi,su_cut)
                     call vgcap( sl, rp1, rp2, rp3, rp4,  
     2                    rp7f(it), rp8f(it), rp9f(it),
     3                    rp10f(it), rp6f(it), su_cut,      
     4                    cp1f(it),cp2f(it),hp, dhp, ireg    )
                     star = (sl-rp1f(it))/(rp2f(it)-rp1f(it))
                     call  vgrlps(2, sl, star, rp3, rp4, rp1, rp2, 
     2                    tol_l, tol_u, rl, drls, rv, drvs )
                     star = (sl-rp11f(it))/(rp12f(it)-rp11f(it))
                     call  vgrlps(2, sl, star, rp13f(it),
     2                    rp14f(it), rp11f(it), rp12f(it),
     3                    tol_l, tol_u, rl1, drls1, rv1, drvs1)
                     permb = akf * porf + akm * (1. - porf)
                     rl=(akf*rl1*porf+akm*rl*(1.0-porf))/permb
                     drls=(akf*drls1*porf+akm*drls*(1.0-porf))/permb
                     rv = 1. - rl
                     drvs = -drls
                  endif
c     upstream weighting of intrinsic permeability

                  if(iupk.gt.0) then
                     permb = permb*darcyf
                     rl = rl*permb
                     rv = rv*permb
                     drls = drls*permb
                     drvs = drvs*permb
                     permb = 1.0/darcyf
c     set saturated permeabilities(assume isotropic)
                     pnx(mi) = permb*1.e+6
                     pnz(mi) = pnx(mi) * zfperm(it)
                     pny(mi) = pnx(mi) * yfperm(it)
                     pnx(mi) = pnx(mi) * xfperm(it)
                  else
c     set saturated permeabilities(assume isotropic)
                     if(porf.ge.0.0) then
                        pnx(mi) = permb*1.e+6
                        pny(mi) = pnx(mi) * yfperm(it)
                        pnz(mi) = pnx(mi) * zfperm(it)
                        pnx(mi) = pnx(mi) * xfperm(it)
                     endif
                  endif
c     set capillary pressures
                  pcp(mi) = 9.8e-3 * hp
                  dpcef(mi) = 9.8e-3 * dhp

               elseif(irpd.eq.11) then
c Brooks-Corey
                  call bcrlp(sl,star,hp,dhp,rl,drls,rv,drvs)
               else if(irpd.eq.21) then
c Thomeer
c  -2 tests linear cap pressure
c               call thomeercap(-2,i,it,sl,star,hp,dhp,rl,drls,rv,drvs)
c call for capillary pressures
                  call thomeercap(2,i,it,sl,star,hp,dhp,rl,drls,rv,drvs)
                  pcp(mi) = hp
                  dpcef(mi) = dhp
c call for relative permeabilities
c  -3 tests relative permeabilities
c               call thomeercap(-3,i,it,sl,star,hp,dhp,rl,drls,rv,drvs)
c call for relative permeabilities
                  call thomeercap(3,i,it,sl,star,hp,dhp,rl,drls,rv,drvs)
         
               endif

            endif
            rlf(mi) = rl
            rvf(mi) = rv
            drlef(mi) = drls
            drvef(mi) = drvs
            drlpf(mi) = 0.0
            drvpf(mi) = 0.0
            if(ireg.eq.100.and.ireg0.ne.100) strd = strd
         enddo
c     
c     check perms if wellbore is enabled
c     
c        call welbor(2)

c
c add residual rlperm if directed
c GAZ 3-19-00
c
         if(irlp_fac.eq.1.and.iz.ne.0.and.rlp_flag.eq.1) then
            do i=1,n
               rlf(i) = rlf(i) + rlp_fac(i)
               rvf(i) = rvf(i) + rlp_fac(n0+i)
            enddo
         endif

      endif

      
      return
      end
