      module comci
!    comci                                                    
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
!***********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 Global include file for array variables and pointers (FEHMN application).
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Revision                    ECD
!D2 Date         Programmer     Number  Comments
!D2
!D2 30-SEP-93    Z. Dash        22      Add prolog.
!D2              G. Zyvoloski           Initial implementation.
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comci.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:28   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:38   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:24   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:22   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:32   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:28 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.2   Tue Jan 09 14:33:18 1996   llt
!D2 increased drc to 10000, so can use with debugger
!D2 
!D2    Rev 1.1   03/18/94 16:22:58   gaz
!D2 Added solve_new and cleaned up memory management.
!D2 
!D2    Rev 1.0   01/20/94 10:21:54   pvcs
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
!D4   Identifier      Type     Block  Description
!D4
!D4   ***** COMMON Block fcc pointers and associated variables *****
!D4   ipakc           POINTER  fcc    Pointer to variable array akc
!D4   ipdanl          POINTER  fcc    Pointer to variable array danl
!D4   ipdanv          POINTER  fcc    Pointer to variable array danv
!D4   ipddva          POINTER  fcc    Pointer to variable array ddva
!D4   ipddvac         POINTER  fcc    Pointer to variable array ddvac
!D4   ipddvae         POINTER  fcc    Pointer to variable array ddvae
!D4   ipddvap         POINTER  fcc    Pointer to variable array ddvap
!D4   ipdeef          POINTER  fcc    Pointer to variable array deef
!D4   ipdelef         POINTER  fcc    Pointer to variable array delef
!D4   ipdelf          POINTER  fcc    Pointer to variable array delf
!D4   ipdenci         POINTER  fcc    Pointer to variable array denci
!D4   ipdenei         POINTER  fcc    Pointer to variable array denei
!D4   ipdeni          POINTER  fcc    Pointer to variable array deni
!D4   ipdepf          POINTER  fcc    Pointer to variable array depf
!D4   ipdeqh          POINTER  fcc    Pointer to variable array deqh
!D4   ipdevef         POINTER  fcc    Pointer to variable array devef
!D4   ipdevf          POINTER  fcc    Pointer to variable array devf
!D4   ipdgle          POINTER  fcc    Pointer to variable array dgle
!D4   ipdglp          POINTER  fcc    Pointer to variable array dglp
!D4   ipdgve          POINTER  fcc    Pointer to variable array dgve
!D4   ipdgvp          POINTER  fcc    Pointer to variable array dgvp
!D4   ipdil           POINTER  fcc    Pointer to variable array dil
!D4   ipdile          POINTER  fcc    Pointer to variable array dile
!D4   ipdilp          POINTER  fcc    Pointer to variable array dilp
!D4   ipdiv           POINTER  fcc    Pointer to variable array div
!D4   ipdive          POINTER  fcc    Pointer to variable array dive
!D4   ipdivp          POINTER  fcc    Pointer to variable array divp
!D4   ipdmef          POINTER  fcc    Pointer to variable array dmef
!D4   ipdmpf          POINTER  fcc    Pointer to variable array dmpf
!D4   ipdpcef         POINTER  fcc    Pointer to variable array dpcef
!D4   ipdq            POINTER  fcc    Pointer to variable array dq
!D4   ipdqh           POINTER  fcc    Pointer to variable array dqh
!D4   ipdrc           POINTER  fcc    Pointer to variable array drc
!D4   ipdstm          POINTER  fcc    Pointer to variable array dstm
!D4   ipdtpa          POINTER  fcc    Pointer to variable array dtpa
!D4   ipdtpae         POINTER  fcc    Pointer to variable array dtpae
!D4   ipenlf          POINTER  fcc    Pointer to variable array enlf
!D4   ipenva          POINTER  fcc    Pointer to variable array enva
!D4   ipenvac         POINTER  fcc    Pointer to variable array envac
!D4   ipenvae         POINTER  fcc    Pointer to variable array envae
!D4   ipenvap         POINTER  fcc    Pointer to variable array envap
!D4   ipenvf          POINTER  fcc    Pointer to variable array envf
!D4   ipineluf        POINTER  fcc    Pointer to variable array ineluf
!D4   ipnsf           POINTER  fcc    Pointer to variable array nsf
!D4   iprolf          POINTER  fcc    Pointer to variable array rolf
!D4   iprovf          POINTER  fcc    Pointer to variable array rovf
!D4 
!D4   akc             REAL*8   fcc    Tracer accumulation term derivative with
!D4                                     respect to total concentration 
!D4   danl            REAL*8   fcc    Derivative of liquid phase concentration
!D4                                     with respect to total concentration 
!D4   danv            REAL*8   fcc    Derivative of vapor phase concentration
!D4                                     with respect to total concentration 
!D4   ddva            REAL*8   fcc    ?
!D4   ddvac           REAL*8   fcc    ?
!D4   ddvae           REAL*8   fcc    ?
!D4   ddvap           REAL*8   fcc    ?
!D4   deef            REAL*8   fcc    Derivative of energy accumulation with
!D4                                     respect to energy variable 
!D4   denci           REAL*8   fcc    Tracer accumulation term 
!D4   denei           REAL*8   fcc    Energy accumulation term 
!D4   deni            REAL*8   fcc    Mass accumulation term 
!D4   depf            REAL*8   fcc    Derivative of energy accumulation term
!D4                                     with respect to pressure  
!D4   deqh            REAL*8   fcc    Derivative of energy source term with
!D4                                     respect to energy variable 
!D4   devef           REAL*8   fcc    Derivative of enthalpy with respect to
!D4                                     energy variable  
!D4   devf            REAL*8   fcc    Derivative of enthalpy with respect to
!D4                                     pressure 
!D4   delef           REAL*8   fcc    Derivative of liquid enthalpy with
!D4                                     respect to energy variable 
!D4   delf            REAL*8   fcc    Derivative of liquid enthalpy with
!D4                                     respect to pressure  
!D4   dgle            REAL*8   fcc    Derivative of liquid mass gravity term
!D4                                     with respect to energy variable 
!D4   dglp            REAL*8   fcc    Derivative of liquid mass gravity term
!D4                                     with respect to pressure  
!D4   dgve            REAL*8   fcc    Derivative of vapor mass gravity term
!D4                                     with respect to energy variable 
!D4   dgvp            REAL*8   fcc    Derivative of vapor mass gravity term
!D4                                     with respect to pressure  
!D4   dil             REAL*8   fcc    Liquid transmissibility 
!D4   dile            REAL*8   fcc    Derivative of liquid transmissibility
!D4                                     with respect to energy variable 
!D4   dilp            REAL*8   fcc    Derivative of liquid transmissibility
!D4                                     with respect to pressure  
!D4   div             REAL*8   fcc    Vapor transmissibility 
!D4   dive            REAL*8   fcc    Derivative of vapor transmissibility
!D4                                     with respect to energy variable 
!D4   divp            REAL*8   fcc    Derivative of vapor transmissibility
!D4                                     with respect to pressure  
!D4   dmef            REAL*8   fcc    Derivative of mass accumulation term
!D4                                     with respect to energy variable 
!D4   dmpf            REAL*8   fcc    Derivative of mass accumulation term
!D4                                     with respect to pressure  
!D4   dpcef           REAL*8   fcc    Derivative of capillary pressure with
!D4                                     respect to the energy variable 
!D4   dq              REAL*8   fcc    Derivative of mass source term with
!D4                                     respect to pressure 
!D4   dqh             REAL*8   fcc    Derivative of energy source term with
!D4                                     respect to pressure  
!D4   drc             REAL*8   fcc    Tracer source term derivative with
!D4                                     respect to total concentration 
!D4   dstm            REAL*8   fcc    Steam mass 
!D4   dtpa            REAL*8   fcc    Derivative of temperature with respect
!D4                                     to pressure  
!D4   dtpae           REAL*8   fcc    Derivative of temperature with respect
!D4                                     to energy variable  
!D4   enlf            REAL*8   fcc    Liquid enthalpy 
!D4   enva            REAL*8   fcc    ?
!D4   envac           REAL*8   fcc    ?
!D4   envae           REAL*8   fcc    ?
!D4   envap           REAL*8   fcc    ?
!D4   envf            REAL*8   fcc    Vapor enthalpy 
!D4   ineluf          INT      fcc    ?
!D4   nsf             INT      fcc    ?
!D4   rolf            REAL*8   fcc    Liquid density 
!D4   rovf            REAL*8   fcc    Vapor density 
!D4
!D4 Global Subprograms
!D4
!D4   None
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
!D9 N/A
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

      integer nbigblock, nstorepiv
      
      real*8, allocatable, target :: bigblock(:)
      real*8, pointer :: b(:)
      real*8, pointer :: gmres(:)
      real*8, pointer :: piv(:)
!      real*8, pointer ::  akc(:)
!      real*8, pointer ::  danl(:) 
!      real*8, pointer ::   danv(:) 
      real*8, pointer ::   ddvac(:)
      real*8, pointer ::   ddvae(:) 
      real*8, pointer ::   ddvap(:) 
      real*8, pointer ::  deef(:)
      real*8, pointer ::   delef(:)
      real*8, pointer ::   delf(:) 
!      real*8, pointer ::   denci(:)
      real*8, pointer ::   denei(:) 
      real*8, pointer ::   deni(:) 
      real*8, pointer ::   denvac(:)
      real*8, pointer ::   denvae(:) 
      real*8, pointer ::   denvap(:) 
      real*8, pointer ::  depf(:)
      real*8, pointer ::   deqh(:) 
      real*8, pointer ::   devef(:)
      real*8, pointer ::   devf(:) 
      real*8, pointer ::   dgle(:)
      real*8, pointer ::   dglp(:) 
      real*8, pointer ::   dgve(:) 
      real*8, pointer ::   dgvp(:)
!      real*8, pointer ::   dil(:) 
      real*8, pointer ::   dile(:) 
      real*8, pointer ::   dilp(:) 
!      real*8, pointer ::   div(:)
      real*8, pointer ::   dive(:) 
      real*8, pointer ::   divp(:) 
      real*8, pointer ::  dmef(:)
      real*8, pointer ::  dmpf(:)
!      real*8, pointer ::   dpcef(:) 
      real*8, pointer ::   dq(:)
      real*8, pointer ::  dqh(:)
      real*8, pointer ::   dqt(:)
      real*8, pointer ::   drc(:) 
      real*8, pointer ::   dstm(:) 
      real*8, pointer ::   dtpa(:) 
      real*8, pointer ::   dtpae(:) 
      real*8, pointer ::   dva(:) 
      real*8, pointer ::   enlf(:) 
      real*8, pointer ::   enva(:)
      real*8, pointer ::   envf(:)
!      real*8, pointer ::   rolf(:) 
      logical cden
      integer ispcden
      real*8 factcden
!      real*8, pointer ::   rovf(:) 

! Rel. perm. shared space with derivative arrays
      real*8, pointer :: rlf(:)
      real*8, pointer :: drlef(:)
      real*8, pointer :: drlpf(:)
      real*8, pointer :: rvf(:)
      real*8, pointer :: drvef(:)
      real*8, pointer :: drvpf(:)

! Not using shared space
      real*8, pointer ::   dvas(:)
      real*8, allocatable ::   denci(:)
      real*8, allocatable ::   rolf(:) 
      real*8, allocatable ::   rovf(:) 
      real*8, allocatable ::   dil(:) 
      real*8, allocatable ::   div(:)
      real*8, allocatable ::   akc(:)
      real*8, allocatable ::   danl(:) 
      real*8, allocatable ::   danv(:) 
      real*8, allocatable ::   dpcef(:)
c gaz 050620  added pure water versions of rolf and dil  
      real*8, allocatable ::   rolf_pure(:)  
      real*8, allocatable ::   dil_pure(:) 
! allocatable storage for saved LU factors
      real*8, allocatable :: bsave(:)
      real*8, allocatable :: pivsave(:)
      real*8,  allocatable ::   den_spatial(:)
      real*8,  allocatable ::   vis_spatial(:)
      real*8,  allocatable ::   comp_spatial(:)
      real*8,  allocatable ::   deng_spatial(:)
      real*8,  allocatable ::   visg_spatial(:)
c gaz 100517
      real*8,  allocatable :: dva_save(:)     
c gaz debug 090511
      real*8, allocatable  ::   denei_ch(:) 
      real*8, allocatable  ::   deni_ch(:) 
      real*8, allocatable  ::   denpci_ch(:) 
      integer, allocatable  ::   ieos_bal(:) 

      end module comci
