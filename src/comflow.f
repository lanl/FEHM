	module comflow
!     module comflow
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
!D1  PURPOSE
!D1
!D1  Include file for array variables and pointers for particle tracking flow.
!D1
!***********************************************************************
!D2
!D2  REVISION HISTORY 
!D2
!D2 $Log:   /pvcs.config/fehm90/src/comflow.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:42:30   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 08:56:48   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:05:34   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:22:32   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 11:56:52   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:39:42 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
!D2 
!D2    Rev 1.4   Wed May 29 14:28:28 1996   hend
!D2 Added variable diffusion with water content
!D2 
!D2    Rev 1.3   Thu Mar 21 13:20:08 1996   hend
!D2 Added variable for trac index ind. of istrw
!D2 
!D2    Rev 1.2   Mon Mar 04 15:58:54 1996   hend
!D2 Removed uneccessary calculations from coneq1 and added trac input option
!D2 
!D2    Rev 1.1   Tue Jan 16 15:55:14 1996   zvd
!D2 Added prolog.
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3  N/A
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4  
!D4  None
!D4  
!***********************************************************************

	real*8, allocatable :: a_axy(:),a_vxy(:),a_wvxy(:),alphaconl(:)
	real*8, allocatable :: alphaconv(:),sehvell(:),sehvelv(:)
	real*8, allocatable :: sehdiff(:),sehdiffv(:)
!       real*8, allocatable :: thetaresid(:,:), thetaresidv(:,:)
        integer dispsame,sehdonevel,hvliquid,hvvapor,sehsize,numflux

c s kelkar 3 July 2014, for calculating heat flow vectors
	logical flag_heat_out
	real*8, allocatable :: e_axy_adv(:), e_axy_cond(:)
	real*8, allocatable :: e_adv_nodal(:,:),e_cond_nodal(:,:)

	end module comflow

