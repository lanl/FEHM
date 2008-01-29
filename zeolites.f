      subroutine zeolites(temperature,
     2     kappa_zeolites, theta_old, theta_new, water_pressure,
     3     dpdt, sk_zeolites,
     4     dsk_zeolites_t, qh_zeolites, dqh_zeolites_t)
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
C***********************************************************************
CD1
CD1  PURPOSE
CD1
CD1  This subroutine computes the energy/water source/sink terms
CD1  due to a zeolite dyhydration reaction
CD1
C***********************************************************************
CD2
CD2  REVISION HISTORY 
CD2
CD2  $Log:   /pvcs.config/fehm90/src/zeolites.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:25:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
!D2 
!D2    Rev 2.3   14 Nov 2001 13:29:30   pvcs
!D2 FEHM Version 2.12, STN 10086-2.12-00
!D2 
!D2    Rev 2.2   06 Jun 2001 13:28:54   pvcs
!D2 FEHM Version 2.11, STN 10086-2.11-00
!D2 
!D2    Rev 2.1   30 Nov 2000 12:13:44   pvcs
!D2 FEHM Version 2.10, STN 10086-2.10-00
!D2 
!D2    Rev 2.0   Fri May 07 14:48:34 1999   pvcs
!D2 FEHM Version 2.0, SC-194 (Fortran 90)
CD2 
CD2    Rev 1.3   Fri Feb 16 13:59:46 1996   zvd
CD2 Added requirement.
CD2 
CD2    Rev 1.2   Thu Jan 11 13:10:40 1996   hend
CD2 Added ECD Number
CD2 
CD2    Rev 1.1   Thu Jan 11 13:06:36 1996   hend
CD2 Added Prolog
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3  2.3.2 Heat- and mass-transfer equations
CD3  2.3.7 Sources and sinks
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C***********************************************************************

      implicit none

      real*8 temperature,kappa_zeolites,theta_old,theta_new
      real*8 sk_zeolites,dsk_zeolites_t,qh_zeolites
      real*8 dqh_zeolites_t,t_kelvin,water_pressure,dpdt
      real*8 dtheta_new_t,df_dtheta,logwp,fterm2,f_theta
      real*8 a_theta, b_theta, c_theta, d_theta
      real*8 a_qh,b_qh,c_qh,twob_theta,twoc_theta,twob_qh,threec_qh
      parameter( a_theta = 3.41674,
     2           b_theta = -0.0703461,
     3           c_theta = 0.000549001,
     4           d_theta = 0.0499952 )
      parameter( a_qh = -0.120480,
     2           b_qh = 0.068535,
     3           c_qh = -0.02498 )
      parameter( twob_theta = 2.*b_theta)
      parameter( twoc_theta = 2.*c_theta)
      parameter( twob_qh = 2.*b_qh)
      parameter( threec_qh = 3.*c_qh)

      logwp = dlog(water_pressure)
      t_kelvin = temperature + 273.16
      theta_new = a_theta +
     2     b_theta*dlog(1./t_kelvin)**2 +
     3     c_theta*logwp**2 +
     4     d_theta*logwp
      sk_zeolites = kappa_zeolites * (theta_new - theta_old)
      dtheta_new_t = twob_theta*dlog(t_kelvin)/t_kelvin +
     2     d_theta*dpdt/water_pressure +
     3     twoc_theta*logwp*dpdt/water_pressure
      dsk_zeolites_t = kappa_zeolites*dtheta_new_t
      fterm2 = theta_new*theta_new
      f_theta = a_qh*(theta_new-theta_old) + 
     2     b_qh*(fterm2 - theta_old*theta_old) + 
     3     c_qh*(fterm2*theta_new -
     4     theta_old*theta_old*theta_old)
      qh_zeolites = kappa_zeolites*f_theta/0.018152
      df_dtheta = a_qh+twob_qh*theta_new+threec_qh*fterm2
      dqh_zeolites_t = df_dtheta*dsk_zeolites_t/0.018152

      return
      end
