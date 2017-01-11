real*8 function cden_correction(i)
  !***********************************************************************
  ! Copyright 2011 Los Alamos National Security, LLC  All rights reserved
  ! Unless otherwise indicated,  this information has been authored by an
  ! employee or employees of the Los Alamos National Security, LLC (LANS),
  ! operator of the  Los  Alamos National  Laboratory  under Contract  No.
  ! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
  ! Government   has   rights  to  use,  reproduce,  and  distribute  this
  ! information.  The  public may copy  and  use this  information without
  ! charge, provided that this  Notice and any statement of authorship are
  ! reproduced on all copies.  Neither  the  Government nor LANS makes any
  ! warranty,   express   or   implied,   or   assumes  any  liability  or
  ! responsibility for the use of this information.      
  !***********************************************************************       
  ! Calculate concentration-dependent density correction
  ! requires cden_flag, node number (i).  requires ispcden to have been identified
  ! cden_flag=0  original formulation, based on C of single tracer (ispcden) and F (input parameter in cden macro)
  !      cden_correction = C*F
  ! cden_flag=1  ( based on molecular weight table provided in trxn macro )
  !      cden_correction = sum over all mw table species: C*mw*700*1.E-3
  ! cden_flag=2  (co2 keyword)
  !      cden_correction = Concentration of species ispcden (used to calculate density in ther_co2_h2o.f)
  !             by default will be Na but could be any species if named in cden macro

  use comci, only : ispcden, factcden
  use comdi, only : icns, nspeci, anl
  use comdti, only : n0
  use comrxni, only : cden_flag,  mw_speci

  implicit none

  integer :: i, i2, iaq, ispeci
  real(8) :: cofdc

  parameter (cofdc = 700.d0)
  
  cden_correction = 0.d0
  iaq = 0
  select case (cden_flag)
  case (0)
     ! Single species density correction (original implementation)
     i2 = i + (ispcden - 1)*n0
     if (anl(i2) .gt. 0.d0) cden_correction = anl(i2) * factcden

  case (1)
     ! All acqueous species, C - moles/kg-water, MW - g/mole)
     do ispeci = 1, nspeci
        select case (icns(ispeci))
        case (1, 2, -2)
           iaq = iaq + 1
           i2 = i + (ispeci - 1)*n0
           if (anl(i2) .gt. 0.d0) cden_correction = cden_correction + anl(i2) * mw_speci(iaq)
        end select
     end do
     ! cofdc - Coefficient of fluid density change drho/dC = 700 kg m^3
     cden_correction = cofdc * cden_correction * 1.d-3
  case (2)
     i2 = i + (ispcden - 1)*n0
     if (anl(i2) .gt. 0.d0) cden_correction = anl(i2)
 
  end select

end function cden_correction
