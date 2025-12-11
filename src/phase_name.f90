function ic(name, ictype)
!*************************************************************************
! Copyright  2015.   Los Alamos National Security, LLC.  This material was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos
! National  Laboratory  (LANL),  which is operated by  Los Alamos National
! Security, LLC  for the U. S. Department of Energy.  The U. S. Government
! has rights to use, reproduce, and distribute this software.  Neither the
! U. S. Government nor Los Alamos National Security, LLC or persons acting
! on their behalf,  make any warranty,  express or implied, or assumes any
! liability for the accuracy, completeness, or usefulness of the software,
! any information pertaining to the software,  or  represents that its use
! would not infringe privately owned rights.

! The  software  being licensed  may  be Export Controlled.  It may not be
! distributed  or  used by individuals  or entities prohibited from having
! access to the software package, pursuant to United States export control
! laws and regulations. An export control review and determination must be
! completed before LANS will provide access to the identified Software.
!*************************************************************************


  use comco2, only: icarb
  use comrlp, only: rlp_type, cap_pos
  implicit none
  character*30 name
  integer ic, ictype(3)
  
!      phase_comp(20,1)=3;phase_comp(20,2)=6  water/co2_liquid
!      phase_comp(21,1)=3;phase_comp(21,2)=4   water/co2_gas
!      phase_comp(22,1)=6;phase_comp(22,2)=4   co2_liquid/co2_gas 
!      phase_comp(23,1)=1;phase_comp(23,2)=7   oil/water
!      phase_comp(24,1)=1;phase_comp(24,2)=8    gas/water
!      phase_comp(25,1)=1;phase_comp(25,2)=2     air/water
!      phase_comp(26,1)=1;phase_comp(26,2)=6;phase_comp(26,3)=4    water/co2_liquid/co2_gas
!      phase_comp(28,1)=2;phase_comp(28,2)=5  water/vapor
  

  ictype = 0
  select case(name)	
  case('water','WATER') 
     if(icarb.eq.0) then
        ic=2
     else
        ic=3
     endif
  case ('air','AIR')
     ic = 1
  case ('h2o_gas','H2O_GAS','vapor','VAPOR')
     ic = 5
  case ('co2_gas','CO2_GAS')
     ic = 4
  case ('co2_liquid','CO2_LIQUID')
     ic = 6
   case ('oil', 'OIL')
     ic = 7
  case ('gas', 'GAS')
     ic = 8   
  case ('co2_sc','CO2_SC')
     ic = 9
  case ('methane_hydrate','METHANE_HYDRATE')
     ic = 11
   case ('oil_water', 'OIL_WATER' , 'water/oil', 'oil/water')
     ic = 23
     ictype(1) = 2
     ictype(2) = 7
  case ('oil_gas', 'OIL_GAS' , 'gas/oil', 'oil/gas')
     ic = 27
     ictype(1) = 7
     ictype(2) = 8
  case ('air/water', 'water/air')
     ic = 25
     ictype(1) = 1
     ictype(2) = 2
  case ('water/co2_liquid', 'co2_liquid/water')
     ic = 20
     ictype(1) = 3
     ictype(2) = 6
  case ('water/co2_gas', 'co2_gas/water')
     ic = 21
     ictype(1) = 3
     ictype(2) = 4
  case ('co2_liquid/co2_gas', 'co2_gas/co2_liquid')
     ic = 22
     ictype(1) = 4
     ictype(2) = 6
  case ('water/vapor', 'vapor/water')
     ic = 28
     ictype(1) = 2
     ictype(2) = 5
  case ('air/vapor', 'vapor/air')
     ic = 29
     ictype(1) = 1
     ictype(2) = 5
  case ('gas/water', 'water/gas')
     ic = 24
     ictype(1) = 2
     ictype(2) = 7
  case ('water/co2_liquid/co2_gas')
     ic=26               		
     ictype(1) = 3
     ictype(2) = 4
     ictype(3) = 6
 case ('liquid/gas', 'gas/liquid')
     ! not sure what this means            
     ic = 10
  case default
     write(*,*) 'phase ',name,' not recognized'
     stop
  end select
  return
end function ic
function pn(k)
  implicit none
  character*30 pn
  integer k
  select case(k)
  case(2,3)
     pn='water'
  case(1)	
     pn='air'
  case(5)
     pn='h2o_gas'
  case(4)
     pn='co2_gas'
  case(6)
     pn='co2_liquid'
  case(9)
     pn='co2_sc'
  case(11)
     pn='methane_hydrate'
  case(7)
     pn='oil'
  case(8)
     pn='gas'
  case(23)
     pn= 'water/oil'
  case(27)
     pn= 'gas/oil'
  case(25)
     pn='air/water'
  case(20)
     pn='water/co2_liquid'
  case(21)
     pn='water/co2_gas'
  case(22)
     pn='co2_liquid/co2_gas'
  case(28)
     pn='water/vapor'
  case(29)
     pn='air/vapor'
  case(24)
     pn='gas/water'
  case(26)
     pn='water/co2_liquid/co2_gas'
  case(10)             		
     pn='liquid/gas'   ! not sure what this means    
  case default
     pn='unrecognized phase'
  end select
  return
end function pn
