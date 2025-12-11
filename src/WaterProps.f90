subroutine WaterProps(Temp,Pres_in,Select,NeglectPres,Answer, &
                      TempDeriv,PresDeriv)
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
! Purpose: Calculate the density, enthalpy, or viscosity of
!          H2O between -30 and 30 C, taking into account
!          the liquid-solid phase transition 

! Inputs: Temp = Temperature (C), Pres = Pressure (MPa),
!         Select = Selection code, 
!         NeglectPres = 1 to neglect the pressure 
!                     = 0 to include pressure dependence
! Output: Answer = density   (kg/m^3) for Select = 0
!                  enthalpy  (MJ/kg)  ...        = 1
!                  viscosity (Pa*s)   ...        = 2
!         TempDeriv = derivative wrt temperature of density
!                     enthalpy, or viscosity
!         PresDeriv = derivative wrt pressure of ...
use comai, only : ierr
implicit none
real*8 :: Temp,Pres,Answer,IceRefDens,LiqRefDens,IceEndTemp, &
   LiqEndTemp,Slope,CpIce,CpLiq,IceRefEnth,Xi,               &
   LiqRefEnth,IceVisc,RefLiqVisc,IceCompress,                &
   LiqCompress,CompressConvert,RefPres,MPaToPa,              &
   LiqExpansion,IceExpansion,LiqEndVisc,IceEndVisc,          &
   LiqDCpDTemp,IceDCpDTemp,Pres_in,LargePos,LargeNeg,        &
   LiqDEnthDTemp,IceDEnthDTemp,                              &
   LiqDViscDTemp,LiqEndDens,IceEndDens,TempDeriv,            &
   PresDeriv,LiqEndEnth,IceEndEnth,Vis_tol
integer :: Select,NeglectPres
parameter( Vis_tol= 1.d-6)

!----------------------DATA SECTION------------------------
! NOTE: 'Ref' in variable names refers to T = phase change 
!       temperature and zero pressure; 'End' refers to 
!       phase change temperatures
IceRefDens = 917.0              ! in kg/m^3
LiqRefDens = 1000.0
IceEndTemp = -0.10              ! in C
LiqEndTemp =  0.10
CpLiq = 4.182                   ! in kJ/kg*C
CpIce = 2.108
IceRefEnth = 0                  ! in kJ/kg
LiqRefEnth = 334.0
IceVisc = 1.E20                 ! in Pa*s
RefLiqVisc = 1.787E-3

if (NeglectPres==0) then
  IceCompress = 1.004E-3        ! in cm^3/(mol*bar)
  LiqCompress = 8.584E-4
  CompressConvert = 1.0/(18*1.E8)
  IceCompress = IceCompress*CompressConvert
  LiqCompress = LiqCompress*CompressConvert
  ! Now IceCompress and LiqCompress are compressibilities in
  ! 1/Pa when multiplied by density
else if (NeglectPres==1) then
  IceCompress = 0.0
  LiqCompress = 0.0
else
  write(*,*) 'NeglectPres is not specified - WaterProps'
  write(*,*) 'NeglectPres: ',NeglectPres
  write(ierr,*) 'NeglectPres is not specified - WaterProps'
  write(ierr,*) 'NeglectPres: ',NeglectPres  
!  STOP
end if
RefPres = 1.E5 ! in Pa
MPaToPa = 1.E6
LiqExpansion = -0.011 ! in kg/(m^3*C)
IceExpansion = -0.26
LiqDCpDTemp = -0.0012 ! in kJ/(kg*K^2)
IceDCpDTemp =  0.0046  
LiqDEnthDTemp = CpLiq + Temp*LiqDCpDTemp
IceDEnthDTemp = CpIce + Temp*IceDCpDTemp
LiqDViscDTemp = -4E-5

! Frozen wall benchmark only !!
!LiqDViscDTemp = 0.0
!RefLiqVisc = 0.0015
!LiqCompress = 0.0
!IceExpansion = 0.0
!LiqExpansion = 0.0
!IceCompress = 0.0
!LiqDCpDTemp = 0.0
!IceDCpDTemp = 0.0
!IceRefDens = 920.0

LargePos = 10.0
LargeNeg = -10.0
!----------------------------------------------------------

if ((Temp<-150.0).or.(Temp>150.0)) then
  ! Return error
  write(*,*) 'Temperature out of bounds --WaterProps'
  write(*,*) 'Temperature: ',Temp
  write(ierr,*) 'Temperature out of bounds --WaterProps'
  write(ierr,*) 'Temperature: ',Temp 
!  STOP
else
  ! Convert Pres to Pa
  Pres = Pres_in*MPaToPa
  if (Select==0) then
    !-----------------CALCULATE DENSITY--------------------
    if (Temp < IceEndTemp) then
      Answer = IceRefDens + &
               IceExpansion*(Temp-IceEndTemp) + &
               IceRefDens**2.0*IceCompress*(Pres-RefPres)
      TempDeriv = IceExpansion
      PresDeriv = IceRefDens**2.0*IceCompress
    elseif ((Temp >= IceEndTemp).and. &
            (Temp <= LiqEndTemp)) then
      LiqEndDens = LiqRefDens + &
        LiqRefDens**2.0*LiqCompress*(Pres-RefPres)
      IceEndDens = IceRefDens + &
        IceRefDens**2.0*IceCompress*(Pres-RefPres)
      Slope = (LargePos-LargeNeg)/(LiqEndTemp-IceEndTemp)
      Xi = Slope*(Temp-IceEndTemp)+LargeNeg
      Answer = (LiqEndDens-IceEndDens)*(tanh(Xi)+1.0)/2.0 + &
        IceEndDens
      TempDeriv = ((LiqEndDens-IceEndDens)/2.0)*Slope* &
        (1.0-tanh(Xi)**2.0)
      PresDeriv = (LiqRefDens**2.0*LiqCompress - &
                   IceRefDens**2.0*IceCompress)* &
                  (tanh(Xi)-1.0)/2.0 + IceRefDens**2.0*IceCompress
    else
      Answer = LiqRefDens + &
               LiqExpansion*(Temp-LiqEndTemp) + &
               LiqRefDens**2.0*LiqCompress*(Pres-RefPres)
      TempDeriv = LiqExpansion
      PresDeriv = LiqRefDens**2.0*LiqCompress
    end if
  elseif (Select==1) then
    !-----------------CALCULATE ENTHALPY-------------------
    if (Temp < IceEndTemp) then
      Answer = IceRefEnth + &
               IceDEnthDTemp*(Temp-IceEndTemp)
      TempDeriv = IceDEnthDTemp
      PresDeriv = 0.0
    elseif ((Temp >= IceEndTemp).and. &
            (Temp <= LiqEndTemp)) then
      LiqEndEnth = LiqRefEnth
      IceEndEnth = IceRefEnth
      Slope = (LargePos-LargeNeg)/ &
              (LiqEndTemp-IceEndTemp)
      Xi = Slope*(Temp-IceEndTemp) + LargeNeg
      Answer = (LiqEndEnth-IceEndEnth)*(tanh(Xi)+1.0)/2.0 + &
        IceEndEnth
      TempDeriv = ((LiqEndEnth-IceEndEnth)/2.0)*Slope* &
        (1.0-tanh(Xi)**2.0)
      PresDeriv = 0.0
    else
      Answer = LiqRefEnth + &
               LiqDEnthDTemp*(Temp-LiqEndTemp)
      TempDeriv = LiqDEnthDTemp
      PresDeriv = 0.0
    end if
    ! Convert to MJ/kg
    Answer = Answer/1000.0  
    TempDeriv = TempDeriv/1000.0
    PresDeriv = PresDeriv/1000.0
  elseif (Select==2) then
    !-----------------CALCULATE VISCOSITY-----------------
    if (Temp < IceEndTemp) then
      Answer = IceVisc
      TempDeriv = 0.0
      PresDeriv = 0.0
    elseif ((Temp >= IceEndTemp).and. &
            (Temp <= LiqEndTemp)) then
      LiqEndVisc = RefLiqVisc
      IceEndVisc = IceVisc
      Slope = (LargePos-LargeNeg)/(LiqEndTemp-IceEndTemp)
      Xi = Slope*(Temp-IceEndTemp)+LargeNeg
      Answer = (LiqEndVisc-IceEndVisc)*(tanh(Xi)+1.0)/2.0 + &
        IceEndVisc
      TempDeriv = ((LiqEndVisc-IceEndVisc)/2.0)*Slope* &
        (1.0-tanh(Xi)**2.0)
      PresDeriv = 0.0
    else
      Answer = RefLiqVisc + LiqDViscDTemp*(Temp-LiqEndTemp)
      TempDeriv = LiqDViscDTemp
      PresDeriv = 0.0
    end if
    if(Answer.lt.Vis_tol) then
     Answer = Vis_tol
     TempDeriv = 0.0
     PresDeriv = 0.0
    endif
  else
    ! Return error
    write(*,*) 'Selection code not set properly --WaterProps'
    write(*,*) 'Selection code: ',Select
    STOP
  end if
end if

end

subroutine WaterProps_Vap(Temp,Pres_in,Select,NeglectPres,Answer, &
                      TempDeriv,PresDeriv)
!  Purpose: Calculate the density, enthalpy, or viscosity of
!          H2O between -30 and 30 C, taking into account
!          the Vapuid-solid phase transition 

! Inputs: Temp = Temperature (C), Pres = Pressure (MPa),
!         Select = Selection code, 
!         NeglectPres = 1 to neglect the pressure 
!                     = 0 to include pressure dependence
! Output: Answer = density   (kg/m^3) for Select = 0
!                  enthalpy  (MJ/kg)  ...        = 1
!                  viscosity (Pa*s)   ...        = 2
!         TempDeriv = derivative wrt temperature of density
!                     enthalpy, or viscosity
!         PresDeriv = derivative wrt pressure of ...
use comai, only : iout, iptty
use comdi, only : LiqEndTemp, VapEndTemp
implicit none
real*8 :: Temp,Pres,Answer,LiqRefDens,VapRefDens,          &
   Slope,CpLiq,CpVap,LiqRefEnth,                           &
   VapRefEnth,LiqVisc,RefLiqVisc,RefVapVisc,LiqCompress,   &
   VapCompress,CompressConvert,RefPres,MPaToPa,            &
   VapExpansion,LiqExpansion,VapEndVisc,LiqEndVisc,        &
   VapDCpDTemp,LiqDCpDTemp,Pres_in,                        &
   VapDEnthDTemp,LiqDEnthDTemp,LiqRefTemp,VapRefTemp,      &
   VapDViscDTemp,VapEndDens,LiqEndDens,TempDeriv,          &
   PresDeriv,VapEndEnth,LiqEndEnth,Vis_tol,                &
   VapTransDens,LiqTransDens,VapTransEnth,LiqTransEnth,    &
   VapTransVisc,LiqTransVisc,LiqViscTerm1,LiqViscTerm2,    &
   VapRefVisc,LiqRefVisc,DslopeP
integer :: Select,NeglectPres
parameter( Vis_tol= 1.d-6)

!----------------------DATA SECTION------------------------
! NOTE: 'Ref' in variable names refers to T = phase change 
!       temperature and zero pressure; 'End' refers to 
!       phase change temperatures
!   Pressure reference is 2.32 Mpa
LiqRefDens = 841.8             ! in kg/m^3 at 20 C
VapRefDens = 11.62             ! at 20 C
! LiqEndTemp = 218.              ! in C
! VapEndTemp =  222.
! gaz 032011 set test temps at 152 and 154
!LiqEndTemp = 151.              ! in C
!VapEndTemp =  155.
LiqRefTemp = LiqEndTemp
VapRefTemp = VapEndTemp
LiqTransDens = 841.8           ! at 218 C
VapTransDens = 11.52           ! at 222 C
LiqTransEnth =  0.934          ! at 218 C
VapTransEnth =  2.806          ! at 222 C
CpVap = 2.8e-3                    ! in MJ/kg*C 220-280 C
CpLiq = 4.28e-3                   !  20-220 C
! Ref values modified for GAZ application
LiqRefEnth = 0.934                 ! in MJ/kg
VapRefEnth = 2.806
LiqViscTerm2 = 0.0193           ! b in visc = a + b/T
LiqRefVisc = 1.19E-4             ! in Pa*s at 218 C
LiqViscTerm1 = LiqRefVisc-LiqViscTerm2/LiqEndTemp ! a in visc = a + b/T
VapRefVisc = 1.65E-6           ! at 222 C
if (NeglectPres==0) then
  LiqCompress = 0.5        ! kg/(m3*Mpa)
  VapCompress = 0.
  ! Now LiqCompress and VapCompress are compressibilities in
  ! 1/Pa when multiplied by density
else if (NeglectPres==1) then
  LiqCompress = 0.0
  VapCompress = 0.0
else
  write(*,*) 'NeglectPres is not specified - WaterProps_Vap'
  write(*,*) 'NeglectPres: ',NeglectPres
  STOP
end if
RefPres = 2.32! in Mpa
MPaToPa = 1.E6
VapExpansion = -0.011 ! in kg/(m^3*C)
LiqExpansion = -0.794
VapDCpDTemp = 0.0
LiqDCpDTemp =  0.0
VapDEnthDTemp = CpLiq + Temp*VapDCpDTemp
LiqDEnthDTemp = CpVap + Temp*LiqDCpDTemp
VapDViscDTemp = 5.e-9 
!----------------------------------------------------------

if (Temp>600) then
  ! Return error
  if(iout.ne.0) then
   write(iout,*) 'Temperature out of bounds --WaterPropsVap'
   write(iout,*) 'Temperature: ',Temp, ' Stopping'
  endif
  if(iptty.ne.0) then
   write(iptty,*) 'Temperature out of bounds --WaterPropsVap'
   write(iptty,*) 'Temperature: ',Temp, ' Stopping'
  endif
  STOP
else
  ! Convert Pres to Pa 
  Pres = Pres_in
  if (Select==0) then
    !-----------------CALCULATE DENSITY--------------------
    If (Temp < 0.0) then
     Answer = LiqRefDens + &
               LiqExpansion*(-LiqRefTemp) + &
               LiqCompress*(Pres-RefPres)
       TempDeriv = 0.0
       PresDeriv = LiqCompress
    elseif (Temp < LiqEndTemp) then
      Answer = LiqRefDens + &
               LiqExpansion*(Temp-LiqRefTemp) + &
               LiqCompress*(Pres-RefPres)
      TempDeriv = LiqExpansion
      PresDeriv = LiqCompress
    elseif ((Temp >= LiqEndTemp).and. &
            (Temp <= VapEndTemp)) then
      VapEndDens =  VapTransDens*(Pres/RefPres)
      LiqEndDens = LiqTransDens + &
        LiqCompress*(Pres-RefPres)
      Slope = (VapEndDens-LiqEndDens)/ &
              (VapEndTemp-LiqEndTemp)
      Answer = Slope*(Temp-LiqEndTemp)+LiqEndDens
      TempDeriv = Slope
      DSlopeP = (1./(VapEndTemp-LiqEndTemp))* &
          (VapEndDens/RefPres-LiqCompress)
      PresDeriv = DSlopeP*(Temp-LiqEndTemp)+ LiqCompress    
    else
!   Perfect gas law (temperature only)    
      Answer = VapTransDens*(VapEndTemp/Temp)*(Pres/RefPres)               
      TempDeriv = -Answer/Temp
      PresDeriv = Answer/Pres
    end if
  elseif (Select==1) then
    !-----------------CALCULATE ENTHALPY-------------------
     If (Temp < 0.0) then
      Answer = LiqRefEnth + &
               LiqDEnthDTemp*(-LiqRefTemp)
       TempDeriv = 0.0
       PresDeriv = 0.0
    elseif (Temp < LiqEndTemp) then
      Answer = LiqRefEnth + &
               LiqDEnthDTemp*(Temp-LiqRefTemp)
      TempDeriv = LiqDEnthDTemp
      PresDeriv = 0.0
    elseif ((Temp >= LiqEndTemp).and. &
            (Temp <= VapEndTemp)) then
      VapEndEnth = VapTransEnth
      LiqEndEnth = LiqTransEnth
      Slope = (VapEndEnth-LiqEndEnth)/ &
              (VapEndTemp-LiqEndTemp)
      Answer = LiqEndEnth + &
               Slope*(Temp-LiqEndTemp)
      TempDeriv = Slope
      PresDeriv = 0.0
    else
      Answer = VapRefEnth + &
               VapDEnthDTemp*(Temp-VapRefTemp)
      TempDeriv = VapDEnthDTemp
      PresDeriv = 0.0
    end if
  elseif (Select==2) then
    !-----------------CALCULATE VISCOSITY-----------------
     if (Temp < 1.) then
      Answer = LiqViscTerm1  + LiqViscTerm2           
      TempDeriv = 0.0
      PresDeriv = 0.0
    elseif (Temp < LiqEndTemp) then
      Answer = LiqViscTerm1  + LiqViscTerm2/Temp            
      TempDeriv = - LiqViscTerm2/Temp**2
      PresDeriv = 0.0
    elseif ((Temp >= LiqEndTemp).and. &
            (Temp <= VapEndTemp)) then
      VapEndVisc = VapRefVisc
      LiqEndVisc = LiqRefVisc
      Slope = (VapEndVisc-LiqEndVisc)/ &
              (VapEndTemp-LiqEndTemp)
      Answer = Slope*(Temp-LiqEndTemp)+LiqEndVisc
      TempDeriv = Slope
      if (Temp==VapEndTemp) then
        Answer = VapEndVisc
      end if
      PresDeriv = 0.0
    else
      Answer = VapRefVisc + VapDViscDTemp*(Temp-VapRefTemp)
      TempDeriv = VapDViscDTemp
      PresDeriv = 0.0
    end if
  else
    ! Return error
    write(*,*) 'Selection code not set properly --WaterProps'
    write(*,*) 'Selection code: ',Select
    STOP
  end if
end if
end
