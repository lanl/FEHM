!     Last change:  JD   12 Oct 2006    4:56 am
module property_interpolate
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


  private

  public   read_interpolation_data,              &
       get_error_message,                    &
       get_property_type,                    &
       get_co2_state,                        &
       get_co2_sat_pressure,                 &
       get_co2_sat_temperature,              &
       get_co2_properties,                   &
       get_co2_sat_properties_pressure,      &
       get_co2_sat_props_temperature,        &
       get_co2_sat_line_props_pressure,      &
       get_co2_sat_line_props_temperat,      &
       interpolation_arrays_deallocate

  ! Reduced name length (<31 max)
  !get_co2_sat_properties_temperature,   &
  !get_co2_sat_line_props_temperature,   &


  ! -- Global variables

  character*1                :: at                       ! "u" if inform grid; "n" if nonuniform

  integer                    :: nt                       ! Number of temperatures in array
  integer                    :: np                       ! Number of pressures in array
  integer                    :: na                       ! Number of arrays
  real*8, allocatable          :: t(:)                     ! Temperatures at array points
  real*8, allocatable          :: p(:)                     ! Pressures at array points

  real*8                       :: t_index_offset           ! Temperature offset
  real*8                       :: t_index_factor           ! Temperature index factor
  real*8                       :: p_index_offset           ! Pressure offset
  real*8                       :: p_index_factor           ! Pressure offset factor
  real*8                       :: t_table_min              ! Minimum of scaled temperatures
  real*8                       :: t_table_max              ! Maximum of scaled temperatures
  real*8                       :: p_table_min              ! Minimum of scaled pressures
  real*8                       :: p_table_max              ! Maximum of scaled pressures

  character*20, allocatable  :: property_type(:)         ! Property type in each array
  real*8, allocatable          :: rarray(:,:,:)            ! Array of properties
  logical, allocatable       :: satline(:,:)             ! Cells intersected by saturation line

  integer                    :: nsat                     ! Number of saturation line vertices
  real*8                       :: tsat_min                 ! Minimum scaled temperature along sat. line
  real*8                       :: tsat_max                 ! Maximum scaled temperature along sat. line
  real*8                       :: psat_min                 ! Minimum scaled pressure along sat. line
  real*8                       :: psat_max                 ! Maximum scaled pressure along sat. line
  real*8, allocatable          :: tsat(:)                  ! Scaled temperatures along saturation line
  real*8, allocatable          :: psat(:)                  ! Scaled pressures along saturation line
  integer, allocatable       :: csat(:)                  ! Array cell number
  integer, allocatable       :: ssat(:)                  ! Type of intersection
  real*8, allocatable          :: msat(:)                  ! Scaled slope of p vs t along saturation line segment
  real*8, allocatable          :: lpsat(:,:)               ! Liquid properties along sat. line for all arrays
  real*8, allocatable          :: gpsat(:,:)               ! Vapour properties along sat. line for all arrays
  integer                    :: sat_index_last_tt=0      ! Used to expedite array searching
  integer                    :: sat_index_last_pp=0      ! Used to expedite array searching
  integer                    :: sat_index_last_ic=0      ! Used to expedite array searching
  integer                    :: ic_last=0                ! Used to expedite array searching
  integer                    :: it_last=0                ! Used to expedite array searching
  integer                    :: ip_last=0                ! Used to expedite array searching

  integer                    :: isatclose                ! Times saturation line is numerically close to grid node
  integer, allocatable       :: satclose(:,:)            ! Saturation line closeness index array

  character*500              :: amessage=' '             ! Error message text

contains


  subroutine get_error_message(astring)

    ! -- Subroutine get_error_message gets an error message.

    implicit none
    character*(*), intent(out) :: astring

    integer              :: l1,l2,l

    l1=len(astring)
    l2=len_trim(amessage)
    if(l2.eq.0)then
       astring=' '
    else
       l=min(l1,l2)
       astring=amessage(1:l)
    end if
    return
  end subroutine get_error_message




  subroutine read_interpolation_data(ifail,infile)

    ! --  Subroutine read_interpolation_data reads an interpolation dataset.

    implicit none
    integer, intent(out)       :: ifail
    character*(*), intent(in)  :: infile

    integer                    :: iunit,ierr,it,ip,ia,isat
    real*8                       :: rtemp
    character*20               :: atemp
    character*200              :: afile

    ifail=0

    iunit=nextunit()
    call addquote(infile,afile)
    open(unit=iunit,file=infile,status='old',iostat=ierr)
    if(ierr.ne.0)then
       write(amessage,20) trim(afile)
20     format('Cannot open file ',a,' to read interpolation data.')
       go to 9890
    end if

    ! -- The grid type is read.

    read(iunit,*) atemp
    call lowcase(atemp)
    if(atemp.eq.'uniform')then
       at='u'
    else if(atemp.eq.'nonuniform')then
       at='n'
    else
       write(amessage,22) trim(afile)
22     format('First line of interpolation data file ',a,' should be "uniform" or "nonuniform".')
       go to 9890
    end if

    ! -- Array table dimensions are read.

    read(iunit,*,iostat=ierr) nt,np,na
    if(ierr.ne.0)then
       write(amessage,40) trim(afile)
40     format('Error reading dimensional information from second line of interpolation ',  &
            'data file ',a,'.')
       go to 9890
    end if
    if((nt.le.0).or.(np.le.0).or.(na.le.0))then
       write(amessage,50) trim(afile)
50     format('Illegal values for one or more dimensions on second line of interpolation ', &
            'data file ',a,'.')
       go to 9890
    end if

    ! -- Temperature and pressure index factors and offsets are read.

    read(iunit,*,iostat=ierr) t_index_factor, t_index_offset
    if(ierr.ne.0)then
       write(amessage,60) trim(afile)
60     format('Error reading temperature factor and/or offset from third line of interpolation ',  &
            'data file ',a,'.')
       go to 9890
    end if
    if(t_index_factor.le.0.0)then
       write(amessage,70) trim(afile)
70     format('Illegal value for temperature factor on third line of interpolation data file ',a,'.')
       go to 9890
    end if

    read(iunit,*,iostat=ierr) p_index_factor, p_index_offset
    if(ierr.ne.0)then
       write(amessage,80) trim(afile)
80     format('Error reading pressure factor and/or offset from fourth line of interpolation ',  &
            'data file ',a,'.')
       go to 9890
    end if
    if(p_index_factor.le.0.0)then
       write(amessage,90) trim(afile)
90     format('Illegal value for pressure factor on fourth line of interpolation data file ',a,'.')
       go to 9890
    end if

    ! -- The saturation line closeness flag is read.

    read(iunit,*,iostat=ierr) isatclose
    if(ierr.ne.0) then
       write(amessage,92) trim(afile)
92     format('Error reading saturation line closeness index from fifth line of interpolation ', &
            'data file ',a,'.')
       go to 9890
    end if
    if(isatclose.lt.0)then
       write(amessage,93) trim(afile)
93     format('Illegal value for saturation line closeness index on fifth line of ',  &
            'interpolation data file ',a,'.')
       go to 9890
    end if

    ! -- The temperature and pressure vectors are read.

    allocate(t(nt),p(np),stat=ierr)
    if(ierr.ne.0) go to 9400
    read(iunit,*,err=9200,end=9250)
    read(iunit,*,err=9200,end=9250) (t(it),it=1,nt)
    read(iunit,*,err=9300,end=9350)
    read(iunit,*,err=9300,end=9350) (p(ip),ip=1,np)

    ! -- The property type held within each array is now read.

    allocate(property_type(na),stat=ierr)
    if(ierr.ne.0) go to 9400
    read(iunit,*,iostat=ierr)
    if(ierr.ne.0)then
       write(amessage,95) trim(afile)
       go to 9890
    end if
    do ia=1,na
       read(iunit,'(a)',iostat=ierr) property_type(ia)
       if(ierr.ne.0)then
          write(amessage,95) trim(afile)
95        format('Error reading property type names from interpolation data file ',a,'.')
          go to 9890
       end if
       property_type(ia)=adjustl(property_type(ia))
       call lowcase(property_type(ia))
    end do

    ! -- The arrays are read.

    allocate(rarray(nt,np,na),satline(nt,np),stat=ierr)
    if(ierr.ne.0) go to 9400
    if(isatclose.gt.0)then
       allocate(satclose(nt,np),stat=ierr)
       if(ierr.ne.0) go to 9400
    end if

    do ia=1,na
       read(iunit,*,err=9100,end=9150)
       do ip=1,np
          read(iunit,*,err=9100,end=9150) (rarray(it,ip,ia),it=1,nt)
       end do
    end do
    read(iunit,*,err=9120,end=9170)
    do ip=1,np
       read(iunit,*,err=9120,end=9170) (satline(it,ip),it=1,nt)
    end do
    if(isatclose.gt.0)then
       read(iunit,*,err=9050,end=9070)
       do ip=1,np
          read(iunit,*,err=9050,end=9070) (satclose(it,ip),it=1,nt)
       end do
    end if

    ! -- Now we read information pertaining to intersections of the saturation line with the table.
    ! -- First the dimension of the intersection table.

    read(iunit,*,iostat=ierr)
    if(ierr.ne.0)then
       write(amessage,97) trim(afile)
       go to 9890
    end if
    read(iunit,*,iostat=ierr) nsat
    if(ierr.ne.0)then
       write(amessage,97) trim(afile)
97     format('Error reading number of saturation line vertices from interpolation data file ',a,'.')
       go to 9890
    end if
    if(nsat.le.0)then
       write(amessage,120) trim(afile)
120    format('Number of saturation line vertices supplied as zero or less in interpolation ',   &
            'data file ',a,'.')
       go to 9890
    end if

    ! -- Saturation Data is read.

    allocate(tsat(nsat),psat(nsat),stat=ierr)
    if(ierr.ne.0) go to 9400
    allocate(csat(nsat),ssat(nsat),msat(nsat),stat=ierr)
    if(ierr.ne.0) go to 9400
    read(iunit,*,err=9450,end=9450)
    do isat=1,nsat
       read(iunit,*,err=9450,end=9450) psat(isat),tsat(isat),msat(isat),csat(isat),ssat(isat)
    end do

    ! -- The extremes are evaluated.

    tsat_min=tsat(1)
    tsat_max=tsat(nsat)
    psat_min=psat(1)
    psat_max=psat(nsat)

    ! -- Liquid properties along the saturation line are now read.

    allocate(lpsat(nsat,na),gpsat(nsat,na),stat=ierr)
    if(ierr.ne.0) go to 9400
    atemp='liquid properties'
    read(iunit,*,err=9500,end=9500)
    read(iunit,*,err=9500,end=9500)
    read(iunit,*,err=9500,end=9500)
    do isat=1,nsat
       read(iunit,*,err=9500,end=9500) (lpsat(isat,ia),ia=1,na)
    end do

    atemp='vapour properties'
    read(iunit,*,err=9500,end=9500)
    read(iunit,*,err=9500,end=9500)
    read(iunit,*,err=9500,end=9500)
    do isat=1,nsat
       read(iunit,*,err=9500,end=9500) (gpsat(isat,ia),ia=1,na)
    end do

    ! -- The coordinates of intersection of the saturation line are now scaled for the uniform case.

    if(at.eq.'u')then
       do isat=1,nsat
          tsat(isat)=(tsat(isat)-t_index_offset)*t_index_factor
       end do
       do isat=1,nsat
          psat(isat)=(psat(isat)-p_index_offset)*p_index_factor
       end do

       ! -- Slopes of saturation line segments are now scaled.

       rtemp=p_index_factor/t_index_factor
       do isat=1,nsat
          msat(isat)=msat(isat)*rtemp
       end do

       ! -- Table coordinates are now scaled for the uniform case.

       do it=1,nt
          t(it)=(t(it)-t_index_offset)*t_index_factor
       end do
       do ip=1,np
          p(ip)=(p(ip)-p_index_offset)*p_index_factor
       end do

    end if

    ! -- The scaled saturation line limits are now calculated.

    tsat_min=tsat(1)
    tsat_max=tsat(nsat)
    psat_min=psat(1)
    psat_max=psat(nsat)

    ! -- The scaled table limits are now calculated

    t_table_min=t(1)
    t_table_max=t(nt)
    p_table_min=p(1)
    p_table_max=p(np)

    close(unit=iunit)

    return


9050 write(amessage,9060) trim(afile)
9060 format('Error reading saturation line closeness array from interpolation data ', &
         'file ',a,'.')
    go to 9890
9070 write(amessage,9080) trim(afile)
9080 format('Premature end encountered to interpolation data file ',a,    &
         ' while reading saturation line closeness array.')
    go to 9890
9100 write(amessage,9010) trim(property_type(ia)),trim(afile)
9010 format('Error reading ',a,' array from interpolation data file ',a,'.')
    go to 9890
9120 write(amessage,9130) trim(afile)
9130 format('Error reading saturation line intersection array from interpolation data ', &
         'file ',a,'.')
    go to 9890
9150 write(amessage,9160) trim(afile),trim(property_type(ia))
9160 format('Premature end encountered to interpolation data file ',a,    &
         ' while reading ',a,' array.')
    go to 9890
9170 write(amessage,9180) trim(afile)
9180 format('Premature end encountered to interpolation data file ',a,    &
         ' while reading saturation line intersection array.')
    go to 9890
9200 write(amessage,9210) trim(afile)
9210 format('Error reading table temperatures from interpolation data file ',a,'.')
    go to 9890
9250 write(amessage,9260) trim(afile)
9260 format('Premature end to file ',a,' encountered while reading table temperatures.')
    go to 9890
9300 write(amessage,9310) trim(afile)
9310 format('Error reading table pressures from interpolation array file ',a,'.')
    go to 9890
9350 write(amessage,9360) trim(afile)
9360 format('Premature end to file ',a,' encountered while reading table pressures.')
    go to 9890
9400 write(amessage,9410)
9410 format('Error in allocating memory for CO2 interpolation data arrays.')
    go to 9890
9450 write(amessage,9460) trim(afile)
9460 format('Error reading saturation line data from interpolation ',  &
         'data file ',a,'.')
    go to 9890
9500 write(amessage,9510) trim(atemp),trim(afile)
9510 format('Error reading saturation line ',a,' from interpolation data file ',a,'.')
    go to 9890



9890 ifail=1
    close(unit=iunit,iostat=ierr)
    return


  end subroutine read_interpolation_data



  subroutine get_property_type(ifail,ind,astring)

    ! -- Subroutine GET_PROPERTY_TYPE retrieves the property type pertaining to an array index.

    implicit none

    integer, intent(out)        :: ifail
    integer, intent(in)         :: ind
    character*(*), intent(out)  :: astring
    
    integer                     :: l1,l2,l

    ifail=0
    if((ind.lt.1).or.(ind.gt.na))then
       write(amessage,10)
10     format('Error in subroutine GET_PROPERTY_TYPE: supplied index out of range.')
       ifail=1
       return
    end if

    l1=len(astring)
    l2=len_trim(property_type(ind))
    l=min(l1,l2)
    astring=property_type(ind)(1:l)

    return
  end subroutine get_property_type



  subroutine get_co2_state(ifail,temperature,pressure,state)

    ! -- Subroutine GET_CO2_STATE returns the state of CO2 at the requested
    !    (unscaled) temperature and pressure.

    ! -- States are as follows:-
    !      1 - liquid
    !      2 - vapour
    !      3 - supercritical
    !      4 - exactly on saturation line

    implicit none
    integer, intent(out)  :: ifail
    real*8, intent(in)      :: temperature, pressure
    integer, intent(out)  :: state

    real*8 tt,pp

    ifail=0
    tt=(temperature-t_index_offset)*t_index_factor
    pp=(pressure-p_index_offset)*p_index_factor

    call get_state(ifail,tt,pp,state)

    return

  end subroutine get_co2_state



  subroutine get_state(ifail,tt,pp,state)

    ! -- Subroutine GET_STATE returns the state of CO2 at the requested
    !    scaled temperature and pressure.

    ! -- States are as follows:-
    !      1 - liquid (new state: liquid, old state liquid)
    !      2 - vapour (new state: 2-phase, old state: vapor)
    !      3 - supercritical (new state: vapor, old state: supercritical)
    !      4 - exactly on saturation line (new state: supercritical, old state: 2-phase)

    implicit none
    integer, intent(out) :: ifail
    real*8, intent(in)     :: tt,pp
    integer, intent(out) :: state

    integer isat
    real*8 p_ref

    ! -- Easy selections are made first.

    ifail=0

    if(tt.gt.tsat_max)then
       if(pp.lt.psat_max)then
          state=3
          return
       else
          state=4
          return
       end if
    else
       if(pp.gt.psat_max)then
          state=1
          return
       end if
    end if

    if((pp.lt.p_table_min).or.(tt.lt.t_table_min))then
       write(amessage,20)
20     format('Error in subroutine GET_STATE. Temperature and/or pressure out of range of ',  &
            'interpolation table.')
       go to 9890
    end if
    ! -- Note that the above test assumes that the table extends to higher temperatures and pressures
    !    than those prevailing at the end of the saturation line.

    ! -- We now determine what column of the table we are in.

    call which_index(tt,nsat,tsat,isat,sat_index_last_tt)

    if (isat .eq. 0) then
       ! Failed to find valid table index
       write(amessage,20)
       go to 9890
    end if

    p_ref=psat(isat)+msat(isat)*(tt-tsat(isat))
    if(pp.gt.p_ref)then
       state=1
    else if(pp.lt.p_ref)then
       state=3
    else
       state=2
    end if

    return

9890  ifail=1
      return

  end subroutine get_state



  subroutine which_sat_index_array(it,ip,isat)
    
    ! -- Subroutine WHICH_SAT_INDEX_ARRAY finds the index of the saturation-line cells
    !    pertaining to a particular cell in the interpolation table. It is assumed that
    !    the saturation line actually passes through this cell. So this is not checked
    !    in order to speed execution time.
    
    implicit none
    integer, intent(in)   :: it,ip
    integer, intent(out)  :: isat

    integer               :: ic,i

    call tp2cell(it,ip,ic)
    if(sat_index_last_ic.eq.0) then
       sat_index_last_ic=1
    end if

    if(ic.eq.ic_last)then
       isat=sat_index_last_ic
       return
    end if
    if(sat_index_last_ic.ne.1)then
       if(ic.eq.csat(sat_index_last_ic-1))then
          isat=sat_index_last_ic-1
          go to 200
       end if
    else if(sat_index_last_ic.ne.nsat)then
       if(ic.eq.csat(sat_index_last_ic+1))then
          isat=sat_index_last_ic+1
          go to 200
       end if
    end if

100 continue
    do i=sat_index_last_ic,nsat-1    ! hopefully the "-1" is ok here.
       if(ic.eq.csat(i))then
          isat=i
          go to 200
       end if
    end do
    do i=sat_index_last_ic-1,1,-1
       if(ic.eq.csat(i))then
          isat=i
          go to 200
       end if
    end do

    write(*,*) 'Error in subroutine WHICH_SAT_INDEX_ARRAY.'
    write(*,*) 'it = ',it
    write(*,*) 'ip = ',ip
    write(*,*) 'ic = ',ic
    stop

200 ic_last=ic
    sat_index_last_ic=isat
    return

  end subroutine which_sat_index_array




  subroutine which_index(rr,nr,r,ind,ind_last)

    ! -- Subroutine WHICH_INDEX calculates the segment of an array of numbers in which a current number lies.
    !    For speed of execution, it is assumed that rr is within range of the search domain (i.e. it is assumed
    !    that this has been previously tested.)

    implicit none
    real*8, intent(in)          :: rr
    integer, intent(in)       :: nr
    real*8, intent(in)          :: r(nr)
    integer, intent(out)      :: ind
    integer, intent(inout)    :: ind_last

    integer                   :: i

    if((ind_last.le.0).or.(ind_last.ge.nr)) then
       ind_last=1
    end if

    if(rr.ge.r(ind_last))then
       do i=ind_last+1,nr
          if(rr.le.r(i))then
             ind=i-1
             ind_last=ind
             return
          end if
       end do
    else
       do i=ind_last-1,1,-1
          if(rr.ge.r(i))then
             ind=i
             ind_last=ind
             return
          end if
       end do
    end if

    ind = 0

  end subroutine which_index



  subroutine get_co2_sat_pressure(ifail,temperature,pressure,dp_dt)

    ! -- Subroutine GET_CO2_SAT_PRESSURE gets the saturation pressure corresponding to a certain
    !    temperature.

    ! -- Note carefully. The return value of ifail is -1 if temperature is above end of saturation line.

    ! -- Question - will FEHM be supplying variables like temperature in double precision.

    implicit none
    integer, intent(out)    :: ifail
    real*8, intent(in)        :: temperature
    real*8, intent(out)       :: pressure
    real*8, intent(out)       :: dp_dt

    real*8                    :: tt,pp,mm

    ifail=0
    tt=(temperature-t_index_offset)*t_index_factor
    call get_sat_pressure(ifail,tt,pp,mm)
    if(ifail.ne.0) return
    pressure=pp/p_index_factor+p_index_offset
    dp_dt=mm*t_index_factor/p_index_factor

    return

  end subroutine get_co2_sat_pressure



  subroutine get_sat_pressure(ifail,tt,pp,mm)

    ! -- Subroutine GET_SAT_PRESSURE gets the scaled saturation pressure corresponding to a certain
    !    scaled temperature.

    implicit none
    integer, intent(out)      :: ifail
    real*8, intent(in)          :: tt
    real*8, intent(out)         :: pp
    real*8, intent(out)         :: mm

    integer                   :: isat

    ifail=0
    if(tt.gt.tsat_max)then
       ifail=-1              ! No error message as this is likely to be a common occurrence.
       return
    else if(tt.lt.tsat_min)then
       ifail=1
       write(amessage,10)
10     format('Error in subroutine GET_SAT_PRESSURE: supplied temperature is out of ',  &
            'range of interpolation table.')
       return
    else
       call which_index(tt,nsat,tsat,isat,sat_index_last_tt)
       pp=psat(isat)+msat(isat)*(tt-tsat(isat))
       mm=msat(isat)
    end if

    return

  end subroutine get_sat_pressure



  subroutine get_co2_sat_temperature(ifail,pressure,temperature,dt_dp)

    ! -- Subroutine GET_CO2_SAT_TEMPERATURE gets the saturation temperature corresponding to a certain
    !    pressure.

    ! -- Note carefully. The return value of ifail is -1 if pressure is above end of saturation line.

    ! -- Question - will FEHM be supplying variables like temperature in double precision.

    implicit none
    integer, intent(out)    :: ifail
    real*8, intent(in)        :: pressure
    real*8, intent(out)       :: temperature
    real*8, intent(out)       :: dt_dp

    real*8                    :: pp,tt,mm

    ifail=0
    pp=(pressure-p_index_offset)*p_index_factor
    call get_sat_temperature(ifail,pp,tt,mm)
    if(ifail.ne.0) return
    temperature=tt/t_index_factor+t_index_offset
    dt_dp=p_index_factor/mm/t_index_factor

    return

  end subroutine get_co2_sat_temperature



  subroutine get_sat_temperature(ifail,pp,tt,mm)

    ! -- Subroutine GET_SAT_TEMPERATURE gets the scaled saturation temperature corresponding to a certain
    !    scaled pressure.

    implicit none
    integer, intent(out)      :: ifail
    real*8, intent(in)          :: pp
    real*8, intent(out)         :: tt
    real*8, intent(out)         :: mm

    integer                   :: isat

    ifail=0
    if(pp.gt.psat_max)then
       ifail=-1              ! No error message as this is likely to be a common occurrence.
       return
    else if(pp.lt.psat_min)then
       ifail=1
       write(amessage,10)
10     format('Error in subroutine GET_SAT_TEMPERATURE: supplied pressure is out of ',  &
            'range of interpolation table.')
       return
    else
       call which_index(pp,nsat,psat,isat,sat_index_last_pp)
       tt=tsat(isat)+(pp-psat(isat))/msat(isat)
       mm=msat(isat)
    end if
    return

  end subroutine get_sat_temperature



  subroutine rectangle_interpolation_factors(x1,x2,y1,y2,xp,yp,fac)

    ! -- Subroutine RECTANGULAR_INTERPOLATION_FACTORS computes factors for rectangular interpolation.

    implicit none
    real*8, intent(in)              :: x1,x2,y1,y2
    real*8, intent(in)              :: xp,yp
    real*8, intent(out)             :: fac(4)

    real*8                          :: y2_yp,x2_xp,xp_x1,yp_y1,den

    y2_yp=y2-yp
    x2_xp=x2-xp
    xp_x1=xp-x1
    yp_y1=yp-y1
    den=1.0/((x2-x1)*(y2-y1))

    fac(1)=y2_yp*x2_xp*den
    fac(2)=y2_yp*xp_x1*den
    fac(3)=yp_y1*x2_xp*den
    fac(4)=yp_y1*xp_x1*den


    return

  end subroutine rectangle_interpolation_factors



  subroutine quad_interpolate(fac,z1,z2,z3,z4,rval)

    ! -- Subroutine QUAD_INTERPOLATE performs interpolation where there are four interpolation factors.

    implicit none
    real*8, intent(in)       :: fac(4)
    real*8, intent(in)       :: z1,z2,z3,z4
    real*8, intent(out)      :: rval

    rval=fac(1)*z1+fac(2)*z2+fac(3)*z3+fac(4)*z4
    return

  end subroutine quad_interpolate




  subroutine triangle_interpolation_factors (x1,x2,x3,y1,y2,y3,xp,yp,fac)

    ! -- Subroutine TRIANGE_INTERPOLATION_FACTORS computes interpolation factors for a triangle.

    implicit none
    real*8, intent(in)     :: x1,x2,x3,y1,y2,y3
    real*8, intent(in)     :: xp,yp
    real*8, intent(out)    :: fac(3)

    real*8                 :: det,x1_xp,x2_xp,x3_xp,y1_yp,y2_yp,y3_yp

    det=1.0/((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
    x1_xp=x1-xp
    x2_xp=x2-xp
    x3_xp=x3-xp
    y1_yp=y1-yp
    y2_yp=y2-yp
    y3_yp=y3-yp
    fac(1)=(x2_xp*y3_yp-x3_xp*y2_yp)*det
    fac(2)=(x3_xp*y1_yp-x1_xp*y3_yp)*det
    fac(3)=(x1_xp*y2_yp-x2_xp*y1_yp)*det
    return

  end subroutine triangle_interpolation_factors



  subroutine triangle_interpolate(fac,z1,z2,z3,rval)

    ! -- Subroutine TRIANGLE_INTERPOLATE performs triangle interpolation on the basis of pre-calculated
    !    interpolation factors.

    implicit none
    real*8, intent(in)       :: fac(3)
    real*8, intent(in)       :: z1,z2,z3
    real*8, intent(out)      :: rval


    rval=fac(1)*z1+fac(2)*z2+fac(3)*z3
    return

  end subroutine triangle_interpolate



  subroutine linear_interpolation_factors(x1,x2,xp,fac)

    ! -- Subroutine LINEAR_INTERPOLATION_FACTORS computes interpolation factors along a line.

    implicit none
    real*8, intent(in)     :: x1,x2
    real*8, intent(in)     :: xp
    real*8, intent(out)    :: fac(2)

    real*8                 :: x2_x1

    x2_x1=x2-x1
    fac(1)=(x2-xp)/x2_x1
    fac(2)=(xp-x1)/x2_x1

    return

  end subroutine linear_interpolation_factors

  !parallelepiped_interpolation_factors
  subroutine parallelepiped_interp_factors(itype,jtype,x1,x2,x3,x4,y1,y2,y3,y4,xp,yp,fac)

    ! -- Subroutine PARALLELEPIPED_INTERPOLATION_FACTORS computes interpolation factors for a pararallelepiped
    !    in which the parallel lines are horizontal (itype=2) or vertical (itype=4). The jtype variable
    !    indicates which of the other sides is non-vertical.

    implicit none
    integer, intent(in) :: itype,jtype
    real*8, intent(in)    :: x1,x2,x3,x4,y1,y2,y3,y4
    real*8, intent(in)    :: xp,yp
    real*8, intent(out)   :: fac(4)

    real*8                :: x5_x1,x2_x1,x4_x1,x2_x5,x2_x3,xp_x1,x2_xp
    real*8                :: w1,w2,x5,y5
    real*8                :: l1,l2,l3,l4,l5,l6
    real*8                :: yp_y1,y3_yp,y3_y1,y3_y2,y3_y5,y4_y1,y5_y1
    real*8                :: a


    if(itype.eq.2)then
       if(jtype.eq.1)then
          x5=x2+(yp-y1)/(y3-y1)*(x4-x2)
          x5_x1=x5-x1
          x2_x1=x2-x1
          x4_x1=x4-x1
          w1=(xp-x1)/x5_x1
          w2=(x5-xp)/x5_x1
          l1=w1*x2_x1
          l2=w2*x2_x1
          l3=w1*x4_x1
          l4=w2*x4_x1
          l5=w1*x5_x1
          l6=w2*x5_x1
       else
          x5=x1+(yp-y1)/(y3-y1)*(x3-x1)
          x2_x5=x2-x5
          x2_x1=x2-x1
          x2_x3=x2-x3
          w1=(xp-x5)/x2_x5
          w2=(x2-xp)/x2_x5
          l1=w1*x2_x1
          l2=w2*x2_x1
          l3=w1*x2_x3
          l4=w2*x2_x3
          l5=w1*x2_x5
          l6=w2*x2_x5
       end if
       yp_y1=yp-y1
       y3_yp=y3-yp
       y3_y1=y3-y1
       a=1.0/(l1+l2+l3+l4)/y3_y1
       fac(1)=(l4+l6)*y3_yp*a
       fac(2)=(l3+l5)*y3_yp*a
       fac(3)=(l2+l6)*yp_y1*a
       fac(4)=(l1+l5)*yp_y1*a
    else if(itype.eq.4)then
       if(jtype.eq.1)then
          y5=y1+(xp-x1)/(x2-x1)*(y2-y1)
          y3_y1=y3-y1
          y3_y2=y3-y2
          y3_y5=y3-y5
          w1=(yp-y5)/y3_y5
          w2=(y3-yp)/y3_y5
          l1=w1*y3_y1
          l2=w2*y3_y1
          l3=w1*y3_y2
          l4=w2*y3_y2
          l5=w1*y3_y5
          l6=w2*y3_y5
       else
          y5=y3+(xp-x1)/(x2-x1)*(y4-y3)
          y3_y1=y3-y1
          y4_y1=y4-y1
          y5_y1=y5-y1
          w1=(yp-y1)/y5_y1
          w2=(y5-yp)/y5_y1
          l1=w1*y3_y1
          l2=w2*y3_y1
          l3=w1*y4_y1
          l4=w2*y4_y1
          l5=w1*y5_y1
          l6=w2*y5_y1
       end if
       xp_x1=xp-x1
       x2_xp=x2-xp
       x2_x1=x2-x1
       a=1.0/(l1+l2+l3+l4)/(x2_x1)
       fac(1)=(l4+l6)*x2_xp*a
       fac(2)=(l2+l6)*xp_x1*a
       fac(3)=(l3+l5)*x2_xp*a
       fac(4)=(l1+l5)*xp_x1*a
    end if

    return

  end subroutine parallelepiped_interp_factors
  ! parallelepiped_interpolation_factors



  subroutine project_to_saturation_line(itype,xc,xa,xb,yc,ya,yb,xp,yp,xi,yi)

    ! -- Subroutine PROJECT_TO_SATURATION_LINE is used in modified rectangular interplation for
    !    cells in which type 1 or type 3 intersection of the saturation line occurs. It calculates
    !    the projection of interpolation points into the saturation line.

    implicit none
    integer, intent(in)     :: itype
    real*8, intent(in)        :: xc,xa,xb,yc,ya,yb
    real*8, intent(in)        :: xp,yp
    real*8, intent(out)       :: xi,yi

    real*8                    :: xc_xp,mp,mi

    xc_xp=xc-xp
    if(abs(xc_xp).lt.1e-4*(xc+xp))then    ! arbitrary
       xi=xc
       if(itype.eq.1)then
          yi=ya
       else
          yi=yb
       end if
       return
    end if
    mp=(yc-yp)/(xc-xp)
    mi=(yb-ya)/(xb-xa)
    xi=(mp*xp-mi*xa-(yp-ya))/(mp-mi)
    yi=mi*(xi-xa)+ya

    return

  end subroutine project_to_saturation_line



  subroutine get_co2_properties(ifail,iphase,ncode,icode,temperature,pressure,value)

    ! -- Subroutine GET_CO2_PROPERTIES obtains CO2 properties through interpolation from a table.

    implicit none

    integer, intent(out)      :: ifail
    integer, intent(in)       :: iphase
    integer, intent(in)       :: ncode
    integer, intent(in)       :: icode(ncode)
    real*8,    intent(in)       :: temperature
    real*8,    intent(in)       :: pressure
    real*8,    intent(out)      :: value(ncode)

    integer                   :: it,ip,isat,itype,istart,i,ii,iit,iip
    real*8                      :: tt,pp,ti,pi,zi
    real*8                      :: z1,z2,z3,z4
    real*8                      :: fac(5),faci(4),facp(4),facl(2),rarr(4)

    ! -- We will not check that NCODE is no greater than NA - it would take too much time to do this
    !    on every occasion that this subroutine is called.

    ifail=0

    if(iphase.eq.4)then
       write(amessage,5)
5      format('Error calling subroutine TABLE_INTERPOLATION; this subroutine should ',   &
            'not be called with IPHASE equal to 4.')
       go to 9890
    end if
    tt=(temperature-t_index_offset)*t_index_factor
    pp=(pressure-p_index_offset)*p_index_factor
    if((tt.lt.t_table_min).or.(tt.gt.t_table_max).or.   &
         (pp.lt.p_table_min).or.(pp.gt.p_table_max))then
       write(amessage,10)
10     format('Temperature or pressure out of interpolation range in call to ',  &
            'subroutine GET_CO2_PROPERTIES.')
       go to 9890
    end if
    if(at.eq.'u')then
       it_last=int(tt)
       ip_last=int(pp)
    end if
    call which_index(tt,nt,t,it,it_last)
    call which_index(pp,np,p,ip,ip_last)

    do i=1,ncode
       if(icode(i).ne.0)then
          istart=i
          go to 50
       end if
    end do
    write(amessage,40)
40  format('Ilegal call to subroutine GET_CO2_PROPERTIES. No interpolation ',  &
         'has been requested.')
    go to 9890
50  continue

    if(.not.satline(it,ip))then

       ! The saturation line does not intersect this cell.
       ! Note that we won't even check that the phase is correct. We will assume for the sake
       !      of speed that it is already so.

       call rectangle_interpolation_factors(t(it),t(it+1),p(ip),p(ip+1),tt,pp,fac)
       do i=istart,ncode
          if(icode(i).ne.0)then
             if(isatclose.eq.0)then
                call quad_interpolate(fac,rarray(it,ip,i),rarray(it+1,ip,i),rarray(it,ip+1,i),   &
                     rarray(it+1,ip+1,i),value(i))
             else
                ii=0
                do iip=ip,ip+1
                   do iit=it,it+1
                      ii=ii+1
                      isat=satclose(iit,iip)
                      if(isat.eq.0)then
                         rarr(ii)=rarray(iit,iip,i)
                      else
                         if(iphase.eq.1)then
                            rarr(ii)=lpsat(isat,i)
                         else
                            rarr(ii)=gpsat(isat,i)
                         end if
                      end if
                   end do
                end do
                call quad_interpolate(fac,rarr(1),rarr(2),rarr(3),rarr(4),value(i))
             end if
          end if
       end do
       return
    else
       call which_sat_index_array(it,ip,isat)
       itype=ssat(isat)

       if(iphase.eq.1)then
          if((itype.eq.1).or.(itype.eq.5).or.(itype.eq.6).or.(itype.eq.8))then
             call triangle_interpolation_factors (tsat(isat),t(it),  tsat(isat+1),   &
                  psat(isat),p(ip+1),psat(isat+1),tt,pp,fac)
             do i=istart,ncode
                call triangle_interpolate(fac,lpsat(isat,i),rarray(it,ip+1,i),lpsat(isat+1,i),value(i))
             end do
             return
          else if(itype.eq.3)then
             call project_to_saturation_line(3,t(it+1),tsat(isat),tsat(isat+1),   &
                  p(ip),psat(isat),psat(isat+1),tt,pp,ti,pi)
             call rectangle_interpolation_factors(t(it),t(it+1),p(ip),p(ip+1),ti,pi,faci)
             call rectangle_interpolation_factors(t(it),t(it+1),p(ip),p(ip+1),tt,pp,facp)
             call linear_interpolation_factors(tsat(isat),tsat(isat+1),ti,facl) ! Assumes sat line not too steep.
             do i=istart,ncode
                if(icode(i).ne.0)then
                   z1=rarray(it,ip,i)
                   z3=rarray(it,ip+1,i)
                   z4=rarray(it+1,ip+1,i)
                   if(abs(faci(2)).gt.1.0e-5)then   ! arbitrary
                      zi=facl(1)*lpsat(isat,i)+facl(2)*lpsat(isat+1,i)
                      z2=(zi-faci(1)*z1-faci(3)*z3-faci(4)*z4)/faci(2)
                   else
                      z2=0.5*(lpsat(isat,i)+lpsat(isat+1,i))
                   end if
                   call quad_interpolate(facp,z1,z2,z3,z4,value(i))
                end if
             end do
             return
          else if((itype.eq.2).or.(itype.eq.9))then
             ! parallelepiped_interpolation_factors
             call parallelepiped_interp_factors(2,1,t(it),tsat(isat),t(it),tsat(isat+1),  &
                  p(ip),psat(isat),p(ip+1),p(isat+1),tt,pp,fac)
             do i=istart,ncode
                if(icode(i).ne.0)then
                   call quad_interpolate(fac,rarray(it,ip,i),lpsat(isat,i),rarray(it,ip+1,i),lpsat(isat+1,i),value(i))
                end if
             end do
             return
          else if((itype.eq.4).or.(itype.eq.7))then
             call parallelepiped_interp_factors(4,1,tsat(isat),tsat(isat+1),t(it),t(it+1),  &
                  psat(isat),psat(isat+1),p(ip+1),p(ip+1),tt,pp,fac)
             do i=istart,ncode
                if(icode(i).ne.0)then
                   call quad_interpolate(fac,lpsat(isat,i),lpsat(isat+1,i),rarray(it,ip+1,i),    &
                        rarray(it+1,ip+1,i),value(i))
                end if
             end do
             return
          end if
       else
          if(itype.eq.1)then
             call project_to_saturation_line(1,t(it),tsat(isat),tsat(isat+1),   &
                  p(ip+1),psat(isat),psat(isat+1),tt,pp,ti,pi)
             call rectangle_interpolation_factors(t(it),t(it+1),p(ip),p(ip+1),ti,pi,faci)
             call rectangle_interpolation_factors(t(it),t(it+1),p(ip),p(ip+1),tt,pp,facp)
             call linear_interpolation_factors(tsat(isat),tsat(isat+1),ti,facl) ! assumes sat line not too steep
             do i=istart,ncode
                if(icode(i).ne.0)then
                   z1=rarray(it,ip,i)
                   z2=rarray(it+1,ip,i)
                   z4=rarray(it+1,ip+1,i)
                   if(abs(faci(3)).gt.1.0e-5)then
                      zi=facl(1)*gpsat(isat,i)+facl(2)*gpsat(isat+1,i)
                      z3=(zi-faci(1)*z1-faci(2)*z2-faci(4)*z4)/faci(3)
                   else
                      z3=0.5*(gpsat(isat,i)+gpsat(isat+1,i))
                   end if
                   call quad_interpolate(facp,z1,z2,z3,z4,value(i))
                end if
             end do
             return
          else if((itype.eq.3).or.(itype.eq.6).or.(itype.eq.7).or.(itype.eq.9))then
             call triangle_interpolation_factors (tsat(isat),tsat(isat+1),t(it+1),              &
                  psat(isat),psat(isat+1),p(ip),tt,pp,fac)
             do i=istart,ncode
                if(icode(i).ne.0)then
                   call triangle_interpolate(fac,gpsat(isat,i),gpsat(isat+1,i),rarray(it+1,ip,i),value(i))
                end if
             end do
             return
          else if((itype.eq.2).or.(itype.eq.5))then
             call parallelepiped_interp_factors(2,2,tsat(isat),t(it+1),tsat(isat+1),t(it+1),  &
                  psat(isat),p(ip),psat(isat+1),p(ip+1),tt,pp,fac)
             do i=istart,ncode
                if(icode(i).ne.0)then
                   call quad_interpolate(fac,gpsat(isat,i),rarray(it+1,ip,i),gpsat(isat+1,i),   &
                        rarray(it+1,ip+1,i),value(i))
                end if
             end do
             return
          else if((itype.eq.4).or.(itype.eq.8))then
             call parallelepiped_interp_factors(4,2,t(it),t(it+1),tsat(isat),tsat(isat+1),  &
                  p(ip),p(ip),psat(isat),psat(isat+1),tt,pp,fac)
             do i=istart,ncode
                if(icode(i).ne.0)then
                   call quad_interpolate(fac,rarray(it,ip,i),rarray(it+1,ip,i),gpsat(isat,i),   &
                        gpsat(isat+1,i),value(i))
                end if
             end do
             return
          end if
       end if
    end if

    return
9890 ifail=1
    return

  end subroutine get_co2_properties



  subroutine get_co2_sat_properties_pressure(ifail,iphase,ncode,icode,pressure,value)

    ! -- Subroutine GET_CO2_SAT_PROPERTIES_PRESSURE gets saturation properties pertaining to a certain
    !    pressure.

    implicit none

    integer, intent(out)      :: ifail
    integer, intent(in)       :: iphase
    integer, intent(in)       :: ncode
    integer, intent(in)       :: icode(ncode)
    real*8,    intent(in)       :: pressure
    real*8,    intent(out)      :: value(ncode)

    integer                   :: isat,i
    real*8                      :: pp
    real*8                      :: fac(2)

    ifail=0

    pp=(pressure-p_index_offset)*p_index_factor
    if(pp.gt.psat_max)then
       ifail=-1
       return
    else if(pp.lt.psat_min)then
       write(amessage,10)
10     format('Error in subroutine GET_CO2_SAT_PROPERTIES_PRESSURE: pressure out of table ', &
            'interpolation range.')
       ifail=1
       return
    else
       call which_index(pp,nsat,psat,isat,sat_index_last_pp)
       call linear_interpolation_factors(psat(isat),psat(isat+1),pp,fac)
       if(iphase.eq.1)then
          do i=1,ncode
             if(icode(i).ne.0)then
                value(i)=lpsat(isat,i)*fac(1)+lpsat(isat+1,i)*fac(2)
             end if
          end do
       else if(iphase.eq.2)then
          do i=1,ncode
             if(icode(i).ne.0)then
                value(i)=gpsat(isat,i)*fac(1)+gpsat(isat+1,i)*fac(2)
             end if
          end do
       else
          write(amessage,20)
20        format('Error in subroutine GET_CO2_SAT_PROPERTIES_PRESSURE: phase code must be "1" or "2".')
          ifail=1
          return
       end if
    end if

    return

  end subroutine get_co2_sat_properties_pressure


  !get_co2_sat_properties_temperature
  subroutine get_co2_sat_props_temperature(ifail,iphase,ncode,icode,temperature,value)

    ! -- Subroutine GET_CO2_SAT_PROPERTIES_TEMPERATURE gets saturation properties pertaining to a certain
    !    temperature.

    implicit none

    integer, intent(out)      :: ifail
    integer, intent(in)       :: iphase
    integer, intent(in)       :: ncode
    integer, intent(in)       :: icode(ncode)
    real*8,    intent(in)       :: temperature
    real*8,    intent(out)      :: value(ncode)

    integer                   :: isat,i
    real*8                      :: tt
    real*8                      :: fac(2)

    ifail=0

    tt=(temperature-t_index_offset)*t_index_factor
    if(tt.gt.tsat_max)then
       ifail=-1
       return
    else if(tt.lt.tsat_min)then
       write(amessage,10)
10     format('Error in subroutine GET_CO2_SAT_PROPERTIES_TEMPERATURE: temperature out of table ', &
            'interpolation range.')
       ifail=1
       return
    else
       call which_index(tt,nsat,tsat,isat,sat_index_last_tt)
       call linear_interpolation_factors(tsat(isat),tsat(isat+1),tt,fac)
       if(iphase.eq.1)then
          do i=1,ncode
             if(icode(i).ne.0)then
                value(i)=lpsat(isat,i)*fac(1)+lpsat(isat+1,i)*fac(2)
             end if
          end do
       else if(iphase.eq.2)then
          do i=1,ncode
             if(icode(i).ne.0)then
                value(i)=gpsat(isat,i)*fac(1)+gpsat(isat+1,i)*fac(2)
             end if
          end do
       else
          write(amessage,20)
20        format('Error in subroutine GET_CO2_SAT_PROPERTIES_TEMPERATURE: phase code must be "1" or "2".')
          ifail=1
          return
       end if
    end if

    return

  end subroutine GET_CO2_SAT_PROPS_TEMPERATURE
  ! get_co2_sat_properties_temperature



  subroutine get_co2_sat_line_props_pressure(ifail,iphase,ncode,icode,pressure,value)

    ! -- Subroutine GET_CO2_SAT_LINE_PROPS_PRESSURE gets saturation properties pertaining
    !    to a certain pressure. However derivatives are taken along the actual saturation line.

    implicit none

    integer, intent(out)      :: ifail
    integer, intent(in)       :: iphase
    integer, intent(in)       :: ncode
    integer, intent(in)       :: icode(ncode)
    real*8,    intent(in)       :: pressure
    real*8,    intent(out)      :: value(ncode)

    integer                   :: isat,i
    real*8                      :: pp
    real*8                      :: delta_p,delta_t
    real*8                      :: fac(2)

    ifail=0

    pp=(pressure-p_index_offset)*p_index_factor
    if(pp.gt.psat_max)then
       ifail=-1
       return
    else if(pp.lt.psat_min)then
       write(amessage,10)
10     format('Error in subroutine GET_CO2_SAT_LINE_PROPS_PRESSURE: ',   &
            'pressure out of table interpolation range.')
       ifail=1
       return
    else
       call which_index(pp,nsat,psat,isat,sat_index_last_pp)
       call linear_interpolation_factors(psat(isat),psat(isat+1),pp,fac)
       delta_p=psat(isat+1)-psat(isat)
       delta_t=tsat(isat+1)-tsat(isat)
       if(iphase.eq.1)then
          do i=1,ncode
             if(icode(i).ne.0)then
                if((i.eq.3).or.(i.eq.6).or.(i.eq.9))then
                   value(i)=(lpsat(isat+1,i-2)-lpsat(isat,i-2))/delta_p
                else if((i.eq.2).or.(i.eq.5).or.(i.eq.8))then
                   value(i)=(lpsat(isat+1,i-1)-lpsat(isat,i-1))/delta_t
                else
                   value(i)=lpsat(isat,i)*fac(1)+lpsat(isat+1,i)*fac(2)
                end if
             end if
          end do
       else if(iphase.eq.2)then
          do i=1,ncode
             if(icode(i).ne.0)then
                if((i.eq.3).or.(i.eq.6).or.(i.eq.9))then
                   value(i)=(gpsat(isat+1,i-2)-gpsat(isat,i-2))/delta_p
                else if((i.eq.2).or.(i.eq.5).or.(i.eq.8))then
                   value(i)=(gpsat(isat+1,i-1)-gpsat(isat,i-1))/delta_t
                else
                   value(i)=gpsat(isat,i)*fac(1)+gpsat(isat+1,i)*fac(2)
                end if
             end if
          end do
       else
          write(amessage,20)
20        format('Error in subroutine GET_CO2_SAT_LINE_PROPS_PRESSURE: ',  &
               'phase code must be "1" or "2".')
          ifail=1
          return
       end if
    end if

    return

  end subroutine get_co2_sat_line_props_pressure


  ! get_co2_sat_line_props_temperature
  subroutine get_co2_sat_line_props_temperat(ifail,iphase,ncode,icode,temperature,value)

    ! -- Subroutine GET_CO2_SAT_LINE_PROPS_TEMPERATURE gets saturation properties pertaining ',
    !    to a certain temperature. However derivatives are taken along the actual saturation line.

    implicit none

    integer, intent(out)      :: ifail
    integer, intent(in)       :: iphase
    integer, intent(in)       :: ncode
    integer, intent(in)       :: icode(ncode)
    real*8,    intent(in)       :: temperature
    real*8,    intent(out)      :: value(ncode)

    integer                   :: isat,i
    real*8                      :: tt
    real*8                      :: delta_p,delta_t
    real*8                      :: fac(2)

    ifail=0

    tt=(temperature-t_index_offset)*t_index_factor
    if(tt.gt.tsat_max)then
       ifail=-1
       return
    else if(tt.lt.tsat_min)then
       write(amessage,10)
10     format('Error in subroutine GET_CO2_SAT_LINE_PROPS_TEMPERATURE: ',  &
            'temperature out of table interpolation range.')
       ifail=1
       return
    else
       call which_index(tt,nsat,tsat,isat,sat_index_last_tt)
       call linear_interpolation_factors(tsat(isat),tsat(isat+1),tt,fac)
       delta_p=psat(isat+1)-psat(isat)
       delta_t=tsat(isat+1)-tsat(isat)
       if(iphase.eq.1)then
          do i=1,ncode
             if(icode(i).ne.0)then
                if((i.eq.3).or.(i.eq.6).or.(i.eq.9))then
                   value(i)=(lpsat(isat+1,i-2)-lpsat(isat,i-2))/delta_p
                else if((i.eq.2).or.(i.eq.5).or.(i.eq.8))then
                   value(i)=(lpsat(isat+1,i-1)-lpsat(isat,i-1))/delta_t
                else
                   value(i)=lpsat(isat,i)*fac(1)+lpsat(isat+1,i)*fac(2)
                end if
             end if
          end do
       else if(iphase.eq.2)then
          do i=1,ncode
             if(icode(i).ne.0)then
                if((i.eq.3).or.(i.eq.6).or.(i.eq.9))then
                   value(i)=(gpsat(isat+1,i-2)-gpsat(isat,i-2))/delta_p
                else if((i.eq.2).or.(i.eq.5).or.(i.eq.8))then
                   value(i)=(gpsat(isat+1,i-1)-gpsat(isat,i-1))/delta_t
                else
                   value(i)=gpsat(isat,i)*fac(1)+gpsat(isat+1,i)*fac(2)
                end if
             end if
          end do
       else
          write(amessage,20)
20        format('Error in subroutine GET_CO2_SAT_LINE PROPS_TEMPERATURE: ',  &
               'phase code must be "1" or "2".')
          ifail=1
          return
       end if
    end if

    return

  end subroutine get_co2_sat_line_props_temperat
  ! get_co2_sat_line_props_temperature



  subroutine interpolation_arrays_deallocate

    implicit none

    integer ierr

    if(allocated(t)) deallocate(t,stat=ierr)
    if(allocated(p)) deallocate(p,stat=ierr)
    if(allocated(property_type)) deallocate(property_type,stat=ierr)
    if(allocated(rarray)) deallocate(rarray,stat=ierr)
    if(allocated(satline)) deallocate(satline,stat=ierr)
    if(allocated(satclose)) deallocate(satclose,stat=ierr)
    if(allocated(tsat)) deallocate(tsat,stat=ierr)
    if(allocated(psat)) deallocate(psat,stat=ierr)
    if(allocated(csat)) deallocate(csat,stat=ierr)
    if(allocated(ssat)) deallocate(ssat,stat=ierr)
    if(allocated(msat)) deallocate(msat,stat=ierr)
    if(allocated(lpsat)) deallocate(lpsat,stat=ierr)
    if(allocated(gpsat)) deallocate(gpsat,stat=ierr)

    return

  end subroutine interpolation_arrays_deallocate




  subroutine addquote(afile,aqfile)

    ! -- Subroutine ADDQUOTE adds quotes to a filename if it has a space in it.

    character (len=*), intent(in)   :: afile
    character (len=*), intent(out)  :: aqfile
    integer nbb

    if(index(trim(afile),' ').eq.0)then
       aqfile=afile
    else
       aqfile(1:1)='"'
       aqfile(2:)=trim(afile)
       nbb=len_trim(aqfile)+1
       aqfile(nbb:nbb)='"'
    end if

    return
  end subroutine addquote


  integer function nextunit()

    ! -- Function nextunit determines the lowest unit number available for
    ! -- opening.

    logical::lopen

    do nextunit=10,300
       inquire(unit=nextunit,opened=lopen)
       if(.not.lopen) return
    end do

  end function nextunit


  subroutine tp2cell(it,ip,ic)

    ! -- Subroutine TP2CELL converts temperature and pressure to cell number.

    implicit none
    integer, intent(in)  :: it,ip
    integer, intent(out) :: ic

    ic=(ip-1)*nt+it

  end subroutine tp2cell




  SUBROUTINE LOWCASE(ASTRNG)

    ! -- Subroutine LOWCASE converts a string to lower case.

    INTEGER I,J
    CHARACTER*(*) ASTRNG

    DO 10 I=1,len_trim(ASTRNG)
       J=ICHAR(ASTRNG(I:I))
       IF((J.GE.65).AND.(J.LE.90)) ASTRNG(I:I)=CHAR(J+32)
10     CONTINUE
    RETURN
  END SUBROUTINE LOWCASE



end module property_interpolate



