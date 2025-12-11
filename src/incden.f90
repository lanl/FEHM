subroutine incden(inpt)
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
  ! Read concentration-dependent density input for all liquid species
  
  use comai, only : ierr, iout, iptty
  use comci, only : cden
  use comchem, only : ncpnt, cpntnam
  use comrxni, only : mw_speci
  implicit none

  integer :: i, inpt, unit, open_file
  integer :: nwds, imsg(1), msg(1)
  real*8 :: xmsg(1), mw
  character(4) :: char4
  character(20) :: name
  character(32) :: cmsg(1)
  character(100) :: table_file, dumstring
  logical :: null1

  ! Read table of molecular weights, only assigning values for aqueous species
  allocate (mw_speci(ncpnt))
  read (inpt, *) char4
  if (char4 .eq. 'file') then
     ! We will read the table from an external file
     read (inpt, '(a100)') table_file
     unit = open_file(table_file, 'old')
     ! Read past header lines (assumes they start with a character)
     do 
        read (unit, '(a4)') char4
        call parse_string(char4,imsg,msg,xmsg,cmsg,nwds)
        if (msg(1) .ne. 3) then
           backspace (unit)
           exit
        end if
     end do
  else
     backspace (inpt)
     unit = inpt
  end if

  mw_speci = 0.
  do
     read (unit, '(a100)', end=1) dumstring
     if (null1(dumstring)) exit
     read(dumstring, *) mw, name
     do i = 1, ncpnt
        if (name .eq.  cpntnam(i)) then
           mw_speci(i) = mw           
           exit
        end if
     end do
  end do

1 cden = .true.
  if (unit .ne. inpt) close (unit)

end subroutine incden
  
