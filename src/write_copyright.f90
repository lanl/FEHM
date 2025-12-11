subroutine write_copyright(termunit)
!***********************************************************************
! Copyright 2015. Los Alamos National Security, LLC.  This material was 
! produced under U.S. Government contract DE-AC52-06NA25396 for Los 
! Alamos National Laboratory (LANL), which is operated by Los Alamos 
! National Security, LLC for the U.S. Department of Energy. The U.S. 
! Government has rights to use, reproduce, and distribute this software.
! Neither the U.S. Government nor Los Alamos National Security, LLC or 
! persons acting on their behalf, make any warranty, express or implied,
! or assumes any liability for the accuracy, completeness, or usefulness
! of the software, any information pertaining to the software, or 
! represents that its use would not infringe privately owned rights.
!
! The software being licensed may be Export Controlled.   It may not be 
! distributed or used by individuals or entities prohibited from having 
! access to the software package, pursuant to United States export 
! control laws and regulations.  An export control review and 
! determination must be completed before LANS will provide access to the
! identified Software.
!***********************************************************************

  use comai, only : verno, jdate, jtime
  implicit none

  integer termunit

  write (termunit, 1) verno, jdate, jtime
  write (termunit, 2)
  write (termunit, 3)

1 format (a30, 3x, a11, 3x, a8, /)

2 format ('Copyright  2015.   Los Alamos National Security, LLC.  This material was', /, & 
       'produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos', /, &
       'National  Laboratory  (LANL),  which is operated by  Los Alamos National', /, & 
       'Security, LLC  for the U. S. Department of Energy.  The U. S. Government', /, & 
       'has rights to use, reproduce, and distribute this software.  Neither the', /, & 
       'U. S. Government nor Los Alamos National Security, LLC or persons acting', /, & 
       'on their behalf,  make any warranty,  express or implied, or assumes any', /, & 
       'liability for the accuracy, completeness, or usefulness of the software,', /, & 
       'any information pertaining to the software,  or  represents that its use', /, & 
       'would not infringe privately owned rights.', /)

3 format ('The  software  being licensed  may  be Export Controlled.  It may not be', /, &  
       'distributed  or  used by individuals  or entities prohibited from having', /, &  
       'access to the software package, pursuant to United States export control', /, &  
       'laws and regulations. An export control review and determination must be', /, &  
       'completed before LANS will provide access to the identified Software.', /)

end subroutine write_copyright
