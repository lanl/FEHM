      subroutine wtsi_column
!***********************************************************************
!  Copyright, 2005,  The  Regents  of the  University of California.
!  This program was prepared by the Regents of the University of 
!  California at Los Alamos National Laboratory (the University) under  
!  contract No. W-7405-ENG-36 with the U.S. Department of Energy (DOE). 
!  All rights in the program are reserved by the DOE and the University. 
!  Permission is granted to the public to copy and use this software 
!  without charge, provided that this Notice and any statement of 
!  authorship are reproduced on all copies. Neither the U.S. Government 
!  nor the University makes any warranty, express or implied, or 
!  assumes any liability or responsibility for the use of this software.
!***********************************************************************!D1
!D1  PURPOSE
!D1
!D1  This subroutine generates or reads node column data for 
!D1  free surface problems
!D1
!***********************************************************************
!D2
!D2 REVISION HISTORY
!D2
!D2 Initial implementation: 17-AUG-05, Programmer: Z. Dash
!D2 
!D2 $Log:   /pvcs.config/fehm90/src/wtsi_column.f_a  $
!D2 
!***********************************************************************
!D3
!D3  REQUIREMENTS TRACEABILITY
!D3
!D3 2.6 Provide Input/Output Data Files
!D3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
!D3
!***********************************************************************
!D4
!D4  SPECIAL COMMENTS AND REFERENCES
!D4
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4
!**********************************************************************
      
      use comai, only : icnl, ierr, igrav, iptty, m
      use combi, only : cord
      use comdi, only : izone_free_nodes, ps
      use comdti, only : n0
      use comwt
      use comxi

      implicit none
      logical null1, exists
      integer iroot, isuffix, nmax2, n_wt_col2
      integer i, ij, im, imm, ip, iswt, j, match, mm, nmax
      real*8, allocatable :: col_temp(:,:), n_col_temp(:)
      character*100 column_root

      allocate (wcol(n0))

! If there isn't a defined column data filename
! Generate column data filename
      if (null1(nmfil(27))) then
         if (nmfil(3) .ne. nmfily(3) .and. nmfil(3) .ne. nmfil(2)
     &        .and. nmfil(3) .ne. ' ') then
! Use root of grid file name
            call file_prefix(nmfil(3), iroot)
            if (iroot .gt. 94) iroot = 94
            column_root(1:iroot) = nmfil(3)(1:iroot)
         else if (.not. null1(root_name)) then
! Use root name if defined
            iroot = len_trim (root_name)
            if (iroot .gt. 94) iroot = 94
            column_root(1:iroot) = root_name(1:iroot)
         else if (nmfil(5) .ne. nmfily(3) .and. nmfil(5) .ne. ' ') then
! Use root of output file name
            call file_prefix(nmfil(5), iroot)
            if (iroot .gt. 94) iroot = 94
            column_root(1:iroot) = nmfil(5)(1:iroot)
         else if (nmfil(2) .ne. nmfily(3) .and. nmfil(2) .ne. ' ') then
! Use root of input file  name
            call file_prefix(nmfil(2), iroot)
            if (iroot .gt. 94) iroot = 94
            column_root(1:iroot) = nmfil(2)(1:iroot)
         else
            write (ierr, *) 'FILE ERROR: ',
     &           ' unable to determine column file prefix'
            stop
         end if
         nmfil(27)(1:iroot) = column_root(1:iroot)
         isuffix = len_trim(suffix(27))
         nmfil(27)(iroot+1:iroot+isuffix) = suffix(27)
      end if
! Check to see if the file exists
      inquire (file = nmfil(27), EXIST=exists)
      iswt = nufilb(27)
      open (unit = iswt, file = nmfil(27), form = cform(27),
     &     status = cstats(27))
      nmax = 0
      n_wt_cols = 0
      if (exists) then
! Make sure the data can be read
         read(iswt,'(2i10)', end = 305) n_wt_cols, nmax
         backspace iswt
 305     if (n_wt_cols .eq. 0) then
! Generate column data
            exists = .FALSE.
         end if
      end if
      if (.not. exists) then
! Generate the column data
c     figure this out and write to this file
!         allocate(col(n0,min(n0,1000)))
! Guess 1/5 of n0 is number of columns, 0.3 root of n0 is nodes in a column
         n_wt_col2 = n0 / 5
         nmax2 = n0**.3
         allocate (col(n_wt_col2, nmax2), n_col(n_wt_col2))
         n_wt_cols=0
         n_col=0
         col=0;column_read=1
         do i=1,n0
!            if(izone_free_nodes(i).ne.0..and.ps(i).gt.0.) then
               do ic=1,n_wt_cols
                  match=0
c     this code assumes 2-d or radial= igrav=2; 3-d= igrav=3
                  if(icnl.ne.0) then
                     if(cord(i,1).eq.cord(col(ic,1),1)) match=1
                  else
                     if(cord(i,1).eq.cord(col(ic,1),1).and.
     &                    cord(i,2).eq.cord(col(ic,1),2)) match=1
                  endif
                  if(match.eq.1) then
                     n_col(ic)=n_col(ic)+1
                     wcol(i)=ic
                     if (n_col(ic) .gt. nmax2) then
! Our col array isn't big enough, we need to resize it
                        allocate (col_temp(n_wt_col2, nmax2))
                        nmax = 1.1 * nmax2 + 1
                        col_temp = col
                        deallocate (col)
                        allocate (col(n_wt_col2, nmax))
                        col = 0
                        col(1:n_wt_col2, 1:nmax2) = col_temp
                        nmax2 = nmax
                        deallocate (col_temp)
                     end if
c     place in array according to igrav
                     do j=1,n_col(ic)-1
                        if(cord(col(ic,j),igrav) .lt. 
     &                       cord(i,igrav)) then
                           do mm=n_col(ic),j+1,-1
                              col(ic,mm)=col(ic,mm-1)
                           end do
                           col(ic,j)=i
                           goto 44
                        endif
                     end do
c     if this is the lowest z node
                     col(ic,n_col(ic))=i
                     goto 44
                  endif
               end do
               n_wt_cols=n_wt_cols+1
               if (n_wt_cols .gt. n_wt_col2) then
! Our arrays aren't big enough, we need to resize them
                  allocate (col_temp(n_wt_col2, nmax2), 
     &                 n_col_temp(n_wt_col2))
                  col_temp = col
                  n_col_temp = n_col
                  deallocate (col, n_col)
                  nmax = 1.1 * n_wt_col2 + 1
                  allocate (col(nmax, nmax2), 
     &                 n_col(nmax))
                  col = 0
                  col(1:n_wt_col2, 1:nmax2) = col_temp
                  n_col(1:n_wt_col2) = n_col_temp
                  n_wt_col2 = nmax
                  deallocate (col_temp, n_col_temp)
               end if
               n_col(n_wt_cols)=1
               col(n_wt_cols,1)=i
               wcol(i)=n_wt_cols
!            endif
 44      end do
c     determine max # of nodes in a column
         nmax=0
         do ij=1,n_wt_cols
            nmax=max(nmax, n_col(ij))
         end do
c     write results to file
         write(iswt,'(2i10)') n_wt_cols,nmax
         do ij=1,n_wt_cols
            do ip=1,n_col(ij)
               im=col(ij,ip)
               write(iswt,304) ij,ip,im,cord(im,1),
     &              cord(im,2),cord(im,3)
            end do
         end do
 304     format(3i15,3f15.3,i10)
         rewind(iswt)
         deallocate (col, n_col)
      end if
! if file exists, read from the file, otherwise 
! we are re-reading so we minimize space allocation
      read(iswt,'(2i10)') n_wt_cols, nmax
      allocate (col(n_wt_cols, nmax), n_col(n_wt_cols))
      do i=1,1000000
         read(iswt,*,end=566) ij,im,imm
         n_col(ij)=im
         col(ij,im)=imm
c     izone_free_nodes(imm)=ij
         wcol(imm)=ij
      end do
 566  n_wt_cols=ij
      
      close(iswt)
      allocate (col_out(m))

      end
