      subroutine wtsictr(iflg)
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
CD1  This subroutine manages initization and boundary conditions 
CD1  that involve free surfaces
CD1
C***********************************************************************
CD2
CD2 REVISION HISTORY
CD2
CD2 Initial implementation: 19-APR-02, Programmer: George Zyvoloski
CD2 
CD2 $Log:   /pvcs.config/fehm90/src/wtsictr.f_a  $
!D2 
!D2    Rev 2.5   06 Jan 2004 10:44:32   pvcs
!D2 FEHM Version 2.21, STN 10086-2.21-00, Qualified October 2003
!D2 
!D2    Rev 2.4   29 Jan 2003 09:25:22   pvcs
!D2 FEHM Version 2.20, STN 10086-2.20-00
CD2
C***********************************************************************
CD3
CD3  REQUIREMENTS TRACEABILITY
CD3
CD3 2.6 Provide Input/Output Data Files
CD3 3.0 INPUT AND OUTPUT REQUIREMENTS                  
CD3
C***********************************************************************
CD4
CD4  SPECIAL COMMENTS AND REFERENCES
CD4
CD4 Requirements from SDN: 10086-RD-2.20-00
CD4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
CD4   FEHM Application Version 2.20
CD4
C**********************************************************************
c iflag
c     -1 sort nodes by column, adjust constant head nodes, source/sinks
c      0 read wtsi macro
c      1 initialize head12 arry, izone_ifree_nodes
c      2 adjust rlxyf, etc.
c      3
c      4
c      5
c      6 initialization:  sort wtsi nodes by column
c      7 check to see if any sources need to be moved
c      9 inititalize dzrg array (cell thickness)
c     11 output water table elevation
c     12 flow_wt macro or a water table output desired
      use comdti
      use comai
      use combi
      use comci
      use comdi
      use comii
      use comfi
      use comwt
      use comxi
      implicit none
c     
      integer iflg,icode, izone, inode,imove,match,im,imm,in,iu,ip
      integer mi,neqp1,i,i1,i2,ii,j,jj,kb,itop,mm,ij,iij, ihead_ck
      integer ishfile, ishfile2, open_file, inodeb
      integer nmove, mmm, iwm, k, itp
      real*8  hnode, hmax, hmin, x1,x2,min_sat,cut,hmid
      real*8  htol, dentol, cord_z_max, dis_ck, cap_val
      real*8 head_ck1, head_ck2, strd_wtsi, hgrid
      real*8 hfac_h,hfac_l,pt,p_dum,wt_elev2
      character*80 form_string
      character*3 dum_wtsi
      logical:: debug_flag = .FALSE.
      logical:: zonewt_flag = .FALSE.
      parameter (htol=1.d-15,dentol=1.d+1)
      parameter (min_sat=1.d-5, strd_wtsi= 1.d00)
      parameter (hfac_h=1.00d0, hfac_l=1.00d0)

! Defined as parameter in comwt
!      sattol = 1.d-1
!      phi_inc = head0*crl(1,1)*(-grav)

c     flag=7 
      cut=0.0
c================================================================
      if(ifree.eq.0) return
c================================================================
      if(iflg.eq.0) then
c     
c     read input       
c     
c     contr option ifree = -1 (heads,sat) , +1 (phi,sat)
         ifree = -1
         wtsi_isot = 0
         isw_term = 0
         zfac_ani = 1.d0
         sattol = 0.1
c
c iad_up_wtsi fixed at 1000 (no explicit component)
c
         iad_up_wtsi = 1000
         read(inpt,'(a80)') wdd1
         do i = 1,80
          if(wdd1(i:i).eq.'a'.or.wdd1(i:i).eq.'A') then
           zonewt_flag = .true.
          endif
         enddo
         if(zonewt_flag) then
          nfree = 1
         else
          backspace inpt
          read(inpt,*) nfree
         endif
         if(.not.allocated(izone_free)) then
            allocate(izone_free(max(1,nfree)))
            allocate(ifreef(max(1,nfree)))
            allocate(ifree_im_ex(1))
            allocate(dry_zone(n0))
            allocate(t91(nn))
            allocate(dfidf(nn))
            allocate(dfid1f(nn))
            allocate(rlzf(n0))
            allocate(rlxyf(n0))
            allocate(drlxyf(n0))
            allocate(drlzf(n0))
c            allocate(dpcef(n0))
            allocate(izone_free_nodes(n0))  
            dry_zone = 0 
         end if   
         allocate(head12(n0,2))
         allocate(dzrg(n0))
         backspace inpt
          if(zonewt_flag) then
           izone_free(1) = -1
           read(inpt, *, end=555, err=555) dum_wtsi,
     &        head_tol,rlptol,zfac_ani,head_id,sattol
          else
           read(inpt, *, end=555, err=555) nfree, 
     &        (izone_free(i),i=1,nfree),head_tol,rlptol,zfac_ani,
     &        head_id,sattol          
          endif

         

         if(abs(zfac_ani).lt.1.0d0) then
	    zfac_ani = 1.0d0
	    wtsi_isot = 0
         else if(zfac_ani.lt.0.d0) then
            wtsi_isot = 2
            zfac_ani = abs(zfac_ani)
         else
	    wtsi_isot = 1
         endif

         go to 556
 555     continue
         head_tol = 1.d-2
         isw_term = 0
         wtsi_isot = 0
 556     continue	     
c         iad_up= min(iad_up,mm)
c         iad_up= max(iad_up,-1)

c     Loop over each zone for determining izone_free array

         izone_free_nodes=0
c     set number of partialy filled cells = 0
         ifree1 = 0
         if(zonewt_flag) then
           do inode = 1, n0          
              izone_free_nodes(inode) = 1              
            end do
         else
          do izone = 1, nfree
            do inode = 1, n0
               if(izonef(inode).eq.izone_free(izone)) then
                  izone_free_nodes(inode) = izone_free(izone)
               end if
            end do
          end do
         endif
    
c     

      else if(iflg.eq.12) then
c     
c     if there's a flow_wt macro or a water table output desired
c     sort wtsi nodes into columns ordered by igrav coordinate
c     
c         if(ishiswt.ne.0.or.move_wt.ne.0.) then
            call wtsi_column
c         endif
c     adjust constant head nodes to ensure that h(node) >= 
c     h(bottom of cell)
c     If flow_wt node, place specified head
c     at appropriate node (where head=z(node)).  If not flow_wt
c     node, release specified head if head < z(node)
c     
c     modified 7-26-05 gaz so seepage face (ka=-3 still works)
         if(l.eq.0) then
            do inode=1,n0
               if(izone_free_nodes(inode).ne.0.and. ka(inode).lt.0) then
                  call headctr(4,inode,head12(inode,1),hmin)
                  call headctr(4,inode,head12(inode,2),hmax)
                  call headctr(4,inode,pflow(inode),p_dum)
                  hnode = cord(inode,abs(igrav))
                  iu=0
c     check to see what kind of node this is
                  if(move_wt.eq.1)then
                     if(move_type(inode).eq.1) iu=1
                  endif
                  if(iu.eq.0) then
c     if this is a normal flow macro, check to see if specified head 
c     is reasonable, if not, release it
                     if(hmin.gt.p_dum) then
                        wellim(inode)=1.d-15
                        s(inode) = rlptol
                        ka(inode) = 0
                        if (debug_flag) then
                           if (iptty .ne. 0) then
                              write(iptty, 101) inode
                              write(iptty, 102) hnode, p_dum
                              write(iptty, 103) cord(inode,1),
     &                             cord(inode,2)
                           end if
                           if (iout .ne. 0) then
                              write(iout, 101) inode
                              write(iout, 102) hnode, p_dum
                              write(iout, 103) cord(inode,1), 
     &                             cord(inode,2)
                           end if
                        end if
                     endif
                  else
c     if this is a flow_wt macro, check to see if specified head is 
c     reasonable, if not, move to correct node
                     if(p_dum.lt.hmin.or.p_dum.gt.hmax) then
                        ic=wcol(inode)
                        do i=1,n_col(ic)
                           in=col(ic,i)
                           call headctr(4,in,head12(in,1),hmin)
                           call headctr(4,in,head12(in,2),hmax)	
                           if(p_dum.le.hmax.and.p_dum.ge.hmin) then
c     
c     found the correct node;set head to desired value
c     
                              wellim(in)=wellim(inode)
                              ka(in)=-1
                              pflow(in)= pflow(inode)
                              esk(in)=esk(inode)
                              move_type(in)=1
                              izoneflxz(in)=izoneflxz(inode)
                              izonef(in)=izonef(inode)
c     
                              if (iptty .ne. 0) write(iptty, 104) 
     &                             inode, in, p_dum
                              if (iout .ne. 0) write(iout, 104) 
     &                             inode, in, p_dum
c     
c     also set all heads below that node - 
c     
                              do ij=i+1,n_col(ic)
                                 in=col(ic,i)
c     call headctr(5,in,pt,pflow(inode))
                                 pho(in)=pflow(inode) + 
     &                                head0*crl(1,1)*(-grav)
                                 phini(in) = pho(in)
                                 phi(in)=pho(in)
                              end do	
c     needs to be called in startup where initial values are set
c     
                              wellim(inode)=1.d-15
                              s(inode) = rlptol
                              ka(inode) = 0	
                              pflow(inode)=0.0
                              goto 803	
                           endif
                        enddo	
                        if (iptty .ne. 0) then
                           write(iptty, 105) inode
                           write(iptty, 106) p_dum
                        end if
                        if (iout .ne. 0) then
                           write(iout, 105) inode
                           write(iout, 106) p_dum
                        end if

c     release head from original node
                        wellim(inode)=1.d-15
                        s(inode) = rlptol
                        ka(inode) = 0
                        pflow(inode)=0.0	
                        goto 803	 
                     endif
                  endif
               endif
 803        enddo
         end if
         if (debug_flag) then
            ishfile = open_file('new-specified-heads','unknown')
            ishfile2 = open_file('new-specified-heads.csv','unknown')
            write(ishfile,'(a4)') 'flow'
            do ii=1,n0
               call headctr(4,ii,pflow(ii),hmin)
               if(ka(ii).lt.0.) write(ishfile,203) ii,ii,1,hmin,1.,
     &              wellim(ii)/1.0e+06
               write(ishfile2,206) (cord(ii,mmm),mmm=1,3),hmin
            end do
            write(ishfile,'(a4)') ' '
            close(ishfile)
            close(ishfile2)
         end if

 101     format (1x,'releasing head in flow macro at node ', i8)
 102     format (1x,'z(node) = ', g16.9, ' applied head was ', g16.9)
 103     format (1x,'x,y = ', 2(g16.9, 1x))
 104     format (1x,'spec head ', i8,' moved to ',i8,' head=',g16.9)
 105     format (1x,'there are no wtsi nodes for spec head ', i8)
 106     format (1x,'head ', g16.9, 'is released')
 203     format(3i10,3f15.3)
 206     format(4(f15.3,','))

      else if(iflg.eq.1) then
c     
c     calculate max and min pressures associated with cell 
c     
c     calculate head values for all nodes (based on pho)
c     
c     first remove specified head nodes from wtsi list
c     
c save this commented out code for later 
c         do inode=1,n0
c            if(ka(inode).eq.-1) izone_free_nodes(inode)=0
c         enddo
c  fill the array head
c
c         rho1grav = crl(1,1)*(9.81d-6)
         call headctr(2,0,0.,0.)
	 head12=0.
         do inode=1,n0

            if(izone_free_nodes(inode).ne.0) then
               i1=nelm(inode)+1
               i2=nelm(inode+1)
               hmax=cord(inode,igrav)
               hmin=cord(inode,igrav)
               hgrid = 0.5d0*dzrg(inode)
               hmax=cord(inode,igrav)+hgrid
               hmin=cord(inode,igrav)-hgrid
               if (rlp_flag .eq. 1) then
                  if(icap(inode).ne.0) then
                     itp= icapt(icap(inode))
                     if(itp.eq.1) then
                        cap_val = cp1f(icap(inode))/rho1grav
                     else
                        cap_val = 0.0d0
                     endif
                     hmin = min(hmin,cord(inode,igrav)-cap_val)
                  endif
               end if
               if(head(inode).lt.hmin) then
c     the cell is dry
c     count as partially filled cell
                  ifree1 = ifree1+1 
                  izone_free_nodes(inode) = 3
                  call headctr(5,inode,pho(inode),hmin)
c                  pho(inode) = crl(4,1) + phi_inc
               else if(head(inode).lt.hmax) then   
c cell is partially full
                  ifree1 = ifree1+1 
                  izone_free_nodes(inode) = 2 
               else
c cell is full
                  izone_free_nodes(inode) = 1                     
               endif
               if(abs(hmax-hmin).le.htol) then
                  izone_free_nodes(inode)=0
               else
c     this call changes head12 from meters to pressure
                  call headctr(5,inode,head12(inode,1),hmin)
                  call headctr(5,inode,head12(inode,2),hmax)
               endif
            endif
         enddo
         phi = pho

c     
      else if(iflg.eq.2) then
c
c call relative perm module for wtsi 
c
         call rlperm_wtsi(1)   
c         
      else if (iflg.eq.7) then
         nmove = 0
         do i=1,n0
            imove=0
            if(ka(i).eq.1.and.move_type(i).gt.0.) then
c     this is a moving water table source
               ij=wcol(i)
               do im=1,n_col(ij)
c we're searching from the top
                  inode=col(ij,im)
                  if(s(inode).gt.rlptol) then
                     if(inode.ne.i) then
                        imove=inode
                        nmove=nmove+1
                        goto 88
                     endif
                     goto 88
		  endif
               end do
!              write(33,*) day,' all wtsi cells in column ',i,' are dry'
			    imove=-1
 88            if(imove.gt.0) then
c     move the source/sink to node imove only
c     if imove is not already a sink/source move source to target
c     (note: need to figure out how to get the source back
c     if water table rises later)
                  if(ka(imove).eq.0) then
                     sk(imove)=sk(i)
                     qc(imove)=qc(i)
                     wellim(imove)=wellim(i)
                     esk(imove)=esk(i)
                     pflow(imove)=pflow(i)
                     ka(imove)=ka(i)
                     move_type(imove)=move_type(i)
c     transfer node in flxz array 
                     izoneflxz(imove)=izoneflxz(i)
                     izonef(imove)=izonef(i)
                  endif
c     remove source from the node 
!                  tcc=tcc+sk(imove)
                  sk(i)=0.
                  wellim(i)=0.
                  esk(i)=0.
                  pflow(i)=0.
                  ka(i)=0.
                  qc(i)=0.
                  move_type(i)=0.
                  izoneflxz(i)=0
                  izonef(i)=0
               endif
            endif
         end do
         continue		

         if (iout .ne. 0) write(iout,*) 'moving sources: ',nmove
         if (iptty .ne. 0) write(iptty,*) 'moving sources: ',nmove
c     
      else if(iflg.eq.8)	then
c     
c     set last time step saturation for phase change information
c     
         so = s

      else if(iflg.eq.9)	then
c     
c     calculate the gridblock length in the gravity direction
c     

c         rho1grav = crl(1,1)*(9.81d-6)
c         rho1grav = 1.
	 do i = 1,neq
            i1=nelm(i)+1
            i2=nelm(i+1)
            hmid=cord(i,igrav)
            hmin=0.
	    hmax=0.
            do ii =i1,i2
               kb=nelm(ii)
               hmax=max(cord(kb,igrav)-hmid,hmax)
               hmin=min(cord(kb,igrav)-hmid,hmin)
            enddo
            if(ivf.eq.-1) then
               dzrg(i) = max(hmax,abs(hmin))
            else
               dzrg(i) = abs(hmax-hmin)/2.	            
            endif      
         enddo
      else if (iflg.eq.11) then
c     write wt output for time history
c         do i=1,n_wt_cols
         iwm=0
         do ij=1,m
            i=wcol(nskw(ij)) 
            do k=1,iwm
               if(i.eq.col_out(k)) goto 566
            end do
            iwm=iwm+1
            col_out(iwm)=i
            do im=n_col(i),1,-1
               inode=col(i,im)
               if(s(inode) .lt. 1. .or. im .eq. 1) then
                  if(im .ne. 1 .or. s(inode) .lt. 1.) then
                     call headctr(4,inode,phi(inode),wt_elev)
                  else
                     wt_elev=cord(inode,3)+dzrg(i)/rho1grav
                  end if
                  if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                     inode=col(i,im+1)
                     call headctr(4,inode,phi(inode),wt_elev)
                  endif
                  if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                     if (iptty .ne. 0) write(iptty, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                     if (iout .ne. 0) write(iout, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                  endif
                  goto 4009
               endif
            end do
 4009       continue
c 	search from top
            do im=1,n_col(i)
               inode=col(i,im)
               if(s(inode).gt.rlptol) then
                  call headctr(4,inode,phi(inode),wt_elev2)
                  if(wt_elev2.eq.0..and.im.ne.n_col(i)) then
                     inode=col(i,im+1)
                     call headctr(4,inode,phi(inode),wt_elev2)
                  endif
                  if(wt_elev2.eq.0..and.im.eq.1) then
                     if (iptty .ne. 0) write(iptty, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                     if (iout .ne. 0) write(iout, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                  endif
                  goto 4007
               endif
            end do
 4007       continue
            if (form_flag .eq. 2) then
               write(ishiswt,4004) days, cord(inode,1),
     &              cord(inode,2), cord(inode,3), wt_elev,
     &              ps(inode), inode, nskw(ij)
c     &              wt_elev2, ps(inode), inode, nskw(ij)
            else
               write(ishiswt,4005) days, cord(inode,1),
     &              cord(inode,2), cord(inode,3), wt_elev,
     &              wt_elev2, ps(inode), inode, nskw(ij)
            end if
 566     end do
 4004    format(1x,7(g16.9,', '),i8,', ',i8)
 4005    format(1x,7(g16.9,1x),2(i8,1x))
 4006    format(1x,'all wtsi nodes are dry ', i8, 1x, 2(g16.9, 1x))

      else if (iflg.eq.14) then
c     write wt output for contours
         if (altc(1:4) .eq. 'avsx') then
            write (form_string, 4015) ' : ', ' : '
            write (isconwt, 4019)
         else if (altc(1:3) .eq. 'sur') then
            write (form_string, 4015) ', ', ', '
            write (isconwt, 4020)
         else if (altc(1:3) .eq. 'avs' .or. altc(1:3) .eq. 'tec') then
            write (form_string, 4015) ' ', ' '
            if (altc(1:3) .eq. 'tec') then
               write (isconwt, 4018) days
            else
               write (isconwt, '("06 1 1 1 1 1 1")')
               write (isconwt, '(a)') 'X coordinate (m), (m)'
               write (isconwt, '(a)') 'Y coordinate (m), (m)'
               write (isconwt, '(a)') 'Z coordinate (m), (m)'
               write (isconwt, '(a)') 'Zone (no dim), (no dim)'
               write (isconwt, '(a)') 'Water table elevation (m), (m)'
               write (isconwt, '(a)') 'Water table elevation2 (m), (m)'
            end if            
         end if
         do i=1,n_wt_cols
            do im=n_col(i),1,-1
               inode=col(i,im)
               if(s(inode).lt.1.or.im.eq.1) then
                  call headctr(4,inode,phi(inode),wt_elev)
                  if(wt_elev.eq.0..and.im.ne.n_col(i)) then
                     inode=col(i,im+1)
                     call headctr(4,inode,phi(inode),wt_elev)
                  endif
                  if(wt_elev.eq.0..and.im.eq.n_col(i)) then
                     if (iptty .ne. 0) write(iptty, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                     if (iout .ne. 0) write(iout, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                  endif
                  goto 4010
               endif
            end do
 4010       continue
            inodeb = inode
c 	search from top
            do im=1,n_col(i)
               inode=col(i,im)
               if(s(inode).gt.rlptol) then
                  call headctr(4,inode,phi(inode),wt_elev2)
                  if(wt_elev2.eq.0..and.im.ne.n_col(i)) then
                     inode=col(i,im+1)
                     call headctr(4,inode,phi(inode),wt_elev2)
                  endif
                  if(wt_elev2.eq.0..and.im.eq.1) then
                     if (iptty .ne. 0) write(iptty, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                     if (iout .ne. 0) write(iout, 4006) inode,
     &                    cord(inode,1), cord(inode,2)
                  endif
                  goto 4011
               endif
            end do
 4011       continue
            write(isconwt, form_string) cord(inode,1), cord(inode,2),
     &           cord(inodeb,3), izonef(inodeb), wt_elev, wt_elev2
         end do

 4015    format("(1x, 3(g16.9, '", a, "'), i4, 2('", a, "', g16.9))")
 4020    format(1x, 'X (m), Y (m), Z (m), Zone, WT elev (m), ', 
     &        'WT elev2 (m)')
 4019    format(1x, 'X (m) : Y (m) : Z (m) : Zone : WT elev (m) : ',
     &        'WT elev2 (m)')           
 4018    format('variables = "X (m)" "Y (m)" "Z (m)" " Zone" ', 
     &        '"WT elev (m) "', '"WT elev2 (m)"', / 
     &        'zone t = "Simulation time ', g16.9, ' days"') 

      else if (iflg.lt.0) then
c     find wt elevation in that column
         in=-1.*iflg
         wt_elev=0.
         i=wcol(in)
         if(i.eq.0) return
         do im=n_col(i),1,-1
            inode=col(i,im)
            if(s(inode).lt.1.) then
               call headctr(4,inode,phi(inode),wt_elev)
               return
            endif
         end do
         wt_elev=cord(inode,3)-dzrg(inode)
      endif

c     
      return
      end                
