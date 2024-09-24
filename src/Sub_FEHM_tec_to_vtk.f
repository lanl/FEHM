       subroutine FEHM_tec_to_vtk(iflg)
c      program gaz
c gaz 013022 changed to subroutine
c program to combine FEHM tecplot files to tecplot 7.x files
c george zyvoloski initial implementation 041514
c modified to add multi-zone files and separate (or none) geometry files    
c 081319 gaz added  velocity vector 
c 091622 gaz added gdkm output      
c
      use avsio, only : avs_root, icall_max, iovelocity,
     &  ioconcentration, ioscalar, iomaterial, iogdkm
      use comai, only : verno, jdate, jtime, iccen
      
      implicit none
      real*8 x0,y0,z0,x1,x2,y1,y2,z1,z2,timez
      real*8 coor_tol
      real*8 aminvx,amaxvx,aminvy,amaxvy,aminvz,amaxvz,vscale
      integer neq, i, j, ii, jj, ij, i1, i2, i_wt,kk, nmax, num_files
      integer ichar1, ichar2, ichar3, idem, ic, num_var, num_tot_i, ns
      integer ichar4, izone, max_var, ifile, inode, iword_len, icnl
      integer maxlines, lines, node_elem, num_elem, total_char
      integer ivtk_elem, izid_position
      integer ivec
c gaz 071522      
      integer icall_max_use
c gaz 020922      
      integer j_vtk, jj_vtk, i_vtk, iq
      character*200 input_data_file
      character*200 temp_file
      character*200 contour_file, appended_file
      character*132 time_text,ele_char
      character*400 dum_zone1, dum_zone2, dum_zone3, dum_var
c gaz 020522 added log file definition      
      character*400 cont_log_file
      character*400 vtk_log_file
      character*145 dum_zone4
      logical file_status, exists

      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8, allocatable :: z(:)
      real*8, allocatable :: temp(:,:)
      integer, allocatable :: node_i(:)
      integer, allocatable :: idum(:,:)
      integer, allocatable :: zid_save(:)
      character*50, allocatable :: title_var(:)
      character*50, allocatable :: title_var0(:)
      character*30 tail_tec1, tail_tec2
      character*200, allocatable :: tec_file_name(:)
      
c gaz 020122 added flag structure      
      integer iflg, j_cont_type , n_cont_types
      parameter (nmax = 20, maxlines = 100000, n_cont_types = 5)
c gaz 070922 changed  iword_len = 15 to iword_len = 50     
      parameter (max_var = 20, total_char = 400, iword_len = 50)
      parameter (coor_tol = 1.d-12)
      parameter (vscale = 1.d00)
   
      ivec = 0
      if(iflg.eq.0) then
c initialize some vtk variables          
      else 
c read existing tec files
c gaz 020522 changing control file to avs log file that includes tecplot info
c initialize velocity id          
      ivec = 0
      cont_log_file = trim(avs_root) 
      iq = scan(cont_log_file,'.',back=.true.)
      write(cont_log_file(iq+1:iq+6),10) '00001'
10    format(a5)        
c remove tec log file if exists
      temp_file(1:iq) = cont_log_file(1:iq)
      write(temp_file(iq+1:iq+7),'(a7)') 'avs_log'
      inquire(file = temp_file, exist = exists)
c 
c create vtk         
      
      vtk_log_file(1:iq) = cont_log_file(1:iq)
      vtk_log_file(iq:iq+13) = '_vtk_log.file'
      open(unit=8,file=vtk_log_file,status='unknown')
      write(8,1) verno, jdate, jtime
1     format(a30,1x,a11,1x,a8,
     &       /,'This program combines FEHM tecplot files',
     &       /,'creates multiple vtk files with geometry attached ',
     &       /,'NOTE: we assume that the coordinate and element ',
     &       /,'data is appended to the first tecplot file',
     &       /,'vtk files are readable by PARAVIEW or SVS')
c     allocate array space
      
      do j_cont_type = 1, n_cont_types
       if(.not.allocated(tec_file_name))
     &   allocate(tec_file_name(icall_max))
       if(j_cont_type.eq.1) then
           if(ioscalar.eq.0) go to 100
           write(8,*)  '>>> Checking for scalar variables <<<'
           icall_max_use = icall_max
       else if(j_cont_type.eq.2) then
         if(iovelocity.eq.0) go to 100
         write(8,*)
         write(8,*)  '>>> Checking for velocity variables <<<'
         icall_max_use = icall_max
       else if(j_cont_type.eq.3) then
           write(8,*)
           if(ioconcentration.eq.0) go to 100
           write(8,*)  '>>> Checking for concentration variables <<<'
           icall_max_use = icall_max
       else if(j_cont_type.eq.4) then
           write(8,*)
           if(iogdkm.eq.0) go to 100
           write(8,*)  '>>> Checking for gdkm variables <<<'
           icall_max_use = icall_max           
       else
           write(8,*)
           if(iomaterial.eq.0) go to 100
           write(8,*)  '>>> Checking for a material file <<<' 
c gaz 071522  only one material file time = 0.0
           timez = 0.0d0
           icall_max_use = 1
       endif
      do j = 1, icall_max_use
       ii = 100000+j
       if(j_cont_type.eq.1) then
         write(tail_tec1,11) ii,'_sca_node.dat'
         write(cont_log_file(iq+1:iq+18),12) tail_tec1(2:19)
       else if(j_cont_type.eq.2) then 
         write(tail_tec1,11) ii,'_vec_node.dat'
         write(cont_log_file(iq+1:iq+18),12) tail_tec1(2:19)
       else if(j_cont_type.eq.3) then 
         write(tail_tec1,11) ii,'_con_node.dat'   
         write(cont_log_file(iq+1:iq+18),12) tail_tec1(2:19)
       else if(j_cont_type.eq.4) then 
         write(tail_tec1,111) ii,'_sca_gdkm_node.dat'   
         write(cont_log_file(iq+1:iq+23),112) tail_tec1(2:24)         
       else if(j_cont_type.eq.5) then 
         write(tail_tec1,'(a)') '_mat_node.dat'  
         write(cont_log_file(iq+1:iq+18),'(a)') tail_tec1(2:19)
       endif
11     format(i6,a13) 
111    format(i6,a18) 
12     format(a18) 
112     format(a23)        
       inquire(file = cont_log_file, exist = exists)
        if(exists.eqv..true.) then
          jj = len_trim(cont_log_file)
          tec_file_name(j)(1:jj) = cont_log_file(1:jj)
        else
          write(8,*) 'file ', cont_log_file, 'does not exist, stopping'
          stop
        endif
      enddo
c process tec files and create VTK files   

      do j_vtk = 1, icall_max_use
       jj_vtk = 0   
       do i2 = 1,196
         if(tec_file_name(j_vtk)(i2:i2+3).eq.'.dat') then
           jj_vtk = i2+3   
           go to 200
         endif 
       enddo    
200    do i2 = 1,200
          write(contour_file(i2:i2),'(a)') ' '
          write(input_data_file(i2:i2),'(a)') ' '
       enddo
c       jj_vtk = len_trim(tec_file_name(j_vtk))
       input_data_file(1:jj_vtk) = tec_file_name(j_vtk)(1:jj_vtk)
       contour_file(1:jj_vtk) = input_data_file(1:jj_vtk)
       i_vtk = scan(contour_file,'.', back = .true.)
       contour_file(i_vtk+1:i_vtk+3) = 'vtk'

      write(8,3) contour_file(1:i_vtk+3)
      
2     format('tec contour file to be converted = ', a80)
3     format(/'output vtk contour file = ', a)
4     format('tec contour files to be processed = ', i6,/)
      if(j_vtk.eq.1) then      
c
c  open initial files  
c
      open(unit=9,file=input_data_file,status='old')
      open(unit=10,file=contour_file,status='unknown')

c
c  read headers from first tec file
c
      read(9,'(a400)') dum_zone1

      read(9,'(a400)') dum_zone1
  
c count variables for vtk output      
      ic = 0
      num_var = 0
      do i = 1,total_char
       if(dum_zone1(i:i).eq.'"') then
        ic = ic + 1
        if(mod(ic,2).ne.0) then
         num_var = num_var+1
        endif
       endif
      enddo 
c
      dum_var(1:total_char) = dum_zone1(1:total_char)
      read(9,'(a400)') dum_zone1
c
       do i = 1, total_char-20
       if(dum_zone1(i:i+3).eq.'time') then
        ichar1 = i +5
        read(dum_zone1(ichar1:ichar1+13),'(1g13.04)') timez
       endif
       if(dum_zone1(i:i+3).eq.'N =') then
        read(dum_zone1(i+4:i+13),*) neq
       endif
       if(dum_zone1(i:i+5).eq.', E = ') then
        read(dum_zone1(i+6:i+13),*) num_elem
       endif
       if(dum_zone1(i:i+3).eq.'DATA') then
        ichar2 = i
       endif
       if(dum_zone1(i:i+1).eq.'FE') then
        ichar3 = i+2
c vtk triangle = 5        
        if(dum_zone1(ichar3:ichar3+2).eq.'TRI') then
         ivtk_elem = 5
         idem = 2
         ns = 3
c vtk quad = 9      
        else if(dum_zone1(ichar3:ichar3+2).eq.'QUA') then
         ivtk_elem = 9
         idem = 2
         ns = 4
c vtk tet = 10      
        else if(dum_zone1(ichar3:ichar3+2).eq.'TET') then
         ivtk_elem = 10
         idem = 3
         ns = 4
c vtk hex (brick) = 12      
        else if(dum_zone1(ichar3:ichar3+2).eq.'BRI') then
         ivtk_elem = 12
         idem = 3
         ns = 8
        endif
       endif
       enddo
c gaz 071522 use timez = 0.0 for material file
c gaz set timez = 0.0 beginning of
      write(8,*)'simulation time for this file = ', timez
      write(8,*)
      write(8,*) 
     & '>>>>>> ', num_var-idem-1,' contour variables found',' <<<<<<<'
c      write(8,*)
c gaz 021022 only search variables on first file      

c
c sort out variables for vtk output  
c
      if(.not.allocated(title_var)) then
        allocate (title_var(max_var))
        allocate (title_var0(max_var))
      endif
      do i = 1, max_var
        title_var(i)= '            '
        title_var0(i)= '            '
      enddo
      ic = 0
      num_var = 0
      izid_position = 0
      do i = 1,total_char-50
       if(dum_var(i:i).eq.'"') then
        ic = ic + 1
        if(mod(ic,2).ne.0) then
         num_var = num_var+1
c gaz 070922 discover variable name and fill in spaces with '_'
         do ii = i+1, i + iword_len
          if(dum_var(ii:ii).eq.' ') dum_var(ii:ii) = '_'
          if(dum_var(ii:ii).eq.'"') then
            i2 = ii-i-1
            title_var0(num_var)(1:i2) = dum_var(i+1:ii-1)  
            go to 101
          endif
         enddo
101       if(num_var.gt.idem+1) then
          jj = num_var-(idem+1)
c gaz 042422 
          title_var(jj)(1:i2) = title_var0(num_var)(1:i2)
c gaz 071522 removed various titles
300       continue
         endif
        endif
       endif
       enddo 
c 
      write(8,*) 'Contour Variables:'
      write(8,'(20(/,a50))') (title_var(i),i = 1, num_var-idem-1)
       if(j_cont_type.eq.2) then
         if(num_var-idem-1.eq.1) ivec = 1 
         if(num_var-idem-1.eq.2) ivec = 2
         if(num_var-idem-1.eq.3) ivec = 3
       else
           ivec = 0
       endif
c end do loop J_vtk = 1   
      else 
c for j_vtk gt 1, read 1 line for simulation time

       open(unit=9,file=input_data_file,status='old')
       open(unit=10,file=contour_file,status='unknown')  

       read(9,'(a400)') dum_zone1
c
       do i = 1, total_char-20
        if(dum_zone1(i:i+3).eq.'time') then
         ichar1 = i +5
         read(dum_zone1(ichar1:ichar1+13),'(1g13.04)') timez
        endif    
       enddo
      endif 
C  VTK format  header    
      write(10,186) timez, neq
186   format('# vtk DataFile Version 2.0',/,'FEHM VTK model, time:',

     &  1x,1p,g14.6,/,'ASCII',/,'DATASET UNSTRUCTURED_GRID',/,
     & 'POINTS',i9,1x,'double')
      
c gaz 021022 just allocate memory and read coor,elem for first tec file
      if(j_vtk.eq.1) then
c
      if(.not.allocated(x)) then
       allocate(x(neq))
       allocate(y(neq))
       allocate(node_i(neq))
       allocate(z(neq))
      endif
c       
c read variable list - FEHM output : node, x, y, z
      if(idem.eq.2) then
       num_tot_i = num_var-3
c       if(izid_position.ne.0) izid_position = izid_position -3
      else
       num_tot_i = num_var-4
c       if(izid_position.ne.0) izid_position = izid_position -4
      endif  
c gaz 031022 reallocate at j = 1 for sca,vel,con   
      if(j_vtk.eq.1) then
       if(allocated(temp)) deallocate(temp)
       allocate(temp(neq,num_tot_i))
      endif
      if(.not.allocated(zid_save).and.izid_position.ne.0)
     & allocate(zid_save(neq))
c read coordinate and variable data from tecplot file      
      do i = 1,neq
       if(idem.eq.2) then
        read(9,*) x(i), y(i), node_i(i), (temp(i,ii), ii = 1, num_tot_i)
       else
        read(9,*) 
     &    x(i), y(i), z(i), node_i(i), (temp(i,ii), ii = 1, num_tot_i)  
       endif
       if(izid_position.ne.0) zid_save(i) = temp(i,izid_position)
      enddo
      if(.not.allocated(idum)) allocate(idum(num_elem,ns))
      do i = 1,num_elem
       read(9,*) (idum(i,j),j = 1,ns)
      enddo
c write coordinates to vtk file
      do i = 1,neq
      if(abs(x(i)).lt.coor_tol) x(i) = 0.0d0
      if(abs(y(i)).lt.coor_tol) y(i) = 0.0d0
       if(idem.eq.2) then
        z(i) = 0.0d0
       else
        if(abs(z(i)).lt.coor_tol) z(i) = 0.0d0
       endif
       write(10,'(1p,3g14.6,i8,15g14.6)') x(i),y(i),z(i)
      enddo
c write element type data to vtk file and the total numbers to 
c define the element list(ns for each element and node numbers for each element      
      write(10,187) num_elem,(ns+1)*num_elem
187   format('CELLS', i9, i9)    
c need to start nodes at "0" not "1"  
      do  i = 1,num_elem
      write(10,'(9(1x,i9))') ns,(idum(i,inode)-1,inode=1,ns)
      enddo
      write(10,188) num_elem
188   format('CELL_TYPES', i9)  
      write(10,189) (ivtk_elem, i = 1,num_elem)
189   format(8(1x,i9)) 
c now write out variables in vtk format
      write(10,190) neq
      if(ivec.eq.0) then
c write out scalar name
      do ii = 1, num_tot_i
      write(10,191) title_var(ii)
c gaz 122722 inlarged the format
      write(10,192)        
        write(10,'(1p,10(1x,g17.6))') (temp(i,ii), i = 1, neq)
      enddo
      else
c write out velocity vectors
      
       write(10,291) 'Velocity_(m/s)'
       if(idem.eq.3)then
        do i = 1, neq 
         write(10,'(1p,3(1x,g16.6))') 
     &      (temp(i,ii), ii = num_tot_i-ivec+1, num_tot_i)
        enddo
       else
         do i = 1, neq 
          write(10,'(1p,3(1x,g16.6))') 
     &      (temp(i,ii), ii = num_tot_i-ivec+1, num_tot_i), 0.d00
        enddo
       endif
      endif
      
190   format('POINT_DATA',1x,i9)   
191   format('SCALARS',1x,a40,1x,'double 1')
291   format('VECTORS',1x,a30,1x,'double')      
192   format('LOOKUP_TABLE default')
c write scalar table with node numbers to all files  
      
c      write(10,190) neq
c write out scalar name
      write(10,191) 'Node_Num'
      write(10,192) 
      write(10,'(1p,10(1x,i14))') (i, i = 1, neq)      
c gaz 021322 added delete for tec files      
      close(9)
      close(10)
      else
c j_vtk gt 1 
c read coordinate and variable data from tecplot file          
c       read(9,'(a25,f13.1)') dum_zone3(1:25),timez (done already)
c       write(*,*)
c     &   'simulation time for this file = ', timez
       write(8,*)
     &   'simulation time for this file = ', timez     
C  VTK format  header    
c       write(10,186) timez, neq (done already)
c write coordinates to vtk file
      do i = 1,neq
      if(abs(x(i)).lt.coor_tol) x(i) = 0.0d0
      if(abs(y(i)).lt.coor_tol) y(i) = 0.0d0
       if(idem.eq.2) then
        z(i) = 0.0d0
       else
        if(abs(z(i)).lt.coor_tol) z(i) = 0.0d0
       endif
       write(10,'(1p,3g14.6,i8,15g14.6)') x(i),y(i),z(i)
      enddo
c write element type data to vtk file     
      write(10,187) num_elem,(ns+1)*num_elem 
c need to start nodes at "0" not "1"  
      do  i = 1,num_elem
      write(10,'(9(1x,i9))') ns,(idum(i,inode)-1,inode=1,ns)
      enddo
      write(10,188) num_elem 
      write(10,189) (ivtk_elem, i = 1,num_elem)
c read variable information 
      if(izid_position.eq.0) then
       do i = 1,neq
        read(9,*) kk, (temp(kk,ii), ii = 1, num_tot_i)
       enddo
      else 
       do i = 1,neq
        read(9,*) kk, (temp(kk,j), j = 1, izid_position-1),
     &   (temp(kk,ii), ii = izid_position+1,num_tot_i)
        temp(kk,izid_position) = zid_save(kk)
       enddo
      endif
c now write out variables in vtk format
      write(10,190) neq
c write out scalar name
      if(ivec.eq.0) then
c gaz 122722 enlarged format
       do ii = 1, num_tot_i
        write(10,191) title_var(ii)
        write(10,192) 
        write(10,'(1p,10(1x,g17.6))') (temp(i,ii), i = 1, neq)
       enddo
      else
c write out velocity vectors
       aminvx = -1.e20
       aminvx = 1.e20
       aminvy = -1.e20
       aminvy = 1.e20
       aminvz = -1.e20
       aminvz = 1.e20
       write(10,291) 'velocity'
       do i = 1, neq 
        if(idem.ne.3) then
        write(10,'(1p,3(1x,g14.6))') 
     &      (vscale*temp(i,j), j = num_tot_i-ivec+1, num_tot_i),0.d00  
        aminvx = min(aminvx,temp(i,num_tot_i-ivec+1))
        amaxvx = max(amaxvx,temp(i,num_tot_i-ivec+1))
        aminvy = min(aminvy,temp(i,num_tot_i-ivec+2))
        amaxvy = max(amaxvy,temp(i,num_tot_i-ivec+2))  
        else if(idem.eq.3) then
         write(10,'(1p,3(1x,g14.6))') 
     &      (vscale*temp(i,j), j = num_tot_i-ivec+1, num_tot_i)  
         aminvz = min(aminvz,temp(i,num_tot_i-ivec+3))
         amaxvz = max(amaxvz,temp(i,num_tot_i-ivec+3))  
        endif
        enddo
      endif  
c      if(ivec.ne.0) then
c       write (*,*) 'file number = ', ifile       
c       write (*,'(1p,6(1x,g14.6))')  aminvx, amaxvx,  aminvy, amaxvy,
c     &  aminvz, amaxvz
c      endif
c write scalar table with node numbers to all files  

c      write(10,190) neq
c write out scalar name
      write(10,191) 'Node_Num'
      write(10,192) 
      write(10,'(1p,10(1x,g14.6))') (i, i = 1, neq)
      close(9)      
      close(10)
c end of loop for vtk output files  
      endif
      enddo
c end of loop for sca,vec,con countour types   
100   continue      
      enddo
c end of iflg if block         
      endif   
      return 
      end
