		! This function should get called ONCE to initialize global silo variables
		! to avoid restepping through same code and creating mulitple arrays
	subroutine init_silo()
		use silo
		use avsio, only : iocord, iogeo, iokd, geoname, iodual, iogdkm,
     &     iogrid
		use comai
		use combi, only : corz      
		use comdi
		use comchem, only : cpntnam
		use comrxni, only : rxn_flag
		use comriv, only : iriver
		use davidi
		use silo
		
		implicit none
		
		include "silo.inc"
					
		integer i, j, lu, ifdual, maxtitle, mout, length, ic1, ic2
		integer il, open_file, nelm(ns_in), iocord_temp
		integer icord1, icord2, icord3, neq_read, neq_write, istart, iend
		parameter(maxtitle = 22)
		character*3 dls
		character*5 char_type, dual_char
		character*15 nform
		character*300 temp_string
		character*42 title(maxtitle), units(maxtitle), pstring
		character*42, allocatable :: title_kd(:)
		character*600 print_title, vstring
		real*8 perm_fac
		parameter (perm_fac=1.d-6)

		! globally allocate size of globally used arrays 
		!( you must re-create the same mesh and connectivity each timestep )
		silo_meshtype = ''
		allocate(x_cord(n0))
      	allocate(y_cord(n0))
      	allocate(z_cord(n0))
      	allocate(silo_cv(nei_in*ns_in))
      	
      	! Set case meshtype(shape)
		select case (ns_in)
		case (5,6,8)
			silo_meshtype = 'brick'
			ndims = 3
			shapetype(1) = DB_ZONETYPE_HEX
			shapesize(1) = 8
		case (4)
		if (icnl .eq. 0) then 
		 	silo_meshtype = 'tetrahedron'
		 	ndims = 3
		 	shapetype(1) = DB_ZONETYPE_TET
		 	shapesize(1) = 4
		else
		 	silo_meshtype = 'quadrilateral'
		 	ndims = 2
			shapetype(1) = DB_ZONETYPE_QUAD
		 	shapesize(1) = 4
		end if
		case (3)
			silo_meshtype = 'triangle'
			ndims = 2
			shapetype(1) = DB_ZONETYPE_TRIANGLE
			shapesize(1) = 3
		case (2)
			silo_meshtype = 'lineseg'
			ndims = 2
			shapetype(1) = DB_ZONETYPE_BEAM
			shapesize(1) = 2
		case (0)
		end select
		
		! Set case cord		
		select case (icnl)
		case (1, 4)
		  icord1 = 1
		  icord2 = 2
		  icord3 = 1
		  x_write = .true.
          y_write = .true.
		case (2, 5)
		  icord1 = 1
		  icord2 = 3
		  icord3 = 2
		  x_write = .true.
          z_write = .true.
		case(3, 6)
		  icord1 = 2
		  icord2 = 3
		  icord3 = 1
		  y_write = .true.
          z_write = .true.
		case default
		  icord1 = 1
		  icord2 = 3
		  icord3 = 1
		  x_write = .true.
          y_write = .true.
          z_write = .true.
		end select
									
		shapecounts(1) = nei_in
		nnodes = ndims * n0
		print *, "SILO TPL initizaled: "
		print *, "DIMENSIONS = ", ndims
		print *, "TOTAL NODES = ", nnodes
		print *, "MESHTYPE = ", silo_meshtype
		print *, ''

		! Set coordinates
		do j = 1, n0		
			do i = icord1, icord2, icord3			   		
				if (i .eq. 1) then					
					x_cord(j) = corz(j, i)
				else if (i .eq. 2) then
					y_cord(j) = corz(j, i)
				else if (i .eq. 3) then
					z_cord(j) = corz(j, i)
				end if
			end do
		end do
			
	return
	end
	
	
	! This function creates silo file depending on parameter (material, scalar, vector, con, etc..)
	! will get called once per time step for time dependant variables
	subroutine create_silo(silo_func)	
		use silo		
		implicit none		
		include "silo.inc"
		integer i
		character*3 silo_func
		
		if(silo_func .eq. 'mat') then
			! Create the Silo mat file		
			silofile = 'materials.silo'
			print *, 'Creating SILO file: ', silofile
			silo_ierr = dbcreate(silofile, 15, DB_CLOBBER, DB_LOCAL, "Unstructured 3d mesh", 20, DB_HDF5, silo_dbfile)
			
		else if(silo_func .eq. 'sca') then
			! Create the Silo sca file
			silofile = ''
			i = silo_cycl
			write(silofile,'(a,i4.4,a)')"scalar",i,".silo"		
			print *, 'Creating SILO file: ', silofile
			silo_ierr = dbcreate(silofile, 15, DB_CLOBBER, DB_LOCAL, "Unstructured 3d mesh", 20, DB_HDF5, silo_dbfile)			
		endif

		if(silo_dbfile.eq.-1) then
			  write (6,*) 'Could not create Silo file!\n'
			  stop
		endif	
	
	return
	end
	
	
	! creates silo file for material properties
	subroutine write_silo_mat(ifdual)
	! Author: Shane McKinney, Dylan Harp
	  use avsio, only : iocord, iogeo, iokd, geoname, iodual, iogdkm,
     &     iogrid
      use comai
      use combi, only : corz      
      use comdi
      use comchem, only : cpntnam
      use comrxni, only : rxn_flag
      use comriv, only : iriver
      use davidi
      use silo

      implicit none
      
      include "silo.inc"
                    
      integer i, j, lu, ifdual, maxtitle, mout, length, ic1, ic2
      integer il, open_file, nelm(ns_in), iocord_temp
      integer icord1, icord2, icord3, neq_read, neq_write, istart, iend
      parameter(maxtitle = 22)
      character*3 dls
      character*5 char_type, dual_char
      character*15 nform
      character*300 temp_string
      character*42 title(maxtitle), units(maxtitle), pstring
      character*42, allocatable :: title_kd(:)
      character*600 print_title, vstring
      real*8 perm_fac
      parameter (perm_fac=1.d-6)
      logical :: xon = .false., yon = .false., zon = .false.      
      !(SILO)	Must dynamically allocate these to fucntion for complex geometry(more than 1 type of polyhedron)
      real*8 x_perm(size(pnx)), y_perm(size(pny)), z_perm(size(pnz))
      real*8 thermal_cond_x(n0), thermal_cond_y(n0), thermal_cond_z(n0), porosity(n0)
      real*8 rel_perm(n0), rock_spec(n0), cap_pres(n0), rock_bulk(n0)      		
      	
      	call create_silo('mat')
      	
		if (ifdual .eq. 1) then
			istart = 1
			iend = neq - neq_primary
			if (icnl .eq. 0) then
			   iocord = 3
			else
			   iocord = 2
			end if
		else
			istart = 1
			iend = neq_primary
		end if
      	  		
		! Set read counter		
		if(iriver.ne.2.and.gdkm_flag.ne.1) then
		  neq_read = neq
		else
		  neq_read = neq_primary
		endif		

		! Set connectivity (must read it from file?)
        if (iogeo .eq. 1 .and. altc(1:4) .eq. 'silo') then
            il = open_file(geoname,'old')
            read(il,*) i
            if (i .ne. neq_primary) backspace il
			do i = 1, neq_read
				read(il,*)
			end do
			icv = 1
			do i = 1, nei_in
				read (il,*) ic1,ic2,char_type,(nelm(j), j=1,ns_in)
				do kk = 1, ns_in
					silo_cv(icv) = nelm(kk)
					icv = icv + 1
				end do
			end do
			close (il)
		end if
		
		! Set arrays that require calculations
		do i = 1, iend
			if (idoff .ne. -1) then
				x_perm(i) = pnx(i)*perm_fac
				y_perm(i) = pny(i)*perm_fac
				z_perm(i) = pnz(i)*perm_fac
			end if  			
		end do
		
! Checking constructed silo arrays		
!		do ii = 1, size(x_cord)		
!			print *, "connectivity array: ", ii, silo_cv(ii)
!			print *, "coordinate array: ", ii, x_cord(ii), y_cord(ii),z_cord(ii)
!		end do
						
		
! New version connectivity silo call 		
		silo_err = dbputzl2(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
     &		nei_in*ns_in, 1, 0, 0, shapetype, shapesize, shapecounts, NSHAPETYPES, DB_F77NULL, silo_ierr)          
     
! Depricated silo connectivity call	
!     	silo_err = dbputzl(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
!     &		ns_in*nei_in, 1, shapesize, shapecounts, NSHAPETYPES, silo_ierr)

! This unstruct mesh write call does not match manual, see -->
! http://visit.ilight.com/svn/visit/trunk/src/tools/DataManualExamples/CreatingCompatible/fucd3d.f
	
     	! option list to add units to variables
		silo_err = dbmkoptlist(100, optlistid)
		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m", 4)
	 	silo_err = dbputum(silo_dbfile, "pointmesh", 9, ndims, x_cord, y_cord, z_cord, "X", 1, "Y", 1, "Z", 1,
     &		DB_FLOAT, iend, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, optlistid, silo_ierr)
 
     
     	! Write permeability    
		if (idoff .ne. -1) then
			silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m**2", 4)    
			if (x_write) then
				silo_ierr = dbputuv1(silo_dbfile, "permeability_x", 14, "pointmesh", 9, x_perm, iend,
     &			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		end if
     		if (y_write) then
				silo_ierr = dbputuv1(silo_dbfile, "permeability_y", 14, "pointmesh", 9, y_perm, iend,
     &			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		end if
     		if (z_write) then
     			silo_ierr = dbputuv1(silo_dbfile, "permeability_z", 14, "pointmesh", 9, z_perm, iend,
     &			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		end if
     	end if     	
     	
     	! Conductivity will be written
		if (ico2 .gt. 0 .or. ice .ne. 0) then
			   silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "W/m*K", 5)
			   if (x_write) then
			   	   silo_ierr = dbputuv1(silo_dbfile, "thermal_conductivity_x", 22, "pointmesh", 9, thx, iend,
     &	    	   DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
			   end if
			   if (y_write) then
			   	   silo_ierr = dbputuv1(silo_dbfile, "thermal_conductivity_y", 22, "pointmesh", 9, thy, iend,
     &	    	   DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
			   end if
			   if (z_write) then
			   	   silo_ierr = dbputuv1(silo_dbfile, "thermal_conductivity_z", 22, "pointmesh", 9, thz, iend,
     &	    	   DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
			   end if
		end if	
		
		! Porosity, bulk density, and specific heat will be written
		silo_ierr = dbputuv1(silo_dbfile, "porosity", 8, "pointmesh", 9, ps, iend,
     &  DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, DB_F77NULL, silo_ierr)
     
     	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/m**3", 7)
     	silo_ierr = dbputuv1(silo_dbfile, "rock_bulk_desnity", 17, "pointmesh", 9, denr, iend,
     &  DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     
     	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "MJ/kg*K", 7)
     	silo_ierr = dbputuv1(silo_dbfile, "rock_specific_heat", 18, "pointmesh", 9, cpr, iend,
     &  DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     
     	! Capillary pressure will be written
		if (rlp_flag .eq. 1) then
			silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "MPa", 3)
			silo_ierr = dbputuv1(silo_dbfile, "capillary_pressure", 18, "pointmesh", 9, pcp, iend,
     &    	DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
		end if
            
! rlp and cap model flags will be written if rlp_flag .eq. 1            
!            if (rlp_flag .eq. 1) then            	 
!                silo_ierr = dbputuv1(silo_dbfile, "relative_permeability_model", 27, "pointmesh", 9, rel_perm, iend,
!     &	    	DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, DB_F77NULL, silo_ierr)
!            end if
!            
! Kd will be written a1adfl(i,itrc(j))
!            if(iccen .eq. 1 .and. iokd .eq. 1 .and. rxn_flag .eq.0) then
!               silo_ierr = dbputuv1(silo_dbfile, "capillary_pressure_model", 24, "pointmesh", 9, a1adfl, iend,
!     &	       DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, DB_F77NULL, silo_ierr)
!            end if


		
! This silo call does not work yet matches input in the manual 4.10.2
!   	    silo_err = dbputum(silo_dbfile, "pointmesh123", 12, ndims, x_cord, y_cord, z_cord, "X", 1, "Y", 1, "Z", 1,
!     &		ns_in*nei_in, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, DB_DOUBLE, DBOPT_PHZONELIST, silo_ierr)	

!This will write mesh without connectivity
!		silo_err = dbputpm(silo_dbfile, "pointmesh", 9, ndims, x_cord, y_cord, z_cord, iend, DB_DOUBLE, DB_F77NULL, silo_ierr)
	
		silo_ierr = dbfreeoptlist(optlistid)
		silo_ierr = dbclose(silo_dbfile)
	return
	end
	
	
	! writes scalar variables to file for each time step
	subroutine write_silo_s(icall, neq, nscalar, ifdual, iriver2)		      
	  use avsio
      use comai, only : altc, days, grav, iadif, icnl, ico2, idof, 
     &     ichead, ihead, nei_in, ns_in, phi_inc, istrs, ivf,
     &     neq_primary, rho1grav, ifdm_elem, igrav
!     &     neq_primary, rho1grav, ifdm_elem, i_subsid, igrav
      use combi, only : corz, izonef, nelm, nelmdg, sx1
      use comci, only : rolf, rovf
      use comdi
      use comfi, only : pci
      use comflow, only : a_axy, a_vxy
      use comfem
      use comii, only : crl
      use comwt, only : sattol, head_id, rlptol
      use comsi
      use davidi
!     RJP 1/12/07 added following
      use comco2
      use comriv, only : npoint_riv, nnelm_riv, nelm_riv
      use silo
      
      implicit none
      
      include "silo.inc"
      
      integer maxscalar
      parameter (maxscalar = 37)
      integer neq,nscalar,lu,ifdual,icall,open_file,offset,iriver2
      integer i,j,iolp,iovp,nout,iz,iendz,il,idz, i1, i2, index, iaxy, k
      integer size_head, size_pcp, istart, iend, ic1, ic2, length, nadd
      integer icord1, icord2, icord3, ns_in0, irivp, iocord_tmp 
      integer, allocatable :: nelm2(:)
      real*8 hdum, sdum, px, py, pz, flxdum
      real*8 pdum, tdum, rolconv, dumconv, dumconv1
      character*80 title(2*maxscalar+3)
      character*150 :: tecstring = ''
      character*150 :: tecstring_riv = ''
      character*500 string
      character*20 vstring
      character*43 tstring
      character*5 char_type
      character*3 dls      
!(SILO)	Must dynamically allocate these to fucntion for complex geometry(more than 1 type of polyhedron) ?
!	  Over allocation below  ?    
      real*8 pressure(n0), silo_temp(n0)

      
      call create_silo('sca')
      
      iocord_tmp = iocord
      if (iogdkm .ne. 0 .and. ifdual .ne. 0) then
      	
      	  !Output for gdkm nodes
         istart = neq_primary + 1
         iend = neq
         offset = 0
         nadd = 0
         if (icnl .eq. 0) then
            iocord = 3
         else
            iocord = 2
         end if
      else if (ifdual .ne. 0)then
         istart = neq + 1
         iend = neq * 2
         nadd = nelm(neq+1)-neq-1
         offset = neq
      else 
         if (iriver2 .ne. 0) then
         	 
         	 ! Output for river/well nodes         
            istart = neq_primary + 1
            iend = neq_primary + npoint_riv
            nadd = 0
            offset = 0
            if(iriver2.eq.2) then
               irivp = 2
               ns_in0 = ns_in
               ns_in = 2
            endif
         else
            istart = 1
            iend = neq_primary
            nadd = 0
            offset = 0
            irivp = 0
         end if
      endif
      
      if (iozone .ne. 0 ) then
         iendz = nsurf
      else 
         iendz = 1
         idz = iozone
      end if
        
      	silo_time = days
		silo_err = dbmkoptlist(100, optlistid)		
		silo_err = dbaddiopt(optlistid, DBOPT_CYCLE, silo_cycl)		
		silo_err = dbaddropt(optlistid, DBOPT_TIME, silo_time)		
				
		! Rewrite mesh and connectivity
		silo_err = dbputzl2(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
     &	nei_in*ns_in, 1, 0, 0, shapetype, shapesize, shapecounts, NSHAPETYPES, optlistid, silo_ierr)
     
     	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m", 1)
	 	silo_err = dbputum(silo_dbfile, "pointmesh", 9, ndims, x_cord, y_cord, z_cord, "X", 1, "Y", 1, "Z", 1,
     &	DB_FLOAT, iend, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, optlistid, silo_ierr)
     
     	! Boolean logic to write correct variables
     		! Pressure liquid
     		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "MPa", 3)
            if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
               if (size(pcp) .ne. 1) then
					do i = 1, n0
						pressure(i) = phi(i)-pcp(i)
					end do               	   
               	    silo_ierr = dbputuv1(silo_dbfile, "Liquid_Pressure", 15, "pointmesh", 9, pressure, n0,
     &  			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)                  
               else
          			do i = 1, n0
						pressure(i) = phi(i)-phi_inc
					end do
               	    silo_ierr = dbputuv1(silo_dbfile, "Liquid_Pressure", 15, "pointmesh", 9, pressure, n0,
     &  			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)                  
               end if
            end if
            
            ! Pressure vapor
            if (iopressure .eq. 1 .and. iovapor .eq. 1) then
            		do i = 1, n0
						pressure(i) = phi(i)
					end do
            		silo_ierr = dbputuv1(silo_dbfile, "Vapor_Pressure", 14, "pointmesh", 9, pressure, n0,
     &  			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
               if (iadif .eq. 1) then
               	   	do i = 1, n0
						pressure(i) = phi(i)-pci(i)
					end do
               	    silo_ierr = dbputuv1(silo_dbfile, "Water_vapor_Pressure", 20, "pointmesh", 9, pressure, n0,
     &  			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
               end if
            end if
            
            ! Temperature            
            if (iotemperature .eq. 1) then
            	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "deg_C", 5)
            	silo_ierr = dbputuv1(silo_dbfile, "Temperature", 20, "pointmesh", 9, t, n0,
     & 			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
            end if
                        
            ! Density Liquid          
            if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
            	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/m**3", 7)
            	silo_ierr = dbputuv1(silo_dbfile, "Liquid_Density", 14, "pointmesh", 9, rolf, n0,
     & 			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)            	
            end if
            
            ! Density Vapor
            if (iodensity .eq. 1 .and. iovapor .eq. 1) then
            	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/m**3", 7)
            	silo_ierr = dbputuv1(silo_dbfile, "Vapor_Density", 13, "pointmesh", 9, rovf, n0,
     & 			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)            	
            end if
            ! Calculate fluxes and set silo_temp array to correct flux value
            if (ioflx .eq. 1) then
				do i = istart, iend
					! Liquid Flux
					if (ioflx .eq. 1 .and. ioliquid .eq. 1) then
					   if (.not. net_flux) then
						  iaxy = nelmdg (i) - (neq + 1) + nadd
						  if (vol_flux) then
							 if (sx1(i) .gt. 0.) then
								flxdum = a_axy(iaxy) / sx1(i)
							 else
								flxdum = 0.d0
							 end if
						  else
							 flxdum = a_axy(iaxy)
						  end if
						  silo_temp(i) = flxdum            	            	
					   else
						  i1 = nelm(i) + 1
						  i2 = nelm(i+1)
						  flxdum = 0.
						  do index = i1, i2
							 iaxy = index - neq - 1 + nadd
							 if(a_axy(iaxy).gt.0.) then
								flxdum = flxdum + a_axy(iaxy)
							 end if
						  end do
						  if (vol_flux) then
							 if (sx1(i) .gt. 0.) then
								flxdum = flxdum / sx1(i)
							 else
								flxdum = 0.d0
							 end if
						  end if
						  silo_temp(j) = flxdum            	           	
					   end if              
					end if 
				end do			
			end if
			! Writing liquid flux values
        	if (vol_flux) then
        		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/m**3", 7)
        	else
        		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/s", 7)
        	end if
        	if (ioflx .eq. 1 .and. .not. net_flux) then
					silo_ierr = dbputuv1(silo_dbfile, "Liquid_Flux", 11, "pointmesh", 9, silo_temp, n0,
     &	 			DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		else if (ioflx .eq. 1 .and. net_flux) then
					silo_ierr = dbputuv1(silo_dbfile, "Net_Liquid_Flux", 15, "pointmesh", 9, silo_temp, n0,
     & 				DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr) 		
            end if            
            ! Calculating vapor flux            						
			do i = istart, iend            							
				if (ioflx .eq. 1 .and. iovapor .eq. 1) then
				   iaxy = nelmdg (i) - (neq + 1) + nadd
				   silo_temp(i) = a_vxy(iaxy)
				end if 
			end do
			! Vapor Flux            
            if (ioflx .eq. 1 .and. iovapor .eq. 1) then
            		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "kg/s", 7)
					silo_ierr = dbputuv1(silo_dbfile, "Vapor_Flux", 10, "pointmesh", 9, silo_temp, n0,
     & 				DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
     		end if
     		! Capillary Pressure     		
            if (iocapillary .eq. 1) then
            		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "MPa", 3)
					silo_ierr = dbputuv1(silo_dbfile, "Capillary_Pressure", 18, "pointmesh", 9, pcp, n0,
     & 				DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, optlistid, silo_ierr)
            end if     		
     		
     		
     	silo_ierr = dbfreeoptlist(optlistid)

     	silo_cycl = silo_cycl + 1
     	silo_ierr = dbclose(silo_dbfile)
	return
	end