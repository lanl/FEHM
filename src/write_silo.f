
	
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

      	allocate(silo_cv(nei_in*ns_in))
      	
		select case (icnl)
		case (1, 4)
			x_write = .true.
			y_write = .true.
		case (2, 5)
			x_write = .true.
			z_write = .true.
		case(3, 6)
			y_write = .true.
			z_write = .true.
		case default
			x_write = .true.
			y_write = .true.
			z_write = .true.
		end select      	
      	
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
									
		shapecounts(1) = nei_in
		nnodes = ndims * n0
		print *, "SILO TPL initizaled: "
		print *, "DIMENSIONS = ", ndims
		print *, "TOTAL NODES = ", nnodes
		if (neq_primary .ne. n0) then
			print *, "GDKM NODES = ", nnodes/ndims-neq_primary
		end if
		print *, "MESHTYPE = ", silo_meshtype
		print *, ''
		
			
	return
	end	
		
	! This function creates silo file depending on parameter (material, scalar, vector, con, etc..)
	! will get called once per time step for time dependant variables
	subroutine create_silo(silo_func)	
		use avsio, only : timec_string
		use comai, only : days
		use silo		
		implicit none		
		include "silo.inc"
		integer i
		character*3 silo_func
		
		if(silo_func .eq. 'mat') then
			! Create the Silo mat file		
			silofile = 'materials.silo'
			print *, 'Creating SILO file: ', trim(silofile)
			silo_ierr = dbcreate(silofile, 14, DB_CLOBBER, DB_LOCAL, "Unstructured 3d mesh", 20, DB_HDF5, silo_dbfile)
			
		else if(silo_func .eq. 'sca') then
			! Create the Silo sca file
			silofile = ''
			i = silo_cycl
			write(silofile,'(a,i4.4,a)')"scalar",i,".silo"		
			print *, 'Creating SILO file: ', trim(silofile), '  at ', trim(timec_string)
			silo_ierr = dbcreate(silofile, 15, DB_CLOBBER, DB_LOCAL, "Unstructured 3d mesh", 20, DB_HDF5, silo_dbfile)			
		
		else if(silo_func .eq. 'con') then
			! Create the Silo con file
			silofile = ''
			i = silo_cycl
			write(silofile,'(a,i4.4,a)')"concen",i,".silo"		
			print *, 'Creating SILO file: ', trim(silofile), '  at ', trim(timec_string)
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
		use avsio, only : iocord, iogeo, iokd, geoname, iodual, iogdkm,
     &  iogrid
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
		

		real*8 x_perm(n0), y_perm(n0), z_perm(n0)
		real*8 thermal_cond_x(n0), thermal_cond_y(n0), thermal_cond_z(n0), porosity(n0)
		real*8 rel_perm(n0), rock_spec(n0), cap_pres(n0), rock_bulk(n0)     		
      	
      	call create_silo('mat')     	      	
      	  		
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
				if(char_type(1:3) .eq. 'tet') then
					do kk = 1, 4						
						silo_cv(icv) = nelm(kk)
						icv = icv + 1					
					end do
				else if(char_type(1:3) .eq. 'tri') then
					do kk = 1, 3						
						silo_cv(icv) = nelm(kk)
						icv = icv + 1					
					end do
				else
					do kk = 1, ns_in						
						silo_cv(icv) = nelm(kk)
						icv = icv + 1					
					end do
				end if
			end do
			close (il)
		end if
		
		! Set arrays that require calculations
		do i = 1, n0
			if (idoff .ne. -1) then
				x_perm(i) = pnx(i)*perm_fac
				y_perm(i) = pny(i)*perm_fac
				z_perm(i) = pnz(i)*perm_fac
			end if  			
		end do
						
		
		! New version connectivity silo call 		
		silo_err = dbputzl2(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
     &	nei_in*ns_in, 1, 0, 0, shapetype, shapesize, shapecounts, NSHAPETYPES, DB_F77NULL, silo_ierr)          
     
	
     	! option list to add units to variables
		silo_err = dbmkoptlist(100, optlistid)
		silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m", 4)
		
		!creates mesh
	 	silo_err = dbputum(silo_dbfile, "pointmesh", 9, ndims, corz(:,1), corz(:,2), corz(:,3), "X", 1, "Y", 1, "Z", 1,
     &	DB_DOUBLE, neq_primary, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, optlistid, silo_ierr)
 

     	! Write permeability 
		if (idoff .ne. -1) then	
		
			if (x_write) then							
				call silo_put_unstruc2("permeability_x", "m**2", x_perm)
     		end if
     		
     		if (y_write) then
    			call silo_put_unstruc2("permeability_y", "m**2", y_perm)
     		end if
     		
     		if (z_write) then
				call silo_put_unstruc2("permeability_z", "m**2", z_perm)
     		end if
     		
     	end if     	
     	
     	! Conductivity will be written
		if (ico2 .gt. 0 .or. ice .ne. 0) then
			
			   if (x_write) then     
					call silo_put_unstruc2("thermal_conductivity_x", "W/m*K", thx)
			   end if
			   
			   if (y_write) then
					call silo_put_unstruc2("permeability_z", "W/m*K", thy)
			   end if
			   
			   if (z_write) then
					call silo_put_unstruc2("permeability_z", "W/m*K", thz)
			   end if
			   
		end if	
		
		! Porosity, bulk density, and specific heat will be written
		call silo_put_unstruc2("porosity", "none", ps)

		call silo_put_unstruc2("rock_bulk_density", "kg/m**3", denr)
     
		call silo_put_unstruc2("rock_specific_heat", "MJ/kg*K", cpr)
     
     	! Capillary pressure will be written
		if (rlp_flag .eq. 1) then
			call silo_put_unstruc2("capillary_pressure", "MPa", pcp)
		end if
 
!MISSING		
! 			rlp and cap model flags will be written if rlp_flag .eq. 1            
!            if (rlp_flag .eq. 1) then            	 
!                silo_ierr = dbputuv1(silo_dbfile, "relative_permeability_model", 27, "pointmesh", 9, rel_perm, iend,
!     &	    	DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, DB_F77NULL, silo_ierr)
!            end if
!            
! 			Kd will be written a1adfl(i,itrc(j))
!            if(iccen .eq. 1 .and. iokd .eq. 1 .and. rxn_flag .eq.0) then
!               silo_ierr = dbputuv1(silo_dbfile, "capillary_pressure_model", 24, "pointmesh", 9, a1adfl, iend,
!     &	       DB_F77NULLSTRING, 0, DB_DOUBLE, DB_BLOCKCENT, DB_F77NULL, silo_ierr)
!            end if


		! Clear option list and close file
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
      allocate( silo_pres(n0), silo_flx(n0), silo_hdum(n0), silo_co2_pgas(n0), silo_du(n0),
     &	silo_dv(n0), silo_dw(n0), silo_sdum(n0), silo_vflx(n0), silo_fw(n0),
     &	silo_fl(n0), silo_fg(n0), silo_yc(n0), silo_co2_prop(n0), silo_ices(n0),
     &	silo_vpres(n0), silo_wvpres(n0) )
   



     	silo_cycl = icall
     	call create_silo('sca')	  
     	
		! Set arrays that require calculations
		do i = 1, n0
			
			! Displacement
			if (iodisp .eq. 1 .and. idisp_rel .ne. 0) then
				silo_du(i) = du(i)-du_ini(i)
				silo_dv(i) = dv(i)-dv_ini(i)
				if (icnl .eq. 0) then
					silo_dw(i) = dw(i)-dw_ini(i)
				end if
			else if (iodisp .eq. 1) then
				silo_du(i) = du(i)
				silo_dv(i) = dv(i)
				if (icnl .eq. 0) then
					silo_dw(i) = dw(i)
				end if 
			end if
			
			! Liquid Pressure
			if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
				if (size(pcp) .ne. 1) then
					silo_pres(i) = phi(i)-pcp(i)
				else
					silo_pres(i) = phi(i)-phi_inc
				end if
			end if
			
			! Vapor Pressure
			if (iopressure .eq. 1 .and. iovapor .eq. 1) then
				silo_vpres(i) = phi(i)
				if (iadif .eq. 1) then
					silo_wvpres(i) = phi(i)-pci(i)
				end if  				
			end if
		end do
		
		! Calculate conditional values and add to silo array (DP)
		if (ioflx .eq. 1 .or. iohead .eq. 1 .or. iosaturation .eq. 1 .or. ioco2 .eq. 1) then
			do i = 1, n0
				
				! Calc Saturations, Density, and phase state zeroes
				if (ioco2 .eq. 1) then
					if (ps(i) .le. 0) then
						sdum = 0.0d0
						silo_fw(i) = sdum
						silo_fl(i) = sdum
						silo_fg(i) = sdum
						silo_yc(i) = sdum
						silo_co2_prop(i) = sdum
						silo_co2_pgas(i) = sdum
						silo_ices(i) = sdum
					else
						silo_fw(i) = fw(i)
						silo_fl(i) = fl(i)
						silo_fg(i) = fg(i)
						silo_yc(i) = yc(i)
						silo_co2_prop(i) = co2_prop(i)
						silo_co2_pgas(i) = co2_prop(9*neq_primary+i)
						silo_ices(i) = ices(i)
					end if
				end if
				
				! Calc Liquid Flux
				if (ioflx .eq. 1) then
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
					   end if              
					end if
					silo_flx(j) = flxdum
				end if
				
			   ! Calc Hydrualic Head (UNTESTED)
			   if (iohead .eq. 1) then
				   if (ichead .eq. 1) then
					  if (ps(i) .le. 0.) then
						 hdum = 0.d0
					  else
						 call headctr(4, i   ,pho(i), hdum)
					  end if					  
				   else if (ihead.eq.1 .and. size_head .ne. 1) then
					  if (ps(i) .le. 0.) then
						 hdum = 0.d0
					  else
						 hdum = head(i) 
						 if (irdof .ne. 13 .or. ifree .ne. 0) then
							if (s(i).lt.sattol+rlptol) hdum = head_id
						 endif 
					  end if					  
				   end if
				   silo_hdum(i) = hdum
				end if
			   
			   ! Calc Saturation
				if (iosaturation .eq. 1) then
					  if (ps(i) .le. 0.) then
						 sdum = 0.d0
					  else if (irdof .ne. 13 .or. ifree .ne. 0) then
						  sdum = s(i)
						  !print *, "sdum:", sdum						  
					  else
						 sdum = 1.d0
					  end if
					  silo_sdum(i) = sdum
					  !print *, "ss", silo_sdum(i)
				end if
				
				! Calculating vapor flux 
				if (ioflx .eq. 1 .and. iovapor .eq. 1) then
				   iaxy = nelmdg (i) - (neq + 1) + nadd
				   silo_vflx(i) = a_vxy(iaxy)
				end if 		
				
			end do			
		end if		
		
		! Set timestep for multiple function calls of write_silo_s()
      	silo_time = days
		silo_err = dbmkoptlist(100, optlistid)		
		silo_err = dbaddiopt(optlistid, DBOPT_CYCLE, silo_cycl)		
		silo_err = dbaddropt(optlistid, DBOPT_DTIME, silo_time)		
				
		! Rewrite mesh and connectivity
		silo_err = dbputzl2(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
     &	nei_in*ns_in, 1, 0, 0, shapetype, shapesize, shapecounts, NSHAPETYPES, optlistid, silo_ierr)
     
     	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m", 1)
	 	silo_err = dbputum(silo_dbfile, "pointmesh", 9, ndims, corz(:,1), corz(:,2), corz(:,3), "X", 1, "Y", 1, "Z", 1,
     &	DB_DOUBLE, neq_primary, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, optlistid, silo_ierr)
     
		! Pressure liquid
		if (iopressure .eq. 1 .and. ioliquid .eq. 1) then
			call silo_put_unstruc2("Liquid_Pressure", "MPa", silo_pres)     
		end if
		
		! Pressure vapor
		if (iopressure .eq. 1 .and. iovapor .eq. 1) then
				call silo_put_unstruc2("Vapor_Pressure", "MPa", silo_vpres)        
		   if (iadif .eq. 1) then
				call silo_put_unstruc2("Water_Vapor_Pressure", "MPa", silo_wvpres)     				
		   end if
		end if
		
		! Temperature            
		if (iotemperature .eq. 1) then
			call silo_put_unstruc2("Temperature", "deg_C", t)
		end if                    

		! Density Liquid          
		if (iodensity .eq. 1 .and. ioliquid .eq. 1) then
			call silo_put_unstruc2("Liquid_Density", "kg/m**3", rolf)
		end if
		
		! Density Vapor
		if (iodensity .eq. 1 .and. iovapor .eq. 1) then
			call silo_put_unstruc2("Vapor_Density", "kg/m**3", rovf)
		end if
		
		! Displacements
		if (iodisp .eq. 1) then   				
		
			call silo_put_unstruc2("Displacement_X", "m", silo_du)
			
			call silo_put_unstruc2("Displacement_Y", "m", silo_dv)
			
			if (icnl .eq. 0) then
				call silo_put_unstruc2("Displacement_Z", "m", silo_dw)     	
			end if
			
		end if            	            	
		
		! Saturation (UNTESTED)           
		if (iosaturation .eq. 1) then        	
			call silo_put_unstruc2("Saturation", "none", silo_sdum)     
		end if
		
		! Hydrualic Head
		if (iohead .eq. 1) then
			call silo_put_unstruc2("Hydrualic_Head", "none", silo_hdum) 
		end if

		! Liquid flux values
		if (vol_flux) then
			silo_unit = "kg/m**3"
		else
			silo_unit = "kg/s"
		end if
		
		if (ioflx .eq. 1 .and. .not. net_flux) then
			call silo_put_unstruc2("Liquid_Flux", silo_unit, silo_flx)
		else if (ioflx .eq. 1 .and. net_flux) then
			call silo_put_unstruc2("Net_Liquid_Flux", silo_unit, silo_flx)     
		end if
		
		! Water Saturation, Super-Critical/Liquid CO2 Saturation, Gaseous CO2 Saturation, Dissolved CO2 Mass Fraction
		! CO2 Liquid Density (kg/m**3), CO2 Gas Density (kg/m**3), CO2 Phase State
		! (UNTESTED)
		if (ioco2 .eq. 1) then         		

			call silo_put_unstruc2("Water_Saturation", "none", silo_fw)
	  
			call silo_put_unstruc2("Super_Critical_Liquid_CO2_Saturation", "none", silo_fl)		 
 
			call silo_put_unstruc2("Gaseous_CO2_Saturation", "none", silo_fg)		 
 
			call silo_put_unstruc2("Dissolved_CO2_Mass_Fraction", "none", silo_yc)		
 
			call silo_put_unstruc2("CO2_Phase_State", "none", silo_ices)      				
			 
			call silo_put_unstruc2("CO2_Liquid_Desnsity", "kg/m**3", silo_co2_prop)	
			
			call silo_put_unstruc2("CO2_Gas_Density", "kg/m**3", silo_co2_pgas)		

		end if                             
		
		! Vapor Flux            
		if (ioflx .eq. 1 .and. iovapor .eq. 1) then
 			call silo_put_unstruc2("Vapor_Flux", "kg/s", silo_vflx)
		end if
		
		! Source (UNTESTED)            
		if (iosource .eq. 1) then
 			call silo_put_unstruc2("Source", "kg/s", sk)
		end if     		
		
		! Capillary Pressure     		
		if (iocapillary .eq. 1) then
 			call silo_put_unstruc2("Capillary_Pressure", "MPa", pcp)
		end if
		
		! Stress, excess shear, and strain (UNTESTED)				
		if (iostress .ne. 0) then
			
			if(iPlastic.eq.1) then
				call silo_put_unstruc2("Plastic_Strain", "none", pstrain)
			endif			
		    
			call silo_put_unstruc2("Stress_X", "MPa", str_x)     
			 
			call silo_put_unstruc2("Stress_Y", "MPa", str_y)     
			
			call silo_put_unstruc2("Stress_XY", "MPa", str_xy)     				
				
			if(icnl.eq.0) then
			 
				call silo_put_unstruc2("Stress_Z", "MPa", str_z)  
				 
				call silo_put_unstruc2("Stress_XZ", "MPa", str_xz)					
							
				call silo_put_unstruc2("Stress_YZ", "MPa", str_yz)
				
			end if  
			
			if(flag_excess_shear.eq.1) then  
			
				call silo_put_unstruc2("Youngs_Mod", "MPa", elastic_mod)     
			
				call silo_put_unstruc2("Excess_Shear", "MPa", excess_shear)          			
			
				call silo_put_unstruc2("Shear_Angle", "deg", shear_angle)     
			
			endif               
		end if
		
		! Volume Strain
		if (iostrain .eq. 1) then
			call silo_put_unstruc2("Volume_Strain", "none", vol_strain)
		end if

		deallocate( silo_pres, silo_flx, silo_hdum, silo_co2_pgas, silo_du,
     &	silo_dv, silo_dw, silo_sdum, silo_vflx, silo_fw,
     &	silo_fl, silo_fg, silo_yc, silo_co2_prop, silo_ices,
     &	silo_vpres, silo_wvpres )
     	! Clear silo variables, close file, and increase the cycle for next time input
     	silo_ierr = dbfreeoptlist(optlistid)

     	
     	silo_ierr = dbclose(silo_dbfile)
     	return
	end
	
	! writes scalar variables to file for each time step
	subroutine write_silo_c(silo_anv, silo_anl, icall,npt,neq,nspeci, ifdual)		      
		use avsio
		use comai, only : altc, days, icnl, jdate, jtime, nei_in,
     &     ns_in, verno, wdd, neq_primary, ivf, ifdm_elem
		use combi, only : corz, izonef, nelm
		use comchem
		use comdi, only : nsurf, izone_surf, izone_surf_nodes, icns,
     &     an, anv
		
		use compart, only : ptrak, pout
		use comrxni
		use comdti
		use silo
		
		implicit none
		
		include "silo.inc"
		
		integer add_dual, maxcon, iz, idz, iendz, il, open_file
		integer neq,nspeci,lu,ifdual,icall,length,i1,i2
		integer icord1, icord2, icord3, iaq, ivap, isolid
		integer npt(*), ip1, ip2
		integer, allocatable ::  nelm2(:)
		parameter (maxcon = 100)
		real*8, allocatable :: an_dum(:,:)
		real*8, allocatable :: anv_dum(:,:)
		real*8, allocatable :: antmp(:,:)
		character*60, allocatable :: title(:)
		character*14 tailstring
		character*8 dual_char
		character*3 dls, snum
		character*5 char_type
		character*60 fstring
		character*35 tmpname
		character*30 cordname(3)
		character*150 :: string = '', tecstring = '', sharestring = ''
		real*8 write_array(maxcon)
		integer i, ic, im, in, iv, ix, istep, j, k, n
		integer t1(maxcon),itotal2,write_total, iocord_temp
		integer irxn_title
		real*8 complex_conc
		real*8 minc, maxc
		character*10 silo_str_i		
		parameter (minc = 1.0d-90, maxc = 1.0d+20)
		real*8 silo_anv(n0, nspeci), silo_anl(n0, nspeci), silo_anr(n0, nspeci)
		

		
		
     	call create_silo('con')		
		
		
		! Set timestep for multiple function calls of write_silo_s()
      	silo_time = days
		silo_err = dbmkoptlist(100, optlistid)		
		silo_err = dbaddiopt(optlistid, DBOPT_CYCLE, silo_cycl)		
		silo_err = dbaddropt(optlistid, DBOPT_DTIME, silo_time)		
				
		! Rewrite mesh and connectivity
		silo_err = dbputzl2(silo_dbfile, "zonelist", 8, nei_in, ndims, silo_cv,
     &	nei_in*ns_in, 1, 0, 0, shapetype, shapesize, shapecounts, NSHAPETYPES, optlistid, silo_ierr)
     
     	silo_err = dbaddcopt(optlistid, DBOPT_UNITS, "m", 1)
	 	silo_err = dbputum(silo_dbfile, "pointmesh", 9, ndims, corz(:,1), corz(:,2), corz(:,3), "X", 1, "Y", 1, "Z", 1,
     &	DB_DOUBLE, neq_primary, nei_in, "zonelist", 8, DB_F77NULLSTRING, 0, optlistid, silo_ierr)		
     
     	
     	! Flip through all con variables, if rxn then make set the title, if not, then make our own general title
     	do i = 1, nspeci
			if(rxn_flag .eq. 0) write(silo_str_i, '(I5)') i			
			
			! Print vapor and liquid (eventually solid - zorb will be added) both for -2 or 2 (Henry) 
			select case(icns(i))
			case(2, -2)
				! Appending "_vapor" or "_liquid" for multiple output varaibles
				silo_conname = trim(cpntnam(i)) // "_Vapor"
				if(rxn_flag .eq. 0) silo_conname = "Vapor_Species" // trim(adjustl(silo_str_i))
				print *, " CON INPUT NAME: ", silo_conname
				call silo_put_unstruc2(silo_conname, "Moles/kg Vapor", silo_anv)				
				
				silo_conname = trim(cpntnam(i)) // "_Liquid"
				if(rxn_flag .eq. 0) silo_conname = "Aqueous_Species" // trim(adjustl(silo_str_i))
				print *, " CON INPUT NAME: ", silo_conname
				call silo_put_unstruc2(silo_conname, "Moles/kg H20", silo_anl)
			case(-1)
				silo_conname = trim(cpntnam(i)) // "_Vapor"
				if(rxn_flag .eq. 0) silo_conname = "Vapor_Species" // trim(adjustl(silo_str_i))
				print *, " CON INPUT NAME: ", silo_conname
				call silo_put_unstruc2(silo_conname, "Moles/kg Vapor", silo_anv)
			case(1)
				silo_conname = trim(cpntnam(i)) // "_Liquid"
				if(rxn_flag .eq. 0) silo_conname = "Aqueous_Species" // trim(adjustl(silo_str_i))
				print *, " CON INPUT NAME: ", silo_conname
				call silo_put_unstruc2(silo_conname, "Moles/kg H20", silo_anl)
			case(0)
				silo_conname = trim(cpntnam(i)) // "_Solid"
				if(rxn_flag .eq. 0) silo_conname = "Solid_Species" // trim(adjustl(silo_str_i))
				print *, " CON INPUT NAME: ", silo_conname
				print *, "NO ROCK CON OUTPUT -coming soon-"
			end select			
     	end do
     	


     	
     	! Clear silo variables, close file, and increase the cycle for next time input     		
     	silo_ierr = dbfreeoptlist(optlistid)

     	silo_cycl = icall
     	silo_ierr = dbclose(silo_dbfile)     			
		return      
	end