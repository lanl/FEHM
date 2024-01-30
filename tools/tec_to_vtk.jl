#===============================================================================
tec_to_vtk.jl

This script converts tecplot formatted fehm cont files into the vtk format for
reading into paraview. Cont files must be output in the tecplot format only.
To run, call this script in the folder where the fehmn.files is stored.
The script uses fehmn.files to point to the grid and other outfiles needed for
node locations.

Command to run: julia path/to/directory/tec_to_vtk.jl        
Dylan R. Harp (dharp@lanl.gov)
================================================================================

USAGE:

    --sca_diff			Sca file to subtract to create delta variables
    --sca_gdkm_diff		GDKM sca file to subtract to create delta variables
    --var_diff			Variables to create delta variables for
    --output_dir, -o		Folder to create and place vtu2 files
    --mesh_file, -m		Mesh file location
    --run_name, -r		Run name: prefix for visualization file names

================================================================================#

using WriteVTK
using Glob
using DataFrames
using ArgParse

function get_tec_time(f)
	# Tecplot to dataframe
    fh = open(f)
    # Read in header
    i = 1
    header = true
	tm = ""
    while header
        l = readline(fh)
        if startswith(l,"TITLE")
			title = strip(split(l,"=")[2])
        elseif startswith(l,"VARIABLES")
			variables = split(split(l,"=")[2],"\" \"")
			variables = [replace(v,"\""=>"") for v in variables]
			variables = [strip(replace(v,r"\(.+\)"=>"")) for v in variables]
        elseif startswith(l,"ZONE")
            tm = split(split(l,"\"")[2])[3]
            if occursin("days\"",tm)
				tm = replace(tm,"days\""=>"")
            elseif occursin("days",tm)
				tm = replace(tm,"days"=>"")
			end
        elseif startswith(l,"FILETYPE")
			continue
        elseif startswith(l,"SOLUTION")
			continue
        elseif occursin(r"^[a-z,\s]+$",split(l)[1])
            println(string("ERROR: Unrecognized header: ",split(l)[1],"; currently looks for TITLE, VARIABLES, and ZONE; might need to add it"))
        else
			header = false
		end
        i = i + 1
	end

	return tm
end

function tec_to_df(f,variables=0)
	# Tecplot to dataframe
    fh = open(f)
    # Read in header
    i = 1
    header = true
	tm = ""
    while header
        l = readline(fh)
        if startswith(l,"TITLE")
			title = strip(split(l,"=")[2])
        elseif startswith(l,"VARIABLES")
			variables = split(split(l,"=")[2],"\" \"")
			variables = [replace(v,"\""=>"") for v in variables]
			variables = [strip(replace(v,r"\(.+\)"=>"")) for v in variables]
        elseif startswith(l,"ZONE")
            tm = split(split(l,"\"")[2])
			if length(tm) <= 2
				i = i + 1
				continue
			else # Probably material file
            	tm = tm[3]
			end
            if occursin("days\"",tm)
				tm = replace(tm,"days\""=>"")
			end
        elseif startswith(l,"FILETYPE")
        	i = i + 1
			continue
        elseif startswith(l,"SOLUTION")
        	i = i + 1
			continue
		elseif isreal(parse(Float64,split(l)[1]))
			header = false
        else
            println(string("ERROR: Unrecognized header: ",split(l)[1],"; currently looks for TITLE, VARIABLES, and ZONE; might need to add it"))
		end
        i = i + 1
	end

	seek(fh,0)
	for i in 1:i-2
		readline(fh)
	end
	lns = readlines(fh)

	if occursin("gdkm",f)
		n_lines = n_gdkm
	else
		n_lines = numnodes
	end

	ncols = length(split(strip(lns[i])))
	# Remove x, y, z variable names they don't exist in tec file 
	# If variables length is 3 more than ncols, assume the first three variables are x,y,z coords and remove
	if length(variables) == ncols + 3
		variables = variables[4:end]
	end
	dtemp = Array{Float64,2}(undef,n_lines,length(variables))
	for i in 1:n_lines
		vs = split(strip(lns[i]))
		for j in 1:length(vs)
			dtemp[i,j] = parse(Float64,vs[j])
		end
	end

	d = Array{Float64,2}(undef,numnodes,length(variables))
	d[:,:] .= 0.
	if occursin("gdkm",f)
		d[gdkm_map[:,1],:] .= dtemp[:,:]
	else
		d[:,:] .= dtemp[:,:]
	end

	#d = readdlm(fh)
	close(fh)
	df = DataFrame(d,:auto)
	var_sym = Array{Symbol,length(variables)}
	variables = [replace(v,"/"=>"_") for v in variables]
	variables = [replace(v,"-"=>"_") for v in variables]
	variables = [replace(v," "=>"_") for v in variables]
	var_sym = [Symbol(v) for v in variables]
	rename!(df,var_sym)

	if (in(:X_coordinate,names(df)))
		delete!(df,:X_coordinate)
	end
	if (in(:Y_coordinate,names(df)))
		delete!(df,:Y_coordinate)
	end
	if (in(:Z_coordinate,names(df)))
		delete!(df,:Z_coordinate)
	end

	if (in(:X,names(df)))
		delete!(df,:X)
	end
	if (in(:Y,names(df)))
		delete!(df,:Y)
	end
	if (in(:Z,names(df)))
		delete!(df,:Z)
	end
	if (in(:Zone,names(df)))
		delete!(df,:Zone)
	end

	return df
end 

function write_vtk(f,coords, elems; diff_df=DataFrame(), variables=0, fi=-1)
	println(f)
	fnlst = split(f,'.')
    root = fnlst[1]
    if occursin("sca_node",fnlst[2])
		lst = split(fnlst[2],"_")
		if (fi==-1)
			fi = lst[1]
		end
        root = string(root,"_sca_node")
    elseif occursin(fnlst[2],"sca_gdkm_node")
		lst = split(fnlst[2],"_")
		if (fi==-1)
			fi = lst[1]
		end
        root = string(root,"_sca_gdkm_node")
    elseif occursin("con_node",fnlst[2])
		lst = split(fnlst[2],"_")
		if (fi==-1)
			fi = lst[1]
		end
        root = string(root,"_con_node")
    elseif occursin("con_gdkm_node",fnlst[2])
		lst = split(fnlst[2],"_")
		if (fi==-1)
			fi = lst[1]
		end
        root = string(root,"_con_gdkm_node")
    elseif occursin("mat_node",fnlst[2])
        root = string(root,"_mat_node")
	end

	df = tec_to_df(f,variables)

	if !isdir(parsed_args["output_dir"])
		mkdir(parsed_args["output_dir"])
	end
	if fi > -1
		fname = string(parsed_args["output_dir"],"/",root,lpad(fi,4,"0"))
	else
		fname = string(parsed_args["output_dir"],"/",root)
	end
	vtkfile = vtk_grid(fname, coords, elems)
	for v in names(df)
		vtk_point_data(vtkfile, df[!,v], string(v))
		if in(v,names(diff_df))
			if length(var_diff)>0
				if in(string(v),var_diff)
					vtk_point_data(vtkfile, df[!,v]-diff_df[!,v], string("Delta ",v))
				end
			else
				vtk_point_data(vtkfile, df[!,v]-diff_df[!,v], string("Delta ",v))
			end
		end
	end

	#outfiles = vtk_save(vtkfile)
	
	return vtkfile

end

s = ArgParseSettings(description = "Convert FEHM gdkm tecplot output to vtu2 files")

@add_arg_table! s begin
    "--sca_diff"
		help = "Sca file to subtract to create delta variables"
		arg_type = String
		default = "none"
    "--sca_gdkm_diff"
		help = "GDKM sca file to subtract to create delta variables"
		arg_type = String
		default = "none"
    "--var_diff"
		help = "Variables to create delta variables for"
		arg_type = String
		default = "none"
    "--output_dir", "-o"         # another option, with short form
		help = "Folder to create and place vtu2 files"
		arg_type = String
		default = "vtu"
    "--mesh_file", "-m"         # another option, with short form
		help = "Mesh file location"
		arg_type = String
		default = "none"
    "--run_name", "-r"         # another option, with short form
		help = "Run name: prefix for visualization file names"
		arg_type = String
		default = "*"
end

parsed_args = parse_args(s) # the result is a Dict{String,Any}
println("Parsed args:")
for (key,val) in parsed_args
    println("  $key  =>  $(repr(val))")
end

if parsed_args["var_diff"] != "none"
	var_diff = split(parsed_args["var_diff"],",")
else
	var_diff = Set()
end

if parsed_args["mesh_file"] == "none"
	fh = open("fehmn.files")
	mesh = ""
	for l in eachline(fh)
		if (startswith(l,"grid"))
			print(l)
			global mesh = strip(split(l,":")[2])
			break
		end
	end
	close(fh)
else
	mesh = parsed_args["mesh_file"]
end

fh = open(mesh)
lns = readlines(fh)
close(fh)
numnodes = parse(Int,strip(lns[2]))

# Points
coords = zeros(Float64,(3,numnodes))
for i in 3:numnodes+2
    vs = split(strip(lns[i]))
    coords[:,i-2] = [parse(Float64,vs[2]), parse(Float64,vs[3]), parse(Float64,vs[4])]
end

# Elements
vs = split(strip(lns[numnodes+5]))

numconns = parse(Int,vs[1])
numelems = parse(Int,vs[2])
elems = Array{MeshCell}(undef,numelems)
# Determine cell type
if ( numconns == 3 )
	celltype = VTKCellTypes.VTK_TRIANGLE
elseif ( numconns == 4 )
	if length(unique(coords[1,:]))>1 && length(unique(coords[2,:]))>1 && length(unique(coords[3,:]))>1
		celltype = VTKCellTypes.VTK_TETRA
	else
		celltype = VTKCellTypes.VTK_QUAD
	end
elseif ( numconns == 8 )
	celltype = VTKCellTypes.VTK_HEXAHEDRON
end
j = 1
for i in numnodes+6:numnodes+5+numelems 
    local vs = split(strip(lns[i]))
	if ( numconns == 3 )
    	elems[j] = MeshCell(celltype, [parse(Int,vs[2]), parse(Int,vs[3]), parse(Int,vs[4])])
	elseif ( numconns == 4 )
		# Check for triangle elements for DFNs
		if parse(Int,vs[5]) == 0
    			elems[j] = MeshCell(VTKCellTypes.VTK_TRIANGLE, [parse(Int,vs[2]), parse(Int,vs[3]), parse(Int,vs[4])])
		else
    			elems[j] = MeshCell(celltype, [parse(Int,vs[2]), parse(Int,vs[3]), parse(Int,vs[4]), parse(Int,vs[5])])
		end
	elseif ( numconns == 8 )
    	elems[j] = MeshCell(celltype, [parse(Int,vs[2]), parse(Int,vs[3]), parse(Int,vs[4]), parse(Int,vs[5]), parse(Int,vs[6]), parse(Int,vs[7]), parse(Int,vs[8]), parse(Int,vs[9])])
	end
	global j = j + 1
end

# Scalar files
fs = Glob.glob(string(parsed_args["run_name"],".","*_sca_node.dat"))
if (length(fs) > 0)
	pvd_sca = paraview_collection(string(parsed_args["output_dir"],"/sca"))
	# Collect variables from first file 
	sca_o = tec_to_df(fs[1])
	#print(names(sca_o))
	variables = [string(nm) for nm in names(sca_o)]
	if occursin("none",parsed_args["sca_diff"])
		sca_zero = DataFrame()
	else
		println(parsed_args["sca_diff"])
		sca_zero = tec_to_df(parsed_args["sca_diff"],variables)
	end
	#fi = Array(Int,len(fs))
	for (fi,f) in enumerate(fs)
	#	fi[fi] = fi
		vfile = write_vtk(f,coords,elems,diff_df=sca_zero,variables=variables,fi=fi);
		tm = get_tec_time(f)
		collection_add_timestep(pvd_sca, vfile, parse(Float64,tm))
	end
	#pmap(write_vtk,fs,fi)
	vtk_save(pvd_sca)
end

# Concentration files
fs = Glob.glob(string(parsed_args["run_name"],".","*_con_node.dat"))
if (length(fs) > 0)
	pvd_con = paraview_collection(string(parsed_args["output_dir"],"/con"))
	# Collect variables from first file 
	con_zero = tec_to_df(fs[1])
	variables = [string(nm) for nm in names(con_zero)]
	con_zero = DataFrame()
	for (fi,f) in enumerate(fs)
		vfile = write_vtk(f,coords,elems,diff_df=con_zero,variables=variables,fi=fi)
		tm = get_tec_time(f)
		collection_add_timestep(pvd_con, vfile, parse(Float64,tm))
	end
	vtk_save(pvd_con)
end

fs = Glob.glob(string(parsed_args["run_name"],".","mat_node.dat"))
if (length(fs) > 0)
	for (fi,f) in enumerate(fs)
		vfile = write_vtk(f,coords,elems,fi=fi)
		vtk_save(vfile)
	end
end

# Write gdkm.vtu
if isfile("gdkm.dat")
	fh = open("gdkm.dat")
	readline(fh)
	n_gdkm = parse(Int,split(strip(readline(fh)))[2])
	d = Array{Float64,(n_gdkm,2)}
	for i in 1:n_gdkm
		local vs = split(strip(readline(fh)))
		d[i,:] = [parse(Float64,vs[2]),parse(Float64,vs[3])]
	end
	readline(fh)
	gdkm_map = Array{Float64,(n_gdkm,2)}
	for i in 1:n_gdkm
		local vs = split(strip(readline(fh)))
		gdkm_map[i,:] = [parse(Int,vs[1]),numnodes+parse(Int,vs[4])]
	end
	close(fh)
	gdkm_dict = Dict(zip(gdkm_map[:,2],gdkm_map[:,1]))
	gdkm_dat = ones(Float64,(numnodes,2))*-9999
	gdkm_dat[gdkm_map[:,1],:] = d[:,:]

	if !isdir(parsed_args["output_dir"])
		mkdir(parsed_args["output_dir"])
	end
	vtkfile = vtk_grid(string(parsed_args["output_dir"],"/gdkm"), coords, elems)
	vtk_point_data(vtkfile, gdkm_dat[:,1], "VFrac Primary")
	vtk_point_data(vtkfile, gdkm_dat[:,2], "GDKM X")
	outfiles = vtk_save(vtkfile)
end


# Write perm.vtu
if isfile("perm.dat")
	perm_dat = ones(Float64,(numnodes,3))*-9999
	fh = open("perm.dat")
	readline(fh)
	for l in eachline(fh)
		local vs = split(strip(l))
		if isempty(strip(l))
			break
		end
		if parse(Int,vs[1]) > 1
			perm_dat[gdkm_dict[parse(Int,vs[1])],:] = [parse(Float64,vs[4]),parse(Float64,vs[5]),parse(Float64,vs[6])]
		end
	end
	close(fh)

	if !isdir(parsed_args["output_dir"])
		mkdir(parsed_args["output_dir"])
	end
	vtkfile = vtk_grid(string(parsed_args["output_dir"],"/perm"), coords, elems)
	vtk_point_data(vtkfile, perm_dat[:,1], "Fracture Perm X")
	vtk_point_data(vtkfile, perm_dat[:,1], "Fracture Perm Y")
	vtk_point_data(vtkfile, perm_dat[:,1], "Fracture Perm Z")
	outfiles = vtk_save(vtkfile)
end

# Scalar gdkm files
fs = Glob.glob(string(parsed_args["run_name"],".","*_sca_gdkm_node.dat"))
if (length(fs) > 0)
	pvd_sca = paraview_collection(string(parsed_args["output_dir"],"/sca_gdkm"))
	# Collect variables from first file 
	sca_o = tec_to_df(fs[1])
	#print(names(sca_o))
	variables = [string(nm) for nm in names(sca_o)]
	if occursin("none",parsed_args["sca_gdkm_diff"])
		sca_zero = DataFrame()
	else
		println(parsed_args["sca_gdkm_diff"])
		sca_zero = tec_to_df(parsed_args["sca_gdkm_diff"])
	end
	#fi = Array(Int,len(fs))
	for (fi,f) in enumerate(fs)
	#	fi[fi] = fi
		vfile = write_vtk(f,coords,elems,diff_df=sca_zero,variables=variables,fi=fi);
		tm = get_tec_time(f)
		collection_add_timestep(pvd_sca, vfile, parse(Float64,tm))
	end
	#pmap(write_vtk,fs,fi)
	vtk_save(pvd_sca)
end

# Concentration files
fs = Glob.glob(string(parsed_args["run_name"],".","*_con_gdkm_node.dat"))
if (length(fs) > 0)
	pvd_con = paraview_collection(string(parsed_args["output_dir"],"/con_gdkm"))
	# Collect variables from first file 
	con_zero = tec_to_df(fs[1])
	variables = [string(nm) for nm in names(con_zero)]
	con_zero = DataFrame()
	for (fi,f) in enumerate(fs)
		vfile = write_vtk(f,coords,elems,diff_df=con_zero,variables=variables,fi=fi)
		tm = get_tec_time(f)
		collection_add_timestep(pvd_con, vfile, parse(Float64,tm))
	end
	vtk_save(pvd_con)
end


