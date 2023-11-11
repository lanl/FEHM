import numpy
import numpy as np
from numpy.lib.recfunctions import append_fields
import os,sys
from os.path import join
from glob import glob
import re
from matk import pest_io
from shutil import rmtree,copyfile
from subprocess import call,STDOUT
import shlex
from scipy import interpolate

def compile_results(dirs, zones=[26]):
    # Initialize header and format
    header_c = ['Time[h] Sample ID [moles 133Xe/kg-air]:']
    header_m = ['Time[h] Sample ID [moles 133Xe/day]:']
    fmt = ['%d']

    for j,dr in enumerate(dirs):
        # Figure out initial time
        with open(join(dr,'tracer','run.out'), 'r') as fh:
           lns = fh.readlines()
           for i,s in enumerate(lns):
                if "Timing Information" in s: break
           ti = float(lns[i+2].split()[1]) - float(lns[i+2].split()[2])
        ######## CONCENTRATIONS ################
        # Get top nodes
        fh = open(join(dr,'tracer','run.his'),'r')
        for i in range(9): fh.readline()
        nds = numpy.genfromtxt(fh)
        nds = nds[numpy.where(nds[:,2]==numpy.max(nds[:,2]))[0],0]
        nds = [str(int(v)) for v in nds]
        fh.close()
        # Create data frame so that columns associated with top nodes can be collected
        fh = open(join(dr,'tracer','run_xenon.dat'),'r')
        fh.readline()
        ln = fh.readline()
        varnames = ln.strip().split('=')[1].split('" "Node')
        varnames = [v.strip() for v in varnames]
        varnames = [v.replace('"','') for v in varnames]
        ln = fh.readline()
        ln = fh.readline()
        c = numpy.genfromtxt(fh,names=varnames)
        fh.close()
        # Find max concentrations along top
        cs = numpy.max(numpy.array(c[nds].tolist()),axis=1)
        # Set spurious negative values to zero
        numpy.place(cs,cs<0,0)
        # Collect times
        ts = c['Time_days']-ti
        if j == 0:
            # Create hourly times
            tinterp = numpy.arange(1./24, ts.max(), 1./24)
            # Initialize data array
            cout = tinterp*24
        # Create interpolator function
        f = interpolate.interp1d(ts, cs)
        # Interpolate concentations to hourly times
        cinterp = f(tinterp)
        # Add results to data array
        cout = numpy.column_stack([cout,cinterp])
        # Append header
        header_c.append(str(dr))
        fmt.append('%g')
        ######## MASS FLUX ################
        d = cflx_recarray(join(dr,'tracer','runxenon.cflx'))
        for i,z in enumerate(zones):
            if i == 0:
                ms = d['Sink'+str(z)]
            else:
                ms += d['Sink'+str(z)]
        #m = numpy.genfromtxt(join(dr,'tracer','runxenon.cflx'),skip_header=2)
        # Collect matrix and fracture (gdkm) NetOut values
        # Set spurious negative values to zero
        #m6 = m[:,4]
        #numpy.place(m6,m6<0,0)
        #m7 = m[:,9]
        #numpy.place(m7,m7<0,0)
        #ms = m6+m7
        # Collect times
        ts = d['Time_days']-ti
        if j == 0:
            # Initialize data array
            mout = tinterp*24
        # Create interpolator function
        f = interpolate.interp1d(ts, ms)
        # Interpolate concentations to hourly times
        minterp = f(tinterp)
        # Add results to data array
        mout = numpy.column_stack([mout,minterp])
        # Append header
        header_m.append(str(dr))
    numpy.savetxt('gas_concs.dat',cout,header=' '.join(header_c),fmt=fmt)
    numpy.savetxt('gas_massflux.dat',mout,header=' '.join(header_m),fmt=fmt)

def cflx_recarray(cflx_file):
    with open(cflx_file,'r') as fh:
        fh.readline()
        zns = [int(v) for v in fh.readline().split(":")[1].split()]
        names = ['Time_days']
        for zn in zns:
            name_str = 'Source%d Sink%d NetIn%d NetOut%d Boundary%d'%(zn,zn,zn,zn,zn)
            names += name_str.split()
        return numpy.genfromtxt(fh,names=names)


def tecplot_dict(tecplot_file,variables=None):
    with open(tecplot_file,'r') as fh:
        header = True
        header_count = 0
        while header:
            ln = fh.readline()
            if not ln.split()[0].isalpha():
                header = False
                break
            hs = ln.split(',')
            for h in hs:
                l = h.split('=')
                if 'VARIABLES' in l[0]:
                    if variables is None:
                        variables = shlex.split(l[1])
                elif 'TITLE' in l[0]:
                    pass
                elif 'Zone T' == l[0].strip():
                    t = float(l[1].split()[2])
                elif 'N' == l[0].strip(): 
                    nn = int(l[1])
                elif 'E' == l[0].strip(): 
                    ne = int(l[1])
            header_count += 1
        fh.seek(0)
        d = numpy.genfromtxt(fh,skip_header=header_count,max_rows=nn,names=variables)
    return d

def finfile_dict(finfile):
    with open(finfile,'r') as fh:
        fh.readline()
        fh.readline()
        t = float(fh.readline().strip())
        nn = int(fh.readline().split()[0])
        d = {}
        d['node'] = numpy.arange(nn) + 1
        while True:
            k = fh.readline().strip()
            if k == 'no fluxes' or k == '':
                break
            d[k] = numpy.fromfile(fh,count=nn,sep=' ')

    return d

def zonefile_dict(zone_file):
    with open(zone_file,'r') as fh:
        if fh.readline().strip() not in ['zone','zonn']:
            print("Error: %s does not appear to be an FEHM zone file"%zone_file)
            return
        d = {}
        while True:
            l = fh.readline()
            if not l.strip() in ['\n','\r\n','stop']:
                k = int(l.strip().split()[0])
                if 'nnum' in fh.readline().strip().split()[0]:
                    n = int(fh.readline().split()[0])
                    d[k] = numpy.fromfile(fh,count=n,dtype=int,sep=' ')
                else:
                    print("Only zone type nnum supported currently!")
                    return
            else: break
    return d

def grid_recarray(mesh_file):
    # Read grid
    fh = open(mesh_file,'r')
    fh.readline()
    numnodes = int(fh.readline().strip())
    dtype = {'names':('node','x','y','z'),
             'formats':('i8','f8','f8','f8')
            }
    grid = numpy.zeros((numnodes,),dtype=dtype)
    for i in range(numnodes): 
        vs = fh.readline().strip().split()
        grid[i] = (int(vs[0]), float(vs[1]), float(vs[2]), float(vs[3]))
    fh.close()
    return grid

def add_damage(grid,avs_file,stor_file,maxfracap,damage_field='dam_map'):
    # Read damage
    fh = open(avs_file,'r')
    ln = fh.readline().strip().split()
    nnodes = int(ln[0])
    if not nnodes == len(grid):
        print("AVS file does not appear to match grid file: node counts are different")
        sys.exit(0)
    nelems = int(ln[1])
    numatts = int(ln[2])
    for a in range(nnodes+nelems+1): fh.readline()
    names = ['id']
    for a in range(numatts): names.append(fh.readline().split(',')[0])
    dam_rc = numpy.genfromtxt(fh,names=names,max_rows=nnodes)
    grid = append_fields(grid,'damage',dam_rc[damage_field])
    fh.close()

    # Read arealist
    fh = open(stor_file,'r')
    fh.readline()
    fh.readline()
    nnodes = int(fh.readline().split()[1])
    if not nnodes == len(grid):
        print("FEHM stor file does not appear to match grid file: node counts are different")
        sys.exit(0)
    grid = append_fields(grid,'area',numpy.fromfile(fh,sep=" ", dtype=float, count=nnodes))

    # Calculate scaled damage (divide by max damage to ensure that max damage is 1
    grid = append_fields(grid,'scaled_damage',maxfracap*grid['damage']/numpy.max(grid['damage']))

    # Calculate fracture perms based on cubic law
    grid = append_fields(grid,'k_f',grid['scaled_damage']**2./12.)

    return grid

def gdkm_nocavchim_macro(grid,matzonefile,vfrac_primary=0.99,gdkm_x=1.e-4,macro_dir='.'):
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir) 

    # gdkm nodes, expects that second zone in file is gdkm
    # third zone is source nodes
    # fourth zone is fault top node
    fh = open(matzonefile,'r')
    for i in range(3): fh.readline()
    n = int(fh.readline().strip())
    t = numpy.fromfile(fh,count=n,sep=' ')
    for i in range(2): fh.readline()
    n = int(fh.readline().strip())
    gdkm_nodes = numpy.fromfile(fh,count=n,sep=' ').astype('int')
    #for i in range(2): fh.readline()
    #n = int(fh.readline().strip())
    #source_nodes = numpy.fromfile(fh,count=n,sep=' ').astype('int')
    #for i in range(2): fh.readline()
    #n = int(fh.readline().strip())
    #faulttop_nodes = numpy.fromfile(fh,count=n,sep=' ').astype('int')

    # Create subset of grid with only gdkm and source nodes
    gdkm_inds = numpy.where(numpy.logical_or.reduce([grid['node'] == x for x in gdkm_nodes])==True)[0]
    #source_inds = numpy.where(numpy.logical_or.reduce([grid['node'] == x for x in source_nodes])==True)[0]
    #faulttop_inds = numpy.where(numpy.logical_or.reduce([grid['node'] == x for x in faulttop_nodes])==True)[0]
    sort_inds = numpy.sort(numpy.concatenate([gdkm_inds]))    
    gdkm = grid[sort_inds]
    gdkm = numpy.sort(gdkm, order='node') 
    #top = grid[faulttop_inds]

    # Add secondary nodes 
    gdkm = append_fields(gdkm,'node2',len(grid) + numpy.arange(1,len(gdkm)+1))

    fh = open(os.path.join(macro_dir,'gdkmzone.zone'),'w')
    fh.write('zone\n102\nnnum\n'+str(len(gdkm))+'\n')
    numpy.savetxt(fh,gdkm['node2'],fmt=' %d')
    fh.write('\nstop\n')
    fh.close()

    #fh = open(os.path.join(macro_dir,'topnodes.zonn'),'w')
    #fh.write('zonn\n106\nnnum\n'+str(len(top))+'\n')
    #numpy.savetxt(fh,top['node'],fmt=' %d')
    #fh.write('\nstop\n')
    #fh.close()

    fh = open(os.path.join(macro_dir,'gdkm.dat'),'w')
    outtable_gdkm = numpy.ones([len(gdkm),3])
    # Shouldn't this be different?
    # VFRAC_PRIMARY could be based on aperture...
    outtable_gdkm[:,1] = vfrac_primary
    # GDKM_X should be smaller for greater damamge...
    outtable_gdkm[:,2] = gdkm_x
    outtable_gdkm2 = numpy.ones([len(gdkm),4])
    outtable_gdkm2[:,0] = gdkm['node']
    outtable_gdkm2[:,1] = gdkm['node']
    outtable_gdkm2[:,2] = 1
    outtable_gdkm2[:,3] = numpy.arange(1,len(gdkm)+1)
    fh.write('gdkm\n1 '+str(len(gdkm))+'\n')
    numpy.savetxt(fh,outtable_gdkm,fmt='%d %16.15g %16.15g')
    fh.write('\n')
    numpy.savetxt(fh,outtable_gdkm2,fmt='%d %d %d %d')
    fh.write('\nstop\n')
    fh.close()

    return gdkm

def gdkm_macro(grid,mat_cavchim,k_gdkm_cutoff,mat_outside,macro_dir='.',
              gdkm_top_nodes=True,mat_cavchim_zoneno=11):

    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir) 

    # Cavity and chimney nodes, expects that zone 11 in file is cavity/chimney 
    zc = zonefile_dict(mat_cavchim)
    #ccc = zccc[2]

    # Get dictionary of outside zones
    zo = zonefile_dict(mat_outside)

    # Create gdkm recarray by removing cavchim nodes
    # Instead try setting all cavchim node perms to max k_f so that they are not removed
    # Modify scaled damage
    gdkm = grid
    ccc_nodes = numpy.where(numpy.logical_or.reduce([grid['node'] == x for x in zc[mat_cavchim_zoneno]])==True)[0]
    for nd in ccc_nodes: 
        gdkm[nd]['k_f'] = numpy.max(gdkm['k_f'])
        gdkm[nd]['scaled_damage'] = numpy.max(gdkm['scaled_damage'])
    # Remove top nodes from gdkm (assume zone 6 is top (lagrit convention))
    if not gdkm_top_nodes:
        not_top_nodes = numpy.where(numpy.logical_or.reduce([grid['node'] == x for x in zo[6]])!=True)[0]
        gdkm = gdkm[not_top_nodes]
    # Then remove nodes with perm less than k_gdkm_cutoff
    gdkm = gdkm[numpy.where(gdkm['k_f']>=k_gdkm_cutoff)[0]]
    
    # Add secondary nodes 
    gdkm = append_fields(gdkm,'node2',len(grid) + numpy.arange(1,len(gdkm)+1))

    #fh = open(os.path.join(macro_dir,'cavchim.zonn'),'w')
    #fh.write('zonn\n101\nnnum\n'+str(len(zc[2]))+'\n')
    #numpy.savetxt(fh,zc[2],fmt=' %d')
    #fh.write('\nstop\n')
    #fh.close()

    fh = open(os.path.join(macro_dir,'gdkmzone.zone'),'w')
    fh.write('zone\n102\nnnum\n'+str(len(gdkm))+'\n')
    numpy.savetxt(fh,gdkm['node2'],fmt=' %d')
    fh.write('\nstop\n')
    fh.close()

    #fh = open(os.path.join(macro_dir,'gdkmtopnodes.zonn'),'w')
    #gdkmtopnodes3 = gdkm['node2'][numpy.where(gdkm['y']==numpy.max(gdkm['y']))[0]] 
    ##print gdkm['y'][numpy.where(gdkm['y']==numpy.max(gdkm['y']))[0]] 
    #fh.write('zonn\n106\nnnum\n'+str(len(gdkmtopnodes3))+'\n')
    #numpy.savetxt(fh,gdkmtopnodes3,fmt=' %d')
    #fh.write('\nstop\n')
    #fh.close()

    fh = open(os.path.join(macro_dir,'gdkm.dat'),'w')
    outtable_gdkm = numpy.ones([len(gdkm),3])
    # VFRAC_PRIMARY; approximated assuming a square finite volume with length sqrt(A) sides
    outtable_gdkm[:,1] = 1 - gdkm['scaled_damage'] / numpy.sqrt(gdkm['area'])
    # GDKM_X should be smaller for greater damamge...
    #outtable_gdkm[:,2] = gdkm['scaled_damage'] / distdiv
    # GDKM_X now calculated to be half of cell width assuming square finite volume
    outtable_gdkm[:,2] = numpy.sqrt(gdkm['area'])/2.
    outtable_gdkm2 = numpy.ones([len(gdkm),4])
    outtable_gdkm2[:,0] = gdkm['node']
    outtable_gdkm2[:,1] = gdkm['node']
    outtable_gdkm2[:,2] = 1
    outtable_gdkm2[:,3] = numpy.arange(1,len(gdkm)+1)
    fh.write('gdkm\n1 '+str(len(gdkm))+'\n')
    numpy.savetxt(fh,outtable_gdkm,fmt='%d %16.15g %16.15g')
    fh.write('\n')
    numpy.savetxt(fh,outtable_gdkm2,fmt='%d %d %d %d')
    fh.write('\nstop\n')
    fh.close()

    ## Write cavchim secondary zone file
    #ccc_secs = numpy.arange(gdkm['node2'][-1]+1,gdkm['node2'][-1]+len(ccc_nodes) + 1)
    #fh = open(os.path.join(macro_dir,'cavchim_gdkm.zonn'),'w')
    #fh.write('zonn\n10\nnnum\n'+str(len(ccc_nodes))+'\n')
    #numpy.savetxt(fh,ccc_secs,fmt=' %d')
    #fh.write('\nstop\n')
    #fh.close()

    return gdkm

def node_macro(grid,macro_dir='.'):
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir) 

    fh = open(os.path.join(macro_dir,'nodeline.dat'),'w')
    # Vertical line of nodes (x=0)
    get = grid[numpy.where(grid['x']==0)[0][0:-1:5]]
    get.sort(order='y')
    get = get[0:-1]['node']
    # Horizontal line of nodes (y max)
    gett = grid[numpy.where(grid['y']==numpy.max(grid['y']))[0]]
    gett.sort(order='x')
    gett = gett['node']
    print_all = numpy.concatenate([gett,get])
    fh.write('node\n'+str(len(print_all))+'\n')
    numpy.savetxt(fh,print_all,fmt=' %d')
    fh.write('\nstop\n')
    fh.close()

def perm_macro(grid, gdkm, k_bg, macro_dir='.',zones=[],zone_perms=[]):
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir) 

    fh = open(os.path.join(macro_dir,'perm.dat'),'w')
    fh.write('perm\n  1 0 0 '+str(k_bg)+' '+str(k_bg)+' '+str(k_bg)+'\n')
    ps = numpy.column_stack([gdkm['node2'],gdkm['node2'],numpy.ones_like(gdkm['node2']),gdkm['k_f'],gdkm['k_f'],gdkm['k_f']])
    numpy.savetxt(fh,ps,fmt=' %d %d %d %16.15g %16.15g %16.15g')
    for zn,k in zip(zones,zone_perms):
        if zn < 0:
            fh.write(' {0} 0 0 {1} {1} {1}\n'.format( int(zn), k ))
        else:
            fh.write(' {0} {0} 1 {1} {1} {1}\n'.format( int(zn), k  ))
    fh.write('\nstop\n')
    fh.close()

def flow_macro(grid,init_dir,macro_dir='.',flow_depth_perc=0.):
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir) 

    fs = sorted(glob(os.path.join(init_dir,'*_sca_node.dat')))

    # Collect variables from first file 
    fh = open(fs[0],'r')
    fh.readline()
    variables = fh.readline().strip().split('=')[1].split('" "')
    variables = [re.sub('"','',v) for v in variables]
    variables = [re.sub('\(.+\)','',v).strip() for v in variables]
    fh.close()

    # Add vapor pressures to grid structure
    with open(fs[-1],'r') as fh:
        for l in fh:
            if l.startswith('ZONE'):
                time = l.split()[4]
                if 'days"' in time: time = re.sub('days"','',time)
                break
        grid = append_fields(grid,'P',numpy.genfromtxt(fh,names=variables)['Vapor_Pressure'])

    xmax = numpy.max(grid['x'])
    ymax = numpy.max(grid['y'])
    max_x = grid[numpy.where(grid['x']>=xmax)]
    max_x_low = max_x[numpy.where(max_x['y']<=ymax*(1-flow_depth_perc/100.))]

    with open(os.path.join(macro_dir,'flow.dat'),'w') as fh:
        fh.write('flow\n')
        for i in range(len(max_x_low)):
            fh.write(' {0:5} {1:5} 1 {2:16} -0.1 1.e3\n'.format(max_x_low['node'][i], max_x_low['node'][i], max_x_low['P'][i]))
        fh.write('\nstop\n')

def run_fehm(p,name,exe=None,tpl_dir='templates',macro_dir='macros',verbose=False):
    #  # Figure out FEHM executable
    #  if exe is None and 'FEHM_EXE' in os.environ: exe = os.environ['FEHM_EXE']
    #  else: exe = 'xfehm'
    # Figure out FEHM executable
    if exe is None and 'FEHM_EXE' in os.environ: exe = os.environ['FEHM_EXE']
    elif exe is not None: exe = exe  #user can provide exe path
    else: exe = 'xfehm'
    # Set verbosity in fehmn.files
    if verbose: p['verbose']='all'
    else: p['verbose']='none'

    # Create directory, delete first if it exists
    if os.path.isdir(name):
        rmtree(name)
    os.mkdir(name)
    # Create FEHM files
    pest_io.tpl_write(p,join(tpl_dir,'run_'+name+'.tpl'),join(name,'run.dat'))
    pest_io.tpl_write(p,join(tpl_dir,'fehmn_'+name+'.tpl'),join(name,'fehmn.files'))
    os.chdir(name)
    #  # Link macros
    #  for f in glob(join('..',macro_dir,'*')):
        #  lnk = f.split(os.sep)[-1]
        #  if not os.path.isfile(lnk):
            #  os.symlink(f,lnk)
    # Copy macros instead of soft linking them  (for posterity)
    for f in glob(join('..',macro_dir,'*')):
        dst = f.split(os.sep)[-1]
        copyfile(f,dst)
    # Run model
    if verbose:
        ierr = call([exe])
    else:
        FNULL = open(os.devnull, 'w')
        ierr = call([exe], stdout=FNULL, stderr=STDOUT)
        FNULL.close()
    os.chdir('..')
    return ierr

#  def write_boun(out_file, t, ys, step="ti_linear", varname='pa', zones=[7,106], comment=None):
    #  '''
    #  Create a .boun file with single or multiple variables for a single model
    #  number (e.g. model1 = just pressure, pressure & temp, etc).
#  
    #  Parameters
    #  ----------
    #  out_file : str
        #  Output file name.
    #  t : ndarray
        #  Time steps [d] for boundary condition changes.
    #  ys : ndarray or list(ndarrays)
        #  Variable(s) boundary condition changes. Use list for multiple.
    #  step : str
    #  Time sequence for changes to follow (see FEHM documentation). Choices are:
        #  "ti", "ti_linear", "cy", "cy_linear".
    #  varname : str or list(str)
        #  Variable(s) names. Use list for multiple.
    #  zones : int or list(int)
        #  Zone numbers to apply BC changes to. Use list for multiple.
    #  '''
    #  fout = open(out_file,'w')
    #  #fout.write('# Mean pressure: '+str(numpy.mean(y))+'\n')
    #  fout.write('boun\n')
    #  ys = numpy.atleast_2d(numpy.array(ys))
    #  if type(varname) == str: varname=[varname]
    #  model = 1
    #  v = 0 #variable count
    #  for y in ys:
        #  # write the timesteps for BC changes
        #  if v==0:
            #  fout.write('model '+str(model)+'\n'+step+'\n')
            #  model += 1
            #  fout.write(str(len(t))+'\n')
            #  cnt = 0
            #  for i in range(len(t)):
                #  #  fout.write(' %12g' % t[i])
                #  #  fout.write(' %12f' % t[i])
                #  fout.write('{: >14.6f}'.format(t[i]))
                #  cnt += 1
                #  if cnt == 5 or i == len(t)-1:
                    #  fout.write('\n')
                    #  cnt = 0
        #  # write the variable BC changes
        #  fout.write(varname[v]+'\n')
        #  cnt = 0
        #  for i in range(len(y)):
            #  fout.write(' %12g' % y[i])
            #  #  fout.write(' %12f' % y[i])
            #  #  fout.write('{: >14.2f}'.format(y[i]))
            #  cnt += 1
            #  if cnt == 5 or i == len(y)-1:
                #  fout.write('\n')
                #  cnt = 0
        #  v+=1
    #  # write the zones numbers
    #  fout.write('\n')
    #  zones = numpy.atleast_2d(numpy.array(zones))
    #  for m,zs in enumerate(zones):
        #  for z in zs:
            #  if z < 0: fout.write(' %d 0 0 %d\n' % (z,m+1))
            #  else: fout.write(' %d %d 1 %d\n' % (z,z,m+1))
    #  fout.write('\nstop\n')
    #  if comment:fout.write('# '+comment)
    #  fout.close()

#  #  #TESTING WRITE_BOUN()
#  #-- (1 model, single var)
#  out_file = 'test.boun'
#  t = np.linspace(0.,10.)
#  temps = [4.]*50; temps = np.array(temps)
#  ys = temps
#  step = "cy_linear"
#  varname = "t"
#  zones = -100
#  comment = None
#  #  #-- (1 model, multiple vars)
#  #  out_file = 'test.boun'
#  #  t = np.linspace(0.,10.)
#  #  temps = [4.]*50; temps = np.array(temps)
#  #  ps = [800.]*50; ps = np.array(ps)
#  #  ys = [temps,ps]
#  #  step = "cy_linear"
#  #  varname = ["t","pw"]
#  #  zones = -100
#  #  comment = None
#  #  #-- (2 models, multiple vars)
#  #  out_file = 'test.boun'
#  #  ts = np.linspace(0.,10.)
#  #  t = [ ts, [0.,1e20] ]
#  #  temps = [4.]*50; temps = np.array(temps)
#  #  ps = [800.]*50; ps = np.array(ps)
#  #  pas = [900.]*2; pas = np.array(pas)
#  #  ys = [ [temps,ps], pas ]
#  #  step = ["cy_linear","ti_linear"]
#  #  varname = [ ["t","pw"], ["pa"] ]
#  #  zones = [ [-100],[-200] ]
#  #  comment = None
#  #  write_boun(out_file, t, ys, step, varname,zones)
#  #  #-- (2 models, multiple vars, model 2 zone specified using JA JB JC)
#  #  out_file = 'test.boun'
#  ts = np.linspace(0.,10.)
#  t = [ ts, [0.,1e20] ]
#  temps = [4.]*50; temps = np.array(temps)
#  ps = [800.]*50; ps = np.array(ps)
#  pas = [900.]*2; pas = np.array(pas)
#  ys = [ [temps,ps], pas ]
#  step = ["cy_linear","ti_linear"]
#  varname = [ ["t","pw"], ["pa"] ]
#  zones = [ [-100],[2101, 2152] ]
#  #  zones = [ [2101, 2152], [-100] ]
#  comment = None
#  #  write_boun(out_file, t, ys, step, varname,zones)

def write_boun(out_file, t, ys, step="ti_linear", varname='pa', zones=[7,106], comment=None, highPrecision=None):
    '''
    Create a .boun file with single or multiple variables for a single model
    number (e.g. model1 = just pressure, pressure & temp, etc).

    Parameters
    ----------
    out_file : str
        Output file name.
    t : ndarray or list(ndarrays)
        Time steps [d] for boundary condition changes.
    ys : ndarray or list(ndarrays) or list(list(ndarrays))
        Variable(s) boundary condition changes. Use list for multiple.
    step : str or list(str)
        Time sequence for changes to follow (see FEHM documentation). Choices are:
        "ti", "ti_linear", "cy", "cy_linear". Number of step entries dictates
        the number of models that will bein the .boun file.
    varname : str or list(str) or list(list(str))
        Variable(s) names.
        Use str for single variable (e.g., "pw").
        Use list for multiple variables (e.g., ["t","pw"] is one model with
        temperature and water pressure changes.
        Use a list of lists for multiple models (e.g., [["t","pw"],["pa"]] has
        model1: temp and water pressure changes, and model2: air pressure
        changes).
    zones : int or list(int) or list(list(int))
        Zone numbers to apply BC changes to. Use list for multiple.
        Use list of lists for multiple models (zones for each model
        should be in a sub-list).
        - e.g., zones = [ [-100],[2101, 2152] ]. This will assign model 1 to
          zone number 100, and model 2 using the JA JB JC format in FEHM
          (default to 1-node spacing: i.e., 2101 2152 1).
        For now, assume you can only have one zone per model per call.
    '''

    err = "ERROR: Problem writing {}\n".format(out_file)
    fout = open(out_file,'w')
    # Check if there are multiple models desired in this .boun
    n_models = 1
    if type(step)==list:
        n_models = len(step)
    zones = numpy.atleast_2d(numpy.array(zones))
    #fout.write('# Mean pressure: '+str(numpy.mean(y))+'\n')
    fout.write('boun\n')
    # Make formatting of variables the same regardless of # of models
    if n_models==1:
        step = [step]
        varname = [varname]  #possible problem ??
        #  if len(varname) == 1:
            #  varname = [varname]  #possible problem ??
        ys = [ys]
        t = [t]
    # Check to make sure there are parameters provided for each model
    if len(t) < n_models: sys.exit(err+"Fewer time arrays than number of models.")
    #  elif len(np.shape(ys)) < n_models: sys.exit(err+"Fewer ys arrays than number of models.")
    elif len(step) < n_models: sys.exit(err+"Fewer step calls than number of models.")
    elif len(varname) < n_models: sys.exit(err+"Fewer varname calls than number of models.")
    #  elif np.shape(zones)[0] < n_models: sys.exit(err+"Fewer zone calls than number of models.")
    elif max(np.shape(zones)) < n_models: sys.exit(err+"Fewer zone calls than number of models.")
    # Loop for each model
    for model in range(1,n_models+1):
        #  ys = numpy.atleast_2d(numpy.array(ys))[model-1]
        #  ys = numpy.atleast_2d(numpy.array(ys[model-1]))
        yss = numpy.array(ys[model-1]) #ys for current model    #possible problem  ??
        # Make formatting of each model's ys the same
        #-- 1st index=# of models, 2nd index=max # of variable types for all models
        if len(np.shape(yss))==1:
            yss = [yss]
        ts = numpy.array(t[model-1])  #times for current model
        v = 0 #variable count

        for y in yss:
            # write the timesteps for BC changes
            if v==0:
                fout.write('model '+str(model)+'\n'+step[model-1]+'\n')
                #  model += 1
                fout.write(str(len(t[model-1]))+'\n')
                cnt = 0
                for i in range(len(ts)):
                    #  fout.write(' %12g' % t[i])
                    #  fout.write(' %12f' % t[i])
                    # Fix issue with decimal accuracy for big numbers
                    if ts[i] > 1e3:
                        #  fout.write('{: >14.6g}'.format(ts[i]))
                        # Use python default for v. large numbers (e.g., 1e20)
                        #  if ts[i] >=1e6: fout.write('{: >14.6g}'.format(ts[i]))
                        if ts[i] >=1e6: fout.write('{: >14.2e}'.format(ts[i]))
                        # Use very high prevission (8 decimals) if specified by user
                        elif highPrecision is not None:
                            print('User triggered highPrecision!')
                            fout.write('{: >14.8f}'.format(ts[i]))
                        # Otherwise, accuracy to 6 decimals (needed for FEHM interp)
                        else: fout.write('{: >14.6f}'.format(ts[i]))
                    else: fout.write('{: >14.6f}'.format(ts[i]))
                    cnt += 1
                    if cnt == 5 or i == len(t[model-1])-1:
                        fout.write('\n')
                        cnt = 0
            # write the variable BC changes
            #  fout.write(varname[model-1][v]+'\n')
            if type(varname[model-1]) is list: fout.write(varname[model-1][v]+'\n')
            else: fout.write(varname[model-1]+'\n')
            cnt = 0
            for i in range(len(y)):
                fout.write(' %12g' % y[i])
                #  fout.write(' %12f' % y[i])
                #  fout.write('{: >14.2f}'.format(y[i]))
                cnt += 1
                if cnt == 5 or i == len(y)-1:
                    fout.write('\n')
                    cnt = 0
            v+=1
    # Write the zones numbers for each model at the end
    #---- Perhaps for now, assume you can only apply 1 model per zone per call ...
    fout.write('\n')
    zones_store = zones
    try: zones[model-1]
    except IndexError: zones = zones_store[0]
    for model in range(1,n_models+1):
        z = zones[model-1]
        # type(z)=list means a mix of ZONEN and JA JB JC
        if type(z) is list:
            # If entry is a ZONEN specification
            if len(z) == 1:
                if z[0] < 0: fout.write(' %d 0 0 %d\n' % (z[0],model))
                else: fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
            # If entry is a JA JB JC specification
            else:
                fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
        # type(z)=numpy.ndarray  means a just ZONEN
        elif type(z) is numpy.ndarray:
            if z[0] < 0: fout.write(' %d 0 0 %d\n' % (z[0],model))
            # If entry is a JA JB JC specification
            else:
                fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
        #  
        #  #  zones = zones_store
        #  #  try: zones[model-1]
        #  #  except IndexError: zones = zones_store[0]
        #  for z in zones[model-1]:
            #  print('zones is')
            #  print(zones)
            #  print('z is')
            #  print(z)
            #  # If z is list, there may be a JA JB JC specified zone 
            #  if type(z) is list:
                #  print('len(z) = {}'.format(len(z)))
                #  print('type(z) is list')
                #  if len(z) == 1:
                    #  print('-- check1 -- ')
                    #  print('z0 is {}'.format(z[0]))
                    #  if z[0] < 0: fout.write(' %d 0 0 %d\n' % (z[0],model))
                    #  else: fout.write(' %d %d 1 %d\n' % (z,z,model))
                #  else:  #JA JB JC format
                    #  print('len(z) = {}'.format(len(z)))
                    #  print('-- check2 --')
                    #  #  fout.write(' %d %d 1 %d\n' % (z[0],z[1],model))
                    #  fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))
            #  # If z is not a list, there must not be a mix of zone specifications 
            #  else:
                #  print('zones is')
                #  print(zones)
                #  print('z is')
                #  print(z)
                #  print('-- check3 --')
                #  if z < 0: fout.write(' %d 0 0 %d\n' % (z,model))
                #  else: fout.write(' %d %d 1 %d\n' % (z,z,model))
    #  #  for m,zs in enumerate(zones):
        #  #  for z in zs:
            #  #  if z < 0: fout.write(' %d 0 0 %d\n' % (z,m+1))
            #  #  else: fout.write(' %d %d 1 %d\n' % (z,z,m+1))
    fout.write('\nstop\n')
    if comment:fout.write('# '+comment)
    fout.close()

#  #DEBUGGING write_userr
#  t = [0.000000, 0.042812, 0.085624, 0.128436, 0.171249]
#  ys = [[4e-5, 1e-4, 3e-4, 5e-5, 1.5e-6],
      #  [1e-4, 2e-4, 3e-4, 4e-4, 5e-4]]
#  zones = [[10201, 10251],
         #  [10150, 10200]]
#  
#  out_file = 'test.dat'
#  varname = 'distcoeffs'
#  
def write_userr(out_file, t, ys, varname='distcoeffs', zones=[7,106]):
    '''
    Create a userr_data.dat file for the userr subroutine in FEHM rxn macro.
    Each 'model' is a different set of distribution coefficents that can be
    applied to different zones/nodes in the domain.
    Uses the same timestep changes for all models.

    Parameters
    ----------
    out_file : str
        Output file name.
    t : ndarray or list(ndarrays)
        Time steps [d] for boundary condition changes.
    ys : ndarray or list(ndarrays) or list(list(ndarrays))
        Variable(s) boundary condition changes. Use list for multiple.
    varname : str or list(str) or list(list(str))
        Variable(s) names.
        Use str for single variable (e.g., "distcoeffs").
        'distcoeffs' currently the only supported option for userr.
    zones : int or list(int) or list(list(int))
        Zone numbers to apply BC changes to. Use list for multiple.
        Use list of lists for multiple models (zones for each model
        should be in a sub-list).
        - e.g., zones = [ [-100],[2101, 2152] ]. This will assign model 1 to
          zone number 100, and model 2 using the JA JB JC format in FEHM
          (default to 1-node spacing: i.e., 2101 2152 1).
        For now, assume you can only have one zone per model per call.
    '''

    err = "ERROR: Problem writing {}\n".format(out_file)
    fout = open(out_file,'w')
    # Check if there are multiple models desired in this .boun
    n_times = 1
    n_models = 1
    #  if type(step)==list:
    if type(ys)==list:
        #  n_models = len(step)
        n_models = len(ys)
        n_times = len(t)
    zones = numpy.atleast_2d(numpy.array(zones))
    # Write num_times num_models to first line
    fout.write('{}   {}\n'.format(n_times,n_models))
    # Make formatting of variables the same regardless of # of models
    t = [t]         # always will only have one set of times
    varname = [varname]*n_models  #repeat distcoeffs for now even though only possibility
    if n_models==1:
        ys = [ys]
    # Check to make sure there are parameters provided for each model
    if np.shape(t)[0]>1: sys.exit(err+"Currently only one time array allowed for all models.")
    elif np.shape(t)[1] != np.shape(ys)[1]: sys.exit(err+"Different number of BC changes from timesteps.")
    #  #  elif len(np.shape(ys)) < n_models: sys.exit(err+"Fewer ys arrays than number of models.")
    #  elif len(step) < n_models: sys.exit(err+"Fewer step calls than number of models.")
    #  elif len(varname) < n_models: sys.exit(err+"Fewer varname calls than number of models.")
    #  #  elif np.shape(zones)[0] < n_models: sys.exit(err+"Fewer zone calls than number of models.")
    #  elif max(np.shape(zones)) < n_models: sys.exit(err+"Fewer zone calls than number of models.")
    # Loop for each model
    for model in range(1,n_models+1):
        yss = numpy.array(ys[model-1]) #ys for current model    #possible problem  ??
        # Make formatting of each model's ys the same
        #-- 1st index=# of models, 2nd index=max # of variable types for all models
        if len(np.shape(yss))==1:
            yss = [yss]
        ts = numpy.array(t[0])  #times for current model (not entirely necessary)
        v = 0 #variable count
        for y in yss:
            # write the timesteps onlyonce for all BC changes
            #  if v==0:
            fout.write('model '+str(model)+'\n')
            if model==1:
                #  fout.write('model '+str(model)+'\n')
                fout.write('times (days)\n')
                #  model += 1
                cnt = 0
                for i in range(len(ts)):
                    # Fix issue with decimal accuracy for big numbers
                    if ts[i] > 1e3:
                        if ts[i] >=1e6: fout.write('{: >14.2e}'.format(ts[i]))
                        # Otherwise, accuracy to 6 decimals (needed for FEHM interp)
                        else: fout.write('{: >14.6f}'.format(ts[i]))
                    else: fout.write('{: >14.6f}'.format(ts[i]))
                    cnt += 1
                    if cnt == 5 or i == n_times-1:
                        fout.write('\n')
                        cnt = 0
            # Write the variable BC changes
            if type(varname[model-1]) is list: fout.write('\n'+varname[model-1][v]+'\n')
            else: fout.write(varname[model-1]+'\n')
            cnt = 0
            for i in range(len(y)):
                fout.write(' %12g' % y[i])
                cnt += 1
                if cnt == 5 or i == len(y)-1:
                    fout.write('\n')
                    cnt = 0
            v+=1
    # Write the zones numbers for each model at the end
    #---- Perhaps for now, assume you can only apply 1 model per zone per call ...
    fout.write('\n')
    zones_store = zones
    try: zones[model-1]
    except IndexError: zones = zones_store[0]
    for model in range(1,n_models+1):
        z = zones[model-1]
        # type(z)=list means a mix of ZONEN and JA JB JC
        if type(z) is list:
            # If entry is a ZONEN specification
            if len(z) == 1:
                if z[0] < 0: fout.write(' %d 0 0 %d\n' % (z[0],model))
                else: fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
            # If entry is a JA JB JC specification
            else:
                fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
        # type(z)=numpy.ndarray  means a just ZONEN
        elif type(z) is numpy.ndarray:
            if z[0] < 0: fout.write(' %d 0 0 %d\n' % (z[0],model))
            # If entry is a JA JB JC specification
            else:
                fout.write(' {} {} 1 {}\n'.format(z[0],z[1],model))#JA JB JC format
    fout.write('\nstop\n')
    fout.close()



def write_pflotran_boun_tpl(out_file, ts, ys):
    ''' From baro_decomp/convert_boun.py'''
    with open(out_file,'w') as fout:
        fout.write('ptf #\n')
        fout.write('TIME_UNITS d\n')
        fout.write('DATA_UNITS Pa\n')
        for t,y in zip(ts,ys):
            fout.write(' %12g #%12g + pressure_offset#\n' % (t,y*1e6))

def write_zone(zonefile, nodes, zone_numbers, zonn=False):
    vs_per_line = 10
    nodes = numpy.atleast_2d(numpy.array(nodes))
    with open(zonefile, 'w') as fout:
        if zonn: fout.write('zonn\n')
        else: fout.write('zone\n')
        for nmbr,ns in zip(zone_numbers,nodes):
            #  fout.write('    %d\n'%nmbr)
            fout.write('    %d\n'%abs(nmbr)) #make sure zone_number is positive
            fout.write('nnum\n')
            fout.write('    %d\n'%len(ns))
            for i in range(0, len(ns), vs_per_line):
                fout.write("    "+" ".join([str(v) for v in ns[i:i+vs_per_line]])+'\n')
        fout.write('\nstop\n')

def fourier_series(periods,ampls,shifts,mean,ts):
    '''
    Generate fourier series

    :param periods: Periods of sine waves
    :type periods: list(floats)
    :param ampls: Amplitudes of sine waves
    :type ampls: list(floats)
    :param shifts: Shifts of sine waves
    :type shifts: list(floats)
    :param mean: Mean value of signal
    :type mean: float
    :param ts: Times at which to evaluate fourier series
    :type ts: list(floats)
    :returns: fourier series values each time in ts
    '''
    y_synth = mean
    for p,a,phi in zip(periods,ampls,shifts):
        w = 2*numpy.pi/p # Frequency
        y_synth -= a*(numpy.sin(w*ts+phi))
    return y_synth


if __name__=="__main__":
    km = 1.e-13
    k_bg = 1.2928e-018 # Background perm
    ke_ccc = 1.e-7 # effective perm of chimney
    distdiv = 10
    maxfracap = 0.0013407
    meshdir = '/scratch/nts/jordan/grid_gen/une_fracnet/mesh/tuff125m/'
    matfile = '/scratch/fwo/jordan/une/fracnetfiles/matlabfiles/need_for_all_tuff125.mat'
    init_dir = '/scratch/er/dharp/projects/gas_seepage/anchorage_smooth/tuff/init'

    g = grid_recarray(meshdir,maxfracap)
    gdkm_data = gdkm_macro(g,matfile,km,ke_ccc,distdiv,macro_dir='macros')
    node_macro(g,macro_dir='macros')
    perm_macro(g,gdkm_data,k_bg,ke_ccc,macro_dir='macros')
    flow_macro(g,init_dir,macro_dir='macros')

    if not os.path.isdir('spinup'):
        os.mkdir('spinup')
    os.chdir('spinup')
    for f in glob('../macros/*'):
        lnk = f.split(os.sep)[-1]
        if not os.path.isfile(lnk):
            os.symlink(f,lnk)
    os.chdir('..')
    if not os.path.isdir('tracer'):
        os.mkdir('tracer')
    os.chdir('tracer')
    for f in glob('../macros/*'):
        lnk = f.split(os.sep)[-1]
        if not os.path.isfile(lnk):
            os.symlink(f,lnk)
    os.chdir('..')
