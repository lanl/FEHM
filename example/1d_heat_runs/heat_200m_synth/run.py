import os
import pickle
import plotting

# reads exe and others from master_model.py
#  exec(open('../master_model.py').read())
exec(open('../heat_master_model.py').read())
# # Alternatively, you can provide an exe file here:
#  exe = '/home/jportiz/software/FEHM/src/xfehm_v3.3.1_mars'

# for testing...
#  if __name__=='__main__':
p = {'mesh'               : 1,         #1D, 0.1x200m mesh
     'sim_time'           : 52000,     #[sols]
     'geothermal_gradient': 0.012908,  #[°C/m]
     'thermal_cond'       : 2.0,       #[W/m/K]
     'T_ref'              : -46.93,    #[°C] avg surf temp @ Gale [Klusman2022]
     'synthetic'          : True,      #create synthetic pressure/temperature record
    }
# 'geothermal_gradient': 0.012908 # [Klusman2023]
# 'geothermal_gradient': 0.004500 # [Jones2011]
model(p,exe=exe,verbose=True)

#-------------------------------------------------- 
#     Run plotting scripts
#-------------------------------------------------- 
# Load saved python variables
with open('v.pkl','rb') as f:
    pardict,geothermal_gradient,T_ref,T_0,T_upshift = pickle.load(f)
# Create output/ directory if doesn't exist
if not os.path.exists('output'): os.makedirs('output')
#-- Plot geothermal gradient from init0 sim
print('Plotting geothermal gradient...')
plotting.plot_geothermal_gradient(geothermal_gradient,T_upshift,T_ref)
print('    Done.')
#-- Plot timeseries of subsurface temperatures
print('Plotting subsurface temperature timeseries...')
plotting.plot_timeseries(T_0=T_0,T_upshift=T_upshift)
print('    Done.')
#-- Plot timeseries of subsurface temperatures
print('Plotting subsurface temperature ranges...')
plotting.plot_ranges(T_0=T_0,T_upshift=T_upshift)
print('    Done.')

