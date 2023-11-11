import os,sys
import pickle
sys.path.append('/project/gas_seepage/jportiz/scripts/')
from mars_tools import sols2days
import plot

# reads exe and others rom master_model.py
exec(open('../master_model.py').read())
# # Alternatively, you can provide an exe file here:
#  exe = '/home/jportiz/software/FEHM/src/xfehm_v3.3.1_mars'

# for testing...
#  if __name__=='__main__':
p = {'mesh': 'depth200_fracDen035',    # levy 50x200m fracture density = 0.035% (@b=1mm)
     'sim_time': sols2days(668*75),#52000,   # [days]
     'methane_prod_rate': 6.67e-16,    # [mol/s/m^3]
     'c_initial': 2.8106e-6,           # [mol/kg vapor]
     'geothermal_gradient': 0.012908,  # [°C/m]   #probably don't need  here
     'thermal_cond'       : 2.0,       # [W/m/K]  #here, it can probably be 0.0 (no heat flow)
     'active_regolith_base'      : 5.0, #50.0,# [m] depth to base of thermally active adsorptive regolith (positive depth)
     'constant_regolith_base'      : 'bottom', #30.0, # [m] depth to base of constant adsorptive regolith (positive depth)
     #  'source_depth_top'   : 5,         # [m]  default is bottom of domain
     'T_ref'              : -46.93,    # [°C] mean surf. temp @ Gale
     'df'                 : 0.001,     # [m] fracture aperture (bigger, 1mm)
     'synthetic'          : True,      # use synthetic P and T?
     #  'noSorp_dir'         : 'noSorp_depth200_b1mm_fracDen035',  #for now, can restart from non-synthetic record sims
     #  'no_adsorption'      : True,      #default is False
    }
# 'geothermal_gradient': 0.012908 # [Klusman2023]
# 'geothermal_gradient': 0.004500 # [Jones2011]
print('tpl_dir = {}'.format(tpl_dir))
model(p,exe=exe,tpl_dir=tpl_dir,verbose=True)
#  model(p,verbose=True)

#---- Run plotting script
with open('v.pkl','rb') as f:
    sourcetype,fehm_upshift,T_upshift = pickle.load(f)
plot.make_plots(fehm_upshift=fehm_upshift, T_upshift=T_upshift)
