"""Functions for FEHM thermodynamic variables calculations."""

"""
Copyright 2013.
Los Alamos National Security, LLC. 
This material was produced under U.S. Government contract DE-AC52-06NA25396 for 
Los Alamos National Laboratory (LANL), which is operated by Los Alamos National 
Security, LLC for the U.S. Department of Energy. The U.S. Government has rights 
to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES 
ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce 
derivative works, such modified software should be clearly marked, so as not to 
confuse it with the version available from LANL.

Additionally, this library is free software; you can redistribute it and/or modify 
it under the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation; either version 2.1 of the License, or (at your option) 
any later version. Accordingly, this library is distributed in the hope that it 
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General 
Public License for more details.
"""

import numpy as np
from ftool import*

from fdflt import*
dflt = fdflt()

YEL=[ 0.25623465e-03, 0.10184405e-02, 0.22554970e-04, 0.34836663e-07, 0.41769866e-02, -0.21244879e-04, 
   0.25493516e-07, 0.89557885e-04, 0.10855046e-06, -0.21720560e-06,]
YEV=[ 0.31290881e+00, -0.10e+01, 0.25748596e-01, 0.38846142e-03, 0.11319298e-01, 0.20966376e-04, 
   0.74228083e-08, 0.19206133e-02, -0.10372453e-03, 0.59104245e-07,]
YDL=[ 0.10000000e+01, 0.17472599e-01, -0.20443098e-04, -0.17442012e-06, 0.49564109e-02, -0.40757664e-04, 
   0.50676664e-07, 0.50330978e-04, 0.33914814e-06, -0.18383009e-06,] 
YDV=[ 0.15089524e-05, 0.10000000e+01, -0.10000000e+01, -0.16676705e-02, 0.40111210e-07, 0.25625316e-10, 
  -0.40479650e-12, 0.43379623e-01, 0.24991800e-02, -0.94755043e-04,]
YVL=[ 0.17409149e-02, 0.18894882e-04, -0.66439332e-07, -0.23122388e-09, -0.31534914e-05, 0.11120716e-07, 
  -0.48576020e-10, 0.28006861e-07, 0.23225035e-09, 0.47180171e-10,] 
YVV=[-0.13920783e-03, 0.98434337e-02, -0.51504232e-03, 0.62554603e-04, 0.27105772e-04, 0.84981906e-05, 
   0.34539757e-07, -0.25524682e-03, 0, 0.12316788e-05,] 
 
ZEL=[ 0.10000000e+01, 0.23513278e-01, 0.48716386e-04, -0.19935046e-08, -0.50770309e-02, 0.57780287e-05, 
   0.90972916e-09, -0.58981537e-04, -0.12990752e-07, 0.45872518e-08,]
ZEV=[ 0.12511319e+00, -0.36061317e+00, 0.58668929e-02, 0.99059715e-04, 0.44331611e-02, 0.50902084e-05, 
  -0.10812602e-08, 0.90918809e-03, -0.26960555e-04, -0.36454880e-06,]
ZDL=[ 0.10009476e-02, 0.16812589e-04, -0.24582622e-07, -0.17014984e-09, 0.48841156e-05, -0.32967985e-07, 
   0.28619380e-10, 0.53249055e-07, 0.30456698e-09, -0.12221899e-09,]
ZDV=[ 0.12636224e+00, -0.30463489e+00, 0.27981880e-02, 0.51132337e-05, 0.59318010e-02, 0.80972509e-05, 
  -0.43798358e-07, 0.53046787e-03, -0.84916607e-05, 0.48444919e-06,]
ZVL=[ 0.10000000e+01, 0.10523153e-01, -0.22658391e-05, -0.31796607e-06, 0.29869141e-01, 0.21844248e-03, 
  -0.87658855e-06, 0.41690362e-03, -0.25147022e-05, 0.22144660e-05,]
ZVV=[ 0.10000000e+01, 0.10000000e+01, -0.10e1 , -0.10e1 , 0.10000000e+01,  0.0000000e+01,
  -0.22934622e-03, 0.10000000e+01, 0 , 0.25834551e-01,]    
  
YSP=[ 0.71725602e-03, 0.22607516e-04, 0.26178556e-05, -0.10516335e-07, 0.63167028e-09,]
YST=[-0.25048121e-05, 0.45249584e-02, 0.33551528e+00, 0.10000000e+01, 0.12254786e+00,]

ZSP=[ 0.10000000e+01, -0.22460012e-02, 0.30234492e-05, -0.32466525e-09, 0.0,]
ZST=[0.20889841e-06, 0.11587544e-03, 0.31934455e-02, 0.45538151e-02, 0.23756593e-03,]  

if WINDOWS:
	co2_interp_path = dflt.co2_interp_path
else:
	co2_interp_path = dflt.co2_interp_path

co2Vars = False
if co2_interp_path != '' and os.path.isfile(co2_interp_path): co2Vars = True

if co2Vars:
	with open(co2_interp_path,'r') as f:
		f.readline()
		line = f.readline()
		tn,pn,na = line.split()[:3]
		tn = int(tn); pn = int(pn); na = int(na)
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		# read in temperature data
		keepReading = True
		T = []
		while keepReading:
			line = f.readline()
			if '>' in line: break
			T.append(line.strip().split())
		T = list(flatten(T))
		Tval = np.array([float(t) for t in T])
		Tdict = dict([(t,i) for i,t in enumerate(Tval)])
		# read in pressure data
		P = []
		while keepReading:
			line = f.readline()
			if '>' in line: break
			P.append(line.strip().split())
		P = list(flatten(P))
		Pval = np.array([float(p) for p in P])	
		Pdict = dict([(p,i) for i,p in enumerate(Pval)])
		# read to array data
		while keepReading:
			line = f.readline()
			if '>' in line: break
		# read in array data
		arraynames = ['density','dddt','dddp','enthalpy','dhdt','dhdp','viscosity','dvdt','dvdp']
		arrays = {}
		for arrayname in arraynames:
			array = []
			while keepReading:
				line = f.readline()
				if '>' in line: break
				array.append(line.strip().split())
			array = list(flatten(array))
			array = np.array([float(a) for a in array])	
			arrays.update(dict(((arrayname,array),)))
	
def dens(P,T,derivative=''):
	"""Return liquid water, vapor water and CO2 density, or derivatives with respect to temperature or pressure, for specified temperature and pressure.
	
	:param P: Pressure (MPa). 
	:type P: fl64
	:param T: Temperature (degC)
	:type T: fl64
	:param derivative: Supply 'T' or 'temperature' for derivatives with respect to temperature, or 'P' or 'pressure' for derivatives with respect to pressure.
	:type T: str
	:returns: Three element tuple containing (liquid, vapor, CO2) density or derivatives if requested.
	
	"""
	P = np.array(P)
	T = np.array(T)
	# calculate water properties
	if not derivative:
		YL0 = (YDL[0] +
		       YDL[1]*P +
		       YDL[2]*P**2 +
		       YDL[3]*P**3 +
		       YDL[4]*T +
		       YDL[5]*T**2 +
		       YDL[6]*T**3 +
		       YDL[7]*P*T +
		       YDL[8]*P**2*T +
		       YDL[9]*P*T**2)
		ZL0 = (ZDL[0] +
		       ZDL[1]*P +
		       ZDL[2]*P**2 +
		       ZDL[3]*P**3 +
		       ZDL[4]*T +
		       ZDL[5]*T**2 +
		       ZDL[6]*T**3 +
		       ZDL[7]*P*T +
		       ZDL[8]*P**2*T +
		       ZDL[9]*P*T**2)
		YV0 = (YDV[0] +
		       YDV[1]*P +
		       YDV[2]*P**2 +
		       YDV[3]*P**3 +
		       YDV[4]*T +
		       YDV[5]*T**2 +
		       YDV[6]*T**3 +
		       YDV[7]*P*T +
		       YDV[8]*P**2*T +
		       YDV[9]*P*T**2)
		ZV0 = (ZDV[0] +
		       ZDV[1]*P +
		       ZDV[2]*P**2 +
		       ZDV[3]*P**3 +
		       ZDV[4]*T +
		       ZDV[5]*T**2 +
		       ZDV[6]*T**3 +
		       ZDV[7]*P*T +
		       ZDV[8]*P**2*T +
		       ZDV[9]*P*T**2)
		dens_l = YL0/ZL0
		dens_v = YV0/ZV0
	elif derivative in ['P','pressure']:
		# terms
		YL0 = (YDL[0] +
		       YDL[1]*P +
		       YDL[2]*P**2 +
		       YDL[3]*P**3 +
		       YDL[4]*T +
		       YDL[5]*T**2 +
		       YDL[6]*T**3 +
		       YDL[7]*P*T +
		       YDL[8]*P**2*T +
		       YDL[9]*P*T**2)
		ZL0 = (ZDL[0] +
		       ZDL[1]*P +
		       ZDL[2]*P**2 +
		       ZDL[3]*P**3 +
		       ZDL[4]*T +
		       ZDL[5]*T**2 +
		       ZDL[6]*T**3 +
		       ZDL[7]*P*T +
		       ZDL[8]*P**2*T +
		       ZDL[9]*P*T**2)
		YV0 = (YDV[0] +
		       YDV[1]*P +
		       YDV[2]*P**2 +
		       YDV[3]*P**3 +
		       YDV[4]*T +
		       YDV[5]*T**2 +
		       YDV[6]*T**3 +
		       YDV[7]*P*T +
		       YDV[8]*P**2*T +
		       YDV[9]*P*T**2)
		ZV0 = (ZDV[0] +
		       ZDV[1]*P +
		       ZDV[2]*P**2 +
		       ZDV[3]*P**3 +
		       ZDV[4]*T +
		       ZDV[5]*T**2 +
		       ZDV[6]*T**3 +
		       ZDV[7]*P*T +
		       ZDV[8]*P**2*T +
		       ZDV[9]*P*T**2)
		# derivatives
		YL1 = (YDL[1] +
		       YDL[2]*P*2 +
		       YDL[3]*P**2*3 +
		       YDL[7]*T +
		       YDL[8]*P*2*T +
		       YDL[9]*T**2)
		ZL1 = (ZDL[1] +
		       ZDL[2]*P*2 +
		       ZDL[3]*P**2*3 +
		       ZDL[7]*T +
		       ZDL[8]*P*2*T +
		       ZDL[9]*T**2)
		YV1 = (YDV[1] +
		       YDV[2]*P*2 +
		       YDV[3]*P**2*3 +
		       YDV[7]*T +
		       YDV[8]*P*2*T +
		       YDV[9]*T**2)
		ZV1 = (ZDV[1] +
		       ZDV[2]*P*2 +
		       ZDV[3]*P**2*3 +
		       ZDV[7]*T +
		       ZDV[8]*P*2*T +
		       ZDV[9]*T**2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	elif derivative in ['T','temperature']:
		# terms
		YL0 = (YDL[0] +
		       YDL[1]*P +
		       YDL[2]*P**2 +
		       YDL[3]*P**3 +
		       YDL[4]*T +
		       YDL[5]*T**2 +
		       YDL[6]*T**3 +
		       YDL[7]*P*T +
		       YDL[8]*P**2*T +
		       YDL[9]*P*T**2)
		ZL0 = (ZDL[0] +
		       ZDL[1]*P +
		       ZDL[2]*P**2 +
		       ZDL[3]*P**3 +
		       ZDL[4]*T +
		       ZDL[5]*T**2 +
		       ZDL[6]*T**3 +
		       ZDL[7]*P*T +
		       ZDL[8]*P**2*T +
		       ZDL[9]*P*T**2)
		YV0 = (YDV[0] +
		       YDV[1]*P +
		       YDV[2]*P**2 +
		       YDV[3]*P**3 +
		       YDV[4]*T +
		       YDV[5]*T**2 +
		       YDV[6]*T**3 +
		       YDV[7]*P*T +
		       YDV[8]*P**2*T +
		       YDV[9]*P*T**2)
		ZV0 = (ZDV[0] +
		       ZDV[1]*P +
		       ZDV[2]*P**2 +
		       ZDV[3]*P**3 +
		       ZDV[4]*T +
		       ZDV[5]*T**2 +
		       ZDV[6]*T**3 +
		       ZDV[7]*P*T +
		       ZDV[8]*P**2*T +
		       ZDV[9]*P*T**2)
		# derivatives
		YL1 = (YDL[4] +
		       YDL[5]*T*2 +
		       YDL[6]*T**2*3 +
		       YDL[7]*P +
		       YDL[8]*P**2 +
		       YDL[9]*P*T*2)
		ZL1 = (ZDL[4] +
		       ZDL[5]*T*2 +
		       ZDL[6]*T**2*3 +
		       ZDL[7]*P +
		       ZDL[8]*P**2 +
		       ZDL[9]*P*T*2)
		YV1 = (YDV[4] +
		       YDV[5]*T*2 +
		       YDV[6]*T**2*3 +
		       YDV[7]*P +
		       YDV[8]*P**2 +
		       YDV[9]*P*T*2)
		ZV1 = (ZDV[4] +
		       ZDV[5]*T*2 +
		       ZDV[6]*T**2*3 +
		       ZDV[7]*P +
		       ZDV[8]*P**2 +
		       ZDV[9]*P*T*2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	else: print 'not a valid derivative'; return
	
	if not co2Vars: return (dens_l,dens_v,np.array([]))
	
	# calculate co2 properties
	if not derivative: k = 'density'
	elif derivative in ['P','pressure']: k = 'dddp'
	elif derivative in ['T','temperature']: k = 'dddt'
	
	if not P.shape: P = np.array([P])
	if not T.shape: T = np.array([T])
	if P.size == 1 and not T.size == 1: P = P*np.ones((1,len(T)))[0]
	elif T.size == 1 and not P.size == 1: T = T*np.ones((1,len(P)))[0]
	# calculate bounding values of P
	dens_c = []
	for Pi, Ti in zip(P,T):
		if Pi<=Pval[0]: p0 = Pval[0]; p1 = Pval[0]
		elif Pi>=Pval[-1]: p0 = Pval[-1]; p1 = Pval[-1]
		else:
			p0 = Pval[0]
			for p in Pval[1:]:
				if Pi<=p: 
					p1 = p; break
				else: 
					p0 = p
		#print Pi,Pval[0],Pval[-1],p0,p1
		# calculate bounding values of T
		if Ti<=Tval[0]: t0 = Tval[0]; t1 = Tval[0]
		elif Ti>=Tval[-1]: t0 = Tval[-1]; t1 = Tval[-1]
		else:
			t0 = Tval[0]
			for t in Tval[1:]:
				if Ti<=t: 
					t1 = t; break
				else: 
					t0 = t
		# calculate four indices
		t0 = Tdict[t0]; t1 = Tdict[t1]; p0 = Pdict[p0]; p1 = Pdict[p1]
		i1 = p0*tn+t0
		i2 = p0*tn+t1
		i3 = p1*tn+t0
		i4 = p1*tn+t1
		# locate value in array
		v1 = arrays[k][i1]
		v2 = arrays[k][i2]
		v3 = arrays[k][i3]
		v4 = arrays[k][i4]
		if p0 == 0 and p1 == 0:
			dens_c.append((0.5*(t1*v3+t0*v4)/(t1+t0) + 0.5*(t1*v1+t0*v2)/(t1+t0)))
		else:
			dens_c.append((p0*(t1*v3+t0*v4)/(t1+t0) + p1*(t1*v1+t0*v2)/(t1+t0))/(p0+p1))
	
	return (dens_l,dens_v,np.array(dens_c))
def enth(P,T,derivative=''):
	"""Return liquid water, vapor water and CO2 enthalpy, or derivatives with respect to temperature or pressure, for specified temperature and pressure.
	
	:param P: Pressure (MPa). 
	:type P: fl64
	:param T: Temperature (degC)
	:type T: fl64
	:param derivative: Supply 'T' or 'temperature' for derivatives with respect to temperature, or 'P' or 'pressure' for derivatives with respect to pressure.
	:type T: str
	:returns: Three element tuple containing (liquid, vapor, CO2) enthalpy or derivatives if requested.
	
	"""
	P = np.array(P)
	T = np.array(T)
	if not derivative:
		YL0 = (YEL[0] +
		       YEL[1]*P +
		       YEL[2]*P**2 +
		       YEL[3]*P**3 +
		       YEL[4]*T +
		       YEL[5]*T**2 +
		       YEL[6]*T**3 +
		       YEL[7]*P*T +
		       YEL[8]*P**2*T +
		       YEL[9]*P*T**2)
		ZL0 = (ZEL[0] +
		       ZEL[1]*P +
		       ZEL[2]*P**2 +
		       ZEL[3]*P**3 +
		       ZEL[4]*T +
		       ZEL[5]*T**2 +
		       ZEL[6]*T**3 +
		       ZEL[7]*P*T +
		       ZEL[8]*P**2*T +
		       ZEL[9]*P*T**2)
		YV0 = (YEV[0] +
		       YEV[1]*P +
		       YEV[2]*P**2 +
		       YEV[3]*P**3 +
		       YEV[4]*T +
		       YEV[5]*T**2 +
		       YEV[6]*T**3 +
		       YEV[7]*P*T +
		       YEV[8]*P**2*T +
		       YEV[9]*P*T**2)
		ZV0 = (ZEV[0] +
		       ZEV[1]*P +
		       ZEV[2]*P**2 +
		       ZEV[3]*P**3 +
		       ZEV[4]*T +
		       ZEV[5]*T**2 +
		       ZEV[6]*T**3 +
		       ZEV[7]*P*T +
		       ZEV[8]*P**2*T +
		       ZEV[9]*P*T**2)
		dens_l = YL0/ZL0
		dens_v = YV0/ZV0
	elif derivative in ['P','pressure']:
		# terms
		YL0 = (YEL[0] +
		       YEL[1]*P +
		       YEL[2]*P**2 +
		       YEL[3]*P**3 +
		       YEL[4]*T +
		       YEL[5]*T**2 +
		       YEL[6]*T**3 +
		       YEL[7]*P*T +
		       YEL[8]*P**2*T +
		       YEL[9]*P*T**2)
		ZL0 = (ZEL[0] +
		       ZEL[1]*P +
		       ZEL[2]*P**2 +
		       ZEL[3]*P**3 +
		       ZEL[4]*T +
		       ZEL[5]*T**2 +
		       ZEL[6]*T**3 +
		       ZEL[7]*P*T +
		       ZEL[8]*P**2*T +
		       ZEL[9]*P*T**2)
		YV0 = (YEV[0] +
		       YEV[1]*P +
		       YEV[2]*P**2 +
		       YEV[3]*P**3 +
		       YEV[4]*T +
		       YEV[5]*T**2 +
		       YEV[6]*T**3 +
		       YEV[7]*P*T +
		       YEV[8]*P**2*T +
		       YEV[9]*P*T**2)
		ZV0 = (ZEV[0] +
		       ZEV[1]*P +
		       ZEV[2]*P**2 +
		       ZEV[3]*P**3 +
		       ZEV[4]*T +
		       ZEV[5]*T**2 +
		       ZEV[6]*T**3 +
		       ZEV[7]*P*T +
		       ZEV[8]*P**2*T +
		       ZEV[9]*P*T**2)
		# derivatives
		YL1 = (YEL[1] +
		       YEL[2]*P*2 +
		       YEL[3]*P**2*3 +
		       YEL[7]*T +
		       YEL[8]*P*2*T +
		       YEL[9]*T**2)
		ZL1 = (ZEL[1] +
		       ZEL[2]*P*2 +
		       ZEL[3]*P**2*3 +
		       ZEL[7]*T +
		       ZEL[8]*P*2*T +
		       ZEL[9]*T**2)
		YV1 = (YEV[1] +
		       YEV[2]*P*2 +
		       YEV[3]*P**2*3 +
		       YEV[7]*T +
		       YEV[8]*P*2*T +
		       YEV[9]*T**2)
		ZV1 = (ZEV[1] +
		       ZEV[2]*P*2 +
		       ZEV[3]*P**2*3 +
		       ZEV[7]*T +
		       ZEV[8]*P*2*T +
		       ZEV[9]*T**2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	elif derivative in ['T','temperature']:
		# terms
		YL0 = (YEL[0] +
		       YEL[1]*P +
		       YEL[2]*P**2 +
		       YEL[3]*P**3 +
		       YEL[4]*T +
		       YEL[5]*T**2 +
		       YEL[6]*T**3 +
		       YEL[7]*P*T +
		       YEL[8]*P**2*T +
		       YEL[9]*P*T**2)
		ZL0 = (ZEL[0] +
		       ZEL[1]*P +
		       ZEL[2]*P**2 +
		       ZEL[3]*P**3 +
		       ZEL[4]*T +
		       ZEL[5]*T**2 +
		       ZEL[6]*T**3 +
		       ZEL[7]*P*T +
		       ZEL[8]*P**2*T +
		       ZEL[9]*P*T**2)
		YV0 = (YEV[0] +
		       YEV[1]*P +
		       YEV[2]*P**2 +
		       YEV[3]*P**3 +
		       YEV[4]*T +
		       YEV[5]*T**2 +
		       YEV[6]*T**3 +
		       YEV[7]*P*T +
		       YEV[8]*P**2*T +
		       YEV[9]*P*T**2)
		ZV0 = (ZEV[0] +
		       ZEV[1]*P +
		       ZEV[2]*P**2 +
		       ZEV[3]*P**3 +
		       ZEV[4]*T +
		       ZEV[5]*T**2 +
		       ZEV[6]*T**3 +
		       ZEV[7]*P*T +
		       ZEV[8]*P**2*T +
		       ZEV[9]*P*T**2)
		# derivatives
		YL1 = (YEL[4] +
		       YEL[5]*T*2 +
		       YEL[6]*T**2*3 +
		       YEL[7]*P +
		       YEL[8]*P**2 +
		       YEL[9]*P*T*2)
		ZL1 = (ZEL[4] +
		       ZEL[5]*T*2 +
		       ZEL[6]*T**2*3 +
		       ZEL[7]*P +
		       ZEL[8]*P**2 +
		       ZEL[9]*P*T*2)
		YV1 = (YEV[4] +
		       YEV[5]*T*2 +
		       YEV[6]*T**2*3 +
		       YEV[7]*P +
		       YEV[8]*P**2 +
		       YEV[9]*P*T*2)
		ZV1 = (ZEV[4] +
		       ZEV[5]*T*2 +
		       ZEV[6]*T**2*3 +
		       ZEV[7]*P +
		       ZEV[8]*P**2 +
		       ZEV[9]*P*T*2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	else: print 'not a valid derivative'; return
	
	if not co2Vars: return (dens_l,dens_v,np.array([]))
		
	# calculate co2 properties
	if not derivative: k = 'enthalpy'
	elif derivative in ['P','pressure']: k = 'dhdp'
	elif derivative in ['T','temperature']: k = 'dhdt'
	
	if not P.shape: P = np.array([P])
	if not T.shape: T = np.array([T])
	if P.size == 1 and not T.size == 1: P = P*np.ones((1,len(T)))[0]
	elif T.size == 1 and not P.size == 1: T = T*np.ones((1,len(P)))[0]
	# calculate bounding values of P
	dens_c = []
	for Pi, Ti in zip(P,T):
		if Pi<=Pval[0]: p0 = Pval[0]; p1 = Pval[0]
		elif Pi>=Pval[-1]: p0 = Pval[-1]; p1 = Pval[-1]
		else:
			p0 = Pval[0]
			for p in Pval[1:]:
				if Pi<=p: 
					p1 = p; break
				else: 
					p0 = p
		# calculate bounding values of T
		if Ti<=Tval[0]: t0 = Tval[0]; t1 = Tval[0]
		elif Ti>=Tval[-1]: t0 = Tval[-1]; t1 = Tval[-1]
		else:
			t0 = Tval[0]
			for t in Tval[1:]:
				if Ti<=t: 
					t1 = t; break
				else: 
					t0 = t
		# calculate four indices
		dt0 = abs(Ti-t0); dt1 = abs(Ti-t1); dp0 = abs(Pi-p0); dp1 = abs(Pi-p1)
		t0 = Tdict[t0]; t1 = Tdict[t1]; p0 = Pdict[p0]; p1 = Pdict[p1]
		i1 = p0*tn+t0
		i2 = p0*tn+t1
		i3 = p1*tn+t0
		i4 = p1*tn+t1
		# locate value in array
		v1 = arrays[k][i1]
		v2 = arrays[k][i2]
		v3 = arrays[k][i3]
		v4 = arrays[k][i4]
		dens_c.append((p0*(t1*v3+t0*v4)/(t1+t0) + p1*(t1*v1+t0*v2)/(t1+t0))/(p0+p1))
	
	return (dens_l,dens_v,np.array(dens_c))
def visc(P,T,derivative=''):
	"""Return liquid water, vapor water and CO2 viscosity, or derivatives with respect to temperature or pressure, for specified temperature and pressure.
	
	:param P: Pressure (MPa). 
	:type P: fl64
	:param T: Temperature (degC)
	:type T: fl64
	:param derivative: Supply 'T' or 'temperature' for derivatives with respect to temperature, or 'P' or 'pressure' for derivatives with respect to pressure.
	:type T: str
	:returns: Three element tuple containing (liquid, vapor, CO2) viscosity or derivatives if requested.
	
	"""
	P = np.array(P)
	T = np.array(T)
	if not derivative:
		YL0 = (YVL[0] +
		       YVL[1]*P +
		       YVL[2]*P**2 +
		       YVL[3]*P**3 +
		       YVL[4]*T +
		       YVL[5]*T**2 +
		       YVL[6]*T**3 +
		       YVL[7]*P*T +
		       YVL[8]*P**2*T +
		       YVL[9]*P*T**2)
		ZL0 = (ZVL[0] +
		       ZVL[1]*P +
		       ZVL[2]*P**2 +
		       ZVL[3]*P**3 +
		       ZVL[4]*T +
		       ZVL[5]*T**2 +
		       ZVL[6]*T**3 +
		       ZVL[7]*P*T +
		       ZVL[8]*P**2*T +
		       ZVL[9]*P*T**2)
		YV0 = (YVV[0] +
		       YVV[1]*P +
		       YVV[2]*P**2 +
		       YVV[3]*P**3 +
		       YVV[4]*T +
		       YVV[5]*T**2 +
		       YVV[6]*T**3 +
		       YVV[7]*P*T +
		       YVV[8]*P**2*T +
		       YVV[9]*P*T**2)
		ZV0 = (ZVV[0] +
		       ZVV[1]*P +
		       ZVV[2]*P**2 +
		       ZVV[3]*P**3 +
		       ZVV[4]*T +
		       ZVV[5]*T**2 +
		       ZVV[6]*T**3 +
		       ZVV[7]*P*T +
		       ZVV[8]*P**2*T +
		       ZVV[9]*P*T**2)
		dens_l = YL0/ZL0
		dens_v = YV0/ZV0
	elif derivative in ['P','pressure']:
		# terms
		YL0 = (YVL[0] +
		       YVL[1]*P +
		       YVL[2]*P**2 +
		       YVL[3]*P**3 +
		       YVL[4]*T +
		       YVL[5]*T**2 +
		       YVL[6]*T**3 +
		       YVL[7]*P*T +
		       YVL[8]*P**2*T +
		       YVL[9]*P*T**2)
		ZL0 = (ZVL[0] +
		       ZVL[1]*P +
		       ZVL[2]*P**2 +
		       ZVL[3]*P**3 +
		       ZVL[4]*T +
		       ZVL[5]*T**2 +
		       ZVL[6]*T**3 +
		       ZVL[7]*P*T +
		       ZVL[8]*P**2*T +
		       ZVL[9]*P*T**2)
		YV0 = (YVV[0] +
		       YVV[1]*P +
		       YVV[2]*P**2 +
		       YVV[3]*P**3 +
		       YVV[4]*T +
		       YVV[5]*T**2 +
		       YVV[6]*T**3 +
		       YVV[7]*P*T +
		       YVV[8]*P**2*T +
		       YVV[9]*P*T**2)
		ZV0 = (ZVV[0] +
		       ZVV[1]*P +
		       ZVV[2]*P**2 +
		       ZVV[3]*P**3 +
		       ZVV[4]*T +
		       ZVV[5]*T**2 +
		       ZVV[6]*T**3 +
		       ZVV[7]*P*T +
		       ZVV[8]*P**2*T +
		       ZVV[9]*P*T**2)
		# derivatives
		YL1 = (YVL[1] +
		       YVL[2]*P*2 +
		       YVL[3]*P**2*3 +
		       YVL[7]*T +
		       YVL[8]*P*2*T +
		       YVL[9]*T**2)
		ZL1 = (ZVL[1] +
		       ZVL[2]*P*2 +
		       ZVL[3]*P**2*3 +
		       ZVL[7]*T +
		       ZVL[8]*P*2*T +
		       ZVL[9]*T**2)
		YV1 = (YVV[1] +
		       YVV[2]*P*2 +
		       YVV[3]*P**2*3 +
		       YVV[7]*T +
		       YVV[8]*P*2*T +
		       YVV[9]*T**2)
		ZV1 = (ZVV[1] +
		       ZVV[2]*P*2 +
		       ZVV[3]*P**2*3 +
		       ZVV[7]*T +
		       ZVV[8]*P*2*T +
		       ZVV[9]*T**2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	elif derivative in ['T','temperature']:
		# terms
		YL0 = (YVL[0] +
		       YVL[1]*P +
		       YVL[2]*P**2 +
		       YVL[3]*P**3 +
		       YVL[4]*T +
		       YVL[5]*T**2 +
		       YVL[6]*T**3 +
		       YVL[7]*P*T +
		       YVL[8]*P**2*T +
		       YVL[9]*P*T**2)
		ZL0 = (ZVL[0] +
		       ZVL[1]*P +
		       ZVL[2]*P**2 +
		       ZVL[3]*P**3 +
		       ZVL[4]*T +
		       ZVL[5]*T**2 +
		       ZVL[6]*T**3 +
		       ZVL[7]*P*T +
		       ZVL[8]*P**2*T +
		       ZVL[9]*P*T**2)
		YV0 = (YVV[0] +
		       YVV[1]*P +
		       YVV[2]*P**2 +
		       YVV[3]*P**3 +
		       YVV[4]*T +
		       YVV[5]*T**2 +
		       YVV[6]*T**3 +
		       YVV[7]*P*T +
		       YVV[8]*P**2*T +
		       YVV[9]*P*T**2)
		ZV0 = (ZVV[0] +
		       ZVV[1]*P +
		       ZVV[2]*P**2 +
		       ZVV[3]*P**3 +
		       ZVV[4]*T +
		       ZVV[5]*T**2 +
		       ZVV[6]*T**3 +
		       ZVV[7]*P*T +
		       ZVV[8]*P**2*T +
		       ZVV[9]*P*T**2)
		# derivatives
		YL1 = (YVL[4] +
		       YVL[5]*T*2 +
		       YVL[6]*T**2*3 +
		       YVL[7]*P +
		       YVL[8]*P**2 +
		       YVL[9]*P*T*2)
		ZL1 = (ZVL[4] +
		       ZVL[5]*T*2 +
		       ZVL[6]*T**2*3 +
		       ZVL[7]*P +
		       ZVL[8]*P**2 +
		       ZVL[9]*P*T*2)
		YV1 = (YVV[4] +
		       YVV[5]*T*2 +
		       YVV[6]*T**2*3 +
		       YVV[7]*P +
		       YVV[8]*P**2 +
		       YVV[9]*P*T*2)
		ZV1 = (ZVV[4] +
		       ZVV[5]*T*2 +
		       ZVV[6]*T**2*3 +
		       ZVV[7]*P +
		       ZVV[8]*P**2 +
		       ZVV[9]*P*T*2)
		dens_l = (ZL0*YL1-YL0*ZL1)/ZL0**2
		dens_v = (ZV0*YV1-YV0*ZV1)/ZV0**2
	else: print 'not a valid derivative'; return
	
	# calculate co2 properties
	if not derivative: k = 'viscosity'
	elif derivative in ['P','pressure']: k = 'dvdp'
	elif derivative in ['T','temperature']: k = 'dvdt'
	
	if not co2Vars: return (dens_l,dens_v,np.array([]))
	
	if not P.shape: P = np.array([P])
	if not T.shape: T = np.array([T])
	if P.size == 1 and not T.size == 1: P = P*np.ones((1,len(T)))[0]
	elif T.size == 1 and not P.size == 1: T = T*np.ones((1,len(P)))[0]
	# calculate bounding values of P
	dens_c = []
	for Pi, Ti in zip(P,T):
		if Pi<=Pval[0]: p0 = Pval[0]; p1 = Pval[0]
		elif Pi>=Pval[-1]: p0 = Pval[-1]; p1 = Pval[-1]
		else:
			p0 = Pval[0]
			for p in Pval[1:]:
				if Pi<=p: 
					p1 = p; break
				else: 
					p0 = p
		# calculate bounding values of T
		if Ti<=Tval[0]: t0 = Tval[0]; t1 = Tval[0]
		elif Ti>=Tval[-1]: t0 = Tval[-1]; t1 = Tval[-1]
		else:
			t0 = Tval[0]
			for t in Tval[1:]:
				if Ti<=t: 
					t1 = t; break
				else: 
					t0 = t
		# calculate four indices
		dt0 = abs(Ti-t0); dt1 = abs(Ti-t1); dp0 = abs(Pi-p0); dp1 = abs(Pi-p1)
		t0 = Tdict[t0]; t1 = Tdict[t1]; p0 = Pdict[p0]; p1 = Pdict[p1]
		i1 = p0*tn+t0
		i2 = p0*tn+t1
		i3 = p1*tn+t0
		i4 = p1*tn+t1
		# locate value in array
		v1 = arrays[k][i1]
		v2 = arrays[k][i2]
		v3 = arrays[k][i3]
		v4 = arrays[k][i4]
		dens_c.append((p0*(t1*v3+t0*v4)/(t1+t0) + p1*(t1*v1+t0*v2)/(t1+t0))/(p0+p1))
	
	return (dens_l,dens_v,np.array(dens_c)*1e-6)
def sat(T):
	"""Return saturation pressure and first derivative for given temperature.
	
	:param T: Temperature (degC)
	:type T: fl64
	:returns: Two element tuple containing (saturation pressure, derivative).
	
	"""
	Y0 = (YSP[0]+
	      YSP[1]*T+
	      YSP[2]*T**2+
	      YSP[3]*T**3+
	      YSP[4]*T**4)
	Z0 = (ZSP[0]+
	      ZSP[1]*T+
	      ZSP[2]*T**2+
	      ZSP[3]*T**3+
	      ZSP[4]*T**4)
	Y1 = (YSP[1]+
	      YSP[2]*T*2+
	      YSP[3]*T**2*3+
	      YSP[4]*T**3*4)
	Z1 = (ZSP[1]+
	      ZSP[2]*T*2+
	      ZSP[3]*T**2*3+
	      ZSP[4]*T**3*4)
	satP = Y0/Z0
	dsatPdT = (Z0*Y1-Y0*Z1)/Z0**2
	return (satP, dsatPdT)
def tsat(P):
	"""Return saturation temperature and first derivative for given pressure.
	
	:param P: Pressure (degC)
	:type P: fl64
	:returns: Two element tuple containing (saturation temperature, derivative).
	
	"""
	Y0 = (YST[0]+
	      YST[1]*P+
	      YST[2]*P**2+
	      YST[3]*P**3+
	      YST[4]*P**4)
	Z0 = (ZST[0]+
	      ZST[1]*P+
	      ZST[2]*P**2+
	      ZST[3]*P**3+
	      ZST[4]*P**4)
	Y1 = (YST[1]+
	      YST[2]*P*2+
	      YST[3]*P**2*3+
	      YST[4]*P**3*4)
	Z1 = (ZST[1]+
	      ZST[2]*P*2+
	      ZST[3]*P**2*3+
	      ZST[4]*P**3*4)
	satT = Y0/Z0
	dsatTdP = (Z0*Y1-Y0*Z1)/Z0**2
	return (satT, dsatTdP)
def fluid_column(z,Tgrad,Tsurf,Psurf,iterations = 3):
	'''Calculate thermodynamic properties of a column of fluid.
	
	:param z: Vector of depths at which to return properties. If z does not begin at 0, this will be prepended.
	:type z: ndarray
	:param Tgrad: Temperature gradient in the column (degC / m).
	:type Tgrad: fl64
	:param Tsurf: Surface temperature (degC).
	:type Tsurf: fl64
	:param Psurf: Surface pressure (MPa).
	:type Psurf:
	:param iterations: Number of times to recalculate column pressure based on updated density.
	:type iterations: int
	:returns: Three element tuple containing (liquid, vapor, CO2) properties. Each contains a three column array corresponding to pressure, temperature, density, enthalpy and viscosity of the fluid.
	
	'''
	z = abs(np.array(z))
	
	if z[-1] < z[0]: z = np.flipud(z)
	if z[0] != 0: z = np.array([0,]+list(z))
	
	if isinstance(Tgrad,str): 		# interpret Tgrad as a down well temperature profile
		if not os.path.isfile(Tgrad): print 'ERROR: cannot find temperature gradient file \''+Tgrad+'\'.'; return
		
		tempfile = open(Tgrad,'r')
		ln = tempfile.readline()
		tempfile.close()
		commaFlag = False; spaceFlag = False
		if len(ln.split(',')) > 1: commaFlag = True
		elif len(ln.split()) > 1: spaceFlag = True
		if not commaFlag and not spaceFlag: print 'ERROR: incorrect formatting for \''+Tgrad+'\'. Expect first column depth (m) and second column temperature (degC), either comma or space separated.'; return
		if commaFlag: tempdat = np.loadtxt(Tgrad,delimiter=',')
		else: tempdat = np.loadtxt(Tgrad)
		zt = tempdat[:,0]; tt = tempdat[:,1]
		
		T = np.interp(z,zt,tt)
		
	else:
		Tgrad = abs(Tgrad)
		T = Tsurf + Tgrad*z
	Pgrad = 800*9.81/1e6
	Phgrad = 1000*9.81/1e6
	if co2Vars:
		Pco2 = Psurf + Pgrad*z
	Ph = Psurf + Phgrad*z
	
	for i in range(iterations):
		if co2Vars:
			rho = dens(Pco2,T)[2]
			Pco2=np.array([(abs(np.trapz(rho[:i+1],z[:i+1]))*9.81/1e6)+Pco2[0] for i in range(len(rho))])
		rho = dens(Ph,T)[0]
		Ph=np.array([(abs(np.trapz(rho[:i+1],z[:i+1]))*9.81/1e6)+Ph[0] for i in range(len(rho))])
	
	if co2Vars:
		rho = dens(Pco2,T)
		H = enth(Pco2,T)
		mu = visc(Pco2,T)
	
	rhoh = dens(Ph,T)
	Hh = enth(Ph,T)
	muh = visc(Ph,T)
	
	if co2Vars:
		return (np.array([Ph,T,rhoh[0],Hh[0],muh[0]]).T,np.array([Ph,T,rhoh[1],Hh[1],muh[1]]).T,np.array([Pco2,T,rho[2],H[2],mu[2]]).T)
	else:
		return (np.array([Ph,T,rhoh[0],Hh[0],muh[0]]).T,np.array([Ph,T,rhoh[1],Hh[1],muh[1]]).T,np.array([]))
	