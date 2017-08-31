# PyTOUGH setup script
from distutils.core import setup

setup(name='PyFEHM',
	version='1.0.4',
	description='Python scripting library for FEHM simulations',
	author='David Dempsey',
	author_email='d.dempsey@lanl.gov',
	url='pyfehm.lanl.gov',
	license='LGPL',
	py_modules=['ftool','fgrid','fdata','fpost','fvars','fdflt','fhelp','ftemp'],
	packages=['pyvtk'],
	scripts = ['scripts/diagnose.py','scripts/fehm_paraview.py']
	)
