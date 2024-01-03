FEHM-mars
=========


Citation
--------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10455952.svg)](https://doi.org/10.5281/zenodo.10455952)


Full model description can be found in:
    <JGR Citation>

Please cite if you use the software or any of these scripts in your work. 


LICENSE
-------

This software is open source software available under the BSD-3 license.

Copyright © 2018. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so. This is open source software; you can redistribute it and/or modify it under the terms of the BSD-3 License. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



Clone and Build FEHM-mars
-------------------------

First, clone the FEHM-mars branch using the following command:

    git clone -b mars git@github.com:lanl/FEHM.git

Then build the FEHM-mars executable by following the instructions:

    https://github.com/lanl/FEHM/tree/mars

FEHM-mars example problem located in:
    <REPO>/example/

Dependencies
------------

Install the MATK (Model Analysis ToolKit; ``http://dharp.github.io/matk/``) for Python by following the instructions on GitHub:

    http://dharp.github.io/matk/

    
Heat Flow Simulation
--------------------

Heat flow simulation must be performed once in order to generate the time
series of subsurface temperatures. These are eventually used to calculate the
temperature-dependent adsorption coefficients. 

1. Navigate to ``1d_heat_runs/heat_200m_synth/`` 
2. Run command:
    python run.py
3. If you have not set the ``RUNDIR`` or ``FEHM_MARS_EXE`` environment variables, the ``heat_master_model.py`` script will prompt you to set them (with instructions). 
4. Once simulation is finished, proceed to next step. 


Tracer Flow & Transport Simulations
-----------------------------------

1. Navigate to ``sim_dir``.
2. Run command:
    python run.py
3. The ``run.py`` script contains variables specific to the type of simulation you are running. 
    
    - Simulations with different parameters can be set up in a different directory and modifying the ``run.py`` file. 
    - Ensure that the ``master_model.py`` script is still one directory above the simulation directory 
4. If you have not set the ``RUNDIR`` or ``FEHM_MARS_EXE`` environment variables, the ``master_model.py`` script will prompt you to set them (with instructions). 
5. Initialization, spinup, and transport simulations will then run sequentially, as described in the manuscript. 
    
    - Please note that the transport simulations can take several days/weeks to complete due to computational time required for solutions to reach a cyclic steady-state condition.

6. Plots will out be placed in ``output`` sub-directory. 


1-D Atmospheric Mixing Simulations
----------------------------------

1. Navigate to ``sim_dir`` and Create a symbolic link to the p``pbl_diffusion.py`` script in the ``example/`` directory:

    cd sim_dir
    ln -s ../pbl_diffusion.py .

2. If you have not already run the script and performed the parameter estimation algorithm, you can run the whole script with the following command:

    python pbl_diffusion.py

3. Plots are generated and placed in ``pbl_output`` sub-directory.


