Running FEHM-mars
=================

Example problem located in:
    <REPO>/example/

Full model description can be found in:
    <JGR Citation>

Please cite if you use any of these scripts in your work. 


Building FEHM-mars
------------------

Build the FEHM-mars executable by following the instructions:

    https://github.com/lanl/FEHM/tree/mars


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


