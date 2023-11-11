Running FEHM-mars
=================

Example problem located in:
    <REPO>/example/

Full model description can be found in:
    <JGR Citation>

Please cite if you use any of these scripts in your work. 


Dependencies
------------

Install the MATK (Model Analysis ToolKit) for Python by following the instructions on GitHub:
    https://github.com/dharp/matk


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


1-D Atmospheric Mixing Simulations
----------------------------------

