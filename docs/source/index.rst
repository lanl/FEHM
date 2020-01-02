===================================================
Software Users Manual (UM) for the FEHM Application
===================================================

This online version of the FEHM User's Manual documents the use of the FEHM application.
It is the most current description of FEHM and its use and is continually updated as code is developed.

The following modeling and functional capabilities lists are additional to the
FEHM UMV3.1 and are provided as an overview and command context for FEHM
users.

Overview of FEHM Capabilities
-----------------------------

* 3-dimensional complex geometries with unstructured grids
* saturated and unsaturated media  
* simulation of production from gas hydrate reservoirs  
* simulation of geothermal reservoirs
* non-isothermal, multi-phase flow of gas, water, oil  
* non-isothermal, multi-phase flow of air, water  
* non-isothermal, multi-phase flow of CO2, water  
* multiple chemically reactive and sorbing tracers  
* preconditioned conjugate gradient solution of coupled linear equations
* fully implicit, fully coupled Newton Raphson solution of nonlinear equations  
* double porosity and double porosity/double permeability capabilities  
* control volume (CV) and finite element method (FE) methods
* coupled geomechanics (THM) problems (fluid flow and heat transfer coupled with stress/deformation) including non-linear elastic and plastic deformation, nonlinear functional dependence of rock properties (e.g.permeability, porosity, Young's modulus) on pressure, temperature and damage/stress

FEHM Documentation
------------------

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   SystemInterface.rst
   DataFiles.rst
   DefinitionsAndAcronyms.rst
   ExamplesAndSamples.rst

.. toctree::
   :maxdepth: 1
   :caption: Program Specifications

   Macros.rst
   InputData.rst
   Output.rst
   ProgramConsiderations.rst

.. toctree::
   :maxdepth: 1
   :caption: References & Support

   References.rst
   Support.rst
   ReleaseNotes/release_notes.rst

