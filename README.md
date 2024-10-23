## FEHM: Finite Element Heat and Mass Transfer Code ##
**LANL Software: LA-CC-2012-083  No. C13022**
**LANL Documents: LA-UR-12-24493**


#### [FEHM Homepage](https://fehm.lanl.gov) • [FEHM Documentation](http://lanl.github.io/FEHM/) • [Latest Release V3.4](https://github.com/lanl/FEHM/releases/tag/v3.4.0)


The numerical background of the FEHM computer code can be traced to the early 1970s when it was used to simulate geothermal and hot dry rock reservoirs. The primary use over a number of years was to assist in the understanding of flow fields and mass transport in the saturated and unsaturated zones below the potential Yucca Mountain repository. Today FEHM is used to simulate groundwater and contaminant flow and transport in deep and shallow, fractured and un-fractured porous media throughout the US DOE complex. FEHM has proved to be a valuable asset on a variety of projects of national interest including Environmental Remediation of the Nevada Test Site, the LANL Groundwater Protection Program, geologic CO2 sequestration, Enhanced Geothermal Energy (EGS) programs, Oil and Gas production, Nuclear Waste Isolation, and Arctic Permafrost. Subsurface physics has ranged from single fluid/single phase fluid flow when simulating basin scale groundwater aquifers to complex multifluid/ multi-phase fluid flow that includes phase change with boiling and condensing in applications such as unsaturated zone surrounding nuclear waste storage facility or leakage of CO2/brine through faults or wellbores. The numerical method used in FEHM is the control volume method (CV) for fluid flow and heat transfer equations which allows FEHM to exactly enforce energy/mass conservation; while an option is available to use the finite element (FE) method for displacement equations to obtain more accurate stress calculations. In addition to these standard methods, an option to use FE for flow is available, as well as a simple Finite Difference scheme.


## License ##

FEHM is distributed as as open-source software under a BSD 3-Clause License. See [Copyright License](LICENSE.md)

## Collaborators and Development ##

FEHM Developers follow instructions at [Developers](develop.md)

External Collaborators must sign a Contribution Agreement. [Contribution Agreement for External Collaborators](CONTRIBUTING.md)

## Build FEHM ##

FEHM V3.6.2 compiles and passes fehmpytests and VV Test Suites using Linux GNU Fortran (GCC) 13.2.0 and Windows Intel® Fortran Compiler Classic 2021.11.0 [Intel(R) 64]. 

> [!NOTE]
> Intel Classic 2021.12.0 has less tolerance for old fortran code. It passes all fehmpytests but the VV Test Suite has incomplete runs for fracture_aperture, reverse_tracking, colloid_stream, stress_3D, stress_3Dbeam, and uz_test. The remaining 78 tests are successful. This issue will be addressed in future versions.* 

Download the Https repo to your current directory by running:

```bash
git clone https://github.com/lanl/FEHM.git
cd FEHM
```
Or Download the ZIP file [HERE](https://github.com/lanl/FEHM/archive/refs/heads/master.zip).

Next, you will need to create the executable. This is different depending on your operating system. For Linux/Mac, build FEHM using the FORTRAN code in `FEHM/src`. For Windows, Build FEHM using Visual Studio with the Fortran extension in the `FEHM` directory.

*	**Linux/Mac**

	This uses a MakeFile to create an executable.

	```bash
	cd src; make
	```
	This will create an executable named `xfehm` in the `/src` directory. *Note. Makefile will overwrite dated.f using dated.template.* 


*	**Windows**

	This requires that you have Visual Studio installed along with Fortran extensions enabled.

	1. Using your File Explorer, navigate to the location where you downloaded `FEHM`.
	2. Double click on the `.vfproj` file. You may need to manually select Visual Studio as the correct program to run it. This will open Visual Studio with the file listed under the Solution Explorer.
	3. Make sure the "Release" and "x64" are selected in the drop downs on the menu bar and that (IFORT) is shown on the Solution Explorer, then click "Start". This will create an executable in the `FEHM/x64/release` directory.


> [!NOTE]
> Windows Users - Using IFORT may give a deprecation warning in Visual Studio. It is safe to ignore this message. FEHM runs best with IFORT and using IFX may give unpredicatble results.

### Test FEHM ###

FEHM provides a set of python driven tests to verify the FEHM installation. These tests require Python be installed.

```bash
cd FEHM/fehmpytests/
python fehmpytests.py <FEHM executable path>
# Example Linux/Mac:
# python fehmpytests.py ../src/xfehm

# Example Windows
# python fehmpytests.py ../x64/release/FEHM3.6_VER2.exe
```

This will run a series of tests that will take several minutes to run and will look similar to this:

```
-----------------------------------------------------------------------------------
Ran 30 tests in 120.000s

OK

```

Detailed testing information is available at [https://lanl.github.io/FEHM/fehmpytest-info/testing](https://lanl.github.io/FEHM/fehmpytest-info/testing)

See **fehmpytests** documentation:
[https://lanl.github.io/FEHM/getting-started/fehmpytests](https://lanl.github.io/FEHM/getting-started/fehmpytests)



## V&V Verification Test Suite ##

This large test suite was created under the YMP QA program and has been the main source for testing and development over the years. The disadvantage of this test suite is its large size, complicated perl scripts, fortran compare executables, and platform dependency. The Windows version of the test suite is most up to date, runs over 80 test problems, and takes over 1.5 hours to run. Most test cases in FEHMPYTESTS are constructed from problems in the Windows V&V test suite.


The V&V FEHM Test Suite documentation: https://www.lanl.gov/orgs/ees/fehm/docs/FEHM_VERIFICATION_V3.3.0.pdf

The V&V FEHM Test Suite can be downloaded from Assets under Releases.
https://github.com/lanl/FEHM/releases
- VERIFICATION_V3.3.0lnx.tar.gz
- VERIFICATION_V3.3.0mac.tar.gz
- VERIFICATION_V3.3.0win.zip 


## FEHM Applications ##

**dfnWorks** is a python driven parallelized computational suite to generate three-dimensional discrete fracture networks (DFN) and simulate flow and transport. It uses LagriT, PETSC, PFLOTRAN, Python, and FEHM. See dfnWorks [Documentation](https://dfnworks.lanl.gov) and [github](https://github.com/lanl/dfnWorks)

**PYFEHM** is a set of classes and methods to enable use of FEHM and auxiliary tasks within the Python scripting environment. This is maintained outside of LANL is is no longer supported. See PYFEHM scripts on GitHub: https://github.com/lanl/PyFEHM


## FEHM Release Versions ##


See Versions and Notes under the Releases tab this repository.

The Most recent distributed release is FEHM V3.4.0 (September 2019) which is the version cloned for this repository. The FEHM software is a continuation of QA work performed for the Yucca Mountain Project (YMP) under Software Configuration Control Request (SCCR) (Software Tracking Numbers STN: 10086-2.21-00 August 2003, V2.22, STN 10086-2.22-01, V2.23, STN 10086-2.23-00, V2.24-01, STN 10086-2.24-01, and V2.25, STN 10086-2.25-00). 
The QA for these codes started under YMP QA and continue under under LANL EES-16 Software QA Policy and Proceedures as outlined in: "EES-16-13-003.SoftwareProcedure.pdf" 

Before distribution of FEHM software, tests are executed and verified as acceptable on LANL computers with operating systems Linux, Mac OSX, and WINDOWS. The overall validation effort for the FEHM software consists of a suite of directories and scripts that test the model whenever possible, against known analytical solutions of the same problem. The test suite was developed under YMP QA for FEHM RD.10086-RD-2.21-00 and is available for download.

