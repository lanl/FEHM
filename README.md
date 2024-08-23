## FEHM: Finite Element Heat and Mass Transfer Code ##
**LANL Software: LA-CC-2012-083  No. C13022**
**LANL Documents: LA-UR-12-24493**

[![PyPI](https://img.shields.io/pypi/l/Django.svg)](https://github.com/lanl/FEHM/blob/master/LICENSE.md)

[![readthedocs](https://img.shields.io/static/v1?label=Documentation&message=Read%20online&color=blue&style=for-the-badge&logo=read-the-docs)](http://lanl.github.io/FEHM/) <br/>
[![readthedocs](https://img.shields.io/static/v1?label=FEHM%20Home&message=Read%20online&color=blue&style=for-the-badge&logo=read-the-docs)](https://fehm.lanl.gov/) <br/>

The numerical background of the FEHM computer code can be traced to the early 1970s when it was used to simulate geothermal and hot dry rock reservoirs. The primary use over a number of years was to assist in the understanding of flow fields and mass transport in the saturated and unsaturated zones below the potential Yucca Mountain repository. Today FEHM is used to simulate groundwater and contaminant flow and transport in deep and shallow, fractured and un-fractured porous media throughout the US DOE complex. FEHM has proved to be a valuable asset on a variety of projects of national interest including Environmental Remediation of the Nevada Test Site, the LANL Groundwater Protection Program, geologic CO2 sequestration, Enhanced Geothermal Energy (EGS) programs, Oil and Gas production, Nuclear Waste Isolation, and Arctic Permafrost. Subsurface physics has ranged from single fluid/single phase fluid flow when simulating basin scale groundwater aquifers to complex multifluid/ multi-phase fluid flow that includes phase change with boiling and condensing in applications such as unsaturated zone surrounding nuclear waste storage facility or leakage of CO2/brine through faults or wellbores. The numerical method used in FEHM is the control volume method (CV) for fluid flow and heat transfer equations which allows FEHM to exactly enforce energy/mass conservation; while an option is available to use the finite element (FE) method for displacement equations to obtain more accurate stress calculations. In addition to these standard methods, an option to use FE for flow is available, as well as a simple Finite Difference scheme.

---

## Quick Start Guide

### Downloading & Building FEHM ###

Download the Https repo to your current directory by running:

```bash
git clone https://github.com/lanl/FEHM.git
cd FEHM/
```

Next, you will need to create the executable. This is different depending on your operating system.

*	**Linux/Mac**

	This uses a MakeFile to create an executable.

	```bash
	cd src/
	```
	This will create an executable named `xfehm` in the `/src` directory.


*	**Windows**

	This requires that you have Visual Studio installed along with Fortran extensions enabled.

	1. Using your File Explorer, navigate to the location where you downloaded FEHM.
	2. Navigate to `FEHM/src/PC/`
	3. Double click on the `.vfproj` file. You may need to manually select Visual Studio as the correct program to run it. This will open Visual Studio with the file listed under the Solution Explorer.
	4. Make sure the "Release" and "x64" are selected in the drop downs on the menu bar and then click "Start". This will create an executable in the `/src/x64/release` directory.

### Testing ###

FEHM provides a small subset of tests to verify the FEHM installation. These tests require Python be installed.

```bash
cd FEHM/fehmpytests/
python fehmpytests.py <FEHM executable path>
# Example Linux/Mac:
# python fehmpytests.py ../src/xfehm

# Example Windows
# python fehmpytests.py ../src/PC/x64/release/FEHM3.4_VER2.exe
```

This will run a series of tests and will look similar to this:

```
-----------------------------------------------------------------------------------
Ran 29 tests in 120.000s

OK

```

Detailed testing information is available at [https://lanl.github.io/FEHM/fehmpytest-info/testing](https://lanl.github.io/FEHM/fehmpytest-info/testing)

See **fehmpytests** documentation:
[https://lanl.github.io/FEHM/getting-started/fehmpytests](https://lanl.github.io/FEHM/getting-started/fehmpytests)


## Advanced for Developers ##

External Collaborators must sign a Contribution Agreement. [Contribution Agreement for External Collaborators](CONTRIBUTING.md)

The following are reminders for FEHM code developers using this repository.

A Git workflow follows these basic steps:

* Make changes to files
* Test changes by adding them to fehmpytests and running to verify compatibility
* Add the files (‘stage’ files)
* ‘Commit’ the staged files
* Push the commit (containing all modified files) to the central repo

To get started:

1. Clone repo and create executable as above.
 
2. Let’s say you’ve done some editing. The next step is to add your new test to fehmpytests (if not already there) and run the test suite to confirm that the code works correctly.

	To run Tests:

	* Add test files under `/fehmpytests`
	* Add Test Case to `fehmpytests.py`
	* Run test as a solo test

	Detailed information on adding tests is available at [https://lanl.github.io/FEHM/fehmpytest-info/newtest](https://lanl.github.io/FEHM/fehmpytest-info/newtest)

	```bash
	python fehmpytests.py ../src/xfehm testcase
	```
	Detailed testing information is available at [https://lanl.github.io/FEHM/fehmpytest-info/testing](https://lanl.github.io/FEHM/fehmpytest-info/testing)

3. If the test is sucessful and you’re ready to push your changes to the FEHM repository.
Run the command

	```
	git add file1 file2 ... fileN
	```
	 
	to add any files you have changed. You can also just run `git add .` if you want to add every changed file.
 
4. Now, run the command
 
	``` 
	git status
	```
	 
	This gives an overview of all tracked and untracked files.
	A tracked file is one that Git considers as part of the repo.
	Untracked files are everything else – think of *.o files, or some test data output generated by an FEHM run.
	 
	Tracked files can be:
	* Unmodified (you haven’t made any changes to it, relative to the last commit)
	* Modified (you have edited the file since the last commit)
	* Staged (the file has been added and is ready to be committed and then pushed)
	 
	Untracked files become tracked by using
	```
	git add filename
	```
 
5. After verifying (with `git status`) that all the files you want to be pushed are properly staged, commit them using

	```
	git commit -m "My first Git commit!"
	```
	 
	Then, push the files onto the GitHub repo with

	```
	git push origin master
	```
 
6. If someone else has made edits, you will need to pull their changes to your local FEHM clone before you can push.
	 
	```
	git pull origin master
	git push origin master
	```

## PYFEHM ##

PYFEHM is a set of classes and methods to enable use of FEHM and auxiliary tasks within the Python scripting environment.

PYFEHM scripts on GitHub: https://github.com/lanl/PyFEHM

Documentation: https://lanl.github.io/PyFEHM/


## FEHM Release Versions ##

See Versions and Notes under the Releases tab this repository.

The Most recent distributed release is FEHM V3.4.0 (November 2019) which is the version cloned for this repository. The FEHM software is a continuation of QA work performed for the Yucca Mountain Project (YMP) under Software Configuration Control Request (SCCR) (Software Tracking Numbers STN: 10086-2.21-00 August 2003, V2.22, STN 10086-2.22-01, V2.23, STN 10086-2.23-00, V2.24-01, STN 10086-2.24-01, and V2.25, STN 10086-2.25-00). 
The QA for these codes started under YMP QA and continue under under LANL EES-16 Software QA Policy and Proceedures as outlined in: "EES-16-13-003.SoftwareProcedure.pdf" 

Before distribution of FEHM software, tests are executed and verified as acceptable on LANL computers with operating systems Linux, Mac OSX, and WINDOWS. The overall validation effort for the FEHM software consists of a suite of directories and scripts that test the model whenever possible, against known analytical solutions of the same problem. The test suite was developed under YMP QA for FEHM RD.10086-RD-2.21-00 and is available for download.

