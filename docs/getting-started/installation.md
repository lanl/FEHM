---
title : Installation Instructions
layout : page_getting-started
permalink: /getting-started/installation
hero_height: is-hidden
---

# Installation Instructions

## Build FEHM ##

FEHM V3.6.2 compiles and passes fehmpytests and VV Test Suites using Linux GNU Fortran (GCC) 13.2.0 and Windows Intel® Fortran Compiler Classic 2021.11.0 [Intel(R) 64]. 

> **_NOTE_**
> 
>Intel Classic 2021.12.0 has less tolerance for old fortran code. It passes all fehmpytests but the VV Test Suite has incomplete runs for fracture_aperture, reverse_tracking, colloid_stream, stress_3D, stress_3Dbeam, and uz_test. The remaining 78 tests are successful. This issue will be addressed in future versions.* 

Download the Https repo to your current directory by running:

```bash
git clone https://github.com/lanl/FEHM.git
cd FEHM
```
Or Download the ZIP file [HERE](https://github.com/lanl/FEHM/archive/refs/heads/master.zip).

Next, you will need to create the executable. This is different depending on your operating system. For Linux/Mac, build FEHM using the FORTRAN code in `FEHM/src`. For Windows, Build FEHM using Visual Studio with the Fortran extension in the `FEHM` directory.

* **Linux/Mac**

  This uses a MakeFile to create an executable.

  ```bash
  cd src; make
  ```
  This will create an executable named `xfehm` in the `/src` directory. *Note. Makefile will overwrite dated.f using dated.template.* 


* **Windows**

  This requires that you have Visual Studio installed along with Fortran extensions enabled.

  1. Using your File Explorer, navigate to the location where you downloaded `FEHM`.
  2. Double click on the `.vfproj` file. You may need to manually select Visual Studio as the correct program to run it. This will open Visual Studio with the file listed under the Solution Explorer.
  3. Make sure the "Release" and "x64" are selected in the drop downs on the menu bar and that (IFORT) is shown on the Solution Explorer, then click "Start". This will create an executable in the `FEHM/x64/release` directory.


> **_NOTE_**
>
> Windows Users - Using IFORT may give a deprecation warning in Visual Studio. It is safe to ignore this message. FEHM runs best with IFORT and using IFX may give unpredictable results.

### Test FEHM ###

FEHM provides a set of python driven tests to verify the FEHM installation. These tests require Python be installed. More information on the simplified test suite can be found [**Here**](fehmpytests)

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


The V&V FEHM Test Suite documentation: [https://www.lanl.gov/orgs/ees/fehm/docs/FEHM_VERIFICATION_V3.3.0.pdf](https://www.lanl.gov/orgs/ees/fehm/docs/FEHM_VERIFICATION_V3.3.0.pdf)

The V&V FEHM Test Suite can be downloaded from Assets under Releases.
[https://github.com/lanl/FEHM/releases](https://github.com/lanl/FEHM/releases)
- VERIFICATION_V3.3.0lnx.tar.gz
- VERIFICATION_V3.3.0mac.tar.gz
- VERIFICATION_V3.3.0win.zip

