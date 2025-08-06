---
title : FehmPyTests Suite
layout : page_getting-started
permalink: /fehmpytest-info/testing
hero_height: is-hidden
---

# FehmPyTests Suite

### Topics

* [Installation](install.md)
* Testing FEHM
* [Creating New Test-Cases](newtest.md)
* [Test-Case Description](testdesc.md)

---



FEHM provides a set of python driven tests to verify the FEHM installation. These tests require Python3 be installed.


## Running the Tests

To run and verify the default set of tests:

Navigate to the folder *fehmpytests* in a terminal.
Type the following command into the terminal:

   ``python fehmpytests.py <fehm-path>``
       
   where ```<fehm-path>``` is the path to the FEHM executable, for instance ```../src/xfehm```


To run a single test case:

   ``python fehmpytests.py <fehm-path> <test-case>``
     
   where ```<fehm-path>``` is the location of the FEHM executable and <test-case> is the name of the test-case to run.
   

## Options

Options can be a -letter or --keyword as shown below. They should be used on the command line before the executable path and name.

    ``python fehmpytests.py [-h] [-a | -d | -f | -s] [-l] [-v {0,1,2,3}] [--clean] exe [testcase]``
    

The following options are available with fehmpytests:


``-h  --help``               
    Show this help message and exit.
    
``-a  --admin``              
    Run the entire test-suite, this is the default mode.
    
``-d  --dev``                
    Run a portion of the test-suite (set in fehmpytests.py by developers).

``-f  --full``              
    Run the entire test-suite including those under development with possible failures.
    
``-l  --log``                
    Create a fail statistics file 'fail_log.txt'.
    
``-v {0,1,2,3}  --verbose {0,1,2,3}``    
    Verbosity level: 0 = No output, 1 = Minimal Output, 2 = Detailed Output, 3 = All Output, for developers
    
``--clean``                  
   Clean up, remove fehm output files and exit.
   

## Example fehmpytests session

```bash
cd FEHM/fehmpytests/
python fehmpytests.py ../src/xfehm

```

This will run a series of tests that will take several minutes to run and will look similar to this:

```
avdonin (__main__.fehmTest.avdonin) ... ok
baro_vel (__main__.fehmTest.baro_vel) ... ok
bodyforce (__main__.fehmTest.bodyforce) ... ok
...
transport3d_validation (__main__.fehmTest.transport3d_validation) ... ok
vapor_extraction (__main__.fehmTest.vapor_extraction) ... ok
wvtest (__main__.fehmTest.wvtest) ... ok
rad_decay (__main__.fehmTest.rad_decay) ... ok
-----------------------------------------------------------------------------------
Ran 30 tests in 120.000s

OK

```


## V&V Verification Test Suite ##

Th large V&V verification test suite was created under the YMP QA program and has been the main source for testing and development over many years. The disadvantage of this test suite is its large size, complicated perl scripts, old fortran compare executables, and platform dependency. The Windows version of the test suite is used to verify modern FEHM simulations. The Windows Test Suite runs over 80 test problems, and can takes over an hour to run. The Linux and Mac versions of the Test Suite are no longer supported. They have perl scripts that no longer work and old input files. Though no longer used for code verification, the V&V problems can be run individually to explore examples and aid development.

Most thee test cases in fehmpytests are constructed from problems in the Windows V&V test suite.

The V&V FEHM Test Suite documentation: https://www.lanl.gov/orgs/ees/fehm/docs/FEHM_VERIFICATION_V3.3.0.pdf

The V&V FEHM Test Suite can be downloaded from Assets under Releases.
https://github.com/lanl/FEHM/releases
- VERIFICATION_V3.3.0lnx.tar.gz
- VERIFICATION_V3.3.0mac.tar.gz
- VERIFICATION_V3.3.0win.zip
- 
   
> ##### WARNING
>
> Developers must run in default or admin modes before commiting new code.
{: .block-warning}





