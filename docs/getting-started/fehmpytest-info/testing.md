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


Use fehmpytests to test changes to FEHM to ensure correct simulations. Please ensure that you are using Python3.


## Running the Tests

To test the default set of tests:

Navigate to the folder *fehmpytests* in a terminal.
Type the following command into the terminal:

   ``python fehmpytests.py <fehm-path>``
       
   where ```<fehm-path>``` is the path to the FEHM executable, for instance ```../src/xfehm```


To run a single test case:

   ``python fehmpytests.py <fehm-path> <test-case>``
     
   where ```<fehm-path>``` is the location of the FEHM executable and <test-case> is the name of the test-case to run.
   

## Options

In addition to mode settings, the following options are available. Option flags should occur on command line before the executable path and name.

    ``python fehmpytests.py [-h] [-a | -d | -f | -s] [-l] [-v {0,1,2,3}] [--clean] exe [testcase]``


**-h  --help**               
    Show this help message and exit.
    
**-a  --admin**              
    Run the entire test-suite, this is the default mode.
    
**-d  --dev**                
    Run a portion of the test-suite (set in fehmpytests.py by developers).

**-f  --full**              
    Run the entire test-suite including those under development with possible failures.
    
**-l  --log**                
    Create a fail statistics file 'fail_log.txt'.
    
**-v {0,1,2,3}  --verbose {0,1,2,3}**     
    Verbosity level: 0 = No output, 1 = Minimal Output, 2 = Detailed Output, 3 = All Output, for developers
    
**--clean**                  
   Clean up, remove fehm output files and exit.
   



   
> ##### WARNING
>
> Developers must run in default or admin modes before commiting new code.
{: .block-warning}





