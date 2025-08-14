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


Use fehmpytests to test changes to FEHM to ensure correct simulation. Fehmpytests can be run in four different modes, default, admin, developer (currently in development), and solo which each run a set of tests. Currently, default and admin mode run the same set of tests and solo mode runs a single test. Developer mode is not implemented yet but will run a subset of admin tests. Please ensure that you are using Python >3.5.

## Testing in Default Mode

To test the default suite:

1. Navigate to the folder *fehmpytests* in a terminal.
2. Type the following command into the terminal:

   ``python fehmpytests.py <fehm-path>``
       
   where ```<fehm-path>``` is the path to the FEHM executable.

## Testing in Admin Mode

To test the admin suite:

1. Navigate to the folder *fehmpytests* in a terminal.
2. Type the following command into the terminal:

   ``python fehmpytests.py --admin <fehm-path>``
   
   where ```<fehm-path>``` is the path to the FEHM executable.
   
## Testing in Developer Mode

(In Development) To test the developer suite:

1. Navigate to the folder *fehmpytests* in a terminal.
2. Type the following command into the terminal:

   ``python fehmpytests.py --dev <fehm-path>``
   
   where ```<fehm-path>``` is the path to the FEHM executable.
                
## Testing in Solo Mode

The process for testing a single test case in solo mode is similar to testing 
a suite in the other modes. There is an additional command line argument needed.
 
To test a singe test-case:

1. Navigate to the folder *fehmpytests* in a terminal.
2. Type the following command into the terminal:

   ``python fehmpytests.py <fehm-path> <test-case>``
     
   where ```<fehm-path>``` is the location of the FEHM executable and <test-case> 
   is the name of the test-case method.
   
> ##### WARNING
>
> Developers must run in default, admin, or developer mode before commiting new code. 
{: .block-warning}


## Creating an Error Log

An error log .txt file can be created to show details about an error and where 
it occurred. To generate an error log, add the switch *log* after 
**fehmpytests.py** and before ```<fehm-path>```. Here is an example:

``python fehmpytests.py --admin --log <fehm-path>``


