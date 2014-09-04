Introduction
============

Fehmpytests is a new test suite for FEHH. Its goal is to enable FEHM developers 
to easily test new code and add new tests for existing or future functionality.
To meet these goals, Fehmpytests uses the Python unittest module and a general 
test method that can be called for each new test case. Currently, there are
thirteen tests that Fehmpytests performs through a command line interface.  
Future plans are to integrate it into the FEHM build process and provide more 
modularity by improving the devloper interface. 

Fehmpytests uses python 2.7 and expects the scipy module to be available.
To verify your python version enter the following in a terminal:

   ``python --version``

If scipy is not available, contact your system adminstrator.

