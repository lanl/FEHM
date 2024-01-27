Introduction
============

Fehmpytests is a new test suite for FEHM. Its goal is to enable FEHM developers 
to easily test new code and add new tests for existing or future functionality.
To meet these goals, Fehmpytests uses the Python "unittest" module and a general 
test method that can be called for each new test case. Currently, there are
twenty tests that Fehmpytests performs through a command line interface.  
Future plans are to integrate it into the FEHM build process and provide more 
modularity by improving the devloper interface. 

Fehmpytests now contains a simplified script from `pyfehm <http://www.github.com/lanl/pyfehm/>`_ (fpost.py) that requires only python>=2.7 without any external packages.
