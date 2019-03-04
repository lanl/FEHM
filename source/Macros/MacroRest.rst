========
``rest``
========

The 'restart' macro controls the content and format of the initial condition file (filename appears after keyword rsti: in the control file) and the final condition file (filename appears after keyword rsto: in the control file).

* Group 1 -	CHDUM 

+----------------+--------------+--------------------------------------------------------------------------------------------+
| Input Variable | Format       | Description                                                                                |
+================+==============+============================================================================================+
| CHDUM          | character*80 | Keyword(s).  Keywords are entered one per line and terminated with ‘end' or a blank line.  |
+----------------+--------------+--------------------------------------------------------------------------------------------+

Keywords are entered one per line and terminated with ‘end' or a blank line.  

Two keywords 'read' and 'write' are followed by a list of variables; all others stand alone. 

Valid keywords (case insensitive) are:

* Keywords that control format:

  - ‘ascii' - both read and write restart files are ascii (default) 

  - ‘binary' - both read and write restart files are unformatted 

  - ‘rbinary' - unformatted read restart file 

  - ‘wbinary' - unformatted write restart file  

  - 'old' use old format (input / output format and content is hard-wired)

  - 'new' use new format  (input / output format and content is controlled by restart macro keywords)

* Keywords that control flux output:

  - 'noflux' - do not output flux

  - 'flux' - output liquid and vapor flux

  - 'lflux' -  output liquid flux

  - 'vflux' - output vapor flux

* Keywords that control the list of variables to read / write
  (default is to read all variables in restart file; write
  all variables in current simulation)

  - 'read' (followed by list of variables, on same line)

  - 'write' (followed by list of variables, on same line)

Possible read/write variables:

* none

* all

* temp

* pres

* poro

* trac

* ptrk

* gasp

* pini

* saturation

* co2

* mass

* disp (disx, disy, disz)

* strs or stre or strs (strx, stry, strz, stxy, stxz, styz)


The following is an example of rest. In this example restart data will be written
to an unformatted file and liquid flux will be output. If a read restart file
is used it will be in ascii format and if liquid flux data is present in the
file it will be read.

+---------+
| rest    |
+---------+
| wbinary |
+---------+
| lflux   |
+---------+
|         |
+---------+
