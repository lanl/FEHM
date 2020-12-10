========
``wflo``
========

Create a new flow macro to represent boundary conditions on an extracted submodel. Alternate submodel boundary output.

* Group 1 - KEYMS1, KEYMS2, KEYMS3, KEYMS3

or, if KEYMS4 = ‘type':

* Group 1 - KEYMS1, KEYMS2, KEYMS3, KEYMS4, ITYPSD

Enter one line per model defined, terminated with a blank line.

* Group 2 - JA, JB, JC, ISUBMD

+----------------+-----------+-------------------------------------------------------------------------------------+
| Input Variable | Format    | Description                                                                         |
+================+===========+=====================================================================================+
| KEYMS1         | character | | Macro SKD type. The letters given in ( ) are sufficient to identify the keyword:  |
|                |           | | (p)ressure                                                                        |
|                |           | | (h)ead                                                                            |
|                |           | | (f)lux                                                                            |
+----------------+-----------+-------------------------------------------------------------------------------------+
| KEYMS2         | character | | Macro ESKD type. The letters given in ( ) are sufficient to identify the keyword: |
|                |           | | (s)aturation                                                                      |
|                |           | | (t)emperature                                                                     |
|                |           | | (e)nthalpy                                                                        |
|                |           | | (w)ater (water only source, output saturation as 1.0)                             |
+----------------+-----------+-------------------------------------------------------------------------------------+
| KEYMS3         | character | | Macro AIPED type:                                                                 |
|                |           | | imph  - Impedance parameter = 1.0                                                 |
|                |           | | impl - Impedance parameter = 1.0e-4                                               |
|                |           | | impn  - Impedance parameter = -1.0                                                |
|                |           | | If KEYSM1 = ‘flux' the impedance parameter = 0.0, otherwise for any other input,  |
|                |           |   the impedance parameter = 1.0e2.                                                  |
+----------------+-----------+-------------------------------------------------------------------------------------+
| KEYMS4         | character | | Flag to indicate submodel type will be input:                                     |
|                |           | | type                                                                              |
+----------------+-----------+-------------------------------------------------------------------------------------+
| ITYPSD         | INTEGER   | | Submodel type:                                                                    |
|                |           | | :math:`ITYPSD = 0`, generate 'flow' macro.                                        |
|                |           | | :math:`ITYPSD \ne 0`, generate 'flo3' macro.                                      |
+----------------+-----------+-------------------------------------------------------------------------------------+
| ISUBMD         | integer   | Submodel assignment.                                                                |
+----------------+-----------+-------------------------------------------------------------------------------------+

The following is an example of the wflo macro. A flow macro will be written for two output zones.
Because, they are specified using two models, the macros will be written to two separate files.
The macro file names will be generated using the root file name (as input or determined from the output
file name) appended with the model number and suffix ``.wflow`` (e.g. ``file.0001.wflow`` and ``file.0002.wflow``).

+------+-----+------+----+
| wflo |     |      |    |
+------+-----+------+----+
| flux | wat | impl | na |
+------+-----+------+----+
| flux | wat | impl | na |
+------+-----+------+----+
|      |     |      |    |
+------+-----+------+----+
| -308 | 0   | 0    | 1  |
+------+-----+------+----+
| -309 | 0   | 0    | 2  |
+------+-----+------+----+
|      |     |      |    |
+------+-----+------+----+

