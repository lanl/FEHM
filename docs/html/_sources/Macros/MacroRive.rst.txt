====================
``rive`` or ``well``
====================

River or implicit well package.

* Group 1 - KEYWORD

  KEYWORD "wellmodel"

  NRIVER, IRIVER

  If IRIVER = 1

  	INRIVERF, ISRIVERF, IFRIVERF,IWSP

  KEYWORD "*end macro*" (required)

The input is terminated with keyword "*end rive*", "*endrive*", "*end well*" or
"*endwell*", where the macro name should match the input macro name that was used.

+-------------------------+-----------+--------------------------------------------------------------------------------------+
| Input Variable          | Format    | Description                                                                          |
+=========================+===========+======================================================================================+
| KEYWORD                 | character |                                                                                      |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| KEYWORD "*wellmodel*"   |           |                                                                                      |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| NRIVER                  | integer   | Number of models defined for this call.                                              |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| IRIVER                  | integer   | | Type of surface/well flow                                                          |
|                         |           | | IRIVER = 0 no routing, ponding, groundwater connection                             |
|                         |           | | IRIVER = 1 simple fluid routing, no groundwater connection                         |
|                         |           | | IRIVER = 2 simple stream definition, no groundwater connection                     |
|                         |           | | IRIVER = 3 simple stream definition, with groundwater connection                   |
|                         |           | | IRIVER = 4 simple stream definition, with groundwater connection (with ponding)    |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| INRIVERF                |           | Section number (id) of ith section.                                                  |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| ISRIVERF                |           | Number of layers for the ith section.                                                |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| IFRIVERF                |           | Number of coordinate points in the ith section.                                      |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| IWSP                    |           |                                                                                      |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
| KEYWORD "*endwell*"     |           |                                                                                      |
+-------------------------+-----------+--------------------------------------------------------------------------------------+
