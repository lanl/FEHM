========
``cgdp``
========

Assign rate-limited gdpm nodes [i.e. make the connection value large for those nodes that serve only to link gdpm nodes (diffusion only) to the flow field].

Group 1 - wdd1

Group 2 - `JA, JB, JC, IGDPM_RATE_NODES <InputData.html#JA>`_

+------------------+-----------+----------------------------------------------+
| Input Variable   | Format    | Description                                  |
+==================+===========+==============================================+
| macro            | character | Process to be modified: 'heat', 'tran'       |
+------------------+-----------+----------------------------------------------+
| IGDPM_RATE_NODES | integer   | GDPM nodes for which rate should be limited. |
+------------------+-----------+----------------------------------------------+
