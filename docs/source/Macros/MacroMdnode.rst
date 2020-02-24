==========
``mdnode``
==========

Enables extra connections to be made to nodes. This is useful for simulating wellbore connections, faults, and flow across internal boundaries.

* Group 1 -	NUM_MD, MAX_CON, IELIM, SX_MULT

* Group 2 -	NODE, IPAR, NPAR (repeated NUM_MD times) 

+----------------+---------+---------+------------------------------------------------------+
| Input Variable | Format  | Default | Description                                          |
+================+=========+=========+======================================================+
| NUM_MD         | integer | 0       | Number of new connections to be entered.             |
+----------------+---------+---------+------------------------------------------------------+
| MDMAX          | integer | 0       | Maximum number of new connections to a given node.   |
|                |         |         | This does not include old connections. Thus,         |
|                |         |         | if a node was already connected to 5 neighboring     |
|                |         |         | nodes and two new connections were added to this     |
|                |         |         | node in this macro statement and this was the        |
|                |         |         | maximum number of connections added in this          |
|                |         |         | macro statement, then MDMAX = 2.                     |
+----------------+---------+---------+------------------------------------------------------+
| I_ELIM         | integer | 0       | IF I_ELIM Š 0, then no action. IF I_ELIM < 0,        |
|                |         |         | then nodal connections are eliminated as needed      |
|                |         |         | if redundant.                                        |
+----------------+---------+---------+------------------------------------------------------+
| SX_MULT        | real*8  | 1.0     | Multiplier for equilibrium conditions.               |
+----------------+---------+---------+------------------------------------------------------+
| NODE           | integer | 0       | Node to which new connection is established.         |
+----------------+---------+---------+------------------------------------------------------+
| IPAR           | integer | 0       | IPAR is not used at present. Its value is ignored.   |
|                |         |         | However the entered number must be an integer.       |
+----------------+---------+---------+------------------------------------------------------+
| NPAR           | integer | 0       | NPAR is the new connected node. If NPAR = NODE,      |
|                |         |         | no new connection is established.                    |
+----------------+---------+---------+------------------------------------------------------+

The following are examples of mdnode. In the first example (top), 3 new connections
are specified, node 10 is connected to node 15, node 100 is connected to node
106, and node 10 is connected to node 320. A maximum of 2 new connections are
allowed per node. The multiplier for equilibrium conditions is set to 10. In the
second example (bottom), 4 new connections are specified, node 1 is connected to
node 16, node 2 is connected to node 1, node 4 is connected to node 1 and node
10 is connected to node 203. A maximum of 3 new connections are allowed per node.
The multiplier for equilibrium conditions is set to 100. 


+--------+---+-----+----+
| mdnode |   |     |    |
+--------+---+-----+----+
| 3      | 2 | 0   | 10 |
+--------+---+-----+----+
| 10     | 0 | 15  |    |
+--------+---+-----+----+
| 100    | 0 | 106 |    |
+--------+---+-----+----+
| 10     | 0 | 320 |    |
+--------+---+-----+----+


+--------+---+-----+-----+
| mdnode |   |     |     |
+--------+---+-----+-----+
| 4      | 3 | 0   | 100 |
+--------+---+-----+-----+
| 1      | 0 | 16  |     |
+--------+---+-----+-----+
| 2      | 0 | 1   |     |
+--------+---+-----+-----+
| 4      | 0 | 1   |     |
+--------+---+-----+-----+
| 10     | 0 | 203 |     |
+--------+---+-----+-----+

