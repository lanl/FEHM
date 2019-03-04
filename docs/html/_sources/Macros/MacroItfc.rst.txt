========
``itfc``
========

Data to define flow and transport parameters at interfaces between pairs of zones.

* Group 1 -	ZONE_PAIR(I,1), ZONE_PAIR(I,2), RED_FACTOR(I)- an arbitrary number of lines of input, terminated by a blank line.

* Group 2 -	(FILTER_FLAG(J), J= 1,NSPECI)

* Group 3 -	ZONEC_PAIR(K,1), ZONEC_PAIR(K,2), FTN_FACTOR(K)- an arbitrary number of lines of input, terminated by a blank line.

  - KEYWORD ‘file'

  - SFILENAME

  - ITFCPORSIZE(I), ITFCPROBSIZE(I)- an arbitrary number of lines of input, terminated by a blank line.

+----------------+--------------+-----------------------------------------------------------+
| Input Variable | Format       | Description                                               |
+================+==============+===========================================================+
| ZONE_PAIR      | integer      | Zone number for the zones for which the code              |
|                |              | identifies the interface connections when                 |
|                |              | applying the permeability reduction factor.               |
+----------------+--------------+-----------------------------------------------------------+
| RED_FACTOR     | real         | Reduction factor multiplying the harmonically             |
|                |              | weighted saturated permeability for all connections       |
|                |              | at the interface identified by ZONE_PAIR                  |
+----------------+--------------+-----------------------------------------------------------+
| FILTER_FLAG    | integer      | FEHM has a provision to apply transport mechanisms        |
|                |              | for size exclusion or filtration at interfaces            |
|                |              | defined in the itfc macro. These provisions can           |
|                |              | be used to simulate conditions in which, for              |
|                |              | example, abrupt changes in properties occur at            |
|                |              | interfaces, or hydrologic conditions not                  |
|                |              | explicitly incorporated in a model (a thin clay           |
|                |              | layer, for example) are thought to be present             |
|                |              | that affect transport across the interface.               |
|                |              | The means for specifying these interface transport        |
|                |              | conditions is the itfc macro. Thus, this parameter        |
|                |              | is a flag used to distinguish whether the size            |
|                |              | exclusion or filtration is to be implemented              |
|                |              | (a value 1) or not (a value 0) for each species           |
|                |              | identified in the trac, ptrk, or mptr macros.             |
|                |              | The default value is 0. See the definition of             |
|                |              | FTN_FACTOR below for details on how to invoke             |
|                |              | the size exclusion or filtration model.                   |
+----------------+--------------+-----------------------------------------------------------+
| ZONEC_PAIR     | integer      | Zone number for the zones for which the code              |
|                |              | identifies the interface connections when applying        |
|                |              | the transport filtration or size exclusion factors.       |
+----------------+--------------+-----------------------------------------------------------+
| FTN_FACTOR     | real         | Filtration or size exclusion factor applied for           |
|                |              | all connections at the interface identified by            |
|                |              | ZONEC_PAIR. For the trac macro, a size exclusion          |
|                |              | model is implemented, where FTN_FACTOR = 0 (size          |
|                |              | exclusion) or 1 (no exclusion) are options. For           |
|                |              | ptrk or mptr, a filtration model is implemented,          |
|                |              | where the parameter is the probability of the             |
|                |              | particle passing through the interface (if 0,             |
|                |              | filtration is guaranteed; if 1, there is no filtration).  |
|                |              | For the particle tracking model, FTN_FACTOR < 0           |
|                |              | denotes that the pore size distribution is being          |
|                |              | used. This option is used with the particle size          |
|                |              | distribution option in ptrk and mptr, so that             |
|                |              | each particle is assigned a size. The cumulative          |
|                |              | pore size distribution is then used as a                  |
|                |              | probability distribution function, and when a             |
|                |              | particle encounters the interface, a pore size is         |
|                |              | randomly selected from the distribution. If the           |
|                |              | particle is larger than the pore, it is filtered.         |
|                |              | Note that filtered particles remain at that location      |
|                |              | in the model and are no longer transported.               |
+----------------+--------------+-----------------------------------------------------------+
| KEYWORD        | character*4  | Optional keyword ‘file' designating that the pore         |
|                |              | size distribution information is being input in           |
|                |              | a separate file. This input is entered only for           |
|                |              | interfaces in which FTN_FACTOR < 0 is used.               |
+----------------+--------------+-----------------------------------------------------------+
| SFILENAME      | character*80 | Optional file name containing the pore size distribution  |
|                |              | table. This input is entered only for interfaces          |
|                |              | in which FTN_FACTOR < 0 is used.                          |
+----------------+--------------+-----------------------------------------------------------+
| ITFCPORSIZE    | real         | Pore size for this entry of the pore size distribution    |
|                |              | table (paired with a value of ITFCPROBSIZE).              |
|                |              | An arbitrary number of entries can be input,              |
|                |              | terminated with a blank line. These entries are           |
|                |              | located in the file SFILENAME if specified, or            |
|                |              | in the itfc input file if the alternate input             |
|                |              | file is not used. The code decides if particles           |
|                |              | are irreversibly filtered by comparing the                |
|                |              | particle size to the randomly selected pore size.         |
|                |              | This input is entered only for interfaces in which        |
|                |              | FTN_FACTOR < 0 is used.                                   |
+----------------+--------------+-----------------------------------------------------------+
| ITFCPROBSIZE   | real         | Cumulative probability for the distribution of            |
|                |              | pore sizes (paired with a value of ITFCPORSIZE).          |
|                |              | See description of ITFCPORSIZE above for details.         |
|                |              | The final entry of the table must have                    |
|                |              | ITFCPROBSIZE = 1, since the distribution is assumed       |
|                |              | to be normalized to unity. This input is entered          |
|                |              | only for interfaces in which FTN_FACTOR < 0 is used.      |
+----------------+--------------+-----------------------------------------------------------+

Note that data for each numbered group must be input. The other input is optional.
If filtration is not implemented for any species, a single blank line is input
for Groups 2 and 3, signaling the end of itfc input.

The following is an example of itfc. In this example, the permeability reduction
factor of 0.1 is applied to all node connections at the interface between zones
6 and 10, or 6 and 11. 

+------+----+-----+
| itfc |    |     |
+------+----+-----+
| 6    | 10 | 0.1 |
+------+----+-----+
| 6    | 11 | 0.1 |
+------+----+-----+
|      |    |     |
+------+----+-----+
|      |    |     |
+------+----+-----+