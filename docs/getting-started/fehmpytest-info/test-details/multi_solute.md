---
title : multi_solute
layout : page_getting-started
hero_height: is-hidden
---

# multi_solute

**Test Multi-Solute Transport with Chemical Reaction**

The test case is constructed from the VV Test Suite multi_solute model problem. It uses the trac and trxn macros and is a comparison of FEHM and PDREACT Solution for Concentration vs Time.
Test compares the generated tracer files to old tracer files known to be correct. All concentraction values are tested.

Test Cases for Multi-Solute Transport with Chemical Reaction - trac macro

  Species Co_aq              
  Species Fe_aq                     
  Species EDTA_aq                    
  Species CoEDTA_aq              
  Species FeEDTA_aq                      
  Species CoEDTA_s                       
  Species Co_s                          
  Species FeEDTA_s   

Multi-Solute Transport with Chemical Reaction - trxn macro

  Species Co_aq               
  Species Fe_aq                   
  Species EDTA_aq                 
  Species CoEDTA_aq                  
  Species FeEDTA_aq                    
  Species CoEDTA_s                  
  Species Co_s                        
  Species FeEDTA_s    

Comparison files:

<pre>
multi_solute_trac_Cobalt[aq].trc	       multi_solute_trxn_Cobalt_s.trc
multi_solute_trac_Cobalt[s].trc            multi_solute_trxn_Cobalt.trc
multi_solute_trac_Co-EDTA[aq].trc	       multi_solute_trxn_Co-EDTA_a.trc
multi_solute_trac_Co-EDTA[s].trc	       multi_solute_trxn_Co-EDTA_s.trc
multi_solute_trac_EDTA[aq].trc             multi_solute_trxn_EDTA.trc
multi_solute_trac_Fe-EDTA[aq].trc          multi_solute_trxn_Fe-EDTA_a.trc
multi_solute_trac_Fe-EDTA[s].trc	       multi_solute_trxn_Fe-EDTA_s.trc
multi_solute_trac_FreeIon_Cobalt[aq].trc   multi_solute_trxn_FreeIon_Cobalt.trc
multi_solute_trac_FreeIon_EDTA[aq].trc     multi_solute_trxn_FreeIon_EDTA.trc
multi_solute_trac_FreeIon_Iron[aq].trc     multi_solute_trxn_FreeIon_Iron.trc
multi_solute_trac_Iron[aq].trc             multi_solute_trxn_Iron.trc        
</pre>

Test Directory: [FEHM/fehmpytests/multi_solute](https://github.com/lanl/FEHM/tree/master/fehmpytests/multi_solute)


### Example File multi_solute.in
<pre>

COMPARE FEHMN and PDREACT: Linear Sorption w/ Surface Exchange  
cond
1 202 1 2.7 2.7 2.7

ctrl
50 1e-6 8
1 202 1 2

1 0 0.5
25 2. 1.e-10 1.e-1
1 0
flow
1 202 101 -0.05556 -25 0
101 202 101 1. -25 -1

init
1. 25 25 0 1000 25 0 0
node
1
202
hist 
concentration
end
perm
1 202 1 5.0e-13 5.0e-30 5.0e-30

rock
1 202 1 1500 1000 0.4

sol
1 -1
time
1.e-6 7.25 1000 10 92 11

# solute 1: Total Cobalt Concentration
# solute 2: Total Iron Concentration
# solute 3: Total EDTA Concentration
# solute 4: CoEDTA adsorbed concentration
# solute 5: Co adsorbed concentration
# solute 6: FeEDTA adsorbed concentration
trac
1.d-80 1.0 1.e-6 0.5
1. 2000 1.0 2000
5 5.0 1.e-6 2.8333e-3
6
1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 3.1623e-5 1.0 4.16667

1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 1.e-13 1.0 4.16667

1
1 0. 0. 1. 1.e-9 .05 1.e-34 1.e-34

1 202 1 1


1 202 101 3.1623e-5 1.0 4.16667

0


0


0


rxn
** NCPLX, NUMRXN
2,4
** Coupling of the aqueous components (dRi/dUj)
2
1 0 1
0 1 0
** IDCPNT(IC),CPNTNAM(IC),IFXCONC(IC),CPNTPRT(IC) (comp,name,cond.; NCPNT rows)
1    Cobalt[aq]     0     0     1.e-9
2      Iron[aq]     0     0     1.e-9
3      EDTA[aq]     0     0     1.e-9
** IDCPLX(IX), CPLXNAM(IX),CPLXPRT(IX) (ID # and name of complex, NCPLX rows)
101   Co-EDTA[aq]   0
102   Fe-EDTA[aq]   0
** IDIMM(IM), IMMNAM(IM),IMMPRT(IM)(ID # and name of immoblie spec, NIMM rows)
1   Co-EDTA[s]     0
2   Fe-EDTA[s]     0
3    Cobalt[s]     0
** IDVAP(IV), VAPNAM(IM), VAPPRT(IV) (ID # and name of vapor spec, NVAP rows)
** Skip nodes
0
** RSDMAX
1.0e-10
**** Chemical reaction information for equilibrium reactions ******
** LOGKEQ (=0 if stability constants are given as K, =1 if given as log(K))
         0
** CKEQ(IX) ,HEQ(IX) (Stability constants and Enthaplys, NCPLX rows)
1.0e+18    0
6.31e+27   0
** STOIC(IX,IC) (Stoichiometric coeff: NCPLX rows, NCPNT columns)
1.0       0.0       1.0
0.0       1.0       1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
      101      1
** Distribution coeffienct (kg water/ kg rock) **  
    0.533
** Mass transfer coefficient (1/hr) **
      1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
        1      3
** Distribution coeffienct (kg rock/ kg water) **  
     5.07
** Mass transfer coefficient (1/hr) **
      1.0
** LINEAR KINETIC REACTION (type 1) **
        1
** Where does the reaction take place? **
 1 0 0

** Aqueous Component/Complex #, Solid Component #
      102      2
** Distribution coeffienct (kg rock/ kg water) **  
    0.427
** Mass transfer coefficient (1/hr) **
      1.0
** GENERAL EXCHANGE REACTION (type 3) **
        3
** Where does the reaction take place? **
 1 0 0

** # of solid, liquid and vapor species **
        3   0   0
** forward and reverse rate constants (1/hr) **
  1.26e-2 0
** Solid Species in reaction **
  1      2      3
** Stoichiometry **
  1.0   -1.0  -1.0
coor   n/a     
   202
         1        0.00000        1.00000        0.00000
         2        0.10000        1.00000        0.00000
         3        0.20000        1.00000        0.00000
         4        0.30000        1.00000        0.00000
         5        0.40000        1.00000        0.00000
         6        0.50000        1.00000        0.00000
         7        0.60000        1.00000        0.00000
         8        0.70000        1.00000        0.00000
         9        0.80000        1.00000        0.00000
        10        0.90000        1.00000        0.00000
        11        1.00000        1.00000        0.00000
        12        1.10000        1.00000        0.00000
        13        1.20000        1.00000        0.00000
        14        1.30000        1.00000        0.00000
        15        1.40000        1.00000        0.00000
        16        1.50000        1.00000        0.00000
        17        1.60000        1.00000        0.00000
        18        1.70000        1.00000        0.00000
        19        1.80000        1.00000        0.00000
        20        1.90000        1.00000        0.00000
        21        2.00000        1.00000        0.00000
        22        2.10000        1.00000        0.00000
        23        2.20000        1.00000        0.00000
        24        2.30000        1.00000        0.00000
        25        2.40000        1.00000        0.00000
        26        2.50000        1.00000        0.00000
        27        2.60000        1.00000        0.00000
        28        2.70000        1.00000        0.00000
        29        2.80000        1.00000        0.00000
        30        2.90000        1.00000        0.00000
        31        3.00000        1.00000        0.00000
        32        3.10000        1.00000        0.00000
        33        3.20000        1.00000        0.00000
        34        3.30000        1.00000        0.00000
        35        3.40000        1.00000        0.00000
        36        3.50000        1.00000        0.00000
        37        3.60000        1.00000        0.00000
        38        3.70000        1.00000        0.00000
        39        3.80000        1.00000        0.00000
        40        3.90000        1.00000        0.00000
        41        4.00000        1.00000        0.00000
        42        4.10000        1.00000        0.00000
        43        4.20000        1.00000        0.00000
        44        4.30000        1.00000        0.00000
        45        4.40000        1.00000        0.00000
        46        4.50000        1.00000        0.00000
        47        4.60000        1.00000        0.00000
        48        4.70000        1.00000        0.00000
        49        4.80000        1.00000        0.00000
        50        4.90000        1.00000        0.00000
        51        5.00000        1.00000        0.00000
        52        5.10000        1.00000        0.00000
        53        5.20000        1.00000        0.00000
        54        5.30000        1.00000        0.00000
        55        5.40000        1.00000        0.00000
        56        5.50000        1.00000        0.00000
        57        5.60000        1.00000        0.00000
        58        5.70000        1.00000        0.00000
        59        5.80000        1.00000        0.00000
        60        5.90000        1.00000        0.00000
        61        6.00000        1.00000        0.00000
        62        6.10000        1.00000        0.00000
        63        6.20000        1.00000        0.00000
        64        6.30000        1.00000        0.00000
        65        6.40000        1.00000        0.00000
        66        6.50000        1.00000        0.00000
        67        6.60000        1.00000        0.00000
        68        6.70000        1.00000        0.00000
        69        6.80000        1.00000        0.00000
        70        6.90000        1.00000        0.00000
        71        7.00000        1.00000        0.00000
        72        7.10000        1.00000        0.00000
        73        7.20000        1.00000        0.00000
        74        7.30000        1.00000        0.00000
        75        7.40000        1.00000        0.00000
        76        7.50000        1.00000        0.00000
        77        7.60000        1.00000        0.00000
        78        7.70000        1.00000        0.00000
        79        7.80000        1.00000        0.00000
        80        7.90000        1.00000        0.00000
        81        8.00000        1.00000        0.00000
        82        8.10000        1.00000        0.00000
        83        8.20000        1.00000        0.00000
        84        8.30000        1.00000        0.00000
        85        8.40000        1.00000        0.00000
        86        8.50000        1.00000        0.00000
        87        8.60000        1.00000        0.00000
        88        8.70000        1.00000        0.00000
        89        8.80000        1.00000        0.00000
        90        8.90000        1.00000        0.00000
        91        9.00000        1.00000        0.00000
        92        9.10000        1.00000        0.00000
        93        9.20000        1.00000        0.00000
        94        9.30000        1.00000        0.00000
        95        9.40000        1.00000        0.00000
        96        9.50000        1.00000        0.00000
        97        9.60000        1.00000        0.00000
        98        9.70000        1.00000        0.00000
        99        9.80000        1.00000        0.00000
       100        9.90000        1.00000        0.00000
       101       10.00000        1.00000        0.00000
       102        0.00000        0.00000        0.00000
       103        0.10000        0.00000        0.00000
       104        0.20000        0.00000        0.00000
       105        0.30000        0.00000        0.00000
       106        0.40000        0.00000        0.00000
       107        0.50000        0.00000        0.00000
       108        0.60000        0.00000        0.00000
       109        0.70000        0.00000        0.00000
       110        0.80000        0.00000        0.00000
       111        0.90000        0.00000        0.00000
       112        1.00000        0.00000        0.00000
       113        1.10000        0.00000        0.00000
       114        1.20000        0.00000        0.00000
       115        1.30000        0.00000        0.00000
       116        1.40000        0.00000        0.00000
       117        1.50000        0.00000        0.00000
       118        1.60000        0.00000        0.00000
       119        1.70000        0.00000        0.00000
       120        1.80000        0.00000        0.00000
       121        1.90000        0.00000        0.00000
       122        2.00000        0.00000        0.00000
       123        2.10000        0.00000        0.00000
       124        2.20000        0.00000        0.00000
       125        2.30000        0.00000        0.00000
       126        2.40000        0.00000        0.00000
       127        2.50000        0.00000        0.00000
       128        2.60000        0.00000        0.00000
       129        2.70000        0.00000        0.00000
       130        2.80000        0.00000        0.00000
       131        2.90000        0.00000        0.00000
       132        3.00000        0.00000        0.00000
       133        3.10000        0.00000        0.00000
       134        3.20000        0.00000        0.00000
       135        3.30000        0.00000        0.00000
       136        3.40000        0.00000        0.00000
       137        3.50000        0.00000        0.00000
       138        3.60000        0.00000        0.00000
       139        3.70000        0.00000        0.00000
       140        3.80000        0.00000        0.00000
       141        3.90000        0.00000        0.00000
       142        4.00000        0.00000        0.00000
       143        4.10000        0.00000        0.00000
       144        4.20000        0.00000        0.00000
       145        4.30000        0.00000        0.00000
       146        4.40000        0.00000        0.00000
       147        4.50000        0.00000        0.00000
       148        4.60000        0.00000        0.00000
       149        4.70000        0.00000        0.00000
       150        4.80000        0.00000        0.00000
       151        4.90000        0.00000        0.00000
       152        5.00000        0.00000        0.00000
       153        5.10000        0.00000        0.00000
       154        5.20000        0.00000        0.00000
       155        5.30000        0.00000        0.00000
       156        5.40000        0.00000        0.00000
       157        5.50000        0.00000        0.00000
       158        5.60000        0.00000        0.00000
       159        5.70000        0.00000        0.00000
       160        5.80000        0.00000        0.00000
       161        5.90000        0.00000        0.00000
       162        6.00000        0.00000        0.00000
       163        6.10000        0.00000        0.00000
       164        6.20000        0.00000        0.00000
       165        6.30000        0.00000        0.00000
       166        6.40000        0.00000        0.00000
       167        6.50000        0.00000        0.00000
       168        6.60000        0.00000        0.00000
       169        6.70000        0.00000        0.00000
       170        6.80000        0.00000        0.00000
       171        6.90000        0.00000        0.00000
       172        7.00000        0.00000        0.00000
       173        7.10000        0.00000        0.00000
       174        7.20000        0.00000        0.00000
       175        7.30000        0.00000        0.00000
       176        7.40000        0.00000        0.00000
       177        7.50000        0.00000        0.00000
       178        7.60000        0.00000        0.00000
       179        7.70000        0.00000        0.00000
       180        7.80000        0.00000        0.00000
       181        7.90000        0.00000        0.00000
       182        8.00000        0.00000        0.00000
       183        8.10000        0.00000        0.00000
       184        8.20000        0.00000        0.00000
       185        8.30000        0.00000        0.00000
       186        8.40000        0.00000        0.00000
       187        8.50000        0.00000        0.00000
       188        8.60000        0.00000        0.00000
       189        8.70000        0.00000        0.00000
       190        8.80000        0.00000        0.00000
       191        8.90000        0.00000        0.00000
       192        9.00000        0.00000        0.00000
       193        9.10000        0.00000        0.00000
       194        9.20000        0.00000        0.00000
       195        9.30000        0.00000        0.00000
       196        9.40000        0.00000        0.00000
       197        9.50000        0.00000        0.00000
       198        9.60000        0.00000        0.00000
       199        9.70000        0.00000        0.00000
       200        9.80000        0.00000        0.00000
       201        9.90000        0.00000        0.00000
       202       10.00000        0.00000        0.00000
         0        0.00000        0.00000        0.00000
elem
   4  100  0
       1     102     103       2       1
       2     103     104       3       2
       3     104     105       4       3
       4     105     106       5       4
       5     106     107       6       5
       6     107     108       7       6
       7     108     109       8       7
       8     109     110       9       8
       9     110     111      10       9
      10     111     112      11      10
      11     112     113      12      11
      12     113     114      13      12
      13     114     115      14      13
      14     115     116      15      14
      15     116     117      16      15
      16     117     118      17      16
      17     118     119      18      17
      18     119     120      19      18
      19     120     121      20      19
      20     121     122      21      20
      21     122     123      22      21
      22     123     124      23      22
      23     124     125      24      23
      24     125     126      25      24
      25     126     127      26      25
      26     127     128      27      26
      27     128     129      28      27
      28     129     130      29      28
      29     130     131      30      29
      30     131     132      31      30
      31     132     133      32      31
      32     133     134      33      32
      33     134     135      34      33
      34     135     136      35      34
      35     136     137      36      35
      36     137     138      37      36
      37     138     139      38      37
      38     139     140      39      38
      39     140     141      40      39
      40     141     142      41      40
      41     142     143      42      41
      42     143     144      43      42
      43     144     145      44      43
      44     145     146      45      44
      45     146     147      46      45
      46     147     148      47      46
      47     148     149      48      47
      48     149     150      49      48
      49     150     151      50      49
      50     151     152      51      50
      51     152     153      52      51
      52     153     154      53      52
      53     154     155      54      53
      54     155     156      55      54
      55     156     157      56      55
      56     157     158      57      56
      57     158     159      58      57
      58     159     160      59      58
      59     160     161      60      59
      60     161     162      61      60
      61     162     163      62      61
      62     163     164      63      62
      63     164     165      64      63
      64     165     166      65      64
      65     166     167      66      65
      66     167     168      67      66
      67     168     169      68      67
      68     169     170      69      68
      69     170     171      70      69
      70     171     172      71      70
      71     172     173      72      71
      72     173     174      73      72
      73     174     175      74      73
      74     175     176      75      74
      75     176     177      76      75
      76     177     178      77      76
      77     178     179      78      77
      78     179     180      79      78
      79     180     181      80      79
      80     181     182      81      80
      81     182     183      82      81
      82     183     184      83      82
      83     184     185      84      83
      84     185     186      85      84
      85     186     187      86      85
      86     187     188      87      86
      87     188     189      88      87
      88     189     190      89      88
      89     190     191      90      89
      90     191     192      91      90
      91     192     193      92      91
      92     193     194      93      92
      93     194     195      94      93
      94     195     196      95      94
      95     196     197      96      95
      96     197     198      97      96
      97     198     199      98      97
      98     199     200      99      98
      99     200     201     100      99
     100     201     202     101     100
       0       0       0       0       0
stop

</pre>
