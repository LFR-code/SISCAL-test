# verbose
0
# nT
44
# nFisheries
6
# nAges
30
# nSplitAge
1
# splitAge
1
# a1F
5
# a2F
8
# baranovIter
1
## Fishery order is:
## LL.3 LL.4 OT.3 OT.4 RV HS
## Selectivity options-----------
## Dome (1=yes, 0=no)
# domeSelectivity
0 0 1 1 1 0
# dM
0.23 0.23 1.26 1.26 0. 0.

## Minimization phases--------------------------
# ph_log_avgR
1
# ph_logit_h
-5
# ph_log_SSB0
-5
# ph_initRecDevs
-2
# ph_recDevs
2
# ph_M
-6
# ph_sel50
5 5 -5 -5 4 4

## Selectivity offSet parameters.
## use < 0 in rows 2,3 for
## asymptotic form
# ph_selOffsets
-5 -5 -5 -5 -5 -5
-7 -7 -7 -7 -7 -7
-7 -7 -8 -8 -8 -7
# ph_selDevs_S50_a
-7
# ph_selDevs_S95_a
-7
## Keep these last devs off all the time
# ph_selDevs_S95_d
-7
# ph_selDevs_S50_d
-7

## Fishing mortality prior to 1970
# ph_Finit
3

## Prior means and std deviations--------------
# priorMean_h
0.75
# priorSD_h
0.01
# prior_initM
0.145
# priorSD_initM
0.001
# priorSD_rwM
0.0005
# priorSD_R
0.25

## Selectivity random-walk std devs--------------
# sd_S50_a
0.1 0.1 0.1 0.1 0.1 0.1
# sd_S95_a
0.1 0.1 0.1 0.1 0.1 0.1
# sd_S50_d
0.1 0.1 0.1 0.1 0.1 0.1
# sd_S95_d
0.1 0.1 0.1 0.1 0.1 0.1

## Maturity--------------
# aMat_50
8.5
# aMat_95
11.5
## Growth parameters--------------
# lInf_m
134
# lInf_f
205
# vonK_m
0.18
# vonK_f
0.12
# wt_a
0.00673
# wt_b
3.12
# cvL_m
0.15
# cvL_f
0.15
