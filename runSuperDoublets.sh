#!/bin/bash

while [ X$# != X0 ] ; do
    case ${1} in
        -n) shift; n=$1; shift;;
        -mm) shift; mm=$1; shift;;
        -D)  shift; D=$1; shift;;
        -us)  shift; us=$1; shift;;
    esac
done

echo -e ".L doubletAnalysis.C++ \n .L loadTTa.C \n\
t9 = new TChain(\"trkTree/tree\"); loadTTa_v2(t9); ScanChainMockSuperDoublets(t9, ${n}, true, ${mm}, ${D}, ${us});\n.qqqq"\
    | root -l -b >& a_ttbar_stats_mm${mm}_sd${D}_us${us}_ntv2.log

