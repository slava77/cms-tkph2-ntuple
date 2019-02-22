#!/bin/bash

while [ X$# != X0 ] ; do
    case ${1} in
        -n) shift; n=$1; shift;;
        -mm) shift; mm=$1; shift;;
        -DB)  shift; DB=$1; shift;;
        -DE)  shift; DE=$1; shift;;
        -us)  shift; us=$1; shift;;
    esac
done

echo -e ".L doubletAnalysis.C++ \n\
t9 = new TChain(\"trackingNtuple/tree\"); t9->Add(\"/Data/cms/tkph2/stublink/CMSSW_9_1_0_pre1-tkNtuple/2023PU140/20634.0_TTbar_14TeV+TTbar_14TeV_TuneCUETP8M1_2023D10PU_GenSimHLBeamSpotFull14+DigiFullTriggerPU_2023D10PU+RecoFullGlobalPU_2023D10PU+HARVESTFullGlobalPU_2023D10PU/trackingNtuple.root\"); ScanChainMockSuperDoublets(t9, ${n}, true, ${mm}, ${DB}, ${DE}, ${us}, false, true, true);\n.qqqq"\
    | root -l -b >& a_ttbar_stats_91X_2023PU140_20634.0_mm${mm}_sd${DB}cm${DE}cm_us${us}.log

