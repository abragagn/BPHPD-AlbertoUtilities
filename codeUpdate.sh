#!/bin/sh
rm ./albertoCode.tar.gz ./albertoCodeNoExamples.tar.gz
tar -czvf ./albertoCode.tar.gz ./PDAnalysis/
tar -czvf ./albertoCodeNoExamples.tar.gz ./PDAnalysis/Ntu/BuildFile.xml ./PDAnalysis/Ntu/bin/BuildFile.xml ./PDAnalysis/Ntu/bin/TopDecayClassifier.cc ./PDAnalysis/Ntu/bin/AlbertoUtil.cc ./PDAnalysis/Ntu/bin/TopDecayMode.h ./PDAnalysis/Ntu/bin/script.sh ./PDAnalysis/Ntu/bin/AlbertoUtil.h ./PDAnalysis/Ntu/bin/PDSoftMuonMvaEstimator.h ./PDAnalysis/Ntu/bin/TopDecayClassifier.h ./PDAnalysis/Ntu/bin/PDSoftMuonMvaEstimator.cc
scp ./albertoCode* abragagn@t2-ui-12:/lustre/cmswork/abragagn/BPH/
