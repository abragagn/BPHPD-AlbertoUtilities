#!/bin/sh
rm ./albertoCode.tar.gz ./albertoCodeNoExamples.tar.gz
tar -czvf ./albertoCode.tar.gz ./PDAnalysis/
tar -czvf ./albertoCodeNoExamples.tar.gz ./PDAnalysis/Ntu/BuildFile.xml ./PDAnalysis/Ntu/bin/BuildFile.xml ./PDAnalysis/Ntu/bin/OSMuonMvaTag.* ./PDAnalysis/Ntu/bin/TopDecayClassifier.* ./PDAnalysis/Ntu/bin/AlbertoUtil.* ./PDAnalysis/Ntu/bin/TopDecayMode.h ./PDAnalysis/Ntu/bin/PDSoftMuonMvaEstimator.h*
scp -2 -P 4444 -oNoHostAuthenticationForLocalhost=yes ./albertoCode*  abragagn@localhost:/lustre/cmswork/abragagn/BPH/

git commit -a -m "Update"
git push
