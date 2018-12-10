#ifndef PDSoftMuonMvaEstimator_H
#define PDSoftMuonMvaEstimator_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDAnalyzerUtil.h"
#include "PDMuonVar.h"
#include "PDAnalysis/Ntu/interface/PDEnumString.h"

#include "TString.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

class PDSoftMuonMvaEstimator:   public virtual PDAnalyzerUtil
,                               public virtual PDMuonVar {

public:
    PDSoftMuonMvaEstimator();
    ~PDSoftMuonMvaEstimator();

    void inizializeMuonMvaReader(TString methodName, TString path);
    float computeMuonMva(int iMuon);
    bool muonPassedPreselection(int iMuon);
    int IPsign_(int iMuon);
    int IPsign_(int iMuon, int iPV);

protected:
    
    TString methodNameBarrel_, methodNameEndcap_;
    TString weightFileBarrel_, weightFileEndcap_;

private:
    void computeMvaVariables(int iMuon);
    void methodSetup(TString methodName, TString path);
    TString methodNameFromWeightName(TString weightname);
    bool useIp(TString methodName);
    bool useIso(TString methodName);

    TMVA::Reader muonMvaIdReader_;

    // MVA Variable
    float muoPt_;
    float absMuoEta_;

    float muoSegmComp_;
    float muoChi2LM_;
    float muoChi2LP_;
    float muoGlbTrackTailProb_;
    float muoIValFrac_;
    float muoLWH_;
    float muoTrkKink_;
    float muoGlbKinkFinderLOG_;
    float muoTimeAtIpInOutErr_;
    float muoOuterChi2_;
    float muoInnerChi2_ ;
    float muoTrkRelChi2_;
    float muoVMuonHitComb_;

    float muoGlbDeltaEtaPhi_;
    float muoStaRelChi2_;
    float muoTimeAtIpInOut_;
    float muoValPixHits_;
    float muoNTrkVHits_;
    float muoGNchi2_;
    float muoVMuHits_;
    float muoNumMatches_;


    float trkDxy_;
    float trkDz_;

    float muoPFiso_;
    float muoQprod_;

    float DUMMY_;

};

#endif
