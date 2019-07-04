#ifndef OSMuonMvaTag_H
#define OSMuonMvaTag_H

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDSoftMuonMvaEstimator.h"
#include "PDAnalyzerUtil.h"
#include "AlbertoUtil.h"

#include "TString.h"
#include "TGraph.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

class OSMuonMvaTag : public virtual PDAnalyzerUtil
,                    public virtual AlbertoUtil
{

public:
    OSMuonMvaTag();
    ~OSMuonMvaTag();

    void    inizializeOsMuonTagVars();

    bool    makeOsMuonTagging();
    int     selectOsMuon(); 

    int     getOsMuon(){ return osMuonIndex_; }
    int     getOsMuonTag(){ return osMuonTagDecision_; }
    float   getOsMuonTagMvaValue(){ return osMuonTagMvaValue_; }
    float   getOsMuonTagMistagProbRaw(){ return osMuonTagMistagProbRaw_; }
    float   getOsMuonTagMistagProbCalProcess(){ return osMuonTagMistagProbCalProcess_; }
    float   getOsMuonTagMistagProbCalProcessBuBs(){ return osMuonTagMistagProbCalProcessBuBs_; }

    void    setVtxOsMuonTag(int iB, int iPV) { ssIndex_ = iB; pvIndex_ = iPV;}
    void    setOsMuonMvaCut(float wp);
    void    setOsMuonDzCut(float dzCut);
    void    inizializeOSMuonMvaReader(TString, TString);
    bool    inizializeOSMuonCalibration(TString process, TString processBuMC, TString processBsMC, TString methodPath);

    int     getNosMuons(){ return nMuonsSel_; }

private:
    
    TString methodNameFromWeightName();
    void    computeOsMuonTagVariables();

    TMVA::Reader osMuonTagReader_;
    TString weightsFile_;
    TString methodName_;
    TString methodPath_ = "/lustre/cmswork/abragagn/mvaWeights/OsMuonTag/";

    int ssIndex_;
    int pvIndex_;
    int osMuonIndex_;
    int osMuonTrackIndex_;
    int osMuonTagDecision_;

    float osMuonTagMvaValue_;
    float osMuonTagMistagProbRaw_;
    float osMuonTagMistagProbCalProcess_;
    float osMuonTagMistagProbCalProcessBuBs_;

    float wp_;
    float dzCut_;
    float PFIsoCut_;

    int nMuonsSel_;

    //MVA Variables
    float muoPt_;
    float muoEta_;
    float muoDxy_;
    float muoExy_;
    float muoDz_;
    float muoEz_;
    float muoSoftMvaValue_;
    float muoDrB_;
    float muoPFIso_;
    float muoConeCleanPt_;
    float muoConeCleanPtRel_;
    float muoConeCleanDr_;
    float muoConeCleanEnergyRatio_;
    float muoConeCleanQ_;
    float muoCharge_;

    float DUMMY_;

    //MISTAG VARIABLES
    TF1 *wCalProcess_;
    TF1 *wCalBuMC_;
    TF1 *wCalBsMC_;
    TF1 *wCalBuBs_;

};

#endif
