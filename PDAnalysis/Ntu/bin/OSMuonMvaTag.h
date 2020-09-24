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
    bool    makeOsMuonTaggingNoMistag();
    int     selectOsMuon();
    void    computeOsMuonTagVariables();

    int     getOsMuon(){ return osMuonIndex_; }
    int     getOsMuonTag(){ return osMuonTagDecision_; }
    double  getOsMuonTagMvaValue(){ return osMuonTagMvaValue_; }
    double  getOsMuonTagMistagProbRaw(){ return osMuonTagMistagProbRaw_; }
    double  getOsMuonTagMistagProbCalProcess(){ return osMuonTagMistagProbCalProcess_; }
    double  getOsMuonTagMistagProbCalProcessBuBs(){ return osMuonTagMistagProbCalProcessBuBs_; }

    void    setVtxOsMuonTag(int iB, int iPV) { ssIndex_ = iB; pvIndex_ = iPV;}
    void    setOsMuonMvaCut(double wp);
    void    setOsMuonDzCut(double dzCut);
    void    inizializeOSMuonMvaReader(TString, TString);
    bool    inizializeOSMuonCalibration(TString process, TString processBuMC, TString processBsMC, TString methodPath);

    int     getNosMuons(){ return nMuonsSel_; }

    float getMuoPt(){ return muoPt_;};
    float getMuoEta(){ return muoEta_;};
    float getMuoDxy(){ return muoDxy_;};
    float getMuoExy(){ return muoExy_;};
    float getMuoDz(){ return muoDz_;};
    float getMuoEz(){ return muoEz_;};
    float getMuoSoftMvaValue(){ return muoSoftMvaValue_;};
    float getMuoDrB(){ return muoDrB_;};
    float getMuoPFIso(){ return muoPFIso_;};
    float getMuoConeCleanPt(){ return muoConeCleanPt_;};
    float getMuoConeCleanPtRel(){ return muoConeCleanPtRel_;};
    float getMuoConeCleanDr(){ return muoConeCleanDr_;};
    float getMuoConeCleanEnergyRatio(){ return muoConeCleanEnergyRatio_;};
    float getMuoConeCleanQ(){ return muoConeCleanQ_;};
    float getMuoCharge(){ return muoCharge_;};

private:
    
    TString methodNameFromWeightName();

    TMVA::Reader osMuonTagReader_;
    TString weightsFile_;
    TString methodName_;
    TString methodPath_ = "/lustre/cmswork/abragagn/mvaWeights/OsMuonTag/";

    int ssIndex_;
    int pvIndex_;
    int osMuonIndex_;
    int osMuonTrackIndex_;
    int osMuonTagDecision_;

    double osMuonTagMvaValue_;
    double osMuonTagMistagProbRaw_;
    double osMuonTagMistagProbCalProcess_;
    double osMuonTagMistagProbCalProcessBuBs_;

    double wp_;
    double dzCut_;
    double PFIsoCut_;

    int nMuonsSel_;

    //MVA Variables (have to be float)
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
