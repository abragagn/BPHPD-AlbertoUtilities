#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

using namespace std;

OSMuonMvaTag::OSMuonMvaTag():
                osMuonTagReader_("!Color:Silent")
,               ssIndex_(-1)
,               pvIndex_(-1)
,               osMuonIndex_(-1)
,               osMuonTrackIndex_(-1)
,               osMuonTagDecision_(0)
,               osMuonTagMvaValue_(-1.)
,               osMuonTagMistagProbRaw_(-1.)
,               osMuonTagMistagProbCalProcess_(-1.)
,               osMuonTagMistagProbCalProcessBuBs_(-1.)
,               wp_(0.21)
,               dzCut_(1.)
,               nMuonsSel_(0)
{}

OSMuonMvaTag::~OSMuonMvaTag() {}

// =====================================================================================
void OSMuonMvaTag::inizializeOsMuonTagVars()
{
    ssIndex_ = -1;
    pvIndex_ = -1;
    osMuonIndex_ = -1;
    osMuonTrackIndex_ = -1;
    osMuonTagDecision_ = 0;
    osMuonTagMvaValue_ = -1;
    osMuonTagMistagProbRaw_ = -1;
    osMuonTagMistagProbCalProcess_ = -1;
    osMuonTagMistagProbCalProcessBuBs_ = -1;
    nMuonsSel_ = 0;
}

void OSMuonMvaTag::setOsMuonMvaCut(float wp = 0.21)
{
    wp_ = wp;
}

void OSMuonMvaTag::setOsMuonDzCut(float dzCut = 1.)
{
    dzCut_ = dzCut;
}

void OSMuonMvaTag::inizializeOSMuonMvaReader(
    TString methodName = "DNNOsMuonHLTJpsiMu"
,   TString methodPath = ""
    )
{

    if(methodPath == "") methodPath = methodPath_;
    TMVA::PyMethodBase::PyInitialize();
    methodName_ = methodName;

    osMuonTagReader_.AddVariable("muoPt", &muoPt_);
    osMuonTagReader_.AddVariable("muoEta", &muoEta_);
    osMuonTagReader_.AddVariable("muoDxy", &muoDxy_);
    osMuonTagReader_.AddVariable("muoExy", &muoExy_);
    osMuonTagReader_.AddVariable("muoDz", &muoDz_);
    osMuonTagReader_.AddVariable("muoEz", &muoEz_);
    osMuonTagReader_.AddVariable("muoSoftMvaValue", &muoSoftMvaValue_);
    osMuonTagReader_.AddVariable("muoDrB", &muoDrB_);
    osMuonTagReader_.AddVariable("muoPFIso", &muoPFIso_);
    osMuonTagReader_.AddVariable("muoConeCleanPt", &muoConeCleanPt_);
    osMuonTagReader_.AddVariable("muoConeCleanPtRel", &muoConeCleanPtRel_);
    osMuonTagReader_.AddVariable("muoConeCleanDr", &muoConeCleanDr_);
    osMuonTagReader_.AddVariable("muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio_);
    osMuonTagReader_.AddVariable("muoConeCleanQ", &muoConeCleanQ_);
    osMuonTagReader_.BookMVA( methodName_, methodPath + "TMVAClassification_" + methodName_ + ".weights.xml" );
}

bool OSMuonMvaTag::inizializeOSMuonCalibration( 
    TString process = "BuJPsiKData2018"
,   TString processBuMC = "BuJPsiKMC2018"
,   TString processBsMC = "BsJPsiPhiMC2018"
,   TString methodPath = ""  
)
{
    if(methodPath == "") methodPath = methodPath_;
    auto *f   = new TFile(methodPath + "OSMuonTaggerCalibration" + process + ".root");
    auto *fBu = new TFile(methodPath + "OSMuonTaggerCalibration" + processBuMC + ".root");
    auto *fBs = new TFile(methodPath + "OSMuonTaggerCalibration" + processBsMC + ".root");

    if(f->IsZombie() || fBu->IsZombie() || fBs->IsZombie()) return false;
    f->cd();
    wCalProcess_ = (TF1*)f->Get("osMuonCal");
    f->Close();

    fBu->cd();
    wCalBuMC_ = (TF1*)f->Get("osMuonCal");
    fBu->Close();

    fBs->cd();
    wCalBsMC_ = (TF1*)f->Get("osMuonCal");
    fBs->Close();

    wCalBuBs_ = new TF1("osMuonCalBuBs","[0]-[1]*[2]/[3]+[2]/[3]*x",0.,1.);
    float qs = wCalBsMC_->GetParameter(0);
    float ms = wCalBsMC_->GetParameter(1);
    float qu = wCalBuMC_->GetParameter(0);
    float mu = wCalBuMC_->GetParameter(1);
    wCalBuBs_->SetParameters(qs, qu, ms, mu);

    delete f;
    delete fBu;
    delete fBs;
    return true;  
}

bool OSMuonMvaTag::makeOsMuonTagging(){
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
    selectOsMuon();
    if(osMuonIndex_ < 0){ osMuonTagDecision_ = 0; return 1;}
    else osMuonTagDecision_ = -1*trkCharge->at(osMuonTrackIndex_); 

    computeOsMuonTagVariables();
    osMuonTagMvaValue_ = osMuonTagReader_.EvaluateMVA(methodName_);
    osMuonTagMistagProbRaw_ = 1 - osMuonTagMvaValue_;
    osMuonTagMistagProbCalProcess_ = wCalProcess_->Eval(osMuonTagMistagProbRaw_);
    osMuonTagMistagProbCalProcessBuBs_ = wCalBuBs_->Eval(osMuonTagMistagProbCalProcess_);

    return 1;
}

int OSMuonMvaTag::selectOsMuon(){
    int iB = ssIndex_;
    int iPV = pvIndex_;

    vector <int> tkSsB = tracksFromSV(iB);
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
    if(pvIndex_ < 0) pvIndex_ = GetBestPV(iB, tB);

    int bestMuIndex = -1;
    float bestMuPt = 2.;
    int bestMuTrack = -1;
    nMuonsSel_ = 0;

    for(int iMuon = 0; iMuon < nMuons; ++iMuon ){

        int itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu<0) continue;

        if(std::find(tkSsB.begin(), tkSsB.end(), itkmu) != tkSsB.end()) continue;

        if(muoPt->at( iMuon ) < 2.) continue;
        if(fabs(muoEta->at( iMuon )) > 2.4) continue;
        if(!IsMvaMuon(iMuon, wp_)) continue;
        if(fabs(dZ(itkmu, iPV)) > dzCut_) continue;
        if(deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon)) < 0.4) continue;
        //if(GetMuoPFiso(iMuon) > PFIsoCut_)  continue;

        nMuonsSel_++;

        if(muoPt->at( iMuon ) > bestMuPt){
            bestMuPt = muoPt->at( iMuon );
            bestMuIndex = iMuon;
            bestMuTrack = itkmu;
        }
    }

    osMuonIndex_ = bestMuIndex;
    osMuonTrackIndex_ = bestMuTrack;

    return bestMuIndex;
}

void OSMuonMvaTag::computeOsMuonTagVariables()
{
    int iB = ssIndex_;
    int iMuon = osMuonIndex_;
    int iPV = pvIndex_;

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    TLorentzVector tB = GetTLorentzVecFromJpsiX(iB);
    vector <int> tkSsB = tracksFromSV(iB);

    TLorentzVector tConeClean(0.,0.,0.,0.), tMu;
    tMu.SetPtEtaPhiM(muoPt->at(iMuon),muoEta->at(iMuon),muoPhi->at(iMuon),MassMu);

    float kappa = 1;
    float drCone = 0.4;

    float qConeClean = 0., ptConeClean = 0.;
    float muoConeCleanNF = 0., muoConeCleanCF = 0;
    int muoConeCleanNCH = 0;

    for(int ipf=0; ipf<nPF; ++ipf){
        float ptpfc = pfcPt->at(ipf);
        float etapfc = pfcEta->at(ipf);
        if( deltaR(etapfc, pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > drCone) continue;
        if(ptpfc < 0.5) continue;
        if(fabs(etapfc) > 3.0) continue;
        if(pfcCharge->at(ipf) == 0) continue;
        if(std::find(tkSsB.begin(), tkSsB.end(), pfcTrk->at(ipf)) != tkSsB.end()) continue;
        if(pfcTrk->at(ipf)<0) continue;
        if(fabs(dZ(pfcTrk->at(ipf), iPV))>=1.0) continue;
  
        TLorentzVector a;
        a.SetPxPyPzE(pfcPx->at(ipf), pfcPy->at(ipf), pfcPz->at(ipf), pfcE->at(ipf));
        tConeClean += a;

        qConeClean += pfcCharge->at(ipf) * pow(ptpfc, kappa);
        ptConeClean += pow(ptpfc, kappa);

        if(pfcCharge->at(ipf)==0) muoConeCleanNF += pfcE->at(ipf);
        if(abs(pfcCharge->at(ipf))==1){
            muoConeCleanNCH++;
            muoConeCleanCF += pfcE->at(ipf);
        }
    }

    if(ptConeClean != 0) qConeClean /= ptConeClean;
    else qConeClean = 1;
    qConeClean *= trkCharge->at(itkmu);
    if(tConeClean.E()!=0){
        muoConeCleanCF /= tConeClean.E();
        muoConeCleanNF /= tConeClean.E();
    }

    muoConeCleanQ_ = qConeClean;
    muoConeCleanPt_ = tConeClean.Pt();
    muoConeCleanDr_ = deltaR(tConeClean.Eta(), tConeClean.Phi(), muoEta->at( iMuon ), muoPhi->at(iMuon));
    if(tConeClean.E()!=0) muoConeCleanEnergyRatio_ = muoE->at(iMuon) / tConeClean.E();
    else muoConeCleanEnergyRatio_ = 1;

    tConeClean -= tMu;
    muoConeCleanPtRel_ = muoPt->at( iMuon ) * (tMu.Vect().Unit() * tConeClean.Vect().Unit());
    tConeClean += tMu; // for IP sign

    muoPt_ = muoPt->at( iMuon );
    muoEta_ = muoEta->at( iMuon );
    muoCharge_ = trkCharge->at(itkmu);
    muoDxy_ = dSign(itkmu, tConeClean.Px(), tConeClean.Py())*abs(trkDxy->at(itkmu));
    muoExy_ = trkExy->at(itkmu);
    muoDz_ = dZ(itkmu, iPV);
    muoEz_ = trkEz->at(itkmu);
    muoSoftMvaValue_ = computeMuonMva(iMuon);
    muoDrB_ = deltaR(tB.Eta(), tB.Phi(), muoEta->at(iMuon), muoPhi->at(iMuon));
    muoPFIso_ = GetMuoPFiso(iMuon);
    
}
