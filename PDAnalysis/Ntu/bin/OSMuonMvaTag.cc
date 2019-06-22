#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

using namespace std;

OSMuonMvaTag::OSMuonMvaTag():
                osMuonTagReader_("!Color:Silent")
,               ssIndex_(-1)
,               pvIndex_(-1)
,               osMuonIndex_(-1)
,               osMuonTrackIndex_(-1)
,               osMuonTagMvaValue_(-1)
,               wp_(0.21)
,               dzCut_(1.)
,               nMuonsSel_(0)
{}

OSMuonMvaTag::~OSMuonMvaTag() {}

// =====================================================================================
void OSMuonMvaTag::inizializeTagVariables()
{
    ssIndex_ = -1;
    pvIndex_ = -1;
    osMuonIndex_ = -1;
    osMuonTrackIndex_ = -1;
    osMuonTagMvaValue_ = -1;
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
    osMuonTagReader_.BookMVA( methodName_, methodPath + methodName_ + ".weights.xml" );
}

bool OSMuonMvaTag::inizializeOSMuonCalibration( 
    TString process = "BuJPsiKData2018"
,   TString methodPath = ""  
)
{
    if(methodPath == "") methodPath = methodPath_;
    auto *f = new TFile(methodPath + "OSMuonTaggerCalibration" + process + ".root");
    if(f->IsZombie()) return false;
    f->cd();
    wCal_ = (TF1*)f->Get("osMuonCal");
    f->Close();
    delete f;
    return true;  
}

int OSMuonMvaTag::getOsMuon()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -2; }

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

void OSMuonMvaTag::computeVariables()
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

int OSMuonMvaTag::getOsMuonTag()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
    if(osMuonIndex_ < 0) getOsMuon();
    if(osMuonIndex_ < 0) return 0;
    return -1*trkCharge->at(osMuonTrackIndex_); 
}

float OSMuonMvaTag::getOsMuonTagMvaValue()
{
    if(ssIndex_ < 0){ cout<<"SS NOT INITIALIZED"<<endl; return -999; }
    if(osMuonIndex_ < 0){ cout<<"WARNING: OS MU NOT INITIALIZED"<<endl; getOsMuon(); }
    
    computeVariables();
    osMuonTagMvaValue_ = osMuonTagReader_.EvaluateMVA(methodName_);
    return osMuonTagMvaValue_;
}

pair<float, float> OSMuonMvaTag::getOsMuonTagMistagProb()
{
    float dnnValue = osMuonTagMvaValue_;
    if(dnnValue == -1) dnnValue = getOsMuonTagMvaValue();

    float wPred = 1. - dnnValue;
    float evtMistagProb = wCal_->Eval(wPred);
    float evtMistagProbError = sqrt(pow(wCal_->GetParError(0),2)+pow((wCal_->GetParError(1))*(wPred),2));

    if(evtMistagProb > 1.) evtMistagProb = 1.;
    if(evtMistagProb < 0.) evtMistagProb = 0.;

    return make_pair(evtMistagProb, evtMistagProbError);
}
