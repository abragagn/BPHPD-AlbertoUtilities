//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      INSTRUCTIONS
//
//      ---Preparations---
//
//      -Include PDMuonVar.cc and PDSoftMuonMvaEstimator.cc in your PDAnalyzer.cc
//      -Include PDMuonVar.h and PDSoftMuonMvaEstimator.h in your PDAnalyzer.h
//      -Add PDMuonVar and PDSoftMuonMvaEstimator as public virtual classes to class PDAnalyzer in PDAnalyzer.h
//
//      ---Definitions---

//      BARREL <-> abs(muoEta)<1.2
//      ENDCAP <-> abs(muoEta)>=1.2
//      BDTs trained with global muons with pT>2 GeV and abs(eta)<2.4 and basic quality cuts
//
//
//      ---How to use the discriminator ---
//
//      0. You can find the weights in /lustre/cmswork/abragagn/weights/
//      1. Initialize the discriminator in PDAnalyzer::beginJob with 'void PDSoftMuonMvaEstimator::setupMuonMvaReader(TString methodName)'
//          -methodName should be in the form of prefix + year + variable flags (w = with, wo = without. The order is "IP" followed by "Iso")
//              --e.g "DNNGlobal2016woIPwIso", but even "DNNGlobalBarrel2016woIPwIso" is accepted
//      2. In PDAnalyzer::analyze compute the needed muon variables for each event with 'void computeMuonVar() 
//          and fill the Cartesian coordinates vectors of muons, tracks, jet and pfcs
//          e.g convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
//      3. Compute the Mva response with 'float PDSoftMuonMvaEstimator::computeMva(int iMuon)'
//          -- Remember that Barrel and Encap use completely different methods
//
//
//
//
//      ---Possible output values---
//
//      -1 = the muon is not a global muon
//      -2 = the muon do not pass preselection
//      [0, 1] = mva discriminator response
//
//
//      Author: Alberto Bragagnolo (alberto.bragagnolo@pd.infn.it)
//      Based on the BMM4 soft muon ID developed by Stephan Wiederkehr (wistepha@phys.ethz.ch)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

PDSoftMuonMvaEstimator::PDSoftMuonMvaEstimator():
    reader_("!Color:Silent")
{
    TMVA::PyMethodBase::PyInitialize();
}

PDSoftMuonMvaEstimator::~PDSoftMuonMvaEstimator() {}

// =====================================================================================
void PDSoftMuonMvaEstimator::setupMuonMvaReader(TString methodName, TString path = "/lustre/cmswork/abragagn/weights/")
{

    methodSetup(methodName, path);

    reader_.AddVariable( "muoPt", &muoPt_ );
    reader_.AddVariable( "abs(muoEta)", &absMuoEta_ );
    reader_.AddVariable( "muoSegmComp", &muoSegmComp_ );
    reader_.AddVariable( "muoChi2LM", &muoChi2LM_ );
    reader_.AddVariable( "muoChi2LP", &muoChi2LP_ );
    reader_.AddVariable( "muoGlbTrackTailProb", &muoGlbTrackTailProb_ );
    reader_.AddVariable( "muoIValFrac", &muoIValFrac_ );
    reader_.AddVariable( "muoLWH", &muoLWH_ );
    reader_.AddVariable( "muoTrkKink", &muoTrkKink_ );
    reader_.AddVariable( "muoGlbKinkFinderLOG", &muoGlbKinkFinderLOG_ );
    reader_.AddVariable( "muoTimeAtIpInOutErr", &muoTimeAtIpInOutErr_ );
    reader_.AddVariable( "muoOuterChi2", &muoOuterChi2_ );
    reader_.AddVariable( "muoInnerChi2", &muoInnerChi2_ );
    reader_.AddVariable( "muoTrkRelChi2", &muoTrkRelChi2_ );
    reader_.AddVariable( "muoVMuonHitComb", &muoVMuonHitComb_ );
    reader_.AddVariable( "muoGlbDeltaEtaPhi", &muoGlbDeltaEtaPhi_ );
    reader_.AddVariable( "muoStaRelChi2", &muoStaRelChi2_ );
    reader_.AddVariable( "muoTimeAtIpInOut", &muoTimeAtIpInOut_ );
    reader_.AddVariable( "muoValPixHits", &muoValPixHits_ );
    reader_.AddVariable( "muoNTrkVHits", &muoNTrkVHits_ );
    reader_.AddVariable( "muoGNchi2", &muoGNchi2_ );
    reader_.AddVariable( "muoVMuHits", &muoVMuHits_ );
    reader_.AddVariable( "muoNumMatches", &muoNumMatches_ );
    reader_.AddVariable( "muoQprod", &muoQprod_ );
    if(useIp(methodNameBarrel_)){
        reader_.AddVariable( "trkDxy/trkExy", &trkDxy_ );
        reader_.AddVariable( "trkDz/trkEz", &trkDz_ );
    }
    if(useIso(methodNameBarrel_)) reader_.AddVariable( "muoPFiso", &muoPFiso_ );

    reader_.AddSpectator( "muoEvt", &DUMMY_ );

    reader_.BookMVA( methodNameBarrel_, weightFileBarrel_ );
    reader_.BookMVA( methodNameEndcap_, weightFileEndcap_ );

    return;

}


// =====================================================================================
void PDSoftMuonMvaEstimator::computeMvaVariables(int iMuon){

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );

    DUMMY_ = -1;

    muoPt_ = muoPt->at(iMuon);
    absMuoEta_ = abs(muoEta->at(iMuon));

    muoSegmComp_ = muoSegmComp->at(iMuon);
    muoChi2LM_ = muoChi2LM->at(iMuon);
    muoChi2LP_ = muoChi2LP->at(iMuon);
    muoGlbTrackTailProb_ = muoGlbTrackTailProb->at(iMuon);
    muoIValFrac_ = muoIValFrac->at(iMuon);
    muoLWH_ = muoLWH->at(iMuon);
    muoTrkKink_ = muoTrkKink->at(iMuon);
    muoGlbKinkFinderLOG_ = muoGlbKinkFinderLOG->at(iMuon);
    muoTimeAtIpInOutErr_ = muoTimeAtIpInOutErr->at(iMuon);
    muoOuterChi2_ = muoOuterChi2->at(iMuon);
    muoInnerChi2_ = muoInnerChi2->at(iMuon);
    muoTrkRelChi2_ = muoTrkRelChi2->at(iMuon);
    muoVMuonHitComb_ = muoVMuonHitComb->at(iMuon);

    muoGlbDeltaEtaPhi_ = muoGlbDeltaEtaPhi->at(iMuon);
    muoStaRelChi2_ = muoStaRelChi2->at(iMuon);
    muoTimeAtIpInOut_ = muoTimeAtIpInOut->at(iMuon);
    muoValPixHits_ = muoValPixHits->at(iMuon);
    muoNTrkVHits_ = muoNTrkVHits->at(iMuon);
    muoGNchi2_ = muoGNchi2->at(iMuon);
    muoVMuHits_ = muoVMuHits->at(iMuon);
    muoNumMatches_ = muoNumMatches->at(iMuon);

    trkDxy_ = (abs(trkDxy->at(itkmu)) * IPsign_(iMuon)) / trkExy->at(itkmu);
    trkDz_ = trkDz->at(itkmu) / trkEz->at(itkmu);

    float PFIso = muoSumCPpt->at(iMuon)/muoPt->at(iMuon);
    float betaCorr = muoSumNHet->at(iMuon) + muoSumPHet->at(iMuon)-0.5*(muoSumPUpt->at(iMuon));
    betaCorr/=muoPt->at(iMuon);
    if(betaCorr>0) PFIso+=betaCorr;

    muoPFiso_ = PFIso;

    muoQprod_ = muoQprod->at(iMuon);

    return;

}


// =====================================================================================
float PDSoftMuonMvaEstimator::computeMva(int iMuon)
{   

    if( !( muoType->at(iMuon) & PDEnumString::global ) )
    {
        return -1;
    }

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );

    if( itkmu < 0 )
    {
        return -1;
    }   

    if( !(( trkQuality->at( itkmu ) >> 2 ) & 1) )
    {
        return -2;
    }

    //VARIABLE EXTRACTION
    computeMvaVariables(iMuon);
    //PRESELECTION
    if(!MuonPassedPreselection(iMuon))
    {
        return -2;
    }

    return (abs(muoEta->at( iMuon ))<1.2) ? reader_.EvaluateMVA(methodNameBarrel_) : reader_.EvaluateMVA(methodNameEndcap_); 

}

// =====================================================================================
bool PDSoftMuonMvaEstimator::MuonPassedPreselection(int iMuon)
{

    if ( muoChi2LM->at( iMuon ) > 5000 ) {return false;}
    if ( muoChi2LP->at( iMuon ) > 2000 ) {return false;}
    if ( muoGlbTrackTailProb->at( iMuon ) > 5000 ) {return false;}
    if ( muoTrkKink->at( iMuon ) > 900 ) {return false;}
    if ( muoGlbKinkFinderLOG->at( iMuon ) > 50 ) {return false;}
    if ( muoTimeAtIpInOutErr->at( iMuon ) > 4 ) {return false;}
    if ( muoOuterChi2->at( iMuon ) > 1000 ) {return false;}
    if ( muoInnerChi2->at( iMuon ) > 10 ) {return false;}
    if ( muoTrkRelChi2->at( iMuon ) > 3 ) {return false;}

    return true;
}

// =====================================================================================
int PDSoftMuonMvaEstimator::IPsign_(int iMuon)
{

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    int ipftkmu = trkPFC->at(itkmu);
    int IPsign = ((double)rand() / (RAND_MAX)) < 0.5 ? -1 : +1; //random value +-1

    int iJet = trkJet->at(itkmu);
    if(iJet<0 && ipftkmu>=0) iJet=pfcJet->at(ipftkmu);

    if(iJet>=0){
        IPsign = dSign(itkmu, jetPx->at( iJet ), jetPy->at( iJet ));
    }else{
        int coneNtrk = 0;
        float pxCone = 0, pyCone = 0;

        for(int ipf = 0; ipf<nPF; ++ipf){

            if( deltaR(pfcEta->at(ipf), pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > 0.4 ) continue;
            //if(std::find(signalTracks.begin(), signalTracks.end(), pfcTrk->at(ipf)) != signalTracks.end()) continue;
            if(pfcTrk->at(ipf) == itkmu) continue;
            if(pfcPt->at(ipf) < 0.2) continue;
            if(abs(pfcEta->at(ipf)) > 2.5) continue;
            //if( !(( trkQuality->at( itk ) >> 2 ) & 1) ) continue;
            ++coneNtrk;
            pxCone += pfcPt->at(ipf)*TMath::Cos(pfcPhi->at(ipf));
            pyCone += pfcPt->at(ipf)*TMath::Sin(pfcPhi->at(ipf));

        }

        if(coneNtrk>=2) IPsign = dSign(itkmu, pxCone, pyCone);
    }

    return IPsign;
}

// =====================================================================================
int PDSoftMuonMvaEstimator::IPsign_(int iMuon, int iPV)
{

    int itkmu = muonTrack( iMuon, PDEnumString::muInner );
    int ipftkmu = trkPFC->at(itkmu);
    int IPsign = ((double)rand() / (RAND_MAX)) < 0.5 ? -1 : +1; //random value +-1

    int iJet = trkJet->at(itkmu);
    if(iJet<0 && ipftkmu>=0) iJet=pfcJet->at(ipftkmu);

    if(iJet>=0){
        IPsign = dSign(itkmu, jetPx->at( iJet ), jetPy->at( iJet ), pvtX->at(iPV), pvtY->at(iPV));
    }else{
        int coneNtrk = 0;
        float pxCone = 0, pyCone = 0;

        for(int ipf = 0; ipf<nPF; ++ipf){

            if( deltaR(pfcEta->at(ipf), pfcPhi->at(ipf), muoEta->at(iMuon), muoPhi->at(iMuon)) > 0.4 ) continue;
            if(pfcPt->at(ipf) < 0.2) continue;
            if(abs(pfcEta->at(ipf)) > 2.5) continue;

            ++coneNtrk;
            pxCone += pfcPx->at(ipf);
            pyCone += pfcPy->at(ipf);

        }
        if(coneNtrk>2) IPsign = dSign(itkmu, pxCone, pyCone, pvtX->at(iPV), pvtY->at(iPV));
    }

    return IPsign;
}

// =====================================================================================
bool PDSoftMuonMvaEstimator::useIp(TString methodName)
{
    return !methodName.Contains("woIP");
}

// =====================================================================================
bool PDSoftMuonMvaEstimator::useIso(TString methodName)
{
    return !methodName.Contains("woIso");
}

// =====================================================================================
TString PDSoftMuonMvaEstimator::methodNameFromWeightName(TString weightsName)
{
    TString prefix = "TMVAClassification_";
    int start = weightsName.Index(prefix) + prefix.Length();
    int length = weightsName.Index(".weights") - start;
    TString name( weightsName(start, length) );
    return name;
}

// =====================================================================================
void PDSoftMuonMvaEstimator::methodSetup(TString methodName, TString path)
{

    TString year = "";
    TString var = "";

    if(methodName.Contains("2016")) year = "2016";
    if(methodName.Contains("2017")) year = "2017";
    if(methodName.Contains("2018")) year = "2018";

    if(useIp(methodName)) var += "wIP"; 
        else var += "woIP";
    if(useIso(methodName)) var += "wIso"; 
        else var += "woIso";

    TString name( methodName( 0, methodName.Index(year) ) );

    if(methodName.Contains("Barrel") || methodName.Contains("Endcap")) name.Remove(name.Length() - 6, 6);

    weightFileBarrel_ = path + year + "/" + "TMVAClassification_" + name + "Barrel" + year + var + ".weights.xml";
    weightFileEndcap_ = path + year + "/" + "TMVAClassification_" + name + "Endcap" + year + var + ".weights.xml";

    methodNameBarrel_ = methodNameFromWeightName(weightFileBarrel_);
    methodNameEndcap_ = methodNameFromWeightName(weightFileEndcap_);

    return;

}