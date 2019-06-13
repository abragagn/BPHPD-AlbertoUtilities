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
//
//      Methods trained with global muons with pT>2 GeV and abs(eta)<2.4 and basic quality cuts
//
//      ---How to use the discriminator ---
//
//      0. You can find the weights in /lustre/cmswork/abragagn/mvaWeights/MvaMuonID/
//      1. Initialize the discriminator in PDAnalyzer::beginJob with 'void PDSoftMuonMvaEstimator::inizializeMuonMvaReader()'
//      2. In PDAnalyzer::analyze compute the needed muon variables for each event with 'void computeMuonVar() 
//          and fill the Cartesian coordinates vectors of muons, tracks, jet and pfcs
//          e.g convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
//      3. Compute the Mva response with 'float PDSoftMuonMvaEstimator::computeMuonMva(int iMuon)'
//
//
//      ---Possible output values---
//
//      -1 = the muon is not a global muon
//      -2 = the muon do not pass preselection
//      [0, 1] = mva discriminator response
//
//
//      ---Preselection efficienty for particle with pT>2 GeV and abs(eta)<2.4---
//
//          mu           --> %
//          K, pi, p, e  --> %
//
//
//      Author: Alberto Bragagnolo (alberto.bragagnolo@cern.ch)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

PDSoftMuonMvaEstimator::PDSoftMuonMvaEstimator():
    muonMvaIdReader_("!Color:Silent")
{
    TMVA::PyMethodBase::PyInitialize();
}

PDSoftMuonMvaEstimator::~PDSoftMuonMvaEstimator() {}

// =====================================================================================
void PDSoftMuonMvaEstimator::inizializeMuonMvaReader(TString methodName = "DNNMuonID", TString path = "/lustre/cmswork/abragagn/mvaWeights/MvaMuonID/")
{
    methodSetup(methodName, path);

    muonMvaIdReader_.AddVariable( "muoPt", &muoPt_ );
    muonMvaIdReader_.AddVariable( "muoEta", &muoEta_ );
    muonMvaIdReader_.AddVariable( "muoSegmComp", &muoSegmComp_ );
    muonMvaIdReader_.AddVariable( "muoChi2LM", &muoChi2LM_ );
    muonMvaIdReader_.AddVariable( "muoChi2LP", &muoChi2LP_ );
    muonMvaIdReader_.AddVariable( "muoGlbTrackTailProb", &muoGlbTrackTailProb_ );
    muonMvaIdReader_.AddVariable( "muoIValFrac", &muoIValFrac_ );
    muonMvaIdReader_.AddVariable( "muoLWH", &muoLWH_ );
    muonMvaIdReader_.AddVariable( "muoTrkKink", &muoTrkKink_ );
    muonMvaIdReader_.AddVariable( "muoGlbKinkFinderLOG", &muoGlbKinkFinderLOG_ );
    muonMvaIdReader_.AddVariable( "muoTimeAtIpInOutErr", &muoTimeAtIpInOutErr_ );
    muonMvaIdReader_.AddVariable( "muoOuterChi2", &muoOuterChi2_ );
    muonMvaIdReader_.AddVariable( "muoInnerChi2", &muoInnerChi2_ );
    muonMvaIdReader_.AddVariable( "muoTrkRelChi2", &muoTrkRelChi2_ );
    muonMvaIdReader_.AddVariable( "muoVMuonHitComb", &muoVMuonHitComb_ );
    muonMvaIdReader_.AddVariable( "muoGlbDeltaEtaPhi", &muoGlbDeltaEtaPhi_ );
    muonMvaIdReader_.AddVariable( "muoStaRelChi2", &muoStaRelChi2_ );
    muonMvaIdReader_.AddVariable( "muoTimeAtIpInOut", &muoTimeAtIpInOut_ );
    muonMvaIdReader_.AddVariable( "muoValPixHits", &muoValPixHits_ );
    muonMvaIdReader_.AddVariable( "muoNTrkVHits", &muoNTrkVHits_ );
    muonMvaIdReader_.AddVariable( "muoGNchi2", &muoGNchi2_ );
    muonMvaIdReader_.AddVariable( "muoVMuHits", &muoVMuHits_ );
    muonMvaIdReader_.AddVariable( "muoNumMatches", &muoNumMatches_ );
    muonMvaIdReader_.AddVariable( "muoQprod", &muoQprod_ );
    if(useIso(methodName_)) muonMvaIdReader_.AddVariable( "muoPFiso", &muoPFiso_ );
    muonMvaIdReader_.AddSpectator( "muoEvt", &DUMMY_ );

    muonMvaIdReader_.BookMVA( methodName_, weightFile_ );

    return;
}


// =====================================================================================
void PDSoftMuonMvaEstimator::computeMvaVariables(int iMuon)
{
    DUMMY_ = -1;

    muoPt_ = muoPt->at(iMuon);
    muoEta_ = muoEta->at(iMuon);

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

    float PFIso = muoSumCPpt->at(iMuon)/muoPt->at(iMuon);
    float betaCorr = muoSumNHet->at(iMuon) + muoSumPHet->at(iMuon)-0.5*(muoSumPUpt->at(iMuon));
    betaCorr/=muoPt->at(iMuon);
    if(betaCorr>0) PFIso+=betaCorr;

    muoPFiso_ = PFIso;
    muoQprod_ = muoQprod->at(iMuon);

    return;
}


// =====================================================================================
float PDSoftMuonMvaEstimator::computeMuonMva(int iMuon)
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
    if(!muonPassedPreselection(iMuon))
    {
        return -2;
    }

    return muonMvaIdReader_.EvaluateMVA(methodName_);
}

// =====================================================================================
bool PDSoftMuonMvaEstimator::muonPassedPreselection(int iMuon)
{
    if ( muoChi2LM->at( iMuon ) > 5000 ) return false;
    if ( muoChi2LP->at( iMuon ) > 2000 ) return false;
    if ( muoGlbTrackTailProb->at( iMuon ) > 5000 ) return false;
    if ( muoTrkKink->at( iMuon ) > 900 ) return false;
    if ( muoGlbKinkFinderLOG->at( iMuon ) > 50 ) return false;
    if ( muoTimeAtIpInOutErr->at( iMuon ) > 4 ) return false;
    if ( muoOuterChi2->at( iMuon ) > 1000 ) return false;
    if ( muoInnerChi2->at( iMuon ) > 10 ) return false;
    if ( muoTrkRelChi2->at( iMuon ) > 3 ) return false;

    return true;
}

// =====================================================================================
bool PDSoftMuonMvaEstimator::useIso(TString methodName)
{
    return !methodName.Contains("woIso");
}

// =====================================================================================
void PDSoftMuonMvaEstimator::methodSetup(TString methodName, TString path)
{
    weightFile_ = path + "TMVAClassification_" + methodName + ".weights.xml";
    methodName_ = methodName;
    return;
}
