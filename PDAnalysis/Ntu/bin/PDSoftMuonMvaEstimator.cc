//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		INSTRUCTIONS
//
//		---Preparations---
//
//		-Include PDMuonVar.cc and PDSoftMuonMvaEstimator.cc in your PDAnalyzer.cc
//		-Include PDMuonVar.h and PDSoftMuonMvaEstimator.h in your PDAnalyzer.h
//		-Add PDMuonVar and PDSoftMuonMvaEstimator as public virtual classes to class PDAnalyzer in PDAnalyzer.h
//
//		---Definitions---

//		BARREL <-> abs(muoEta)<1.2
//		ENDCAP <-> abs(muoEta)>=1.2
//		BDTs trained with global muons with pT>2 GeV and abs(eta)<2.4 and basic quality cuts
//
//
//		---How to use the discriminator (example for barrel, endcap is the same)---
//
//		0. Change isoFlag_ and ipFlag_ with the corresponding function if you want to veto these variables
//		1. Initialize the discriminator in PDAnalyzer::beginJob with 'void PDSoftMuonMvaEstimator::setupReaderBarrel(std::string weightFileBarrel)'
//			-You can find the weights in /lustre/cmswork/abragagn/weights/
//		2. Compute the needed muon variables for each event with 'void computeMuonVar()' (in PDAnalyzer::analyze)
//		3. Compute the Mva response with 'float PDSoftMuonMvaEstimator::computeMvaBarrel(int iMuon)'
//			-Make sure to feed muons with the right pseudorapidity range otherwise you'll get invalid response
//		4. For endcap is the same, change every name accordingly
//			-a barrel ID and a endcap ID can be initialize at the same time
//			4.1 The function 'float AlbertoUtil::getMvaMuonValue(int iMuon)' takes care automatically of the pseudorapidity range
//
//
//
//
//		---Possible output values---

//		-1 = the muon is not a global muon
//		-2 = the muon do not pass basic preselection
//		-3 = the muon is not in the pseudorapidity range of the discriminator
//		[0, 1] = mva discriminator response
//
//
//		Author: Alberto Bragagnolo (alberto.bragagnolo@pd.infn.it)
//		Based on the BMM4 soft muon ID developed by Stephan Wiederkehr (wistepha@phys.ethz.ch)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"

PDSoftMuonMvaEstimator::PDSoftMuonMvaEstimator():
	readerBarrel("!Color:Silent"),
	readerEndcap("!Color:Silent")
{
	TMVA::PyMethodBase::PyInitialize();
}

PDSoftMuonMvaEstimator::~PDSoftMuonMvaEstimator() {}

// =====================================================================================
void PDSoftMuonMvaEstimator::setupReader(TMVA::Reader &reader, TString weightFile, TString methodName)
{

	reader.AddVariable( "muoPt", &muoPt_ );
	reader.AddVariable( "abs(muoEta)", &absMuoEta_ );
	reader.AddVariable( "muoSegmComp", &muoSegmComp_ );
	reader.AddVariable( "muoChi2LM", &muoChi2LM_ );
	reader.AddVariable( "muoChi2LP", &muoChi2LP_ );
	reader.AddVariable( "muoGlbTrackTailProb", &muoGlbTrackTailProb_ );
	reader.AddVariable( "muoIValFrac", &muoIValFrac_ );
	reader.AddVariable( "muoLWH", &muoLWH_ );
	reader.AddVariable( "muoTrkKink", &muoTrkKink_ );
	reader.AddVariable( "muoGlbKinkFinderLOG", &muoGlbKinkFinderLOG_ );
	reader.AddVariable( "muoTimeAtIpInOutErr", &muoTimeAtIpInOutErr_ );
	reader.AddVariable( "muoOuterChi2", &muoOuterChi2_ );
	reader.AddVariable( "muoInnerChi2", &muoInnerChi2_ );  
	reader.AddVariable( "muoTrkRelChi2", &muoTrkRelChi2_ );
	reader.AddVariable( "muoVMuonHitComb", &muoVMuonHitComb_ );
	reader.AddVariable( "muoGlbDeltaEtaPhi", &muoGlbDeltaEtaPhi_ );
	reader.AddVariable( "muoStaRelChi2", &muoStaRelChi2_ );
	reader.AddVariable( "muoTimeAtIpInOut", &muoTimeAtIpInOut_ );
	reader.AddVariable( "muoValPixHits", &muoValPixHits_ );
	reader.AddVariable( "muoNTrkVHits", &muoNTrkVHits_ );
	reader.AddVariable( "muoGNchi2", &muoGNchi2_ );
	reader.AddVariable( "muoVMuHits", &muoVMuHits_ );
	reader.AddVariable( "muoNumMatches", &muoNumMatches_ );
	if(ipFlag_) reader.AddVariable( "trkDxy/trkExy", &trkDxy_ );
	if(ipFlag_) reader.AddVariable( "trkDz/trkEz", &trkDz_ );	
	if(isoFlag_) reader.AddVariable( "muoPFiso", &muoPFiso_ );
	reader.AddVariable( "muoQprod", &muoQprod_ );

	reader.AddSpectator( "muoEvt", &DUMMY_ );

	reader.BookMVA( methodName, weightFile );

	return;

}

// =====================================================================================
void PDSoftMuonMvaEstimator::setupReaderBarrel(TString weightFileBarrel, TString methodNameBarrel)
{

	methodNameBarrel_ = methodNameBarrel;
	setupReader(readerBarrel, weightFileBarrel, methodNameBarrel);

	return;

}

// =====================================================================================
void PDSoftMuonMvaEstimator::setupReaderEndcap(TString weightFileEndcap, TString methodNameEndcap)
{

	methodNameEndcap_ = methodNameEndcap;
	setupReader(readerEndcap, weightFileEndcap, methodNameEndcap);

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
	muoInnerChi2_ =  muoInnerChi2->at(iMuon);
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
float PDSoftMuonMvaEstimator::computeMvaBarrel(int iMuon)
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

	if( abs(muoEta->at(iMuon))>=1.2)
	{
		return -3;
	}

	//VARIABLE EXTRACTION
	computeMvaVariables(iMuon);

  //PRESELECTION
	if(!MuonPassedPreselection(iMuon))
	{
		return -2;
	}

	return readerBarrel.EvaluateMVA(methodNameBarrel_);

}

// =====================================================================================
float PDSoftMuonMvaEstimator::computeMvaEndcap(int iMuon)
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

	if( abs(muoEta->at(iMuon))<1.2)
	{
		return -3;
	}

	//VARIABLE EXTRACTION
	computeMvaVariables(iMuon);

	//PRESELECTION
	if(!MuonPassedPreselection(iMuon))
	{
		return -2;
	}

	return readerEndcap.EvaluateMVA(methodNameEndcap_);

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
    IPsign = dSign(itkmu, iJet);
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
void PDSoftMuonMvaEstimator::ipFlag(bool newValue)
{
	ipFlag_ = newValue;
	return;
}

// =====================================================================================
void PDSoftMuonMvaEstimator::isoFlag(bool newValue)
{
	isoFlag_ = newValue;
	return;
}
