#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"

// additional features
#include "PDSecondNtupleWriter.h"
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

using namespace std;

/* EXAMPLE OF USAGE
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -v useHLT t -v writeVars f -n 2000000
*/
PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "false" );
    setUserParameter( "writeVars", "true" );

    setUserParameter( "outputFile", "ntu.root" );

    setUserParameter( "muonIdWpBarrel", "0.8910" ); 
    setUserParameter( "muonIdWpEndcap", "0.8925" ); 

    setUserParameter( "muonMvaMethod",      "BDTMuonID2017woIPwIso" ); 
    setUserParameter( "osMuonTagMvaMethod", "DNNOsMuonHLTJpsiMu_test241" ); 

    setUserParameter( "ptCut", "40.0" ); //needed for paolo's code (no influence in the code whatsoever)

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

    PDAnalyzerUtil::beginJob();

    // user parameters are retrieved as strings by using their names;
    // numeric parameters ( int, float or whatever ) can be directly set
    // by passing the corresponding variable,
    // e.g. getUserParameter( "name", x )

    getUserParameter( "verbose", verbose );

    getUserParameter( "process", process );
    getUserParameter( "useHLT", useHLT );
    getUserParameter( "writeVars", writeVars );

    getUserParameter( "outputFile", outputFile );

    getUserParameter( "muonIdWpBarrel", muonIdWpBarrel ); 
    getUserParameter( "muonIdWpEndcap", muonIdWpEndcap ); 

    getUserParameter( "muonMvaMethod", muonMvaMethod );
    getUserParameter( "osMuonTagMvaMethod", osMuonTagMvaMethod );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code (no influence in the code whatsoever)

//  additional features
    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    setOsMuonCuts(muonIdWpBarrel, muonIdWpEndcap, 1. ); //set wp for muonID and Dz cut

    inizializeMuonMvaReader( muonMvaMethod );           //initialize mva muon id
    inizializeOSMuonMvaTagReader( osMuonTagMvaMethod ); //inizialize mva os muon method
    bool osInit = inizializeOSMuonMvaMistagMethods();   //read os muon method for per-event-mistag
    if(!osInit) cout<<"METHOD NOT INIZIALIZATED. ABORT"<<endl;

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.50);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.50);

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility

    float min = 5.0;
    float max = 5.5;
    float nbin = 250;


    autoSavedObject =
    hmass_ssB       = new TH1D( "hmass_ssB", "hmass_ssB", nbin, min, max );

    return;

}


void PDAnalyzer::reset() {

    autoReset();
    return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

    if ( verbose ) {
        cout << " +++++++++++++++++++++++++++ " << endl;
        cout << "entry: "
             << entry << " " << event_file << " " << event_tot << endl;
        cout << "run: " << runNumber << " , "
             << "evt: " << eventNumber << endl;
    }
    else {
        if ( (!(event_tot%10) && event_tot<100 ) || 
        (!(event_tot %100) && event_tot<1000 ) || 
        (!(event_tot %1000)&& event_tot<10000 ) || 
        (!(event_tot %10000) && event_tot<100000 ) || 
        (!(event_tot %100000) && event_tot<1000000 ) || 
        (!(event_tot %1000000) && event_tot<10000000 ) )
            cout << " == at event " << event_file << " " << event_tot << endl;
    }

// additional features
    computeMuonVar();   //compute muon variable for soft id
    inizializeTagVariables(); //initialiaze some variable for tagging
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    if( !((process=="BsJPsiPhi")||(process=="BuJPsiK")) ) {
        cout<<"!$!#$@$% PROCESS NAME WRONG"<<endl;
        return false;
    }

//------------------------------------------------HLT---------------------------------------

    bool jpsimu = false;
    bool jpsitktk = false;
    bool jpsitk = false;

    if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;

    if( !jpsimu ) return false; //currently only jpsimu is allowed
    if( jpsimu ) SetJpsiMuCut();    //set analysis cut for jpsimu
    if( !jpsimu ) SetJpsiTrkTrkCut();   //set analysis cut for jpsitktk

    if(useHLT && process=="BsJPsiPhi" && !(jpsimu || jpsitktk)) return false;
    if(useHLT && process=="BuJPsiK" && !(jpsimu || jpsitk)) return false;

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int iSsB = GetCandidate(process); //get best svt
    if(iSsB<0) return false;

    bool isTight = false;
    int iSsBtight = GetTightCandidate(process); //get best svt that pass analysis selection
    if(iSsBtight>=0){
        isTight = true;
        iSsB = iSsBtight;
    }else return false; //apply analysis selection to the event

    TLorentzVector tB = GetTLorentzVecFromJpsiX(iSsB);
    int iSsPV = GetBestPV(iSsB, tB); //select PV
    if(iSsPV < 0) return false;

    setSsForTag(iSsB, iSsPV); //set PV and SVT for tagging class

    hmass_ssB->Fill(svtMass->at(iSsB));
    (tWriter->ssbMass) = svtMass->at(iSsB);
    
//-----------------------------------------TAG-----------------------------------------

    int bestMuIndex = getOsMuon(); //get OS muon
    int tagDecision = getOsMuonTag();   //get tag decision 

    if( tagDecision == 0 ){
        cout<<"no os muon founded"<<endl;
        (tWriter->osMuonTag) = 0;
        return true;
    }

    pair<float,float> osMuonTagMistag = getOsMuonTagMistagProb(2); //.first =  mistag, .second =  error (still not implemented). Argument define what method to use

     (tWriter->osMuonTag) = tagDecision;
     (tWriter->osMuonTagMistag) = osMuonTagMistag.first;

    cout<<"os muon "<<bestMuIndex<<" founded with decision "<<tagDecision<<" and mistag "<<osMuonTagMistag.first<<endl;

    return true;

}

void PDAnalyzer::endJob() {

// additional features
    tWriter->close();   // second ntuple
    return;
}


void PDAnalyzer::save() {
#   if UTIL_USE == FULL
    // explicit saving not necessary for "autoSavedObjects"
    autoSave();
#elif UTIL_USE == BARE
    // explicit save histos when not using the full utility

#endif

    return;
}
