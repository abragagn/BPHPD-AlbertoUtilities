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

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // default values can be set in the analyzer class contructor

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "outputFile", "ntu.root" );
    setUserParameter( "ptCut", "40.0" ); //needed for paolo's code for unknow reasons

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

    PDAnalyzerUtil::beginJob();

    // user parameters are retrieved as strings by using their names;
    // numeric parameters ( int, float or whatever ) can be directly set
    // by passing the corresponding variable,
    // e.g. getUserParameter( "name", x )

    getUserParameter( "process", process );
    getUserParameter( "outputFile", outputFile );
    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

//  additional features
//  tWriter = new PDSecondNtupleWriter; // second ntuple
//  tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    inizializeMuonMvaReader(); // initialize TMVA methods for muon ID
    inizializeOSMuonMvaReader(); // initialize TMVA methods for muon tagger
    bool osInit = inizializeOSMuonCalibration(); // initialize calibration methods for muon tagger
    if(!osInit) cout<<endl<<"!!! FAILED TO INIZIALIZED TAG CALIBRATION"<<endl<<endl;

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.65);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.65);

    return;

}


void PDAnalyzer::book() {


    float min = 5.15;
    float max = 5.65;
    float nbin = 250;


    autoSavedObject =
    hmass_ssB       = new TH1D( "hmass_ssB", "hmass_ssB", nbin, min, max );

    autoSavedObject =
    hmass_ssB_os    = new TH1D( "hmass_ssB_os", "hmass_ssB_os", nbin, min, max );

    return;
}


void PDAnalyzer::reset() {
    autoReset();
    return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

    if ( (!(event_tot%10) && event_tot<100 ) || 
    (!(event_tot %100) && event_tot<1000 ) || 
    (!(event_tot %1000)&& event_tot<10000 ) || 
    (!(event_tot %10000) && event_tot<100000 ) || 
    (!(event_tot %100000) && event_tot<1000000 ) || 
    (!(event_tot %1000000) && event_tot<10000000 ) )
        cout << " == at event " << event_file << " " << event_tot << endl;

// additional features
    computeMuonVar(); // compute variable needed for the muon ID
    inizializeTagVariables(); // initialize variables for muon tagger 
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz); // needed for the methods
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz); // needed for the methods
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz); // needed for the methods
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz); // needed for the methods

    if( !((process=="BsJPsiPhi")||(process=="BuJPsiK")) ) {
        cout<<"!$!#$@$% PROCESS NAME WRONG"<<endl;
        return false;
    }

//------------------------------------------------HLT---------------------------------------

    bool jpsimu = false;
    bool jpsitktk = false;

    if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;

    if( !jpsimu ) return false; // hlt veto
    SetJpsiMuCut(); //set selection for jpsimu

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int ssbSVT = GetCandidate(process);
    if(ssbSVT<0) return false;

    bool isTight = false;
    int ssbSVTtight = GetTightCandidate(process);
    if(ssbSVTtight>=0){
        isTight = true;
        ssbSVT = ssbSVTtight;
    }

    int iJPsi = (subVtxFromSV(ssbSVT)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(ssbSVT);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(ssbSVT);

    int ssbPVT = GetBestPV(ssbSVT, tB);
    if(ssbPVT < 0) return false;

    setVtxForTag(ssbSVT, ssbPVT); // set vertices index dor muon tagger

    hmass_ssB->Fill(svtMass->at(ssbSVT));
    
//-----------------------------------------OPPOSITE SIDE-----------------------------------------

    int bestMuIndex = getOsMuon(); // get OS muon index
    int tagDecision = getOsMuonTag(); get Tag decision // 1*trkCharge->at(osMuonTrackIndex_), 0 -> no muon.

    if( tagDecision == 0 ){
        return true;
    }

    hmass_ssB_os->Fill(svtMass->at(ssbSVT));

    float osMuonTagMvaValue = getOsMuonTagMvaValue();
    pair<float,float> osMuonTagMistag = getOsMuonTagMistagProb();

    cout<<"muon "<<bestMuIndex<<", tag "<<tagDecision<<", mistag "<<osMuonTagMistag.first<<" +- "<<osMuonTagMistag.second<<endl;

    return true;

}

void PDAnalyzer::endJob() {

// additional features
//    tWriter->close();   // second ntuple

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
