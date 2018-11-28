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
#include "PDSecondNtupleWriter.h"   // second ntuple
//#include "DataSetFilter.cc"       // dataset filter
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"

using namespace std;

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2016Lists/BsToJpsiPhi_BMuonFilter_2016_DCAP.list hist.root -v histoMode RECREATE -v use_gen t -n 100000
*/

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // user parameters are set as names associated to a string, 
    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );
    setUserParameter( "minPtMuon", "2." );
    setUserParameter( "maxEtaMuon", "2.4" );
    setUserParameter( "outputFile", "ntu.root" );

    setUserParameter( "mvaMethod", "DNNGlobal2016woIPwIso" ); 

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

    getUserParameter( "verbose", verbose );
    getUserParameter( "minPtMuon", minPtMuon );
    getUserParameter( "maxEtaMuon", maxEtaMuon );

    getUserParameter( "mvaMethod", mvaMethod );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

// to skim the N-tuple "uncomment" the following lines
// // dropBranch( "*tau*" ); // drop some branch if required
// initWSkim( new TFile( "skim.root", "RECREATE" ) );


// additional features
//  DataSetFilter::beginJob();
    tWriter = new PDSecondNtupleWriter;
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    setupMuonMvaReader( mvaMethod );

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility

    autoSavedObject =
    hmass_JPsi              = new TH1D( "hmass_JPsi", "hmass_JPsi", 400, 2.5, 3.5 ); 

    autoSavedObject =
    hmass_Bs                = new TH1D( "hmass_Bs", "hmass_Bs", 100, 5.25, 5.5 ); 


    return;

}


void PDAnalyzer::reset() {
// automatic reset
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
     (!(event_tot %1000) && event_tot<10000 ) || 
     (!(event_tot %10000) && event_tot<100000 ) || 
     (!(event_tot %100000) && event_tot<1000000 ) || 
     (!(event_tot %1000000) && event_tot<10000000 ) )
            cout << " == at event " << event_file << " " << event_tot << endl;
    }

// additional features
    computeMuonVar();
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    // flag to be set "true" or "false" for events accepted or rejected

    nselMu=0;

 // generation information
    vector <int> ListLongLivedB;
        
    for( unsigned int i=0 ; i<genId->size() ; i++ ){
        unsigned int Code = abs(genId->at(i));
        if( Code == 511 || Code == 521 || Code ==531 || Code == 541 || Code == 5122 || Code == 5132 || Code ==5232 ){
            ListLongLivedB.push_back(i);
            continue;
        }
    }

    //search for muons

    for ( int iMuon = 0; iMuon < nMuons; ++iMuon ){

        //SELECTION

        if(muoPt->at( iMuon )<minPtMuon) continue;
        if(abs(muoEta->at( iMuon ))>maxEtaMuon) continue;

        nselMu++;

        int itkmu = muonTrack( iMuon, PDEnumString::muInner );
        if(itkmu < 0) continue;

        // gen info

        int genMuIndex=-1, muoLund=0, muoAncestor=-1;

        genMuIndex = GetClosestGen( muoEta->at(iMuon), muoPhi->at(iMuon), muoPt->at(iMuon) );

        if( genMuIndex >= 0 ) {
            muoLund = genId->at(genMuIndex);
            muoAncestor = GetAncestor( genMuIndex, &ListLongLivedB );
        }

        //TWRITER FILLING

        (tWriter->muoPt)->push_back( muoPt->at( iMuon ) );
        (tWriter->muoEta)->push_back( muoEta->at( iMuon ) );
        (tWriter->muoPhi)->push_back( muoPhi->at(iMuon) );

        (tWriter->trkDxy)->push_back( abs(trkDxy->at(itkmu)) * IPsign(iMuon) );
        (tWriter->trkDz)->push_back( trkDz->at(itkmu) );
        (tWriter->trkExy)->push_back( trkExy->at(itkmu) );
        (tWriter->trkEz)->push_back( trkEz->at(itkmu) );

        (tWriter->muoPFiso)->push_back(GetMuoPFiso(iMuon));

        (tWriter->muoSoftMvaValue)->push_back( computeMva(iMuon) );

        (tWriter->muoLund)->push_back( muoLund );
        (tWriter->muoAncestor)->push_back( muoAncestor ); 

    }

    if( (nselMu) == 0 ) return false;

    tWriter->fill();


// to skim the N-tuple "uncomment" the following line
// if ( flag ) fillSkim();

    return true;

}


void PDAnalyzer::endJob() {
// to skim the N-tuple "uncomment" the following line
//  closeSkim();

// additional features
//  DataSetFilter::endJob();
    tWriter->close();
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


// to plot some histogram immediately after the ntuple loop
// "uncomment" the following lines
/*
void PDAnalyzer::plot() {
    TCanvas* can = new TCanvas( "muoPt", "muoPt", 800, 600 );
    can->cd();
    can->Divide( 1, 2 );
    can->cd( 1 );
    hptmumax->Draw();
    hptmu2nd->Draw();
    return;
}
*/


// ======MY FUNCTIONS===============================================================================