#include "PDMuonVar.h"

#include <iostream>
#include <sstream>
#include <stdio.h>

PDMuonVar::PDMuonVar() {
  muoInnerChi2  = new std::vector<number>;
  muoValPixHits = new std::vector<number>;
  muoNTrkVHits  = new std::vector<number>;
  muoOuterChi2  = new std::vector<number>;
  muoGNchi2     = new std::vector<number>;
  muoVMuHits    = new std::vector<number>;
  muoQprod      = new std::vector<number>;
  muoLWH        = new std::vector<number>;
}


PDMuonVar::~PDMuonVar() {
}


void PDMuonVar::computeMuonVar() {
  muoInnerChi2 ->clear();
  muoValPixHits->clear();
  muoNTrkVHits ->clear();
  muoOuterChi2 ->clear();
  muoGNchi2    ->clear();
  muoVMuHits   ->clear();
  muoQprod     ->clear();
  muoLWH       ->clear();
  muoInnerChi2 ->resize( nMuons, -1 );
  muoValPixHits->resize( nMuons, -1 );
  muoNTrkVHits ->resize( nMuons, -1 );
  muoOuterChi2 ->resize( nMuons, -1 );
  muoGNchi2    ->resize( nMuons, -1 );
  muoVMuHits   ->resize( nMuons, -1 );
  muoQprod     ->resize( nMuons, -1 );
  muoLWH       ->resize( nMuons, -1 );
  const std::vector<int>& iTracks = muonTracks( PDEnumString::muInner      );
  const std::vector<int>& oTracks = muonTracks( PDEnumString::muStandalone );
  const std::vector<int>& gTracks = muonTracks( PDEnumString::muGlobal     );
  int iMu;
  for ( iMu = 0; iMu < nMuons; ++iMu ) {
    int iTk = iTracks[iMu];
    int oTk = oTracks[iMu];
    int gTk = gTracks[iMu];
    if ( iTk >= 0 ) {
      muoInnerChi2 ->at( iMu ) = trkNormChi2->at( iTk );
      int hp = trkHitPattern->at( iTk ) / 100;
      int lp = trkLayPattern->at( iTk );
      muoNTrkVHits ->at( iMu ) = hp % 100;
      hp /= 100;
      muoValPixHits->at( iMu ) = hp % 100;
      muoLWH       ->at( iMu ) = lp % 100;
    }
    if ( oTk >= 0 ) {
      muoOuterChi2 ->at( iMu ) = trkNormChi2->at( oTk );
      if ( iTk >= 0 )
      muoQprod     ->at( iMu ) = trkCharge->at( iTk ) * trkCharge->at( oTk );
    }
    if ( gTk >= 0 ) {
      muoGNchi2    ->at( iMu ) = trkNormChi2->at( gTk );
      muoVMuHits   ->at( iMu ) = trkHitPattern->at( gTk ) / 1000000;
    }
  }
}

