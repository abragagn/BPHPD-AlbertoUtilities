#ifndef TopDecayClassifier_h
#define TopDecayClassifier_h

#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "TopDecayMode.h"
#include "string"
#include "map"

class TopDecayClassifier: public virtual PDGenHandler {

// public:

 protected:

  TopDecayClassifier();
  virtual ~TopDecayClassifier();

  TopDecayMode::type decayClass( int iGen );
  bool oscillatedLeptonic( int iGen );
  bool checkBHadron( int hGen, int bGen ); // check quark-hadron charge
  bool checkHadron ( int hGen, int qGen ); // check quark-hadron charge
  bool checkBDecay ( int lGen, int bGen ); // check quark-lepton decay chain
  bool checkDecay  ( int lGen, int qGen ); // check quark-lepton decay chain
  static const std::string& decayName( TopDecayMode::type mode );

 private:

  void newEvent();
  unsigned int lastRun;
  unsigned int lastLumi;
  unsigned int lastEvent;
  std::map<int,TopDecayMode::type> decayCache;

  int n_twl ;
  int n_twtl;
  int n_tbl ;
//  int n_txbl;
  int n_tbtl;
//  int n_tbbl;
  int n_tbcl;
  int n_tbol;
  int n_tbxl;
  int n_twcl;
  int n_twxl;
  int n_bl  ;
  int n_btl ;
//  int n_bbl ;
  int n_bcl ;
  int n_bol ;
  int n_bxl ;
  int n_cl  ;
  int n_xl  ;
  int n_tb  ;
  int n_twq ;
//  int n_xb  ;

  const PDGenHandler::TypeSelect** s_twl ;
  const PDGenHandler::TypeSelect** s_twtl;
  const PDGenHandler::TypeSelect** s_tbl ;
//  const PDGenHandler::TypeSelect** s_txbl;
  const PDGenHandler::TypeSelect** s_tbtl;
//  const PDGenHandler::TypeSelect** s_tbbl;
  const PDGenHandler::TypeSelect** s_tbcl;
  const PDGenHandler::TypeSelect** s_tbol;
  const PDGenHandler::TypeSelect** s_tbxl;
  const PDGenHandler::TypeSelect** s_twcl;
  const PDGenHandler::TypeSelect** s_twxl;
  const PDGenHandler::TypeSelect** s_bl  ;
  const PDGenHandler::TypeSelect** s_btl ;
//  const PDGenHandler::TypeSelect** s_bbl ;
  const PDGenHandler::TypeSelect** s_bcl ;
  const PDGenHandler::TypeSelect** s_bol ;
  const PDGenHandler::TypeSelect** s_bxl ;
  const PDGenHandler::TypeSelect** s_cl  ;
  const PDGenHandler::TypeSelect** s_xl  ;
  const PDGenHandler::TypeSelect** s_tb  ;
  const PDGenHandler::TypeSelect** s_twq ;

  TopDecayClassifier           ( const TopDecayClassifier& x );
  TopDecayClassifier& operator=( const TopDecayClassifier& x );

};

#endif

