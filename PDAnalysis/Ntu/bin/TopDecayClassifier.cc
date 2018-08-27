#include "TopDecayClassifier.h"
#include <vector>
#include <iostream>


TopDecayClassifier::TopDecayClassifier() {

  n_twl  = 2;
  s_twl  = new const PDGenHandler::TypeSelect*[n_twl ];
  s_twl [0] = &isDirectW;
  s_twl [1] = &isDirectTQuark;

  n_twtl = 3;
  s_twtl = new const PDGenHandler::TypeSelect*[n_twtl];
  s_twtl[0] = &isDirectTau;
  s_twtl[1] = &isDirectW;
  s_twtl[2] = &isDirectTQuark;

  n_tbl  = 3;
  s_tbl  = new const PDGenHandler::TypeSelect*[n_tbl ];
  s_tbl [0] = &isDirectBHadron;
  s_tbl [1] = &isBQuark;
  s_tbl [2] = &isDirectTQuark;

//  n_txbl = 4;
//  s_txbl = new const PDGenHandler::TypeSelect*[n_txbl];
//  s_txbl[0] = &isDirectBHadron;
//  s_txbl[1] = &isBQuark;
//  s_txbl[2] = &isNotTop;
//  s_txbl[3] = &isTQuark;

  n_tbtl = 4;
  s_tbtl = new const PDGenHandler::TypeSelect*[n_tbtl];
  s_tbtl[0] = &isDirectTau;
  s_tbtl[1] = &isDirectBHadron;
  s_tbtl[2] = &isBQuark;
  s_tbtl[3] = &isDirectTQuark;

  n_tbcl = 4;
  s_tbcl = new const PDGenHandler::TypeSelect*[n_tbcl];
  s_tbcl[0] = &isDirectCHNonCharmonium;
  s_tbcl[1] = &isBHadron;
  s_tbcl[2] = &isBQuark;
  s_tbcl[3] = &isDirectTQuark;

  n_tbol = 4;
  s_tbol = new const PDGenHandler::TypeSelect*[n_tbol];
  s_tbol[0] = &isDirectCharmonium;
  s_tbol[1] = &isBHadron;
  s_tbol[2] = &isBQuark;
  s_tbol[3] = &isDirectTQuark;

  n_tbxl = 4;
  s_tbxl = new const PDGenHandler::TypeSelect*[n_tbxl];
  s_tbxl[0] = &isDirectLHadron;
  s_tbxl[1] = &isBHadron;
  s_tbxl[2] = &isBQuark;
  s_tbxl[3] = &isDirectTQuark;

  n_twcl = 3;
  s_twcl = new const PDGenHandler::TypeSelect*[n_twcl];
  s_twcl[0] = &isDirectCHadron;
  s_twcl[1] = &isW;
  s_twcl[2] = &isDirectTQuark;

  n_twxl = 3;
  s_twxl = new const PDGenHandler::TypeSelect*[n_twxl];
  s_twxl[0] = &isDirectLHadron;
  s_twxl[1] = &isW;
  s_twxl[2] = &isDirectTQuark;

  n_bl   = 1;
  s_bl   = new const PDGenHandler::TypeSelect*[n_bl  ];
  s_bl  [0] = &isDirectBHadron;

  n_btl  = 3;
  s_btl  = new const PDGenHandler::TypeSelect*[n_btl ];
  s_btl[0] = &isDirectTau;
  s_btl[1] = &isDirectBHadron;
  s_btl[2] = &isBQuark;

  n_bcl  = 2;
  s_bcl  = new const PDGenHandler::TypeSelect*[n_bcl ];
  s_bcl [0] = &isDirectCHNonCharmonium;
  s_bcl [1] = &isBHadron;

  n_bol  = 2;
  s_bol  = new const PDGenHandler::TypeSelect*[n_bol ];
  s_bol [0] = &isDirectCharmonium;
  s_bol [1] = &isBHadron;

  n_bxl  = 2;
  s_bxl  = new const PDGenHandler::TypeSelect*[n_bxl ];
  s_bxl [0] = &isDirectLHadron;
  s_bxl [1] = &isBHadron;

  n_cl   = 1;
  s_cl   = new const PDGenHandler::TypeSelect*[n_cl  ];
  s_cl  [0] = &isDirectCHadron;

  n_xl   = 1;
  s_xl   = new const PDGenHandler::TypeSelect*[n_xl  ];
  s_xl  [0] = &isDirectLHadron;

  n_tb   = 1;
  s_tb   = new const PDGenHandler::TypeSelect*[n_tb  ];
  s_tb  [0] = &isDirectTQuark;

  n_twq  = 1;
  s_twq  = new const PDGenHandler::TypeSelect*[n_twq ];
  s_twq [0] = &isDirectW;
  s_twq [1] = &isDirectTQuark;

}


TopDecayClassifier::~TopDecayClassifier() {
}


TopDecayMode::type TopDecayClassifier::decayClass( int iGen ) {

  if ( iGen < 0 ) return TopDecayMode::unknown;
  newEvent();
  std::map<int,TopDecayMode::type>::iterator iter =
           decayCache.find( iGen );
  std::map<int,TopDecayMode::type>::iterator iend =
           decayCache.end();
  if ( iter != iend ) return iter->second;
  TopDecayMode::type& dm = decayCache[iGen];

  if ( isLepton( iGen ) ) {
    iGen = firstMother( iGen );

    bool found = false;
    int jGen = iGen;
    while ( true ) {
      std::cout << jGen << ' ' << genId->at( jGen );
//      if ( isDirectLHadron( jGen ) ) std::cout << " DLH";
      if (       isLHadron( jGen ) ) std::cout << "  LH";
//      if ( isDirectCHadron( jGen ) ) std::cout << " DCH";
      if (       isCHadron( jGen ) ) std::cout << "  CH";
 //     if ( isDirectBHadron( jGen ) ) std::cout << " DBH";
      if (       isBHadron( jGen ) ) std::cout << "  BH";
      std::cout << std::endl;
      if (       isBHadron( jGen ) ) found = true;
      else
      if ( found                   ) break;
      if (         isQuark( jGen ) ) break;
      if (         isGluon( jGen ) ) break;
      if (        isString( jGen ) ) break;
      if (       isCluster( jGen ) ) break;
      jGen = genMother->at( jGen );
    }

    bool  twlFlag = PDGenHandler::hasAncestors( iGen, n_twl , s_twl  );
    bool twtlFlag = PDGenHandler::hasAncestors( iGen, n_twtl, s_twtl );
    bool  tblFlag = PDGenHandler::hasAncestors( iGen, n_tbl , s_tbl  );
    bool tbtlFlag = PDGenHandler::hasAncestors( iGen, n_tbtl, s_tbtl );
//    bool tbblFlag = PDGenHandler::hasAncestors( iGen, n_tbbl, s_tbbl );
    bool tbblFlag = false;
    bool tbclFlag = PDGenHandler::hasAncestors( iGen, n_tbcl, s_tbcl );
    bool tbolFlag = PDGenHandler::hasAncestors( iGen, n_tbol, s_tbol );
    bool tbxlFlag = PDGenHandler::hasAncestors( iGen, n_tbxl, s_tbxl );
    bool twclFlag = PDGenHandler::hasAncestors( iGen, n_twcl, s_twcl );
    bool twxlFlag = PDGenHandler::hasAncestors( iGen, n_twxl, s_twxl );
    bool   blFlag = PDGenHandler::hasAncestors( iGen, n_bl  , s_bl   );
    bool  btlFlag = PDGenHandler::hasAncestors( iGen, n_btl , s_btl  );
//    bool  bblFlag = PDGenHandler::hasAncestors( iGen, n_bbl , s_bbl  );
    bool  bblFlag = false;
    bool  bclFlag = PDGenHandler::hasAncestors( iGen, n_bcl , s_bcl  );
    bool  bolFlag = PDGenHandler::hasAncestors( iGen, n_bol , s_bol  );
    bool  bxlFlag = PDGenHandler::hasAncestors( iGen, n_bxl , s_bxl  );
    bool   clFlag = PDGenHandler::hasAncestors( iGen, n_cl  , s_cl   );
    bool   xlFlag = PDGenHandler::hasAncestors( iGen, n_xl  , s_xl   );
    if ( tblFlag || tbtlFlag || tbblFlag || tbclFlag || tbxlFlag || tbolFlag ) {
      int prev = iGen;
      int hGen = -1;
      while ( true ) {
	std::multimap<int,int>& hMap = findAncestor( prev, isBHadron );
        if ( hMap.size() ) prev = hGen = hMap.begin()->second;
        else               break;
      }
//      std::multimap<int,int>& hMap = findAncestor( iGen, isBHadron );
//      int hGen = ( hMap.size() ? hMap.begin()->second : -1 );
      std::multimap<int,int>& bMap = findAncestor( hGen, isBQuark  );
      std::multimap<int,int>::const_iterator iter = bMap.begin();
      std::multimap<int,int>::const_iterator iend = bMap.end();
      int bGen;
      bool wrongChain = true;
      while ( iter != iend ) {
        bGen = iter++->second;
        if ( !findAncestor( bGen, isTQuark ).size() ) continue;
        bGen = lastDaughter( bGen );
        if ( checkBDecay( iGen, bGen ) ) wrongChain = false;
      }
      if ( wrongChain ) tblFlag = tbtlFlag = //tbblFlag =
                       tbclFlag = tbolFlag = tbxlFlag = false;

/*
//      int bGen = ( bMap.size() ? bMap.begin()->second : -2 );
      int bGen = -2;
      int qGen;
      std::multimap<int,int>::const_iterator iter = bMap.begin();
      std::multimap<int,int>::const_iterator iend = bMap.end();
      const TypeSelect* tSel = &isTQuark;
      while ( iter != iend ) {
        if ( hasAncestors( qGen = iter++->second, 1, &tSel ) ) bGen = qGen;
      }
//      std::cout << iGen << " " << hGen << " " << bGen << std::endl;
      int iDau;
      int nDau;
      int sGen = -1;
      int fGen = -1;
      int jGen;
      bGen = lastDaughter( bGen );
//      std::cout << iGen << " " << hGen << " " << bGen << std::endl;
      const std::vector<int>& qDau = allDaughters( bGen );
      nDau = qDau.size();
      for ( iDau = 0; iDau < nDau; ++iDau ) {
        jGen = qDau[iDau];
//	std::cout << jGen << "(" << genId->at( jGen ) << ")" << std::endl;
        if ( isCluster( jGen ) ) sGen = jGen;
        if ( isString ( jGen ) ) sGen = jGen;
        if ( isBHadron( jGen ) ) fGen = jGen;
//        std::cout << "1-> " << sGen << " " << fGen << std::endl;
//        if ( isString( iDau ) ) {
//          sGen = iDau;
//          break;
//        }
      }
      sGen = lastDaughter( sGen );
      const std::vector<int>& fDau = allDaughters( sGen );
      nDau = ( fGen < 0 ? fDau.size() : 0 );
      for ( iDau = 0; iDau < nDau; ++iDau ) {
        jGen = fDau[iDau];
//	std::cout << jGen << "(" << genId->at( jGen ) << ")" << std::endl;
        if ( isBHadron( jGen ) ) {
          fGen = jGen;
          break;
        }
//        std::cout << "2-> " << sGen << " " << fGen << std::endl;
      }
*/
/*
      int code = ( fGen >= 0 ? genId->at( fGen ) : 0 );
      int sign = ( code > 0 ? 1 : -1 );
      if ( code < 0 ) code = -code;
      code %= 10000;
      if ( code < 1000 ) sign = -sign;
      int prod = ( bGen >= 0 ? genId->at( bGen ) * sign : -1 );
*/
/*
      if ( !checkBHadron( fGen, bGen ) ) std::cout << iGen << " " << tblFlag << tbclFlag << tbxlFlag << " WRONG DECAY CHAIN " << runNumber << " " << eventNumber << " " << hGen << " " << bGen << " " << sGen << " " << fGen << std::endl;
      if ( !checkBHadron( fGen, bGen ) ) tblFlag = tbclFlag = tbxlFlag = false;
//      std::multimap<int,int>& sMap = findDescendant( bGen, isString  );
//      int sGen = ( sMap.size() ? sMap.begin()->second : -1 );
//      std::multimap<int,int>& fMap = findDescendant( sGen, isBHadron );
//      int fGen = ( fMap.size() ? fMap.begin()->second : -1 );
//      if ( hGen * bGen > 0 ) std::cout << hGen << " " << bGen 
*/
    }
    bool cBar = false;
//    if ( tbblFlag || tbclFlag || bblFlag || bclFlag ) {
    if ( tbclFlag || bclFlag ) {
      int cGen = genMother->at( iGen );
      int cHad = genId->at( cGen );
      int bGen = genMother->at( cGen );
      int bHad = genId->at( bGen );
      cBar = ( ( cHad * bHad ) > 0 );
    }
    if ( cBar ) {
      if ( tbclFlag ) {
        tbblFlag = true;
        tbclFlag = false;
      }
      if ( bclFlag ) {
        bblFlag = true;
        bclFlag = false;
      }
    }
    if (  twlFlag                           ) return dm = TopDecayMode::twl;
    if ( twtlFlag                           ) return dm = TopDecayMode::twtl;
    if (  tblFlag                           ) return dm = TopDecayMode::tbl;
    if ( tbtlFlag                           ) return dm = TopDecayMode::tbtl;
    if ( tbblFlag                           ) return dm = TopDecayMode::tbbl;
    if ( tbclFlag                           ) return dm = TopDecayMode::tbcl;
    if ( tbolFlag                           ) return dm = TopDecayMode::tbol;
    if ( tbxlFlag                           ) return dm = TopDecayMode::tbxl;
    if ( twclFlag                           ) return dm = TopDecayMode::twcl;
    if ( twxlFlag                           ) return dm = TopDecayMode::twxl;
    if (   blFlag && ! tblFlag              ) return dm = TopDecayMode::bl;
    if (  btlFlag && !tbtlFlag              ) return dm = TopDecayMode::btl;
    if (  bblFlag && !tbblFlag              ) return dm = TopDecayMode::bbl;
    if (  bclFlag && !tbclFlag              ) return dm = TopDecayMode::bcl;
    if (  bolFlag && !tbolFlag              ) return dm = TopDecayMode::bol;
    if (  bxlFlag && !tbxlFlag              ) return dm = TopDecayMode::bxl;
    if (   clFlag && ! bclFlag && !tbclFlag ) return dm = TopDecayMode::cl;
    if (   xlFlag && ! bxlFlag && !tbxlFlag ) return dm = TopDecayMode::xl;
    return dm = TopDecayMode::other;
  }

  if ( isQuark( iGen ) ) {
    bool   tbFlag = PDGenHandler::hasAncestors( iGen, n_tb  , s_tb   );
    bool  twqFlag = PDGenHandler::hasAncestors( iGen, n_twq , s_twq  );
    if ( isBQuark( iGen ) ) return dm = (  tbFlag ? TopDecayMode::tb  :
                                                    TopDecayMode::xb );
    else                    return dm = ( twqFlag ? TopDecayMode::twq :
                                                    TopDecayMode::hadronic );
  }

  return TopDecayMode::unknown;

}


bool TopDecayClassifier::oscillatedLeptonic( int iGen ) {
  TopDecayMode::type mode = decayClass( iGen );
  switch ( mode ) {
  case TopDecayMode::tbl :
  case TopDecayMode::tbtl:
  case TopDecayMode::tbbl:
  case TopDecayMode::tbcl:
  case TopDecayMode::tbol:
  case TopDecayMode::tbxl:
  case TopDecayMode::bl  :
  case TopDecayMode::btl :
  case TopDecayMode::bbl :
  case TopDecayMode::bcl :
  case TopDecayMode::bol :
  case TopDecayMode::bxl :
    {
      std::multimap<int,int>& bmap = findAncestor( iGen, isBHadron );
      std::multimap<int,int>::const_iterator iter = bmap.begin();
      std::multimap<int,int>::const_iterator iend = bmap.end();
      if ( iter == iend ) return false;
      std::cout << "--- " << iGen << ' ' << iter->second << ' ' << hasOscillated( iter->second ) << std::endl;
      return hasOscillated( iter->second ) > 0;
    }
    break;
  default:
    return false;
  }
  return false;
}


// check quark-hadron charge
bool TopDecayClassifier::checkBHadron( int hGen, int bGen ) {
  if ( !isBHadron( hGen ) ) return false;
  if ( !isBQuark ( bGen ) ) return false;
  int code = chkId( hGen );
  int sign = ( code > 0 ? 1 : -1 );
  if ( code < 0 ) code = -code;
  code %= 10000;
  int ctmp = code / 10;
  int nbqk = 0;
  while ( ctmp ) {
    if ( ( ctmp  % 10 ) == 5 ) ++nbqk;
    ctmp /= 10;
  }
  if ( nbqk >= 2 ) return true;
  if ( code < 1000 ) sign = -sign;
  int prod = ( bGen >= 0 ? genId->at( bGen ) * sign : -1 );
  return ( prod > 0 );
}


// check quark-hadron charge
bool TopDecayClassifier::checkHadron( int hGen, int qGen ) {
  if ( !isHadron( hGen ) ) return false;
  if ( !isQuark ( qGen ) ) return false;
  int hCode = chkId( hGen );
  int qCode = chkId( qGen );
  int hType = ( hCode > 0 ? hCode : -hCode );
  int qType = ( qCode > 0 ? qCode : -qCode );
  hType %= 10000;
//  int qHadr = hType / 100;
//  if ( qHadr > 10 ) qHadr /= 10;
  int hTtmp = hType / 10;
  int qCtmp = ( hCode > 0 ? qCode : -qCode );
  if ( hTtmp == 21 ) {
    if ( qCtmp ==  2 ) return true;
    if ( qCtmp == -1 ) return true;
    return false;
  }
  if ( ( hTtmp == 11  ) ||
       ( hTtmp == 22  ) ||
       ( hType == 331 ) ) return ( qType <= 2 );
  if ( ( hType == 130 ) ||
       ( hType == 310 ) ) return ( qType == 3 );
//  if ( ( hTtmp == 21 ) && ( hCode > 0 ) && ( qCode ==  2 ) ) return true;
//  if ( ( hTtmp == 21 ) && ( hCode > 0 ) && ( qCode == -1 ) ) return true;
//  if ( ( hTtmp == 21 ) && ( hCode < 0 ) && ( qCode == -2 ) ) return true;
//  if ( ( hTtmp == 21 ) && ( hCode < 0 ) && ( qCode ==  1 ) ) return true;
//  if ( ( hTtmp == 11 ) && ( qCode < 3 ) ) return true;
//  if ( ( hTtmp == 22 ) && ( qCode < 3 ) ) return true;
  int qHadr = 0;
  int qHtmp;
  bool qqBar = false;
  while ( hTtmp ) {
    qHtmp = hTtmp % 10;
    if ( qHtmp == qHadr )   qqBar = true;
    if ( qHtmp >  qHadr ) { qqBar = false; qHadr = qHtmp; }
    hTtmp /= 10;
  }
  if ( qHadr != qType ) return false;
  if ( qqBar ) return true;
  int sign = ( ( hType < 1000 ) && ( qType % 2 ) ? -1 : 1 );
  int prod = hCode * qCode * sign;
  //std::cout << "checkHadron " << hCode << " " << qCode << " " << sign << std::endl;
  return ( prod > 0 );
}


// check quark-lepton decay chain
bool TopDecayClassifier::checkBDecay( int lGen, int bGen ) {

  if ( !isBQuark( bGen ) ) return false;
  if ( !isLepton( lGen ) ) return false;

  bGen = lastDaughter( bGen );
  std::vector<int> hList;
  std::vector<int> sList;
  hList.reserve( 10 );
  sList.reserve( 10 );
  int iList;
  int nList;

  const std::vector<int>& qDau = allDaughters( bGen );
  int iDau;
  int nDau;
  int iGen;
  nDau = qDau.size();
  for ( iDau = 0; iDau < nDau; ++iDau ) {
    iGen = qDau[iDau];
    if ( isCluster( iGen ) ) sList.push_back( iGen );
    if ( isString ( iGen ) ) sList.push_back( iGen );
    if ( isBHadron( iGen ) ) hList.push_back( iGen );
  }
  nList = sList.size();
  for ( iList = 0; iList < nList; ++iList ) {
    iGen = sList[iList];
    const std::vector<int>& fDau = allDaughters( iGen );
    nDau = fDau.size();
    for ( iDau = 0; iDau < nDau; ++iDau ) {
      iGen = fDau[iDau];
      if ( isBHadron( iGen ) ) hList.push_back( iGen );
    }
  }

  nList = hList.size();
  for ( iList = 0; iList < nList; ++iList ) {
    iGen = hList[iList];
//    std::cout << "checkBHadron " << iGen << std::endl;
    if ( !checkBHadron( iGen, bGen ) ) continue;
//    std::cout << "check index " << iGen << std::endl;
    TypeSelect hSel( iGen, 'p' );
    if ( findAncestor( lGen, hSel ).size() ) return true;
  }

  return false;

}


// check quark-lepton decay chain
bool TopDecayClassifier::checkDecay( int lGen, int qGen ) {

  if ( !isQuark ( qGen ) ) return false;
  if ( !isLepton( lGen ) ) return false;

  qGen = lastDaughter( qGen );
  std::vector<int> hList;
  std::vector<int> sList;
  hList.reserve( 10 );
  sList.reserve( 10 );
  int iList;
  int nList;

  const std::vector<int>& qDau = allDaughters( qGen );
  int iDau;
  int nDau;
  int iGen;
  nDau = qDau.size();
  for ( iDau = 0; iDau < nDau; ++iDau ) {
    iGen = qDau[iDau];
    if ( isCluster( iGen ) ) sList.push_back( iGen );
    if ( isString ( iGen ) ) sList.push_back( iGen );
    if ( isHadron ( iGen ) ) hList.push_back( iGen );
  }
  nList = sList.size();
  for ( iList = 0; iList < nList; ++iList ) {
    iGen = sList[iList];
    const std::vector<int>& fDau = allDaughters( iGen );
    nDau = fDau.size();
    for ( iDau = 0; iDau < nDau; ++iDau ) {
      iGen = fDau[iDau];
      if ( isHadron( iGen ) ) hList.push_back( iGen );
    }
  }

  nList = hList.size();
  for ( iList = 0; iList < nList; ++iList ) {
    iGen = hList[iList];
    if ( !checkHadron( iGen, qGen ) ) continue;
    TypeSelect hSel( iGen, 'p' );
    if ( findAncestor( lGen, hSel ).size() ) return true;
  }

  return false;

}


const std::string& TopDecayClassifier::decayName(
                                       TopDecayMode::type mode ) {
  static std::vector<std::string> dNames( 1000, "" );
  static bool init = true;
  if ( init ) {
    init = false;
    dNames[TopDecayMode::twl     ] = "twl";
    dNames[TopDecayMode::twtl    ] = "twtl";
    dNames[TopDecayMode::tbl     ] = "tbl";
    dNames[TopDecayMode::tbtl    ] = "tbtl";
    dNames[TopDecayMode::tbbl    ] = "tbbl";
    dNames[TopDecayMode::tbcl    ] = "tbcl";
    dNames[TopDecayMode::tbol    ] = "tbol";
    dNames[TopDecayMode::tbxl    ] = "tbxl";
    dNames[TopDecayMode::twcl    ] = "twcl";
    dNames[TopDecayMode::twxl    ] = "twxl";
    dNames[TopDecayMode::bl      ] = "bl";
    dNames[TopDecayMode::btl     ] = "btl";
    dNames[TopDecayMode::bbl     ] = "bbl";
    dNames[TopDecayMode::bcl     ] = "bcl";
    dNames[TopDecayMode::bol     ] = "bol";
    dNames[TopDecayMode::bxl     ] = "bxl";
    dNames[TopDecayMode::cl      ] = "cl";
    dNames[TopDecayMode::xl      ] = "xl";
    dNames[TopDecayMode::other   ] = "other";
    dNames[TopDecayMode::tb      ] = "tb";
    dNames[TopDecayMode::twq     ] = "twq";
    dNames[TopDecayMode::xb      ] = "xb";
    dNames[TopDecayMode::hadronic] = "hadronic";
    dNames[TopDecayMode::unknown ] = "unknown";
  }
  return dNames[mode];
}


void TopDecayClassifier::newEvent() {
  if ( ( lastRun   != runNumber   ) || 
       ( lastLumi  != lumiSection ) || 
       ( lastEvent != eventNumber ) ) {
    decayCache.clear();
    lastRun   = runNumber;
    lastLumi  = lumiSection;
    lastEvent = eventNumber;
  }
  return;
}

