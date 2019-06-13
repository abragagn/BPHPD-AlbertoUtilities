#ifndef PDMuonVar_h
#define PDMuonVar_h

#include <vector>
#include <set>
#include <map>
#include <string>

#include "PDAnalysis/Ntu/bin/PDAnalyzerUtil.h"
#include "PDAnalysis/Ntu/interface/PDEnumString.h"

class PDMuonVar: public virtual PDAnalyzerUtil {

 public:

  PDMuonVar();
  ~PDMuonVar();

  void computeMuonVar();

  std::vector<number>* muoInnerChi2;
  std::vector<number>* muoValPixHits;
  std::vector<number>* muoNTrkVHits;
  std::vector<number>* muoOuterChi2;
  std::vector<number>* muoGNchi2;
  std::vector<number>* muoVMuHits;
  std::vector<number>* muoQprod;
  std::vector<number>* muoLWH;

 private:

  PDMuonVar           ( const PDMuonVar& a );
  PDMuonVar& operator=( const PDMuonVar& a );

};

#endif

