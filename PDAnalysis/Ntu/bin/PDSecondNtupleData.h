#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset()  { autoReset(); }

    PDSecondNtupleData() {

}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "ssbMass", &ssbMass, "ssbMass/F", &b_ssbMass );
    setBranch( "osMuonTag", &osMuonTag, "osMuonTag/I", &b_osMuonTag );
    setBranch( "osMuonTagMistag", &osMuonTagMistag, "osMuonTagMistag/F", &b_osMuonTagMistag );

}

int osMuonTag;
float ssbMass, osMuonTagMistag;

TBranch *b_osMuonTag, *b_ssbMass, *b_osMuonTagMistag;
private:

    PDSecondNtupleData         ( const PDSecondNtupleData& a );
    PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

