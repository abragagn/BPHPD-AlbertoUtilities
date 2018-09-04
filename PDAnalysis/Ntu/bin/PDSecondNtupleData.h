#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset()  { autoReset(); }

    PDSecondNtupleData() {

        muoPt           = new vector <float>; 
        muoEta          = new vector <float>;
        muoPhi          = new vector <float>;
        trkDxy          = new vector <float>;
        trkExy          = new vector <float>;
        trkDz           = new vector <float>;
        trkEz           = new vector <float>;

        muoPFiso        = new vector <float>;

        muoLund         = new vector <int>;
        muoAncestor     = new vector <int>;

        muoSoftMvaValue = new vector <float>;


}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "muoPt", &muoPt , 8192, 99, &b_muoPt );
    setBranch( "muoEta", &muoEta , 8192, 99, &b_muoEta );
    setBranch( "muoPhi", &muoPhi , 8192, 99, &b_muoPhi );
    setBranch( "trkDxy", &trkDxy , 8192, 99, &b_trkDxy );
    setBranch( "trkExy", &trkExy , 8192, 99, &b_trkExy );
    setBranch( "trkDz", &trkDz , 8192, 99, &b_trkDz );
    setBranch( "trkEz", &trkEz , 8192, 99, &b_trkEz );

    setBranch( "muoPFiso", &muoPFiso , 8192, 99, &b_muoPFiso );

    setBranch( "muoLund", &muoLund , 8192, 99, &b_muoLund );
    setBranch( "muoAncestor", &muoAncestor , 8192, 99, &b_muoAncestor );

    setBranch( "muoSoftMvaValue", &muoSoftMvaValue , 8192, 99, &b_muoSoftMvaValue );

}

vector <int> *muoLund, *muoAncestor;
vector <float> *muoSoftMvaValue;
vector <float> *muoPt, *muoEta, *muoPhi, *trkDxy, *trkExy, *trkDz, *trkEz, *muoPFiso;

TBranch *b_muoPt, *b_muoEta, *b_muoPhi, *b_trkDxy, *b_trkExy, *b_trkDz, *b_trkEz, *b_muoPFiso;
TBranch *b_muoLund, *b_muoAncestor, *b_muoSoftMvaValue;

private:

    PDSecondNtupleData         ( const PDSecondNtupleData& a );
    PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

