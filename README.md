# OS Muon Tagger for BsJpsiPhi CMS Run-2 Analysis

This repository contains the class needed to use the OS Muon code. 

As today it has been developed and calibrated only for the JpsiMuon trigger.

The needed weight/calibration files are stored locally at 

```
/lustre/cmswork/abragagn/mvaWeights/
```

and remotely at

https://github.com/abragagn/PDMvaMethods-WeightFiles

The calibration files contained in 
```
mvaWeights/OsMuonTag
```
store:
* the calibration function for the given process
* the TGraphAsymmetricErrors where the function was actually fitted
* the TFitResultPtr of the fit

## Installation

Extract the code in your .../workingArea/src/

```
cd .../workingArea/src/
tar -xvzf albertoCode.tar.gz
```
Edit PDAnalyzer.h to include in the header:

```
#include "PDMuonVar.h"
#include "PDSoftMuonMvaEstimator.h"
#include "AlbertoUtil.h"
#include "OSMuonMvaTag.h"
```

and in the list of base classes:

```
,                 public virtual PDMuonVar
,                 public virtual PDSoftMuonMvaEstimator
,                 public virtual AlbertoUtil
,                 public virtual OSMuonMvaTag
```

Edit PDAnalyzer.cc ro include in the header
```
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"
```

May be necessary to include additional ROOT classes such as TLorenzVector, TFile, etc, follow your compiler complains
```
#include "TLorentzVector.h"
#include "TFile.h"
```

Compile:

```
cd .../workingArea/src/
scram b clean
scram b
```

## Usage

> A set of PDAnalyzer.* examples can be found in PDAnalysis/Ntu/bin/examples/


**The following lines of code need to be included in the code to correcly initialize the various methods:**

* In PDAnalyzer::beginJob() 
```
inizializeMuonMvaReader(); // initialize TMVA methods for muon ID
inizializeOSMuonMvaReader(); // initialize TMVA methods for muon tagger
bool osInit = inizializeOSMuonCalibration("BuJPsiKData2018", "BuJPsiKMC2018", "BsJPsiPhiMC2018"); 
// initialize calibration methods for muon tagger, the parameters are respectively:
// the calibration process, and the two MC processes you want to construct the Bu->Bs transformation
```

* In PDAnalyzer::analyze() at the top
```
computeMuonVar(); // compute variable needed for the muon ID
inizializeOsMuonTagVars(); // initialize variables for muon tagger 
convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz); // needed for the methods
convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz); // needed for the methods
convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz); // needed for the methods
convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz); // needed for the methods
```

* In PDAnalyzer::analyze() after the SVT and PVT selection
```
setVtxOsMuonTag(ssbSVT, ssbPVT); // set vertices index dor muon tagger
```

**All tagging variables are computed with the following function:**

```
makeOsMuonTagging();
```
**And can be retrieved with the following set of function:**
```
int     getOsMuon(); // os muon index
int     getOsMuonTag(); // Tag decision // 1*trkCharge->at(osMuonTrackIndex_), 0 -> no muon.
float   getOsMuonTagMvaValue();
float   getOsMuonTagMistagProbRaw(); // 1- dnn score
float   getOsMuonTagMistagProbCalProcess(); // mistag calibrated to the calibration choosen as a first argument of in inizializeOSMuonCalibration()
float   getOsMuonTagMistagProbCalProcessBuBs(); // mistag calibrated with the Bu->Bs correction using the calibration files specified as second and third arguments of in inizializeOSMuonCalibration()


```

See OSMuonMvaTag.cc/.h for more informations. 

### Addedum: usage outside PD network

The code is flexible but all the default paths to the various TMVA weights are related to the PD network.

To use the code outside PD please look at the function declarations in the source code. Usually the path to the weiths is one of the optional arguments (if not please tell me).

**Unavoidable (to my knowledge) nuisance**

Keras pyTMVA .xlm weights files have hardcoded the path to the Keras .h5 weight files (both of them are needed) in the following line:
```
 <Option name="FilenameTrainedModel" modified="No">/path/</Option>
```

e.g.

```
 <Option name="FilenameTrainedModel" modified="No">/lustre/cmswork/abragagn/BPH/BTag/osMuonV13/src/PDAnalysis/Ntu/bin/mvaTraining/dataset/weights/TrainedModel_DNNOsMuonHLTJpsiMu.h5</Option>
```

You need to edit this line to reflect your new path in the following files:
```
TMVAClassification_DNNOsMuonHLTJpsiMu.weights
TMVAClassification_DNNMuonID.weights
```


## Contacts

Alberto Bragagnolo – alberto.bragagnolo@cern.ch
