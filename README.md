# OS Muon Tagger for BsJpsiPhi CMS Run-2 Analysis

This repository contains the class needed to use the OS Muon code. 

As today it has been developed and calibrated only for the JpsiMuon trigger.

The needed weights are stored locally at 

```
/lustre/cmswork/abragagn/mvaWeights/
```

and remotely at

https://github.com/abragagn/PDMvaMethods-WeightFiles

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

Edit PDAnalyzed.cc ro include in the header
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


## Usage example



## Meta

Alberto Bragagnolo – alberto.bragagnolo@cern.ch
