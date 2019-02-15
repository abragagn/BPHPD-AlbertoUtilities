#include "AlbertoUtil.h"

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

AlbertoUtil::AlbertoUtil() {}

AlbertoUtil::~AlbertoUtil() {}


// =====================================================================================
bool AlbertoUtil::IsB( uint genindex ) 
{

    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundBmesons); ++i ) if( genCode == listLundBmesons[i] ) return true;
    for( uint i=0; i<ARRAY_SIZE(listLundBbaryons); ++i ) if( genCode == listLundBbaryons[i] ) return true;
    return false;

}

// =================================================================================================
bool AlbertoUtil::IsBottomium(uint genindex)
{

    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundBottonium); ++i ) if( genCode == listLundBottonium[i] ) return true;
    return false;

}

// =====================================================================================
bool AlbertoUtil::IsC( uint genindex ) 
{

    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundCmesons); ++i ) if( genCode == listLundCmesons[i] ) return true;
    for( uint i=0; i<ARRAY_SIZE(listLundCbaryons); ++i ) if( genCode == listLundCbaryons[i] ) return true;
    return false;

}

// =================================================================================================
bool AlbertoUtil::IsCharmonium(uint genindex)
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundCharmonium); ++i ) if( genCode == listLundCharmonium[i] ) return true;
    return false;

}

// =====================================================================================
bool AlbertoUtil::IsLongLived( uint genindex ) 
{

    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(LongLivedList); ++i ) if( genCode == LongLivedList[i] ) return true;
    return false;

}

// =====================================================================================
int AlbertoUtil::GetClosestGen( float eta, float phi, float pt ) 
{
    double drb = 0.12;
    double dpb = 0.3; 
    int best = -1;
    
    for( uint i=0; i<genId->size(); ++i ){
       if( !IsLongLived(i) ) continue;
       float dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       float dpt = fabs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return best;
}
// =====================================================================================
int AlbertoUtil::GetOverlappedTrack( int trk, vector <int> *List )
{
    double drb = 0.1;
    double dpb = 0.1; 
    int best = -1;
    
    for(int it:*List){

       float dr = deltaR(trkEta->at(trk), trkPhi->at(trk), trkEta->at(it), trkPhi->at(it));
       float dpt = fabs(trkPt->at(it) - trkPt->at(trk))/trkPt->at(it);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = it;
       drb = dr;
    }

    return best;
}


// =====================================================================================
int AlbertoUtil::GetClosestGenLongLivedB( float eta, float phi, float pt, vector <int> *GenList ) {

    double drb = 0.4;
    double dpb = 0.4; 
    int best = -1;
    
    for(int it:*GenList){

       float dr = deltaR(eta, phi, genEta->at(it), genPhi->at(it));
       float dpt = fabs(genPt->at(it) - pt)/genPt->at(it);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = it;
       drb = dr;
    } 

    return best;
}

// =====================================================================================
int AlbertoUtil::GetAncestor( uint iGen, vector <int> *GenList ) 
{
    const vector <int>* aM = &allMothers(iGen);
    while( aM->size()>0 ){ 
       int a = aM->at(0);
       for( uint i=0; i<GenList->size(); ++i ) if( GenList->at(i) == a ) return GenList->at(i);
       aM = &allMothers( aM->at(0) );
    }
    return -1;
}

// =====================================================================================
int AlbertoUtil::WhichMuon(int trk)
{
    for( int iMuon = 0; iMuon<nMuons; ++iMuon){
       if(muonTrack( iMuon, PDEnumString::muInner ) == trk) return iMuon;
    }
    return -1;
}
// =====================================================================================
float AlbertoUtil::GetGenCT( uint genIndex ) 
{

    const vector <int>& aD = allDaughters(genIndex);
    if( aD.size() == 0 ) return -1 ;

    uint mthIndex = aD[0];

    if( genId->at( genIndex ) == - genId->at(genMother->at(genIndex)) ) mthIndex = genMother->at(genIndex ); 
 
    TLorentzVector pGen ;
    pGen.SetPtEtaPhiM( (double) genPt->at(genIndex), (double) genEta->at(genIndex),
         (double) genPhi->at(genIndex), (double) genMass->at(genIndex) );

    float dx = genVx->at(genIndex)-genVx->at(mthIndex);
    float dy = genVy->at(genIndex)-genVy->at(mthIndex);
    float dz = genVz->at(genIndex)-genVz->at(mthIndex);

    return sqrt( dx*dx+dy*dy+dz*dz )/pGen.Beta()/pGen.Gamma();
 
}

// ========================================================================================
int AlbertoUtil::GetBestBstrange()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBsJPsiPhi) ) continue;
       if( svtMass->at(iB)<BsMassRange[0] || svtMass->at(iB)>BsMassRange[1] ) continue;

       if( svtChi2->at(iB)>bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBup()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBuJPsiK) ) continue;
       if( svtMass->at(iB)<BuMassRange[0] || svtMass->at(iB)>BuMassRange[1] ) continue;

       if( svtChi2->at(iB) > bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBdown()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBdJPsiKx) ) continue;
       if( svtMass->at(iB)<BdMassRange[0] || svtMass->at(iB)>BdMassRange[1] ) continue;
       if( svtChi2->at(iB) > bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
bool AlbertoUtil::IsTightJPsi(int iJPsi)
{
    if(fabs(svtMass->at(iJPsi) - MassJPsi) > 0.15 ) return false;

    vector <int> tkJpsi = tracksFromSV(iJPsi);
    TLorentzVector tJPsi(0,0,0,0);

    for( uint i=0; i<tkJpsi.size(); ++i ){
        int j = tkJpsi[i];
        if(trkPt->at(j) < bMuPtCut) return false;
        if(fabs(trkEta->at(j)) > bMuEtaCut) return false;
        TLorentzVector a;
        a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), MassMu );
        tJPsi += a;
    }

    //if(tJPsi.Pt() < 7.0) return false;

    return true;
}

// ========================================================================================
bool AlbertoUtil::IsTightPhi(int iPhi)
{
    if(fabs(svtMass->at(iPhi) - MassPhi) > 0.01 ) return false;
    
    vector <int> tkPhi = tracksFromSV(iPhi);
    for( uint i=0; i<tkPhi.size(); ++i ){
        int j = tkPhi[i];
        if(trkPt->at(j) < bKPtCut) return false;
        if(fabs(trkEta->at(j)) > bKEtaCut) return false;

        int K_Hits = trkHitPattern->at(tkPhi[i]);
        K_Hits = ( int(K_Hits) / 100 ) % 10000;
        K_Hits = K_Hits % 100;
        if(K_Hits<4 ) return false;

    }
    return true;
}


// ========================================================================================
int AlbertoUtil::GetBestBstrangeTight()
{
    int index = -1;
    float best = 0.;
    SetBsMassRange(5.2, 5.65);

    for( int iB=0; iB<nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBsJPsiPhi) ) continue;

        int iJPsi = (subVtxFromSV(iB)).at(0);
        int iPhi = (subVtxFromSV(iB)).at(1);

        vector <int> tkSsB = tracksFromSV(iB);
        vector <int> tkJpsi = tracksFromSV(iJPsi);
        vector <int> tkPhi = tracksFromSV(iPhi);

        //JPSI
        if(!IsTightJPsi(iJPsi)) continue;

        //PHI
        if(!IsTightPhi(iPhi)) continue;

        TLorentzVector tB(0,0,0,0);

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            float m = MassK;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MassMu;
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
        }

        //BS
        float bVprob = ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BsMassRange[0] || svtMass->at(iB)>BsMassRange[1] ) continue;
        if( bVprob < bVprobCut ) continue;
        if(tB.Pt() < bPtCut) continue;

        int PV = GetBestPV(iB, tB);
        if(PV<0) continue;
        if(GetCt2D(tB, iB) < bCtCut) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBupTight()
{
    int index = -1;
    float best = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBuJPsiK) ) continue;

        int iJPsi = (subVtxFromSV(iB)).at(0);
        if(!IsTightJPsi(iJPsi)) continue;

        vector <int> tkJpsi = tracksFromSV(iJPsi);
        vector <int> tkSsB = tracksFromSV(iB);

        TLorentzVector tB(0,0,0,0);
        float KaonPt = 0;

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            float m = MassK;
            if( j == tkJpsi[0] || j == tkJpsi[1] ){ m = MassMu; }else{ KaonPt = trkPt->at(j); }
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
       }

       //Bu
        float bVprob = ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BuMassRange[0] || svtMass->at(iB)>BuMassRange[1] ) continue;
        if( bVprob < bVprobCut ) continue;
        if(tB.Pt() < bPtCut) continue;
        if(KaonPt < 1.6) continue;
        int PV = GetBestPV(iB, tB);
        if(PV<0) continue;
        if(GetCt2D(tB, iB) < bCtCut) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBdownTight()
{
    int index = -1;
    float best = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBdJPsiKx) ) continue;

        int iJPsi = (subVtxFromSV(iB)).at(0);
        if(!IsTightJPsi(iJPsi)) continue;

        vector <int> tkJpsi = tracksFromSV(iJPsi);
        vector <int> tkSsB = tracksFromSV(iB);

        TLorentzVector tB(0,0,0,0);

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            float m = MassKx;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MassMu;
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;

        }

        //Bd
        float bVprob = ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BdMassRange[0] || svtMass->at(iB)>BdMassRange[1] ) continue;
        if( bVprob < bVprobCut ) continue;
        if(tB.Pt() < bPtCut) continue;
        int PV = GetBestPV(iB, tB);
        if(PV<0) continue;
        if(GetCt2D(tB, iB) < bCtCut) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;
    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestJpsi()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( int i=0; i<nSVertices; ++i ){

        if((svtType->at(i)!=PDEnumString::svtJPsi) ) continue;
        if( fabs(svtMass->at(i)-MassJPsi) > MassRangeJPsi) continue;

        if( svtChi2->at(i)>bestChi2 ) continue;
        index = i;
        bestChi2 = svtChi2->at(i);

    }
    return index;
}

// ========================================================================================
float AlbertoUtil::GetInvMass(int i1, int i2, float mass1, float mass2)
{

    float px1=trkPx->at(i1);
    float px2=trkPx->at(i2);
    float py1=trkPy->at(i1);
    float py2=trkPy->at(i2);
    float pz1=trkPz->at(i1);
    float pz2=trkPz->at(i2);

    float E1=sqrt( pow(mass1,2) + pow(px1,2)+pow(py1,2)+pow(pz1,2) );
    float E2=sqrt( pow(mass2,2) + pow(px2,2)+pow(py2,2)+pow(pz2,2) );

    float m=sqrt( pow( (E1+E2),2) - ( pow(px1+px2,2) + pow(py1+py2,2) + pow(pz1+pz2,2) ) );

    return m;
}

// =====================================================================================
int AlbertoUtil::TagMixStatus( uint genindex )
{
    int Code = genId->at( genindex );

    const vector <int>& aD = allDaughters(genindex);
    if( aD.size()>0 && genId->at(aD[0]) == -Code ) return 2;

    const vector <int>& aM = allMothers(genindex); 
    if( aM.size()>0 && genId->at(aM[0]) == -Code ) return 1;

    return 0;
 
}

// ========================================================================================
float AlbertoUtil::GetMuoPFiso (int iMuon)
{

    float PFIso = muoSumCPpt->at(iMuon)/muoPt->at(iMuon);
    float betaCorr = muoSumNHet->at(iMuon) + muoSumPHet->at(iMuon)-0.5*(muoSumPUpt->at(iMuon));
    betaCorr/=muoPt->at(iMuon);
    if(betaCorr>0) PFIso+=betaCorr;

    return PFIso;
}

// ========================================================================================
bool AlbertoUtil::isMvaMuon(int iMuon, float wpB, float wpE)
{

    if((fabs(muoEta->at( iMuon ))<1.2)&&(computeMuonMva(iMuon)>=wpB)) return true;
    if((fabs(muoEta->at( iMuon ))>=1.2)&&(computeMuonMva(iMuon)>=wpE)) return true;

    return false;

}

// =====================================================================================
float AlbertoUtil::GetJetCharge(int iJet, float kappa)
{

    float QJet = 0;
    float ptJet = 0;

    vector <int> list = pfCandFromJet( iJet );

    for(int it:list){

       float pt = pfcPt->at(it);
       float eta = pfcEta->at(it);

       if(pt<0.2) continue;
       if(fabs(eta)>2.5) continue;

       QJet += pfcCharge->at(it) * pow(pt, kappa);
       ptJet += pow(pt, kappa);

    }

    QJet /= ptJet;

    return QJet; 

}

// =====================================================================================
int AlbertoUtil::IPsign(int iMuon)
{
    return IPsign_(iMuon);
}
// =====================================================================================
int AlbertoUtil::IPsign(int iMuon, int iPV)
{

    return IPsign_(iMuon, iPV);
}

// =====================================================================================
float AlbertoUtil::GetJetProbb(int iJet)
{
    float jetdfprob = 0;
    int iTagType = 0;
    for(int iTag=0; iTag<nTags; ++iTag){
       if(tagJet->at(iTag) != iJet) continue;
       ++iTagType;
       if((iTagType == 4) || (iTagType == 5)) jetdfprob += tagProb->at(iTag); 
    }
    return jetdfprob;
}
// =====================================================================================
float AlbertoUtil::CountEventsWithFit(TH1 *hist, TString process){

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    float mean=MassBs;
    if(process=="BsJPsiPhi")   mean=MassBs;
    if(process=="BuJPsiK")     mean=MassBu;
    if(process=="BdJPsiKx")    mean=MassBd;
    if(process=="BdKxMuMu")    mean=MassBd;

    float sigma = 0.015;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [4], true)";
    sgnDef +=       "+[2]*TMath::Gaus(x, [0], [5], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [6], true)";
    TString bkgDef = "[7]+[8]*TMath::Erfc([9]*(x-[10]))";

    if(hist->GetEntries()<=250) bkgDef = "[7]";

    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, hist->GetBinLowEdge(1), hist->GetBinLowEdge(hist->GetNbinsX()));

    //SIGNAL
    float limit = hist->GetEntries()*hist->GetBinWidth(1);

    func->SetParameter(0, mean);
    func->SetParameter(1, limit/3);
    func->SetParameter(2, limit/3);
    func->SetParameter(3, limit/3);
    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);
    func->SetParLimits(0, mean-sigma, mean+sigma);
    func->SetParLimits(1, 0, limit);
    func->SetParLimits(2, 0, limit);
    func->SetParLimits(3, 0, limit);
    func->SetParLimits(4, sigma/2, sigma*2);
    func->SetParLimits(5, sigma/2, sigma*2);
    func->SetParLimits(6, sigma/2, sigma*2);

    //BKG
    float bkgHeight = 0.;
    int nBinsBkgEst = 7;
    for(int i=0; i<nBinsBkgEst; i++)
        bkgHeight += hist->GetBinContent(hist->GetNbinsX()-i);
    bkgHeight /= nBinsBkgEst;

    func->SetParameter(7, bkgHeight);
    if(hist->GetEntries()>250){
        func->SetParameter(8, hist->GetBinContent(1)/2);
        func->SetParameter(9, 10);
        func->SetParameter(10, 5);
        func->SetParLimits(7, 0, bkgHeight*2);
        func->SetParLimits(8, 0, hist->GetBinContent(1));
    }
    

    hist->Fit("func","MRLQ");

    TF1 *fit = hist->GetFunction("func");

    float nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);

    nEvt/=hist->GetBinWidth(0);

    return nEvt;

}
// =====================================================================================
int AlbertoUtil::GetBestPV(int isvt, TLorentzVector t)
{
    int ssPV = -1;
    float bestCos = -1;

    TVector3 vB(t.Px(),t.Py(),t.Pz());
    TVector3 vSVT( svtX->at(isvt), svtY->at(isvt), svtZ->at(isvt) );

    for( int i=0; i<nPVertices; ++i ){

       if(fabs(svtZ->at(isvt) - pvtZ->at( i )) > 0.5 ) continue;

       TVector3 vPV(pvtX->at( i ), pvtY->at( i ), pvtZ->at( i ) );
       TVector3 vPointing;
       vPointing = vSVT - vPV;
       float cos = vPointing.Unit() * vB.Unit();

       if(cos > bestCos ){
         bestCos = cos;
         ssPV = i;
       }

    }

    return ssPV;
}
// =====================================================================================
TLorentzVector AlbertoUtil::GetTLorentzVecFromJpsiX(int iSvt)
{
    int iJPsi = (subVtxFromSV(iSvt)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi) ;
    vector <int> tkSsB = tracksFromSV(iSvt);
    TLorentzVector t(0,0,0,0);

    for( uint i=0; i<tkSsB.size(); ++i ){
        int j = tkSsB.at(i);
        float m = MassK;
        if( (j==tkJpsi.at(0)) || (j==tkJpsi.at(1)) ) m = MassMu;
        TLorentzVector a;
        a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
        t += a;
    }

    return t;
}
// =====================================================================================
int AlbertoUtil::GetTightCandidate(TString process)
{
    if(process=="BsJPsiPhi") return GetBestBstrangeTight();

    if(process=="BuJPsiK") return GetBestBupTight();

    return -1;
}
// =====================================================================================
int AlbertoUtil::GetCandidate(TString process)
{
    if(process=="BsJPsiPhi") return GetBestBstrange();

    if(process=="BuJPsiK") return GetBestBup();

    return -1;
}
// =====================================================================================
float AlbertoUtil::dZ(int itk, int iPV)
{
    return PDAnalyzerUtil::dZ(itk, pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
}
// =====================================================================================
float AlbertoUtil::dXYjet(int itk, int iPV, int iJet)
{
    return fabs(dXY( itk, pvtX->at(iPV), pvtY->at(iPV) ))*dSign( itk, iJet, pvtX->at(iPV), pvtY->at(iPV) );
}
// =====================================================================================
float AlbertoUtil::GetMuonSignedDxy(int iMuon, int iPV)
{
    float dxy = dXY( muonTrack( iMuon, PDEnumString::muInner ), pvtX->at(iPV), pvtY->at(iPV) );
    int sign = IPsign(iMuon, iPV);
    return fabs(dxy)*sign;
}
// =================================================================================================
void AlbertoUtil::printMotherChain(int iGen)
{
    cout<<genId->at(iGen)<<" << ";
    const vector <int>& vM = allMothers(iGen);
    uint nmot = vM.size();
    if(nmot>1) for(uint im=0; im<nmot; ++im) if(genId->at(vM[im])!=21) cout<<genId->at(vM[im])<<" ";
    if(nmot==1) printMotherChain(vM[0]);
    return;
}
// =================================================================================================
void AlbertoUtil::printDaughterTree(int iGen, const string & pre)
{
    cout<<genId->at(iGen)<<endl;

    const vector <int>& vD = allDaughters(iGen);
    uint ndau = vD.size();
    if(ndau == 0) return;

    bool lastLevel = true;
    for(uint id =0; id<ndau; ++id){
        if ( hasDaughter( vD[id] ) ) {
            lastLevel = false;
            break;
        }
    }

    if( lastLevel ){
        cout << pre << "+-> ";
        for( uint id=0; id<ndau; ++id ) {
            int d = vD[id];
            cout<< genId->at( d ) <<" ";
        }
        cout << endl;
        return;
    }

  for( uint id=0; id<ndau; ++id ) {
    int d = vD[id];
    cout << pre << "+-> ";
    string prepre( pre );
    if ( id == ndau - 1 ) prepre += "    ";
    else prepre += "|   ";
    printDaughterTree( d, prepre );
  }
}

// =================================================================================================
bool AlbertoUtil::hasDaughter(int iGen)
{
    const vector <int>& vD = allDaughters(iGen);
    return vD.size()>0 ? true : false;
}
// =================================================================================================
float AlbertoUtil::GetCt2D(TLorentzVector t, int iSV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), 0. );
  TVector3 vPV( bsX, bsY, 0. );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  return MassBs/t.Pt() * (vPointing * vBs.Unit());

}

// =================================================================================================
float AlbertoUtil::GetCt2D(TLorentzVector t, int iSV, int iPV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), 0. );
  TVector3 vPV( pvtX->at(iPV), pvtY->at(iPV), 0. );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  return MassBs/t.Pt() * (vPointing * vBs.Unit());
  
}

// =================================================================================================
float AlbertoUtil::GetCt3D(TLorentzVector t, int iSV, int iPV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );
  TVector3 vPV( pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV) );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  return MassBs/t.P() * vPointing.Dot(vBs)/vBs.Mag();

}
// =================================================================================================
/*float AlbertoUtil::GetCt2DErr(TLorentzVector t, int iSV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), 0. );
  TVector3 vPV( bsX, bsY, 0. );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  TMatrixF covSV(3,3);
  float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
                      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
                      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
  covSV.SetMatrixArray(covSVArray);

  TMatrixF covPV(3,3);
  float covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
                      pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
                      pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
  covPV.SetMatrixArray(covPVArray);

  TMatrixF covTot= covSV+covPV;

  float distArray2D[]={float(vPointing.X()),float(vPointing.Y()),0.};
  TVectorF diff2D(3,distArray2D);

  if (diff2D.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 

  return MassBs/t.Pt() * sqrt(covTot.Similarity(diff2D)) / sqrt(diff2D.Norm2Sqr()); 
   
}
*/
// =================================================================================================
float AlbertoUtil::GetCt2DErr(TLorentzVector t, int iSV, int iPV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), 0. );
  TVector3 vPV( pvtX->at(iPV), pvtY->at(iPV), 0. );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  TMatrixF covSV(3,3);
  float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
                      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
                      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
  covSV.SetMatrixArray(covSVArray);

  TMatrixF covPV(3,3);
  float covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
                      pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
                      pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
  covPV.SetMatrixArray(covPVArray);

  TMatrixF covTot= covSV+covPV;

  float distArray2D[]={float(vPointing.X()),float(vPointing.Y()),0.};
  TVectorF diff2D(3,distArray2D);

  if (diff2D.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 

  return MassBs/t.Pt() * sqrt(covTot.Similarity(diff2D)) / sqrt(diff2D.Norm2Sqr()); 
   
}

// =================================================================================================
float AlbertoUtil::GetCt3DErr(TLorentzVector t, int iSV, int iPV)
{

  TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );
  TVector3 vPV( pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV) );
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = t.Vect();

  TMatrixF covSV(3,3);
  float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
                      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
                      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
  covSV.SetMatrixArray(covSVArray);

  TMatrixF covPV(3,3);
  float covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
                      pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
                      pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
  covPV.SetMatrixArray(covPVArray);

  TMatrixF covTot= covSV+covPV;

  float distArray[]={float(vPointing.X()),float(vPointing.Y()),float(vPointing.Z())};
  TVectorF diff(3,distArray);

  if ( diff.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 
 
  return MassBs/t.P() * sqrt(covTot.Similarity(diff)) / sqrt(diff.Norm2Sqr()); 

}
// =================================================================================================
void AlbertoUtil::SetJpsiMuCut()
{
    SetBctCut(0.007);
    SetBptCut(11.);
    SetBmuptCut(3.5);
    SetBkptCut(1.2);
}
// =================================================================================================
void AlbertoUtil::SetJpsiTrkTrkCut()
{
    SetBctCut(0.02);
    SetBptCut(10.);
    SetBmuptCut(3.6);
    SetBkptCut(1.0);
}