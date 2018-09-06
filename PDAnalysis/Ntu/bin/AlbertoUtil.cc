#include "AlbertoUtil.h"

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

AlbertoUtil::AlbertoUtil() {}

AlbertoUtil::~AlbertoUtil() {}


// =====================================================================================
bool AlbertoUtil::IsB( unsigned int genindex ) 
{

    unsigned int genCode = abs( genId->at(genindex) );
    for( unsigned int i=0; i<ARRAY_SIZE(listLundBmesons); ++i ) if( genCode == listLundBmesons[i] ) return true;
    for( unsigned int i=0; i<ARRAY_SIZE(listLundBbaryons); ++i ) if( genCode == listLundBbaryons[i] ) return true;
    return false;

}

// =================================================================================================
bool AlbertoUtil::IsBottomium(unsigned int genindex)
{

    unsigned int genCode = abs( genId->at(genindex) );
    for( unsigned int i=0; i<ARRAY_SIZE(listLundBottonium); ++i ) if( genCode == listLundBottonium[i] ) return true;
    return false;

}

// =====================================================================================
bool AlbertoUtil::IsC( unsigned int genindex ) 
{

    unsigned int genCode = abs( genId->at(genindex) );
    for( unsigned int i=0; i<ARRAY_SIZE(listLundCmesons); ++i ) if( genCode == listLundCmesons[i] ) return true;
    for( unsigned int i=0; i<ARRAY_SIZE(listLundCbaryons); ++i ) if( genCode == listLundCbaryons[i] ) return true;
    return false;

}

// =================================================================================================
bool AlbertoUtil::IsCharmonium(unsigned int genindex)
{
    unsigned int genCode = abs( genId->at(genindex) );
    for( unsigned int i=0; i<ARRAY_SIZE(listLundCharmonium); ++i ) if( genCode == listLundCharmonium[i] ) return true;
    return false;

}

// =====================================================================================
bool AlbertoUtil::IsLongLived( unsigned int genindex ) 
{

    unsigned int genCode = abs( genId->at(genindex) );
    for( unsigned int i=0; i<ARRAY_SIZE(LongLivedList); ++i ) if( genCode == LongLivedList[i] ) return true;
    return false;

}

// =====================================================================================
int AlbertoUtil::GetClosestGen( float eta, float phi, float pt ) 
{
    double drb = 0.12;
    double dpb = 0.3; 
    int best = -1;
    
    for( unsigned int i=0; i< genId->size(); ++i ){
       if( !IsLongLived(i) ) continue;
       float dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       float dpt = abs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return best;
}

// =====================================================================================
int AlbertoUtil::GetClosestGenLongLivedB( float eta, float phi, float pt, vector <int> *GenList ) {

    double drb = 0.4;
    double dpb = 0.4; 
    int best = -1;
    
    for(std::vector<int>::iterator it = GenList->begin(); it != GenList->end(); ++it){

       float dr = deltaR(eta, phi, genEta->at(*it), genPhi->at(*it));
       float dpt = abs(genPt->at(*it) - pt)/genPt->at(*it);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = *it;
       drb = dr;
    } 

    return best;
}

// =====================================================================================
int AlbertoUtil::GetAncestor( unsigned int iGen, vector <int> *GenList ) 
{
    const vector <int>* aM = &allMothers(iGen);
    while( aM->size()>0 ){ 
       int a = aM->at(0);
       for( unsigned int i=0; i<GenList->size(); ++i ) if( GenList->at(i) == a ) return i;
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
float AlbertoUtil::GetCT( unsigned int genIndex ) 
{

    const vector <int>& aD = allDaughters(genIndex);
    if( aD.size() == 0 ) return -1 ;

    unsigned int mthIndex = aD[0];

    if( genId->at( genIndex ) == - genId->at(genMother->at(genIndex)) ) mthIndex = genMother->at(genIndex ); 
 
    TLorentzVector pGen ;
    pGen.SetPtEtaPhiM( (double) genPt->at(genIndex), (double) genEta->at(genIndex),
         (double) genPhi->at(genIndex), (double) genMass->at(genIndex) );

    float dx = genVx->at(genIndex)-genVx->at(mthIndex);
    float dy = genVy->at(genIndex)-genVy->at(mthIndex);
    float dz = genVz->at(genIndex)-genVz->at(mthIndex);

    float ct = sqrt( dx*dx+dy*dy+dz*dz )/pGen.Beta()/pGen.Gamma();
    return ct;
 
}

// ========================================================================================
int AlbertoUtil::GetBestBstrange()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

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
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

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
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

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

    for( unsigned int i=0; i<tkJpsi.size(); ++i ){
       int j = tkJpsi[i];
       if(trkPt->at(j) < 4.0) return false;
       if(fabs(trkEta->at(j)) > 2.2) return false;
       TLorentzVector a;
       a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), MassMu );
       tJPsi += a;
    }

    if(tJPsi.Pt() < 8.0) return false;

    return true;
}

// ========================================================================================
int AlbertoUtil::GetBestBstrangeTight()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBsJPsiPhi) ) continue;
       if( svtMass->at(iB)<BsMassRange[0] || svtMass->at(iB)>BsMassRange[1] ) continue;

       if((svtDist2D->at(iB)/svtSigma2D->at(iB)) < 3) continue;
       if( ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) ) < 0.10 ) continue;

       int iJPsi = (subVtxFromSV(iB)).at(0);
       if(!IsTightJPsi(iJPsi)) continue;

       int iPhi = (subVtxFromSV(iB)).at(1);

       vector <int> tkSsB = tracksFromSV(iB);
       vector <int> tkJpsi = tracksFromSV(iJPsi);
       vector <int> tkPhi = tracksFromSV(iPhi);

       if(trkPt->at(tkPhi[0]) < 0.7) continue;
       if(trkPt->at(tkPhi[1]) < 0.7) continue;
       if(fabs(svtMass->at(iPhi) - MassPhi) > 0.01 ) continue;

       TLorentzVector tB(0,0,0,0);

       for( unsigned int i=0; i<tkSsB.size(); ++i ){

         int j = tkSsB[i];

         float m = MassK;

         if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MassMu;

         TLorentzVector a;
         a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
         tB += a;

       }

       if(tB.Pt() < 10.0) continue;

       if( svtChi2->at(iB)>bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBupTight()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBuJPsiK) ) continue;
       if( svtMass->at(iB)<BuMassRange[0] || svtMass->at(iB)>BuMassRange[1] ) continue;

       if((svtDist2D->at(iB)/svtSigma2D->at(iB)) < 3) continue;
       if( ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) ) < 0.10 ) continue;

       int iJPsi = (subVtxFromSV(iB)).at(0);
       if(!IsTightJPsi(iJPsi)) continue;

       vector <int> tkJpsi = tracksFromSV(iJPsi);
       vector <int> tkSsB = tracksFromSV(iB);

       TLorentzVector tB(0,0,0,0);
       float KaonPt = 0;

       for( unsigned int i=0; i<tkSsB.size(); ++i ){

         int j = tkSsB[i];

         float m = MassK;

         if( j == tkJpsi[0] || j == tkJpsi[1] ){ m = MassMu; }else{ KaonPt = trkPt->at(j); }

         TLorentzVector a;
         a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
         tB += a;

       }

       if(tB.Pt() < 10.0) continue;
       if(KaonPt < 1.6) continue;

       if( svtChi2->at(iB)>bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestBdownTight()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( unsigned short int iB=0; iB<nSVertices; ++iB ){

       if((svtType->at(iB)!=PDEnumString::svtBdJPsiKx) ) continue;
       if( svtMass->at(iB)<BdMassRange[0] || svtMass->at(iB)>BdMassRange[1] ) continue;

       if((svtDist2D->at(iB)/svtSigma2D->at(iB)) < 3) continue;
       if( ChiSquaredProbability( svtChi2->at(iB), svtNDOF->at(iB) ) < 0.10 ) continue;

       int iJPsi = (subVtxFromSV(iB)).at(0);
       if(!IsTightJPsi(iJPsi)) continue;
       if(fabs(svtMass->at(iJPsi) - MassJPsi) > 0.10 ) continue;

       vector <int> tkJpsi = tracksFromSV(iJPsi);
       vector <int> tkSsB = tracksFromSV(iB);

       TLorentzVector tB(0,0,0,0);

       for( unsigned int i=0; i<tkSsB.size(); ++i ){

         int j = tkSsB[i];

         float m = MassKx;

         if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MassMu;

         TLorentzVector a;
         a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
         tB += a;

       }

       if(tB.Pt() < 10.0) continue;

       if( svtChi2->at(iB)>bestChi2 ) continue;
       index = iB;
       bestChi2 = svtChi2->at(iB);

    }
    return index;
}

// ========================================================================================
int AlbertoUtil::GetBestJpsi()
{
    int index = -1;
    float bestChi2 = 1e9;
    for( unsigned short int i=0; i<nSVertices; ++i ){

       if((svtType->at(i)!=PDEnumString::svtJPsi) ) continue;
       if( abs(svtMass->at(i)-MassJPsi) > MassRangeJPsi) continue;

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
unsigned short int AlbertoUtil::TagMixStatus( unsigned int genindex )
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

    if((abs(muoEta->at( iMuon ))<1.2)&&(computeMvaBarrel(iMuon)>=wpB)) return true;
    if((abs(muoEta->at( iMuon ))>=1.2)&&(computeMvaEndcap(iMuon)>=wpE)) return true;

    return false;

}

// =====================================================================================
float AlbertoUtil::GetJetCharge(int iJet, float kappa)
{

    float QJet = 0;
    float ptJet = 0;

    vector <int> list = pfCandFromJet( iJet );

    for(std::vector<int>::iterator it = list.begin(); it != list.end(); ++it){

       float pt = pfcPt->at(*it);
       float eta = pfcEta->at(*it);

       if(pt<0.2) continue;
       if(abs(eta)>2.5) continue;

       QJet += pfcCharge->at(*it) * pow(pt, kappa);
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
float AlbertoUtil::GetMvaMuonValue(int iMuon)
{ 
    return (abs(muoEta->at( iMuon ))<1.2) ? computeMvaBarrel(iMuon) : computeMvaEndcap(iMuon); 
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

    float mean=MassBs;
    if(process=="BsJPsiPhi")   mean=MassBs;
    if(process=="BuJPsiK")     mean=MassBp;
    if(process=="BdJPsiKx")       mean=MassB0;
    if(process=="BdKxMuMu")       mean=MassB0;

    float sigma = 0.015;

    TString funcDef = "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)+[5]*exp([6]*x)";

    TF1 *func = new TF1("func", funcDef, 5.0, 5.5);

    func->FixParameter(1, mean);
    func->SetParameter(2, sigma);
    func->SetParameter(4, sigma);

    func->SetParameter(5, 1);
    func->SetParameter(6, -1);

    func->SetParLimits(0, 0, 1e6);
    func->SetParLimits(3, 0, 1e6);

    func->SetParLimits(1, mean-sigma, mean+sigma);
    func->SetParLimits(2, sigma/3, sigma*3);
    func->SetParLimits(4, sigma/3, sigma*3);

    func->SetParLimits(5, 0, 1e9);
    func->SetParLimits(6, -1e2, 0);

    hist->Fit("func","MRLQ");

    TF1 *fit = hist->GetFunction("func");

    float a1 = fit->GetParameter(0);
    float a2 = fit->GetParameter(3);
    float sigma1 = fit->GetParameter(2);
    float sigma2 = fit->GetParameter(4);
    float nEvt = TMath::Sqrt(2*TMath::Pi())*(a1*sigma1+a2*sigma2);
    nEvt/=hist->GetBinWidth(0);

    return nEvt;

}
// =====================================================================================
int AlbertoUtil::GetBestPV(int isvt, TLorentzVector t)
{
    int ssPV = -1;
    float bestCos = -1;

    TVector3 vSsBp(t.Px(),t.Py(),t.Pz());
    TVector3 vSVT( svtX->at(isvt), svtY->at(isvt), svtZ->at(isvt) );

    for( int i=0; i<nPVertices; ++i ){

       if(fabs(svtZ->at(isvt) - pvtZ->at( i )) > 1.0 ) continue;

       TVector3 vPV(pvtX->at( i ), pvtY->at( i ), pvtZ->at( i ) );
       TVector3 vPointing;
       vPointing = vSVT - vPV;
       float cos = vPointing.Unit() * vSsBp.Unit();

       if(cos > bestCos ){
         bestCos = cos;
         ssPV = i;
       }

    }

    return ssPV;
}
// =====================================================================================
float AlbertoUtil::GetSignedDxy(int iMuon, int iPV)
{

    float dxy = dXY( muonTrack( iMuon, PDEnumString::muInner ), pvtX->at(iPV), pvtY->at(iPV) );
    int sign = IPsign(iMuon, iPV);

    return dxy*sign;

}
// =====================================================================================
TLorentzVector AlbertoUtil::GetTLorentzVecFromJpsiX(int iSvt)
{
    int iJPsi = (subVtxFromSV(iSvt)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi) ;
    vector <int> tkSsB = tracksFromSV(iSvt);

    TLorentzVector t(0,0,0,0);

    for( unsigned int i=0; i<tkSsB.size(); ++i ){

       int j = tkSsB[i];

       float m = MassK;

       if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MassMu;

       TLorentzVector a ;
       a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m ) ;
       t += a ;

    }

    return t;
}
// =====================================================================================
int AlbertoUtil::GetCandidate(TString process, bool useTightSel)
{
    if(process=="BsJPsiPhi"){
       if(useTightSel){
         return GetBestBstrangeTight();
       }else{
         return GetBestBstrange();
       }
    }

    if(process=="BuJPsiK"){
       if(useTightSel){
         return GetBestBupTight();
       }else{
         return GetBestBup();
       }
    }

    return -1;
}
// =====================================================================================
float AlbertoUtil::dZ(int itk, int iPV)
{
    return PDAnalyzerUtil::dZ(itk, pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
}
    
