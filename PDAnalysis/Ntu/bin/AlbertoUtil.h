#ifndef AlbertoUtil_H
#define AlbertoUtil_H

#include <vector>
#include <string>
#include "PDAnalysis/Ntu/interface/PDEnumString.h"
#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "PDAnalyzerUtil.h"
#include "PDSoftMuonMvaEstimator.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "Math/MinimizerOptions.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;


class AlbertoUtil:  public virtual PDAnalyzerUtil
,                   public virtual PDSoftMuonMvaEstimator
,                   public virtual PDGenHandler{

public:

    AlbertoUtil();
    virtual ~AlbertoUtil();

    const double MUMASS  = 0.1056583745;
    const double EMASS   = 0.00051099895;
    const double PMASS   = 0.93827208816;

    const double PIMASS  = 0.13957039;
    const double PI0MASS = 0.1349768;

    const double KMASS   = 0.493677;
    const double K0MASS  = 0.497611;
    const double KXMASS  = 0.89166;

    const double BSMASS  = 5.36688;
    const double B0MASS  = 5.27965;
    const double BUMASS  = 5.27934;
    const double BCMASS  = 6.2749;

    const double JPSIMASS = 3.0969;
    const double PHIMASS  = 1.019461;

    const double LAMBDAMASS = 5.61960;

    const unsigned int listLundBmesons[24] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 
                                            20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 
                                            10541, 543, 10543, 20543, 545};

    const unsigned int listLundBbaryons[35] = {5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 
                                            5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 
                                            5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 
                                            5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554};


    const unsigned int listLundBottonium[29] = {551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 
                                            20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 
                                            220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 
                                            110555, 120555, 200555, 557, 100557};

    const unsigned int LongLivedList[7] = {11,13,211,321,2212}; //e, mu, pi, K, p   

    const unsigned int listLundCmesons[18] = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 
                                            415, 425, 431, 10431, 433, 10433, 20433, 435};

    const unsigned int listLundCbaryons[22] = {4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 
                                            4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

    const unsigned int listLundCharmonium[13] = {441, 10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443,
                                            9010443, 9020443, 445, 100445};

    double MassRangeJPsi = 0.15;

    double BsMassRange[2] = {5.0, 6.0};
    double BuMassRange[2] = {5.0, 6.0};
    double BdMassRange[2] = {5.0, 6.0};

    //B tight cuts JpsiMu / JpsiTrkTrk
    double bPtCut    = 11;    //10
    double bCtCut    = 0.007; //0.02
    double bVprobCut = 0.001;
    double bMuPtCut  = 3.5;   //3.6
    double bMuEtaCut = 2.4;
    double bKPtCut   = 1.2;   //1.0
    double bKEtaCut  = 2.5;

    void SetBctCut(double newValue = 0.007){ bCtCut = newValue; }
    void SetBptCut(double newValue = 10){ bPtCut = newValue; }
    void SetBmuptCut(double newValue = 3.5){ bMuPtCut = newValue; }
    void SetBkptCut(double newValue = 1.2){ bKPtCut = newValue; }

    void SetJpsiMuCut();
    void SetJpsiTrkTrkCut();

    void SetMassRangeJPsi(double newValue){ MassRangeJPsi = newValue; }
    void SetBsMassRange(double lower, double upper) { BsMassRange[0] = lower; BsMassRange[1] = upper; }
    void SetBuMassRange(double lower, double upper) { BuMassRange[0] = lower; BuMassRange[1] = upper; }
    void SetBdMassRange(double lower, double upper) { BdMassRange[0] = lower; BdMassRange[1] = upper; }

    bool IsLongLived( unsigned int i );
    bool IsB( unsigned int i ) ;
    bool IsC( unsigned int i ) ;
    bool IsCharmonium( unsigned int i ) ;
    bool IsBottomium( unsigned int i ) ;

    int GetClosestGen( double eta, double phi, double pt );
    int GetClosestGenLongLivedB( double eta, double phi, double pt, std::vector <int> *GenList );
    int GetOverlappedTrack( int trk, std::vector <int> *List );
    bool AreOverlapped( double pt1, double eta1, double phi1, double pt2, double eta2, double phi2 );
    int GetAncestor( unsigned int iGen, std::vector <int> *GenList );
    int MuonFromTrack(int trk);
    double GetGenCT( unsigned int genIndex );

    int GetBestBstrange();
    int GetBestBdown();
    int GetBestBup();
    int GetBestBstrangeTight();
    int GetBestBdownTight();
    int GetBestBupTight();
    bool IsTightJPsi(int iJPsi);
    bool IsTightPhi(int iPhi);
    int GetBestJpsi();
    int GetTightCandidate(TString process);
    int GetCandidate(TString process);
    
    double GetInvMass(int i1, int i2, double mass1, double mass2);
    int TagMixStatus( unsigned int genIndex );
    double GetMuoPFiso (int iMuon);
    double GetJetCharge(int iJet, double kappa);
    double GetListCharge(std::vector <int> *list, double kappa);
    bool IsMvaMuon(int iMuon, double wp);
    int IPsign(int iMuon, int iPV);
    double GetJetProbb(int iJet);
    double CountEventsWithFit(TH1 *hist, TString process);
    int GetBestPV(int isvt, TLorentzVector t);
    TLorentzVector GetTLorentzVecFromJpsiX(int iSvt);
    double dZ(int itk, int iPV);
    double dXYjet(int itk, int iPV, int iJet);
    void PrintMotherChain(int iGen);
    void PrintDaughterTree(int iGen, const std::string & pre);
    void PrintDaughterTreePt(int iGen, const std::string & pre);
    bool HasDaughter(int iGen);

    double GetCt2D(TLorentzVector t, int iSV);
    double GetCt2D(TLorentzVector t, int iSV, int iPV);
    double GetCt3D(TLorentzVector t, int iSV, int iPV);

    double GetCt2DErr(TLorentzVector t, int iSV);
    double GetCt2DErr(TLorentzVector t, int iSV, int iPV);
    double GetCt3DErr(TLorentzVector t, int iSV, int iPV);

    bool isTrkHighPurity(int itk){ return (( trkQuality->at( itk ) >> 2 ) & 1); }


protected:


};

#endif