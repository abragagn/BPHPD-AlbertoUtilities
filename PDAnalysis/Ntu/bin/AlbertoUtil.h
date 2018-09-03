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

class AlbertoUtil: 	public virtual PDAnalyzerUtil
,										public virtual PDSoftMuonMvaEstimator
,										public virtual PDGenHandler{

public:

	AlbertoUtil();
	~AlbertoUtil();

	const float MassJPsi = 3.097;
	const float MassMu = 0.1057;
	const float MassK = 0.4937;
	const float MassKx = 0.8917;
	const float MassBs = 5.3663;
	const float MassBp = 5.2793;
	const float MassB0 = 5.2796;
	const float MassKs = 0.4977;
	const float MassPi = 0.1396;
	const float MassPhi = 1.019;
	const float MassProton = 0.9382;
	const float MassLambda = 1.116;

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

	float MassRangeJPsi = 0.15;

	float BsMassRange[2] = {5.0, 6.0};
	float BuMassRange[2] = {5.0, 6.0};
	float BdMassRange[2] = {5.0, 6.0};

	TString weightFileBarrel = "/lustre/cmswork/abragagn/weights/2016/TMVAClassification_DNNGlobalBarrel2016wIPwIso.weights.xml";
	TString weightFileEndcap = "/lustre/cmswork/abragagn/weights/2016/TMVAClassification_DNNGlobalEndcap2016wIPwIso.weights.xml";

	void SetMassRangeJPsi(float newValue){ MassRangeJPsi = newValue; }
	void SetBsMassRange(float lower, float upper) {	BsMassRange[0] = lower; BsMassRange[1] = upper; }
	void SetBpMassRange(float lower, float upper) {	BuMassRange[0] = lower; BuMassRange[1] = upper; }
	void SetBdMassRange(float lower, float upper) {	BdMassRange[0] = lower; BdMassRange[1] = upper; }

	bool IsLongLived( unsigned int i );
	bool IsB( unsigned int i ) ;
	bool IsC( unsigned int i ) ;
	bool IsCharmonium( unsigned int i ) ;
	bool IsBottomium( unsigned int i ) ;

	int GetClosestGen( float eta, float phi, float pt );
	int GetClosestGenLongLivedB( float eta, float phi, float pt, std::vector <int> *GenList );
	int GetAncestor( unsigned int iGen, std::vector <int> *GenList );
	int WhichMuon(int trk);
	float GetCT( unsigned int genIndex );

	int GetBestBstrange();
	int GetBestBdown();
	int GetBestBup();
	int GetBestBstrangeTight();
	int GetBestBdownTight();
	int GetBestBupTight();
	bool IsTightJPsi(int iJPsi);
	int GetBestJpsi();
	
	float GetInvMass(int i1, int i2, float mass1, float mass2);
	unsigned short int TagMixStatus( unsigned int genIndex );
	float GetMuoPFiso (int iMuon);
	float GetJetCharge(int iJet, float kappa);
	bool isMvaMuon(int iMuon, float wpB, float wpE);
	float GetMvaMuonValue(int iMuon);
	int IPsign(int iMuon);
	int IPsign(int iMuon, int iPV);
	float GetJetProbb(int iJet);
	float CountEventsWithFit(TH1 *hist, TString process);
	int GetBestPV(int isvt, TLorentzVector t);
	float GetSignedDxy(int iMuon, int iPV);


protected:


};

#endif