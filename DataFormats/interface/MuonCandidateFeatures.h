#ifndef LLPReco_DataFormats_MuonCandidateFeatures_h
#define LLPReco_DataFormats_MuonCandidateFeatures_h

 
namespace llpdnnx {
 
struct MuonCandidateFeatures
{ 

    int  mu_jetIdx ; 
    bool  mu_isGlobal ; 
    bool  mu_isTight ; 
    bool  mu_isMedium ; 
    bool  mu_isLoose ; 
    bool  mu_isStandAlone ;

    
    float mu_ptrel ;
    float mu_deta;
    float mu_dphi;
    float mu_charge ; 
    float mu_energy;
    float mu_et ;
    float mu_jetDeltaR ; 
    float mu_numberOfMatchedStations ;

    float mu_2dIp ; 
    float mu_2dIpSig ;
    float mu_3dIp ; 
    float mu_3dIpSig ; 

    float mu_EtaRel;
    float mu_dxy ; 
    float mu_dxyError ; 
    float mu_dxySig ; 
    float mu_dz ; 
    float mu_dzError ; 
    float mu_numberOfValidPixelHits; 
    float mu_numberOfpixelLayersWithMeasurement ; 
    float mu_numberOfstripLayersWithMeasurement ; //that does not help. needs to be discussed.

    float mu_chi2 ; 
    float mu_ndof ; 

    float mu_caloIso ; 
    float mu_ecalIso ; 
    float mu_hcalIso ;

    float mu_sumPfChHadronPt ; 
    float mu_sumPfNeuHadronEt ; 
    float mu_Pfpileup ; 
    float mu_sumPfPhotonEt ; 

    float mu_sumPfChHadronPt03 ; 
    float mu_sumPfNeuHadronEt03 ; 
    float mu_Pfpileup03 ; 
    float mu_sumPfPhotonEt03 ; 

    float mu_sumChHadronPt ; 
    float mu_sumNeuHadronEt ; 
    float mu_pileup ; 
    float mu_sumPhotonEt ; 

    float mu_timeAtIpInOut ; 
    float mu_timeAtIpInOutErr ; 
    float mu_timeAtIpOutIn ; 

};

 
}

#endif
