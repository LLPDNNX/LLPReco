#ifndef LLPReco_DataFormats_MuonCandidateFeatures_h
#define LLPReco_DataFormats_MuonCandidateFeatures_h

namespace llpdnnx {

struct MuonCandidateFeatures {

    int mu_jetIdx;
    int mu_isGlobal;
    int mu_isTight;
    int mu_isMedium;
    int mu_isLoose;
    int mu_isStandAlone;

    float mu_ptrel;
    float mu_deta;
    float mu_dphi;
    float mu_charge;
    float mu_energy;
    float mu_et;
    float mu_jetDeltaR;
    int mu_numberOfMatchedStations;

    float mu_2dIP;
    float mu_2dIPSig;
    float mu_3dIP;
    float mu_3dIPSig;

    float mu_EtaRel;
    float mu_dxy;
    float mu_dxyError;
    float mu_dxySig;
    float mu_dz;
    float mu_dzError;
    int mu_numberOfValidPixelHits;
    int mu_numberOfpixelLayersWithMeasurement;
    int mu_numberOfstripLayersWithMeasurement; //that does not help. needs to be discussed.

    float mu_chi2;
    int    mu_ndof;

    float mu_caloIso;
    float mu_ecalIso;
    float mu_hcalIso;

    float mu_sumPfChHadronPt;
    float mu_sumPfNeuHadronEt;
    float mu_Pfpileup;
    float mu_sumPfPhotonEt;

    float mu_sumPfChHadronPt03;
    float mu_sumPfNeuHadronEt03;
    float mu_Pfpileup03;
    float mu_sumPfPhotonEt03;


    float mu_timeAtIpInOut;
    float mu_timeAtIpInOutErr;
    float mu_timeAtIpOutIn;


};


}

#endif
