#ifndef LLPReco_DataFormats_ElectronCandidateFeatures_h
#define LLPReco_DataFormats_ElectronCandidateFeatures_h

namespace llpdnnx {
 
class ElectronCandidateFeatures { 

  public:
    int elec_jetIdx;
    float elec_ptrel;
    float elec_deltaR;
    float elec_deta;
    float elec_dphi;
    float elec_px;
    float elec_py;
    float elec_pz;
    float elec_charge;
    float elec_energy;
    float elec_EtFromCaloEn;
    int elec_isEB; 
    int elec_isEE; 
    float elec_ecalEnergy;
    int elec_isPassConversionVeto;
    float elec_convDist;
    int   elec_convFlags;
    float elec_convRadius;
    float elec_hadronicOverEm;
    float elec_ecalDrivenSeed;
    float elec_2dIP;
    float elec_2dIPSig;
    float elec_3dIP;
    float elec_3dIPSig;

    float elecSC_energy;
    float elecSC_deta;
    float elecSC_dphi;
    float elecSC_et;
    float elecSC_eSuperClusterOverP;
    float elec_scPixCharge;
    float elec_superClusterFbrem;

    float elec_eSeedClusterOverP;
    float elec_eSeedClusterOverPout;
    float elec_eSuperClusterOverP;

    // shower shape
    float elec_sigmaEtaEta;
    float elec_sigmaIetaIeta;
    float elec_sigmaIphiIphi;
    float elec_e5x5;
    float elec_e5x5Rel;
    float elec_e1x5Overe5x5;
    float elec_e2x5MaxOvere5x5;
    float elec_r9;
    float elec_hcalOverEcal;
    float elec_hcalDepth1OverEcal;
    float elec_hcalDepth2OverEcal;

    // Track-Cluster Matching Attributes
    float elec_deltaEtaEleClusterTrackAtCalo;
    float elec_deltaEtaSeedClusterTrackAtCalo;
    float elec_deltaPhiSeedClusterTrackAtCalo; 
    float elec_deltaEtaSeedClusterTrackAtVtx;
    float elec_deltaEtaSuperClusterTrackAtVtx;
    float elec_deltaPhiEleClusterTrackAtCalo;
    float elec_deltaPhiSuperClusterTrackAtVtx;

    float elec_sCseedEta;

    // electron gsf variables. 
    float elec_EtaRel;
    float elec_dxy;
    float elec_dxyError;
    float elec_dxySig;
    float elec_dz;
    float elec_dzError;
    float elec_dzSig;
    int elec_nbOfMissingHits;
    float elec_gsfCharge;
    int   elec_ndof;
    float elec_chi2;
    int   elec_numberOfBrems;
    float elec_trackFbrem;
    float elec_fbrem;

    // Isolation block
    float elec_neutralHadronIso;
    float elec_particleIso;
    float elec_photonIso;
    float elec_puChargedHadronIso;
    float elec_trackIso;
    float elec_ecalPFClusterIso;
    float elec_hcalPFClusterIso;
    
    float elec_pfSumPhotonEt;
    float elec_pfSumChargedHadronPt; 
    float elec_pfSumNeutralHadronEt;
    float elec_pfSumPUPt;

    float elec_dr04TkSumPt;
    float elec_dr04EcalRecHitSumEt;
    float elec_dr04HcalDepth1TowerSumEt;
    float elec_dr04HcalDepth1TowerSumEtBc;
    float elec_dr04HcalDepth2TowerSumEt;
    float elec_dr04HcalDepth2TowerSumEtBc;
    float elec_dr04HcalTowerSumEt;
    float elec_dr04HcalTowerSumEtBc;
};

}

#endif
