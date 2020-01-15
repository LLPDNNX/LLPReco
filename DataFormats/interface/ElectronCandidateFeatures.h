#ifndef LLPReco_DataFormats_ElectronCandidateFeatures_h
#define LLPReco_DataFormats_ElectronCandidateFeatures_h

 
namespace llpdnnx {
 
class ElectronCandidateFeatures { 

  public:

  float elec_pt ;
  float elec_jetPtRatio ;
  float elec_jetDeltaR ; 
  float elec_p ;
  float elec_eta;
  float elec_phi;
  float elec_charge ; 
  float elec_energy;
  float elec_EtFromCaloEn ;
  float elec_isEB ; 
  float elec_isEE ; 
  float elec_ecalEnergy ; 
  float elec_isPassConversionVeto ;
  float elec_convDist ; 
  int   elec_convFlags ; 
  float elec_convRadius ; 
  float elec_hadronicOverEm ;
  float elec_ecalDrivenSeed;
 
// superCluster block
 
  float elecSC_energy ; 
  float elecSC_eta ; 
  float elecSC_phi ;
  float elecSC_et ;
  float elecSC_eSuperClusterOverP ; 
  float elec_scE1x5 ; 
  float elec_scE2x5Max ; 
  float elec_scE5x5 ; 
  float elec_scPixCharge ; 
  float elec_scSigmaEtaEta ;
  float elec_scSigmaIEtaIEta ;  
  float elec_superClusterFbrem ; 

  float elec_2dIP ; 
  float elec_2dIPSig ;
  float elec_3dIP ; 
  float elec_3dIPSig ; 
  float eSeedClusterOverP ;
  float eSeedClusterOverPout;
  float eSuperClusterOverP;
  float eTop; 
 
  float elec_deltaEtaEleClusterTrackAtCalo ; 
  float elec_deltaEtaSeedClusterTrackAtCalo ;
  float elec_deltaPhiSeedClusterTrackAtCalo ; 
  float elec_deltaEtaSeedClusterTrackAtVtx ; 
  float elec_deltaEtaSuperClusterTrackAtVtx ;
  float elec_deltaPhiEleClusterTrackAtCalo ; 
  float elec_deltaPhiSuperClusterTrackAtVtx ;
  float elec_sCseedEta ;  

// electron gsf variables. 
  float elec_EtaRel ; 
  float elec_dxy ; 
  float elec_dz ;
  float elec_nbOfMissingHits ; 
  float elec_gsfCharge ;


// 5*5 cells cluster 
//

  float elec_full5x5_sigmaIetaIeta ;
  float elec_full5x5_e1x5 ;
  float elec_full5x5_e2x5Bottom ;
  float elec_full5x5_e2x5Left ;
  float elec_full5x5_e2x5Max ;
  float elec_full5x5_e2x5Right ;
  float elec_full5x5_e2x5Top ;
  float elec_full5x5_e5x5 ;
  float elec_full5x5_eBottom ;
  float elec_full5x5_eLeft;
  float elec_full5x5_eRight;
  float elec_full5x5_eTop;
  float elec_full5x5_hcalDepth1OverEcal ;
  float elec_full5x5_hcalDepth1OverEcalBc ;
  float elec_full5x5_hcalDepth2OverEcal;
  float elec_full5x5_hcalDepth2OverEcalBc ;
  float elec_full5x5_hcalOverEcal ;
  float elec_full5x5_hcalOverEcalBc;   
  float elec_full5x5_r9 ;
 
  int   elec_numberOfBrems ;
  float elec_trackFbrem ; 
  float elec_fbrem ; 
  float elec_e2x5Max ; 
  float elec_e1x5 ; 
  float elec_e5x5 ;
 
// Isolation block
  float elec_neutralHadronIso; 
  float elec_particleIso  ; 
  float elec_photonIso ;
  float elec_puChargedHadronIso ; 
  float elec_trackIso ;  
  float elec_hcalDepth1OverEcal ; 
  float elec_hcalDepth2OverEcal ; 
  float elec_ecalPFClusterIso ;
  float elec_hcalPFClusterIso ;  
  float elec_dr03TkSumPt ; 
  float elec_dr03HcalDepth2TowerSumEt ; 
  float elec_dr03EcalRecHitSumEt ; 
  float elec_dr03HcalDepth1TowerSumEt ;  
  float elec_dr03HcalDepth1TowerSumEtBc ; 
  float elec_pfSumPhotonEt ; 
  float elec_pfSumChargedHadronPt ; 
  float elec_pfSumNeutralHadronEt ; 
  float elec_pfSumPUPt ;



// cone 0.4 
//

float elec_dr04EcalRecHitSumEt ;  
float elec_dr04HcalDepth1TowerSumEt ;  
float elec_dr04HcalDepth1TowerSumEtBc ;
float elec_dr04HcalDepth2TowerSumEt ; 
float elec_dr04HcalDepth2TowerSumEtBc  ;
float elec_dr04HcalTowerSumEt  ;
float elec_dr04HcalTowerSumEtBc  ;


 
 
};

 
}

#endif
