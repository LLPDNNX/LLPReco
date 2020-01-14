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
  float elec_hadronicOverEm ;
 
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
 
  float elec_deltaEtaEleClusterTrackAtCalo ; 
  float elec_deltaEtaSeedClusterTrackAtCalo ; 
  float elec_deltaEtaSeedClusterTrackAtVtx ; 
  float elec_deltaEtaSuperClusterTrackAtVtx ; 
  float elec_deltaPhiEleClusterTrackAtCalo ; 
  float elec_deltaPhiSeedClusterTrackAtCalo ; 
  float elec_deltaPhiSuperClusterTrackAtVtx ; 

// electron gsf variables. 
//
  float elec_dxy ; 
  float elec_dz ;
  float elec_nbOfMissingHits ; 
  float elec_gsfCharge ;
  float elec_full5x5_sigmaIetaIeta ; 
 
  int elec_numberOfBrems ; 
  float elec_fbrem ; 
  float elec_e2x5Max ; 
  float elec_e1x5 ; 
  float elec_e5x5 ;
 
// Isolation block 
  float elec_hcalPFClusterIso ;  
  float elec_dr03TkSumPt ; 
  float elec_hcalDepth1OverEcal ; 
  float elec_hcalDepth2OverEcal ; 
  float elec_dr03HcalDepth2TowerSumEt ; 
  float elec_hcalDepth2TowerSumEtNoVeto ; 
  float elec_hcalDepth1TowerSumEtNoVeto ; 
  float elec_dr03EcalRecHitSumEt ; 
  float elec_dr03HcalDepth1TowerSumEt ;  
  float elec_dr03HcalDepth1TowerSumEtBc ; 
  float elec_pfSumPhotonEt ; 
  float elec_pfSumChargedHadronPt ; 
  float elec_pfSumNeutralHadronEt ; 
  float elec_pfSumPUPt ; 
 
};

 
}

#endif
