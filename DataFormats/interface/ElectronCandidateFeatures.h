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
  float elec_et ;
  float elec_EtFromCaloEn ; 
  float elecSC_energy ; 
  float elecSC_eta ; 
  float elecSC_phi ;

  float elec_isEB ; 
  float elec_isEE ; 
  float elec_eSuperClusterOverP ; 
  float elec_ecalEnergy ; 
  float elec_isPassConversionVeto ; 

  float elec_2dIP ; 
  float elec_2dIPSig ;
  float elec_3dIP ; 
  float elec_3dIPSig ; 
 
  float elec_CutVeto ; 
  float elec_CutLoose ; 
  float elec_CutMedium ; 
  float elec_CutTight ;

// electron gsf variables. 
//
  float elec_dxy ; 
  float elec_dz ;
  float elec_nbOfMissingHits ; 
  float elec_gsfCharge ; 
 
//float elec_absdxyError ; 
//float elec_absdxySig ; 
//float elec_absdzError ; 
  
  float elec_fbrem ; 
  float elec_e2x5Max ; 
  float elec_e1x5 ; 
  float elec_e5x5 ; 

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
