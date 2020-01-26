#ifndef LLPReco_DataFormats_ElectronCandidateFeatures_h
#define LLPReco_DataFormats_ElectronCandidateFeatures_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
namespace llpdnnx {
 
class ElectronCandidateFeatures { 

  public:

 

    int elec_jetIdx;
    float elec_ptrel;
    float elec_jetDeltaR ; 
    float elec_deta;
    float elec_dphi;
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
    float elecSC_deta ; 
    float elecSC_dphi ;
    float elecSC_et ;
    float elecSC_eSuperClusterOverP ; 
    float elec_scE1x5Overe5x5 ; 
    float elec_scE2x5MaxOvere5x5 ; 
    float elec_scE5x5 ; 
    float elec_scE5x5Rel ; 
    float elec_scPixCharge ; 
    float elec_scSigmaEtaEta ;
    float elec_scSigmaIEtaIEta ;  
    float elec_superClusterFbrem ; 

    float elec_2dIP ; 
    float elec_2dIPSig ;
    float elec_3dIP ; 
    float elec_3dIPSig ; 
    float elec_eSeedClusterOverP;
    float elec_eSeedClusterOverPout;
    float elec_eSuperClusterOverP;
    float elec_eTopOvere5x5; 

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


    int   elec_numberOfBrems ;
    float elec_trackFbrem ; 
    float elec_fbrem ; 
    float elec_e5x5 ; 
    float elec_e5x5Rel ; 
    float elec_e1x5Overe5x5 ; 
    float elec_e2x5MaxOvere5x5;

    // 5*5 cells cluster 
    //

    float elec_full5x5_e5x5 ;
    float elec_full5x5_e5x5Rel ; 
    float elec_full5x5_sigmaIetaIeta ;
    float elec_full5x5_e1x5Overe5x5 ;
    float elec_full5x5_e2x5BottomOvere5x5 ;
    float elec_full5x5_e2x5LeftOvere5x5 ;
    float elec_full5x5_e2x5MaxOvere5x5 ;
    float elec_full5x5_e2x5RightOvere5x5 ;
    float elec_full5x5_e2x5TopOvere5x5 ;



    float elec_full5x5_eBottomOvere5x5 ;
    float elec_full5x5_eLeftOvere5x5;
    float elec_full5x5_eRightOvere5x5;
    float elec_full5x5_eTopOvere5x5;
    float elec_full5x5_hcalDepth1OverEcal ;
    float elec_full5x5_hcalDepth1OverEcalBc ;
    float elec_full5x5_hcalDepth2OverEcal;
    float elec_full5x5_hcalDepth2OverEcalBc ;
    float elec_full5x5_hcalOverEcal ;
    float elec_full5x5_hcalOverEcalBc;   
    float elec_full5x5_r9 ;



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

    float elec_dr03EcalRecHitSumEt ; 
    float elec_dr03HcalDepth1TowerSumEt ;  
    float elec_dr03HcalDepth1TowerSumEtBc ; 
    float elec_dr03HcalDepth2TowerSumEt ; 
    float elec_dr03HcalDepth2TowerSumEtBc ; 
    float elec_pfSumPhotonEt ; 
    float elec_pfSumChargedHadronPt ; 
    float elec_pfSumNeutralHadronEt ; 
    float elec_pfSumPUPt ;

    float elec_dr04EcalRecHitSumEt ;  
    float elec_dr04HcalDepth1TowerSumEt ;  
    float elec_dr04HcalDepth1TowerSumEtBc ;
    float elec_dr04HcalDepth2TowerSumEt ; 
    float elec_dr04HcalDepth2TowerSumEtBc  ;
    float elec_dr04HcalTowerSumEt  ;
    float elec_dr04HcalTowerSumEtBc  ;

void electronFeatures(const pat::Electron& electron , const pat::Jet& jet , const reco::Vertex& pv ) ;  
 
};

 
}

#endif
