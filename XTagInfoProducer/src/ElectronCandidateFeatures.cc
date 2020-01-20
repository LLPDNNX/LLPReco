 
#include "LLPReco/DataFormats/interface/ElectronCandidateFeatures.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"


void llpdnnx::ElectronCandidateFeatures::electronFeatures(const pat::Electron& electron  , const pat::Jet& jet , const reco::Vertex& pv ) {


 
                elec_ptrel = electron.pt()/jet.pt() ; 
                elec_eta = electron.eta() ; 
                elec_phi = electron.phi() ; 
                elec_charge = electron.charge() ; 
                elec_energy = electron.energy()/electron.pt(); 
                elec_jetDeltaR = reco::deltaR(electron , jet) ; 
                elec_EtFromCaloEn = electron.caloEnergy() * sin(electron.p4().theta());
                elec_ecalDrivenSeed = electron.ecalDrivenSeed() ;

                elec_isEB = electron.isEB() ;  
                elec_isEE  = electron.isEE();
                elec_ecalEnergy  = electron.ecalEnergy()/electron.pt();
                elec_isPassConversionVeto = electron.passConversionVeto();

                elec_convDist = electron.convDist() ; 
                elec_convFlags = electron.convFlags() ; 
                elec_convRadius = electron.convRadius() ; 


                elec_3dIP = electron.dB(pat::Electron::PV3D) ; 
                elec_3dIPSig = electron.dB(pat::Electron::PV3D); 
                elec_2dIP = electron.dB() ; 
                elec_2dIPSig = electron.dB()/electron.edB() ; 
                elec_sCseedEta = electron.superCluster()->seed()->eta();


                elec_numberOfBrems  = electron.numberOfBrems () ;
                elec_trackFbrem  = electron.trackFbrem() ; 
                elec_fbrem = electron.fbrem() ; 


                elec_e5x5 = electron.e5x5() ;
                elec_e5x5Rel = electron.e5x5()/jet.pt() ;
                elec_e1x5Overe5x5 = electron.e1x5()/electron.e5x5() ; 
                elec_e2x5MaxOvere5x5 = electron.e2x5Max()/electron.e5x5() ; 

                elec_eSeedClusterOverP = electron.eSeedClusterOverP() ;
                elec_eSeedClusterOverPout = electron.eSeedClusterOverPout() ; 
                elec_eSuperClusterOverP = electron.eSuperClusterOverP() ; 
                elec_eTopOvere5x5 = electron.eTop()/electron.e5x5() ;  

                elec_hadronicOverEm = electron.hadronicOverEm() ;  
                elec_full5x5_sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();

                elec_full5x5_e5x5  = electron.full5x5_e5x5() ;
                elec_full5x5_e5x5Rel  = electron.full5x5_e5x5()/jet.pt() ;

                elec_full5x5_e1x5Overe5x5  = electron.full5x5_e1x5()/electron.full5x5_e5x5() ;
                elec_full5x5_e2x5BottomOvere5x5  = electron.full5x5_e2x5Bottom()/ electron.full5x5_e5x5();
                elec_full5x5_e2x5LeftOvere5x5  = electron.full5x5_e2x5Left()/ electron.full5x5_e5x5();
                elec_full5x5_e2x5MaxOvere5x5  = electron.full5x5_e2x5Max()/ electron.full5x5_e5x5();
                elec_full5x5_e2x5RightOvere5x5  = electron.full5x5_e2x5Right()/ electron.full5x5_e5x5();
                elec_full5x5_e2x5TopOvere5x5  = electron.full5x5_e2x5Top()/ electron.full5x5_e5x5();
                elec_full5x5_eBottomOvere5x5  = electron.full5x5_eBottom()/ electron.full5x5_e5x5();
                elec_full5x5_eLeftOvere5x5 = electron.full5x5_eLeft()/ electron.full5x5_e5x5();
                elec_full5x5_eRightOvere5x5 = electron.full5x5_eRight()/ electron.full5x5_e5x5();
                elec_full5x5_eTopOvere5x5 = electron.full5x5_eTop()/ electron.full5x5_e5x5();
                elec_full5x5_hcalDepth1OverEcal  = electron.full5x5_hcalDepth1OverEcal() ;
                elec_full5x5_hcalDepth1OverEcalBc  = electron.full5x5_hcalDepth1OverEcalBc() ;
                elec_full5x5_hcalDepth2OverEcal = electron.full5x5_hcalDepth2OverEcal() ;
                elec_full5x5_hcalDepth2OverEcalBc  = electron.full5x5_hcalDepth2OverEcalBc() ;
                elec_full5x5_hcalOverEcal  = electron.full5x5_hcalOverEcal() ;
                elec_full5x5_hcalOverEcalBc = electron.full5x5_hcalOverEcalBc() ;   
                elec_full5x5_r9  = electron.full5x5_r9() ;


                elec_deltaEtaEleClusterTrackAtCalo  = electron.deltaEtaEleClusterTrackAtCalo();
                elec_deltaPhiEleClusterTrackAtCalo = electron.deltaPhiEleClusterTrackAtCalo() ;

                elec_deltaEtaSeedClusterTrackAtCalo = electron.deltaEtaSeedClusterTrackAtCalo () ; 
                elec_deltaPhiSeedClusterTrackAtCalo = electron.deltaPhiSeedClusterTrackAtCalo() ; 

                elec_deltaEtaSeedClusterTrackAtVtx = electron.deltaEtaSeedClusterTrackAtVtx(); 
                elec_deltaEtaSuperClusterTrackAtVtx = electron.deltaEtaSuperClusterTrackAtVtx() ;  
                elec_deltaPhiSuperClusterTrackAtVtx = electron.deltaPhiSuperClusterTrackAtVtx ()  ;


                reco::Candidate::Vector electronMom = electron.gsfTrack()->momentum();

                elec_EtaRel =reco::btau::etaRel(jet.momentum().Unit(), electronMom); 
                elec_dxy = electron.gsfTrack()->dxy(pv.position()) ; 
                elec_dz = electron.gsfTrack()->dz(pv.position()) ;
                elec_nbOfMissingHits = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) ; 
                elec_gsfCharge = electron.gsfTrack()->charge() ;

                elecSC_energy = electron.superCluster()->energy()/electron.pt(); 
                elecSC_deta = std::fabs(electron.superCluster()->eta()-electron.gsfTrack()->eta());
                elecSC_dphi = reco::deltaPhi(electron.superCluster()->phi(),electron.gsfTrack()->phi());
                elecSC_et = electron.superCluster()->energy() * sin(electron.p4().theta())/electron.pt();
                elec_scPixCharge = electron.scPixCharge() ; 
                elec_scSigmaEtaEta = electron.scSigmaEtaEta() ; 
                elec_scSigmaIEtaIEta = electron.scSigmaIEtaIEta() ;
                elec_superClusterFbrem = electron.superClusterFbrem() ; 

                elec_scE5x5 = electron.scE5x5 () ; 
                elec_scE5x5Rel = electron.scE5x5()/jet.pt() ; 
                elec_scE1x5Overe5x5 = electron.scE1x5 ()/electron.scE5x5 () ; 
                elec_scE2x5MaxOvere5x5  = electron.scE2x5Max()/electron.scE5x5 () ; 
                elecSC_eSuperClusterOverP  = electron.eSuperClusterOverP();
                
                elec_particleIso  = electron.particleIso()/electron.pt(); 
                elec_neutralHadronIso  = electron.neutralHadronIso()/electron.pt();
                elec_photonIso = electron.photonIso()/electron.pt(); 
                elec_puChargedHadronIso = electron.puChargedHadronIso()/electron.pt(); 

                elec_trackIso = electron.trackIso()/electron.pt();
                elec_hcalDepth1OverEcal = electron.hcalDepth1OverEcal() ; 
                elec_hcalDepth2OverEcal = electron.hcalDepth2OverEcal() ;  
                elec_ecalPFClusterIso = electron.ecalPFClusterIso()/electron.pt(); 
                elec_hcalPFClusterIso = electron.hcalPFClusterIso()/electron.pt(); 

                elec_dr03TkSumPt = electron.dr03TkSumPt()/electron.pt();
                
                elec_dr03EcalRecHitSumEt = electron.dr03EcalRecHitSumEt()/electron.pt(); 
                
                elec_dr03HcalDepth1TowerSumEt = electron.dr03HcalDepth1TowerSumEt()/electron.pt(); 
                elec_dr03HcalDepth1TowerSumEtBc = electron.dr03HcalDepth1TowerSumEtBc()/electron.pt();
                elec_dr03HcalDepth2TowerSumEt = electron.dr03HcalDepth2TowerSumEt()/electron.pt(); 
                elec_dr03HcalDepth2TowerSumEtBc = electron.dr03HcalDepth2TowerSumEtBc()/electron.pt(); 

                elec_pfSumPhotonEt = electron.pfIsolationVariables().sumPhotonEt/electron.pt(); 
                elec_pfSumChargedHadronPt = electron.pfIsolationVariables().sumChargedHadronPt/electron.pt(); 
                elec_pfSumNeutralHadronEt = electron.pfIsolationVariables().sumNeutralHadronEt/electron.pt(); 
                elec_pfSumPUPt = electron.pfIsolationVariables().sumPUPt/electron.pt(); 


                elec_dr04EcalRecHitSumEt = electron.dr04EcalRecHitSumEt()/electron.pt(); 
                elec_dr04HcalDepth1TowerSumEt = electron.dr04HcalDepth1TowerSumEt()/electron.pt(); 
                elec_dr04HcalDepth1TowerSumEtBc = electron.dr04HcalDepth1TowerSumEtBc()/electron.pt(); 
                elec_dr04HcalDepth2TowerSumEt = electron.dr04HcalDepth2TowerSumEt()/electron.pt(); 
                elec_dr04HcalDepth2TowerSumEtBc = electron.dr04HcalDepth2TowerSumEtBc()/electron.pt();

                elec_dr04HcalTowerSumEt = electron.dr04HcalTowerSumEt()/electron.pt();
                elec_dr04HcalTowerSumEtBc = electron.dr04HcalTowerSumEtBc()/electron.pt();



}


