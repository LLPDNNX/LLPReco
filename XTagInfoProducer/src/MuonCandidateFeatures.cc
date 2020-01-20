 
#include "LLPReco/DataFormats/interface/MuonCandidateFeatures.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"


void llpdnnx::MuonCandidateFeatures::muonFeatures(const pat::Muon& muon  , const pat::Jet& jet , const reco::Vertex& pv ) {

                mu_isGlobal = muon.isGlobalMuon() ;                                   
                mu_isTight = muon.isTightMuon(pv);                                     
                mu_isMedium = muon.isMediumMuon();       
                mu_isLoose = muon.isLooseMuon() ; 
                mu_isStandAlone = muon.isStandAloneMuon() ; 

                mu_ptrel = muon.pt()/jet.pt() ; 
                mu_eta = muon.eta();                                                 
                mu_phi = muon.phi();                                                 
                mu_charge = muon.charge();        
                mu_energy = muon.energy()/muon.pt();                                           
                mu_et = muon.et();   
                mu_jetDeltaR = reco::deltaR(muon ,jet) ; 
                mu_numberOfMatchedStations = muon.numberOfMatchedStations();

                mu_2dIp = muon.dB() ; 
                mu_2dIpSig = muon.dB()/muon.edB() ; 
                mu_3dIp = muon.dB(pat::Muon::PV3D) ; 
                mu_3dIpSig = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D) ;


                reco::Candidate::Vector muonMom = muon.bestTrack()->momentum();

                mu_EtaRel =reco::btau::etaRel(jet.momentum().Unit(), muonMom);
                mu_dxy = muon.bestTrack()->dxy(pv.position());
                mu_dxyError = muon.bestTrack()->dxyError() ; 
                mu_dxySig = muon.bestTrack()->dxy(pv.position())/muon.bestTrack()->dxyError(); 
                mu_dz = muon.bestTrack()->dz(pv.position()) ; 
                mu_dzError = muon.bestTrack()->dzError() ;
                mu_numberOfValidPixelHits = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
                mu_numberOfpixelLayersWithMeasurement = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement() ;
                mu_numberOfstripLayersWithMeasurement = muon.bestTrack()->hitPattern().stripLayersWithMeasurement() ;


                mu_chi2 = muon.bestTrack()->chi2() ;  
                mu_ndof = muon.bestTrack()->ndof() ;


                mu_caloIso =  muon.caloIso()/muon.pt() ; 
                mu_ecalIso =  muon.ecalIso()/muon.pt() ; 
                mu_hcalIso =  muon.hcalIso()/muon.pt() ;     


                mu_sumPfChHadronPt  = muon.pfIsolationR04().sumChargedHadronPt/muon.pt();
                mu_sumPfNeuHadronEt  = muon.pfIsolationR04().sumNeutralHadronEt/muon.pt();
                mu_Pfpileup  = muon.pfIsolationR04().sumPUPt/muon.pt();
                mu_sumPfPhotonEt = muon.pfIsolationR04().sumPhotonEt/muon.pt();



                mu_sumPfChHadronPt03  = muon.pfIsolationR03().sumChargedHadronPt/muon.pt();
                mu_sumPfNeuHadronEt03  = muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
                mu_Pfpileup03  = muon.pfIsolationR03().sumPUPt/muon.pt();
                mu_sumPfPhotonEt03 = muon.pfIsolationR03().sumPhotonEt/muon.pt();       


                mu_timeAtIpInOut = muon.time().timeAtIpInOut ;  
                mu_timeAtIpInOutErr = muon.time().timeAtIpInOutErr ; 
                mu_timeAtIpOutIn = muon.time().timeAtIpOutIn ;  


}


