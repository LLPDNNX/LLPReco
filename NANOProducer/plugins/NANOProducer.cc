// -*- C++ -*-
//
// Package:    LLPReco/NANOProducer
// Class:      NANOProducer
// 
/**\class NANOProducer NANOProducer.cc LLPReco/NANOProducer/plugins/NANOProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vilius Cepaitis
//         Created:  Fri, 10 Jan 2020 15:20:19 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"
#include "LLPReco/DataFormats/interface/LLPLabel.h"
#include "LLPReco/DataFormats/interface/LLPLabelInfo.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"




//
// class declaration
//

class NANOProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOProducer(const edm::ParameterSet&);
      ~NANOProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _tag_src;
      const edm::EDGetTokenT<std::vector<reco::LLPLabelInfo>> _label_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOProducer::NANOProducer(const edm::ParameterSet& iConfig) :
    _tag_src(consumes<std::vector<reco::XTagInfo>>(iConfig.getParameter<edm::InputTag>("srcTags"))),
    _label_src(consumes<std::vector<reco::LLPLabelInfo>>(iConfig.getParameter<edm::InputTag>("srcLabels")))
{
    produces<nanoaod::FlatTable>("global");
    produces<nanoaod::FlatTable>("csv");
    produces<nanoaod::FlatTable>("cpf");
    produces<nanoaod::FlatTable>("npf");
    produces<nanoaod::FlatTable>("sv");
    produces<nanoaod::FlatTable>("length");
    produces<nanoaod::FlatTable>("jetorigin");
// New : 
    produces<nanoaod::FlatTable>("muon") ;
    produces<nanoaod::FlatTable>("electron") ; 
}


NANOProducer::~NANOProducer()
{
}


// ------------ method called to produce the data  ------------
void
NANOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
    
    edm::Handle<std::vector<reco::XTagInfo>> tag_infos;
    iEvent.getByToken(_tag_src, tag_infos);

    edm::Handle<std::vector<reco::LLPLabelInfo>> label_infos;
    iEvent.getByToken(_label_src, label_infos);

    unsigned int ntags = tag_infos->size();
    unsigned int ntruth = label_infos->size();
    if (ntags != ntruth) throw cms::Exception("NANOProducer:: number of jet tags is not equal to the number of labelled jets!");
    auto jetOriginTable = std::make_unique<nanoaod::FlatTable>(ntags, "jetorigin", false, false);

    std::vector<int> isPU;
    std::vector<int> isB;
    std::vector<int> isBB;
    std::vector<int> isGBB;
    std::vector<int> isLeptonic_B;
    std::vector<int> isLeptonic_C;
    std::vector<int> isC;
    std::vector<int> isCC;
    std::vector<int> isGCC;
    std::vector<int> isS;
    std::vector<int> isUD;
    std::vector<int> isG;

    std::vector<int> isLLP_RAD; //no flavour match (likely from wide angle radiation)
    std::vector<int> isLLP_MU; //prompt lepton
    std::vector<int> isLLP_E; //prompt lepton
    std::vector<int> isLLP_Q; //single light quark
    std::vector<int> isLLP_QMU; //single light quark + prompt lepton
    std::vector<int> isLLP_QE; //single light quark + prompt lepton
    std::vector<int> isLLP_QQ; //double light quark
    std::vector<int> isLLP_QQMU; //double light quark + prompt lepton
    std::vector<int> isLLP_QQE; //double light quark + prompt lepton
    std::vector<int> isLLP_B; //single b/c quark
    std::vector<int> isLLP_BMU; //single b/c quark + prompt lepton
    std::vector<int> isLLP_BE; //single b/c quark + prompt lepton
    std::vector<int> isLLP_BB; //double b/c quark
    std::vector<int> isLLP_BBMU; //double b/c quark + prompt lepton
    std::vector<int> isLLP_BBE; //double b/c quark + prompt lepton
    std::vector<int> isUndefined;
        
    std::vector<int> partonFlavor;
    std::vector<int> hadronFlavor;
    std::vector<int> llpId;

    std::vector<float> displacement;
    std::vector<float> displacement_xy;
    std::vector<float> displacement_z;
    std::vector<float> decay_angle;
    std::vector<float> betagamma;
        
    auto lengthTable = std::make_unique<nanoaod::FlatTable>(ntags, "length", false, false);
    std::vector<int> cpf_length;
    std::vector<int> npf_length;
    std::vector<int> sv_length;
    
    auto globalTable = std::make_unique<nanoaod::FlatTable>(ntags, "global", false, false);
    std::vector<float> pt;
    std::vector<float> eta;

    auto csvTable = std::make_unique<nanoaod::FlatTable>(ntags, "csv", false, false);
    std::vector<float> trackSumJetEtRatio;
    std::vector<float> trackSumJetDeltaR;
    std::vector<float> vertexCategory;
    std::vector<float> trackSip2dValAboveCharm;
    std::vector<float> trackSip2dSigAboveCharm;
    std::vector<float> trackSip3dValAboveCharm;
    std::vector<float> trackSip3dSigAboveCharm;
    std::vector<float> jetNSelectedTracks;
    std::vector<float> jetNTracksEtaRel;

    std::vector<float> cpf_trackEtaRel;
    std::vector<float> cpf_trackPtRel;
    std::vector<float> cpf_trackPPar;
    std::vector<float> cpf_trackDeltaR;
    std::vector<float> cpf_trackPParRatio;
    std::vector<float> cpf_trackPtRatio;
    std::vector<float> cpf_trackSip2dVal;
    std::vector<float> cpf_trackSip2dSig;
    std::vector<float> cpf_trackSip3dVal;
    std::vector<float> cpf_trackSip3dSig;
    std::vector<float> cpf_trackJetDistVal;
    std::vector<float> cpf_trackJetDistSig;
    std::vector<float> cpf_ptrel;
    std::vector<float> cpf_drminsv;
    std::vector<float> cpf_vertex_association;
    std::vector<float> cpf_fromPV;
    std::vector<float> cpf_puppi_weight;
    std::vector<float> cpf_track_chi2;
    std::vector<float> cpf_track_quality;
    std::vector<float> cpf_track_ndof;
    
    std::vector<float> npf_ptrel;
    std::vector<float> npf_deltaR;
    std::vector<float> npf_isGamma;
    std::vector<float> npf_hcal_fraction;
    std::vector<float> npf_drminsv;
    std::vector<float> npf_puppi_weight;
    
    std::vector<float> sv_pt;
    std::vector<float> sv_deltaR;
    std::vector<float> sv_mass;
    std::vector<float> sv_ntracks;
    std::vector<float> sv_chi2;
    std::vector<float> sv_ndof;
    std::vector<float> sv_dxy;
    std::vector<float> sv_dxysig;
    std::vector<float> sv_d3d;
    std::vector<float> sv_d3dsig;
    std::vector<float> sv_costhetasvpv;
    std::vector<float> sv_enratio;




    std::vector<int>  mu_isGlobal ; 
    std::vector<int>  mu_isTight ; 
    std::vector<int>  mu_isMedium ; 
    std::vector<int>  mu_isLoose ; 
    std::vector<int>  mu_isStandAlone ;

    std::vector<float> mu_pt ;
    std::vector<float> mu_p ;
    std::vector<float> mu_jetPtRel ;
    std::vector<float> mu_EtaRel; 
    std::vector<float> mu_jetPtRel2 ; 
    std::vector<float> mu_eta;
    std::vector<float> mu_phi;
    std::vector<float> mu_charge ; 
    std::vector<float> mu_energy;
    std::vector<float> mu_et ;
    std::vector<float> mu_jetDeltaR ; 
    std::vector<float> mu_numberOfMatchedStations ;
  
    std::vector<float> mu_2dIp ; 
    std::vector<float> mu_2dIpSig ;
    std::vector<float> mu_3dIp ; 
    std::vector<float> mu_3dIpSig ; 

    std::vector<float> mu_dxy ; 
    std::vector<float> mu_dxyError ; 
    std::vector<float> mu_dxySig ; 
    std::vector<float> mu_dz ; 
    std::vector<float> mu_dzError ; 
    std::vector<float> mu_numberOfValidPixelHits; 
    std::vector<float> mu_numberOfpixelLayersWithMeasurement ; 
    //std::vector<float> mu_numberOfstripLayersWithMeasurement ; 

    std::vector<float> mu_chi2 ; 
    std::vector<float> mu_ndof ; 
    
    std::vector<float> mu_caloIso ; 
    std::vector<float> mu_ecalIso ; 
    std::vector<float> mu_hcalIso ;
  
    std::vector<float> mu_sumPfChHadronPt ; 
    std::vector<float> mu_sumPfNeuHadronEt ; 
    std::vector<float> mu_Pfpileup ; 
    std::vector<float> mu_sumPfPhotonEt ; 
  
    std::vector<float> mu_sumPfChHadronPt03 ; 
    std::vector<float> mu_sumPfNeuHadronEt03 ; 
    std::vector<float> mu_Pfpileup03 ; 
    std::vector<float> mu_sumPfPhotonEt03 ; 
 
    std::vector<float> mu_sumChHadronPt ; 
    std::vector<float> mu_sumNeuHadronEt ; 
    std::vector<float> mu_pileup ; 
    std::vector<float> mu_sumPhotonEt ; 
  
    //std::vector<float> mu_sumChHadronPt03 ; 
    //std::vector<float> mu_sumNeuHadronEt03 ; 
    //std::vector<float> mu_pileup03 ; 
    //std::vector<float> mu_sumPhotonEt03 ;
    
   std::vector<float>  mu_timeAtIpInOut ; 
   std::vector<float>  mu_timeAtIpInOutErr ; 
   std::vector<float>  mu_timeAtIpOutIn ; 
 

// Electron Block 
  
    std::vector<float> elec_CanDr ;
    std::vector<float> elec_pt ;
    std::vector<float> elec_jetPtRatio ;
    std::vector<float> elec_jetDeltaR ; 
    std::vector<float> elec_p ;
    std::vector<float> elec_eta;
    std::vector<float> elec_phi;
    std::vector<float> elec_charge ; 
    std::vector<float> elec_energy;
    std::vector<float> elec_EtFromCaloEn ;
    std::vector<float> elec_isEB ; 
    std::vector<float> elec_isEE ; 
    std::vector<float> elec_ecalEnergy ; 
    std::vector<float> elec_isPassConversionVeto ;
    std::vector<float> elec_convDist ;
    std::vector<int>   elec_convFlags ;

    std::vector<float> elec_convRadius ; 
    std::vector<float> elec_hadronicOverEm ;
    std::vector<float> elec_ecalDrivenSeed;
	
     
   std::vector<float> elecSC_energy ; 
   std::vector<float> elecSC_eta ; 
   std::vector<float> elecSC_phi ;
   std::vector<float> elecSC_et ;
   std::vector<float> elecSC_eSuperClusterOverP ; 
   std::vector<float> elec_scE1x5Overe5x5 ; 
   std::vector<float> elec_scE2x5MaxOvere5x5 ; 
   std::vector<float> elec_scE5x5 ; 
   std::vector<float> elec_scE5x5Rel ; 
   std::vector<float> elec_scPixCharge ; 
   std::vector<float> elec_scSigmaEtaEta ;
   std::vector<float> elec_scSigmaIEtaIEta ;  
   std::vector<float> elec_superClusterFbrem ; 

   std::vector<float> elec_2dIP ; 
   std::vector<float> elec_2dIPSig ;
   std::vector<float> elec_3dIP ; 
   std::vector<float> elec_3dIPSig ; 
   std::vector<float> elec_eSeedClusterOverP ;
   std::vector<float> elec_eSeedClusterOverPout;
   std::vector<float> elec_eSuperClusterOverP;
   std::vector<float> elec_eTopOvere5x5; 
 
   std::vector<float> elec_deltaEtaEleClusterTrackAtCalo ; 
   std::vector<float> elec_deltaEtaSeedClusterTrackAtCalo ;
   std::vector<float> elec_deltaPhiSeedClusterTrackAtCalo ; 
   std::vector<float> elec_deltaEtaSeedClusterTrackAtVtx ; 
   std::vector<float> elec_deltaEtaSuperClusterTrackAtVtx ;
   std::vector<float> elec_deltaPhiEleClusterTrackAtCalo ; 
   std::vector<float> elec_deltaPhiSuperClusterTrackAtVtx ;
   std::vector<float> elec_sCseedEta ;  
///////
   std::vector<float> elec_EtaRel ; 
   std::vector<float> elec_dxy ; 
   std::vector<float> elec_dz ;
   std::vector<float> elec_nbOfMissingHits ; 
   std::vector<float> elec_gsfCharge ;

   std::vector<float> elec_e2x5MaxOvere5x5 ; 
   std::vector<float> elec_e1x5Overe5x5 ; 
   std::vector<float> elec_e5x5 ;
   std::vector<float> elec_e5x5Rel ;
   std::vector<float> elec_full5x5_sigmaIetaIeta ;
   std::vector<float> elec_full5x5_e1x5Overe5x5 ;
   std::vector<float> elec_full5x5_e2x5BottomOvere5x5 ;
   std::vector<float> elec_full5x5_e2x5LeftOvere5x5 ;
   std::vector<float> elec_full5x5_e2x5MaxOvere5x5 ;
   std::vector<float> elec_full5x5_e2x5RightOvere5x5 ;
   std::vector<float> elec_full5x5_e2x5TopOvere5x5 ;
   std::vector<float> elec_full5x5_e5x5 ;
   std::vector<float> elec_full5x5_e5x5Rel ;
   std::vector<float> elec_full5x5_eBottomOvere5x5 ;
   std::vector<float> elec_full5x5_eLeftOvere5x5;
   std::vector<float> elec_full5x5_eRightOvere5x5;
   std::vector<float> elec_full5x5_eTopOvere5x5;

   std::vector<float> elec_full5x5_hcalDepth1OverEcal ;
   std::vector<float> elec_full5x5_hcalDepth1OverEcalBc ;
   std::vector<float> elec_full5x5_hcalDepth2OverEcal;
   std::vector<float> elec_full5x5_hcalDepth2OverEcalBc ;
   std::vector<float> elec_full5x5_hcalOverEcal ;
   std::vector<float> elec_full5x5_hcalOverEcalBc;   
   std::vector<float> elec_full5x5_r9 ;
   std::vector<int>   elec_numberOfBrems ; 
   std::vector<float> elec_trackFbrem ; 
   std::vector<float> elec_fbrem ; 
   std::vector<float> elec_neutralHadronIso; 
   std::vector<float> elec_particleIso  ; 
   std::vector<float> elec_photonIso ;
   std::vector<float> elec_puChargedHadronIso ; 
   std::vector<float> elec_trackIso ;  
   std::vector<float> elec_hcalDepth1OverEcal ; 
   std::vector<float> elec_hcalDepth2OverEcal ; 
   std::vector<float> elec_ecalPFClusterIso ;
   std::vector<float> elec_hcalPFClusterIso ;  
   std::vector<float> elec_dr03TkSumPt ; 
   std::vector<float> elec_dr03HcalDepth2TowerSumEt ; 
   std::vector<float> elec_dr03EcalRecHitSumEt ; 
   std::vector<float> elec_dr03HcalDepth1TowerSumEt ;  
   std::vector<float> elec_dr03HcalDepth1TowerSumEtBc ; 
   std::vector<float> elec_pfSumPhotonEt ; 
   std::vector<float> elec_pfSumChargedHadronPt ; 
   std::vector<float> elec_pfSumNeutralHadronEt ; 
   std::vector<float> elec_pfSumPUPt ;
   std::vector<float> elec_dr04EcalRecHitSumEt ;  
   std::vector<float> elec_dr04HcalDepth1TowerSumEt ;  
   std::vector<float> elec_dr04HcalDepth1TowerSumEtBc ;
   std::vector<float> elec_dr04HcalDepth2TowerSumEt ; 
   std::vector<float> elec_dr04HcalDepth2TowerSumEtBc  ;
   std::vector<float> elec_dr04HcalTowerSumEt  ;
   std::vector<float> elec_dr04HcalTowerSumEtBc  ;

    unsigned int nmu_total = 0;
    unsigned int nelec_total = 0;
    unsigned int ncpf_total = 0;
    unsigned int nnpf_total = 0;
    unsigned int nsv_total = 0;

    for (unsigned int itag= 0; itag < ntags; itag++) {
         const auto& features = tag_infos->at(itag).features();
         const auto& tag_info_features = features.tag_info_features;

         unsigned int nmu = features.mu_features.size();
         unsigned int nelec = features.elec_features.size();
         unsigned int ncpf = features.cpf_features.size();
         unsigned int nnpf = features.npf_features.size();
         unsigned int nsv = features.sv_features.size();
	 nmu_total += nmu; 
	 nelec_total += nelec;
         ncpf_total += ncpf;
         nnpf_total += nnpf;
         nsv_total += nsv;
         cpf_length.push_back(ncpf);
         npf_length.push_back(nnpf);
         sv_length.push_back(nsv);

         pt.push_back(features.jet_features.pt);
         eta.push_back(features.jet_features.eta);

         trackSumJetEtRatio.push_back(tag_info_features.csv_trackSumJetEtRatio);
         trackSumJetDeltaR.push_back(tag_info_features.csv_trackSumJetDeltaR);
         trackSip2dValAboveCharm.push_back(tag_info_features.csv_trackSip2dValAboveCharm);
         trackSip2dSigAboveCharm.push_back(tag_info_features.csv_trackSip2dSigAboveCharm);
         trackSip3dValAboveCharm.push_back(tag_info_features.csv_trackSip3dValAboveCharm);
         trackSip3dSigAboveCharm.push_back(tag_info_features.csv_trackSip3dSigAboveCharm);
         jetNSelectedTracks.push_back(tag_info_features.csv_jetNSelectedTracks);
         jetNTracksEtaRel.push_back(tag_info_features.csv_jetNTracksEtaRel);
         vertexCategory.push_back(tag_info_features.csv_vertexCategory);


         const auto& labels = label_infos->at(itag).features();

         unsigned int _isPU{0};
         unsigned int _isB{0};
         unsigned int _isBB{0};
         unsigned int _isGBB{0};
         unsigned int _isLeptonic_B{0};
         unsigned int _isLeptonic_C{0};
         unsigned int _isC{0};
         unsigned int _isCC{0};
         unsigned int _isGCC{0};
         unsigned int _isS{0};
         unsigned int _isUD{0};
         unsigned int _isG{0};
         unsigned int _isLLP_RAD{0};
         unsigned int _isLLP_MU{0};
         unsigned int _isLLP_E{0};
         unsigned int _isLLP_Q{0};
         unsigned int _isLLP_QMU{0};
         unsigned int _isLLP_QE{0};
         unsigned int _isLLP_QQ{0};
         unsigned int _isLLP_QQMU{0};
         unsigned int _isLLP_QQE{0};
         unsigned int _isLLP_B{0};
         unsigned int _isLLP_BMU{0};
         unsigned int _isLLP_BE{0};
         unsigned int _isLLP_BB{0};
         unsigned int _isLLP_BBMU{0};
         unsigned int _isLLP_BBE{0};
         unsigned int _isUndefined{0};

         if (labels.type == llpdnnx::LLPLabel::Type::isPU) _isPU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isB) _isB = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isBB) _isBB = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isGBB) _isGBB = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLeptonic_B) _isLeptonic_B = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLeptonic_C) _isLeptonic_C = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isC) _isC = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isCC) _isCC = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isGCC) _isGCC = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isS) _isS = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isUD) _isUD = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isG) _isG = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_RAD) _isLLP_RAD = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_MU) _isLLP_MU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_E) _isLLP_E = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_Q) _isLLP_Q = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QMU) _isLLP_QMU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QE) _isLLP_QE = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQ) _isLLP_QQ = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQMU) _isLLP_QQMU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_QQE) _isLLP_QQE = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_B) _isLLP_B = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BMU) _isLLP_BMU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BE) _isLLP_BE = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BB) _isLLP_BB = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BBMU) _isLLP_BBMU = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isLLP_BBE) _isLLP_BBE = 1;
         if (labels.type == llpdnnx::LLPLabel::Type::isUndefined) _isUndefined = 1;




         isPU.push_back(_isPU);
         isB.push_back(_isB);
         isBB.push_back(_isBB);
         isGBB.push_back(_isGBB);
         isLeptonic_B.push_back(_isLeptonic_B);
         isLeptonic_C.push_back(_isLeptonic_C);
         isC.push_back(_isC);
         isCC.push_back(_isCC);
         isGCC.push_back(_isGCC);
         isS.push_back(_isS);
         isUD.push_back(_isUD);
         isG.push_back(_isG);
         isLLP_RAD.push_back(_isLLP_RAD);
         isLLP_MU.push_back(_isLLP_MU);
         isLLP_E.push_back(_isLLP_E);
         isLLP_Q.push_back(_isLLP_Q);
         isLLP_QMU.push_back(_isLLP_QMU);
         isLLP_QE.push_back(_isLLP_QE);
         isLLP_QQ.push_back(_isLLP_QQ);
         isLLP_QQMU.push_back(_isLLP_QQMU);
         isLLP_QQE.push_back(_isLLP_QQE);
         isLLP_B.push_back(_isLLP_B);
         isLLP_BMU.push_back(_isLLP_BMU);
         isLLP_BE.push_back(_isLLP_BE);
         isLLP_BB.push_back(_isLLP_BB);
         isLLP_BBMU.push_back(_isLLP_BBMU);
         isLLP_BBE.push_back(_isLLP_BBE);
         isUndefined.push_back(_isUndefined);

         partonFlavor.push_back(labels.partonFlavor);
         hadronFlavor.push_back(labels.hadronFlavor);
         llpId.push_back(labels.llpId);
         displacement.push_back(labels.displacement);
         displacement_xy.push_back(labels.displacement_xy);
         displacement_z.push_back(labels.displacement_z);
         decay_angle.push_back(labels.decay_angle);
         betagamma.push_back(labels.betagamma);
 
    }

    auto muonTable = std::make_unique<nanoaod::FlatTable>(nmu_total, "muon", false, false);
    auto electronTable = std::make_unique<nanoaod::FlatTable>(nelec_total, "electron", false, false);
    auto cpfTable = std::make_unique<nanoaod::FlatTable>(ncpf_total, "cpf", false, false);
    auto npfTable = std::make_unique<nanoaod::FlatTable>(nnpf_total, "npf", false, false);
    auto svTable = std::make_unique<nanoaod::FlatTable>(nsv_total, "sv", false, false);

    for (unsigned int itag= 0; itag < ntags; itag++) {
         const auto& features = tag_infos->at(itag).features();
	 auto mu  = features.mu_features ; 
	 auto elec = features.elec_features ; 
         auto cpf = features.cpf_features;
         auto npf = features.npf_features;
         auto sv = features.sv_features ;

	 unsigned int nmu = features.mu_features.size();
         unsigned int nelec = features.elec_features.size();

         unsigned int ncpf = features.cpf_features.size();
         unsigned int nnpf = features.npf_features.size();
         unsigned int nsv = features.sv_features.size();

         for (unsigned int i = 0; i < ncpf; i++){
                const auto& cpf_features = cpf.at(i);
                cpf_trackEtaRel.push_back(cpf_features.cpf_trackEtaRel);
                cpf_trackPtRel.push_back(cpf_features.cpf_trackPtRel);
                cpf_trackPPar.push_back(cpf_features.cpf_trackPPar);
                cpf_trackDeltaR.push_back(cpf_features.cpf_trackDeltaR);
                cpf_trackPParRatio.push_back(cpf_features.cpf_trackPParRatio);
                cpf_trackPtRatio.push_back(cpf_features.cpf_trackPtRatio);
                cpf_trackSip2dVal.push_back(cpf_features.cpf_trackSip2dVal);
                cpf_trackSip2dSig.push_back(cpf_features.cpf_trackSip2dSig);
                cpf_trackSip3dVal.push_back(cpf_features.cpf_trackSip3dVal);
                cpf_trackSip3dSig.push_back(cpf_features.cpf_trackSip3dSig);
                cpf_trackJetDistVal.push_back(cpf_features.cpf_trackJetDistVal);
                cpf_trackJetDistSig.push_back(cpf_features.cpf_trackJetDistSig);
                cpf_ptrel.push_back(cpf_features.cpf_ptrel);
                cpf_drminsv.push_back(cpf_features.cpf_drminsv);
                cpf_vertex_association.push_back(cpf_features.cpf_vertex_association);
                cpf_fromPV.push_back(cpf_features.cpf_fromPV);
                cpf_puppi_weight.push_back(cpf_features.cpf_puppi_weight);
                cpf_track_chi2.push_back(cpf_features.cpf_track_chi2);
                cpf_track_quality.push_back(cpf_features.cpf_track_quality);
                cpf_track_ndof.push_back(cpf_features.cpf_track_ndof);
         }
        
        for (unsigned int i = 0; i < nnpf; i++){
            const auto& npf_features = npf.at(i);
            npf_ptrel.push_back(npf_features.npf_ptrel);
            npf_deltaR.push_back(npf_features.npf_deltaR);
            npf_isGamma.push_back(npf_features.npf_isGamma);
            npf_hcal_fraction.push_back(npf_features.npf_hcal_fraction);
            npf_drminsv.push_back(npf_features.npf_drminsv);
            npf_puppi_weight.push_back(npf_features.npf_puppi_weight);
        }
        
        for (unsigned int i = 0; i < nsv; i++){
            const auto& sv_features = sv.at(i);
            sv_pt.push_back(sv_features.sv_pt);
            sv_deltaR.push_back(sv_features.sv_deltaR);
            sv_mass.push_back(sv_features.sv_mass);
            sv_ntracks.push_back(sv_features.sv_ntracks);
            sv_chi2.push_back(sv_features.sv_chi2);
            sv_ndof.push_back(sv_features.sv_ndof);
            sv_dxy.push_back(sv_features.sv_dxy);
            sv_dxysig.push_back(sv_features.sv_dxysig);
            sv_d3d.push_back(sv_features.sv_d3d);
            sv_d3dsig.push_back(sv_features.sv_d3dsig);
            sv_costhetasvpv.push_back(sv_features.sv_costhetasvpv);
            sv_enratio.push_back(sv_features.sv_enratio);
        }



       for(unsigned int i = 0; i < nmu; i++){

	const auto& mu_features = mu.at(i);
	
        mu_isGlobal.push_back(mu_features.mu_isGlobal) ; 
        mu_isTight.push_back(mu_features.mu_isTight) ; 
        mu_isMedium.push_back(mu_features.mu_isMedium) ; 
        mu_isLoose.push_back(mu_features.mu_isLoose) ; 
        mu_isStandAlone.push_back(mu_features.mu_isStandAlone) ;


        mu_pt.push_back(mu_features.mu_pt);
        mu_p.push_back(mu_features.mu_p);
        mu_jetPtRel.push_back(mu_features.mu_jetPtRel);
        mu_EtaRel.push_back(mu_features.mu_EtaRel); 
        mu_jetPtRel2.push_back(mu_features.mu_jetPtRel2); 
        mu_eta.push_back(mu_features.mu_eta);
        mu_phi.push_back(mu_features.mu_phi);
        mu_charge.push_back(mu_features.mu_charge); 
        mu_energy.push_back(mu_features.mu_energy);
        mu_et.push_back(mu_features.mu_et);
        mu_jetDeltaR.push_back(mu_features.mu_jetDeltaR); 
        mu_numberOfMatchedStations.push_back(mu_features.mu_numberOfMatchedStations);
  
        mu_2dIp.push_back(mu_features.mu_2dIp); 
        mu_2dIpSig.push_back(mu_features.mu_2dIpSig);
        mu_3dIp.push_back(mu_features.mu_3dIp); 
        mu_3dIpSig.push_back(mu_features.mu_3dIpSig); 

        mu_dxy.push_back(mu_features.mu_dxy); 
        mu_dxyError.push_back(mu_features.mu_dxyError); 
        mu_dxySig.push_back(mu_features.mu_dxySig); 
        mu_dz.push_back(mu_features.mu_dz); 
        mu_dzError.push_back(mu_features.mu_dzError); 
        mu_numberOfValidPixelHits.push_back(mu_features.mu_numberOfValidPixelHits); 
        mu_numberOfpixelLayersWithMeasurement.push_back(mu_features.mu_numberOfpixelLayersWithMeasurement);

     
        mu_chi2.push_back(mu_features.mu_chi2); 
        mu_ndof.push_back(mu_features.mu_ndof); 
    
        mu_caloIso.push_back(mu_features.mu_caloIso); 
        mu_ecalIso.push_back(mu_features.mu_ecalIso); 
        mu_hcalIso.push_back(mu_features.mu_hcalIso);
  
        mu_sumPfChHadronPt.push_back(mu_features.mu_sumPfChHadronPt); 
        mu_sumPfNeuHadronEt.push_back(mu_features.mu_sumPfNeuHadronEt); 
        mu_Pfpileup.push_back(mu_features.mu_Pfpileup); 
        mu_sumPfPhotonEt.push_back(mu_features.mu_sumPfPhotonEt); 
  
        mu_sumPfChHadronPt03.push_back(mu_features.mu_sumPfChHadronPt03); 
        mu_sumPfNeuHadronEt03.push_back(mu_features.mu_sumPfNeuHadronEt03); 
        mu_Pfpileup03.push_back(mu_features.mu_Pfpileup03); 
        mu_sumPfPhotonEt03.push_back(mu_features.mu_sumPfPhotonEt03); 
 
        mu_sumChHadronPt.push_back(mu_features.mu_sumChHadronPt); 
        mu_sumNeuHadronEt.push_back(mu_features.mu_sumNeuHadronEt); 
        mu_pileup.push_back(mu_features.mu_pileup); 
        mu_sumPhotonEt.push_back(mu_features.mu_sumPhotonEt);
	 
	mu_timeAtIpInOut.push_back(mu_features.mu_timeAtIpInOut) ; 
        mu_timeAtIpInOutErr.push_back(mu_features.mu_timeAtIpInOutErr) ;  
        mu_timeAtIpOutIn.push_back(mu_features.mu_timeAtIpOutIn) ; 
 

       }

       for(unsigned int i = 0; i < nelec ; i++){

         const auto& elec_features = elec.at(i);


      elec_CanDr.push_back(elec_features.elec_CanDr );
      elec_pt.push_back(elec_features.elec_pt );
      elec_jetPtRatio.push_back(elec_features.elec_jetPtRatio );
      elec_jetDeltaR.push_back(elec_features.elec_jetDeltaR ); 
      elec_p.push_back(elec_features.elec_p );
      elec_eta.push_back(elec_features.elec_eta );
      elec_phi.push_back(elec_features.elec_phi );
      elec_charge.push_back(elec_features.elec_charge );
 
      elec_energy.push_back(elec_features.elec_energy );
      elec_EtFromCaloEn.push_back(elec_features.elec_EtFromCaloEn );
      elec_isEB.push_back(elec_features.elec_isEB ); 
      elec_isEE.push_back(elec_features.elec_isEE ); 
      elec_ecalEnergy.push_back(elec_features.elec_ecalEnergy ); 
      elec_isPassConversionVeto.push_back(elec_features.elec_isPassConversionVeto );
      elec_convDist.push_back(elec_features.elec_convDist );
      elec_convFlags.push_back(elec_features.elec_convFlags );

      elec_convRadius.push_back(elec_features.elec_convRadius ); 
      elec_hadronicOverEm.push_back(elec_features.elec_hadronicOverEm );
      elec_ecalDrivenSeed.push_back(elec_features.elec_ecalDrivenSeed );
	
     
     elecSC_energy.push_back(elec_features.elecSC_energy ); 
     elecSC_eta.push_back(elec_features.elecSC_eta ); 
     elecSC_phi.push_back(elec_features.elecSC_phi );
     elecSC_et.push_back(elec_features.elecSC_et );
     elecSC_eSuperClusterOverP.push_back(elec_features.elecSC_eSuperClusterOverP );
 
     elec_scE1x5Overe5x5.push_back(elec_features.elec_scE1x5Overe5x5 ); 
     elec_scE2x5MaxOvere5x5.push_back(elec_features.elec_scE2x5MaxOvere5x5 ); 
     elec_scE5x5.push_back(elec_features.elec_scE5x5 ); 
     elec_scE5x5Rel.push_back(elec_features.elec_scE5x5Rel ); 

     elec_scPixCharge.push_back(elec_features.elec_scPixCharge ); 
     elec_scSigmaEtaEta.push_back(elec_features.elec_scSigmaEtaEta );
     elec_scSigmaIEtaIEta.push_back(elec_features.elec_scSigmaIEtaIEta );  
     elec_superClusterFbrem.push_back(elec_features.elec_superClusterFbrem ); 

     elec_2dIP.push_back(elec_features.elec_2dIP ); 
     elec_2dIPSig.push_back(elec_features.elec_2dIPSig );
     elec_3dIP.push_back(elec_features.elec_3dIP ); 
     elec_3dIPSig.push_back(elec_features.elec_3dIPSig ); 
     elec_eSeedClusterOverP.push_back(elec_features.elec_eSeedClusterOverP );
     elec_eSeedClusterOverPout.push_back(elec_features.elec_eSeedClusterOverPout );
     elec_eSuperClusterOverP.push_back(elec_features.elec_eSuperClusterOverP );
     elec_eTopOvere5x5.push_back(elec_features.elec_eTopOvere5x5 ); 
 
     elec_deltaEtaEleClusterTrackAtCalo.push_back(elec_features.elec_deltaEtaEleClusterTrackAtCalo ); 
     elec_deltaEtaSeedClusterTrackAtCalo.push_back(elec_features.elec_deltaEtaSeedClusterTrackAtCalo );
     elec_deltaPhiSeedClusterTrackAtCalo.push_back(elec_features.elec_deltaPhiSeedClusterTrackAtCalo ); 
     elec_deltaEtaSeedClusterTrackAtVtx.push_back(elec_features.elec_deltaEtaSeedClusterTrackAtVtx ); 
     elec_deltaEtaSuperClusterTrackAtVtx.push_back(elec_features.elec_deltaEtaSuperClusterTrackAtVtx );
     elec_deltaPhiEleClusterTrackAtCalo.push_back(elec_features.elec_deltaPhiEleClusterTrackAtCalo ); 
     elec_deltaPhiSuperClusterTrackAtVtx.push_back(elec_features.elec_deltaPhiSuperClusterTrackAtVtx );
     elec_sCseedEta.push_back(elec_features.elec_sCseedEta );  
	

    elec_EtaRel.push_back(elec_features.elec_EtaRel ); 
    elec_dxy.push_back(elec_features.elec_dxy ); 
    elec_dz.push_back(elec_features.elec_dz );
    elec_nbOfMissingHits.push_back(elec_features.elec_nbOfMissingHits ); 
    elec_gsfCharge.push_back(elec_features.elec_gsfCharge );

    elec_full5x5_sigmaIetaIeta.push_back(elec_features.elec_full5x5_sigmaIetaIeta );
    elec_full5x5_e1x5Overe5x5.push_back(elec_features.elec_full5x5_e1x5Overe5x5 );
    elec_full5x5_e2x5BottomOvere5x5.push_back(elec_features.elec_full5x5_e2x5BottomOvere5x5 );
    elec_full5x5_e2x5LeftOvere5x5.push_back(elec_features.elec_full5x5_e2x5LeftOvere5x5 );
    elec_full5x5_e2x5MaxOvere5x5.push_back(elec_features.elec_full5x5_e2x5MaxOvere5x5 );
    elec_full5x5_e2x5RightOvere5x5.push_back(elec_features.elec_full5x5_e2x5RightOvere5x5 );
    elec_full5x5_e2x5TopOvere5x5.push_back(elec_features.elec_full5x5_e2x5TopOvere5x5 );
    elec_full5x5_e5x5.push_back(elec_features.elec_full5x5_e5x5 );
    elec_full5x5_e5x5Rel.push_back(elec_features.elec_full5x5_e5x5Rel );

    elec_full5x5_eBottomOvere5x5.push_back(elec_features.elec_full5x5_eBottomOvere5x5 );
    elec_full5x5_eLeftOvere5x5.push_back(elec_features.elec_full5x5_eLeftOvere5x5 );
    elec_full5x5_eRightOvere5x5.push_back(elec_features.elec_full5x5_eRightOvere5x5 );
    elec_full5x5_eTopOvere5x5.push_back(elec_features.elec_full5x5_eTopOvere5x5 );
    elec_full5x5_hcalDepth1OverEcal.push_back(elec_features.elec_full5x5_hcalDepth1OverEcal );
    elec_full5x5_hcalDepth1OverEcalBc.push_back(elec_features.elec_full5x5_hcalDepth1OverEcalBc );
    elec_full5x5_hcalDepth2OverEcal.push_back(elec_features.elec_full5x5_hcalDepth2OverEcal );
    elec_full5x5_hcalDepth2OverEcalBc.push_back(elec_features.elec_full5x5_hcalDepth2OverEcalBc );
    elec_full5x5_hcalOverEcal.push_back(elec_features.elec_full5x5_hcalOverEcal );
    elec_full5x5_hcalOverEcalBc.push_back(elec_features.elec_full5x5_hcalOverEcalBc );   
    elec_full5x5_r9.push_back(elec_features.elec_full5x5_r9 );
    elec_numberOfBrems.push_back(elec_features.elec_numberOfBrems ); 
    elec_trackFbrem.push_back(elec_features.elec_trackFbrem ); 
    elec_fbrem.push_back(elec_features.elec_fbrem ); 
    elec_e2x5MaxOvere5x5.push_back(elec_features.elec_e2x5MaxOvere5x5 ); 
    elec_e1x5Overe5x5.push_back(elec_features.elec_e1x5Overe5x5 ); 
    elec_e5x5.push_back(elec_features.elec_e5x5 );
    elec_e5x5Rel.push_back(elec_features.elec_e5x5Rel );
    elec_neutralHadronIso.push_back(elec_features.elec_neutralHadronIso ); 
    elec_particleIso .push_back(elec_features.elec_particleIso ); 
    elec_photonIso.push_back(elec_features.elec_photonIso );
    elec_puChargedHadronIso.push_back(elec_features.elec_puChargedHadronIso ); 
    elec_trackIso.push_back(elec_features.elec_trackIso );  
    elec_hcalDepth1OverEcal.push_back(elec_features.elec_hcalDepth1OverEcal ); 
    elec_hcalDepth2OverEcal.push_back(elec_features.elec_hcalDepth2OverEcal ); 
    elec_ecalPFClusterIso.push_back(elec_features.elec_ecalPFClusterIso );
    elec_hcalPFClusterIso.push_back(elec_features.elec_hcalPFClusterIso );  
    elec_dr03TkSumPt.push_back(elec_features.elec_dr03TkSumPt ); 
    elec_dr03HcalDepth2TowerSumEt.push_back(elec_features.elec_dr03HcalDepth2TowerSumEt ); 
    elec_dr03EcalRecHitSumEt.push_back(elec_features.elec_dr03EcalRecHitSumEt ); 
    elec_dr03HcalDepth1TowerSumEt.push_back(elec_features.elec_dr03HcalDepth1TowerSumEt );  
    elec_dr03HcalDepth1TowerSumEtBc.push_back(elec_features.elec_dr03HcalDepth1TowerSumEtBc ); 
    elec_pfSumPhotonEt.push_back(elec_features.elec_pfSumPhotonEt ); 
    elec_pfSumChargedHadronPt.push_back(elec_features.elec_pfSumChargedHadronPt ); 
    elec_pfSumNeutralHadronEt.push_back(elec_features.elec_pfSumNeutralHadronEt ); 
    elec_pfSumPUPt.push_back(elec_features.elec_pfSumPUPt );
    elec_dr04EcalRecHitSumEt.push_back(elec_features.elec_dr04EcalRecHitSumEt );  
    elec_dr04HcalDepth1TowerSumEt.push_back(elec_features.elec_dr04HcalDepth1TowerSumEt );  
    elec_dr04HcalDepth1TowerSumEtBc.push_back(elec_features.elec_dr04HcalDepth1TowerSumEtBc );
    elec_dr04HcalDepth2TowerSumEt.push_back(elec_features.elec_dr04HcalDepth2TowerSumEt ); 
    elec_dr04HcalDepth2TowerSumEtBc.push_back(elec_features.elec_dr04HcalDepth2TowerSumEtBc );
    elec_dr04HcalTowerSumEt.push_back(elec_features.elec_dr04HcalTowerSumEt );
    elec_dr04HcalTowerSumEtBc.push_back(elec_features.elec_dr04HcalTowerSumEtBc );





       }
    }

    globalTable->addColumn<float>("pt", pt, "global jet pt (log 10, uncorrected)", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("eta", eta, "global jet eta", nanoaod::FlatTable::FloatColumn);

    csvTable->addColumn<float>("trackSumJetEtRatio", trackSumJetEtRatio, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSumJetDeltaR", trackSumJetDeltaR, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("vertexCategory", vertexCategory, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dValAboveCharm", trackSip2dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dSigAboveCharm", trackSip2dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dValAboveCharm", trackSip3dValAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dSigAboveCharm", trackSip3dSigAboveCharm, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("jetNSelectedTracks", jetNSelectedTracks, "doc", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("jetNTracksEtaRel", jetNTracksEtaRel, "doc", nanoaod::FlatTable::FloatColumn);
    
    cpfTable->addColumn<float>("trackEtaRel", cpf_trackEtaRel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRel", cpf_trackPtRel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPPar", cpf_trackPPar, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackDeltaR", cpf_trackDeltaR, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPParRatio", cpf_trackPParRatio, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRatio", cpf_trackPtRatio, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dVal", cpf_trackSip2dVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dSig", cpf_trackSip2dSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dVal", cpf_trackSip3dVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dSig", cpf_trackSip3dSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistVal", cpf_trackJetDistVal, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistSig", cpf_trackJetDistSig, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("ptrel", cpf_ptrel, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("drminsv", cpf_drminsv, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("vertex_association", cpf_vertex_association, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("fromPV", cpf_fromPV, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("puppi_weight", cpf_puppi_weight, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_chi2", cpf_track_chi2, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_quality", cpf_track_quality, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_ndof", cpf_track_ndof, "doc", nanoaod::FlatTable::FloatColumn);
    
    npfTable->addColumn<float>("ptrel", npf_ptrel, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("deltaR", npf_deltaR, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("isGamma", npf_isGamma, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("hcal_fraction", npf_hcal_fraction, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("drminsv", npf_drminsv, "doc", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("puppi_weight", npf_puppi_weight, "doc", nanoaod::FlatTable::FloatColumn);
    
    svTable->addColumn<float>("sv_pt", sv_pt, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_deltaR", sv_deltaR, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_mass", sv_mass, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_ntracks", sv_ntracks, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_chi2", sv_chi2,  "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_ndof", sv_ndof, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_dxy", sv_dxy, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_dxysig", sv_dxysig, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_d3d", sv_d3d, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_d3dsig", sv_d3dsig, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_costhetasvpv", sv_costhetasvpv, "doc", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("sv_enratio", sv_enratio, "doc", nanoaod::FlatTable::FloatColumn);
    
    lengthTable->addColumn<int>("cpf", cpf_length, "cpf offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("npf", npf_length, "npf offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("sv", sv_length, "sv offset", nanoaod::FlatTable::IntColumn);

    jetOriginTable->addColumn<int>("isPU", isPU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isB", isB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isBB", isBB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isGBB", isGBB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_B", isLeptonic_B, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLeptonic_C", isLeptonic_C, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isC", isC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isCC", isCC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isGCC", isGCC, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isS", isS, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isUD", isUD, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isG", isG, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_RAD", isLLP_RAD, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_MU", isLLP_MU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_E", isLLP_E, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_Q", isLLP_Q, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QMU", isLLP_QMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QE", isLLP_QE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQ", isLLP_QQ, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQMU", isLLP_QQMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_QQE", isLLP_QQE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_B", isLLP_B, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BMU", isLLP_BMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BE", isLLP_BE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BB", isLLP_BB, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBMU", isLLP_BBMU, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isLLP_BBE", isLLP_BBE, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("isUndefined", isUndefined, "doc", nanoaod::FlatTable::IntColumn);

    jetOriginTable->addColumn<int>("partonFlavor", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("hadronFlavor", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    jetOriginTable->addColumn<int>("llpId", partonFlavor, "doc", nanoaod::FlatTable::IntColumn);
    
    jetOriginTable->addColumn<float>("displacement", displacement, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_xy", displacement_xy, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("displacement_z", displacement_z, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("decay_angle", decay_angle, "doc", nanoaod::FlatTable::FloatColumn);
    jetOriginTable->addColumn<float>("betagamma", betagamma, "doc", nanoaod::FlatTable::FloatColumn);
     

     muonTable->addColumn<int>("isGlobal",mu_isGlobal,"doc",nanoaod::FlatTable::IntColumn); 
     muonTable->addColumn<int>("isTight",mu_isTight,"doc",nanoaod::FlatTable::IntColumn); 
     muonTable->addColumn<int>("isMedium",mu_isMedium,"doc",nanoaod::FlatTable::IntColumn); 
     muonTable->addColumn<int>("isLoose",mu_isLoose,"doc",nanoaod::FlatTable::IntColumn); 
     muonTable->addColumn<int>("isStandAlone",mu_isStandAlone,"doc", nanoaod::FlatTable::IntColumn);

     muonTable->addColumn<float>("pt",mu_pt,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("p",mu_p,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("jetPtRel",mu_jetPtRel,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("EtaRel",mu_EtaRel,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("jetPtRel2",mu_jetPtRel2,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("phi",mu_phi,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("charge",mu_charge,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("energy",mu_energy,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("jetDeltaR",mu_jetDeltaR,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("numberOfMatchedStations",mu_numberOfMatchedStations,"doc",nanoaod::FlatTable::FloatColumn);
  
     muonTable->addColumn<float>("2dIp",mu_2dIp,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("2dIpSig",mu_2dIpSig,"doc",nanoaod::FlatTable::FloatColumn);
     muonTable->addColumn<float>("3dIp",mu_3dIp,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("3dIpSig",mu_3dIpSig,"doc",nanoaod::FlatTable::FloatColumn); 

     muonTable->addColumn<float>("dxy",mu_dxy,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("dxyError",mu_dxyError,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("dxySig",mu_dxySig,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("dz",mu_dz,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("dzError",mu_dzError,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("numberOfValidPixelHits",mu_numberOfValidPixelHits,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("numberOfpixelLayersWithMeasurement",mu_numberOfpixelLayersWithMeasurement,"doc",nanoaod::FlatTable::FloatColumn); 
     //muonTable->addColumn<float>("numberOfstripLayersWithMeasurement",mu_numberOfstripLayersWithMeasurement,"doc",nanoaod::FlatTable::FloatColumn); 

     muonTable->addColumn<float>("chi2",mu_chi2,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("ndof",mu_ndof,"doc",nanoaod::FlatTable::FloatColumn); 
    
     muonTable->addColumn<float>("caloIso",mu_caloIso,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("ecalIso",mu_ecalIso,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("hcalIso",mu_hcalIso,"doc",nanoaod::FlatTable::FloatColumn);
  
     muonTable->addColumn<float>("sumPfChHadronPt",mu_sumPfChHadronPt,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumPfNeuHadronEt",mu_sumPfNeuHadronEt,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("Pfpileup",mu_Pfpileup,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumPfPhotonEt",mu_sumPfPhotonEt,"doc",nanoaod::FlatTable::FloatColumn); 
  
     muonTable->addColumn<float>("sumPfChHadronPt03",mu_sumPfChHadronPt03,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumPfNeuHadronEt03",mu_sumPfNeuHadronEt03,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("Pfpileup03",mu_Pfpileup03,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumPfPhotonEt03",mu_sumPfPhotonEt03,"doc",nanoaod::FlatTable::FloatColumn); 
 
     muonTable->addColumn<float>("sumChHadronPt",mu_sumChHadronPt,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumNeuHadronEt",mu_sumNeuHadronEt,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("pileup",mu_pileup,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("sumPhotonEt",mu_sumPhotonEt,"doc",nanoaod::FlatTable::FloatColumn); 
  
     //muonTable->addColumn<float>("sumChHadronPt03",mu_sumChHadronPt03,"doc",nanoaod::FlatTable::FloatColumn); 
     //muonTable->addColumn<float>("sumNeuHadronEt03",mu_sumPhotonEt03,"doc",nanoaod::FlatTable::FloatColumn); 
     //muonTable->addColumn<float>("pileup03",mu_sumPhotonEt03,"doc",nanoaod::FlatTable::FloatColumn); 
     //muonTable->addColumn<float>("sumPhotonEt03",mu_sumPhotonEt03,"doc",nanoaod::FlatTable::FloatColumn);
     //

     muonTable->addColumn<float>("timeAtIpInOut",mu_timeAtIpInOut,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("timeAtIpInOutErr",mu_timeAtIpInOutErr,"doc",nanoaod::FlatTable::FloatColumn); 
     muonTable->addColumn<float>("timeAtIpOutIn",mu_timeAtIpOutIn,"doc",nanoaod::FlatTable::FloatColumn); 
 

// Electron block : 


       electronTable->addColumn<float>("CanDr",elec_CanDr,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("pt",elec_pt,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("jetPtRatio",elec_jetPtRatio,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("jetDeltaR",elec_jetDeltaR,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("p",elec_p,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("eta",elec_eta,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("phi",elec_phi,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("charge",elec_charge,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("energy",elec_energy,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("EtFromCaloEn",elec_EtFromCaloEn,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("isEB",elec_isEB,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("isEE",elec_isEE,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("ecalEnergy",elec_ecalEnergy,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("isPassConversionVeto",elec_isPassConversionVeto,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("convDist",elec_convDist,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<int>("convFlags",elec_convFlags,"doc",nanoaod::FlatTable::IntColumn);

       electronTable->addColumn<float>("convRadius",elec_convRadius,"doc",nanoaod::FlatTable::FloatColumn); 
       electronTable->addColumn<float>("hadronicOverEm",elec_hadronicOverEm,"doc",nanoaod::FlatTable::FloatColumn);
       electronTable->addColumn<float>("ecalDrivenSeed",elec_ecalDrivenSeed,"doc",nanoaod::FlatTable::FloatColumn);
	
     
      electronTable->addColumn<float>("SC_energy",elecSC_energy,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("SC_eta",elecSC_eta,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("SC_phi",elecSC_phi,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("SC_et",elecSC_et,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("SC_eSuperClusterOverP",elecSC_eSuperClusterOverP,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scE1x5Overe5x5",elec_scE1x5Overe5x5,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scE2x5MaxOvere5x5",elec_scE2x5MaxOvere5x5,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scE5x5",elec_scE5x5,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scE5x5Rel",elec_scE5x5Rel,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scPixCharge",elec_scPixCharge,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("scSigmaEtaEta",elec_scSigmaEtaEta,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("scSigmaIEtaIEta",elec_scSigmaIEtaIEta,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("superClusterFbrem",elec_superClusterFbrem,"doc",nanoaod::FlatTable::FloatColumn); 

      electronTable->addColumn<float>("2dIP",elec_2dIP,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("2dIPSig",elec_2dIPSig,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("3dIP",elec_3dIP,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("3dIPSig",elec_3dIPSig,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("eSeedClusterOverP",elec_eSeedClusterOverP,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("eSeedClusterOverPout",elec_eSeedClusterOverPout,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("eSuperClusterOverP",elec_eSuperClusterOverP,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("eTopOvere5x5",elec_eTopOvere5x5,"doc",nanoaod::FlatTable::FloatColumn); 
 
      electronTable->addColumn<float>("deltaEtaEleClusterTrackAtCalo",elec_deltaEtaEleClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("deltaEtaSeedClusterTrackAtCalo",elec_deltaEtaSeedClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("deltaPhiSeedClusterTrackAtCalo",elec_deltaPhiSeedClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("deltaEtaSeedClusterTrackAtVtx",elec_deltaEtaSeedClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("deltaEtaSuperClusterTrackAtVtx",elec_deltaEtaSuperClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("deltaPhiEleClusterTrackAtCalo",elec_deltaPhiEleClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("deltaPhiSuperClusterTrackAtVtx",elec_deltaPhiSuperClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("sCseedEta",elec_sCseedEta,"doc",nanoaod::FlatTable::FloatColumn); 
 
      electronTable->addColumn<float>("EtaRel",elec_EtaRel,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dxy",elec_dxy,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dz",elec_dz,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("nbOfMissingHits",elec_nbOfMissingHits,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("gsfCharge",elec_gsfCharge,"doc",nanoaod::FlatTable::FloatColumn);

      electronTable->addColumn<int>("numberOfBrems",elec_numberOfBrems,"doc",nanoaod::FlatTable::IntColumn); 
      electronTable->addColumn<float>("trackFbrem",elec_trackFbrem,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("fbrem",elec_fbrem,"doc",nanoaod::FlatTable::FloatColumn); 

      electronTable->addColumn<float>("e5x5",elec_e5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("e5x5Rel",elec_e5x5Rel,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("e2x5MaxOvere5x5",elec_e2x5MaxOvere5x5,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("e1x5Overe5x5",elec_e1x5Overe5x5,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("full5x5_sigmaIetaIeta",elec_full5x5_sigmaIetaIeta,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e1x5Overe5x5",elec_full5x5_e1x5Overe5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e2x5BottomOvere5x5",elec_full5x5_e2x5BottomOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e2x5LeftOvere5x5",elec_full5x5_e2x5LeftOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e2x5MaxOvere5x5",elec_full5x5_e2x5MaxOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e2x5RightOvere5x5",elec_full5x5_e2x5RightOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e2x5TopOvere5x5",elec_full5x5_e2x5TopOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e5x5",elec_full5x5_e5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_e5x5Rel",elec_full5x5_e5x5Rel,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_eBottomOvere5x5",elec_full5x5_eBottomOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_eLeftOvere5x5",elec_full5x5_eLeftOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_eRightOvere5x5",elec_full5x5_eRightOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_eTopOvere5x5",elec_full5x5_eTopOvere5x5,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalDepth1OverEcal",elec_full5x5_hcalDepth1OverEcal,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalDepth1OverEcalBc",elec_full5x5_hcalDepth1OverEcalBc,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalDepth2OverEcal",elec_full5x5_hcalDepth2OverEcal,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalDepth2OverEcalBc",elec_full5x5_hcalDepth2OverEcalBc,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalOverEcal",elec_full5x5_hcalOverEcal,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("full5x5_hcalOverEcalBc",elec_full5x5_hcalOverEcalBc,"doc",nanoaod::FlatTable::FloatColumn);   
      electronTable->addColumn<float>("full5x5_r9",elec_full5x5_r9,"doc",nanoaod::FlatTable::FloatColumn);

      electronTable->addColumn<float>("neutralHadronIso",elec_neutralHadronIso,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("particleIso",elec_particleIso,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("photonIso",elec_photonIso,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("puChargedHadronIso",elec_puChargedHadronIso,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("trackIso",elec_trackIso,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("hcalDepth1OverEcal",elec_hcalDepth1OverEcal,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("hcalDepth2OverEcal",elec_hcalDepth2OverEcal,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("ecalPFClusterIso",elec_ecalPFClusterIso,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("hcalPFClusterIso",elec_hcalPFClusterIso,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("dr03TkSumPt",elec_dr03TkSumPt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dr03HcalDepth2TowerSumEt",elec_dr03HcalDepth2TowerSumEt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dr03EcalRecHitSumEt",elec_dr03EcalRecHitSumEt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dr03HcalDepth1TowerSumEt",elec_dr03HcalDepth1TowerSumEt,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("dr03HcalDepth1TowerSumEtBc",elec_dr03HcalDepth1TowerSumEtBc,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("pfSumPhotonEt",elec_pfSumPhotonEt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("pfSumChargedHadronPt",elec_pfSumChargedHadronPt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("pfSumNeutralHadronEt",elec_pfSumNeutralHadronEt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("pfSumPUPt",elec_pfSumPUPt,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("dr04EcalRecHitSumEt",elec_dr04EcalRecHitSumEt,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("dr04HcalDepth1TowerSumEt",elec_dr04HcalDepth1TowerSumEt,"doc",nanoaod::FlatTable::FloatColumn);  
      electronTable->addColumn<float>("dr04HcalDepth1TowerSumEtBc",elec_dr04HcalDepth1TowerSumEtBc,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("dr04HcalDepth2TowerSumEt",elec_dr04HcalDepth2TowerSumEt,"doc",nanoaod::FlatTable::FloatColumn); 
      electronTable->addColumn<float>("dr04HcalDepth2TowerSumEtBc",elec_dr04HcalDepth2TowerSumEtBc,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("dr04HcalTowerSumEt",elec_dr04HcalTowerSumEt,"doc",nanoaod::FlatTable::FloatColumn);
      electronTable->addColumn<float>("dr04HcalTowerSumEtBc",elec_dr04HcalTowerSumEtBc,"doc",nanoaod::FlatTable::FloatColumn);



    iEvent.put(std::move(globalTable), "global");
    iEvent.put(std::move(csvTable), "csv");
    iEvent.put(std::move(cpfTable), "cpf");
    iEvent.put(std::move(npfTable), "npf");
    iEvent.put(std::move(svTable), "sv");
    iEvent.put(std::move(lengthTable), "length");
    iEvent.put(std::move(jetOriginTable), "jetorigin");
    iEvent.put(std::move(muonTable), "muon");
    iEvent.put(std::move(electronTable), "electron");
}

void
NANOProducer::beginStream(edm::StreamID)
{
}

void
NANOProducer::endStream() {
}

 
void
NANOProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NANOProducer);
