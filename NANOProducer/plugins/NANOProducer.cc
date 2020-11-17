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
#include <Math/Vector4D.h>

#include "FlatTableFiller.h"

//
// class declaration
//

class NANOProducer : public edm::stream::EDProducer<> {
    public:
        explicit NANOProducer(const edm::ParameterSet&);
        ~NANOProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        PropertyList<llpdnnx::JetFeatures> globalProperties;
        PropertyList<llpdnnx::ShallowTagInfoFeatures> csvProperties;
        PropertyList<llpdnnx::ChargedCandidateFeatures> cpfProperties;
        PropertyList<llpdnnx::NeutralCandidateFeatures> npfProperties;
        PropertyList<llpdnnx::SecondaryVertexFeatures> svProperties;
        PropertyList<llpdnnx::MuonCandidateFeatures> muonProperties;
        PropertyList<llpdnnx::ElectronCandidateFeatures> electronProperties;
        
   private:
      const edm::EDGetTokenT<edm::View<pat::Jet>> _jet_src;
      const edm::EDGetTokenT<std::vector<reco::XTagInfo>> _tag_src;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
};

NANOProducer::NANOProducer(const edm::ParameterSet& iConfig) :
    _jet_src(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJets"))),
    _tag_src(consumes<std::vector<reco::XTagInfo>>(iConfig.getParameter<edm::InputTag>("srcTags")))
{
    globalProperties = {
        PROPERTY(llpdnnx::JetFeatures, pt, "doc"),
        PROPERTY(llpdnnx::JetFeatures, eta, "doc"),
        PROPERTY(llpdnnx::JetFeatures, phi, "doc"),
        PROPERTY(llpdnnx::JetFeatures, mass, "doc"),
        PROPERTY(llpdnnx::JetFeatures, energy, "doc"),

        PROPERTY(llpdnnx::JetFeatures, area, "doc"),

        PROPERTY(llpdnnx::JetFeatures, beta, "doc"),
        PROPERTY(llpdnnx::JetFeatures, dR2Mean, "doc"),
        PROPERTY(llpdnnx::JetFeatures, frac01, "doc"),
        PROPERTY(llpdnnx::JetFeatures, frac02, "doc"),
        PROPERTY(llpdnnx::JetFeatures, frac03, "doc"),
        PROPERTY(llpdnnx::JetFeatures, frac04, "doc"),

        PROPERTY(llpdnnx::JetFeatures, jetR, "doc"),
        PROPERTY(llpdnnx::JetFeatures, jetRchg, "doc"),

        PROPERTY(llpdnnx::JetFeatures, n60, "doc"),
        PROPERTY(llpdnnx::JetFeatures, n90, "doc"),

        PROPERTY(llpdnnx::JetFeatures, chargedEmEnergyFraction, "doc"),
        PROPERTY(llpdnnx::JetFeatures, chargedHadronEnergyFraction, "doc"),
        PROPERTY(llpdnnx::JetFeatures, chargedMuEnergyFraction, "doc"),
        PROPERTY(llpdnnx::JetFeatures, electronEnergyFraction, "doc"),

        PROPERTY(llpdnnx::JetFeatures, tau1, "doc"),
        PROPERTY(llpdnnx::JetFeatures, tau2, "doc"),
        PROPERTY(llpdnnx::JetFeatures, tau3, "doc"),

        PROPERTY(llpdnnx::JetFeatures, relMassDropMassAK, "doc"),
        PROPERTY(llpdnnx::JetFeatures, relMassDropMassCA, "doc"),
        PROPERTY(llpdnnx::JetFeatures, relSoftDropMassAK, "doc"),
        PROPERTY(llpdnnx::JetFeatures, relSoftDropMassCA, "doc"),

        PROPERTY(llpdnnx::JetFeatures, thrust, "doc"),
        PROPERTY(llpdnnx::JetFeatures, sphericity, "doc"),
        PROPERTY(llpdnnx::JetFeatures, circularity, "doc"),
        PROPERTY(llpdnnx::JetFeatures, isotropy, "doc"),
        PROPERTY(llpdnnx::JetFeatures, eventShapeC, "doc"),
        PROPERTY(llpdnnx::JetFeatures, eventShapeD, "doc"),

        PROPERTY(llpdnnx::JetFeatures, numberCpf, "doc"),
        PROPERTY(llpdnnx::JetFeatures, numberNpf, "doc"),
        PROPERTY(llpdnnx::JetFeatures, numberSv, "doc"),
        PROPERTY(llpdnnx::JetFeatures, numberSvAdapted, "doc"),
        PROPERTY(llpdnnx::JetFeatures, numberMuon, "doc"),
        PROPERTY(llpdnnx::JetFeatures, numberElectron, "doc")
    
    };
    
    csvProperties = {
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSumJetEtRatio, "ratio of track sum transverse energy over jet energy"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSumJetDeltaR, "pseudoangular distance between jet axis and track fourvector sum"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, vertexCategory, "category of secondary vertex (Reco, Pseudo, No)"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSip2dValAboveCharm, "track 2D signed impact parameter of first track lifting mass above charm"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSip2dSigAboveCharm, "track 2D signed impact parameter significance of first track lifting mass above charm"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSip3dValAboveCharm, "track 3D signed impact parameter of first track lifting mass above charm"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, trackSip3dSigAboveCharm, "track 3D signed impact parameter significance of first track lifting mass above charm"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, jetNTracksEtaRel, "tracks associated to jet for which trackEtaRel is calculated"),
        PROPERTY(llpdnnx::ShallowTagInfoFeatures, jetNSelectedTracks, "doc")
    };
    
    cpfProperties = {
        PROPERTY(llpdnnx::ChargedCandidateFeatures, ptrel, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, deta, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, dphi, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, deltaR, "doc"),

        PROPERTY(llpdnnx::ChargedCandidateFeatures, px, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, py, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, pz, "doc"),

        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackEtaRel, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackPtRel, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackPPar, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackDeltaR, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackPParRatio, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackPtRatio, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dVal, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dSig, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dVal, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dSig, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackJetDistVal, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackJetDistSig, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, drminsv, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, vertex_association, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, fromPV, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, puppi_weight, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_chi2, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_quality, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_numberOfValidPixelHits, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_pixelLayersWithMeasurement, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_numberOfValidStripHits, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_stripLayersWithMeasurement, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, relmassdrop, "doc"),
        
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dValSV, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dSigSV, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dValSV, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dSigSV, "doc"), 

        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dValSV_adapted, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip2dSigSV_adapted, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dValSV_adapted, "doc"), 
        PROPERTY(llpdnnx::ChargedCandidateFeatures, trackSip3dSigSV_adapted, "doc"), 

        PROPERTY(llpdnnx::ChargedCandidateFeatures, matchedMuon, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, matchedElectron, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, matchedSV, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, matchedSV_adapted, "doc"),
        
        PROPERTY(llpdnnx::ChargedCandidateFeatures, track_ndof, "doc"),
        PROPERTY(llpdnnx::ChargedCandidateFeatures, dZmin, "doc")
    };
    
    
    npfProperties = {
        PROPERTY(llpdnnx::NeutralCandidateFeatures, ptrel, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, deta, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, dphi, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, deltaR, "doc"),

        PROPERTY(llpdnnx::NeutralCandidateFeatures, px, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, py, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, pz, "doc"),

        PROPERTY(llpdnnx::NeutralCandidateFeatures, isGamma, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, hcal_fraction, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, drminsv, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, puppi_weight, "doc"),
        PROPERTY(llpdnnx::NeutralCandidateFeatures, relmassdrop, "doc")
    };
    
    svProperties = {
        PROPERTY(llpdnnx::SecondaryVertexFeatures, ptrel, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, deta, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, dphi, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, deltaR, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, mass, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, ntracks, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, chi2, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, ndof, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, dxy, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, dxysig, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, d3d, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, d3dsig, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, costhetasvpv, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, enratio, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, vx, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, vy, "doc"),
        PROPERTY(llpdnnx::SecondaryVertexFeatures, vz, "doc")
    };
    
    muonProperties = {
        PROPERTY(llpdnnx::MuonCandidateFeatures, isGlobal, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, isTight, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, isMedium, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, isLoose, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, isStandAlone, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, ptrel, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, deta, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dphi, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, px, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, py, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, pz, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, charge, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, energy, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, et, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, deltaR, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, numberOfMatchedStations, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, IP2d, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, IP2dSig, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, IP3d, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, IP3dSig, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, EtaRel, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dxy, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dxyError, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dxySig, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dz, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dzError, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, dzSig, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, numberOfValidPixelHits, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, numberOfpixelLayersWithMeasurement, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, numberOfstripLayersWithMeasurement, "doc"), //that does not help. needs to be discussed.

        PROPERTY(llpdnnx::MuonCandidateFeatures, chi2, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, ndof, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, caloIso, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, ecalIso, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, hcalIso, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfChHadronPt, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfNeuHadronEt, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, Pfpileup, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfPhotonEt, "doc"),

        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfChHadronPt03, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfNeuHadronEt03, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, Pfpileup03, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, sumPfPhotonEt03, "doc"),
    
        PROPERTY(llpdnnx::MuonCandidateFeatures, timeAtIpInOut, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, timeAtIpInOutErr, "doc"),
        PROPERTY(llpdnnx::MuonCandidateFeatures, timeAtIpOutIn, "doc")
    };
    
    
    electronProperties = {
        PROPERTY(llpdnnx::ElectronCandidateFeatures, ptrel, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures, deltaR, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deta, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dphi, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,px, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,py, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,pz, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,charge, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,energy, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,EtFromCaloEn, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures, isEB, "doc"), 
        PROPERTY(llpdnnx::ElectronCandidateFeatures, isEE, "doc"), 
        PROPERTY(llpdnnx::ElectronCandidateFeatures,ecalEnergy, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures, isPassConversionVeto, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,convDist, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,   convFlags, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,convRadius, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,hadronicOverEm, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,ecalDrivenSeed, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,IP2d, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,IP2dSig, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,IP3d, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,IP3dSig, "doc"),

        PROPERTY(llpdnnx::ElectronCandidateFeatures,elecSC_energy, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,elecSC_deta, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,elecSC_dphi, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,elecSC_et, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,elecSC_eSuperClusterOverP, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,scPixCharge, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,superClusterFbrem, "doc"),

        PROPERTY(llpdnnx::ElectronCandidateFeatures,eSeedClusterOverP, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,eSeedClusterOverPout, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,eSuperClusterOverP, "doc"),

        // shower shape
        PROPERTY(llpdnnx::ElectronCandidateFeatures,sigmaEtaEta, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,sigmaIetaIeta, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,sigmaIphiIphi, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,e5x5, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,e5x5Rel, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,e1x5Overe5x5, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,e2x5MaxOvere5x5, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,r9, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,hcalOverEcal, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,hcalDepth1OverEcal, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,hcalDepth2OverEcal, "doc"),

        // Track-Cluster Matching Attributes
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaEtaEleClusterTrackAtCalo, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaEtaSeedClusterTrackAtCalo, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaPhiSeedClusterTrackAtCalo, "doc"), 
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaEtaSeedClusterTrackAtVtx, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaEtaSuperClusterTrackAtVtx, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaPhiEleClusterTrackAtCalo, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,deltaPhiSuperClusterTrackAtVtx, "doc"),

        PROPERTY(llpdnnx::ElectronCandidateFeatures,sCseedEta, "doc"),

        // electron gsf variables. 
        PROPERTY(llpdnnx::ElectronCandidateFeatures,EtaRel, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dxy, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dxyError, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dxySig, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dz, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dzError, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dzSig, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures, nbOfMissingHits, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,gsfCharge, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,   ndof, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,chi2, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,   numberOfBrems, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,fbrem, "doc"),

        // Isolation block
        PROPERTY(llpdnnx::ElectronCandidateFeatures,neutralHadronIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,particleIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,photonIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,puChargedHadronIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,trackIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,ecalPFClusterIso, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,hcalPFClusterIso, "doc"),

        PROPERTY(llpdnnx::ElectronCandidateFeatures,pfSumPhotonEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,pfSumChargedHadronPt, "doc"), 
        PROPERTY(llpdnnx::ElectronCandidateFeatures,pfSumNeutralHadronEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,pfSumPUPt, "doc"),

        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04TkSumPt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04EcalRecHitSumEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalDepth1TowerSumEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalDepth1TowerSumEtBc, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalDepth2TowerSumEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalDepth2TowerSumEtBc, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalTowerSumEt, "doc"),
        PROPERTY(llpdnnx::ElectronCandidateFeatures,dr04HcalTowerSumEtBc, "doc")
    
    };


    produces<nanoaod::FlatTable>("global");
    produces<nanoaod::FlatTable>("csv");
    produces<nanoaod::FlatTable>("cpf");
    produces<nanoaod::FlatTable>("npf");
    produces<nanoaod::FlatTable>("sv");
    produces<nanoaod::FlatTable>("svAdapted");
    produces<nanoaod::FlatTable>("length");
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
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(_jet_src, jets);
    edm::Handle<std::vector<reco::XTagInfo>> tag_infos;
    iEvent.getByToken(_tag_src, tag_infos);

    unsigned int ntags = tag_infos->size();

    std::vector<int> global_jetIdx;
    std::vector<int> csv_jetIdx;
    std::vector<int> cpf_jetIdx;
    std::vector<int> npf_jetIdx;
    std::vector<int> sv_jetIdx;
    std::vector<int> sv_adapted_jetIdx;
    std::vector<int> mu_jetIdx;
    std::vector<int> elec_jetIdx;

    auto lengthTable = std::make_unique<nanoaod::FlatTable>(ntags, "length", false, false);
    std::vector<int> cpf_length;
    std::vector<int> npf_length;
    std::vector<int> sv_length;
    std::vector<int> sv_adapted_length;

    std::vector<int> elec_length;
    std::vector<int> mu_length;

    auto globalTable = std::make_unique<nanoaod::FlatTable>(ntags, "global", false, false);
    auto csvTable = std::make_unique<nanoaod::FlatTable>(ntags, "csv", false, false);

    

    FlatTableFillerList<llpdnnx::JetFeatures> globalFillerList(globalProperties);
    FlatTableFillerList<llpdnnx::ShallowTagInfoFeatures> csvFillerList(csvProperties);
    
    FlatTableFillerList<llpdnnx::ChargedCandidateFeatures> cpfFillerList(cpfProperties);
    FlatTableFillerList<llpdnnx::NeutralCandidateFeatures> npfFillerList(npfProperties);
    FlatTableFillerList<llpdnnx::SecondaryVertexFeatures> svFillerList(svProperties);
    FlatTableFillerList<llpdnnx::SecondaryVertexFeatures> svAdaptedFillerList(svProperties);
    FlatTableFillerList<llpdnnx::MuonCandidateFeatures> muonFillerList(muonProperties);
    FlatTableFillerList<llpdnnx::ElectronCandidateFeatures> electronFillerList(electronProperties);

    unsigned int nmu_total = 0;
    unsigned int nelec_total = 0;
    unsigned int ncpf_total = 0;
    unsigned int nnpf_total = 0;
    unsigned int nsv_total = 0;
    unsigned int nsv_total_adapted = 0;

    for (unsigned int itag= 0; itag < ntags; itag++){
        const auto& features = tag_infos->at(itag).features();

        unsigned int nmu = features.mu_features.size();
        unsigned int nelec = features.elec_features.size();
        unsigned int ncpf = features.cpf_features.size();
        unsigned int nnpf = features.npf_features.size();
        unsigned int nsv = features.sv_features.size();
        unsigned int nsv_adapted = features.sv_adapted_features.size();
    
        nmu_total += nmu;
        nelec_total += nelec;
        ncpf_total += ncpf;
        nnpf_total += nnpf;
        nsv_total += nsv;
        nsv_total_adapted += nsv_adapted;

        cpf_length.push_back(ncpf);
        npf_length.push_back(nnpf);
        sv_length.push_back(nsv);
        sv_adapted_length.push_back(nsv_adapted);
        elec_length.push_back(nelec);
        mu_length.push_back(nmu);

        int jetIdx = -1;
        auto base_jet_ref = tag_infos->at(itag).jet();
        if (base_jet_ref.isAvailable() and base_jet_ref.isNonnull()){
            const auto& base_jet = base_jet_ref.get();
            for (std::size_t ijet = 0; ijet < jets->size(); ijet++) {
                auto jet = jets->at(ijet);
                if (reco::deltaR(base_jet->p4(),jet.p4()) < 1e-4){
                    jetIdx = ijet;
                    break;
                }
            }
        } 
        
        global_jetIdx.push_back(jetIdx);
        csv_jetIdx.push_back(jetIdx);
        
        globalFillerList.push_back(features.jet_features);
        csvFillerList.push_back(features.tag_info_features);
    }

    auto muonTable = std::make_unique<nanoaod::FlatTable>(nmu_total, "muon", false, false);
    auto electronTable = std::make_unique<nanoaod::FlatTable>(nelec_total, "electron", false, false);
    auto cpfTable = std::make_unique<nanoaod::FlatTable>(ncpf_total, "cpf", false, false);
    auto npfTable = std::make_unique<nanoaod::FlatTable>(nnpf_total, "npf", false, false);
    auto svTable = std::make_unique<nanoaod::FlatTable>(nsv_total, "sv", false, false);
    auto svTable_adapted = std::make_unique<nanoaod::FlatTable>(nsv_total_adapted, "svAdapted", false, false);

    for (unsigned int itag= 0; itag < ntags; itag++)
    {
        const auto& features = tag_infos->at(itag).features();
        auto mu  = features.mu_features;
        auto elec = features.elec_features;
        auto cpf = features.cpf_features;
        auto npf = features.npf_features;
        auto sv = features.sv_features;
        auto sv_adapted = features.sv_adapted_features;

        unsigned int nmu = features.mu_features.size();
        unsigned int nelec = features.elec_features.size();
        unsigned int ncpf = features.cpf_features.size();
        unsigned int nnpf = features.npf_features.size();
        unsigned int nsv = features.sv_features.size();
        unsigned int nsv_adapted = features.sv_adapted_features.size();

        int jetIdx = global_jetIdx[itag];

        for (unsigned int i = 0; i < ncpf; i++)
        {
            cpf_jetIdx.push_back(jetIdx);
            cpfFillerList.push_back(cpf.at(i));
        }

        for (unsigned int i = 0; i < nnpf; i++)
        {
            npf_jetIdx.push_back(jetIdx);          
            npfFillerList.push_back(npf.at(i));
        }

        for (unsigned int i = 0; i < nsv; i++)
        {
            sv_jetIdx.push_back(jetIdx);
            svFillerList.push_back(sv.at(i));
        }

        for (unsigned int i = 0; i < nsv_adapted; i++)
        {
            sv_adapted_jetIdx.push_back(jetIdx);
            svAdaptedFillerList.push_back(sv_adapted.at(i));
        }

        for(unsigned int i = 0; i < nmu; i++)
        {
            mu_jetIdx.push_back(jetIdx);
            muonFillerList.push_back(mu.at(i));
        }

        for(unsigned int i = 0; i < nelec ; i++)
        {
            elec_jetIdx.push_back(jetIdx);
            electronFillerList.push_back(elec.at(i));
        }

    }

 
    lengthTable->addColumn<int>("cpf", cpf_length, "charged PF candidate track offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("npf", npf_length, "neutral PF candidate offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("sv", sv_length, "secondary vertex (SV) offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("svAdapted", sv_adapted_length, "secondary vertex (SV) offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("mu", mu_length, "muon offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("ele", elec_length, "electron offset", nanoaod::FlatTable::IntColumn);
    
    
    globalTable->addColumn<int>("jetIdx", global_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    globalFillerList.fill(globalTable);
    
    csvTable->addColumn<int>("jetIdx", csv_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    csvFillerList.fill(csvTable);
   
    cpfTable->addColumn<int>("jetIdx", cpf_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    cpfFillerList.fill(cpfTable);

    npfTable->addColumn<int>("jetIdx", npf_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    npfFillerList.fill(npfTable);

    svTable->addColumn<int>("jetIdx", sv_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    svFillerList.fill(svTable);
    
    svTable_adapted->addColumn<int>("jetIdx", sv_adapted_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    svAdaptedFillerList.fill(svTable_adapted);

    muonTable->addColumn<int>("jetIdx", mu_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    muonFillerList.fill(muonTable);

    electronTable->addColumn<int>("jetIdx", elec_jetIdx, "linked jet Idx", nanoaod::FlatTable::IntColumn);
    electronFillerList.fill(electronTable);
    
    iEvent.put(std::move(globalTable), "global");
    iEvent.put(std::move(csvTable), "csv");
    iEvent.put(std::move(cpfTable), "cpf");
    iEvent.put(std::move(npfTable), "npf");
    iEvent.put(std::move(svTable), "sv");
    iEvent.put(std::move(svTable_adapted), "svAdapted");
    iEvent.put(std::move(lengthTable), "length");
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
