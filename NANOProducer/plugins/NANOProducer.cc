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



//
// class declaration
//

class NANOProducer : public edm::stream::EDProducer<> {
   public:
      explicit NANOProducer(const edm::ParameterSet&);
      ~NANOProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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
    produces<nanoaod::FlatTable>("global");
    produces<nanoaod::FlatTable>("csv");
    produces<nanoaod::FlatTable>("cpf");
    produces<nanoaod::FlatTable>("npf");
    produces<nanoaod::FlatTable>("sv");
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
    std::vector<int> cpf_jetIdx;
    std::vector<int> npf_jetIdx;
    std::vector<int> sv_jetIdx;
    std::vector<int> mu_jetIdx;
    std::vector<int> elec_jetIdx;

    auto lengthTable = std::make_unique<nanoaod::FlatTable>(ntags, "length", false, false);
    std::vector<int> cpf_length;
    std::vector<int> npf_length;
    std::vector<int> sv_length;

    std::vector<int> elec_length;
    std::vector<int> mu_length;

    auto globalTable = std::make_unique<nanoaod::FlatTable>(ntags, "global", false, false);
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

    std::vector<float> area;

    std::vector<int> n60;
    std::vector<int> n90;

    std::vector<float> beta;
    std::vector<float> dR2Mean;
    std::vector<float> frac01;
    std::vector<float> frac02;
    std::vector<float> frac03;
    std::vector<float> frac04;

    std::vector<float> jetR;
    std::vector<float> jetRchg;

    std::vector<float> chargedEmEnergyFraction;
    std::vector<float> chargedHadronEnergyFraction;
    std::vector<float> chargedMuEnergyFraction;
    std::vector<float> electronEnergyFraction;

    std::vector<float> tau1;
    std::vector<float> tau2;
    std::vector<float> tau3;

    std::vector<float> relMassDropMassAK;
    std::vector<float> relMassDropMassCA;
    std::vector<float> relSoftDropMassAK;
    std::vector<float> relSoftDropMassCA;

    std::vector<float> thrust;
    std::vector<float> sphericity;
    std::vector<float> circularity;
    std::vector<float> isotropy;
    std::vector<float> eventShapeC;
    std::vector<float> eventShapeD;


    auto csvTable = std::make_unique<nanoaod::FlatTable>(ntags, "csv", false, false);
    std::vector<float> trackSumJetEtRatio;
    std::vector<float> trackSumJetDeltaR;
    std::vector<float> vertexCategory;
    std::vector<float> trackSip2dValAboveCharm;
    std::vector<float> trackSip2dSigAboveCharm;
    std::vector<float> trackSip3dValAboveCharm;
    std::vector<float> trackSip3dSigAboveCharm;
    std::vector<int> jetNSelectedTracks;
    std::vector<int> jetNTracksEtaRel;

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
    std::vector<float> cpf_deta;
    std::vector<float> cpf_dphi;
    std::vector<float> cpf_drminsv;
    std::vector<float> cpf_vertex_association;
    std::vector<float> cpf_fromPV;
    std::vector<float> cpf_puppi_weight;
    std::vector<float> cpf_track_chi2;
    std::vector<float> cpf_track_quality;
    std::vector<float> cpf_relmassdrop;
    std::vector<int> cpf_track_ndof;
    std::vector<int> cpf_matchedMuon;
    std::vector<int> cpf_matchedElectron;
    std::vector<int> cpf_matchedSV;
    std::vector<int> cpf_track_numberOfValidPixelHits;
    std::vector<int> cpf_track_pixelLayersWithMeasurement;
    std::vector<int> cpf_track_numberOfValidStripHits; 
    std::vector<int> cpf_track_stripLayersWithMeasurement;
    std::vector<float> cpf_dZmin;

    std::vector<float> npf_ptrel;
    std::vector<float> npf_deta;
    std::vector<float> npf_dphi;
    std::vector<float> npf_deltaR;
    std::vector<float> npf_isGamma;
    std::vector<float> npf_hcal_fraction;
    std::vector<float> npf_drminsv;
    std::vector<float> npf_puppi_weight;
    std::vector<float> npf_relmassdrop;

    std::vector<float> sv_ptrel;
    std::vector<float> sv_deta;
    std::vector<float> sv_dphi;
    std::vector<float> sv_deltaR;
    std::vector<float> sv_mass;
    std::vector<float> sv_ntracks;
    std::vector<float> sv_chi2;
    std::vector<int> sv_ndof;
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


    std::vector<float> mu_ptrel;
    std::vector<float> mu_EtaRel;
    std::vector<float> mu_deta;
    std::vector<float> mu_dphi;
    std::vector<float> mu_charge;
    std::vector<float> mu_energy;
    std::vector<float> mu_et;
    std::vector<float> mu_jetDeltaR;
    std::vector<float> mu_numberOfMatchedStations;

    std::vector<float> mu_2dIP;
    std::vector<float> mu_2dIPSig;
    std::vector<float> mu_3dIP;
    std::vector<float> mu_3dIPSig;

    std::vector<float> mu_dxy;
    std::vector<float> mu_dxyError;
    std::vector<float> mu_dxySig;
    std::vector<float> mu_dz;
    std::vector<float> mu_dzError;
    std::vector<int> mu_numberOfValidPixelHits;
    std::vector<int> mu_numberOfpixelLayersWithMeasurement;
    std::vector<int> mu_numberOfstripLayersWithMeasurement ;

    std::vector<float> mu_chi2;
    std::vector<int> mu_ndof;

    std::vector<float> mu_caloIso;
    std::vector<float> mu_ecalIso;
    std::vector<float> mu_hcalIso;

    std::vector<float> mu_sumPfChHadronPt;
    std::vector<float> mu_sumPfNeuHadronEt;
    std::vector<float> mu_Pfpileup;
    std::vector<float> mu_sumPfPhotonEt;

    std::vector<float> mu_sumPfChHadronPt03;
    std::vector<float> mu_sumPfNeuHadronEt03;
    std::vector<float> mu_Pfpileup03;
    std::vector<float> mu_sumPfPhotonEt03;



    std::vector<float>  mu_timeAtIpInOut;
    std::vector<float>  mu_timeAtIpInOutErr;
    std::vector<float>  mu_timeAtIpOutIn;


    // Electron Block

    std::vector<float> elec_ptrel;
    std::vector<float> elec_jetDeltaR;
    std::vector<float> elec_deta;
    std::vector<float> elec_dphi;
    std::vector<float> elec_charge;
    std::vector<float> elec_energy;
    std::vector<float> elec_EtFromCaloEn ;
    std::vector<float> elec_isEB;
    std::vector<float> elec_isEE;
    std::vector<float> elec_ecalEnergy;
    std::vector<float> elec_isPassConversionVeto;
    std::vector<float> elec_convDist;
    std::vector<int>   elec_convFlags;

    std::vector<float> elec_convRadius;
    std::vector<float> elec_hadronicOverEm;
    std::vector<float> elec_ecalDrivenSeed;

    std::vector<float> elecSC_energy;
    std::vector<float> elecSC_deta;
    std::vector<float> elecSC_dphi;
    std::vector<float> elecSC_et;
    std::vector<float> elecSC_eSuperClusterOverP;
    std::vector<float> elec_scE1x5Overe5x5;
    std::vector<float> elec_scE2x5MaxOvere5x5;
    std::vector<float> elec_scE5x5;
    std::vector<float> elec_scE5x5Rel;
    std::vector<float> elec_scPixCharge;
    std::vector<float> elec_scSigmaEtaEta;
    std::vector<float> elec_scSigmaIEtaIEta;
    std::vector<float> elec_superClusterFbrem;

    std::vector<float> elec_2dIP;
    std::vector<float> elec_2dIPSig;
    std::vector<float> elec_3dIP;
    std::vector<float> elec_3dIPSig;
    std::vector<float> elec_eSeedClusterOverP;
    std::vector<float> elec_eSeedClusterOverPout;
    std::vector<float> elec_eSuperClusterOverP;
    std::vector<float> elec_eTopOvere5x5;

    std::vector<float> elec_deltaEtaEleClusterTrackAtCalo;
    std::vector<float> elec_deltaEtaSeedClusterTrackAtCalo;
    std::vector<float> elec_deltaPhiSeedClusterTrackAtCalo;
    std::vector<float> elec_deltaEtaSeedClusterTrackAtVtx;
    std::vector<float> elec_deltaEtaSuperClusterTrackAtVtx;
    std::vector<float> elec_deltaPhiEleClusterTrackAtCalo;
    std::vector<float> elec_deltaPhiSuperClusterTrackAtVtx;
    std::vector<float> elec_sCseedEta;
    ///////
    std::vector<float> elec_EtaRel;
    std::vector<float> elec_dxy;
    std::vector<float> elec_dz;
    std::vector<float> elec_nbOfMissingHits;
    std::vector<float> elec_gsfCharge;
    std::vector<int> elec_ndof;
    std::vector<float> elec_chi2;

    std::vector<float> elec_e2x5MaxOvere5x5;
    std::vector<float> elec_e1x5Overe5x5;
    std::vector<float> elec_e5x5;
    std::vector<float> elec_e5x5Rel;
    std::vector<float> elec_full5x5_sigmaIetaIeta;
    std::vector<float> elec_full5x5_e1x5Overe5x5;
    std::vector<float> elec_full5x5_e2x5BottomOvere5x5;
    std::vector<float> elec_full5x5_e2x5LeftOvere5x5;
    std::vector<float> elec_full5x5_e2x5MaxOvere5x5;
    std::vector<float> elec_full5x5_e2x5RightOvere5x5;
    std::vector<float> elec_full5x5_e2x5TopOvere5x5;
    std::vector<float> elec_full5x5_e5x5;
    std::vector<float> elec_full5x5_e5x5Rel;
    std::vector<float> elec_full5x5_eBottomOvere5x5;
    std::vector<float> elec_full5x5_eLeftOvere5x5;
    std::vector<float> elec_full5x5_eRightOvere5x5;
    std::vector<float> elec_full5x5_eTopOvere5x5;

    std::vector<float> elec_full5x5_hcalDepth1OverEcal;
    std::vector<float> elec_full5x5_hcalDepth1OverEcalBc;
    std::vector<float> elec_full5x5_hcalDepth2OverEcal;
    std::vector<float> elec_full5x5_hcalDepth2OverEcalBc;
    std::vector<float> elec_full5x5_hcalOverEcal;
    std::vector<float> elec_full5x5_hcalOverEcalBc;
    std::vector<float> elec_full5x5_r9;
    std::vector<int>   elec_numberOfBrems;
    std::vector<float> elec_fbrem;
    std::vector<float> elec_neutralHadronIso;
    std::vector<float> elec_particleIso;
    std::vector<float> elec_photonIso;
    std::vector<float> elec_puChargedHadronIso;
    std::vector<float> elec_trackIso;
    std::vector<float> elec_hcalDepth1OverEcal;
    std::vector<float> elec_hcalDepth2OverEcal;
    std::vector<float> elec_ecalPFClusterIso;
    std::vector<float> elec_hcalPFClusterIso;
    std::vector<float> elec_dr03TkSumPt;
    std::vector<float> elec_dr03EcalRecHitSumEt;
    std::vector<float> elec_dr03HcalDepth1TowerSumEt;
    std::vector<float> elec_dr03HcalDepth1TowerSumEtBc;
    std::vector<float> elec_dr03HcalDepth2TowerSumEt;
    std::vector<float> elec_dr03HcalDepth2TowerSumEtBc;
    std::vector<float> elec_pfSumPhotonEt;
    std::vector<float> elec_pfSumChargedHadronPt;
    std::vector<float> elec_pfSumNeutralHadronEt;
    std::vector<float> elec_pfSumPUPt;
    std::vector<float> elec_dr04EcalRecHitSumEt;
    std::vector<float> elec_dr04HcalDepth1TowerSumEt;
    std::vector<float> elec_dr04HcalDepth1TowerSumEtBc;
    std::vector<float> elec_dr04HcalDepth2TowerSumEt;
    std::vector<float> elec_dr04HcalDepth2TowerSumEtBc;
    std::vector<float> elec_dr04HcalTowerSumEt;
    std::vector<float> elec_dr04HcalTowerSumEtBc;

    unsigned int nmu_total = 0;
    unsigned int nelec_total = 0;
    unsigned int ncpf_total = 0;
    unsigned int nnpf_total = 0;
    unsigned int nsv_total = 0;

    for (unsigned int itag= 0; itag < ntags; itag++){
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
        pt.push_back(features.jet_features.pt);
        eta.push_back(features.jet_features.eta);
        phi.push_back(features.jet_features.phi);
        mass.push_back(features.jet_features.mass);

        area.push_back(features.jet_features.area);

        n60.push_back(features.jet_features.n60);
        n90.push_back(features.jet_features.n90);

        beta.push_back(features.jet_features.beta);
        dR2Mean.push_back(features.jet_features.dR2Mean);
        frac01.push_back(features.jet_features.frac01);
        frac02.push_back(features.jet_features.frac02);
        frac03.push_back(features.jet_features.frac03);
        frac04.push_back(features.jet_features.frac04);
    
        jetR.push_back(features.jet_features.jetR);
        jetRchg.push_back(features.jet_features.jetRchg);

        chargedEmEnergyFraction.push_back(features.jet_features.chargedEmEnergyFraction);
        chargedHadronEnergyFraction.push_back(features.jet_features.chargedHadronEnergyFraction);
        chargedMuEnergyFraction.push_back(features.jet_features.chargedMuEnergyFraction);
        electronEnergyFraction.push_back(features.jet_features.electronEnergyFraction);

        tau1.push_back(features.jet_features.tau1);
        tau2.push_back(features.jet_features.tau2);
        tau3.push_back(features.jet_features.tau3);

        relMassDropMassAK.push_back(features.jet_features.relMassDropMassAK);
        relMassDropMassCA.push_back(features.jet_features.relMassDropMassCA);
        relSoftDropMassAK.push_back(features.jet_features.relSoftDropMassAK);
        relSoftDropMassCA.push_back(features.jet_features.relSoftDropMassCA);

        thrust.push_back(features.jet_features.thrust);
        sphericity.push_back(features.jet_features.sphericity);
        circularity.push_back(features.jet_features.circularity);
        isotropy.push_back(features.jet_features.isotropy);
        eventShapeC.push_back(features.jet_features.eventShapeC);
        eventShapeD.push_back(features.jet_features.eventShapeD);

        trackSumJetEtRatio.push_back(tag_info_features.csv_trackSumJetEtRatio);
        trackSumJetDeltaR.push_back(tag_info_features.csv_trackSumJetDeltaR);
        trackSip2dValAboveCharm.push_back(tag_info_features.csv_trackSip2dValAboveCharm);
        trackSip2dSigAboveCharm.push_back(tag_info_features.csv_trackSip2dSigAboveCharm);
        trackSip3dValAboveCharm.push_back(tag_info_features.csv_trackSip3dValAboveCharm);
        trackSip3dSigAboveCharm.push_back(tag_info_features.csv_trackSip3dSigAboveCharm);
        jetNSelectedTracks.push_back(tag_info_features.csv_jetNSelectedTracks);
        jetNTracksEtaRel.push_back(tag_info_features.csv_jetNTracksEtaRel);
        vertexCategory.push_back(tag_info_features.csv_vertexCategory);
    }

    auto muonTable = std::make_unique<nanoaod::FlatTable>(nmu_total, "muon", false, false);
    auto electronTable = std::make_unique<nanoaod::FlatTable>(nelec_total, "electron", false, false);
    auto cpfTable = std::make_unique<nanoaod::FlatTable>(ncpf_total, "cpf", false, false);
    auto npfTable = std::make_unique<nanoaod::FlatTable>(nnpf_total, "npf", false, false);
    auto svTable = std::make_unique<nanoaod::FlatTable>(nsv_total, "sv", false, false);

    for (unsigned int itag= 0; itag < ntags; itag++)
    {
        const auto& features = tag_infos->at(itag).features();
        auto mu  = features.mu_features;
        auto elec = features.elec_features;
        auto cpf = features.cpf_features;
        auto npf = features.npf_features;
        auto sv = features.sv_features;

        unsigned int nmu = features.mu_features.size();
        unsigned int nelec = features.elec_features.size();
        unsigned int ncpf = features.cpf_features.size();
        unsigned int nnpf = features.npf_features.size();
        unsigned int nsv = features.sv_features.size();

        int jetIdx = global_jetIdx[itag];

        for (unsigned int i = 0; i < ncpf; i++)
        {
            const auto& cpf_features = cpf.at(i);
            cpf_jetIdx.push_back(jetIdx);
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
            cpf_deta.push_back(cpf_features.cpf_deta);
            cpf_dphi.push_back(cpf_features.cpf_dphi);
            cpf_drminsv.push_back(cpf_features.cpf_drminsv);
            cpf_vertex_association.push_back(cpf_features.cpf_vertex_association);
            cpf_fromPV.push_back(cpf_features.cpf_fromPV);
            cpf_puppi_weight.push_back(cpf_features.cpf_puppi_weight);
            cpf_track_chi2.push_back(cpf_features.cpf_track_chi2);
            cpf_track_quality.push_back(cpf_features.cpf_track_quality);
            cpf_relmassdrop.push_back(cpf_features.cpf_relmassdrop);
            cpf_track_ndof.push_back(cpf_features.cpf_track_ndof);
            cpf_matchedMuon.push_back(cpf_features.cpf_matchedMuon);
            cpf_matchedElectron.push_back(cpf_features.cpf_matchedElectron);
            cpf_matchedSV.push_back(cpf_features.cpf_matchedSV);
	        cpf_track_numberOfValidPixelHits.push_back(cpf_features.cpf_track_numberOfValidPixelHits);
            cpf_track_pixelLayersWithMeasurement.push_back(cpf_features.cpf_track_pixelLayersWithMeasurement);
	        cpf_track_numberOfValidStripHits.push_back(cpf_features.cpf_track_numberOfValidStripHits);
	        cpf_track_stripLayersWithMeasurement.push_back(cpf_features.cpf_track_stripLayersWithMeasurement);
            cpf_dZmin.push_back(cpf_features.cpf_dZmin);
        }

        for (unsigned int i = 0; i < nnpf; i++)
        {
            const auto& npf_features = npf.at(i);
            npf_jetIdx.push_back(jetIdx);
            npf_ptrel.push_back(npf_features.npf_ptrel);
            npf_deta.push_back(npf_features.npf_deta);
            npf_dphi.push_back(npf_features.npf_dphi);
            npf_deltaR.push_back(npf_features.npf_deltaR);
            npf_isGamma.push_back(npf_features.npf_isGamma);
            npf_hcal_fraction.push_back(npf_features.npf_hcal_fraction);
            npf_drminsv.push_back(npf_features.npf_drminsv);
            npf_puppi_weight.push_back(npf_features.npf_puppi_weight);
            npf_relmassdrop.push_back(npf_features.npf_relmassdrop);
        }

        for (unsigned int i = 0; i < nsv; i++)
        {
            const auto& sv_features = sv.at(i);
            sv_jetIdx.push_back(jetIdx);
            sv_ptrel.push_back(sv_features.sv_ptrel);
            sv_deta.push_back(sv_features.sv_deta);
            sv_dphi.push_back(sv_features.sv_dphi);
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


        for(unsigned int i = 0; i < nmu; i++)
        {
            const auto& mu_features = mu.at(i);
            mu_jetIdx.push_back(jetIdx);
            mu_isGlobal.push_back(mu_features.mu_isGlobal) ;
            mu_isTight.push_back(mu_features.mu_isTight) ;
            mu_isMedium.push_back(mu_features.mu_isMedium) ;
            mu_isLoose.push_back(mu_features.mu_isLoose) ;
            mu_isStandAlone.push_back(mu_features.mu_isStandAlone) ;


            mu_ptrel.push_back(mu_features.mu_ptrel);
            mu_EtaRel.push_back(mu_features.mu_EtaRel);
            mu_deta.push_back(mu_features.mu_deta);
            mu_dphi.push_back(mu_features.mu_dphi);
            mu_charge.push_back(mu_features.mu_charge);
            mu_energy.push_back(mu_features.mu_energy);
            mu_et.push_back(mu_features.mu_et);
            mu_jetDeltaR.push_back(mu_features.mu_jetDeltaR);
            mu_numberOfMatchedStations.push_back(mu_features.mu_numberOfMatchedStations);

            mu_2dIP.push_back(mu_features.mu_2dIP);
            mu_2dIPSig.push_back(mu_features.mu_2dIPSig);
            mu_3dIP.push_back(mu_features.mu_3dIP);
            mu_3dIPSig.push_back(mu_features.mu_3dIPSig);

            mu_dxy.push_back(mu_features.mu_dxy);
            mu_dxyError.push_back(mu_features.mu_dxyError);
            mu_dxySig.push_back(mu_features.mu_dxySig);
            mu_dz.push_back(mu_features.mu_dz);
            mu_dzError.push_back(mu_features.mu_dzError);
            mu_numberOfValidPixelHits.push_back(mu_features.mu_numberOfValidPixelHits);
            mu_numberOfpixelLayersWithMeasurement.push_back(mu_features.mu_numberOfpixelLayersWithMeasurement);
            mu_numberOfstripLayersWithMeasurement.push_back(mu_features.mu_numberOfstripLayersWithMeasurement);


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


            mu_timeAtIpInOut.push_back(mu_features.mu_timeAtIpInOut) ;
            mu_timeAtIpInOutErr.push_back(mu_features.mu_timeAtIpInOutErr) ;
            mu_timeAtIpOutIn.push_back(mu_features.mu_timeAtIpOutIn) ;


        }

        for(unsigned int i = 0; i < nelec ; i++)
        {

            const auto& elec_features = elec.at(i);

            elec_jetIdx.push_back(jetIdx);

            elec_ptrel.push_back(elec_features.elec_ptrel );
            elec_jetDeltaR.push_back(elec_features.elec_jetDeltaR );
            elec_deta.push_back(elec_features.elec_deta );
            elec_dphi.push_back(elec_features.elec_dphi );
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
            elecSC_deta.push_back(elec_features.elecSC_deta );
            elecSC_dphi.push_back(elec_features.elecSC_dphi );
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
            elec_ndof.push_back(elec_features.elec_ndof );
            elec_chi2.push_back(elec_features.elec_chi2 );

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
            elec_dr03EcalRecHitSumEt.push_back(elec_features.elec_dr03EcalRecHitSumEt );
            elec_dr03HcalDepth1TowerSumEt.push_back(elec_features.elec_dr03HcalDepth1TowerSumEt );
            elec_dr03HcalDepth1TowerSumEtBc.push_back(elec_features.elec_dr03HcalDepth1TowerSumEtBc );
            elec_dr03HcalDepth2TowerSumEt.push_back(elec_features.elec_dr03HcalDepth2TowerSumEt );
            elec_dr03HcalDepth2TowerSumEtBc.push_back(elec_features.elec_dr03HcalDepth2TowerSumEtBc );
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
    globalTable->addColumn<int>("jetIdx", global_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);
    globalTable->addColumn<float>("pt", pt, "global jet pt (uncorrected)", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("eta", eta, "global jet eta", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("phi", phi, "global jet phi", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("mass", mass, "global jet mass", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("area", area, "global jet area", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<int>("n60", n60, "n60", nanoaod::FlatTable::IntColumn);
    globalTable->addColumn<int>("n90", n90, "n90", nanoaod::FlatTable::IntColumn);
    globalTable->addColumn<float>("chargedEmEnergyFraction", chargedEmEnergyFraction, "chargedEmEnergyFraction", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("chargedHadronEnergyFraction", chargedHadronEnergyFraction, "chargedHadronEnergyFraction", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("chargedMuEnergyFraction", chargedMuEnergyFraction, "chargedMuEnergyFraction", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("electronEnergyFraction", electronEnergyFraction, "electronEnergyFraction", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("tau1", tau1, "nsubjettiness 1", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("tau2", tau2, "nsubjettiness 2", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("tau3", tau3, "nsubjettiness 3", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("relMassDropMassAK", relMassDropMassAK, "mass drop mass with anti-kT", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("relMassDropMassCA", relMassDropMassCA, "mass drop mass with Cambridge/Aachen", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("relSoftDropMassAK", relSoftDropMassAK, "soft drop mass with anti-kT", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("relSoftDropMassCA", relSoftDropMassCA, "soft drop mass with Cambridge/Aachen", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("thrust", thrust, "thrust", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("sphericity", sphericity, "sphericity", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("circularity", circularity, "circularity", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("isotropy", isotropy, "isotropy", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("eventShapeC", eventShapeC, "eventShapeC", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("eventShapeD", eventShapeD, "eventShapeD", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("beta", beta, "fraction of tracks associated to the PV", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("dR2Mean", dR2Mean, "pt2 weighted average square distance of jet constituents from the jet axis", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("frac01", frac01, "fraction of constituents in the region dR < 0.1 around the jet axis", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("frac02", frac02, "fraction of constituents in the region 0.1 < dR < 0.2 around the jet axis", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("frac03", frac03, "fraction of constituents in the region 0.2 < dR < 0.3 around the jet axis", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("frac04", frac04, "fraction of constituents in the region 0.3 < dR < 0.4 around the jet axis", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("jetR", jetR, "fraction of jet pt carried by lead constituent", nanoaod::FlatTable::FloatColumn);
    globalTable->addColumn<float>("jetRchg", jetRchg, "fraction of jet pt carried by lead jet constituent", nanoaod::FlatTable::FloatColumn);  

    csvTable->addColumn<float>("trackSumJetEtRatio", trackSumJetEtRatio, "ratio of track sum transverse energy over jet energy", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSumJetDeltaR", trackSumJetDeltaR, "pseudoangular distance between jet axis and track fourvector sum", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("vertexCategory", vertexCategory, "category of secondary vertex (Reco, Pseudo, No)", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dValAboveCharm", trackSip2dValAboveCharm, "track 2D signed impact parameter of first track lifting mass above charm", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip2dSigAboveCharm", trackSip2dSigAboveCharm, "track 2D signed impact parameter significance of first track lifting mass above charm", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dValAboveCharm", trackSip3dValAboveCharm, "track 3D signed impact parameter of first track lifting mass above charm", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<float>("trackSip3dSigAboveCharm", trackSip3dSigAboveCharm, "track 3D signed impact parameter significance of first track lifting mass above charm", nanoaod::FlatTable::FloatColumn);
    csvTable->addColumn<int>("jetNSelectedTracks", jetNSelectedTracks, "tracks associated to jet", nanoaod::FlatTable::IntColumn);
    csvTable->addColumn<int>("jetNTracksEtaRel", jetNTracksEtaRel, "number of tracks for which etaRel is computed", nanoaod::FlatTable::IntColumn);

    cpfTable->addColumn<int>("jetIdx", cpf_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<float>("trackEtaRel", cpf_trackEtaRel, "track pseudorapidity, relative to the jet axis", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRel", cpf_trackPtRel, "track transverse momentum, relative to the jet axis", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPPar", cpf_trackPPar, "track parallel momentum, along the jet axis", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackDeltaR", cpf_trackDeltaR, "charged PF candidate track pseudoangular distance (Delta R) from the jet axis", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPParRatio", cpf_trackPParRatio, "track parallel momentum, along the jet axis, normalized to its energy", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackPtRatio", cpf_trackPtRatio, "track transverse momentum, relative to the jet axis, normalized to its energy", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dVal", cpf_trackSip2dVal, "track 2D signed impact parameter", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip2dSig", cpf_trackSip2dSig, "track 2D signed impact parameter significance", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dVal", cpf_trackSip3dVal, "track 3D signed impact parameter", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackSip3dSig", cpf_trackSip3dSig, "track 3D signed impact parameter significance", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistVal", cpf_trackJetDistVal, "minimum track approach distance to jet axis", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("trackJetDistSig", cpf_trackJetDistSig, "minimum track approach distance to jet axis significance", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("ptrel", cpf_ptrel, "Charged PF candidate track transverse momentum over jet transverse momentum (uncorrelated)", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("drminsv", cpf_drminsv, "pseudoangular distance (Delta R) between the charged PF candidate track and closest secondary vertex (SV) within the jet", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("vertex_association", cpf_vertex_association, "Indicates whether the charged PF candidate track is used in the primary vertex fit", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("fromPV", cpf_fromPV, "Indicates whether the charged PF candidate track is associated to the primary vertex (PV)", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("puppi_weight", cpf_puppi_weight, "weight assigned by the PUPPI algorithm to neutral PF candidate track", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_chi2", cpf_track_chi2, "chi^2 of a charged PF candidate track", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("track_quality", cpf_track_quality, "Indicates charged particle track reconstruction quality", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<int>("track_ndof", cpf_track_ndof, "Track fit number of degrees of freedom of charged PF candidate track", nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("matchedMuon", cpf_matchedMuon, "flag to specify whether the track is matched to a PF muon", nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("matchedElectron", cpf_matchedElectron, "flag to specify whether the track is matched to a PF electron", nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("matchedSV", cpf_matchedSV, "flag to specify whether the track is matched to a PF secondary vertex", nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("numberOfValidPixelHits", cpf_track_numberOfValidPixelHits , "number of valid pixel hits " , nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("pixelLayersWithMeasurement", cpf_track_pixelLayersWithMeasurement , "pixel layers with measurment ", nanoaod::FlatTable::IntColumn);  
    cpfTable->addColumn<int>("numberOfValidStripHits" , cpf_track_numberOfValidStripHits , "nb of valid strip hits " , nanoaod::FlatTable::IntColumn);
    cpfTable->addColumn<int>("stripLayersWithMeasurement" , cpf_track_stripLayersWithMeasurement , "nb of strip layers with measurement ", nanoaod::FlatTable::IntColumn); 
    cpfTable->addColumn<float>("relmassdrop", cpf_relmassdrop, "doc", nanoaod::FlatTable::FloatColumn);
    cpfTable->addColumn<float>("dzMin", cpf_dZmin, "min distance to other PV", nanoaod::FlatTable::FloatColumn);

    npfTable->addColumn<int>("jetIdx", npf_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);
    npfTable->addColumn<float>("ptrel", npf_ptrel, "neutral PF candidate transverse momentum over jet transverse momentum (uncorrelated)", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("deta", npf_deta, "absolute difference between the neutral PF candidate eta and jet eta", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("dphi", npf_dphi, "absolute difference between the neutral PF candidate phi and jet phi", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("deltaR", npf_deltaR, "neutral PF candidate pseudoangular distance (Delta R) from the jet axis", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("isGamma", npf_isGamma, "flag for passing loose photon identification requirements", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("hcal_fraction", npf_hcal_fraction, "fraction of energy deposited in the hadronic calorimeter", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("drminsv", npf_drminsv, "pseudoangular distance (Delta R) between neutral PF candidate and closest secondary vertex (SV) within the jet", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("puppi_weight", npf_puppi_weight, "weight assigned by the PUPPI algorithm to neutral PF candidate", nanoaod::FlatTable::FloatColumn);
    npfTable->addColumn<float>("relmassdrop", npf_relmassdrop, "neutral PF candidate mass drop normalized relative to the jet", nanoaod::FlatTable::FloatColumn);

    svTable->addColumn<int>("jetIdx", sv_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);
    svTable->addColumn<float>("ptrel", sv_ptrel, "secondary vertex (SV) transverse momentum relative to the uncorrected jet energy pt", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("deta", sv_deta, "absolute difference between the secondary vertex (SV) eta and jet eta", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("dphi", sv_dphi, "absolute difference between the secondary vertex (SV) phi and jet phi", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("deltaR", sv_deltaR, "pseudoangular distance (Delta R) between jet axis and SV flight direction", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("mass", sv_mass, "invariant mass of reconstructed secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("ntracks", sv_ntracks, "number of tracks associated to the secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("chi2", sv_chi2,  "chi^2 of secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<int>("ndof", sv_ndof, "number of degrees of freedom of secondary vertex (SV) fit", nanoaod::FlatTable::IntColumn);
    svTable->addColumn<float>("dxy", sv_dxy, "transverse impact parameter of the secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("dxysig", sv_dxysig, "transverse impact parameter significance of the secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("d3d", sv_d3d, "3D impact parameter of secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("d3dsig", sv_d3dsig, "3D impact parameter significance of secondary vertex (SV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("costhetasvpv", sv_costhetasvpv, "cos(theta) of the secondary vertex (SV) with respect to the primary vertex (PV)", nanoaod::FlatTable::FloatColumn);
    svTable->addColumn<float>("enratio", sv_enratio, "ratio of secondary vertex (SV) energy with respect to the jet", nanoaod::FlatTable::FloatColumn);

    lengthTable->addColumn<int>("cpf", cpf_length, "charged PF candidate track offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("npf", npf_length, "neutral PF candidate offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("sv", sv_length, "secondary vertex (SV) offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("mu", mu_length, "muon offset", nanoaod::FlatTable::IntColumn);
    lengthTable->addColumn<int>("ele", elec_length, "electron offset", nanoaod::FlatTable::IntColumn);

    muonTable->addColumn<int>("jetIdx", mu_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);

    muonTable->addColumn<int>("isGlobal",mu_isGlobal,"muon fitted from the tracker and muon stations",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("isTight",mu_isTight,"global muon with additional muon-quality requirements.Tight Muon ID selects a subset of the Particle-Flow muons",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("isMedium",mu_isMedium,"loose muon (i.e. PF muon that is either a global or an arbitrated tracker muon) with additional track-quality and muon-quality requirements.",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("isLoose",mu_isLoose,"particle identified as a muon by the Particle-Flow event reconstruction, and that is also reconstructed either as a global-muon or as an arbitrated tracker-muon.",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("isStandAlone",mu_isStandAlone,"particle identified as a muon by the Particle-Flow event reconstruction reconstructed only in the muon chambers", nanoaod::FlatTable::IntColumn);

    muonTable->addColumn<float>("ptrel",mu_ptrel,"muon candidate transverse momentum relative to the uncorrected jet energy pt",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("EtaRel",mu_EtaRel,"muon candidate pseudorapidity, relative to the jet axis",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("dphi",mu_dphi,"absolute difference between the muon candidate phi and jet phi",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("deta",mu_deta,"absolute difference between the muon candidate eta and jet eta",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("charge",mu_charge,"muon candidate charge (-1 or +1)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("energy",mu_energy,"muon candidate energy",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("jetDeltaR",mu_jetDeltaR,"pseudoangular distance between jet axis and muon track fourvector",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("numberOfMatchedStations",mu_numberOfMatchedStations,"number of activated muon stations",nanoaod::FlatTable::FloatColumn);

    muonTable->addColumn<float>("2dIP",mu_2dIP,"muon inner track transverse impact parameter relative to the primary vertex in transverse plane (2D), absolute value",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("2dIPSig",mu_2dIPSig,"muon inner track transverse impact parameter relative to the primary vertex relative to its uncertainty in transverse plane (2D)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("3dIP",mu_3dIP,"muon inner track transverse impact parameter relative to the primary vertex in 3D",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("3dIPSig",mu_3dIPSig,"muon inner track transverse impact parameter relative to the primary vertex relative to its uncertainty in 3D",nanoaod::FlatTable::FloatColumn);

// Adding new documentation
    muonTable->addColumn<float>("dxy",mu_dxy,"transverse impact parameter of the best reconstructed muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("dxyError",mu_dxyError,"error on transverse impact parameter of the best reconstructed muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("dxySig",mu_dxySig,"significance of the transverse impact parameter of the best reconstructed muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("dz",mu_dz,"longitudinal impact parameter of the best reconstructed muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("dzError",mu_dzError,"error on the longitudinal impact parameter of the best reconstructed muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<int>("numberOfValidPixelHits",mu_numberOfValidPixelHits,"number of pixel hits by the muon track to further suppress muons from decays in flight",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("numberOfpixelLayersWithMeasurement",mu_numberOfpixelLayersWithMeasurement,"number of pixel layers with valid hits by the muon track",nanoaod::FlatTable::IntColumn);
    muonTable->addColumn<int>("numberOfstripLayersWithMeasurement",mu_numberOfstripLayersWithMeasurement,"doc",nanoaod::FlatTable::IntColumn);

    muonTable->addColumn<float>("chi2",mu_chi2,"chi^2 of muon track",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<int>("ndof",mu_ndof,"number of degrees of freedom of muon track fit",nanoaod::FlatTable::IntColumn);

    muonTable->addColumn<float>("caloIso",mu_caloIso,"returns the fraction of the summed transverse energy of the muon track recHits in the ECAL & HCAL in a cone of deltaR<0.3 normalized against its transverse momentum",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("ecalIso",mu_ecalIso,"returns the fraction of the transverse energy of the muon track recHits in the ECAL normalized against its transverse momentum",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("hcalIso",mu_hcalIso,"returns the fraction of the transverse energy of the muon track recHits in the HCAL normalized against its transverse momentum",nanoaod::FlatTable::FloatColumn);

    muonTable->addColumn<float>("sumPfChHadronPt",mu_sumPfChHadronPt,"summed pt of charged hadron normalized against the muon's pt (PF based isolation of cone of 0.4 - Suggested)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("sumPfNeuHadronEt",mu_sumPfNeuHadronEt,"summed pt of neutral hadron normalized against the muon's pt (PF based isolation of cone of 0.4 - Suggested)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("Pfpileup",mu_Pfpileup,"summed pt of charged particles not from PV (for PU corrections) normalized against the muon's pt (PF based isolation of cone of 0.4 - Suggested)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("sumPfPhotonEt",mu_sumPfPhotonEt,"summed pt of PF photons normalized against the muon's pt (PF based isolation of cone of 0.4 - Suggested)",nanoaod::FlatTable::FloatColumn);

    muonTable->addColumn<float>("sumPfChHadronPt03",mu_sumPfChHadronPt03,"summed pt of charged hadron normalized against the muon's pt (PF based isolation of cone of 0.3)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("sumPfNeuHadronEt03",mu_sumPfNeuHadronEt03,"summed pt of neutral hadron normalized against the muon's pt (PF based isolation of cone of 0.3)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("Pfpileup03",mu_Pfpileup03,"summed pt of charged particles not from PV (for PU corrections) normalized against the muon's pt (PF based isolation of cone of 0.3)",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("sumPfPhotonEt03",mu_sumPfPhotonEt03,"summed pt of PF photons normalized against the muon's pt (PF based isolation of cone of 0.3)",nanoaod::FlatTable::FloatColumn);

    muonTable->addColumn<float>("timeAtIpInOut",mu_timeAtIpInOut,"the time at the interaction point for muons moving inside-out",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("timeAtIpInOutErr",mu_timeAtIpInOutErr,"error on the time at the interaction point for muons moving inside-out",nanoaod::FlatTable::FloatColumn);
    muonTable->addColumn<float>("timeAtIpOutIn",mu_timeAtIpOutIn,"The time at the interaction point for muons moving outside-in",nanoaod::FlatTable::FloatColumn);
    
    electronTable->addColumn<int>("jetIdx", elec_jetIdx, "linked jet Id", nanoaod::FlatTable::IntColumn);
    electronTable->addColumn<float>("ptrel",elec_ptrel,"electron candidate transverse momentum relative to the uncorrected jet energy pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("jetDeltaR",elec_jetDeltaR,"pseudoangular distance between jet axis and electron track fourvector",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deta",elec_deta,"absolute difference between the electron candidate eta and jet eta",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dphi",elec_dphi,"absolute difference between the electron candidate phi and jet phi",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("charge",elec_charge,"electron candidate charge (-1 or +1)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("energy",elec_energy,"electron candidate energy",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("EtFromCaloEn",elec_EtFromCaloEn,"electron transverse energy from calorimeters",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("isEB",elec_isEB,"electron energy deposited in the ECAL barrel",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("isEE",elec_isEE,"electron energy deposited in the ECAL endcap",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("ecalEnergy",elec_ecalEnergy,"electron energy deposited in the ECAL normalized against the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("isPassConversionVeto",elec_isPassConversionVeto,"vertex fit combined with missing number of hits method",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("convDist",elec_convDist,"minimum distance between conversion tracks",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<int>("convFlags",elec_convFlags,"flag for electron coming from photon conversion",nanoaod::FlatTable::IntColumn);

    electronTable->addColumn<float>("convRadius",elec_convRadius,"electron conversion radius",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("hadronicOverEm",elec_hadronicOverEm,"sum of hcal over ecal seed cluster energy using 1st and 2nd hcal depth (using hcal towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("ecalDrivenSeed",elec_ecalDrivenSeed,"ECAL driven seed",nanoaod::FlatTable::FloatColumn);


    electronTable->addColumn<float>("SC_energy",elecSC_energy,"electron supercluster energy normalized against the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("SC_deta",elecSC_deta,"absolute difference between the electron supercluster eta and gsfTrack electron eta",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("SC_dphi",elecSC_dphi,"absolute difference between the electron supercluster phi and gsfTrack electron phi",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("SC_et",elecSC_et,"electron supercluster transverse energy normalized against the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("SC_eSuperClusterOverP",elecSC_eSuperClusterOverP,"the electron supercluster energy / track momentum at the PCA to the beam spot",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scE1x5Overe5x5",elec_scE1x5Overe5x5,"electron supercluster energy ratio of highest energy 1x5 ECAL block containing seed over total energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scE2x5MaxOvere5x5",elec_scE2x5MaxOvere5x5,"electron supercluster energy ratio of highest energy 2x5 ECAL block containing seed over total energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scE5x5",elec_scE5x5,"electron supercluster energy in 5x5 ECAL block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scE5x5Rel",elec_scE5x5Rel,"electron supercluster energy in 5x5 ECAL block on seed relative to the jet pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scPixCharge",elec_scPixCharge,"electron supercluster pixel charge",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scSigmaEtaEta",elec_scSigmaEtaEta,"electron supercluster Sigma Eta Eta variable (quite unintuitive)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("scSigmaIEtaIEta",elec_scSigmaIEtaIEta,"electron supercluster Sigma iEta iEta variable (improved variable against Sigma Eta Eta to resolv issues in cracks)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("superClusterFbrem",elec_superClusterFbrem,"the brem fraction from supercluster: (supercluster energy - electron cluster energy) / supercluster energy",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<float>("2dIP",elec_2dIP,"electron inner track transverse impact parameter relative to the primary vertex in transverse plane (2D), absolute value",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("2dIPSig",elec_2dIPSig,"electron inner track transverse impact parameter relative to the primary vertex relative to its uncertainty in transverse plane (2D)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("3dIP",elec_3dIP,"electron inner track transverse impact parameter relative to the primary vertex in 3D",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("3dIPSig",elec_3dIPSig,"electron inner track transverse impact parameter relative to the primary vertex relative to its uncertainty in 3D",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("eSeedClusterOverP",elec_eSeedClusterOverP,"the seed cluster energy / track momentum at the PCA to the beam spot",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("eSeedClusterOverPout",elec_eSeedClusterOverPout,"the seed cluster energy / track momentum at calo extrapolated from the outermost track state",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("eSuperClusterOverP",elec_eSuperClusterOverP,"the supercluster energy / track momentum at the PCA to the beam spot",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("eTopOvere5x5",elec_eTopOvere5x5,"electron.full5x5_eTop()/ electron.full5x5_e5x5()",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<float>("deltaEtaEleClusterTrackAtCalo",elec_deltaEtaEleClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaEtaSeedClusterTrackAtCalo",elec_deltaEtaSeedClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaPhiSeedClusterTrackAtCalo",elec_deltaPhiSeedClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaEtaSeedClusterTrackAtVtx",elec_deltaEtaSeedClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaEtaSuperClusterTrackAtVtx",elec_deltaEtaSuperClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaPhiEleClusterTrackAtCalo",elec_deltaPhiEleClusterTrackAtCalo,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("deltaPhiSuperClusterTrackAtVtx",elec_deltaPhiSuperClusterTrackAtVtx,"doc",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("sCseedEta",elec_sCseedEta,"doc",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<float>("EtaRel",elec_EtaRel,"electron pseudorapidity, relative to the jet axis",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dxy",elec_dxy,"transverse impact parameter of the electron",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dz",elec_dz,"longitudinal impact parameter of the best reconstructed electron track",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("nbOfMissingHits",elec_nbOfMissingHits,"number of missing electron hits in its hit pattern",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("gsfCharge",elec_gsfCharge,"gsf electron charge",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<int>("ndof",elec_ndof,"gsf electron number degree of freedom",nanoaod::FlatTable::IntColumn);
    electronTable->addColumn<float>("chi2",elec_chi2,"chi2 of the fit ",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<int>("numberOfBrems",elec_numberOfBrems,"number of Bremsstrahlung electrons in a shower",nanoaod::FlatTable::IntColumn);
    electronTable->addColumn<float>("fbrem",elec_fbrem,"doc",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<float>("e5x5",elec_e5x5,"electron energy in 5x5 ECAL block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("e5x5Rel",elec_e5x5Rel,"electron energy in 5x5 ECAL block on seed relative to the jet pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("e2x5MaxOvere5x5",elec_e2x5MaxOvere5x5,"electron energy ratio of highest energy 2x5 ECAL block containing seed over total energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("e1x5Overe5x5",elec_e1x5Overe5x5,"electron energy ratio of highest energy 1x5 ECAL block containing seed over total energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_sigmaIetaIeta",elec_full5x5_sigmaIetaIeta,"electron full shower Sigma iEta iEta variable (improved variable against Sigma Eta Eta to resolv issues in cracks)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e1x5Overe5x5",elec_full5x5_e1x5Overe5x5,"electron full shower energy ratio of highest energy 1x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e2x5BottomOvere5x5",elec_full5x5_e2x5BottomOvere5x5,"electron full shower energy ratio of bottom 2x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e2x5LeftOvere5x5",elec_full5x5_e2x5LeftOvere5x5,"electron full shower energy ratio of left 2x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e2x5MaxOvere5x5",elec_full5x5_e2x5MaxOvere5x5,"electron full shower energy ratio of highest energy 2x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e2x5RightOvere5x5",elec_full5x5_e2x5RightOvere5x5,"electron full shower energy ratio of right 2x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e2x5TopOvere5x5",elec_full5x5_e2x5TopOvere5x5,"electron full shower energy ratio of top 2x5 ECAL block containing seed over total full shower energy in 5x5 block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e5x5",elec_full5x5_e5x5,"electron full shower energy in 5x5 full shower ECAL block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_e5x5Rel",elec_full5x5_e5x5Rel,"electron full shower energy in 5x5 full shower ECAL block on seed relative to the jet pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_eBottomOvere5x5",elec_full5x5_eBottomOvere5x5,"electron full shower energy in bottom ECAL blocks ratio of 5x5 full shower block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_eLeftOvere5x5",elec_full5x5_eLeftOvere5x5,"electron full shower energy in left ECAL blocks ratio of 5x5 full shower block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_eRightOvere5x5",elec_full5x5_eRightOvere5x5,"electron full shower energy in right ECAL blocks ratio of 5x5 full shower block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_eTopOvere5x5",elec_full5x5_eTopOvere5x5,"electron full shower energy in top ECAL blocks ratio of 5x5 full shower block on seed",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalDepth1OverEcal",elec_full5x5_hcalDepth1OverEcal,"HCAL over ECAL seed cluster full shower energy using 1st HCAL depth (using HCAL towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalDepth1OverEcalBc",elec_full5x5_hcalDepth1OverEcalBc,"HCAL over ECAL seed cluster full shower energy using 1st HCAL depth (using HCAL towers behind clusters)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalDepth2OverEcal",elec_full5x5_hcalDepth2OverEcal,"HCAL over ECAL seed cluster full shower energy using 2nd HCAL depth (using HCAL towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalDepth2OverEcalBc",elec_full5x5_hcalDepth2OverEcalBc,"HCAL over ECAL seed cluster full shower energy using 2nd HCAL depth (using HCAL towers behind clusters)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalOverEcal",elec_full5x5_hcalOverEcal,"HCAL over ECAL seed cluster full shower energy summing the 1st and 2nd HCAL depth (using HCAL towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_hcalOverEcalBc",elec_full5x5_hcalOverEcalBc,"HCAL over ECAL seed cluster full shower energy summing the 1st and 2nd HCAL depth (using HCAL towers behind clusters)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("full5x5_r9",elec_full5x5_r9,"ratio of the 3x3 full shower energy and supercluster energy",nanoaod::FlatTable::FloatColumn);

    electronTable->addColumn<float>("neutralHadronIso",elec_neutralHadronIso,"returns the isolation calculated with only the neutral hadron PFCandidates relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("particleIso",elec_particleIso,"returns the isolation calculated with all the PFCandidates relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("photonIso",elec_photonIso,"returns the isolation calculated with only the gamma PFCandidates relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("puChargedHadronIso",elec_puChargedHadronIso,"returns the isolation calculated with only the pile-up charged hadron PFCandidates relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("trackIso",elec_trackIso,"returns the tracker isolation variable that was stored in this object when produced, or -1.0 if there is none, relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("hcalDepth1OverEcal",elec_hcalDepth1OverEcal,"HCAL over ECAL seed cluster energy using 1st HCAL depth (using HCAL towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("hcalDepth2OverEcal",elec_hcalDepth2OverEcal,"HCAL over ECAL seed cluster energy using 2nd HCAL depth (using HCAL towers within a cone)",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("ecalPFClusterIso",elec_ecalPFClusterIso,"sum pt of isolated ECAL clusters, vetoing clusters part of electron",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("hcalPFClusterIso",elec_hcalPFClusterIso,"sum pt of isolated HCAL clusters, vetoing clusters part of electron",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03TkSumPt",elec_dr03TkSumPt,"track iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03EcalRecHitSumEt",elec_dr03EcalRecHitSumEt,"ECAL iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03HcalDepth1TowerSumEt",elec_dr03HcalDepth1TowerSumEt,"HCAL depth 1 iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03HcalDepth1TowerSumEtBc",elec_dr03HcalDepth1TowerSumEtBc,"HCAL depth 1 iso deposit without towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03HcalDepth2TowerSumEt",elec_dr03HcalDepth2TowerSumEt,"HCAL depth 2 iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr03HcalDepth2TowerSumEtBc",elec_dr03HcalDepth2TowerSumEtBc,"HCAL depth 2 iso deposit without towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("pfSumPhotonEt",elec_pfSumPhotonEt,"sum pt of PF photons // old float photonIso, relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("pfSumChargedHadronPt",elec_pfSumChargedHadronPt,"sum pt of charged hadron // old float chargedHadronIso, relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("pfSumNeutralHadronEt",elec_pfSumNeutralHadronEt,"sum transverse energy of neutral hadron // old float neutralHadronIso, relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("pfSumPUPt",elec_pfSumPUPt,"sum pt of charged particles not from PV (for PU corrections) ",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04EcalRecHitSumEt",elec_dr04EcalRecHitSumEt,"ECAL iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalDepth1TowerSumEt",elec_dr04HcalDepth1TowerSumEt,"HCAL depth 1 iso deposit with electron footprint removed relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalDepth1TowerSumEtBc",elec_dr04HcalDepth1TowerSumEtBc,"HCAL depth 1 iso deposit without towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalDepth2TowerSumEt",elec_dr04HcalDepth2TowerSumEt,"HCAL depth 2 iso deposit with towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalDepth2TowerSumEtBc",elec_dr04HcalDepth2TowerSumEtBc,"HCAL depth 2 iso deposit without towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalTowerSumEt",elec_dr04HcalTowerSumEt,"HCAL depth 1+2 iso deposit with towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);
    electronTable->addColumn<float>("dr04HcalTowerSumEtBc",elec_dr04HcalTowerSumEtBc,"HCAL depth 1+2 iso deposit without towers behind cluster relative to the electron pt",nanoaod::FlatTable::FloatColumn);

    iEvent.put(std::move(globalTable), "global");
    iEvent.put(std::move(csvTable), "csv");
    iEvent.put(std::move(cpfTable), "cpf");
    iEvent.put(std::move(npfTable), "npf");
    iEvent.put(std::move(svTable), "sv");
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
