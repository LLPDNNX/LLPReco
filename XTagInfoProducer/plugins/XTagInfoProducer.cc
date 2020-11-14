//Description: Produces and fill in LLPDNNX features

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "LLPReco/DataFormats/interface/XTagInfo.h"

#include "RecoBTag/FeatureTools/interface/JetConverter.h"
#include "RecoBTag/FeatureTools/interface/ShallowTagInfoConverter.h"
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"
#include "RecoBTag/FeatureTools/interface/NeutralCandidateConverter.h"
#include "RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "LLPReco/XTagInfoProducer/interface/JetSubstructure.h"

#include "TVector3.h"


llpdnnx::SecondaryVertexFeatures fillSVFeatures(const reco::VertexCompositePtrCandidate& sv, const reco::Vertex& pv, const pat::Jet& jet){
    float uncorrectedPt = jet.correctedP4("Uncorrected").pt();
    llpdnnx::SecondaryVertexFeatures sv_features;

    sv_features.ptrel = sv.pt()/uncorrectedPt;
    sv_features.deta = sv.eta()-jet.eta();
    sv_features.dphi = reco::deltaPhi(sv.phi(),jet.phi());
    sv_features.deltaR = reco::deltaR(sv,jet);
    sv_features.mass = sv.mass();
    sv_features.ntracks = sv.numberOfDaughters();
    sv_features.chi2 = sv.vertexChi2();
    sv_features.ndof = sv.vertexNdof();


    reco::Vertex::CovarianceMatrix covsv; 
    sv.fillVertexCovariance(covsv);
    reco::Vertex svtx(sv.vertex(), covsv);

    VertexDistanceXY distXY;
    Measurement1D distanceXY = distXY.distance(svtx, pv);
    sv_features.dxy = distanceXY.value();
    sv_features.dxysig = distanceXY.value()/distanceXY.error();

    VertexDistance3D dist3D;
    Measurement1D distance3D = dist3D.distance(svtx, pv);
    sv_features.d3d = distance3D.value();
    sv_features.d3dsig = distance3D.value()/distance3D.error();

    if (std::isnan(sv_features.dxysig) || std::isnan(sv_features.d3dsig))
    {
        sv_features.dxysig = 0.;
        sv_features.d3dsig = 0.;
    }

    reco::Candidate::Vector distance(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    sv_features.costhetasvpv = sv.momentum().Unit().Dot(distance.Unit());
    sv_features.enratio = sv.energy()/jet.pt();

    return sv_features;
}


class XTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit XTagInfoProducer(const edm::ParameterSet&);
    ~XTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    struct CandidateHash
    {
        long operator() (const reco::CandidatePtr& cand) const 
        {
            return cand.id().id() * 100000 + cand.key();
        }
    };
    
    private:
        virtual void beginStream(edm::StreamID) override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endStream() override;


        edm::EDGetTokenT<edm::View<pat::Jet>> jet_token_;
        edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;
        edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_adapted_token_;
        edm::EDGetTokenT<edm::View<reco::ShallowTagInfo>> shallow_tag_info_token_;
        edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
        typedef std::vector<reco::XTagInfo> XTagInfoCollection;

        edm::EDGetTokenT< pat::MuonCollection > muonsMiniAODToken_;
        edm::EDGetTokenT< pat::ElectronCollection > electronsMiniAODToken_;
};

XTagInfoProducer::XTagInfoProducer(const edm::ParameterSet& iConfig) :
    jet_token_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
    sv_adapted_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices_adapted"))),
    shallow_tag_info_token_(consumes<edm::View<reco::ShallowTagInfo>>(iConfig.getParameter<edm::InputTag>("shallow_tag_infos"))),
    muonsMiniAODToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    electronsMiniAODToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc")))
{
    produces<XTagInfoCollection>();
}


XTagInfoProducer::~XTagInfoProducer(){ }
void XTagInfoProducer::beginStream(edm::StreamID) { }
// ------------ method called to produce the data  ------------
    void
XTagInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    auto output_tag_infos = std::make_unique<XTagInfoCollection>();
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jet_token_, jets);

    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByToken(vtx_token_, vtxs);

    if (vtxs->empty()) {
        // produce empty TagInfos in case no primary vertex
        iEvent.put(std::move(output_tag_infos));
        return;  // exit event
    }
    const auto& pv = vtxs->at(0);
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

    edm::Handle<edm::View<reco::ShallowTagInfo>> shallow_tag_infos;
    iEvent.getByToken(shallow_tag_info_token_, shallow_tag_infos);

    edm::Handle<reco::VertexCompositePtrCandidateCollection> svs;
    iEvent.getByToken(sv_token_, svs);

    edm::Handle<reco::VertexCompositePtrCandidateCollection> svs_adapted;
    iEvent.getByToken(sv_adapted_token_, svs_adapted);

    edm::Handle< pat::MuonCollection > muons;
    iEvent.getByToken(muonsMiniAODToken_, muons);

    edm::Handle<pat::ElectronCollection > electrons;
    iEvent.getByToken(electronsMiniAODToken_, electrons);

    std::unordered_map<reco::CandidatePtr, const pat::Muon*, CandidateHash> muonMap;
    std::unordered_map<reco::CandidatePtr, const pat::Electron*, CandidateHash> electronMap;

    for (const pat::Muon& muon: *muons)
    {
        for (unsigned int i = 0 ; i < muon.numberOfSourceCandidatePtrs(); ++i )
        {
            muonMap[muon.sourceCandidatePtr(i)] = &muon;
        }
    }


    for (const pat::Electron& electron: *electrons)
    {
        for (unsigned int i = 0 ; i < electron.numberOfSourceCandidatePtrs(); ++i )
        {
            electronMap[electron.sourceCandidatePtr(i)] = &electron;
        }
    }

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) 
    {
        const pat::Jet& jet = jets->at(ijet);

        edm::RefToBase<reco::Jet> jet_ref(jets->refAt(ijet)); //upcast

        std::unordered_set<reco::CandidatePtr, CandidateHash> jetConstituentSet;

        bool badJet = 0;
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            jetConstituentSet.insert(jet.daughterPtr(idaughter));
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if (std::isinf(constituent->pt())){
                edm::LogWarning("BadTransverseMomentum") << "dropping jet with input candidate with inf";
                badJet = 1;
            }
        }

        // Cut on eta
        if (std::abs(jet.eta()) > 2.5 or badJet) continue;

        /*
        float NHF  = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF  = jet.chargedHadronEnergyFraction();
        float MUF  = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        int NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        int NumNeutralParticles =jet.neutralMultiplicity();
        float CHM      = jet.chargedMultiplicity();

        int looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jet.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet.eta())>2.4) && abs(jet.eta())<=2.7;
        */

        // create data containing structure
        llpdnnx::XTagFeatures features;
        

        // Start with global jet features
        float uncorrectedPt = jet.correctedP4("Uncorrected").pt();
        
        features.jet_features.pt = uncorrectedPt;  // uncorrected
        features.jet_features.eta = jet.eta();

        features.jet_features.phi = jet.phi();
        
        features.jet_features.energy = jet.energy();
        features.jet_features.area = jet.jetArea();
        
        features.jet_features.n60 = jet.n60();
        features.jet_features.n90 = jet.n90();
        
        features.jet_features.chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
        features.jet_features.chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
        features.jet_features.chargedMuEnergyFraction = jet.chargedMuEnergyFraction();
        features.jet_features.electronEnergyFraction = jet.electronEnergyFraction();

        features.npv = vtxs->size();
        
        llpdnnx::JetSubstructure jetSubstructure(jet);
        
        //mass calculated from constituents seems to be different from the stored jet mass (likely due to reduced numerical precision)
        features.jet_features.mass = jetSubstructure.massFromConstituents();

        features.jet_features.relMassDropMassAK = jetSubstructure.relMassDropMass(llpdnnx::JetSubstructure::ClusterType::AK);
        features.jet_features.relMassDropMassCA = jetSubstructure.relMassDropMass(llpdnnx::JetSubstructure::ClusterType::CA);
        features.jet_features.relSoftDropMassAK = jetSubstructure.relSoftDropMass(llpdnnx::JetSubstructure::ClusterType::AK);
        features.jet_features.relSoftDropMassCA = jetSubstructure.relSoftDropMass(llpdnnx::JetSubstructure::ClusterType::CA);

        // Still need to fix bug!

        features.jet_features.tau1 = jetSubstructure.nSubjettiness(1);
        features.jet_features.tau2 = jetSubstructure.nSubjettiness(2);
        features.jet_features.tau3 = jetSubstructure.nSubjettiness(3);
   
        if (jetSubstructure.nConstituents()>3)
        {
            auto eventShapes = jetSubstructure.eventShapeVariables();
            features.jet_features.thrust = jetSubstructure.thrust();
            features.jet_features.sphericity = eventShapes.sphericity();
            features.jet_features.circularity = eventShapes.circularity();
            features.jet_features.isotropy = eventShapes.isotropy();
            features.jet_features.eventShapeC = eventShapes.C();
            features.jet_features.eventShapeD = eventShapes.D();
        }
        
       
        // Add CSV variables
        const edm::View<reco::ShallowTagInfo>& taginfos = *shallow_tag_infos;
        edm::Ptr<reco::ShallowTagInfo> match;


        for (auto it = taginfos.begin(); it != taginfos.end(); ++it) {
            float dR = reco::deltaR(it->jet()->p4(),jet.p4());
            if (dR<0.01) {
                match = taginfos.ptrAt(it - taginfos.begin());
                break;
            }
        }
        reco::ShallowTagInfo tag_info;
        if (match.isNonnull()) {
            tag_info = *match;
        }  // will be default values otherwise

        reco::TaggingVariableList vars = tag_info.taggingVariables();
        features.tag_info_features.trackSumJetEtRatio = vars.get(reco::btau::trackSumJetEtRatio, -1);
        features.tag_info_features.trackSumJetDeltaR = vars.get(reco::btau::trackSumJetDeltaR, -1);
        features.tag_info_features.vertexCategory = vars.get(reco::btau::vertexCategory, -1);
        features.tag_info_features.trackSip2dValAboveCharm = vars.get(reco::btau::trackSip2dValAboveCharm, -1);
        features.tag_info_features.trackSip2dSigAboveCharm = vars.get(reco::btau::trackSip2dSigAboveCharm, -1);
        features.tag_info_features.trackSip3dValAboveCharm = vars.get(reco::btau::trackSip3dValAboveCharm, -1);
        features.tag_info_features.trackSip3dSigAboveCharm = vars.get(reco::btau::trackSip3dSigAboveCharm, -1);
        features.tag_info_features.jetNTracksEtaRel = vars.get(reco::btau::jetNTracksEtaRel, -1);
        features.tag_info_features.jetNSelectedTracks = vars.get(reco::btau::jetNSelectedTracks, -1);


        std::unordered_set<reco::CandidatePtr, CandidateHash> candidatesMatchedToSV;
        std::unordered_set<reco::CandidatePtr, CandidateHash> candidatesMatchedToSVAdapted;

        // fill features from secondary vertices  
        
        for (unsigned int isv = 0; isv < svs->size(); ++isv)
        {
            const reco::VertexCompositePtrCandidate& sv = svs->at(isv);
            
            if (reco::deltaR(sv,jet)>0.4)
            {
                continue;
            }
            bool matchingTrack = false;
            for (auto const& candidateFromVertex: sv.daughterPtrVector())
            {
                if (jetConstituentSet.find(candidateFromVertex)!=jetConstituentSet.end())
                {
                    candidatesMatchedToSV.insert(candidateFromVertex);
                    matchingTrack = true;
                }
            }
            if (not matchingTrack) continue;



            llpdnnx::SecondaryVertexFeatures sv_features = fillSVFeatures(sv, pv, jet);
            features.sv_features.emplace_back(sv_features);
        }

        std::stable_sort(features.sv_features.begin(),features.sv_features.end());


        // fill features from adapted secondary vertices  
        
        for (unsigned int isv = 0; isv < svs_adapted->size(); ++isv)
        {
            const reco::VertexCompositePtrCandidate& sv_adapted = svs_adapted->at(isv);
            
            if (reco::deltaR(sv_adapted,jet)>0.4)
            {
                continue;
            }
            bool matchingTrack = false;
            for (auto const& candidateFromVertex: sv_adapted.daughterPtrVector())
            {
                if (jetConstituentSet.find(candidateFromVertex)!=jetConstituentSet.end())
                {
                    candidatesMatchedToSVAdapted.insert(candidateFromVertex);
                    matchingTrack = true;
                }
            }
            if (not matchingTrack) continue;

            llpdnnx::SecondaryVertexFeatures sv_adapted_features = fillSVFeatures(sv_adapted, pv, jet);
            features.sv_adapted_features.emplace_back(sv_adapted_features);
        }

        std::stable_sort(features.sv_adapted_features.begin(),features.sv_adapted_features.end());


        // Fill cpf info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()==0 or (not constituent->hasTrackDetails()))
            {
                continue;
            }

            if (constituent->pt() < 1e-10){
                edm::LogWarning("NullTransverseMomentum") << "dropping input candidate with pt<1e-10";
                continue;
            }

            llpdnnx::ChargedCandidateFeatures cpf_features;

            cpf_features.ptrel = constituent->pt()/uncorrectedPt;
            cpf_features.deta = constituent->eta()-jet.eta();
            cpf_features.dphi = reco::deltaPhi(constituent->phi(),jet.phi());
            
            cpf_features.px = constituent->px();
            cpf_features.py = constituent->py();
            cpf_features.pz = constituent->pz();

            cpf_features.drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                cpf_features.drminsv = std::min(cpf_features.drminsv,dR);
            }


            float dZmin = 100;
            for (size_t i = 0; i < vtxs->size(); i++){
                auto vtx = vtxs->at(i);
                if (vtx.isFake() || vtx.ndof() < 4) {
                    continue;
                }
                if ((vtx.position()-pv.position()).mag2()<1e-3) continue; //skip PV
                dZmin = std::min(dZmin, std::abs(constituent->dz(vtx.position())));
            }

            cpf_features.dZmin = dZmin;
            cpf_features.vertex_association = constituent->pvAssociationQuality();
            cpf_features.fromPV = constituent->fromPV();
            cpf_features.puppi_weight = constituent->puppiWeight();
            cpf_features.track_chi2 = constituent->pseudoTrack().chi2();
            cpf_features.track_ndof = constituent->pseudoTrack().ndof();
            cpf_features.track_quality = constituent->pseudoTrack().qualityMask();
    	    cpf_features.track_numberOfValidPixelHits = constituent->pseudoTrack().hitPattern().numberOfValidPixelHits();
    	    cpf_features.track_pixelLayersWithMeasurement  = constituent->pseudoTrack().hitPattern().pixelLayersWithMeasurement();
    	    cpf_features.track_numberOfValidStripHits = constituent->pseudoTrack().hitPattern().numberOfValidStripHits();
    	    cpf_features.track_stripLayersWithMeasurement = constituent->pseudoTrack().hitPattern().stripLayersWithMeasurement();
		

            if (jet.mass()<1e-10)
            {
                cpf_features.relmassdrop = -1;
            }
            else
            {
                cpf_features.relmassdrop = (jet.p4()-constituent->p4()).mass()/jet.mass();
            }
            
            reco::TransientTrack transientTrack = builder->build(constituent->pseudoTrack());
            reco::Candidate::Vector jetDir = jet.momentum().Unit();
            GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());

            Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, jetRefTrackDir, pv).second;
            Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, jetRefTrackDir, pv).second;
            reco::Candidate::Vector trackMom = constituent->pseudoTrack().momentum();
            double trackMag = std::sqrt(trackMom.Mag2());
            TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
            TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());

            cpf_features.trackEtaRel=reco::btau::etaRel(jetDir, trackMom);
            cpf_features.trackPtRel=trackMom3.Perp(jetDir3);
            cpf_features.trackPPar=jetDir.Dot(trackMom);
            cpf_features.trackDeltaR=reco::deltaR(trackMom, jetDir);
            cpf_features.trackPtRatio=cpf_features.trackPtRel / trackMag;
            cpf_features.trackPParRatio=cpf_features.trackPPar / trackMag;

            cpf_features.trackSip2dVal=meas_ip2d.value();
            cpf_features.trackSip2dSig=meas_ip2d.significance();
            cpf_features.trackSip3dVal=meas_ip3d.value();
            cpf_features.trackSip3dSig=meas_ip3d.significance();
            if (std::isnan(cpf_features.trackSip2dSig) || std::isnan(cpf_features.trackSip3dSig))
            {
                cpf_features.trackSip2dSig=0.;
                cpf_features.trackSip3dSig=0.;
            }

            cpf_features.trackJetDistVal = jetdist.value();
            cpf_features.trackJetDistSig = jetdist.significance();

            cpf_features.matchedMuon = 0;
            cpf_features.matchedElectron = 0;
            
            if (candidatesMatchedToSV.find(jet.daughterPtr(idaughter))!=candidatesMatchedToSV.end())
            {
                cpf_features.matchedSV = 1;
            }
            else
            {
                cpf_features.matchedSV = 0;
            }

            if (candidatesMatchedToSVAdapted.find(jet.daughterPtr(idaughter))!=candidatesMatchedToSVAdapted.end())
            {
                cpf_features.matchedSV_adapted = 1;
            }
            else
            {
                cpf_features.matchedSV_adapted = 0;
            }

            
            //find matching muons
            auto findMuon = muonMap.find(jet.daughterPtr(idaughter));  
            if (findMuon!=muonMap.end())
            {
                llpdnnx::MuonCandidateFeatures mu_features; 
                const pat::Muon & muon = *findMuon->second;

                if (not muon.isGlobalMuon() || reco::deltaR(muon, jet) > 0.4) continue;
                cpf_features.matchedMuon = 1;
                mu_features.isGlobal = muon.isGlobalMuon();                                   
                mu_features.isTight = muon.isTightMuon(pv);                                     
                mu_features.isMedium = muon.isMediumMuon();
                mu_features.isLoose = muon.isLooseMuon();
                mu_features.isStandAlone = muon.isStandAloneMuon();

                mu_features.ptrel = muon.pt()/uncorrectedPt;
                mu_features.deta = muon.eta()-jet.eta();                                      
                mu_features.dphi = reco::deltaPhi(muon.phi(),jet.phi());     
                mu_features.px = muon.px();
                mu_features.py = muon.py();
                mu_features.pz = muon.pz();                                          
                mu_features.charge = muon.charge();
                mu_features.energy = muon.energy()/muon.pt();                                   
                mu_features.et = muon.et();
                mu_features.deltaR = reco::deltaR(muon, jet);
                mu_features.numberOfMatchedStations = muon.numberOfMatchedStations();

                mu_features.IP2d = muon.dB();
                mu_features.IP2dSig = muon.dB()/muon.edB();
                mu_features.IP3d = muon.dB(pat::Muon::PV3D);
                mu_features.IP3dSig = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D);

                if (std::isnan(mu_features.IP2dSig) || std::isnan(mu_features.IP3dSig))
                {
                    mu_features.IP2dSig = 0.;
                    mu_features.IP3dSig = 0.;
                }


                reco::Candidate::Vector muonMom = muon.bestTrack()->momentum();

                mu_features.EtaRel =reco::btau::etaRel(jetDir, muonMom);
                mu_features.dxy = muon.bestTrack()->dxy(pv.position());
                mu_features.dxyError = muon.bestTrack()->dxyError();
                mu_features.dxySig = muon.bestTrack()->dxy(pv.position())/(1e-10+std::fabs(muon.bestTrack()->dxyError())); 
                mu_features.dz = muon.bestTrack()->dz(pv.position());
                mu_features.dzError = muon.bestTrack()->dzError();
                mu_features.dzSig = muon.bestTrack()->dz(pv.position())/(1e-10+std::fabs(muon.bestTrack()->dzError()));
                mu_features.numberOfValidPixelHits = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
                mu_features.numberOfpixelLayersWithMeasurement = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement();
                mu_features.numberOfstripLayersWithMeasurement = muon.bestTrack()->hitPattern().stripLayersWithMeasurement();
	

                mu_features.chi2 = muon.bestTrack()->chi2();
                mu_features.ndof = muon.bestTrack()->ndof();

                mu_features.caloIso =  muon.caloIso()/muon.pt();
                mu_features.ecalIso =  muon.ecalIso()/muon.pt(); 
                mu_features.hcalIso =  muon.hcalIso()/muon.pt();     


                mu_features.sumPfChHadronPt  = muon.pfIsolationR04().sumChargedHadronPt/muon.pt();
                mu_features.sumPfNeuHadronEt  = muon.pfIsolationR04().sumNeutralHadronEt/muon.pt();
                mu_features.Pfpileup  = muon.pfIsolationR04().sumPUPt/muon.pt();
                mu_features.sumPfPhotonEt = muon.pfIsolationR04().sumPhotonEt/muon.pt();


                mu_features.sumPfChHadronPt03  = muon.pfIsolationR03().sumChargedHadronPt/muon.pt();
                mu_features.sumPfNeuHadronEt03  = muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
                mu_features.Pfpileup03  = muon.pfIsolationR03().sumPUPt/muon.pt();
                mu_features.sumPfPhotonEt03 = muon.pfIsolationR03().sumPhotonEt/muon.pt();       


                mu_features.timeAtIpInOut = muon.time().timeAtIpInOut;
                mu_features.timeAtIpInOutErr = muon.time().timeAtIpInOutErr;
                mu_features.timeAtIpOutIn = muon.time().timeAtIpOutIn; 

                features.mu_features.emplace_back(mu_features);
            }

            std::stable_sort(features.mu_features.begin(),features.mu_features.end());


            //find matching electrons
            auto findElectron = electronMap.find(jet.daughterPtr(idaughter));  
            if(findElectron!=electronMap.end())
            {
                llpdnnx::ElectronCandidateFeatures elec_features;
                const pat::Electron & electron = *findElectron->second;
                cpf_features.matchedElectron = 1;
		        if (reco::deltaR(electron, jet) > 0.4) continue; 

                elec_features.ptrel = electron.pt()/uncorrectedPt;
                elec_features.deta = electron.eta()-jet.eta();
                elec_features.dphi = reco::deltaPhi(electron.phi(),jet.phi()); 
                elec_features.charge = electron.charge(); 
                elec_features.px = electron.px();
                elec_features.py = electron.py();
                elec_features.pz = electron.pz();
                
                elec_features.energy = electron.energy()/electron.pt(); 
                elec_features.deltaR = reco::deltaR(electron,jet); 
                elec_features.EtFromCaloEn = electron.caloEnergy() * sin(electron.p4().theta())/ electron.pt();
                elec_features.ecalDrivenSeed = electron.ecalDrivenSeed();

                elec_features.isEB = electron.isEB();
                elec_features.isEE  = electron.isEE();
                elec_features.ecalEnergy  = electron.ecalEnergy()/electron.pt();
                elec_features.isPassConversionVeto = electron.passConversionVeto();
        		if(electron.convDist() >= 0.)
        		{
                    elec_features.convDist = electron.convDist(); 
                    elec_features.convFlags = electron.convFlags(); 
                    elec_features.convRadius = electron.convRadius();
         		}
        		else 
        		{
                    elec_features.convDist = -1.; 
                    elec_features.convFlags = -1.; 
                    elec_features.convRadius = -1.;
        		}


                elec_features.IP3d = electron.dB(pat::Electron::PV3D); 
                elec_features.IP3dSig = electron.dB(pat::Electron::PV3D)/electron.edB(pat::Electron::PV3D); 
                elec_features.IP2d = electron.dB();
                elec_features.IP2dSig = electron.dB()/electron.edB();

                if (std::isnan(elec_features.IP2dSig) || std::isnan(elec_features.IP3dSig))
                {
                    elec_features.IP2dSig = 0.;
                    elec_features.IP3dSig = 0.;
                }

                elec_features.sCseedEta = electron.superCluster()->seed()->eta();


                elec_features.eSeedClusterOverP = electron.eSeedClusterOverP();
                elec_features.eSeedClusterOverPout = electron.eSeedClusterOverPout();
                elec_features.eSuperClusterOverP = electron.eSuperClusterOverP();
                elec_features.hadronicOverEm = electron.hadronicOverEm();


                elec_features.deltaEtaEleClusterTrackAtCalo  = electron.deltaEtaEleClusterTrackAtCalo();
                elec_features.deltaPhiEleClusterTrackAtCalo = electron.deltaPhiEleClusterTrackAtCalo();

                elec_features.deltaEtaSeedClusterTrackAtCalo = electron.deltaEtaSeedClusterTrackAtCalo(); 
                elec_features.deltaPhiSeedClusterTrackAtCalo = electron.deltaPhiSeedClusterTrackAtCalo(); 

                elec_features.deltaEtaSeedClusterTrackAtVtx = electron.deltaEtaSeedClusterTrackAtVtx(); 
                elec_features.deltaEtaSuperClusterTrackAtVtx = electron.deltaEtaSuperClusterTrackAtVtx();  
                elec_features.deltaPhiSuperClusterTrackAtVtx = electron.deltaPhiSuperClusterTrackAtVtx();


                reco::Candidate::Vector electronMom = electron.gsfTrack()->momentum();

                elec_features.EtaRel = reco::btau::etaRel(jetDir, electronMom); 
                elec_features.dxy = electron.gsfTrack()->dxy(pv.position());
                elec_features.dxyError = electron.gsfTrack()->dxyError();
                elec_features.dxySig = electron.gsfTrack()->dxy(pv.position())/(1e-10+std::fabs(electron.gsfTrack()->dxyError()));
                elec_features.dz = electron.gsfTrack()->dz(pv.position());
                elec_features.dzError = electron.gsfTrack()->dzError();
                elec_features.dzSig = electron.gsfTrack()->dz(pv.position())/(1e-10+std::fabs(electron.gsfTrack()->dzError()));
                
                elec_features.nbOfMissingHits = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
                elec_features.gsfCharge = electron.gsfTrack()->charge();
        		elec_features.ndof = electron.gsfTrack()->ndof(); 
        		elec_features.chi2 = electron.gsfTrack()->chi2();


                elec_features.elecSC_energy = electron.superCluster()->energy()/electron.pt(); 
                elec_features.elecSC_deta = electron.superCluster()->eta()-electron.gsfTrack()->eta();
                elec_features.elecSC_dphi = reco::deltaPhi(electron.superCluster()->phi(),electron.gsfTrack()->phi());
                elec_features.elecSC_et = electron.superCluster()->energy() * sin(electron.p4().theta())/electron.pt();
                elec_features.scPixCharge = electron.scPixCharge();


                elec_features.numberOfBrems  = electron.numberOfBrems();
		        if(electron.pt() >= 5.){
                    elec_features.fbrem = electron.fbrem();
                    elec_features.sigmaEtaEta = electron.sigmaEtaEta();
                    elec_features.sigmaIetaIeta = electron.sigmaIetaIeta();
                    elec_features.sigmaIphiIphi = electron.sigmaIphiIphi();
                    elec_features.r9 = electron.r9();
                    elec_features.superClusterFbrem = electron.superClusterFbrem();

		        }
        		else 
        		{
        	        elec_features.fbrem = -1.;
        			elec_features.sigmaEtaEta = -1.;
        			elec_features.sigmaIetaIeta = -1.;
                    elec_features.sigmaIphiIphi = -1;
        			elec_features.superClusterFbrem = -1.;
        		}
                elec_features.e5x5 = electron.e5x5();

                elec_features.e5x5Rel = electron.e5x5()/jet.pt();
                elec_features.e1x5Overe5x5 = electron.e1x5()/electron.e5x5();
                elec_features.e2x5MaxOvere5x5 = electron.e2x5Max()/electron.e5x5();

                if (electron.e5x5() == 0){
                    elec_features.e1x5Overe5x5 = -1.;
                    elec_features.e2x5MaxOvere5x5 = -1.;
                }
                elec_features.hcalOverEcal = electron.hcalOverEcal();
                elec_features.hcalDepth1OverEcal = electron.hcalDepth1OverEcal();
                elec_features.hcalDepth2OverEcal = electron.hcalDepth2OverEcal();
	
                elec_features.elecSC_eSuperClusterOverP  = electron.eSuperClusterOverP();
 
                elec_features.neutralHadronIso  = electron.neutralHadronIso()/electron.pt();
                elec_features.photonIso = electron.photonIso()/electron.pt(); 
                elec_features.puChargedHadronIso = electron.puChargedHadronIso()/electron.pt(); 

                elec_features.trackIso = electron.trackIso()/electron.pt();
                elec_features.hcalDepth1OverEcal = electron.hcalDepth1OverEcal(); 
                elec_features.hcalDepth2OverEcal = electron.hcalDepth2OverEcal();  
                elec_features.ecalPFClusterIso = electron.ecalPFClusterIso()/electron.pt(); 
                elec_features.hcalPFClusterIso = electron.hcalPFClusterIso()/electron.pt(); 

                


                elec_features.pfSumPhotonEt = electron.pfIsolationVariables().sumPhotonEt/electron.pt(); 
                elec_features.pfSumChargedHadronPt = electron.pfIsolationVariables().sumChargedHadronPt/electron.pt(); 
                elec_features.pfSumNeutralHadronEt = electron.pfIsolationVariables().sumNeutralHadronEt/electron.pt(); 
                elec_features.pfSumPUPt = electron.pfIsolationVariables().sumPUPt/electron.pt(); 

                // isolation
                elec_features.dr04TkSumPt = electron.dr04TkSumPt()/electron.pt();
                elec_features.dr04EcalRecHitSumEt = electron.dr04EcalRecHitSumEt()/electron.pt(); 
                elec_features.dr04HcalDepth1TowerSumEt = electron.dr04HcalDepth1TowerSumEt()/electron.pt(); 
                elec_features.dr04HcalDepth1TowerSumEtBc = electron.dr04HcalDepth1TowerSumEtBc()/electron.pt(); 
                elec_features.dr04HcalDepth2TowerSumEt = electron.dr04HcalDepth2TowerSumEt()/electron.pt(); 
                elec_features.dr04HcalDepth2TowerSumEtBc = electron.dr04HcalDepth2TowerSumEtBc()/electron.pt();

                elec_features.dr04HcalTowerSumEt = electron.dr04HcalTowerSumEt()/electron.pt();
                elec_features.dr04HcalTowerSumEtBc = electron.dr04HcalTowerSumEtBc()/electron.pt();

                features.elec_features.emplace_back(elec_features);
            }
            
            std::stable_sort(features.elec_features.begin(),features.elec_features.end());


            features.cpf_features.emplace_back(cpf_features);

        } //end loop over charged consistuents
        
    
        std::stable_sort(features.cpf_features.begin(),features.cpf_features.end());
        
        
        // Fill neutral hadron info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()!=0)
            {
                continue;
            }
            llpdnnx::NeutralCandidateFeatures npf_features;

            npf_features.ptrel = constituent->pt()/uncorrectedPt;
            npf_features.deta = constituent->eta()-jet.eta();
            npf_features.dphi = reco::deltaPhi(constituent->phi(),jet.phi());
            
            npf_features.px = constituent->px();
            npf_features.py = constituent->py();
            npf_features.pz = constituent->pz();
            
            npf_features.puppi_weight = constituent->puppiWeight();
            npf_features.deltaR = reco::deltaR(*constituent,jet);
            npf_features.isGamma = abs(constituent->pdgId())==22;
            npf_features.hcal_fraction = constituent->hcalFraction();

            npf_features.drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                npf_features.drminsv = std::min(npf_features.drminsv,dR);
            }

            if (jet.mass()<1e-10) 
            {
                npf_features.relmassdrop = -1;
            }
            else
            {
                npf_features.relmassdrop = (jet.p4()- constituent->p4()).mass()/jet.mass();
            }
            features.npf_features.emplace_back(npf_features);
            
        }
        std::stable_sort(features.npf_features.begin(),features.npf_features.end());

        float jetRchg(-1), jetRntr(-1);
        if (features.cpf_features.size() > 0){
            jetRchg = features.cpf_features.at(0).ptrel;
        }
        
        if (features.npf_features.size() > 0){
            jetRntr = features.npf_features.at(0).ptrel;
        }

        float jetR = std::max(jetRchg, jetRntr);

        features.jet_features.jetRchg = jetRchg;
        features.jet_features.jetR = jetR;

        int beta = 0;
        int frac01 = 0;
        int frac02 = 0;
        int frac03 = 0;
        int frac04 = 0;
        float dR2Mean = 0;
        float pt2Sum = 0;


        for (size_t i = 0; i < features.cpf_features.size(); i++){
            llpdnnx::ChargedCandidateFeatures cpf = features.cpf_features.at(i);
            beta += cpf.fromPV;
            dR2Mean += (cpf.ptrel*cpf.trackDeltaR) * (cpf.ptrel*cpf.trackDeltaR);
            pt2Sum += (cpf.ptrel) * (cpf.ptrel);
            if (cpf.trackDeltaR < 0.1) frac01 ++;
            else if (cpf.trackDeltaR < 0.2) frac02 ++;
            else if (cpf.trackDeltaR < 0.3) frac03 ++;
            else if (cpf.trackDeltaR < 0.4) frac04 ++;
        }

        if (features.cpf_features.size() > 0)
        {
            features.jet_features.beta = (float)beta/(float)features.cpf_features.size();

        }



        for (size_t i = 0; i < features.npf_features.size(); i++){
            llpdnnx::NeutralCandidateFeatures npf = features.npf_features.at(i);
            dR2Mean += (npf.ptrel*npf.deltaR) * (npf.ptrel*npf.deltaR);
            pt2Sum += (npf.ptrel) * (npf.ptrel);
            if (npf.deltaR < 0.1) frac01 ++;
            else if (npf.deltaR < 0.2) frac02 ++;
            else if (npf.deltaR < 0.3) frac03 ++;
            else if (npf.deltaR < 0.4) frac04 ++;
        }

        int nCandidates = features.cpf_features.size()+features.npf_features.size();

        if (nCandidates > 0) {
            features.jet_features.frac01 = 1.*frac01/nCandidates;
            features.jet_features.frac02 = 1.*frac02/nCandidates;
            features.jet_features.frac03 = 1.*frac03/nCandidates;
            features.jet_features.frac04 = 1.*frac04/nCandidates;
            features.jet_features.dR2Mean = dR2Mean/pt2Sum;
        }
        
        features.jet_features.numberCpf = features.cpf_features.size();
        features.jet_features.numberNpf = features.npf_features.size();
        features.jet_features.numberSv = features.sv_features.size();
        features.jet_features.numberSvAdapted = features.sv_adapted_features.size();
        features.jet_features.numberMuon = features.mu_features.size();
        features.jet_features.numberElectron = features.elec_features.size();
        
        output_tag_infos->emplace_back(features, jet_ref);
    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void XTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}
void XTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(XTagInfoProducer);
