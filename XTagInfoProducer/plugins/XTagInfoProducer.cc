/**\class XTagInfoProducer XTagInfoProducer.cc RecoLLP/XTagInfoProducer/plugins/XTagInfoProducer.cc

Description: Produces and fill in LLPDNNX features

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Vilius Cepaitis, Matthias Komm
//         Created:  Fri, 06 Dec 2019 12:09:48 GMT

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

    edm::Handle< pat::MuonCollection > muons;
    iEvent.getByToken(muonsMiniAODToken_, muons);

    edm::Handle<pat::ElectronCollection > electrons;
    iEvent.getByToken(electronsMiniAODToken_, electrons);

    std::unordered_map< reco::CandidatePtr, const  pat::Muon*, CandidateHash> muonMap ;
    std::unordered_map< reco::CandidatePtr, const  pat::Electron*, CandidateHash> electronMap ;

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
    
        std::unordered_set<reco::CandidatePtr, CandidateHash> jetConsituentSet;
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            jetConsituentSet.insert(jet.daughterPtr(idaughter));
        }
    
        // create data containing structure
        llpdnnx::XTagFeatures features;
        

        // Start with global jet features
        features.jet_features.pt = jet.pt();  // uncorrected
        features.jet_features.eta = jet.eta();
        features.jet_features.mass = jet.mass();
        features.jet_features.energy = jet.energy();
        
        features.jet_features.n60 = jet.n60();
        features.jet_features.n90 = jet.n90();
        
        features.jet_features.chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
        features.jet_features.chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
        features.jet_features.chargedMuEnergyFraction = jet.chargedMuEnergyFraction();
        features.jet_features.electronEnergyFraction = jet.electronEnergyFraction();


        features.jet_features.jetIdx = jet_ref.key();

        features.npv = vtxs->size();
        
        llpdnnx::JetSubstructure jetSubstructure(jet);

        features.jet_features.tau1 = jetSubstructure.nSubjettiness(1);
        features.jet_features.tau2 = jetSubstructure.nSubjettiness(2);
        features.jet_features.tau3 = jetSubstructure.nSubjettiness(3);
        

        features.jet_features.relMassDropMassAK = jetSubstructure.massDropMass(llpdnnx::JetSubstructure::ClusterType::AK)/features.jet_features.mass;
        features.jet_features.relMassDropMassCA = jetSubstructure.massDropMass(llpdnnx::JetSubstructure::ClusterType::CA)/features.jet_features.mass;
        features.jet_features.relSoftDropMassAK = jetSubstructure.softDropMass(llpdnnx::JetSubstructure::ClusterType::AK)/features.jet_features.mass;
        features.jet_features.relSoftDropMassCA = jetSubstructure.softDropMass(llpdnnx::JetSubstructure::ClusterType::CA)/features.jet_features.mass;
       
        auto eventShapes = jetSubstructure.eventShapeVariables();
        features.jet_features.thrust = jetSubstructure.thrust();
        features.jet_features.sphericity = eventShapes.sphericity();
        features.jet_features.circularity = eventShapes.circularity();
        features.jet_features.isotropy = eventShapes.isotropy();
        features.jet_features.eventShapeC = eventShapes.C();
        features.jet_features.eventShapeD = eventShapes.D();
        

        // Add CSV variables
        const edm::View<reco::ShallowTagInfo>& taginfos = *shallow_tag_infos;
        edm::Ptr<reco::ShallowTagInfo> match;
        // Try first by 'same index'
        if ((ijet < taginfos.size()) && (taginfos[ijet].jet() == jet_ref))
        {
            match = taginfos.ptrAt(ijet);
        } else {
            // otherwise fail back to a simple search
            for (auto itTI = taginfos.begin(), edTI = taginfos.end(); itTI != edTI; ++itTI)
            {
                if (itTI->jet() == jet_ref) {
                    match = taginfos.ptrAt(itTI - taginfos.begin());
                    //break;
                }
            }
        }
        reco::ShallowTagInfo tag_info;
        if (match.isNonnull()) {
            tag_info = *match;
        }  // will be default values otherwise

        reco::TaggingVariableList vars = tag_info.taggingVariables();
        features.tag_info_features.csv_jetIdx = jet_ref.key();
        features.tag_info_features.csv_trackSumJetEtRatio = vars.get(reco::btau::trackSumJetEtRatio, -1);
        features.tag_info_features.csv_trackSumJetDeltaR = vars.get(reco::btau::trackSumJetDeltaR, -1);
        features.tag_info_features.csv_vertexCategory = vars.get(reco::btau::vertexCategory, -1);
        features.tag_info_features.csv_trackSip2dValAboveCharm = vars.get(reco::btau::trackSip2dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip2dSigAboveCharm = vars.get(reco::btau::trackSip2dSigAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dValAboveCharm = vars.get(reco::btau::trackSip3dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dSigAboveCharm = vars.get(reco::btau::trackSip3dSigAboveCharm, -1);
        features.tag_info_features.csv_jetNTracksEtaRel = vars.get(reco::btau::jetNTracksEtaRel, -1);
        features.tag_info_features.csv_jetNSelectedTracks = vars.get(reco::btau::jetNSelectedTracks, -1);


        std::unordered_set<reco::CandidatePtr, CandidateHash> candidatesMatchedToSV;
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
                if (jetConsituentSet.find(candidateFromVertex)!=jetConsituentSet.end())
                {
                    candidatesMatchedToSV.insert(candidateFromVertex);
                    matchingTrack = true;
                }
            }
            if (not matchingTrack) continue;
            
            llpdnnx::SecondaryVertexFeatures sv_features;

            sv_features.sv_jetIdx = jet_ref.key();
            sv_features.sv_pt = sv.pt(); 
            sv_features.sv_deltaR = reco::deltaR(sv,jet);
            sv_features.sv_mass = sv.mass();
            sv_features.sv_ntracks = sv.numberOfDaughters();
            sv_features.sv_chi2 = sv.vertexChi2();
            sv_features.sv_ndof = sv.vertexNdof();


            reco::Vertex::CovarianceMatrix covsv; 
            sv.fillVertexCovariance(covsv);
            reco::Vertex svtx(sv.vertex(), covsv);

            VertexDistanceXY distXY;
            Measurement1D distanceXY = distXY.distance(svtx, pv);
            sv_features.sv_dxy = distanceXY.value();
            sv_features.sv_dxysig = distanceXY.value()/distanceXY.error();

            VertexDistance3D dist3D;
            Measurement1D distance3D = dist3D.distance(svtx, pv);
            sv_features.sv_d3d = distance3D.value();
            sv_features.sv_d3dsig = distance3D.value()/distance3D.error();

            reco::Candidate::Vector distance(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
            sv_features.sv_costhetasvpv = sv.momentum().Unit().Dot(distance.Unit());
            sv_features.sv_enratio = sv.energy()/jet.pt();


            features.sv_features.emplace_back(sv_features);
        }

        std::stable_sort(features.sv_features.begin(),features.sv_features.end(),[](const auto& d1, const auto& d2)
        {
            return d1.sv_dxysig>d2.sv_dxysig; //sort decreasing
        });


        // Fill cpf info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()==0 or (not constituent->hasTrackDetails()))
            {
                continue;
            }

            llpdnnx::ChargedCandidateFeatures cpf_features;

            cpf_features.cpf_ptrel = constituent->pt()/jet.pt();

            cpf_features.cpf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                cpf_features.cpf_drminsv = std::min(cpf_features.cpf_drminsv,dR);
            }

            cpf_features.cpf_jetIdx = jet_ref.key();
            cpf_features.cpf_vertex_association = constituent->pvAssociationQuality();
            cpf_features.cpf_fromPV = constituent->fromPV();
            cpf_features.cpf_puppi_weight = constituent->puppiWeight();
            cpf_features.cpf_track_chi2 = constituent->pseudoTrack().chi2();
            cpf_features.cpf_track_ndof = constituent->pseudoTrack().ndof();
            cpf_features.cpf_track_quality = constituent->pseudoTrack().qualityMask();

            if (jet.mass()<1e-10)
            {
                cpf_features.cpf_relmassdrop = -1;
            }
            else
            {
                cpf_features.cpf_relmassdrop = (jet.p4()-constituent->p4()).mass()/jet.mass();
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

            cpf_features.cpf_trackEtaRel=reco::btau::etaRel(jetDir, trackMom);
            cpf_features.cpf_trackPtRel=trackMom3.Perp(jetDir3);
            cpf_features.cpf_trackPPar=jetDir.Dot(trackMom);
            cpf_features.cpf_trackDeltaR=reco::deltaR(trackMom, jetDir);
            cpf_features.cpf_trackPtRatio=cpf_features.cpf_trackPtRel / trackMag;
            cpf_features.cpf_trackPParRatio=cpf_features.cpf_trackPPar / trackMag;

            cpf_features.cpf_trackSip2dVal=std::abs(meas_ip2d.value());
            cpf_features.cpf_trackSip2dSig=std::abs(meas_ip2d.significance());
            cpf_features.cpf_trackSip3dVal=std::abs(meas_ip3d.value());
            cpf_features.cpf_trackSip3dSig=std::abs(meas_ip3d.significance());
            if (std::isnan(cpf_features.cpf_trackSip2dSig) || std::isnan(cpf_features.cpf_trackSip3dSig))
            {
                cpf_features.cpf_trackSip2dSig=-1.;
                cpf_features.cpf_trackSip3dSig=-1.;
            }

            cpf_features.cpf_trackJetDistVal = jetdist.value();
            cpf_features.cpf_trackJetDistSig = jetdist.significance();

            cpf_features.cpf_matchedMuon = 0;
            cpf_features.cpf_matchedElectron = 0;
            
            if (candidatesMatchedToSV.find(jet.daughterPtr(idaughter))!=candidatesMatchedToSV.end())
            {
                cpf_features.cpf_matchedSV = 1;
            }
            else
            {
                cpf_features.cpf_matchedSV = 0;
            }

            
            
            //find matching muons
            llpdnnx::MuonCandidateFeatures mu_features; 
            auto findMuon = muonMap.find(jet.daughterPtr(idaughter));  
            if (findMuon!=muonMap.end())
            {
                const pat::Muon & muon = *findMuon->second;

                if (not muon.isGlobalMuon()) continue ;
                cpf_features.cpf_matchedMuon = 1;

                mu_features.mu_isGlobal = muon.isGlobalMuon() ;                                   
                mu_features.mu_isTight = muon.isTightMuon(pv);                                     
                mu_features.mu_isMedium = muon.isMediumMuon();       
                mu_features.mu_isLoose = muon.isLooseMuon() ; 
                mu_features.mu_isStandAlone = muon.isStandAloneMuon() ; 

                mu_features.mu_ptrel = muon.pt()/jet.pt() ; 
                mu_features.mu_eta = muon.eta();                                                 
                mu_features.mu_phi = muon.phi();                                                 
                mu_features.mu_charge = muon.charge();        
                mu_features.mu_energy = muon.energy()/muon.pt();                                           
                mu_features.mu_et = muon.et();   
                mu_features.mu_jetDeltaR = reco::deltaR(muon ,jet) ; 
                mu_features.mu_numberOfMatchedStations = muon.numberOfMatchedStations();

                mu_features.mu_2dIp = muon.dB() ; 
                mu_features.mu_2dIpSig = muon.dB()/muon.edB() ; 
                mu_features.mu_3dIp = muon.dB(pat::Muon::PV3D) ; 
                mu_features.mu_3dIpSig = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D) ;


                reco::Candidate::Vector muonMom = muon.bestTrack()->momentum();

                mu_features.mu_EtaRel =reco::btau::etaRel(jetDir, muonMom);
                mu_features.mu_dxy = muon.bestTrack()->dxy(pv.position());
                mu_features.mu_dxyError = muon.bestTrack()->dxyError() ; 
                mu_features.mu_dxySig = muon.bestTrack()->dxy(pv.position())/muon.bestTrack()->dxyError(); 
                mu_features.mu_dz = muon.bestTrack()->dz(pv.position()) ; 
                mu_features.mu_dzError = muon.bestTrack()->dzError() ;
                mu_features.mu_numberOfValidPixelHits = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
                mu_features.mu_numberOfpixelLayersWithMeasurement = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement() ;
                mu_features.mu_numberOfstripLayersWithMeasurement = muon.bestTrack()->hitPattern().stripLayersWithMeasurement() ;


                mu_features.mu_chi2 = muon.bestTrack()->chi2() ;  
                mu_features.mu_ndof = muon.bestTrack()->ndof() ;


                mu_features.mu_caloIso =  muon.caloIso()/muon.pt() ; 
                mu_features.mu_ecalIso =  muon.ecalIso()/muon.pt() ; 
                mu_features.mu_hcalIso =  muon.hcalIso()/muon.pt() ;     


                mu_features.mu_sumPfChHadronPt  = muon.pfIsolationR04().sumChargedHadronPt/muon.pt();
                mu_features.mu_sumPfNeuHadronEt  = muon.pfIsolationR04().sumNeutralHadronEt/muon.pt();
                mu_features.mu_Pfpileup  = muon.pfIsolationR04().sumPUPt/muon.pt();
                mu_features.mu_sumPfPhotonEt = muon.pfIsolationR04().sumPhotonEt/muon.pt();



                mu_features.mu_sumPfChHadronPt03  = muon.pfIsolationR03().sumChargedHadronPt/muon.pt();
                mu_features.mu_sumPfNeuHadronEt03  = muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
                mu_features.mu_Pfpileup03  = muon.pfIsolationR03().sumPUPt/muon.pt();
                mu_features.mu_sumPfPhotonEt03 = muon.pfIsolationR03().sumPhotonEt/muon.pt();       


                mu_features.mu_timeAtIpInOut = muon.time().timeAtIpInOut ;  
                mu_features.mu_timeAtIpInOutErr = muon.time().timeAtIpInOutErr ; 
                mu_features.mu_timeAtIpOutIn = muon.time().timeAtIpOutIn ;  

                features.mu_features.emplace_back(mu_features);
            }

            std::stable_sort(features.mu_features.begin(),features.mu_features.end(),[](const auto& d1, const auto& d2)
            {
                if (d1.mu_2dIpSig>0 and d2.mu_2dIpSig>0)
                {
                    if (std::fabs(d1.mu_2dIpSig-d2.mu_2dIpSig)>std::numeric_limits<float>::epsilon())
                    {
                        return std::fabs(d1.mu_2dIpSig)>std::fabs(d2.mu_2dIpSig); //sort decreasing
                    }
                }
                return d1.mu_ptrel>d2.mu_ptrel; //sort decreasing
            });


            //find matching electrons
            llpdnnx::ElectronCandidateFeatures elec_features;
            auto findElectron = electronMap.find(jet.daughterPtr(idaughter));  
            if(findElectron!=electronMap.end())
            {
                const pat::Electron & electron = *findElectron->second;

                cpf_features.cpf_matchedElectron = 1;
		elec_features.electronFeatures(electron , jet , pv );

                features.elec_features.emplace_back(elec_features);
            }
            
            std::stable_sort(features.elec_features.begin(),features.elec_features.end(),[](const auto& d1, const auto& d2)


            {
                if (d1.elec_2dIPSig>0 and d2.elec_2dIPSig>0)
                {
                    if (std::fabs(d1.elec_2dIPSig-d2.elec_2dIPSig)>std::numeric_limits<float>::epsilon())
                    {
                        return std::fabs(d1.elec_2dIPSig)>std::fabs(d2.elec_2dIPSig); //sort decreasing
                    }
                }
                return d1.elec_ptrel>d2.elec_ptrel; //sort decreasing
            });
            
            
            features.cpf_features.emplace_back(cpf_features);

        } //end loop over charged consistuents
        
    
        std::stable_sort(features.cpf_features.begin(),features.cpf_features.end(),[](const auto& d1, const auto& d2)
        {
            if (d1.cpf_trackSip2dSig>0 and d2.cpf_trackSip2dSig>0)
            {
                return std::fabs(d1.cpf_trackSip2dSig)>std::fabs(d2.cpf_trackSip2dSig); //sort decreasing
            }
            else if (d1.cpf_trackSip2dSig<0 and d2.cpf_trackSip2dSig>0)
            {
                return false;
            }
            else if (d1.cpf_trackSip2dSig>0 and d2.cpf_trackSip2dSig<0)
            {
                return true;
            }
            else if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
            {
                return d1.cpf_drminsv<d2.cpf_drminsv; //sort increasing
            }
            else
            {
                return d1.cpf_ptrel>d2.cpf_ptrel;  //sort decreasing
            }
            
            return false;
        });
        
        
        // Fill neutral hadron info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter)
        {
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()!=0)
            {
                continue;
            }
            llpdnnx::NeutralCandidateFeatures npf_features;

            npf_features.npf_ptrel = constituent->pt()/jet.pt();
            npf_features.npf_jetIdx = jet_ref.key();
            npf_features.npf_puppi_weight = constituent->puppiWeight();
            npf_features.npf_deltaR = reco::deltaR(*constituent,jet);
            npf_features.npf_isGamma = fabs(constituent->pdgId())==22;
            npf_features.npf_hcal_fraction = constituent->hcalFraction();

            npf_features.npf_drminsv = 0.4;
            for (const auto& sv: *svs.product())
            {
                float dR = reco::deltaR(sv,*constituent);
                npf_features.npf_drminsv = std::min(npf_features.npf_drminsv,dR);
            }

            if (jet.mass()<1e-10) 
            {
                npf_features.npf_relmassdrop = -1;
            }
            else
            {
                npf_features.npf_relmassdrop = (jet.p4()- constituent->p4()).mass()/jet.mass();
            }
            features.npf_features.emplace_back(npf_features);
            
        }
        std::stable_sort(features.npf_features.begin(),features.npf_features.end(),[](const auto& d1, const auto& d2)
        {

            if (std::fabs(d1.npf_drminsv-d2.npf_drminsv)>std::numeric_limits<float>::epsilon())
            {
                return d1.npf_drminsv<d2.npf_drminsv; //sort increasing
            }
            else
            {
                return d1.npf_ptrel>d2.npf_ptrel; //sort decreasing
            }
            return false;
        });
        output_tag_infos->emplace_back(features, jet_ref);
    }

    iEvent.put(std::move(output_tag_infos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void XTagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("ak4PFJetsCHS"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
    desc.add<edm::InputTag>("shallow_tag_infos", edm::InputTag("pfDeepCSVTagInfos"));
}
void XTagInfoProducer::endStream() {};

//define this as a plug-in
DEFINE_FWK_MODULE(XTagInfoProducer);
