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

#include "TVector3.h"


class XTagInfoProducer : public edm::stream::EDProducer<> {
public:
    explicit XTagInfoProducer(const edm::ParameterSet&);
    ~XTagInfoProducer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
    edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> sv_token_;
    edm::EDGetTokenT<edm::View<reco::ShallowTagInfo>> shallow_tag_info_token_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
    typedef std::vector<reco::XTagInfo> XTagInfoCollection;

};

XTagInfoProducer::XTagInfoProducer(const edm::ParameterSet& iConfig) :
jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
vtx_token_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
sv_token_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
shallow_tag_info_token_(
    consumes<edm::View<reco::ShallowTagInfo>>(iConfig.getParameter<edm::InputTag>("shallow_tag_infos")))
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
    edm::Handle<edm::View<reco::Jet>> jets;
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

    for (std::size_t ijet = 0; ijet < jets->size(); ijet++) {
        // create data containing structure
        llpdnnx::XTagFeatures features;
        const pat::Jet& jet = jets->at(ijet); 
        edm::RefToBase<reco::Jet> jet_ref(jets, ijet);



        // Start with global jet features
        features.jet_features.pt = std::log10(jet.pt());  // uncorrected
        features.jet_features.eta = jet.eta();
        features.jet_features.mass = jet.mass();
        features.jet_features.energy = jet.energy();

        features.npv = vtxs->size();

        // Add CSV variables
        const edm::View<reco::ShallowTagInfo>& taginfos = *shallow_tag_infos;
        edm::Ptr<reco::ShallowTagInfo> match;
        // Try first by 'same index'
        if ((ijet < taginfos.size()) && (taginfos[ijet].jet() == jet_ref)) {
        match = taginfos.ptrAt(ijet);
        } else {
        // otherwise fail back to a simple search
        for (auto itTI = taginfos.begin(), edTI = taginfos.end(); itTI != edTI; ++itTI) {
            if (itTI->jet() == jet_ref) {
            match = taginfos.ptrAt(itTI - taginfos.begin());
            break;
                }
            }
        }
        reco::ShallowTagInfo tag_info;
        if (match.isNonnull()) {
            tag_info = *match;
        }  // will be default values otherwise

        reco::TaggingVariableList vars = tag_info.taggingVariables();
        features.tag_info_features.csv_trackSumJetEtRatio = vars.get(reco::btau::trackSumJetEtRatio, -1);
        features.tag_info_features.csv_trackSumJetDeltaR = vars.get(reco::btau::trackSumJetDeltaR, -1);
        features.tag_info_features.csv_vertexCategory = vars.get(reco::btau::vertexCategory, -1);
        features.tag_info_features.csv_trackSip2dValAboveCharm = vars.get(reco::btau::trackSip2dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip2dSigAboveCharm = vars.get(reco::btau::trackSip2dSigAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dValAboveCharm = vars.get(reco::btau::trackSip3dValAboveCharm, -1);
        features.tag_info_features.csv_trackSip3dSigAboveCharm = vars.get(reco::btau::trackSip3dSigAboveCharm, -1);
        features.tag_info_features.csv_jetNTracksEtaRel = vars.get(reco::btau::jetNTracksEtaRel, -1);
        features.tag_info_features.csv_jetNSelectedTracks = vars.get(reco::btau::jetNSelectedTracks, -1);
                

        // fill features from secondary vertices
        for (unsigned int isv = 0; isv < svs->size(); ++isv)
        {
            const reco::VertexCompositePtrCandidate& sv = svs->at(isv);
            if (reco::deltaR(sv,jet)>0.4)
            {
                continue;
            }
            llpdnnx::SecondaryVertexFeatures sv_features;
        
            sv_features.sv_pt = std::log10(sv.pt()); 
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

        std::stable_sort(features.sv_features.begin(),features.sv_features.end(),[&pv](const auto& d1, const auto& d2)
        {
            return d1.sv_dxysig>d2.sv_dxysig; //sort decreasing
        });


        // Fill track info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter){
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
                
            cpf_features.cpf_vertex_association = constituent->pvAssociationQuality();
            cpf_features.cpf_fromPV = constituent->fromPV();
            cpf_features.cpf_puppi_weight = constituent->puppiWeight();
            cpf_features.cpf_track_chi2 = constituent->pseudoTrack().chi2();
            cpf_features.cpf_track_ndof = constituent->pseudoTrack().ndof();
            cpf_features.cpf_track_quality = constituent->pseudoTrack().qualityMask();

            //if (jet.mass()<1e-10) cpf_features.cpf_jetmassdroprel = 0;
            //else cpf_features.jetmassdroprel = (jet.p4()-constituent->p4()).mass()/jet.mass();
            
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
            //cpf_features.cpf_trackPtRatio=cpf_features.cpf_trackPtRel / trackMag;
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
                
            float sumPt = 0.;
            for (unsigned int jdaughter = 0; jdaughter < jet.numberOfDaughters(); ++jdaughter)
            {
                const pat::PackedCandidate* other = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(jdaughter));
                if (other and other!=constituent and reco::deltaR(*other,*constituent)<0.1)
                {
                    sumPt += other->pt();
                }
            }
            //cpf_features.cpf_relIso01 = sumPt/constituent->pt();
            
            //cpf_features.cpf_lostInnerHits = constituent->lostInnerHits(); //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_4/doc/html/d8/d79/classpat_1_1PackedCandidate.html#ab9ef9a12f92e02fa61653ba77ee34274

            //int pdgId = std::abs(constituent->pdgId());
            //cpf_features.cpf_isLepton = pdgId==11 or pdgId==13;
        
            features.cpf_features.emplace_back(cpf_features);
        }
        std::stable_sort(features.cpf_features.begin(),features.cpf_features.end(),[](const auto& d1, const auto& d2)
            {
                if (!std::isnan(d1.cpf_trackSip2dSig) and !std::isinf(d1.cpf_trackSip2dSig))
                {
                    if (!std::isnan(d2.cpf_trackSip2dSig) and !std::isinf(d2.cpf_trackSip2dSig))
                    {
                        if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return std::fabs(d1.cpf_trackSip2dSig)>std::fabs(d2.cpf_trackSip2dSig); //sort decreasing
                        }
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.cpf_trackSip2dSig) and !std::isinf(d2.cpf_trackSip2dSig))
                {
                    return false;
                }
                if (!std::isnan(d1.cpf_drminsv) and !std::isinf(d1.cpf_drminsv))
                {
                    if (!std::isnan(d2.cpf_drminsv) and !std::isinf(d2.cpf_drminsv))
                    {
                        if (std::fabs(d1.cpf_drminsv-d2.cpf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return d1.cpf_drminsv<d2.cpf_drminsv; //sort increasing
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else if (!std::isnan(d2.cpf_drminsv) and !std::isinf(d2.cpf_drminsv))
                {
                    return true;
                }
                
                if (!std::isnan(d1.cpf_ptrel) and !std::isinf(d1.cpf_ptrel))
                {
                    if (!std::isnan(d2.cpf_ptrel) and !std::isinf(d2.cpf_ptrel))
                    {
                        return d1.cpf_ptrel>d2.cpf_ptrel; //sort decreasing
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.cpf_ptrel) and !std::isinf(d2.cpf_ptrel))
                {
                    return false;
                }
                return false;
            });
        // Fill neutral hadron info
        for (unsigned int idaughter = 0; idaughter < jet.numberOfDaughters(); ++idaughter){
            const pat::PackedCandidate* constituent = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(idaughter));
            if ((not constituent) or constituent->charge()!=0)
            {
                continue;
            }
            llpdnnx::NeutralCandidateFeatures npf_features;
        
            npf_features.npf_ptrel = constituent->pt()/jet.pt();
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
            
            //if (jet.mass()<1e-10) npf_features.npf._etmassdroprel = -9;
            //else npf_features.npf_jetmassdroprel = (jet.p4()- constituent->p4()).mass()/jet.mass();

            /*
            float sumPt = 0.;
            for (unsigned int jdaughter = 0; jdaughter < jet.numberOfDaughters(); ++jdaughter)
            {
                const pat::PackedCandidate* other = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(jdaughter));
                if (other and other!=constituent and reco::deltaR(*other,*constituent)<0.1)
                {
                    sumPt += other->pt();
                }
            }
            npf_features.npf_relIso01 = sumPt/constituent->pt();
            */
            features.npf_features.emplace_back(npf_features);
        }
        std::stable_sort(features.npf_features.begin(),features.npf_features.end(),[](const auto& d1, const auto& d2)
            {
                if (!std::isnan(d1.npf_drminsv) and !std::isinf(d1.npf_drminsv))
                {
                    if (!std::isnan(d2.npf_drminsv) and !std::isinf(d2.npf_drminsv))
                    {
                        if (std::fabs(d1.npf_drminsv-d2.npf_drminsv)>std::numeric_limits<float>::epsilon())
                        {
                            return d1.npf_drminsv<d2.npf_drminsv; //sort increasing
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else if (!std::isnan(d2.npf_drminsv) and !std::isinf(d2.npf_drminsv))
                {
                    return true;
                }
                    
                if (!std::isnan(d1.npf_ptrel) and !std::isinf(d1.npf_ptrel))
                    {
                    if (!std::isnan(d2.npf_ptrel) and !std::isinf(d2.npf_ptrel))
                    {
                        return d1.npf_ptrel>d2.npf_ptrel; //sort decreasing
                    }
                    else
                    {
                        return true;
                    }
                }
                else if (!std::isnan(d2.npf_ptrel) and !std::isinf(d2.npf_ptrel))
                {
                    return false;
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
