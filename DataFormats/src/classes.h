#include "DataFormats/Common/interface/Wrapper.h"
#include "LLPReco/DataFormats/interface/XTagFeatures.h"
#include "LLPReco/DataFormats/interface/XTagInfo.h"
#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"
#include "LLPReco/DataFormats/interface/LLPLabel.h"

namespace {

    struct dictionary
    {
        std::vector<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>> dummy0;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>>> dummy1;
        reco::FeaturesTagInfo<llpdnnx::XTagFeatures> dummy2;
        edm::Wrapper<reco::FeaturesTagInfo<llpdnnx::XTagFeatures>> dummy3;


        llpdnnx::XTagFeatures dummy4;
        llpdnnx::JetFeatures dummy5;
        llpdnnx::SecondaryVertexFeatures dummy6;
        llpdnnx::ChargedCandidateFeatures dummy7;
        llpdnnx::NeutralCandidateFeatures dummy8;
        

        llpdnnx::DisplacedGenVertexCollection dummy9;
        edm::Wrapper<llpdnnx::DisplacedGenVertexCollection> dummy10;
        
        edm::Ptr<llpdnnx::DisplacedGenVertex> dummy11;
        edm::Wrapper<edm::Ptr<llpdnnx::DisplacedGenVertex>> dummy12;
        
        edm::Ptr<llpdnnx::DisplacedGenVertexCollection> dummy13;
        edm::Wrapper<edm::Ptr<llpdnnx::DisplacedGenVertexCollection>> dummy14;

        edm::PtrVector<llpdnnx::DisplacedGenVertexCollection> dumm15;
        edm::Wrapper<edm::PtrVector<llpdnnx::DisplacedGenVertexCollection>> dummy16;
        
        edm::PtrVector<reco::GenParticle> dummy17;
        edm::Wrapper<edm::PtrVector<reco::GenParticle>> dummy18;
        
        llpdnnx::LLPLabel dummy19;
        std::vector<reco::FeaturesTagInfo<llpdnnx::LLPLabel>> dummy20;
        edm::Wrapper<std::vector<reco::FeaturesTagInfo<llpdnnx::LLPLabel>>> dummy21;

    };

}
