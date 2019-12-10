#include "DataFormats/Common/interface/Wrapper.h"
#include "LLPReco/DataFormats/interface/XTagFeatures.h"
#include "LLPReco/DataFormats/interface/XTagInfo.h"

namespace {

    struct disctionary_xtag
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

    };

}
