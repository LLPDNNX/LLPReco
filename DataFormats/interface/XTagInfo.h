#ifndef LLPReco_DataFormats_DeepFlavourTagInfo_h
#define LLPReco_DataFormats_DeepFlavourTagInfo_h

#include "LLPReco/DataFormats/interface/XTagFeatures.h"
#include "DataFormats/BTauReco/interface/FeaturesTagInfo.h"

namespace reco {

    typedef  FeaturesTagInfo<llpdnnx::XTagFeatures> XTagInfo;

    DECLARE_EDM_REFS( XTagInfo )

}

#endif

